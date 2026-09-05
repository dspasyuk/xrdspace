#!/usr/bin/env node
// Copyright (c) 2026 Denis Spasyuk. MIT License.
// Build the offline PDB unit-cell / space-group lookup table for `xrdspace
// --valid`.
//
// Downloads the unit-cell parameters and space group of every PDB entry
// directly from RCSB (no mmCIF files involved):
//   1. page the RCSB Search API for all entry ids,
//   2. fetch cell + symmetry in batches of 1000 through the RCSB GraphQL API,
//   3. keep the diffraction entries (a finite unit cell + space-group number),
//   4. store each entry's cell and Niggli-reduced cell, sorted by reduced a.
//
// The result is written to data/pdb-cells.json (git-ignored; ~40 MB for the
// full database). Re-run after the PDB grows to refresh it.
//
// Usage:
//   node scripts/build-pdb-table.js              # full database
//   node scripts/build-pdb-table.js --limit 5000 # first N entries (testing)
//   node scripts/build-pdb-table.js --out out.json

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

import { buildPdbLookupTable } from '../src/pdb-lookup.js';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.join(__dirname, '..');
const DEFAULT_OUT = path.join(ROOT, 'data', 'pdb-cells.json');

const RCSB_SEARCH = 'https://search.rcsb.org/rcsbsearch/v2/query';
const RCSB_GRAPHQL = 'https://data.rcsb.org/graphql';

const ROWS = 10000;          // search API page size (RCSB caps this at 10000)
const CHUNK = 1000;          // GraphQL entry_ids limit
const CONCURRENCY = 8;
const MAX_ATTEMPTS = 3;

async function fetchJson(url, options = {}, timeoutMs = 60000) {
    let lastErr;
    for (let attempt = 1; attempt <= MAX_ATTEMPTS; attempt++) {
        const ctrl = new AbortController();
        const timer = setTimeout(() => ctrl.abort(), timeoutMs);
        try {
            const res = await fetch(url, { ...options, signal: ctrl.signal });
            if (res.status === 204) return null;
            if (!res.ok) throw new Error(`HTTP ${res.status} for ${url}`);
            return await res.json();
        } catch (e) {
            lastErr = e;
            if (attempt < MAX_ATTEMPTS) await new Promise(r => setTimeout(r, 500 * attempt));
        } finally {
            clearTimeout(timer);
        }
    }
    throw lastErr;
}

async function runConcurrent(items, concurrency, fn) {
    const out = new Array(items.length);
    let next = 0;
    const worker = async () => {
        while (next < items.length) {
            const i = next++;
            try { out[i] = { ok: true, value: await fn(items[i], i) }; }
            catch (e) { out[i] = { ok: false, error: e.message }; }
        }
    };
    await Promise.all(Array.from({ length: Math.min(concurrency, items.length) }, worker));
    return out;
}

// 1. All PDB entry ids via the search API (paginated).
async function collectIds(progress) {
    const ids = [];
    let start = 0;
    const body = {
        query: { type: 'terminal', service: 'full_text', parameters: { value: '*' } },
        return_type: 'entry',
        request_options: { paginate: { start: 0, rows: ROWS } },
    };
    for (;;) {
        body.request_options.paginate.start = start;
        const j = await fetchJson(RCSB_SEARCH, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(body),
        });
        if (j === null) break;
        const page = (j.result_set || []).map(r => r.identifier);
        if (!page.length) break;
        ids.push(...page);
        start += ROWS;
        if (start >= (j.total_count != null ? j.total_count : start)) break;
        if (progress) progress(`collected ${ids.length} entry ids`);
    }
    return ids;
}

// 2. Cell + symmetry for a batch of ids via GraphQL.
function graphqlQuery(ids) {
    const fields = `rcsb_id cell { length_a length_b length_c angle_alpha angle_beta angle_gamma } `
        + `symmetry { space_group_name_H_M Int_Tables_number }`;
    return `{ entries(entry_ids: [${ids.map(i => JSON.stringify(i)).join(',')}]) { ${fields} } }`;
}

async function fetchBatch(ids) {
    const j = await fetchJson(RCSB_GRAPHQL, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ query: graphqlQuery(ids) }),
    });
    return (j && j.data && j.data.entries) || [];
}

async function main() {
    const argv = process.argv.slice(2);
    let out = DEFAULT_OUT;
    let limit = Infinity;
    for (let i = 0; i < argv.length; i++) {
        if (argv[i] === '--out') out = path.resolve(argv[++i]);
        else if (argv[i] === '--limit') limit = parseInt(argv[++i], 10);
    }

    let t0 = Date.now();
    const log = (msg) => fs.writeSync(2, msg + '\n');
    log(`xrdspace: building PDB lookup table`);
    log(`  output      : ${out}`);

    const ids = await collectIds((m) => log(m));
    log(`  total ids   : ${ids.length}`);
    const targets = ids.slice(0, limit);

    // Chunk ids and fetch concurrently.
    const chunks = [];
    for (let i = 0; i < targets.length; i += CHUNK) chunks.push(targets.slice(i, i + CHUNK));
    log(`  fetching ${chunks.length} GraphQL batches (${CHUNK} ids each, concurrency ${CONCURRENCY}) ...`);

    const fetched = await runConcurrent(chunks, CONCURRENCY, (chunk) => fetchBatch(chunk));
    let done = 0;
    const rows = [];
    let errors = 0;
    for (const r of fetched) {
        if (!r.ok) { errors++; continue; }
        for (const e of r.value) {
            const c = e.cell;
            if (!c) continue;
            const sym = e.symmetry || {};
            const cell = {
                a: c.length_a, b: c.length_b, c: c.length_c,
                alpha: c.angle_alpha, beta: c.angle_beta, gamma: c.angle_gamma,
            };
            const sg = sym.Int_Tables_number;
            if (!Number.isFinite(sg) || sg < 1 || sg > 230) continue;
            rows.push({ id: e.rcsb_id, cell, sg, hm: sym.space_group_name_H_M || null });
        }
        done++;
        if (done % 25 === 0 || done === fetched.length) {
            log(`  fetched ${done}/${fetched.length} batches (${rows.length} entries with a cell)`);
        }
    }
    if (errors) log(`  WARNING: ${errors} batch(es) failed`);
    log(`  entries with cell+space group : ${rows.length}`);

    // 3. Reduce + sort.
    log(`  computing Niggli-reduced cells ...`);
    const entries = buildPdbLookupTable(rows);
    log(`  kept ${entries.length} entries`);

    // 4. Write.
    fs.mkdirSync(path.dirname(out), { recursive: true });
    const payload = {
        version: 1,
        source: 'RCSB Protein Data Bank (data.rcsb.org, Search + GraphQL APIs)',
        generated: new Date().toISOString(),
        count: entries.length,
        fields: ['id', 'cell(a,b,c,alpha,beta,gamma)', 'red(niggli-reduced)', 'sg(Int Tables number)', 'hm'],
        entries,
    };
    const json = JSON.stringify(payload);
    fs.writeFileSync(out, json, 'utf8');
    const mb = (json.length / 1024 / 1024).toFixed(1);
    log(`  wrote ${out} (${mb} MB, ${entries.length} entries) in ${((Date.now() - t0) / 1000).toFixed(0)}s`);
}

main();
