#!/usr/bin/env node
// Copyright (c) 2026 Denis Spasyuk. MIT License.
// Head-to-head comparison of xrdspace against Bruker XPREP on the cached COD
// datasets used by tests/xrdspace-cod.js.
//
// For every pick (tests/cod-picks.json):
//   1. read the cached HKLs/cod/<id>.hkl,
//   2. write a minimal SHELX .fcf (list code 3) in tests/xprep-work/,
//   3. run XPREP through scripts/xprep-run.py (pexpect) and read the space
//      group it chooses from the .prp log,
//   4. run xrdspace (analyzeHkl) on the same reflections,
//   5. compare both with each other and with the published COD space group.
//
// Usage:
//   node tests/xrdspace-xprep.js                 # all picks
//   node tests/xrdspace-xprep.js --limit 50      # first N
//   node tests/xrdspace-xprep.js --id 1501632    # one entry
//   node tests/xrdspace-xprep.js --sg 14         # one published space group
//   node tests/xrdspace-xprep.js --concurrency 4
//
// Requires XPREP (XPREP_BIN, default xdsgo/executables/xprep) and python3+pexpect.

import fs from 'node:fs';
import path from 'node:path';
import { spawn } from 'node:child_process';
import { fileURLToPath } from 'node:url';

import { parseHkl } from '../src/hkl-parser.js';
import { analyzeHkl } from '../src/index.js';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.join(__dirname, '..');
const CACHE_DIR = path.join(ROOT, 'HKLs', 'cod');
const WORK_DIR = path.join(__dirname, 'xprep-work');
const DRIVER = path.join(ROOT, 'scripts', 'xprep-run.py');
const REPORT_JSON = path.join(__dirname, 'xrdspace-xprep-report.json');

const picks = JSON.parse(fs.readFileSync(path.join(__dirname, 'cod-picks.json'), 'utf8'));

function cellFromPick(p) {
    const [a, b, c, alpha, beta, gamma] = p.cell;
    return { a, b, c, alpha, beta, gamma };
}

function fmt(x, d = 4) {
    return Number(x).toFixed(d);
}

// Build a minimal SHELX .fcf (HKLF list code 3: h k l F^2 sigma) that XPREP
// accepts. The published unit cell from the pick is written explicitly; the
// reflection list comes from xrdspace's own COD parser so both programs see
// exactly the same data.
function buildFcf(id, cell, reflections) {
    const lines = [];
    lines.push('#');
    lines.push(`# COD ${id} — generated for the xdsgo/xrdspace XPREP comparison`);
    lines.push('#');
    lines.push(`data_${id}`);
    lines.push("_shelx_title ' xrdspace-xprep comparison'");
    lines.push('_shelx_refln_list_code          3');
    lines.push(`_cell_length_a    ${fmt(cell.a)}`);
    lines.push(`_cell_length_b    ${fmt(cell.b)}`);
    lines.push(`_cell_length_c    ${fmt(cell.c)}`);
    lines.push(`_cell_angle_alpha  ${fmt(cell.alpha, 3)}`);
    lines.push(`_cell_angle_beta   ${fmt(cell.beta, 3)}`);
    lines.push(`_cell_angle_gamma  ${fmt(cell.gamma, 3)}`);
    lines.push('');
    lines.push('loop_');
    lines.push(' _refln_index_h');
    lines.push(' _refln_index_k');
    lines.push(' _refln_index_l');
    lines.push(' _refln_F_squared_meas');
    lines.push(' _refln_F_squared_sigma');
    for (const r of reflections) {
        lines.push(` ${r.h} ${r.k} ${r.l} ${r.I.toFixed(2)} ${(r.sig || 0).toFixed(2)}`);
    }
    return lines.join('\n') + '\n';
}

function runXprep(id, fcf) {
    return new Promise(resolve => {
        const child = spawn('python3', [DRIVER, fcf], {
            stdio: ['ignore', 'pipe', 'pipe'],
        });
        let out = '';
        const timer = setTimeout(() => {
            try { child.kill('SIGKILL'); } catch { /* ignore */ }
            resolve({ ok: false, error: 'driver timeout' });
        }, 120000);
        child.stdout.on('data', d => { out += d; });
        child.on('error', e => {
            clearTimeout(timer);
            resolve({ ok: false, error: String(e.message || e) });
        });
        child.on('close', () => {
            clearTimeout(timer);
            const lines = out.trim().split(/\r?\n/).filter(Boolean);
            if (!lines.length) return resolve({ ok: false, error: 'no driver output' });
            try {
                resolve(JSON.parse(lines[lines.length - 1]));
            } catch {
                resolve({ ok: false, error: `bad driver JSON: ${lines[lines.length - 1].slice(0, 120)}` });
            }
        });
    });
}

async function runOne(pick) {
    const cache = path.join(CACHE_DIR, `${pick.id}.hkl`);
    if (!fs.existsSync(cache)) return { pick, skip: 'not cached' };
    const text = fs.readFileSync(cache, 'utf8');

    let reflections = [];
    try {
        reflections = parseHkl(text).reflections;
    } catch (e) {
        return { pick, skip: e.message };
    }
    if (!reflections.length) return { pick, skip: 'no reflections parsed' };

    fs.mkdirSync(WORK_DIR, { recursive: true });
    const fcf = path.join(WORK_DIR, `${pick.id}.fcf`);
    fs.writeFileSync(fcf, buildFcf(pick.id, cellFromPick(pick), reflections));

    const xprep = await runXprep(pick.id, fcf);

    let xrd = null;
    try {
        xrd = await analyzeHkl(text, { cell: cellFromPick(pick) });
    } catch (e) {
        xrd = { ok: false, error: e.message };
    }

    const xrdNum = xrd && xrd.ok && xrd.best ? xrd.best.id : null;
    const xrdHm = xrd && xrd.ok && xrd.best ? xrd.best.hm : null;
    const xrdNear = !!(xrd && xrd.ok && xrd.candidates &&
        xrd.candidates.some(c => c.id === pick.sgNumber && c.violations === 0));
    const xpNum = xprep && xprep.ok ? xprep.sgNumber : null;
    const xpNear = !!(xprep && xprep.ok && xprep.candidates &&
        xprep.candidates.some(c => c.number === pick.sgNumber));

    return {
        pick, reflections: reflections.length,
        xrdNum, xrdHm, xrdNear,
        xprepNum: xpNum, xprepHm: xprep.hm, xprepNear: xpNear,
        xprepError: xprep && !xprep.ok ? xprep.error : null,
    };
}

function classify(r) {
    const pub = r.pick.sgNumber;
    const xpPass = r.xprepNum === pub;
    const xrdPass = r.xrdNum === pub;
    const agree = r.xrdNum !== null && r.xprepNum !== null && r.xrdNum === r.xprepNum;
    return { xpPass, xrdPass, agree };
}

function systemOf(n) {
    if (n <= 2) return 'triclinic';
    if (n <= 15) return 'monoclinic';
    if (n <= 74) return 'orthorhombic';
    if (n <= 142) return 'tetragonal';
    if (n <= 167) return 'trigonal';
    if (n <= 194) return 'hexagonal';
    return 'cubic';
}

async function runAll(targets, concurrency) {
    const results = new Array(targets.length);
    let next = 0, done = 0;
    const t0 = Date.now();
    const worker = async () => {
        while (true) {
            const i = next++;
            if (i >= targets.length) break;
            try {
                results[i] = await runOne(targets[i]);
            } catch (e) {
                results[i] = { pick: targets[i], xprepError: e.message };
            }
            done++;
            if (done % 5 === 0 || done === targets.length) {
                const el = ((Date.now() - t0) / 1000).toFixed(0);
                fs.writeSync(1, `\r  ${done}/${targets.length}  ${el}s elapsed`);
            }
        }
    };
    await Promise.all(Array.from({ length: Math.min(concurrency, targets.length) }, worker));
    fs.writeSync(1, '\n');
    return results;
}

function main() {
    const argv = process.argv.slice(2);
    let only = null, limit = null, sgFilter = null, concurrency = 4;
    for (let i = 0; i < argv.length; i++) {
        const a = argv[i];
        if (a === '--limit') { limit = parseInt(argv[++i], 10); continue; }
        if (a === '--sg') { sgFilter = parseInt(argv[++i], 10); continue; }
        if (a === '--id') { only = argv[++i]; continue; }
        if (a === '--concurrency') { concurrency = parseInt(argv[++i], 10); continue; }
        if (!a.startsWith('--')) only = a;
    }

    let targets = picks;
    if (only) targets = picks.filter(p => p.id === only);
    if (sgFilter) targets = picks.filter(p => p.sgNumber === sgFilter);
    if (limit) targets = targets.slice(0, limit);
    if (!targets.length) {
        console.error('No COD entries match the given filter.');
        process.exit(1);
    }

    console.log(`xrdspace vs XPREP — ${targets.length} COD entries`);
    console.log('='.repeat(72));

    runAll(targets, concurrency).then(results => {
        const summary = {
            total: results.length, assessed: 0, skip: 0,
            xrdspace: { pass: 0, near: 0, fail: 0 },
            xprep: { pass: 0, near: 0, fail: 0, error: 0 },
            agree: 0, bothPass: 0, xprepOnly: 0, xrdspaceOnly: 0,
            bothFail: 0, disagree: 0,
            bySystem: {},
        };
        const entries = [];

        for (const r of results) {
            const p = r.pick;
            const sys = systemOf(p.sgNumber);
            summary.bySystem[sys] = summary.bySystem[sys] ||
                { total: 0, xrdPass: 0, xpPass: 0, agree: 0 };

            if (r.skip) {
                summary.skip++;
                summary.bySystem[sys].total++;
                console.log(`  [SKIP] ${p.id}  ${r.skip}`);
                entries.push({ id: p.id, published: p.sgNumber, status: 'SKIP', skip: r.skip });
                continue;
            }
            summary.assessed++;
            const { xpPass, xrdPass, agree } = classify(r);
            if (xrdPass) summary.xrdspace.pass++;
            else if (r.xrdNear) summary.xrdspace.near++;
            else summary.xrdspace.fail++;
            if (r.xprepError) summary.xprep.error++;
            else if (xpPass) summary.xprep.pass++;
            else if (r.xprepNear) summary.xprep.near++;
            else summary.xprep.fail++;
            if (agree) summary.agree++;
            if (xpPass && xrdPass) summary.bothPass++;
            else if (xpPass && !xrdPass) summary.xprepOnly++;
            else if (!xpPass && xrdPass) summary.xrdspaceOnly++;
            else if (agree) summary.bothFail++;
            else summary.disagree++;

            const bs = summary.bySystem[sys];
            bs.total++;
            if (xrdPass) bs.xrdPass++;
            if (xpPass) bs.xpPass++;
            if (agree) bs.agree++;

            const mark = xpPass && xrdPass ? '=' : xpPass ? 'X' : xrdPass ? 'R' : (agree ? '~' : '!');
            console.log(`  [${mark}] ${p.id}  published ${p.sgNumber} (${p.sg})  ` +
                `xrdspace ${r.xrdNum ?? '-'}${r.xrdHm ? ' ' + r.xrdHm : ''}  ` +
                `xprep ${r.xprepNum ?? '-'}${r.xprepHm ? ' ' + r.xprepHm : ''}` +
                (r.xprepError ? `  (xprep: ${r.xprepError})` : ''));

            entries.push({
                id: p.id, published: p.sgNumber, publishedHm: p.sg,
                xrdspace: r.xrdNum, xrdspaceHm: r.xrdHm, xrdspaceNear: r.xrdNear,
                xprep: r.xprepNum, xprepHm: r.xprepHm, xprepNear: r.xprepNear,
                xprepError: r.xprepError, reflections: r.reflections,
                agree, xpPass, xrdPass,
            });
        }

        console.log('='.repeat(72));
        console.log('By crystal system:');
        for (const [s, v] of Object.entries(summary.bySystem)) {
            console.log(`  ${s.padEnd(12)} n=${String(v.total).padStart(3)}  ` +
                `xrdspace PASS ${String(v.xrdPass).padStart(3)}  ` +
                `xprep PASS ${String(v.xpPass).padStart(3)}  agree ${String(v.agree).padStart(3)}`);
        }
        console.log('='.repeat(72));
        const a = summary.assessed || 1;
        console.log(`Assessed: ${summary.assessed}  (skipped ${summary.skip})`);
        console.log(`xrdspace : PASS ${summary.xrdspace.pass}  NEAR ${summary.xrdspace.near}  FAIL ${summary.xrdspace.fail}  ` +
            `(${((summary.xrdspace.pass + summary.xrdspace.near) / a * 100).toFixed(1)}% recovered)`);
        console.log(`xprep    : PASS ${summary.xprep.pass}  NEAR ${summary.xprep.near}  FAIL ${summary.xprep.fail}  ` +
            `ERROR ${summary.xprep.error}  ` +
            `(${((summary.xprep.pass + summary.xprep.near) / a * 100).toFixed(1)}% recovered)`);
        console.log(`head-to-head: agree ${summary.agree}  both PASS ${summary.bothPass}  ` +
            `xprep only ${summary.xprepOnly}  xrdspace only ${summary.xrdspaceOnly}  ` +
            `both wrong (same) ${summary.bothFail}  disagree ${summary.disagree}`);
        console.log(`xrdspace exact-match rate ${(summary.xrdspace.pass / a * 100).toFixed(1)}%  |  ` +
            `xprep exact-match rate ${(summary.xprep.pass / a * 100).toFixed(1)}%`);

        const report = { ...summary, entries };
        fs.writeFileSync(REPORT_JSON, JSON.stringify(report, null, 2));
        console.log(`Report: ${REPORT_JSON}`);
    });
}

main();
