#!/usr/bin/env node
// Copyright (c) 2026 Denis Spasyuk. MIT License.
//
// Wide (10k) head-to-head validation of xrdspace against Bruker XPREP on a
// fresh, diverse set of COD single-crystal entries (tests/cod-picks-10k.json),
// disjoint from the original 2000-entry set.
//
// Two resumable passes, each appending to its own JSONL:
//   --xrd    run xrdspace (analyzeHkl) on every pick: space-group accuracy
//            (PASS/NEAR/FAIL vs the published COD SG) plus output data-quality
//            metrics (Rmerge, Rmeas, R(pim), completeness, d(I/sig=1),
//            d(CC1/2=0.30), mean I/sig, multiplicity, #unique).
//   --xprep  run XPREP (scripts/xprep-run.py) on a pick: the SG it chooses,
//            its R(sym), CFOM and candidate list.
//   --report merge both JSONL files with the picks and print a summary table
//            plus the failure lists (where each program is wrong).
//
// Usage:
//   node tests/xrdspace-compare-10k.js --xrd [--limit N] [--concurrency C]
//   node tests/xrdspace-compare-10k.js --xprep [--limit N] [--concurrency C]
//   node tests/xrdspace-compare-10k.js --report
//   node tests/xrdspace-compare-10k.js --xrd --only-ids a,b,c

import fs from 'node:fs';
import path from 'node:path';
import { spawn } from 'node:child_process';
import { fileURLToPath } from 'node:url';

import { parseHkl } from '../src/hkl-parser.js';
import { analyzeHkl } from '../src/index.js';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.join(__dirname, '..');
const CACHE_DIR = path.join(ROOT, 'HKLs', 'cod');
const WORK_DIR = path.join(__dirname, 'xprep-work-10k');
const DRIVER = path.join(ROOT, 'scripts', 'xprep-run.py');
const PICKS = path.join(__dirname, 'cod-picks-10k.json');
const XRD_JSONL = path.join(__dirname, 'compare-10k-xrd.jsonl');
const XP_JSONL = path.join(__dirname, 'compare-10k-xprep.jsonl');

function systemOf(n) {
    if (n <= 2) return 'triclinic';
    if (n <= 15) return 'monoclinic';
    if (n <= 74) return 'orthorhombic';
    if (n <= 142) return 'tetragonal';
    if (n <= 167) return 'trigonal';
    if (n <= 194) return 'hexagonal';
    return 'cubic';
}

function cellFromPick(p) {
    const [a, b, c, alpha, beta, gamma] = p.cell;
    return { a, b, c, alpha, beta, gamma };
}

// --- resumable JSONL helpers ---

function loadJsonlDone(file) {
    const done = new Set();
    if (!fs.existsSync(file)) return done;
    for (const line of fs.readFileSync(file, 'utf8').split('\n')) {
        if (!line.trim()) continue;
        try { done.add(JSON.parse(line).id); } catch { /* ignore */ }
    }
    return done;
}

function appendJsonl(file, obj) {
    fs.appendFileSync(file, JSON.stringify(obj) + '\n');
}

// --- xrdspace pass ---

async function runXrdOne(pick) {
    const cache = path.join(CACHE_DIR, `${pick.id}.hkl`);
    if (!fs.existsSync(cache)) return { id: pick.id, status: 'MISSING' };
    const text = fs.readFileSync(cache, 'utf8');
    let reflections;
    try {
        reflections = parseHkl(text).reflections;
    } catch (e) {
        return { id: pick.id, status: 'PARSE', error: e.message };
    }
    if (!reflections.length) return { id: pick.id, status: 'EMPTY' };

    let r;
    try {
        r = await analyzeHkl(text, { cell: cellFromPick(pick), quality: true });
    } catch (e) {
        return { id: pick.id, status: 'ERROR', error: e.message };
    }
    if (!r.ok) return { id: pick.id, status: 'ERROR', error: r.error };

    const num = r.best ? r.best.id : null;
    const near = r.candidates.some(c => c.id === pick.sgNumber && c.violations === 0);
    const st = r.merge.statistics;
    const rec = {
        id: pick.id,
        status: num === pick.sgNumber ? 'PASS' : (near ? 'NEAR' : 'FAIL'),
        nRefl: r.summary.nReflections,
        xrd: {
            num,
            hm: r.best ? r.best.hm : null,
            near,
            violations: r.best ? (r.candidates.find(c => c.id === num) || {}).violations : null,
            laue: r.summary.laueClass,
            centering: r.summary.centering,
            system: r.summary.crystalSystem,
            centric: r.summary.centricity,
            chiral: r.summary.chiral,
            rMerge: st.rMerge,
            rMeas: st.rMeas,
            rPim: st.rPim,
            completeness: st.completeness,
            dIsig1: st.dIsig1,
            dCC30: st.dCC30,
            meanIsig: st.meanIsig,
            nUnique: st.nUnique,
            mult: st.meanMultiplicity,
        },
    };
    return rec;
}

async function passXrd(targets, concurrency) {
    const done = loadJsonlDone(XRD_JSONL);
    const todo = targets.filter(p => !done.has(p.id));
    console.log(`xrdspace pass: ${targets.length} picks, ${todo.length} to do, ${concurrency} workers`);
    const t0 = Date.now();
    let next = 0, cnt = 0;
    const worker = async () => {
        while (true) {
            const i = next++;
            if (i >= todo.length) break;
            const p = todo[i];
            let rec;
            try { rec = await runXrdOne(p); }
            catch (e) { rec = { id: p.id, status: 'EXC', error: e.message }; }
            appendJsonl(XRD_JSONL, rec);
            cnt++;
            if (cnt % 50 === 0 || cnt === todo.length) {
                const el = ((Date.now() - t0) / 1000).toFixed(0);
                fs.writeSync(1, `\r  xrd ${cnt}/${todo.length}  ${el}s`);
            }
        }
    };
    await Promise.all(Array.from({ length: Math.min(concurrency, todo.length) }, worker));
    fs.writeSync(1, '\n');
    console.log(`xrdspace pass done -> ${XRD_JSONL}`);
}

// --- XPREP pass ---

function fmt(x, d = 4) { return Number(x).toFixed(d); }

function buildFcf(id, cell, reflections) {
    const lines = [];
    lines.push('#');
    lines.push(`# COD ${id} — xrdspace 10k XPREP comparison`);
    lines.push('#');
    lines.push(`data_${id}`);
    lines.push("_shelx_title ' xrdspace-10k-xprep'");
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
        const child = spawn('python3', [DRIVER, fcf], { stdio: ['ignore', 'pipe', 'pipe'] });
        let out = '';
        const timer = setTimeout(() => {
            try { child.kill('SIGKILL'); } catch { /* ignore */ }
            resolve({ ok: false, error: 'driver timeout' });
        }, 60000);
        child.stdout.on('data', d => { out += d; });
        child.on('error', e => { clearTimeout(timer); resolve({ ok: false, error: String(e.message || e) }); });
        child.on('close', () => {
            clearTimeout(timer);
            const lines = out.trim().split(/\r?\n/).filter(Boolean);
            if (!lines.length) return resolve({ ok: false, error: 'no driver output' });
            try { resolve(JSON.parse(lines[lines.length - 1])); }
            catch { resolve({ ok: false, error: `bad driver JSON` }); }
        });
    });
}

async function runXprepOne(pick) {
    const cache = path.join(CACHE_DIR, `${pick.id}.hkl`);
    if (!fs.existsSync(cache)) return { id: pick.id, xprep: { error: 'not cached' } };
    const text = fs.readFileSync(cache, 'utf8');
    let reflections;
    try { reflections = parseHkl(text).reflections; }
    catch (e) { return { id: pick.id, xprep: { error: e.message } }; }
    if (!reflections.length) return { id: pick.id, xprep: { error: 'no reflections' } };

    fs.mkdirSync(WORK_DIR, { recursive: true });
    const fcf = path.join(WORK_DIR, `${pick.id}.fcf`);
    fs.writeFileSync(fcf, buildFcf(pick.id, cellFromPick(pick), reflections));
    const x = await runXprep(pick.id, fcf);
    return {
        id: pick.id,
        xprep: {
            ok: x.ok,
            num: x.sgNumber ?? null,
            hm: x.hm ?? null,
            rsym: (x.candidates || []).find(c => c.opt === x.chosen)?.rsym ?? null,
            cfom: (x.candidates || []).find(c => c.opt === x.chosen)?.cfom ?? null,
            nCand: (x.candidates || []).length,
            error: x.ok ? null : x.error,
        },
    };
}

async function passXprep(targets, concurrency) {
    const done = loadJsonlDone(XP_JSONL);
    const todo = targets.filter(p => !done.has(p.id));
    console.log(`XPREP pass: ${targets.length} picks, ${todo.length} to do, ${concurrency} workers`);
    const t0 = Date.now();
    let next = 0, cnt = 0;
    const worker = async () => {
        while (true) {
            const i = next++;
            if (i >= todo.length) break;
            const p = todo[i];
            let rec;
            try { rec = await runXprepOne(p); }
            catch (e) { rec = { id: p.id, xprep: { error: e.message } }; }
            appendJsonl(XP_JSONL, rec);
            cnt++;
            if (cnt % 10 === 0 || cnt === todo.length) {
                const el = ((Date.now() - t0) / 1000).toFixed(0);
                fs.writeSync(1, `\r  xprep ${cnt}/${todo.length}  ${el}s`);
            }
        }
    };
    await Promise.all(Array.from({ length: Math.min(concurrency, todo.length) }, worker));
    fs.writeSync(1, '\n');
    console.log(`XPREP pass done -> ${XP_JSONL}`);
}

// --- report ---

// Build a PASS/NEAR/FAIL bar chart (same style as xrdspace-cod.js) for the
// 10 000-entry wide-set validation.
function build10kSvg(bySys, totals) {
    const W = 900, H = 420, PL = 90, PR = 30, PT = 40, PB = 70;
    const colors = { PASS: '#2e9e4f', NEAR: '#e0a800', FAIL: '#d64545' };
    const systems = ['triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'trigonal', 'hexagonal', 'cubic'];
    const maxCount = Math.max(...systems.map(s => (bySys[s] ? bySys[s].total : 0)), 1);
    const innerW = W - PL - PR, innerH = H - PT - PB;
    const groupW = innerW / systems.length;
    const barW = groupW * 0.22;

    const parts = [];
    parts.push(`<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="sans-serif">`);
    parts.push(`<rect x="0" y="0" width="${W}" height="${H}" fill="#ffffff"/>`);
    parts.push(`<text x="${W/2}" y="24" text-anchor="middle" font-size="18" font-weight="bold" fill="#222">xrdspace space-group determination vs COD (10 000 structures, wide set)</text>`);

    for (let g = 0; g <= 4; g++) {
        const v = maxCount * g / 4;
        const y = PT + innerH - innerH * g / 4;
        parts.push(`<line x1="${PL}" y1="${y}" x2="${W-PR}" y2="${y}" stroke="#e5e5e5" stroke-width="1"/>`);
        parts.push(`<text x="${PL-8}" y="${y+4}" text-anchor="end" font-size="11" fill="#666">${Math.round(v)}</text>`);
    }

    systems.forEach((sys, si) => {
        const b = bySys[sys] || { total: 0, pass: 0, near: 0, fail: 0 };
        const gx = PL + si * groupW + (groupW - 3 * barW) / 2;
        [['PASS', b.pass], ['NEAR', b.near], ['FAIL', b.fail]].forEach(([label, val], ci) => {
            const h = val / maxCount * innerH;
            const x = gx + ci * (barW + 2);
            const y = PT + innerH - h;
            parts.push(`<rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${barW.toFixed(1)}" height="${h.toFixed(1)}" fill="${colors[label]}" rx="2"/>`);
            if (val > 0) {
                parts.push(`<text x="${(x + barW/2).toFixed(1)}" y="${(y - 4).toFixed(1)}" text-anchor="middle" font-size="11" fill="#333" font-weight="bold">${val}</text>`);
            }
        });
        const label = sys.charAt(0).toUpperCase() + sys.slice(1);
        parts.push(`<text x="${(gx + groupW/2).toFixed(1)}" y="${H-PB+16}" text-anchor="middle" font-size="12" fill="#333">${label}</text>`);
        parts.push(`<text x="${(gx + groupW/2).toFixed(1)}" y="${H-PB+31}" text-anchor="middle" font-size="11" fill="#888">n=${b.total}</text>`);
    });

    const lx = W/2 - 90, ly = H - 22;
    ['PASS', 'NEAR', 'FAIL'].forEach((label, i) => {
        const x = lx + i * 130;
        parts.push(`<rect x="${x}" y="${ly-10}" width="14" height="14" fill="${colors[label]}" rx="2"/>`);
        parts.push(`<text x="${x+20}" y="${ly+1}" font-size="12" fill="#333">${label}</text>`);
    });

    parts.push(`</svg>`);
    return parts.join('\n');
}

function report(picksFile = PICKS) {
    const picks = JSON.parse(fs.readFileSync(picksFile, 'utf8'));
    const byId = new Map(picks.map(p => [p.id, p]));
    const xrd = new Map();
    if (fs.existsSync(XRD_JSONL)) {
        for (const line of fs.readFileSync(XRD_JSONL, 'utf8').split('\n')) {
            if (!line.trim()) continue;
            const r = JSON.parse(line); xrd.set(r.id, r);
        }
    }
    const xp = new Map();
    if (fs.existsSync(XP_JSONL)) {
        for (const line of fs.readFileSync(XP_JSONL, 'utf8').split('\n')) {
            if (!line.trim()) continue;
            const r = JSON.parse(line); xp.set(r.id, r);
        }
    }
    console.log(`picks=${picks.length}  xrd records=${xrd.size}  xprep records=${xp.size}`);

    // xrdspace accuracy vs published truth
    const sys = {};
    let pass = 0, near = 0, fail = 0, skip = 0;
    for (const p of picks) {
        const r = xrd.get(p.id);
        const s = systemOf(p.sgNumber);
        sys[s] = sys[s] || { total: 0, pass: 0, near: 0, fail: 0, skip: 0 };
        if (!r) { sys[s].skip++; skip++; continue; }
        if (['MISSING', 'PARSE', 'EMPTY', 'ERROR', 'EXC'].includes(r.status)) { sys[s].skip++; skip++; continue; }
        sys[s].total++;
        if (r.status === 'PASS') { sys[s].pass++; pass++; }
        else if (r.status === 'NEAR') { sys[s].near++; near++; }
        else { sys[s].fail++; fail++; }
    }
    const assessed = pass + near + fail;
    console.log('\n=== xrdspace vs published COD space group ===');
    console.log('system         total   PASS   NEAR   FAIL   SKIP');
    for (const s of ['triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'trigonal', 'hexagonal', 'cubic']) {
        const v = sys[s] || { total: 0, pass: 0, near: 0, fail: 0, skip: 0 };
        console.log(`${s.padEnd(13)} ${String(v.total).padStart(5)} ${String(v.pass).padStart(6)} ${String(v.near).padStart(6)} ${String(v.fail).padStart(6)} ${String(v.skip).padStart(6)}`);
    }
    console.log(`${'TOTAL'.padEnd(13)} ${String(assessed).padStart(5)} ${String(pass).padStart(6)} ${String(near).padStart(6)} ${String(fail).padStart(6)} ${String(skip).padStart(6)}`);
    console.log(`xrdspace recovered (PASS+NEAR): ${((pass + near) / (assessed || 1) * 100).toFixed(1)}%   exact PASS: ${(pass / (assessed || 1) * 100).toFixed(1)}%`);

    // head-to-head where both have a record
    const both = picks.filter(p => xrd.get(p.id) && xp.get(p.id) && xp.get(p.id).xprep.ok);
    let agree = 0, xrdOnly = 0, xpOnly = 0, bothWrong = 0, disagree = 0;
    console.log(`\n=== head-to-head (both answered: ${both.length}) ===`);
    for (const p of both) {
        const xr = xrd.get(p.id).xrd.num;
        const xpnum = xp.get(p.id).xprep.num;
        const pub = p.sgNumber;
        if (xr === xpnum) agree++;
        if (xr === pub && xpnum !== pub) xrdOnly++;
        else if (xpnum === pub && xr !== pub) xpOnly++;
        else if (xr !== pub && xpnum !== pub && xr === xpnum) bothWrong++;
        else if (xr !== xpnum) disagree++;
    }
    console.log(`agree ${agree}  xrdspace-only-correct ${xrdOnly}  xprep-only-correct ${xpOnly}  both-wrong-same ${bothWrong}  disagree ${disagree}`);

    // data-quality: mean Rmerge / completeness for PASS entries
    const q = [...xrd.values()].filter(r => r.status === 'PASS' && r.xrd);
    if (q.length) {
        const mean = a => a.reduce((s, v) => s + v, 0) / a.length;
        console.log(`\n=== xrdspace output quality (PASS entries, n=${q.length}) ===`);
        console.log(`  mean Rmerge ${(mean(q.map(r => r.xrd.rMerge)) * 100).toFixed(2)}%   ` +
            `mean completeness ${(mean(q.map(r => r.xrd.completeness)) * 100).toFixed(1)}%   ` +
            `mean <I/sig> ${mean(q.map(r => r.xrd.meanIsig)).toFixed(1)}`);
    }

    // failure lists
    const fails = picks.filter(p => { const r = xrd.get(p.id); return r && r.status === 'FAIL'; });
    console.log(`\n=== xrdspace FAIL list (${fails.length}) — top 40 ===`);
    for (const p of fails.slice(0, 40)) {
        const r = xrd.get(p.id);
        console.log(`  ${p.id}  pub ${p.sgNumber} (${p.sg})  got ${r.xrd.num} (${r.xrd.hm})  [${systemOf(p.sgNumber)}]`);
    }

    // save a merged report
    const merged = picks.map(p => {
        const r = xrd.get(p.id); const x = xp.get(p.id);
        return {
            id: p.id, published: p.sgNumber, publishedHm: p.sg, system: systemOf(p.sgNumber),
            xrd: r ? r.xrd : null, xrdStatus: r ? r.status : null,
            xprep: x ? x.xprep : null,
        };
    });
    const out = path.join(__dirname, 'compare-10k-report.json');
    fs.writeFileSync(out, JSON.stringify({
        totals: { assessed, pass, near, fail, skip, bySystem: sys },
        headToHead: { both: both.length, agree, xrdOnly, xpOnly, bothWrong, disagree },
        entries: merged,
    }, null, 1));
    console.log(`\nreport -> ${out}`);

    // PASS/NEAR/FAIL chart for the wide set.
    const chartPath = path.join(__dirname, 'xrdspace-report-10k.svg');
    fs.writeFileSync(chartPath, build10kSvg(sys, { total: assessed }));
    console.log(`chart  -> ${chartPath}`);
}

// --- main ---

function main() {
    const argv = process.argv.slice(2);
    let mode = null, limit = null, concurrency = 8, onlyIds = null, picksFile = PICKS;
    for (let i = 0; i < argv.length; i++) {
        const a = argv[i];
        if (a === '--xrd' || a === '--xprep' || a === '--report') { mode = a.slice(2); continue; }
        if (a === '--limit') { limit = parseInt(argv[++i], 10); continue; }
        if (a === '--concurrency') { concurrency = parseInt(argv[++i], 10); continue; }
        if (a === '--only-ids') { onlyIds = new Set(argv[++i].split(',')); continue; }
        if (a === '--picks') { picksFile = argv[++i]; continue; }
    }
    if (!mode) {
        console.error('usage: --xrd | --xprep | --report  [--limit N] [--concurrency C] [--only-ids a,b,c] [--picks file]');
        process.exit(1);
    }
    if (mode === 'report') { report(picksFile); return; }

    let picks = JSON.parse(fs.readFileSync(picksFile, 'utf8'));
    if (onlyIds) picks = picks.filter(p => onlyIds.has(p.id));
    if (limit) picks = picks.slice(0, limit);
    if (mode === 'xrd') passXrd(picks, concurrency);
    else passXprep(picks, concurrency);
}

main();
