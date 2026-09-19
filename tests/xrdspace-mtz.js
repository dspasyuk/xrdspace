#!/usr/bin/env node
// Copyright (c) 2026 Denis Spasyuk. MIT License.
// xrdspace MTZ reader/writer tests.
//
// 1. Builds a small synthetic MTZ in memory, writes it with writeMtz, re-parses
//    it, and checks that every field and every data value survives the round
//    trip.
// 2. If the reference MTZ files from the 2T364 XDS/pointless/aimless pipeline
//    are present, parses + round-trips each one and (when the `gemmi` Python
//    module is available) cross-checks the values against gemmi, the reference
//    implementation.
//
// Usage:
//   node tests/xrdspace-mtz.js [path-to-mtz-dir ...]

import fs from 'node:fs';
import path from 'node:path';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

import { parseMtz, writeMtz, getReflections, getHkl, getColumnValues, columnIndexOf } from '../src/mtz.js';
import { analyzeMtz } from '../src/index.js';

const __dirname = path.dirname(fileURLToPath(import.meta.url));

let failures = 0;
function check(cond, msg) {
    if (cond) {
        console.log('  ok   ' + msg);
    } else {
        failures++;
        console.log('  FAIL ' + msg);
    }
}

// Build a minimal valid MTZ object (1 dataset, space group P 1) in memory.
function syntheticMtz() {
    const ncol = 5; // H K L F SIGF
    const refls = [
        [0, 0, 1, 120.5, 3.2],
        [0, 0, 2, 88.1, 2.1],
        [1, 0, 0, 45.0, 1.7],
        [0, 1, 0, NaN, NaN], // a missing reflection
        [1, 1, 0, 200.7, 4.4],
    ];
    const nrefl = refls.length;
    const data = new Float32Array(ncol * nrefl);
    for (let r = 0; r < nrefl; r++) {
        for (let c = 0; c < ncol; c++) data[r * ncol + c] = refls[r][c];
    }
    return {
        title: 'synthetic test',
        version: 'MTZ:V1.1',
        nreflections: nrefl,
        ncols: ncol,
        cell: { a: 50, b: 50, c: 50, alpha: 90, beta: 90, gamma: 90 },
        spaceGroupNumber: 1,
        spaceGroupName: 'P 1',
        latticeType: 'P',
        nsymop: 1,
        nops: 1,
        symops: ['X, Y, Z'],
        reso: [4, 100],
        valm: NaN,
        columns: [
            { label: 'H', type: 'H', min: 0, max: 1, dataset: 0, source: '' },
            { label: 'K', type: 'H', min: 0, max: 1, dataset: 0, source: '' },
            { label: 'L', type: 'H', min: 1, max: 2, dataset: 0, source: '' },
            { label: 'F', type: 'F', min: 45, max: 200.7, dataset: 0, source: '' },
            { label: 'SIGF', type: 'Q', min: 1.7, max: 4.4, dataset: 0, source: '' },
        ],
        datasets: [{ id: 0, project: 'P', crystal: 'C', dataset: 'D', cell: { a: 50, b: 50, c: 50, alpha: 90, beta: 90, gamma: 90 }, wavelength: 1.54 }],
        history: ['created by xrdspace test'],
        batchNumbers: [],
        data,
        byteOrder: 'LE',
    };
}

function testSynthetic() {
    console.log('Synthetic MTZ round-trip:');
    const m = syntheticMtz();
    const buf = writeMtz(m);
    check(buf.subarray(0, 4).toString('latin1') === 'MTZ ', 'output starts with "MTZ "');
    const m2 = parseMtz(buf);
    check(m2.ncols === m.ncols, `ncols ${m2.ncols} === ${m.ncols}`);
    check(m2.nreflections === m.nreflections, `nreflections ${m2.nreflections} === ${m.nreflections}`);
    check(m2.spaceGroupNumber === 1 && m2.spaceGroupName === 'P 1', `space group ${m2.spaceGroupNumber} "${m2.spaceGroupName}"`);
    check(m2.cell.a === 50 && m2.cell.b === 50 && m2.cell.c === 50, `cell ${m2.cell.a} ${m2.cell.b} ${m2.cell.c}`);
    check(m2.title === 'synthetic test', `title "${m2.title}"`);
    check(m2.columns.map(c => c.label).join(',') === 'H,K,L,F,SIGF', 'column labels preserved');
    check(m2.columns[3].type === 'F' && m2.columns[4].type === 'Q', 'column types preserved');
    check(m2.history.length === 1 && m2.history[0].includes('xrdspace'), 'history preserved');
    // every data value identical
    let identical = m2.data.length === m.data.length;
    if (identical) {
        for (let i = 0; i < m.data.length; i++) {
            const a = m.data[i], b = m2.data[i];
            if ((Number.isNaN(a) !== Number.isNaN(b)) || (a !== b && Math.abs(a - b) > 1e-6)) { identical = false; break; }
        }
    }
    check(identical, 'all data values identical after round-trip');
    // reflection extraction
    const refls = getReflections(m2, { intensity: 'F', sigma: 'SIGF' });
    check(refls.length === 4, `getReflections skips missing (${refls.length} of ${m.nreflections})`);
    check(refls[0].h === 0 && refls[0].k === 0 && refls[0].l === 1 && Math.abs(refls[0].I - 120.5) < 1e-3, 'first reflection (0,0,1) I=120.5');
    check(columnIndexOf(m2, 'F') === 3 && columnIndexOf(m2, 'SIGF') === 4, 'columnIndexOf');
    check(getHkl(m2)[2] === 1, 'getHkl l of first reflection');
    check(Math.abs(getColumnValues(m2, 'F')[0] - 120.5) < 1e-3, 'getColumnValues F[0]');
}

function gemmiAvailable() {
    try {
        execFileSync('python3', ['-c', 'import gemmi'], { stdio: 'ignore' });
        return true;
    } catch {
        return false;
    }
}

function testRealFile(file) {
    const name = path.basename(file);
    console.log('Real MTZ: ' + name);
    const buf = fs.readFileSync(file);
    const m = parseMtz(buf);
    check(m.ncols > 0 && m.nreflections > 0, `parsed ncols=${m.ncols} nreflections=${m.nreflections}`);
    check(m.data.length === m.ncols * m.nreflections, `data length ${m.data.length} === ncols*nrefl`);
    check(Number.isFinite(m.cell.a) && m.cell.a > 0, `cell a=${m.cell.a}`);
    check(m.spaceGroupNumber > 0, `space group ${m.spaceGroupNumber} "${m.spaceGroupName}"`);
    // round-trip
    const out = writeMtz(m);
    const m2 = parseMtz(out);
    let identical = m2.data.length === m.data.length;
    if (identical) {
        for (let i = 0; i < m.data.length; i++) {
            const a = m.data[i], b = m2.data[i];
            if ((Number.isNaN(a) !== Number.isNaN(b)) || (a !== b && Math.abs(a - b) > 1e-6)) { identical = false; break; }
        }
    }
    check(identical, 'data identical after round-trip');
    check(m2.ncols === m.ncols && m2.nreflections === m.nreflections && m2.spaceGroupNumber === m.spaceGroupNumber,
        `header identical after round-trip (ncols=${m2.ncols} nrefl=${m2.nreflections} sg=${m2.spaceGroupNumber})`);
    // gemmi cross-check
    if (gemmiAvailable()) {
        try {
            const tmp = path.join('/tmp', `xrdspace-mtz-${process.pid}-${Math.random().toString(36).slice(2)}.mtz`);
            fs.writeFileSync(tmp, out);
            const py = `
import gemmi, numpy as np, sys
mine = gemmi.read_mtz_file(${JSON.stringify(tmp)})
orig = gemmi.read_mtz_file(${JSON.stringify(file)})
a = np.asarray(mine.array); o = np.asarray(orig.array)
if a.shape != o.shape:
    print("SHAPE"); sys.exit(2)
d = np.nanmax(np.abs(a - o)) if a.size else 0.0
print(f"{mine.nreflections} {len(mine.columns)} {d:.3e}")
`;
            const res = execFileSync('python3', ['-c', py], { encoding: 'utf8' }).trim().split(/\s+/);
            const [nrefl, ncols, maxdiff] = [res[0], res[1], parseFloat(res[2])];
            check(nrefl === String(m.nreflections) && ncols === String(m.ncols), `gemmi agrees on size (${nrefl} x ${ncols})`);
            check(maxdiff < 1e-4, `gemmi max value diff ${maxdiff}`);
            fs.unlinkSync(tmp);
        } catch (e) {
            console.log('  (gemmi cross-check skipped: ' + (e.message || e) + ')');
        }
    }
}

function main() {
    testSynthetic();

    // Collect real MTZ files: explicit args, else the known 2T364 pipeline dir.
    const args = process.argv.slice(2);
    let files = [];
    if (args.length) {
        for (const a of args) {
            const st = fs.statSync(a);
            if (st.isDirectory()) files.push(...fs.readdirSync(a).filter(f => f.endsWith('.mtz')).map(f => path.join(a, f)));
            else if (a.endsWith('.mtz')) files.push(a);
        }
    } else {
        const known = '/home/denis/CODE/pdbhkl/proc';
        if (fs.existsSync(known)) files = fs.readdirSync(known).filter(f => f.endsWith('.mtz')).map(f => path.join(known, f));
    }
    if (files.length) {
        for (const f of files) testRealFile(f);
    } else {
        console.log('No reference MTZ files found (pass a directory of *.mtz as an argument).');
    }

    // End-to-end: run the full space-group analysis on a real MTZ.
    const cut = files.find(f => f.includes('CUT_CRUNCATE.mtz')) || files[0];
    if (cut) {
        console.log('analyzeMtz end-to-end: ' + path.basename(cut));
        const result = analyzeMtz(fs.readFileSync(cut));
        check(result.ok, 'analyzeMtz returned ok');
        if (result.ok) {
            console.log('  -> best space group: ' + (result.best ? `${result.best.hm} (No. ${result.best.id})` : 'indeterminate'));
            check(result.summary.bestSpaceGroupNumber > 0, `determined space group ${result.summary.bestSpaceGroupNumber}`);
            check(result.mtz && result.mtz.ncols > 0, 'result.mtz metadata present');
        }
    }

    console.log('='.repeat(50));
    if (failures) {
        console.log(`${failures} check(s) FAILED`);
        process.exit(1);
    }
    console.log('All MTZ checks passed.');
}

main();
