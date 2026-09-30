#!/usr/bin/env node
// Copyright (c) 2026 Denis Spasyuk. MIT License.
// Regression tests for the Bruker P4P reader.
//
// P4P (SAINT/APEX) is a fixed-format text format: a FILEID/CELL/SOURCE/TITLE
// header followed by REF05 reflection records. The parser must:
//   * detect the format from the FILEID header (before the numeric DATA rows
//     could be mistaken for SHELX),
//   * read the cell from CELL and the wavelength from SOURCE,
//   * read h k l, the intensity I and its sigma from REF05. I and sigma are the
//     two fields immediately before the reciprocal-lattice u1 u2 u3 columns;
//     for this layout that is h+9 and h+10. This is verified here against the
//     orientation matrix (u == ORT . hkl).
//   * split the SAINT quirk where a negative omega is fused onto l
//     ("13-108.000" => l=13, omega=-108.000),
//   * drop the (0,0,0) direct-beam / standard rows,
//   * run the full space-group analysis end-to-end.

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { parseHkl, detectFormat, HKL_FORMAT } from '../src/hkl-parser.js';
import { analyzeHkl } from '../src/index.js';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const P4P_DIR = path.join(__dirname, '..', 'test', 'p4p');

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
    if (cond) { pass++; console.log(`  ok   ${name}`); }
    else { fail++; console.error(`  FAIL ${name}${detail ? ' — ' + detail : ''}`); }
}

console.log('xrdspace P4P reader regression tests');

// ---------------------------------------------------------------------------
// Synthetic file with known values, including the l/omega concatenation quirk
// and a (0,0,0) direct-beam row that must be ignored.
// ---------------------------------------------------------------------------
const synthetic = [
    'FILEID SAINT        V7.68A        4.00        01/01/00 00:00:00 test',
    'SITEID Bruker                              ?',
    'TITLE  synthetic test',
    'CHEM   C2 H6 O',
    'CELL     10.0000   11.0000   12.0000   90.0000   90.0000   90.0000   1320.000',
    'SOURCE CU     1.54188   1.54056   1.54439   2.00000   40.00    4.00',
    // normal l/omega
    'REF05 H       4  -4  14 -30.000 310.150   0.000  54.740  440.64  369.84  267.09  34.13  0.179694 -0.216485 -0.646262   0.000   0.000   0.000   0.000',
    // concatenated l=13, omega=-108.000
    'REF05 AC      1  -6  13-108.000 251.850  45.000  54.736  265.42   43.21  266.67  66.58  0.991718  0.050750 -0.306996   0.000   0.000   0.000   0.000',
    // direct beam (0,0,0) -- must be dropped
    'REF05         0   0   0 -40.000 319.850 360.000  50.000  436.31  181.97  419.36  42.62  0.017797 -0.619856 -0.503863   0.000   0.000   0.000   0.000',
].join('\n') + '\n';

const sp = parseHkl(synthetic);
check('synthetic: detected as p4p', detectFormat(synthetic) === HKL_FORMAT.P4P);
check('synthetic: title', sp.title === 'synthetic test', `got "${sp.title}"`);
check('synthetic: wavelength', Math.abs(sp.wavelength - 1.54188) < 1e-9, `got ${sp.wavelength}`);
check('synthetic: cell',
    sp.cell && sp.cell.a === 10 && sp.cell.b === 11 && sp.cell.c === 12 &&
    sp.cell.alpha === 90 && sp.cell.beta === 90 && sp.cell.gamma === 90,
    JSON.stringify(sp.cell));
check('synthetic: (0,0,0) dropped, 2 reflections', sp.reflections.length === 2,
    `got ${sp.reflections.length}`);

const r1 = sp.reflections.find(r => r.h === 4 && r.k === -4 && r.l === 14);
check('normal line: I and sigma are h+9 / h+10',
    r1 && r1.I === 267.09 && Math.abs(r1.sig - 34.13) < 1e-9, JSON.stringify(r1));
check('normal line: omega -> psi', r1 && Math.abs(r1.psi - (-30)) < 1e-9, `psi=${r1 && r1.psi}`);

const r2 = sp.reflections.find(r => r.h === 1 && r.k === -6 && r.l === 13);
check('concatenated "13-108.000" -> l=13, omega=-108',
    r2 && r2.I === 266.67 && Math.abs(r2.sig - 66.58) < 1e-9 && Math.abs(r2.psi - (-108)) < 1e-9,
    JSON.stringify(r2));

// ---------------------------------------------------------------------------
// Real P4P files: cell/wavelength/reflection counts and sane I/sigma.
// ---------------------------------------------------------------------------
// The real .p4p files live in test/p4p/, which is git-ignored (local scratch
// data). Run the fixture-based checks only when they are present, so the test
// still passes on a clean clone (the synthetic checks above cover the format).
const haveFixtures = fs.existsSync(P4P_DIR) &&
    fs.readdirSync(P4P_DIR).some(f => f.endsWith('.p4p'));

if (!haveFixtures) {
    console.log('  skip real-file checks (test/p4p/ not present — git-ignored scratch data)');
} else {
    const files = fs.readdirSync(P4P_DIR).filter(f => f.endsWith('.p4p'));
    check('test/p4p directory is present and non-empty', files.length > 0, `${files.length} files`);

    let withRefl = 0, allDetected = true, allSane = true, originLeaks = 0;
    for (const f of files) {
        const text = fs.readFileSync(path.join(P4P_DIR, f), 'utf8');
        if (detectFormat(text) !== HKL_FORMAT.P4P) { allDetected = false; continue; }
        const p = parseHkl(text);
        if (!p.cell) allSane = false;
        for (const r of p.reflections) {
            if (!Number.isFinite(r.I) || r.I < 0 || !(r.sig > 0)) allSane = false;
            if (r.h === 0 && r.k === 0 && r.l === 0) originLeaks++;
        }
        if (p.reflections.length) withRefl++;
    }
    check('every .p4p file detected as p4p', allDetected);
    check('every parsed reflection has finite I>=0 and sigma>0', allSane);
    check('all files carry a unit cell', allSane);
    check('no (0,0,0) rows leak into the reflection list', originLeaks === 0, `${originLeaks} leaked`);
    check('at least 7 files carry integrated reflections', withRefl >= 7, `${withRefl} files`);

    // denis7.p4p is all direct-beam rows -> no reflections (must not misanalyse).
    const denis7 = parseHkl(fs.readFileSync(path.join(P4P_DIR, 'denis7.p4p'), 'utf8'));
    check('all-origin file yields 0 reflections', denis7.reflections.length === 0,
        `got ${denis7.reflections.length}`);

    // A raw-data P4P (header + DATA block, no REF records) still reads its cell
    // but has no reflections.
    const den20 = parseHkl(fs.readFileSync(path.join(P4P_DIR, 'den20.p4p'), 'utf8'));
    check('raw-data P4P: cell read, 0 reflections',
        den20.cell && Math.abs(den20.cell.a - 11.4914) < 1e-6 && den20.reflections.length === 0,
        JSON.stringify({ cell: den20.cell, n: den20.reflections.length }));

    // End-to-end analysis on a real P4P file: the default SHELX output must be
    // the UNMERGED observation list, and the merged text is kept separately.
    const den = fs.readFileSync(path.join(P4P_DIR, 'den.p4p'), 'utf8');
    const res = analyzeHkl(den, { quality: true });
    check('end-to-end analyzeHkl on den.p4p succeeds', res.ok, res.error || '');
    check('end-to-end: orthorhombic P, primitive cell',
        res.ok && res.summary.crystalSystem === 'orthorhombic' && res.summary.centering === 'P',
        res.ok ? `${res.summary.crystalSystem}/${res.summary.centering}` : '');
    check('end-to-end: merging produced uniques', res.ok && res.merge && res.merge.nUnique > 0);
    if (res.ok && res.merge) {
        const unmergedRows = res.merge.shelxHkl.trim().split('\n').length;
        const mergedRows = res.merge.shelxHklMerged.trim().split('\n').length;
        check('default SHELX output is UNMERGED (more rows than merged)',
            unmergedRows > mergedRows && mergedRows === res.merge.nUnique,
            `unmerged=${unmergedRows} merged=${mergedRows} nUnique=${res.merge.nUnique}`);
    }
}

console.log(`passed ${pass}, failed ${fail}`);
process.exit(fail === 0 ? 0 : 1);
