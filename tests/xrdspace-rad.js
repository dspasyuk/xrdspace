#!/usr/bin/env node
// Copyright (c) 2026 Denis Spasyuk. MIT License.
// Tests for the RADDOSE-style beam-damage analysis (src/raddose.js) and the
// per-observation rotation-angle (PSI) parsing in the XDS_ASCII reader.
//
// Two parts:
//   1. A synthetic XDS_ASCII scan whose intensities decay by a KNOWN constant
//      k_true = I0 * exp(-k_true * dose). The analysis must recover k (and the
//      overall <I>_late / <I>_early ratio) close to the ground truth, and must
//      report the damage as present.
//   2. The real taurine XDS_ASCII.HKL scan: the parser must capture PSI and the
//      scan geometry, and the analysis must be usable with sensible values.
//   3. Graceful degradation: input with no PSI (SHELX five-column) reports
//      "not analysed" instead of throwing.

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

import { parseHkl } from '../src/hkl-parser.js';
import { analyzeBeamDamage } from '../src/raddose.js';
import { analyzeHkl } from '../src/index.js';

const __dirname = path.dirname(fileURLToPath(import.meta.url));

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
    if (cond) { pass++; console.log(`  ok   ${name}`); }
    else { fail++; console.error(`  FAIL ${name}${detail ? ' — ' + detail : ''}`); }
}

const CELL = { a: 10.0, b: 11.0, c: 12.0, alpha: 90, beta: 90, gamma: 90 };

// Deterministic PRNG (mulberry32) so the synthetic data is reproducible.
function rng(seed) {
    let a = seed >>> 0;
    return () => {
        a |= 0; a = (a + 0x6D2B79F5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

// Build a synthetic unmerged XDS_ASCII scan. Each unique reflection is observed
// at several rotation angles (PSI); its intensity at dose d is
// I0 * exp(-k_true * d) * (1 + noise). `k_true` is the ground-truth decay
// constant (per degree).
function makeSyntheticXds({ nRefl = 400, nObs = 6, kTrue = 0.004, total = 360, seed = 42 } = {}) {
    const rand = rng(seed);
    const lines = [
        '!FORMAT=XDS_ASCII    MERGE=FALSE    FRIEDEL\'S_LAW=TRUE',
        '!OUTPUT_FILE=synthetic.hkl',
        '!DATA_RANGE=       1    1800',
        '!ROTATION_AXIS=  1.0 0.0 0.0',
        '!OSCILLATION_RANGE=  0.2',
        '!STARTING_ANGLE=     0.000',
        '!STARTING_FRAME=       1',
        '!X-RAY_WAVELENGTH= 0.71073',
        '!UNIT_CELL_CONSTANTS= 10.0 11.0 12.0 90.0 90.0 90.0',
        '!NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=12',
        '!ITEM_H=1', '!ITEM_K=2', '!ITEM_L=3', '!ITEM_IOBS=4', '!ITEM_SIGMA(IOBS)=5',
        '!ITEM_XD=6', '!ITEM_YD=7', '!ITEM_ZD=8', '!ITEM_RLP=9', '!ITEM_PEAK=10',
        '!ITEM_CORR=11', '!ITEM_PSI=12',
        '!END_OF_HEADER',
    ];
    for (let r = 0; r < nRefl; r++) {
        const h = Math.floor(rand() * 21) - 10;
        const k = Math.floor(rand() * 21) - 10;
        const l = Math.floor(rand() * 21) - 10;
        if (h === 0 && k === 0 && l === 0) continue;
        const I0 = 100 + rand() * 900;
        for (let o = 0; o < nObs; o++) {
            // PSI is reported centred on 0 (as in real XDS files); the dose
            // position is the same unwrapping the analysis uses:
            // dose = (psi - start) mod total, negatives wrapped by +total.
            const psi = (rand() * total - total / 2);
            let dose = (psi - 0) % total;
            if (dose < 0) dose += total;
            const I = I0 * Math.exp(-kTrue * dose) * (1 + (rand() - 0.5) * 0.1);
            const sig = Math.max(1, 0.02 * I0);
            lines.push(`${h} ${k} ${l} ${I.toFixed(3)} ${sig.toFixed(3)} 0 0 0 0.1 100 0 ${psi.toFixed(3)}`);
        }
    }
    return lines.join('\n') + '\n';
}

console.log('xrdspace beam-damage (RADDOSE-style) tests');
console.log('');

// --- 1. synthetic scan with a known decay constant ---
console.log('synthetic scan (known decay)');
{
    const kTrue = 0.004;
    const text = makeSyntheticXds({ kTrue, nRefl: 400, nObs: 8, seed: 7 });
    const parsed = parseHkl(text);
    check('parser captures PSI on every observation',
        parsed.reflections.every(r => Number.isFinite(r.psi)),
        `${parsed.reflections.filter(r => Number.isFinite(r.psi)).length}/${parsed.reflections.length}`);
    check('parser captures scan geometry (360 deg)',
        parsed.geometry && parsed.geometry.totalRotation === 360,
        JSON.stringify(parsed.geometry));

    const res = analyzeBeamDamage(parsed.reflections, parsed.cell, parsed.geometry, { minIsig: 1, shells: 8 });
    check('analysis is usable', res.usable === true, res.reason || '');
    check('damage detected (ratio < 1)', res.decay === 'yes', `decay=${res.decay}, ratio=${res.overall.ratio}`);

    // Ground-truth overall ratio over the early (first 25%) vs late (last 25%)
    // windows of a full 360-degree rotation.
    const expectedRatio = Math.exp(-kTrue * (0.875 - 0.125) * 360); // late mid - early mid
    const ratio = res.overall.ratio;
    check('overall <I>_late/<I>_early matches ground truth (~10%)',
        ratio > 0.5 * expectedRatio && ratio < 2 * expectedRatio,
        `got ${ratio.toFixed(3)}, expected ~${expectedRatio.toFixed(3)}`);

    // Recovered decay constant should be the same order as k_true.
    const kRec = res.overall.k;
    check('recovered k is the same order as k_true',
        Number.isFinite(kRec) && kRec > 0.3 * kTrue && kRec < 3 * kTrue,
        `got ${kRec}, k_true=${kTrue}`);

    check('per-shell table present (>= 3 shells)', res.shells.length >= 3, `${res.shells.length} shells`);
    check('dose range spans the full rotation',
        res.doseMin < 5 && res.doseMax > 355,
        `[${res.doseMin.toFixed(1)}, ${res.doseMax.toFixed(1)}]`);
}
console.log('');

// --- 1b. a scan with NO decay must not be flagged as damaged ---
console.log('synthetic scan (no decay)');
{
    const text = makeSyntheticXds({ kTrue: 0, nRefl: 300, nObs: 8, seed: 11 });
    const parsed = parseHkl(text);
    const res = analyzeBeamDamage(parsed.reflections, parsed.cell, parsed.geometry, { minIsig: 1, shells: 8 });
    check('analysis is usable', res.usable === true, res.reason || '');
    check('no decay -> ratio ~ 1 (within 3%)',
        Number.isFinite(res.overall.ratio) && Math.abs(res.overall.ratio - 1) < 0.03,
        `ratio=${res.overall.ratio}`);
    check('no decay -> not flagged as damaged', res.decay !== 'yes', `decay=${res.decay}`);
}
console.log('');

// --- 2. real taurine scan ---
console.log('real taurine XDS_ASCII.HKL');
{
    const file = path.join(__dirname, '..', 'test', 'XDS_ASCII.HKL');
    if (fs.existsSync(file)) {
        const parsed = parseHkl(fs.readFileSync(file, 'utf8'));
        check('parser captures PSI', parsed.reflections.every(r => Number.isFinite(r.psi)));
        check('parser captures scan geometry',
            parsed.geometry && Number.isFinite(parsed.geometry.totalRotation),
            JSON.stringify(parsed.geometry));

        const res = analyzeBeamDamage(parsed.reflections, parsed.cell, parsed.geometry, { minIsig: 2, shells: 10 });
        check('analysis is usable', res.usable === true, res.reason || '');
        check('nObs > 0 and nRefl > 0', res.nObs > 0 && res.nRefl > 0,
            `nObs=${res.nObs}, nRefl=${res.nRefl}`);
        check('overall ratio is sensible (0 < R < 2)',
            Number.isFinite(res.overall.ratio) && res.overall.ratio > 0 && res.overall.ratio < 2,
            `ratio=${res.overall.ratio}`);
        check('shells cover the resolution range', res.shells.length >= 5, `${res.shells.length} shells`);

        // End-to-end through analyzeHkl with --rad.
        const full = analyzeHkl(fs.readFileSync(file, 'utf8'), { rad: true, quality: true });
        check('analyzeHkl({rad:true}) attaches merge.beamDamage',
            full.ok && full.merge && full.merge.beamDamage && full.merge.beamDamage.usable === true);
    } else {
        check('taurine file present', false, file + ' not found');
    }
}
console.log('');

// --- 3. graceful degradation on input without PSI ---
console.log('input without PSI (graceful)');
{
    const shelx = [
        '0 0 1 100.0 5.0', '0 0 -1 98.0 5.0', '1 0 0 90.0 4.0',
        '0 1 0 85.0 4.0', '-1 0 0 80.0 4.0', '0 -1 0 78.0 4.0',
    ];
    for (let i = 0; i < 200; i++) shelx.push(`${i % 7} ${(i * 3) % 7} ${(i * 5) % 7} ${50 + i} 3.0`);
    const parsed = parseHkl(shelx.join('\n'));
    check('SHELX input has no PSI', parsed.reflections.every(r => r.psi === undefined));
    const res = analyzeBeamDamage(parsed.reflections, CELL, null, {});
    check('reports not analysed (no PSI)', res.usable === false && /rotation angle/.test(res.reason || ''),
        JSON.stringify({ usable: res.usable, reason: res.reason }));

    // And through analyzeHkl it must not throw.
    const full = analyzeHkl(shelx.join('\n'), { cell: CELL, rad: true });
    check('analyzeHkl({rad:true}) on PSI-less input does not throw',
        full.ok === true && full.merge.beamDamage.usable === false);
}

console.log('');
console.log(`passed ${pass}, failed ${fail}`);
process.exit(fail === 0 ? 0 : 1);
