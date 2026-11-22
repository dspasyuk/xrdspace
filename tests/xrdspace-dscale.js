#!/usr/bin/env node
// Copyright (c) 2026 Denis Spasyuk. MIT License.
// Tests for multi-dataset merging: relative dataset scaling (src/dataset-scale.js),
// the MTZ/text auto-detecting loader (parseReflectionFile), the analyzeParsed
// integration (merge.datasets summary), and — when the reference data is present —
// an end-to-end CLI merge of the two-kappa XDS_ASCII set.
//
// Usage:
//   node tests/xrdspace-dscale.js [path-to-dir-with-two-XDS_ASCII-HKL-files]

import fs from 'node:fs';
import path from 'node:path';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

import { compareDatasetCells, datasetCount, scaleDatasets } from '../src/dataset-scale.js';
import { parseReflectionFile, analyzeParsed } from '../src/index.js';
import { writeMtz } from '../src/mtz.js';
import { dSpacing as dSpacingLocal } from '../src/merge.js';

const __dirname = path.dirname(fileURLToPath(import.meta.url));

let failures = 0;
function check(name, cond) {
    if (cond) console.log(`  PASS  ${name}`);
    else { console.log(`  FAIL  ${name}`); failures++; }
}
const near = (a, b, tol) => Number.isFinite(a) && Math.abs(a - b) <= tol;

// Deterministic RNG (LCG) + Box-Muller gaussian.
let seed = 12345;
function rnd() { seed = (seed * 1103515245 + 12345) >>> 0; return seed / 4294967296; }
function gauss() {
    const u = Math.max(rnd(), 1e-12), v = rnd();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

const laueInv = [
    { M: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: [0, 0, 0], opString: 'x, y, z' },
    { M: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]], t: [0, 0, 0], opString: '-x, -y, -z' },
];
const cell = { a: 8, b: 10, c: 12, alpha: 95, beta: 100, gamma: 105 };

// ---------------------------------------------------------------------------
console.log('compareDatasetCells:');
{
    const ref = { a: 10, b: 10, c: 10, alpha: 90, beta: 90, gamma: 90 };
    const r0 = compareDatasetCells(ref, [ref, { ...ref }]);
    check('identical cells are ok', r0.ok && r0.warnings.length === 0);
    check('identical cells have zero deviation', r0.rows[1].maxLenPct === 0 && r0.rows[1].maxAngDeg === 0);

    const r1 = compareDatasetCells(ref, [ref, { a: 10.05, b: 10, c: 10, alpha: 90, beta: 90, gamma: 90 }]);
    check('0.5% length deviation within default tolerance', r1.ok && near(r1.rows[1].maxLenPct, 0.5, 1e-9));

    const r2 = compareDatasetCells(ref, [ref, { a: 10.2, b: 10, c: 10, alpha: 90, beta: 90, gamma: 90 }]);
    check('2% length deviation is flagged', !r2.ok && r2.rows[1].bad && r2.warnings.length === 1);

    const r3 = compareDatasetCells(ref, [ref, { a: 10, b: 10, c: 10, alpha: 92, beta: 90, gamma: 90 }]);
    check('2 deg angle deviation is flagged', !r3.ok && near(r3.rows[1].maxAngDeg, 2, 1e-9));

    const r4 = compareDatasetCells(ref, [ref, null]);
    check('missing cell tolerated (null row, ok)', r4.ok && r4.rows[1] === null);

    const r5 = compareDatasetCells(ref, [ref, { a: 10.15, b: 10, c: 10, alpha: 90, beta: 90, gamma: 90 }], { tolLen: 0.02 });
    check('custom tolerance 2% passes 1.5% deviation', r5.ok);
}

// ---------------------------------------------------------------------------
console.log('\nscaleDatasets on synthetic two-dataset P-1 data:');
// Truth model: I(s) = 2000*exp(-1.2 s) + 200, s = 1/d^2. Dataset 1 carries a
// known scale K_T and temperature factor B_T on top (plus 2% gaussian noise on
// both datasets). The fit must recover the RELATIVE correction:
//   K_1/K_0 ~= 1/K_T   and   B_1 - B_0 ~= -B_T.
const K_T = 1.05, B_T = 0.5, NOISE = 0.02;

function makeSynthetic() {
    const orbits = [];
    const seen = new Set();
    for (let h = -5; h <= 5; h++) {
        for (let k = -6; k <= 6; k++) {
            for (let l = -7; l <= 7; l++) {
                if (!h && !k && !l) continue;
                const d = dSpacingLocal(h, k, l, cell);
                if (!(d >= 0.9 && d <= 8)) continue;
                const rep = h > 0 || (h === 0 && (k > 0 || (k === 0 && l > 0))) ? [h, k, l] : [-h, -k, -l];
                const key = rep.join(',');
                if (seen.has(key)) continue;
                seen.add(key);
                orbits.push({ rep, s: 1 / (d * d), I: 2000 * Math.exp(-1.2 / (d * d)) + 200 });
            }
        }
    }
    const refl = [];
    for (const o of orbits) {
        // Dataset 0: both Friedel mates with 2% noise.
        for (const sgn of [1, -1]) {
            const I = o.I * (1 + NOISE * gauss());
            refl.push({ h: sgn * o.rep[0], k: sgn * o.rep[1], l: sgn * o.rep[2], I, sig: 0.03 * o.I + Math.sqrt(Math.abs(I)), dataset: 0 });
        }
        // Dataset 1: known scale + B on top.
        const g = K_T * Math.exp(-B_T * o.s / 2);
        for (const sgn of [1, -1]) {
            const I = o.I * g * (1 + NOISE * gauss());
            refl.push({ h: sgn * o.rep[0], k: sgn * o.rep[1], l: sgn * o.rep[2], I, sig: 0.03 * o.I * g + Math.sqrt(Math.abs(I)), dataset: 1 });
        }
    }
    return { orbits, refl };
}

const synth = makeSynthetic();
console.log(`  (synthetic: ${synth.orbits.length} orbits, ${synth.refl.length} observations)`);
{
    check('datasetCount reads the tags', datasetCount(synth.refl) === 2);

    const res = scaleDatasets(synth.refl, laueInv, cell, { cycles: 10 });
    check('scaling succeeded', res.ok === true);
    check('converged', res.converged === true);
    const [d0, d1] = res.datasets;
    const relK = d1.K / d0.K;
    const relB = d1.B - d0.B;
    check(`relative K recovers 1/K_T=${(1 / K_T).toFixed(4)} (got ${relK.toFixed(4)})`, near(relK, 1 / K_T, 0.01));
    check(`relative B recovers -B_T=${(-B_T).toFixed(3)} (got ${relB.toFixed(3)})`, near(relB, -B_T, 0.05));
    check('anchor: obs-weighted mean lnK = 0',
        near(d0.nObs * Math.log(d0.K) + d1.nObs * Math.log(d1.K), 0, 1e-6));
    check('rAfter < rBefore for dataset 0', d0.rAfter < d0.rBefore);
    check('rAfter < rBefore for dataset 1', d1.rAfter < d1.rBefore);
    check('rAfter at the noise level (< 60% of rBefore)', d0.rAfter < 0.6 * d0.rBefore && d1.rAfter < 0.6 * d1.rBefore);
    check('all observations returned', res.reflections.length === synth.refl.length);

    // Each corrected intensity equals I * K_d * exp(-B_d s / 2).
    let maxDev = 0;
    for (let i = 0; i < synth.refl.length; i++) {
        const r0 = synth.refl[i], r1 = res.reflections[i];
        // order may differ only if passthroughs exist; none here (all d finite).
        const ds = r0.dataset;
        const d = dSpacingLocal(r0.h, r0.k, r0.l, cell);
        const c = res.datasets[ds].K * Math.exp(-res.datasets[ds].B / (2 * d * d));
        if (r1.h !== r0.h || r1.k !== r0.k || r1.l !== r0.l) { maxDev = Infinity; break; }
        maxDev = Math.max(maxDev, Math.abs(r1.I - r0.I * c));
        maxDev = Math.max(maxDev, Math.abs(r1.sig - Math.abs(r0.sig) * c));
    }
    check('corrected I and sig carry the fitted correction', maxDev < 1e-9);

    // K-only fit: a constant cannot absorb the B_T gradient (it just shifts the
    // mismatch across resolution), so K+B must give clearly better agreement.
    const resK = scaleDatasets(synth.refl, laueInv, cell, { cycles: 10, withB: false });
    check('K-only mode: B stays 0', resK.datasets[0].B === 0 && resK.datasets[1].B === 0);
    check('K+B agrees better than K-only (B term matters)',
        resK.datasets[0].rAfter + resK.datasets[1].rAfter
        > res.datasets[0].rAfter + res.datasets[1].rAfter);

    // Statistics-only mode (cycles = 0): identity corrections, R unchanged.
    const res0 = scaleDatasets(synth.refl, laueInv, cell, { cycles: 0 });
    check('cycles=0 keeps K=1, B=0', res0.datasets.every(d => d.K === 1 && d.B === 0));
    check('cycles=0 rAfter == rBefore', res0.datasets.every(d => d.rAfter === d.rBefore));
    check('cycles=0 returns identity copies',
        res0.reflections.every((r, i) => r.I === synth.refl[i].I && r.dataset === synth.refl[i].dataset));

    // Error path: a single dataset.
    const single = synth.refl.map(r => ({ ...r, dataset: 0 }));
    const resS = scaleDatasets(single, laueInv, cell, {});
    check('single dataset -> ok:false', resS.ok === false && /fewer than two/i.test(resS.reason));
}

// ---------------------------------------------------------------------------
console.log('\nparseReflectionFile (MTZ auto-detection):');
function syntheticMtz(iLabel, sLabel) {
    const ncol = 5;
    const refls = [
        [0, 0, 1, 120.5, 3.2],
        [0, 0, 2, 88.1, 2.1],
        [1, 0, 0, 45.0, 1.7],
        [0, 1, 0, NaN, NaN], // missing reflection -> skipped
        [1, 1, 0, 200.7, 4.4],
    ];
    const data = new Float32Array(ncol * refls.length);
    for (let r = 0; r < refls.length; r++) for (let c = 0; c < ncol; c++) data[r * ncol + c] = refls[r][c];
    return {
        title: 'dscale test',
        version: 'MTZ:V1.1',
        nreflections: refls.length,
        ncols: ncol,
        cell: { a: 50, b: 50, c: 50, alpha: 90, beta: 90, gamma: 90 },
        spaceGroupNumber: 1,
        spaceGroupName: 'P 1',
        latticeType: 'P',
        nsymop: 1, nops: 1, symops: ['X, Y, Z'],
        reso: [4, 100], valm: NaN,
        columns: [
            { label: 'H', type: 'H', min: 0, max: 1, dataset: 0, source: '' },
            { label: 'K', type: 'H', min: 0, max: 1, dataset: 0, source: '' },
            { label: 'L', type: 'H', min: 0, max: 2, dataset: 0, source: '' },
            { label: iLabel, type: 'J', min: 45, max: 200.7, dataset: 0, source: '' },
            { label: sLabel, type: 'Q', min: 1.7, max: 4.4, dataset: 0, source: '' },
        ],
        datasets: [{ id: 0, project: 'P', crystal: 'C', dataset: 'D', cell: { a: 50, b: 50, c: 50, alpha: 90, beta: 90, gamma: 90 }, wavelength: 0.6888 }],
        history: ['xrdspace dscale test'],
        batchNumbers: [],
        data,
        byteOrder: 'LE',
    };
}
{
    const bufMerged = writeMtz(syntheticMtz('IMEAN', 'SIGIMEAN'));
    const p = parseReflectionFile(bufMerged);
    check('MTZ magic detected -> format mtz', p.format === 'mtz');
    check('MTZ cell parsed', p.cell.a === 50 && p.cell.alpha === 90);
    check('MTZ wavelength from dataset record', near(p.wavelength, 0.6888, 1e-6));
    check('IMEAN columns -> merge=true', p.merge === true);
    check('NaN row skipped, 4 reflections', p.reflections.length === 4);
    check('I/sig values read', near(p.reflections[0].I, 120.5, 1e-3) && near(p.reflections[0].sig, 3.2, 1e-3));

    const bufUnmerged = writeMtz(syntheticMtz('I', 'SIGI'));
    const p2 = parseReflectionFile(bufUnmerged);
    check('I/SIGI columns -> merge=false (unmerged MTZ)', p2.merge === false);
    check('unmerged MTZ reflections read', p2.reflections.length === 4 && near(p2.reflections[1].I, 88.1, 1e-3));

    // Explicit column selection overrides the label preference.
    const p3 = parseReflectionFile(bufMerged, { intensity: 'SIGIMEAN', sigma: 'IMEAN' });
    check('--intensity/--sigma style column override', near(p3.reflections[0].I, 3.2, 1e-3) && near(p3.reflections[0].sig, 120.5, 1e-3));

    // Text passthrough: the bundled XDS_ASCII test file.
    const txtPath = path.join(__dirname, '..', 'test', 'XDS_ASCII.HKL');
    if (fs.existsSync(txtPath)) {
        const p4 = parseReflectionFile(fs.readFileSync(txtPath));
        check('text file -> xds_ascii', p4.format === 'xds_ascii');
        check('text reflections parsed', p4.reflections.length > 1000);
    } else {
        console.log('  SKIP  test/XDS_ASCII.HKL not present');
    }
}

// ---------------------------------------------------------------------------
console.log('\nanalyzeParsed multi-dataset integration:');
{
    const meta = [
        { name: 'ds0.hkl', format: 'xds_ascii', cell, nObs: synth.refl.filter(r => r.dataset === 0).length },
        { name: 'ds1.hkl', format: 'xds_ascii', cell: { ...cell, a: cell.a * 1.002 }, nObs: synth.refl.filter(r => r.dataset === 1).length },
    ];
    const parsed = {
        format: 'xds_ascii', title: 'synthetic', cell, wavelength: 0.7,
        merge: false, friedelsLaw: true, reflections: synth.refl,
    };
    const res = analyzeParsed(parsed, { spaceGroup: 2, quality: false, datasets: meta });
    check('analysis ok', res.ok === true);
    const ds = res.merge && res.merge.datasets;
    check('merge.datasets summary present', !!ds && ds.count === 2);
    check('scaling applied by default', /applied/.test(ds.scaling));
    check('per-dataset K/B attached', ds.items.every(it => it.K != null && Number.isFinite(it.B)));
    check('cell deviation attached (0.2% a-axis)', near(ds.items[1].cellDev.maxLenPct, 0.2, 0.01));
    check('no cell warnings at 0.2%', ds.cellCheck.ok === true);
    check('merged intensities use the scaled data', res.merge.statistics.nObs === synth.refl.length);

    const res2 = analyzeParsed(parsed, { spaceGroup: 2, quality: false, datasets: meta, dscale: false });
    check('--no-dscale: state says disabled', /disabled/.test(res2.merge.datasets.scaling));
    check('--no-dscale: K=1 identity', res2.merge.datasets.items.every(it => it.K === 1 && it.B === 0));

    // A bad cell (3% off) must raise a warning in the summary.
    const metaBad = [meta[0], { ...meta[1], cell: { ...cell, b: cell.b * 1.03 } }];
    const res3 = analyzeParsed(parsed, { spaceGroup: 2, quality: false, datasets: metaBad });
    check('3% cell deviation flagged', res3.merge.datasets.cellCheck.ok === false
        && res3.merge.datasets.cellCheck.warnings.length === 1);
}

// ---------------------------------------------------------------------------
// End-to-end CLI run on the reference two-kappa dataset (when present).
console.log('\nEnd-to-end CLI on reference kappa-pair data:');
{
    const dir = process.argv[2] || '/home/denis/X-ray/UofS/MiguelSoto';
    const f1 = path.join(dir, 'XDS_ASCII-1.HKL');
    const f2 = path.join(dir, 'XDS_ASCII-2.HKL');
    if (!fs.existsSync(f1) || !fs.existsSync(f2)) {
        console.log(`  SKIP  ${f1} / ${f2} not present`);
    } else {
        const cli = path.join(__dirname, '..', 'src', 'xrdspace.js');
        const t0 = Date.now();
        const run = spawnSync(process.execPath, [cli, '--hklin', f1, f2, '--no-write'], { encoding: 'utf8', timeout: 600000 });
        check('CLI exits 0', run.status === 0);
        const out = run.stdout || '';
        check('report has MULTI-DATASET MERGE (2 DATASETS)', out.includes('MULTI-DATASET MERGE (2 DATASETS)'));
        check('scaling converged', /Relative scaling.*converged/.test(out) && !out.includes('NOT converged'));
        check('cells compatible', /Cell agreement.*compatible/.test(out));
        // Both datasets must show improved inter-dataset agreement.
        const rows = out.split('\n').filter(l => /XDS_ASCII-[12]\.HKL/.test(l) && /->/.test(l));
        check('per-dataset R-vs-others rows present', rows.length === 2);
        let improved = rows.length === 2;
        for (const row of rows) {
            const m = row.match(/([\d.]+)%\s*->\s*([\d.]+)%/);
            if (!m || !(parseFloat(m[2]) < parseFloat(m[1]))) improved = false;
        }
        check('R-vs-others improved after scaling for both datasets', improved);
        check('space group P -1 determined', /Best space group\s*:\s*P -1/.test(out));
        console.log(`  (CLI run took ${((Date.now() - t0) / 1000).toFixed(1)} s)`);
    }
}

console.log('\n' + '='.repeat(50));
if (failures) {
    console.log(`${failures} check(s) FAILED`);
    process.exit(1);
}
console.log('All multi-dataset merge checks passed.');
