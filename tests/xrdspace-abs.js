// Copyright (c) 2026 Denis Spasyuk. MIT License.
// Test for multi-scan absorption correction (src/abs.js), AIMLESS method.
//
// Verifies: the spherical-harmonic basis (term count, orthonormality on an
// equal-area grid), the full multi-scan correction on synthetic data with a
// known absorption gradient and per-scan crystal orientations, error paths,
// and the class API.

import { runAbsorptionCorrection, buildAbsorptionReport, MultiScanAbsorptionCorr, SphericalHarmonic } from '../src/abs.js';

let failures = 0;
function check(name, cond) {
    if (cond) console.log(`  PASS  ${name}`);
    else { console.log(`  FAIL  ${name}`); failures++; }
}

const laueP1 = [{ M: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: [0, 0, 0], opString: 'x, y, z' }];
const laueInv = [
    { M: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: [0, 0, 0], opString: 'x, y, z' },
    { M: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]], t: [0, 0, 0], opString: '-x, -y, -z' },
];
const cell = { a: 10, b: 10, c: 10, alpha: 90, beta: 90, gamma: 90 };

// Rotation about the z axis by `deg` (crystal->lab orientation).
function rotZ(deg) {
    const t = deg * Math.PI / 180, c = Math.cos(t), s = Math.sin(t);
    return [[c, -s, 0], [s, c, 0], [0, 0, 1]];
}
// Rotation about the x axis by `deg`.
function rotX(deg) {
    const t = deg * Math.PI / 180, c = Math.cos(t), s = Math.sin(t);
    return [[1, 0, 0], [0, c, -s], [0, s, c]];
}
// Rotation about the y axis by `deg`.
function rotY(deg) {
    const t = deg * Math.PI / 180, c = Math.cos(t), s = Math.sin(t);
    return [[c, 0, s], [0, 1, 0], [-s, 0, c]];
}

// --- 1. Spherical-harmonic basis -------------------------------------------
console.log('Spherical harmonics:');
{
    // AIMLESS defaults: lmaxEven=6, lmaxOdd=5 -> 1 + sum_{l=1..5}(2l+1) = 49.
    const sh = new SphericalHarmonic(6, 5);
    check('nTerms() = 49 for (6,5)', sh.nTerms() === 49);

    // Y_00 is the constant 1/sqrt(4 pi).
    check('Y_00 = 1/sqrt(4*pi)', Math.abs(sh.ylm(1, 0, 0)[0] - 1 / Math.sqrt(4 * Math.PI)) < 1e-12);

    // Orthonormality on an equal-area (midpoint) grid:
    //   (1/4pi) int Y_k Y_j dOmega ~= delta_kj
    // A latitude/longitude midpoint grid with the sin(theta) cell weight is an
    // equal-area partition, so a plain average of the products converges to
    // delta_kj as the grid refines.
    // Orthonormality: int Y_k Y_j dOmega ~ delta_kj. A midpoint lat/lon grid
    // with cell weight sin(theta) dtheta dphi is an equal-area partition, so
    // summing f_i * sin(theta_i) dtheta dphi over the grid estimates the
    // surface integral directly (no extra 1/4pi).
    const nLat = 90, nLon = 180;
    const dth = Math.PI / nLat, dph = (2 * Math.PI) / nLon;
    const dirs = [];
    const wts = [];
    for (let i = 0; i < nLat; i++) {
        const th = (i + 0.5) * dth;
        for (let j = 0; j < nLon; j++) {
            const ph = j * dph;
            dirs.push([Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph), Math.cos(th)]);
            wts.push(Math.sin(th) * dth * dph);
        }
    }
    const n = sh.nTerms();
    let maxOff = 0, maxDiag = 0, minDiag = Infinity;
    for (let k = 0; k < n; k++) {
        const vk = new Float64Array(nLat * nLon);
        for (let i = 0; i < dirs.length; i++) vk[i] = sh.ylm(dirs[i][0], dirs[i][1], dirs[i][2])[k];
        for (let j = 0; j < n; j++) {
            let s = 0;
            for (let i = 0; i < dirs.length; i++) {
                s += vk[i] * sh.ylm(dirs[i][0], dirs[i][1], dirs[i][2])[j] * wts[i];
            }
            if (k === j) { if (s > maxDiag) maxDiag = s; if (s < minDiag) minDiag = s; }
            else maxOff = Math.max(maxOff, Math.abs(s));
        }
    }
    check(`basis orthonormal (diag in [${minDiag.toFixed(3)},${maxDiag.toFixed(3)}], off-diag max ${maxOff.toFixed(4)})`,
        maxDiag < 1.05 && minDiag > 0.95 && maxOff < 0.02);

    // l=1, m=0 term is an odd function of z (proportional to cos theta).
    // Basis order: [0]=Y00, [1]=Y11cos, [2]=Y11sin, [3]=Y10, ...
    const sh1 = new SphericalHarmonic(2, 1);
    const v = sh1.ylm(0, 0, 1); // +z
    check('Y_10 odd in z (Y_10(+z) = -Y_10(-z))',
        v[3] > 0 && sh1.ylm(0, 0, -1)[3] < 0 && Math.abs(v[3] + sh1.ylm(0, 0, -1)[3]) < 1e-9);
}

// --- 2. Synthetic multi-scan correction ------------------------------------
console.log('Multi-scan correction:');
{
    // 3 scans x 400 reflections, P1 Laue (hkl and -hkl are distinct).
    // True absorption depends on the LAB-frame z of the scattering vector:
    //   T = 1 + 0.30 * z_lab
    // Each scan has a different crystal orientation, so the same hkl is
    // measured along a different beam path in each scan.
    let seed = 42;
    const rand = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
    // Shared set of CRYSTAL hkl vectors. The same crystal reflection is
    // observed in every scan (the crystal is the same, just re-oriented), so
    // under P1 each (h,k,l) forms a multi-scan equivalent set. In scan s the
    // lab-frame beam path is O_s * (h,k,l), so the absorption (and hence the
    // observed intensity) differs from scan to scan -- exactly the redundancy
    // multi-scan absorption correction exploits.
    // Three distinct orientations about different axes so the lab-frame
    // direction (and thus the absorption) genuinely differs per scan.
    const orientations = [rotZ(0), rotX(60), rotY(60)];
    const shared = [];
    for (let i = 0; i < 300; i++) {
        const h = Math.floor(rand() * 30) - 15;
        const k = Math.floor(rand() * 30) - 15;
        const l = Math.floor(rand() * 30) - 15;
        if (h === 0 && k === 0 && l === 0) continue;
        shared.push([h, k, l]);
    }
    const scans = [];
    for (let s = 0; s < 3; s++) {
        const O = orientations[s];
        const refls = [];
        for (const [h, k, l] of shared) {
            const mag = Math.sqrt(h * h + k * k + l * l);
            const zlab = (O[2][0] * h + O[2][1] * k + O[2][2] * l) / mag;
            const T = 1 + 0.30 * zlab;
            if (!(T > 0.2)) continue;
            refls.push({ h, k, l, I: 1000 / T, sig: 10 });
        }
        scans.push({ reflections: refls, orientation: O });
    }

    const corr = new MultiScanAbsorptionCorr({
        lmaxEven: 4, lmaxOdd: 3,
        maxIterations: 40, convergenceThreshold: 1e-3,
        minScans: 2, minMembers: 2, minGroups: 5, useSD: false,
    });
    corr.setCell(cell);
    scans.forEach((sc, i) => corr.addScan(i, sc.reflections, sc.orientation));
    const result = corr.run(laueP1);

    check('returns ok', result.ok === true);
    if (result.ok) {
        check('nTerms for (4,3) is 25', result.nTerms === 25);
        check('converges', result.converged === true);
        check('overall scale finite and positive', Number.isFinite(result.overallScale) && result.overallScale > 0);
        check('absorption + run-offset coefficients (25 + 2)',
            Array.isArray(result.absorptionCoeffs) && result.absorptionCoeffs.length === 27);
        check('two run offsets (runs 1,2; run 0 reference)', Object.keys(result.runOffsets).length === 2);
        check('run offsets ~ 0 (all scans share the same true intensities)',
            Object.values(result.runOffsets).every(o => Math.abs(o) < 0.02));
        check('corrected reflections returned (3 scans x 300)',
            Array.isArray(corr.getCorrectedReflections()) && corr.getCorrectedReflections().length === 900);
        check('corrected values finite',
            corr.getCorrectedReflections().every(r => Number.isFinite(r.I) && Number.isFinite(r.sig)));
        check('stats present', result.stats && result.stats.nReflections === 900);
        // The fitted l=1 terms should recover the true z-gradient (0.30):
        // l1mag * Y10-normalisation (sqrt(3/4pi) = 0.4886) ~ 0.30.
        const c = result.absorptionCoeffs;
        const l1mag = Math.hypot(c[1], c[2], c[3]);
        const grad = l1mag * Math.sqrt(3 / (4 * Math.PI));
        check(`l=1 terms recover the 0.30 gradient (got ${grad.toFixed(3)})`,
            Math.abs(grad - 0.30) < 0.05);
        console.log(`  info  overallScale=${result.overallScale.toFixed(4)} l1mag=${l1mag.toFixed(4)} ` +
            `converged=${result.converged} iters=${result.iterations} groups=${result.nEquivalentSets}`);
    } else {
        console.log('  error:', result.error, result.detail);
    }
    console.log('');
    console.log(buildAbsorptionReport(result, cell));
}

// --- 2b. Higher-order (quadratic) absorption gradient ----------------------
console.log('Quadratic gradient (l=2):');
{
    // T = 1 + 0.30*z + 0.25*(z^2 - 1/3)  (z = lab z-direction). The z^2 term is
    // an l=2 harmonic, so the l=2 coefficients must pick it up.
    let seed = 7;
    const rand = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
    const orientations = [rotZ(0), rotX(50), rotY(70)];
    const shared = [];
    for (let i = 0; i < 400; i++) {
        const h = Math.floor(rand() * 24) - 12;
        const k = Math.floor(rand() * 24) - 12;
        const l = Math.floor(rand() * 24) - 12;
        if (h || k || l) shared.push([h, k, l]);
    }
    const scans = [];
    for (let s = 0; s < 3; s++) {
        const O = orientations[s];
        const refls = shared.map(([h, k, l]) => {
            const mag = Math.sqrt(h * h + k * k + l * l);
            const z = (O[2][0] * h + O[2][1] * k + O[2][2] * l) / mag;
            return { h, k, l, I: 1000 / (1 + 0.30 * z + 0.25 * (z * z - 1 / 3)), sig: 10 };
        });
        scans.push({ reflections: refls, orientation: O });
    }
    const corr = new MultiScanAbsorptionCorr({
        lmaxEven: 4, lmaxOdd: 3, maxIterations: 40, convergenceThreshold: 1e-3,
        minGroups: 5, useSD: false,
    });
    corr.setCell(cell);
    scans.forEach((sc, i) => corr.addScan(i, sc.reflections, sc.orientation));
    const result = corr.run(laueP1);
    check('quadratic fit ok', result.ok === true);
    if (result.ok) {
        const c = result.absorptionCoeffs;
        // The l=2, m=0 coefficient of -log(T) for
        //   T = 1 + 0.30 z + 0.25 (z^2 - 1/3)
        // is ~ -0.217 (the P2(cos) component of -log T). l=2 m=0 sits at
        // index 8 (Y00, l1: 1,2,3, l2: 4,5,6,7,8).
        const l2m0 = c[8];
        check(`l=2 quadratic term recovered (l2m0=${l2m0.toFixed(3)}, theory ~-0.217)`,
            Math.abs(l2m0 - (-0.217)) < 0.06);
        console.log(`  info  quadratic fit converged=${result.converged} iters=${result.iterations}`);
    }
}

// --- 2c. Outlier rejection -------------------------------------------------
console.log('Outlier rejection:');
{
    let seed = 99;
    const rand = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
    // Box-Muller Gaussian from the seeded uniform PRNG.
    const gauss = () => {
        const u1 = Math.max(1e-12, rand()), u2 = rand();
        return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    };
    const orientations = [rotZ(0), rotX(60), rotY(60)];
    const shared = [];
    for (let i = 0; i < 300; i++) {
        const h = Math.floor(rand() * 30) - 15;
        const k = Math.floor(rand() * 30) - 15;
        const l = Math.floor(rand() * 30) - 15;
        if (h || k || l) shared.push([h, k, l]);
    }
    const scans = [];
    let nBad = 0;
    for (let s = 0; s < 3; s++) {
        const O = orientations[s];
        const refls = shared.map(([h, k, l]) => {
            const mag = Math.sqrt(h * h + k * k + l * l);
            const z = (O[2][0] * h + O[2][1] * k + O[2][2] * l) / mag;
            // Realistic Poisson-like noise (sigma = 10) plus ~4% gross spikes.
            let I = 1000 / (1 + 0.30 * z) + gauss() * 10;
            if (rand() < 0.04) { I *= 6; nBad++; }
            return { h, k, l, I, sig: 10 };
        });
        scans.push({ reflections: refls, orientation: O });
    }
    const corr = new MultiScanAbsorptionCorr({
        lmaxEven: 4, lmaxOdd: 3, maxIterations: 40, convergenceThreshold: 1e-3,
        minGroups: 5, useSD: false, outlierSigma: 5,
    });
    corr.setCell(cell);
    scans.forEach((sc, i) => corr.addScan(i, sc.reflections, sc.orientation));
    const result = corr.run(laueP1);
    check('outlier fit ok', result.ok === true);
    if (result.ok) {
        // The gross spikes (~5 sigma in log space) are caught. Rejection is
        // conservative: a spike briefly biases its group's mean, so some of its
        // clean equivalent-mates are rejected in the same pass too (they are
        // re-admitted on later passes, but a few may remain flagged). What
        // matters is that most injected outliers are caught and the fit is
        // unharmed.
        check(`caught most gross outliers (rejected=${result.nRejected} of ~${nBad} injected)`,
            result.nRejected >= 0.5 * nBad && result.nRejected <= 5 * nBad);
        const c = result.absorptionCoeffs;
        const grad = Math.hypot(c[1], c[2], c[3]) * Math.sqrt(3 / (4 * Math.PI));
        check(`gradient still recovered despite outliers (got ${grad.toFixed(3)})`,
            Math.abs(grad - 0.30) < 0.06);
    }
}

// --- 3. Insufficient scans ---------------------------------------------------
console.log('Insufficient scans:');
{
    const result = runAbsorptionCorrection(
        [[{ h: 1, k: 0, l: 0, I: 100, sig: 5 }]], cell, laueP1);
    check('rejects single scan', result.ok === false && result.error === 'NEED_MORE_SCANS');
}

// --- 4. Object-form scans + inversion Laue ----------------------------------
console.log('Object-form scans:');
{
    const scansObj = {
        0: [
            { h: 1, k: 2, l: 3, I: 500, sig: 10 },
            { h: 4, k: 0, l: 0, I: 800, sig: 10 },
            { h: 0, k: 5, l: 0, I: 600, sig: 10 },
        ],
        1: [
            { h: 1, k: 2, l: 3, I: 480, sig: 10 },
            { h: -4, k: 0, l: 0, I: 790, sig: 10 },
            { h: 0, k: -5, l: 0, I: 610, sig: 10 },
        ],
    };
    const result = runAbsorptionCorrection(scansObj, cell, laueInv, { minGroups: 3, useSD: false });
    check('object-form scans accepted', result.ok === true);
    if (result.ok) check('3 equivalent sets found', result.nEquivalentSets === 3);
}

// --- 4b. P4P ORT orientation -----------------------------------------------
console.log('P4P ORT orientation:');
{
    const fs = await import('node:fs');
    const path = 'test/p4p/ds74.p4p';
    if (fs.existsSync(path)) {
        const { parseHkl } = await import('../src/hkl-parser.js');
        const { dSpacing } = await import('../src/merge.js');
        const parsed = parseHkl(fs.readFileSync(path, 'utf8'));
        check('P4P orientation (ORT matrix) extracted',
            Array.isArray(parsed.orientation) && parsed.orientation.length === 3);
        // |ORT . (h,k,l)| must equal 1/d(hkl) -- ORT is the reciprocal
        // orientation matrix mapping Miller indices to the lab scattering vector.
        const O = parsed.orientation;
        let maxErr = 0, n = 0;
        for (const r of parsed.reflections) {
            if (!r.h && !r.k && !r.l) continue;
            const qx = O[0][0] * r.h + O[0][1] * r.k + O[0][2] * r.l;
            const qy = O[1][0] * r.h + O[1][1] * r.k + O[1][2] * r.l;
            const qz = O[2][0] * r.h + O[2][1] * r.k + O[2][2] * r.l;
            const mag = Math.sqrt(qx * qx + qy * qy + qz * qz);
            const d = dSpacing(r.h, r.k, r.l, parsed.cell);
            if (!(d > 0)) continue;
            maxErr = Math.max(maxErr, Math.abs(mag - 1 / d));
            n++;
        }
        check(`|ORT.(h,k,l)| = 1/d over ${n} reflections (max err ${maxErr.toExponential(1)})`, maxErr < 1e-5);
    } else {
        console.log('  SKIP  P4P sample file not found');
    }
}

// --- 5. Direct class API -----------------------------------------------------
console.log('Direct class API:');
{
    const corr = new MultiScanAbsorptionCorr({ lmaxEven: 2, lmaxOdd: 1, useSD: false, minGroups: 3 });
    corr.setCell(cell);
    corr.addScan(0, [
        { h: 1, k: 0, l: 0, I: 1000, sig: 10 },
        { h: 0, k: 1, l: 0, I: 900, sig: 10 },
        { h: 0, k: 0, l: 1, I: 1100, sig: 10 },
        { h: 2, k: 1, l: 0, I: 500, sig: 10 },
    ]);
    corr.addScan(1, [
        { h: -1, k: 0, l: 0, I: 1050, sig: 10 },
        { h: 0, k: -1, l: 0, I: 870, sig: 10 },
        { h: 0, k: 0, l: -1, I: 1080, sig: 10 },
        { h: -2, k: -1, l: 0, I: 520, sig: 10 },
    ]);
    const result = corr.run(laueInv);
    check('class run ok', result.ok === true);
    if (result.ok) {
        const corrected = corr.getCorrectedReflections();
        check('getCorrectedReflections returns all 8', corrected.length === 8);
        check('corrected intensities are finite',
            corrected.every(r => Number.isFinite(r.I) && Number.isFinite(r.sig)));
        check('nTerms for (2,1) is 9', result.nTerms === 9);
    }
}

console.log('');
if (failures) {
    console.log(`${failures} test(s) FAILED`);
    process.exit(1);
} else {
    console.log('All tests passed.');
}
