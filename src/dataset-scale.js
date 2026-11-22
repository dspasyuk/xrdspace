// Copyright (c) 2026 Denis Spasyuk. MIT License.
// Relative scaling of multiple datasets before merging (XSCALE / AIMLESS style).
//
// When several scans (e.g. different kappa settings of one crystal, or several
// crystals of the same compound) are combined into one merged dataset, each
// input file carries its own arbitrary scale and a slightly different overall
// temperature factor (beam attenuation, absorption, crystal decay). Merging
// the raw intensities adds that systematic disagreement on top of the random
// errors. Following the classical SCALEPACK / AIMLESS approach, every dataset
// is put onto the common scale of the others with a smooth two-parameter
// correction
//
//     I' = I * K * exp(-B / (2 d^2))        (K = scale, B = A^2)
//
// fitted by weighted least squares on ln(<I_all> / <I_d>) versus s = 1/d^2,
// where <I_all> is the inverse-variance weighted mean over ALL datasets (the
// AIMLESS / XSCALE iteration; fitting against the mean of the others only and
// updating simultaneously would diverge for two datasets). Only well-measured
// reflections that occur in at least two datasets enter the fit, which is
// iterated to convergence. The anchor is the observation-weighted mean over
// the datasets, so the merged scale is preserved; per-dataset K and B are
// reported together with the inter-dataset agreement R (each dataset versus
// the mean of the OTHERS) before and after scaling.
//
// A full multi-scan absorption correction (--abs, spherical harmonics) is the
// better model when the scans differ in crystal orientation; this module
// handles the plain relative scale + decay between datasets.

import { canonicalRep } from './op-math.js';
import { reciprocalCell, dSpacing } from './merge.js';

/**
 * Compare each dataset's unit cell against a reference cell.
 * @param {object} refCell {a,b,c,alpha,beta,gamma}
 * @param {Array<object|null>} cells one cell per dataset (null when unknown)
 * @param {object} [options] { tolLen = 0.01 (1%), tolAng = 1.0 (deg) }
 * @returns {{ rows: Array<object|null>, warnings: string[], ok: boolean }}
 */
export function compareDatasetCells(refCell, cells, options = {}) {
    const tolLen = options.tolLen != null ? options.tolLen : 0.01;
    const tolAng = options.tolAng != null ? options.tolAng : 1.0;
    const rows = [];
    const warnings = [];
    for (let i = 0; i < cells.length; i++) {
        const c = cells[i];
        if (!c || !refCell || !Number.isFinite(c.a)) { rows.push(null); continue; }
        const lenPct = {
            a: 100 * Math.abs(c.a - refCell.a) / refCell.a,
            b: 100 * Math.abs(c.b - refCell.b) / refCell.b,
            c: 100 * Math.abs(c.c - refCell.c) / refCell.c,
        };
        const angDeg = {
            alpha: Math.abs(c.alpha - refCell.alpha),
            beta: Math.abs(c.beta - refCell.beta),
            gamma: Math.abs(c.gamma - refCell.gamma),
        };
        const maxLenPct = Math.max(lenPct.a, lenPct.b, lenPct.c);
        const maxAngDeg = Math.max(angDeg.alpha, angDeg.beta, angDeg.gamma);
        const bad = maxLenPct > 100 * tolLen || maxAngDeg > tolAng;
        rows.push({ lenPct, angDeg, maxLenPct, maxAngDeg, bad });
        if (bad) {
            warnings.push(
                `dataset ${i + 1}: cell deviates from the reference by ` +
                `${maxLenPct.toFixed(2)}% in length / ${maxAngDeg.toFixed(2)} deg in angle ` +
                `(tolerance ${100 * tolLen}% / ${tolAng} deg) - check the indexing`);
        }
    }
    return { rows, warnings, ok: warnings.length === 0 };
}

// Number of distinct datasets in a reflection list (by the .dataset tag).
export function datasetCount(reflections) {
    let n = 0;
    for (const r of reflections) {
        const d = Number.isFinite(r.dataset) ? r.dataset : 0;
        if (d + 1 > n) n = d + 1;
    }
    return Math.max(1, n);
}

/**
 * Fit per-dataset relative scales and merge-ready corrected reflections.
 * @param {Array<object>} reflections non-rejected observations with .dataset tags
 * @param {Array<object>} matrices Laue operators ({M, t}) used to group equivalents
 * @param {object} cell reference unit cell
 * @param {object} [options]
 *   cycles {number}      refinement cycles (default 5; 0 = statistics only)
 *   withB {boolean}      refine a B factor per dataset (default true)
 *   maxRelSigma {number} relative-sigma cutoff for fit/R statistics (default 0.25)
 * @returns {{
 *   ok: boolean, reason?: string,
 *   nDatasets: number, cyclesRun: number, converged: boolean, withB: boolean,
 *   datasets: Array<{index, nObs, nCommon, nFit, K, B, rBefore, rAfter, dmin, dmax}>,
 *   reflections: Array<object>   // corrected copies (identity when cycles = 0)
 * }}
 */
export function scaleDatasets(reflections, matrices, cell, options = {}) {
    const cycles = options.cycles != null ? options.cycles : 5;
    const withB = options.withB !== false;
    const maxRelSigma = options.maxRelSigma != null ? options.maxRelSigma : 0.25;
    const LN_K_CLAMP = Math.log(5);   // K in [1/5, 5]
    const B_CLAMP = 10;               // |B| <= 10 A^2

    const recip = reciprocalCell(cell);
    const nDs = datasetCount(reflections);
    if (nDs < 2) return { ok: false, reason: 'fewer than two datasets' };
    if (!recip) return { ok: false, reason: 'bad unit cell' };

    // Pre-compute s = 1/d^2 and the dataset tag per usable reflection.
    const refl = [];
    const sArr = [];
    const dsArr = [];
    const passthrough = [];   // reflections without a finite d (kept unmodified)
    for (const r of reflections) {
        const ds = Number.isFinite(r.dataset) ? r.dataset : 0;
        const d = dSpacing(r.h, r.k, r.l, cell, recip);
        if (!Number.isFinite(d) || !(d > 0) || !Number.isFinite(r.I)) {
            passthrough.push({ h: r.h, k: r.k, l: r.l, I: r.I, sig: r.sig, psi: r.psi, dataset: ds });
            continue;
        }
        refl.push(r);
        sArr.push(1 / (d * d));
        dsArr.push(ds);
    }
    if (!refl.length) return { ok: false, reason: 'no usable reflections' };

    // Group reflections by Laue-equivalent representative (done once; the
    // grouping does not depend on the scales).
    const groups = [];
    const keyMap = new Map();
    for (let i = 0; i < refl.length; i++) {
        const r = refl[i];
        const { rep } = canonicalRep([r.h, r.k, r.l], matrices);
        const key = rep[0] + ',' + rep[1] + ',' + rep[2];
        let gi = keyMap.get(key);
        if (gi === undefined) { gi = groups.length; keyMap.set(key, gi); groups.push([]); }
        groups[gi].push(i);
    }

    // Per-dataset bookkeeping.
    const nObs = new Int32Array(nDs);
    const dmin = new Float64Array(nDs).fill(Infinity);
    const dmax = new Float64Array(nDs).fill(0);
    for (let i = 0; i < refl.length; i++) {
        const ds = dsArr[i];
        nObs[ds]++;
        const d = 1 / Math.sqrt(sArr[i]);
        if (d < dmin[ds]) dmin[ds] = d;
        if (d > dmax[ds]) dmax[ds] = d;
    }

    // Correction parameters: corr(ds, s) = exp(lnK[ds] - B[ds]*s/2).
    const lnK = new Float64Array(nDs);
    const B = new Float64Array(nDs);

    // Per-group per-dataset means of the CURRENTLY corrected intensities.
    // acc[g][ds] = { n, sumI, sumSig2, hasSig } or null.
    function groupMeans() {
        const acc = new Array(groups.length);
        for (let g = 0; g < groups.length; g++) {
            const idxs = groups[g];
            let per = null;
            for (let j = 0; j < idxs.length; j++) {
                const i = idxs[j];
                const ds = dsArr[i];
                if (!per) { per = new Array(nDs).fill(null); }
                let a = per[ds];
                if (!a) { a = { n: 0, sumI: 0, sumSig2: 0, hasSig: true }; per[ds] = a; }
                const r = refl[i];
                const c = Math.exp(lnK[ds] - B[ds] * sArr[i] / 2);
                a.n++;
                a.sumI += r.I * c;
                const sig = Number.isFinite(r.sig) ? Math.abs(r.sig) * c : 0;
                if (sig > 0) a.sumSig2 += sig * sig; else a.hasSig = false;
            }
            acc[g] = per;
        }
        return acc;
    }

    // Mean and sigma-of-the-mean for one dataset slot.
    function meanOf(a) {
        const m = a.sumI / a.n;
        const sig = a.hasSig ? Math.sqrt(a.sumSig2) / a.n : 0;
        return { m, sig };
    }

    // Reference for the FIT of dataset ds in group g: inverse-variance weighted
    // mean over ALL datasets present (including ds itself). Scaling each
    // dataset onto the common mean is the AIMLESS / XSCALE iteration; it
    // converges in a couple of cycles. (Fitting against the mean of the OTHERS
    // and updating all datasets simultaneously diverges for two datasets: each
    // chases the other and the corrections double every cycle.)
    function allMean(per) {
        let wsum = 0, w = 0, ssum = 0, n = 0;
        for (let e = 0; e < per.length; e++) {
            if (!per[e]) continue;
            const { m, sig } = meanOf(per[e]);
            n++;
            if (sig > 0) { const wi = 1 / (sig * sig); wsum += m * wi; w += wi; }
            ssum += m;
        }
        if (!n) return null;
        const m = w > 0 ? wsum / w : ssum / n;
        const sig = w > 0 ? 1 / Math.sqrt(w) : 0;
        return { m, sig, n };
    }

    // Reference for the inter-dataset AGREEMENT statistics: weighted mean of
    // the OTHER datasets only (an honest cross-comparison). Null when the group
    // has no other dataset.
    function otherMean(per, ds) {
        let wsum = 0, w = 0, ssum = 0, n = 0;
        for (let e = 0; e < per.length; e++) {
            if (e === ds || !per[e]) continue;
            const { m, sig } = meanOf(per[e]);
            n++;
            if (sig > 0) { const wi = 1 / (sig * sig); wsum += m * wi; w += wi; }
            ssum += m;
        }
        if (!n) return null;
        const m = w > 0 ? wsum / w : ssum / n;
        const sig = w > 0 ? 1 / Math.sqrt(w) : 0;
        return { m, sig, n };
    }

    // Fit residual corrections of dataset ds against the overall mean:
    // ln(m_all / m_d) = dlnK - (dB/2) * s. Weighted 2-parameter linear LSQ
    // on well-measured reflections (both means positive, rel sigma <= cutoff).
    // Groups where ds is the only dataset give y = 0 and drop out.
    function fitDataset(acc, ds) {
        let w = 0, wx = 0, wx2 = 0, wy = 0, wxy = 0;
        let nFit = 0;
        for (let g = 0; g < acc.length; g++) {
            const per = acc[g];
            if (!per || !per[ds]) continue;
            const ref = allMean(per);
            if (!ref || !ref.n || ref.n < 2) continue;
            const { m, sig } = meanOf(per[ds]);
            if (!(m > 0) || !(ref.m > 0) || !(sig > 0) || !(ref.sig > 0)) continue;
            const rel = Math.sqrt((sig / m) * (sig / m) + (ref.sig / ref.m) * (ref.sig / ref.m));
            if (!(rel > 0) || rel > maxRelSigma) continue;
            // s of this group (all members share hkl equivalence -> same d).
            const s = sArr[groups[g][0]];
            const wi = 1 / (rel * rel);
            const x = withB ? s : 0;
            const y = Math.log(ref.m / m);
            w += wi; wx += wi * x; wx2 += wi * x * x; wy += wi * y; wxy += wi * x * y;
            nFit++;
        }
        if (!w) return { dlnK: 0, dB: 0, nFit: 0 };
        if (!withB) return { dlnK: wy / w, dB: 0, nFit };
        const det = w * wx2 - wx * wx;
        if (!(Math.abs(det) > 1e-12 * w * wx2)) return { dlnK: wy / w, dB: 0, nFit };
        const c0 = (wx2 * wy - wx * wxy) / det;
        const c1 = (w * wxy - wx * wy) / det;
        return { dlnK: c0, dB: -2 * c1, nFit };
    }

    // Per-dataset R against the mean of the others (strong common groups):
    // R = sum |m_d * corr - m_ref| / sum |m_ref|. `identity` skips the corr.
    function rVsOthers(acc, ds, identity) {
        let num = 0, den = 0, n = 0;
        for (let g = 0; g < acc.length; g++) {
            const per = acc[g];
            if (!per || !per[ds]) continue;
            const ref = otherMean(per, ds);
            if (!ref) continue;
            const { m, sig } = meanOf(per[ds]);
            if (!(m > 0) || !(ref.m > 0) || !(sig > 0) || !(ref.sig > 0)) continue;
            const rel = Math.sqrt((sig / m) * (sig / m) + (ref.sig / ref.m) * (ref.sig / ref.m));
            if (!(rel > 0) || rel > maxRelSigma) continue;
            const s = sArr[groups[g][0]];
            const c = identity ? 1 : Math.exp(lnK[ds] - B[ds] * s / 2);
            num += Math.abs(m * c - ref.m);
            den += Math.abs(ref.m);
            n++;
        }
        return { r: den > 0 ? num / den : null, n };
    }

    // "Before" statistics on the unscaled data.
    const accRaw = groupMeans();   // lnK = B = 0 here
    const before = [];
    const nCommon = new Int32Array(nDs);
    for (let ds = 0; ds < nDs; ds++) {
        const { r, n } = rVsOthers(accRaw, ds, true);
        before.push(r);
        nCommon[ds] = n;
    }

    // Iterative refinement of the per-dataset corrections. Only the RELATIVE
    // scales are determined by the data: the common mode of the updates drifts
    // (the fit's quality filter and the inverse-variance weighting of the
    // common mean do not select identical reflection sets), so the weighted
    // mean of each cycle's updates is removed before applying them. The
    // convergence check then measures the relative component only.
    let cyclesRun = 0;
    let converged = cycles === 0;
    const nFitLast = new Int32Array(nDs);
    let totN = 0;
    for (let ds = 0; ds < nDs; ds++) totN += nObs[ds];
    while (cyclesRun < cycles) {
        const acc = groupMeans();
        const dK = new Float64Array(nDs);
        const dBv = new Float64Array(nDs);
        let mdK = 0, mdB = 0;
        for (let ds = 0; ds < nDs; ds++) {
            const f = fitDataset(acc, ds);
            dK[ds] = f.dlnK;
            dBv[ds] = f.dB;
            nFitLast[ds] = f.nFit;
            mdK += nObs[ds] * f.dlnK;
            mdB += nObs[ds] * f.dB;
        }
        if (totN) { mdK /= totN; mdB /= totN; }
        let maxDK = 0, maxDB = 0;
        for (let ds = 0; ds < nDs; ds++) {
            lnK[ds] += dK[ds] - mdK;
            B[ds] += dBv[ds] - mdB;
            // Sanity clamps: a dataset needing more than this is not mergeable.
            if (lnK[ds] > LN_K_CLAMP) lnK[ds] = LN_K_CLAMP;
            if (lnK[ds] < -LN_K_CLAMP) lnK[ds] = -LN_K_CLAMP;
            if (B[ds] > B_CLAMP) B[ds] = B_CLAMP;
            if (B[ds] < -B_CLAMP) B[ds] = -B_CLAMP;
            if (Math.abs(dK[ds] - mdK) > maxDK) maxDK = Math.abs(dK[ds] - mdK);
            if (Math.abs(dBv[ds] - mdB) > maxDB) maxDB = Math.abs(dBv[ds] - mdB);
        }
        cyclesRun++;
        if (options.trace) options.trace(cyclesRun, maxDK, maxDB);
        if (maxDK < 1e-6 && maxDB < 1e-4) { converged = true; break; }
    }

    // Final anchor on the observation-weighted mean so the merged scale is
    // preserved (the per-cycle de-meaning keeps this a no-op up to rounding).
    let sumLnK = 0, sumB = 0;
    for (let ds = 0; ds < nDs; ds++) { sumLnK += nObs[ds] * lnK[ds]; sumB += nObs[ds] * B[ds]; }
    const lnK0 = totN ? sumLnK / totN : 0;
    const B0 = totN ? sumB / totN : 0;
    for (let ds = 0; ds < nDs; ds++) { lnK[ds] -= lnK0; B[ds] -= B0; }

    // "After" statistics with the anchored corrections.
    const accFin = groupMeans();
    const datasets = [];
    for (let ds = 0; ds < nDs; ds++) {
        const { r, n } = rVsOthers(accFin, ds, false);
        datasets.push({
            index: ds,
            nObs: nObs[ds],
            nCommon: nCommon[ds],
            nFit: nFitLast[ds],
            K: Math.exp(lnK[ds]),
            B: B[ds],
            rBefore: before[ds],
            rAfter: r,
            dmin: Number.isFinite(dmin[ds]) ? dmin[ds] : null,
            dmax: Number.isFinite(dmax[ds]) ? dmax[ds] : null,
        });
    }

    // Corrected reflections (new objects; the originals stay untouched).
    const out = [];
    for (let i = 0; i < refl.length; i++) {
        const r = refl[i];
        const ds = dsArr[i];
        const c = Math.exp(lnK[ds] - B[ds] * sArr[i] / 2);
        const o = { h: r.h, k: r.k, l: r.l, I: r.I * c, dataset: ds };
        if (Number.isFinite(r.sig)) o.sig = Math.abs(r.sig) * c;
        if (Number.isFinite(r.psi)) o.psi = r.psi;
        out.push(o);
    }
    for (const r of passthrough) out.push(r);

    return {
        ok: true,
        nDatasets: nDs,
        cyclesRun,
        converged,
        withB,
        maxRelSigma,
        datasets,
        reflections: out,
    };
}
