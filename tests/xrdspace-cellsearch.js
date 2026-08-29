#!/usr/bin/env node
// Copyright (c) 2026 Denis Spasyuk. MIT License.
// Unit tests for the unit-cell database search machinery (COD + PDB).
// These are offline tests — no network requests are made.

import {
    niggliReduce,
    cellSettings,
    transformCell,
    cellSimilarity,
    cellToleranceWindows,
} from '../src/cell-search.js';
import { cellVolume } from '../src/analyze.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
    if (cond) { pass++; console.log(`  ok   ${name}`); }
    else { fail++; console.error(`  FAIL ${name}${detail ? ' — ' + detail : ''}`); }
}
const close = (a, b, tol = 1e-3) => Math.abs(a - b) <= tol;
const near = (name, got, want, tol = 1e-2) => {
    const ok = close(got, want, tol);
    if (ok) pass++;
    else fail++;
    console.log(`  ${ok ? 'ok' : 'FAIL'}   ${name}: got ${got}, want ${want}${ok ? '' : ` (tol ${tol})`}`);
};

console.log('xrdspace cell-search unit tests');

// --- Niggli reduction ---

console.log('Niggli reduction');
{
    const r = niggliReduce({ a: 5, b: 7, c: 9, alpha: 80, beta: 85, gamma: 75 });
    near('triclinic a', r.a, 5);
    near('triclinic b', r.b, 7);
    near('triclinic c', r.c, 9);
    near('triclinic alpha', r.alpha, 80);
    near('triclinic beta', r.beta, 85);
    near('triclinic gamma', r.gamma, 75);
}
{
    const r = niggliReduce({ a: 79.98, b: 79.98, c: 123.95, alpha: 90, beta: 90, gamma: 120 });
    near('hexagonal a', r.a, 79.98);
    near('hexagonal b', r.b, 79.98);
    near('hexagonal c', r.c, 123.95);
    near('hexagonal gamma', r.gamma, 120);
}
{
    const r = niggliReduce({ a: 6, b: 6, c: 6, alpha: 60, beta: 60, gamma: 60 });
    near('rhombohedral a', r.a, 6);
    near('rhombohedral alpha', r.alpha, 60);
}

// --- same lattice, different settings -> identical reduced cell ---

console.log('Setting-invariant reduction');
{
    const m1 = niggliReduce({ a: 10.2, b: 12.3, c: 14.5, alpha: 90, beta: 105.5, gamma: 90 });
    const m2 = niggliReduce({ a: 10.2, b: 12.3, c: 14.5, alpha: 90, beta: 74.5, gamma: 90 });
    check('monoclinic axis flip -> same reduced cell',
        close(m1.a, m2.a) && close(m1.b, m2.b) && close(m1.c, m2.c)
        && close(m1.alpha, m2.alpha) && close(m1.beta, m2.beta) && close(m1.gamma, m2.gamma));
}
{
    const m1 = niggliReduce({ a: 8, b: 11, c: 15, alpha: 90, beta: 90, gamma: 90 });
    const m2 = niggliReduce({ a: 15, b: 8, c: 11, alpha: 90, beta: 90, gamma: 90 });
    check('orthorhombic axis permutation -> same reduced cell',
        close(m1.a, m2.a) && close(m1.b, m2.b) && close(m1.c, m2.c));
}

// --- cellSimilarity ---

console.log('cellSimilarity');
{
    const s = cellSimilarity(
        { a: 79.98, b: 79.98, c: 123.95, alpha: 90, beta: 90, gamma: 120 },
        { a: 80.466, b: 80.466, c: 124.929, alpha: 90, beta: 90, gamma: 120 });
    check('PDB 3UQH-like cell scores ~99%', s.match > 98, `match ${s.match}`);
    check('same lattice (setting permuted) scores 100%',
        cellSimilarity(
            { a: 8, b: 11, c: 15, alpha: 90, beta: 90, gamma: 90 },
            { a: 15, b: 8, c: 11, alpha: 90, beta: 90, gamma: 90 }).match > 99.9);
    const unrelated = cellSimilarity(
        { a: 5, b: 5, c: 5, alpha: 90, beta: 90, gamma: 90 },
        { a: 50, b: 50, c: 50, alpha: 90, beta: 90, gamma: 90 });
    check('unrelated cells score low', unrelated.match < 10, `match ${unrelated.match}`);
    const identical = cellSimilarity(
        { a: 4.7606, b: 4.7606, c: 12.994, alpha: 90, beta: 90, gamma: 120 },
        { a: 4.7606, b: 4.7606, c: 12.994, alpha: 90, beta: 90, gamma: 120 });
    check('identical cells score 100%', identical.match === 100);
}

// --- settings enumeration ---

console.log('cellSettings');
{
    const hex = cellSettings({ a: 79.98, b: 79.98, c: 123.95, alpha: 90, beta: 90, gamma: 120 });
    check('hexagonal has 6 settings', hex.length === 6, `got ${hex.length}`);
    const mono = cellSettings({ a: 10.86, b: 8.7, c: 7.76, alpha: 90, beta: 102.9, gamma: 90 });
    check('monoclinic has 12 settings', mono.length === 12, `got ${mono.length}`);
    const ortho = cellSettings({ a: 8, b: 11, c: 15, alpha: 90, beta: 90, gamma: 90 });
    check('orthorhombic has 6 settings', ortho.length === 6, `got ${ortho.length}`);
    const cubic = cellSettings({ a: 10, b: 10, c: 10, alpha: 90, beta: 90, gamma: 90 });
    check('cubic has 1 setting', cubic.length === 1, `got ${cubic.length}`);
    const tri = cellSettings({ a: 5, b: 7, c: 9, alpha: 80, beta: 85, gamma: 75 });
    check('triclinic has 24 settings', tri.length === 24, `got ${tri.length}`);
}

// All settings of a cell must have the same volume (they describe the same lattice).
{
    const vol0 = cellVolume({ a: 10.86, b: 8.7, c: 7.76, alpha: 90, beta: 102.9, gamma: 90 });
    let maxDev = 0;
    for (const s of cellSettings({ a: 10.86, b: 8.7, c: 7.76, alpha: 90, beta: 102.9, gamma: 90 })) {
        maxDev = Math.max(maxDev, Math.abs(cellVolume(s) - vol0) / vol0);
    }
    check('settings preserve cell volume', maxDev < 1e-6, `max rel. dev ${maxDev}`);
}

// --- transformCell ---

console.log('transformCell');
{
    // Swapping two axes is an involution.
    const cell = { a: 8, b: 11, c: 15, alpha: 90, beta: 90, gamma: 90 };
    const swap = [[0, 1, 0], [1, 0, 0], [0, 0, 1]];
    const t1 = transformCell(cell, swap);
    const t2 = transformCell(t1, swap);
    check('double axis swap returns to the original cell',
        close(t2.a, 8) && close(t2.b, 11) && close(t2.c, 15));
}

// --- cellToleranceWindows ---

console.log('cellToleranceWindows');
{
    const w = cellToleranceWindows([{ a: 10, b: 20, c: 30, alpha: 90, beta: 105, gamma: 90 }], { tolLen: 0.02, tolAng: 2 });
    check('window aMin/aMax', close(w[0].aMin, 9.8) && close(w[0].aMax, 10.2));
    check('window betaMin/betaMax', close(w[0].betaMin, 103) && close(w[0].betaMax, 107));
    check('window gammaMin/gammaMax', close(w[0].gammaMin, 88) && close(w[0].gammaMax, 92));
}

console.log('');
console.log(`passed ${pass}, failed ${fail}`);
process.exit(fail ? 1 : 0);
