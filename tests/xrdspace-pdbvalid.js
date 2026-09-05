#!/usr/bin/env node
// Copyright (c) 2026 Denis Spasyuk. MIT License.
// Unit tests for the offline PDB space-group lookup (--valid).
// Offline tests — no network requests are made. The lookup table is built
// from hand-written rows with buildPdbLookupTable.

import { buildPdbLookupTable, searchPdbLookup, validateSpaceGroupAgainstPdb } from '../src/pdb-lookup.js';
import { niggliReduce } from '../src/cell-search.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
    if (cond) { pass++; console.log(`  ok   ${name}`); }
    else { fail++; console.error(`  FAIL ${name}${detail ? ' — ' + detail : ''}`); }
}
const close = (a, b, tol = 1e-3) => Math.abs(a - b) <= tol;

console.log('xrdspace PDB lookup unit tests');

// --- buildPdbLookupTable: sorting, rounding, filtering ---

console.log('buildPdbLookupTable');
const table = buildPdbLookupTable([
    // Monoclinic protein-like cells (a ~ 40-80 A).
    { id: '4hhb', cell: { a: 63.15, b: 83.59, c: 53.8, alpha: 90, beta: 99.34, gamma: 90 }, sg: 4, hm: 'P 1 21 1' },
    { id: '1a3n', cell: { a: 80.1, b: 80.1, c: 38.0, alpha: 90, beta: 90, gamma: 90 }, sg: 96, hm: 'P 43 21 2' },
    { id: '2mev', cell: { a: 31.9, b: 31.9, c: 31.9, alpha: 90, beta: 90, gamma: 90 }, sg: 205, hm: 'P A -3' },
    { id: '3bep', cell: { a: 40.203, b: 70.08, c: 73.872, alpha: 113.284, beta: 92.073, gamma: 99.361 }, sg: 1, hm: 'P 1' },
    { id: '9zzw', cell: { a: 4.7606, b: 4.7606, c: 12.994, alpha: 90, beta: 90, gamma: 120 }, sg: 176, hm: 'P 63/m' },
    // Bad rows: no cell, bad sg number, non-finite cell -> all dropped.
    { id: 'skip1', cell: null, sg: 1, hm: 'P 1' },
    { id: 'skip2', cell: { a: 10, b: 10, c: 10, alpha: 90, beta: 90, gamma: 90 }, sg: 0, hm: '?' },
    { id: 'skip3', cell: { a: NaN, b: 10, c: 10, alpha: 90, beta: 90, gamma: 90 }, sg: 2, hm: '?' },
]);
check('table keeps 5 valid entries', table.length === 5, `got ${table.length}`);
check('table sorted by reduced a-length', table.every((e, i) => i === 0 || table[i - 1].red[0] <= e.red[0]));
check('reduced cell is sorted a<=b<=c', table.every(e => e.red[0] <= e.red[1] && e.red[1] <= e.red[2]));
{
    const hex = table.find(e => e.id === '9zzw');
    check('hexagonal red cell preserved', close(hex.red[0], 4.7606) && close(hex.red[5], 120));
    const mono = table.find(e => e.id === '4hhb');
    check('monoclinic red cell recorded', close(mono.red[0], 53.8) && close(mono.red[2], 83.59) && close(mono.red[5], 99.34));
}

// --- searchPdbLookup: matches and non-matches ---

console.log('searchPdbLookup');
{
    // Same monoclinic cell, different setting (unique axis c instead of b).
    const s = searchPdbLookup(table, { a: 53.8, b: 63.15, c: 83.59, alpha: 90, beta: 90, gamma: 99.34 });
    check('permuted monoclinic setting found', s.some(m => m.id === '4hhb'), `ids ${s.map(m => m.id)}`);
}
{
    const s = searchPdbLookup(table, { a: 80.1, b: 80.1, c: 38.0, alpha: 90, beta: 90, gamma: 90 });
    check('tetragonal cell found', s.some(m => m.id === '1a3n'));
}
{
    // Slightly different cell -> still within 1% tolerance.
    const s = searchPdbLookup(table, { a: 79.7, b: 79.9, c: 37.9, alpha: 90, beta: 90, gamma: 90 });
    check('small cell drift still matches', s.some(m => m.id === '1a3n'));
}
{
    // Unrelated cell -> nothing.
    const s = searchPdbLookup(table, { a: 5, b: 5, c: 5, alpha: 90, beta: 90, gamma: 90 });
    check('unrelated cell finds nothing', s.length === 0, `got ${s.length}`);
}
{
    // Tight tolerance rejects a drifted cell.
    const s = searchPdbLookup(table, { a: 79.7, b: 79.9, c: 37.9, alpha: 90, beta: 90, gamma: 90 }, { tolLen: 0.002, tolAng: 0.5 });
    check('tight tolerance rejects drift', s.length === 0, `got ${s.length}`);
}

// --- validateSpaceGroupAgainstPdb ---

console.log('validateSpaceGroupAgainstPdb');
{
    const v = validateSpaceGroupAgainstPdb(table, { a: 63.15, b: 83.59, c: 53.8, alpha: 90, beta: 99.34, gamma: 90 }, 4);
    check('correct SG -> verified', v.verdict === 'verified', `verdict ${v.verdict}`);
    check('sgNumbers contains 4', v.sgNumbers.includes(4));
}
{
    const v = validateSpaceGroupAgainstPdb(table, { a: 63.15, b: 83.59, c: 53.8, alpha: 90, beta: 99.34, gamma: 90 }, 14);
    check('wrong SG -> mismatch', v.verdict === 'mismatch', `verdict ${v.verdict}`);
}
{
    // Enantiomorph pair: 96 (P 43 21 2) reported, but 92 (P 41 21 2) determined.
    const v = validateSpaceGroupAgainstPdb(table, { a: 80.1, b: 80.1, c: 38.0, alpha: 90, beta: 90, gamma: 90 }, 92);
    check('enantiomorph -> ambiguous', v.verdict === 'enantiomorph', `verdict ${v.verdict}`);
    check('enantiomorphOf is 96', v.enantiomorphOf === 96, `got ${v.enantiomorphOf}`);
}
{
    const v = validateSpaceGroupAgainstPdb(table, { a: 5, b: 5, c: 5, alpha: 90, beta: 90, gamma: 90 }, 1);
    check('no match -> none', v.verdict === 'none', `verdict ${v.verdict}`);
}

// --- raw entries array also accepted (mirrors loadPdbLookup table shape) ---

console.log('table shape compatibility');
{
    const wrapped = { entries: table, meta: { count: table.length } };
    const s = searchPdbLookup(wrapped, { a: 40.203, b: 70.08, c: 73.872, alpha: 113.284, beta: 92.073, gamma: 99.361 });
    check('wrapped {entries} object works', s.some(m => m.id === '3bep'));
}
{
    // The reduced cell stored must equal a fresh reduction of the original cell.
    const e = table.find(m => m.id === '3bep');
    const red = niggliReduce({ a: 40.203, b: 70.08, c: 73.872, alpha: 113.284, beta: 92.073, gamma: 99.361 });
    check('stored reduced cell matches fresh reduction',
        Math.abs(e.red[0] - red.a) < 1e-3 && Math.abs(e.red[1] - red.b) < 1e-3 && Math.abs(e.red[2] - red.c) < 1e-3);
}

console.log('');
console.log(`passed ${pass}, failed ${fail}`);
process.exit(fail ? 1 : 0);
