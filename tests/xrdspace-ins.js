#!/usr/bin/env node
// Copyright (c) 2026 Denis Spasyuk. MIT License.
// Regression tests for the SHELX .ins generator (LATT / SYMM).
//
// Guards against three historical bugs:
//   1. scrambled LATT lattice-type codes (1=P,2=I,3=R,4=F,5=A,6=B,7=C),
//   2. inverted LATT sign (negative = NON-centrosymmetric),
//   3. SYMM lines that carried lattice-centering translations (SHELXT then
//      reports "Cannot identify Laue group from SYMM instructions").
//
// The strongest check generates an .ins for every space-group setting, reads it
// back with the SHELX semantics (LATT positive adds inversion, LATT magnitude
// selects the centering), and verifies the reconstructed general positions
// match the dictionary exactly.

import { loadSpaceGroups, writeShelxIns, centeringOf } from '../src/index.js';
import { isCentrosymmetric } from '../src/analyze.js';
import { opsFromLattSymm } from '../src/sg-model.js';
import { parseOperation, LATT_CENTERING } from '../src/op-math.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
    if (cond) { pass++; console.log(`  ok   ${name}`); }
    else { fail++; console.error(`  FAIL ${name}${detail ? ' — ' + detail : ''}`); }
}

// Rotation-only key: some centrosymmetric settings store the origin on a
// different inversion centre, so the same group may appear with shifted
// translations while the rotation parts (and order) are unchanged.
const rotKey = (op) => op.R.map(r => r.join(',')).join('|');

const CELL = { a: 10, b: 11, c: 12, alpha: 90, beta: 90, gamma: 90 };
const sgData = loadSpaceGroups();

console.log('xrdspace .ins generator regression tests');
console.log(`  (${sgData.length} space-group settings)`);

// --- known LATT / SYMM values for representative groups ---
const known = [
    { id: 14, latt: 1, symm: 1 },   // P 21/c, centrosymmetric P
    { id: 2, latt: 1, symm: 0 },    // P -1
    { id: 1, latt: -1, symm: 0 },   // P 1, acentric
    { id: 4, latt: -1, symm: 1 },   // P 21, acentric
    { id: 15, latt: 7, symm: 1 },   // C 2/c, C-centred
    { id: 5, latt: -7, symm: 1 },   // C 2, acentric C-centred
    { id: 139, latt: 2, symm: 7 },  // I 4/m m m
    { id: 167, latt: 3, symm: 5 },  // R -3 c (hexagonal axes)
];
for (const k of known) {
    const sg = sgData.find(g => g.id === k.id);
    const txt = writeShelxIns(sg, CELL, { sfac: ['C', 'H'] });
    const latt = parseInt(txt.match(/^LATT\s+(-?\d+)/m)[1], 10);
    const symm = [...txt.matchAll(/^SYMM/gm)].length;
    check(`#${k.id} ${sg.hm}: LATT ${latt}, ${symm} SYMM`, latt === k.latt && symm === k.symm,
        `want LATT ${k.latt}, ${k.symm} SYMM`);
}

// --- every setting round-trips through LATT/SYMM ---
const failures = [];
for (const sg of sgData) {
    const txt = writeShelxIns(sg, CELL, { sfac: ['C', 'H'] });
    const latt = parseInt(txt.match(/^LATT\s+(-?\d+)/m)[1], 10);
    const symms = [...txt.matchAll(/^SYMM (.*)$/gm)].map(m => m[1]);

    const wantLatt = (isCentrosymmetric(sg) ? 1 : -1) * (LATT_CENTERING[centeringOf(sg)] || 1);
    if (latt !== wantLatt) {
        failures.push(`#${sg.id} ${sg.hm}: LATT ${latt} (want ${wantLatt})`);
        continue;
    }
    const rebuilt = opsFromLattSymm(latt, symms);
    const expected = (sg.s || []).map(parseOperation).filter(Boolean);
    const rebuiltRot = rebuilt.map(rotKey).sort();
    const expectedRot = expected.map(rotKey).sort();
    if (rebuilt.length !== expected.length || JSON.stringify(rebuiltRot) !== JSON.stringify(expectedRot)) {
        failures.push(`#${sg.id} ${sg.hm}: ${rebuilt.length} general positions (want ${expected.length})`);
    }
}
check('LATT sign+centering and SYMM reconstruct the full group for every setting',
    failures.length === 0, failures.slice(0, 5).join('; '));

console.log(`passed ${pass}, failed ${fail}`);
process.exit(fail === 0 ? 0 : 1);
