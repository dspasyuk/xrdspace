#!/usr/bin/env node
// Copyright (c) 2026 Denis Spasyuk. MIT License.
// xrdspace validation against real macromolecular (MX) XDS_ASCII.HKL data.
//
// Each dataset has a POINTLESS run (pointless.xml) whose <BestSolution
// Type="spacegroup"> gives the ground-truth space group. We run xrdspace on
// the XDS_ASCII.HKL and compare. These are large protein cells (V >> 64000 A^3)
// so the chiral (Sohncke) restriction is active by default.
//
// Usage:
//   node tests/xrdspace-mx.js            # all datasets in mx/hkl
//   node tests/xrdspace-mx.js <tag>      # a specific dataset tag

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

import { analyzeHkl, isSohncke, cellVolume, loadSpaceGroups } from '../src/index.js';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
// Data location: set XRDSPACE_MX_DIR to a directory containing hkl/ and xml/
// subdirectories (the files are large and not committed to the repo).
const MX_DIR = process.env.XRDSPACE_MX_DIR || path.join(__dirname, '..', 'mx');
const HKL_DIR = path.join(MX_DIR, 'hkl');
const XML_DIR = path.join(MX_DIR, 'xml');

// Parse the POINTLESS ground-truth space group number from pointless.xml.
function parsePointlessSg(xmlText) {
    const m = xmlText.match(/<BestSolution[^>]*Type="spacegroup"[^>]*>[\s\S]*?<SGnumber>\s*(\d+)\s*<\/SGnumber>/);
    if (!m) return null;
    return parseInt(m[1], 10);
}

function runOne(tag) {
    const hklPath = path.join(HKL_DIR, tag + '.HKL');
    const xmlPath = path.join(XML_DIR, tag + '.xml');
    if (!fs.existsSync(hklPath) || !fs.existsSync(xmlPath)) return { tag, error: 'missing files' };
    const text = fs.readFileSync(hklPath, 'utf8');
    const xml = fs.readFileSync(xmlPath, 'utf8');
    const expected = parsePointlessSg(xml);
    const result = analyzeHkl(text, {});
    if (!result.ok) return { tag, expected, error: result.error };
    const determined = result.best ? result.best.id : null;
    const near = result.candidates.some(c => c.id === expected && c.violations === 0);
    const vol = cellVolume(result.cell);
    // result.best is a compact {id,hm,hs}; resolve the full dictionary entry
    // (with the .s operation list) to test chirality.
    const sgData = loadSpaceGroups();
    const fullBest = result.best ? sgData.find(g => g.id === result.best.id) : null;
    return {
        tag, expected, determined, near,
        bestHm: result.best ? result.best.hm : null,
        chiral: result.summary.chiral,
        volume: Math.round(vol),
        bestIsSohncke: fullBest ? isSohncke(fullBest) : null,
    };
}

function main() {
    const only = process.argv[2] || null;
    const tags = fs.readdirSync(HKL_DIR).filter(f => f.endsWith('.HKL')).map(f => f.replace(/\.HKL$/, ''));
    const targets = only ? tags.filter(t => t === only) : tags;
    if (!targets.length) { console.error('No MX datasets found in ' + HKL_DIR); process.exit(1); }

    console.log(`xrdspace validation against MX data — ${targets.length} datasets`);
    console.log('='.repeat(70));
    let pass = 0, near = 0, fail = 0, skip = 0;
    const rows = [];
    for (const tag of targets) {
        const r = runOne(tag);
        let status;
        if (r.error) { status = 'SKIP'; skip++; }
        else if (r.determined === r.expected) { status = 'PASS'; pass++; }
        else if (r.near) { status = 'NEAR'; near++; }
        else { status = 'FAIL'; fail++; }
        const detail = r.error || (r.bestHm ? `${r.bestHm} (No. ${r.determined})` : 'none');
        const chiralMark = r.chiral ? 'chiral' : '      ';
        console.log(`  [${status}] ${tag}  expected ${r.expected}  got ${detail}  [${chiralMark}]`);
        rows.push({ ...r, status });
    }

    console.log('='.repeat(70));
    const assessed = pass + near + fail;
    const rate = assessed ? ((pass + near) / assessed * 100) : 0;
    console.log(`PASS: ${pass}  NEAR: ${near}  FAIL: ${fail}  SKIP: ${skip}  (of ${targets.length})`);
    console.log(`Correct (PASS) or among zero-violation candidates (NEAR): ${rate.toFixed(1)}% of ${assessed}`);

    const report = { total: targets.length, pass, near, fail, skip, rate: +rate.toFixed(1), datasets: rows };
    const reportPath = path.join(__dirname, 'xrdspace-mx-report.json');
    fs.writeFileSync(reportPath, JSON.stringify(report, null, 2));
    console.log(`Report: ${reportPath}`);
    process.exit(rate < 90 ? 1 : 0);
}

main();
