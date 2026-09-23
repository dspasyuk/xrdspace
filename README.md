# xrdspace

**xrdspace** is a JavaScript (Node.js) tool for **space-group determination and
reflection merging** of X-ray diffraction data — an XPREP/POINTLESS-style
analysis of HKL files. It reads reflection data, determines the crystal system,
Laue class, lattice centering and space group, merges the data under the
determined symmetry, and writes out files ready for structure solution with
SHELXD / SHELXT / SHELXS.

- Pure JavaScript (ES modules), **zero runtime dependencies** — runs on any
  Node.js 18+ installation.
- Works both as a **command-line tool** (POINTLESS-style arguments) and as a
  **library** (import `analyzeHkl` from your own code or server).
- Ships with a complete dictionary of all **230 space groups** (every setting,
  with full general-position symmetry operations).
- **Unit-cell database search**: search the **Crystallography Open Database
  (COD)** and the **RCSB Protein Data Bank (PDB)** for structures matching a
  query cell (`--codsearch / --pdbsearch / --search`), ranked by a
  Niggli-reduced-cell match score.
- Validated against **2000** and **10 000 real structures from the
  Crystallography Open Database (COD)**: the published space group is recovered
  exactly (PASS) or appears among the zero-violation candidates (NEAR) in
  **98.1 %** / **99.2 %** of assessed entries (exact match **81.0 %** on the
  10 000) — and against **real macromolecular (protein) data**, where every
  determined space group is chiral (Sohncke). On a head-to-head subset it
  recovers the published group for **99.0 %** of entries versus XPREP's 85.9 %
  exact.

---

## Features

| Step | What xrdspace does |
|---|---|
| 1. Parse | Reads **XDS_ASCII**, **SHELX five-column** and **COD `.hkl` (CIF)** files, extracting reflections, unit cell, wavelength, title |
| 2. Crystal system | From the unit-cell metric (length/angle tolerances), with automatic fallback when the data demands lower symmetry than the metric suggests (pseudo-symmetry) |
| 3. Laue class | R(sym) merge test over all 11 Laue classes (all settings of 2/m tried); the highest-symmetry metric-compatible class whose R(sym) is close to the intrinsic (−1) merge is chosen |
| 4. Centering | Bravais lattice (P/A/B/C/I/F/R) from reflection parity (systematic absences of the centering conditions), picking the most restrictive centering with no significant violations |
| 5. Space group | All candidates of the crystal system + centering are scored by their **systematic-absence conditions** (screw axes and glide planes, op by op, including absences that are simply missing from the data). Ranking: fewest violations → most confirmed absences → Laue-class match → Wilson centricity match. For macromolecular cells (volume > 64 000 Å³, ≈ 40×40×40) candidates are restricted to the **65 chiral (Sohncke) space groups** — see below |
| 6. Centricity | Wilson-style test on \|E²−1\| (centric ≈ 0.968, acentric ≈ 0.736) used as a tie-breaker |
| 7. Merge | Reflections merged under the chosen Laue class with 1/σ² weights; merged σ combines the weighted-mean error with the sample scatter |
| 8. Output | Merged **SHELX** HKL, merged **XDS_ASCII** HKL, **unmerged XDS_ASCII**, a **SHELX `.ins`** instruction file, and a consolidated `xrdspace.log` report: merging statistics (incl. resolution/completeness at I/σ = 1 and CC(1/2) = 0.30), a per-resolution-shell quality table, and artifact flags (outliers / ice rings / anisotropy) |
| 9. Cell search | Search the **Crystallography Open Database (COD)** and the **RCSB Protein Data Bank (PDB)** for structures whose unit cell matches a query cell. Matching is done in the **Niggli-reduced cell**, so different settings of the same lattice (axis permutations, unique-axis choices, obtuse/acute angle conventions) are recognised automatically and ranked by match score |
| 10. PDB validation | `--valid` checks the determined space group against an **offline PDB unit-cell/space-group lookup table** (`data/pdb-cells.json`, built once from RCSB). No network access at validation time: reports **VERIFIED / MISMATCH / AMBIGUOUS (enantiomorph) / INDETERMINATE**, with the space groups PDB assigns to matching cells |

---

## Requirements

- **Node.js 18 or newer** (uses built-in `fetch` in the test harness; the
  library itself only needs `node:fs`, `node:vm`, `node:path`, `node:url`).
- No `npm install` needed — there are no dependencies.

```sh
git clone https://github.com/dspasyuk/xrdspace.git
cd xrdspace
```

---

## Command-line usage

```
node src/xrdspace.js --hklin <file.hkl> [options]
node src/xrdspace.js --search --cell "a b c alpha beta gamma" [options]
```

Bare POINTLESS-style keywords (`hklin`, `hklout`, `spacegroup`, `cell`, …) are
accepted as well, and a bare filename as the first argument is treated as the
input file:

```sh
node src/xrdspace.js hklin data.hkl hklout merged.hkl spacegroup "P 21/c"
```

### Options

| Option | Description |
|---|---|
| `--hklin <file>` | Input HKL file (XDS_ASCII, SHELX five-column, or COD `.hkl`) |
| `--hklout <file>` | Output merged HKL file in **SHELX format** (default: `<input>_merged.hkl`) |
| `--xdsout <file>` | Output merged HKL file in **XDS_ASCII format** (default: `<input>_xds.hkl`) |
| `--unmergedout <file>` | Output **unmerged** HKL file in **XDS_ASCII format** (`MERGE=FALSE`), keeping every observation (default: `<input>_unmerged.hkl`) |
| `--log <file>` | Path for the consolidated, formatted report. The report is **always** written; default `xrdspace.log` next to the input file |
| `--spacegroup <sg>` | **Force** a specific space group — number (`14`) or Hermann–Mauguin symbol (`"P 21/c"`, `"P-1"`, `"P 21 21 21"`). Used for merging/output and checked for consistency with the data |
| `--laue <group>` | **Force** a Laue class for merging (e.g. `-1`, `2/m`, `mmm`, `4/mmm`, `-3m`, `6/mmm`, `m-3m`) |
| `--cell "a b c alpha beta gamma"` | Unit cell, used when the file does not carry cell parameters (skips the interactive prompt) |
| `--resolution "lo hi"` | Restrict the analysis to a resolution range in Å (low = large d, high = small d) |
| `--sigthreshold <n>` | I/σ(I) significance threshold for systematic-absence tests (default `5`) |
| `--sfac "C H N O"` | Expected elements — or a formula such as `"C12 H16 N2 O4"` — written into the SHELX `.ins` `SFAC`/`UNIT` lines for SHELXT |
| `--chiral` | Restrict candidates to the 65 **chiral (Sohncke)** space groups. This is the **default for macromolecular cells** (volume > 64 000 Å³, ≈ 40×40×40 Å) |
| `--no-chiral` | Allow non-chiral (centrosymmetric / mirror) space groups even for large cells |
| `--search` | **Search** both the Crystallography Open Database (COD) and the RCSB Protein Data Bank (PDB) for structures with a unit cell matching `--cell`. No HKL file is needed |
| `--codsearch` | Search only the COD |
| `--pdbsearch` | Search only the PDB |
| `--tol <pct>` | Relative length tolerance in **%** for the cell match (default `1.0`). COD reports esds on the cell parameters (shown in the results); PDB does not, so a few % is usually appropriate for protein data |
| `--tol-angle <deg>` | Angle tolerance in **degrees** (default `1.5`) |
| `--limit <n>` | Maximum number of matches to report (default `20`) |
| `--valid` | **Validate** the determined space group against an offline PDB unit-cell/space-group lookup table (default `data/pdb-cells.json`). No network access: reports VERIFIED / MISMATCH / AMBIGUOUS (enantiomorph) / INDETERMINATE, with the space groups the PDB assigns to cells matching the query cell (uses `--tol` / `--tol-angle`) |
| `--pdb-table <file>` | Path to the PDB lookup table used by `--valid` (default: `data/pdb-cells.json`). Build it once with `node scripts/build-pdb-table.js` |
| `--help`, `-h` | Show help |
| `--version`, `-v` | Show version |

If the input file has no unit-cell parameters and `--cell` is not given,
xrdspace **prompts interactively** for `a b c alpha beta gamma`.

### Example

```sh
node src/xrdspace.js --hklin data.hkl \
    --spacegroup "P 21/c" \
    --sfac "C12 H16 N2 O4" \
    --resolution 50 1.2
```

### Example output

```
==============================================================================
  xrdspace  —  space-group determination and reflection merging
==============================================================================

  Input file         : data.hkl
  Format             : xds_ascii
  Wavelength         : 0.71073 A
  Reflections        : 123456

------------------------------------------------------------------------------
  UNIT CELL
------------------------------------------------------------------------------
           a         b         c     alpha      beta     gamma
     10.500    10.500    14.000    90.000    90.000    90.000
  Volume             : 1544 A^3
  Crystal system     : tetragonal
  Lattice centering  : P

------------------------------------------------------------------------------
  SPACE-GROUP DETERMINATION
------------------------------------------------------------------------------
  Best space group   : P 42/m  (No. 84)
  Laue class         : 4/mmm   R(sym) = 1.82 %
  Centrosymmetric    : centric   (<|E^2-1|> = 0.971)
  Data consistency   : consistent with data

  R(sym) by Laue class:
    -1      order  2   R(sym) =   3.41 %
    2/m     order  4   R(sym) =   3.38 %
    mmm     order  8   R(sym) =   2.95 %
    4/m     order  8   R(sym) =   2.10 %
    4/mmm   order 16   R(sym) =   1.82 %  <-- chosen
    ...

  Space-group candidates (systematic absences):
    No.  HM                          violations
     88  I 41/a                             0
     84  P 42/m                             0  <-- best

------------------------------------------------------------------------------
  MERGING STATISTICS
------------------------------------------------------------------------------
  Resolution range    : 50.00 - 1.20 A
  Resolution (I/σ=1)  : 1.35 A
  Resolution (CC1/2)  : 1.30 A   (CC1/2 = 0.30)
  Observations        : 123456
  Unique reflections  : 7712
  Mean multiplicity   : 16.0
  Completeness        : 98.7 %
  Completeness (I/σ=1): 99.5 %
  Completeness (CC1/2): 99.0 %
  R(merge)            : 1.82 %
  R(meas)             : 2.35 %
  R(pim)              : 0.46 %
  Mean I/sigma(I)     : 14.8   (to CC1/2 = 0.30)

------------------------------------------------------------------------------
  RESOLUTION SHELLS
------------------------------------------------------------------------------
  d range (A)      #obs   #uniq   Compl   Mult   Rmerge    Rmeas    Rpim   <I/σ>  CC(1/2)   %neg
  50.00-3.50      32000    4200   99.8%    7.6     3.1%     3.6%     1.2%   22.1    0.999     1%
  3.50-2.50       38000    5200   99.5%    7.3     8.4%     9.8%     3.2%   11.5    0.995     3%
  2.50-2.00       28000    4100   98.9%    6.8    24.6%    28.5%     9.1%    5.2    0.982     8%
  2.00-1.60       18000    3000   97.2%    6.0    78.3%    90.6%    28.4%    2.3    0.921    18%
  1.60-1.20       12000    2500   94.1%    4.8   210.5%   243.1%    76.2%    1.1    0.702    33%
  ...

------------------------------------------------------------------------------
  QUALITY FLAGS / ARTIFACTS
------------------------------------------------------------------------------
  Outliers |ΔI|/σ>10  : 63 / 123456 (0.05%)
     h    k    l      d(A)         I    sigma   |ΔI|/σ     n
       4    2   -9     3.12       88.40     3.10     48.7      6
      ...
  Ice-ring scan (mean |I| in a 0.05 A band vs local background):
    none detected
    3.90 A   n=  1905   ratio 1.08
    3.67 A   n=  2349   ratio 1.16
    ...
  Anisotropy (resolution at I/σ=2 by reciprocal axis):
    a*   1.55 A   (n=41234)
    b*   1.62 A   (n=39001)
    c*   1.58 A   (n=43022)
    max/min = 1.05   (isotropic)

------------------------------------------------------------------------------
  OUTPUT FILES
------------------------------------------------------------------------------
  data_merged.hkl  (SHELX format, ready for SHELXD/SHELXT)
  data_xds.hkl  (merged XDS_ASCII)
  data_unmerged.hkl  (UNMERGED XDS_ASCII, all observations)
  data_merged.ins  (SHELX instructions, matching cell/space group)
  xrdspace.log  (this report)

==============================================================================
```

The same consolidated report is written to **`xrdspace.log`** (next to the input
file) on every analysis run; use `--log <file>` to choose another path.

### Unit-cell database search (COD + PDB)

Instead of analyzing an HKL file you can search the Crystallography Open
Database (COD) and the RCSB Protein Data Bank (PDB) for structures whose unit
cell matches a query cell — e.g. to identify an unknown phase from the cell,
or to find isostructural (cell-isomorphous) compounds:

```sh
node src/xrdspace.js --pdbsearch --cell "40.203 70.080 73.872 113.284 92.073 99.361"
node src/xrdspace.js --codsearch --cell "4.7606 4.7606 12.994 90 90 120"
node src/xrdspace.js --search    --cell "10.86 8.70 7.76 90 102.9 90" --tol 1.5 --tol-angle 2 --limit 10
```

How the matching works:

- The query cell is expanded into its **standard settings** (axis permutations
  and sign/unique-axis conventions), so a query reported in one setting finds
  entries deposited in another — e.g. a monoclinic cell queried with unique
  axis *b* finds structures published with the equivalent setting.
- Both databases are searched within the requested tolerance windows
  (`--tol` length %, `--tol-angle` degrees; defaults 1 % and 1.5°). PDB does
  not store esds on cell parameters, so a wider tolerance is usually needed for
  protein data; COD esds are shown in the results when present.
- Every hit is reduced to its **Niggli-reduced cell** (Krivý–Gruber 1976) and
  scored against the query: `match = 100 − (max relative length deviation in % +
  max angle deviation in degrees)`. 100 % means an identical lattice. Results
  are ranked by this score and the "best solution" is listed first.

Example output (`--pdbsearch` for the triclinic protein cell above):

```
==============================================
  xrdspace  —  unit-cell database search
==============================================
  Query cell         : 40.203 70.080 73.872  113.3 92.1 99.4
  Databases          : PDB
  Length tolerance   : 1%   angle tolerance: 1.5 deg
  Match score        : 100 - (max rel. length dev % + max angle dev deg)
----------------------------------------------
  Standard settings  : 24
----------------------------------------------
  Rank   Match  DB    ID      Space group      a       b       c    alpha  beta  gamma
  ------------------------------------------------------------------------------------------------
     1  100.0%  PDB  3BEP     P 1             40.20 70.08 73.87  113.3 92.1 99.4
       Structure of a sliding clamp on DNA
       space group #1
==============================================
```

### PDB space-group validation (`--valid`)

Space-group determination sometimes makes mistakes (pseudo-symmetry, weak data).
To cross-check a determination against the PDB **without querying RCSB on every
run**, xrdspace ships an offline lookup: a single JSON table with the unit cell
and space group of every PDB entry.

Build the table once (downloads cell + space-group data for the whole PDB via
the RCSB Search + GraphQL APIs — no mmCIF files, typically under a minute):

```sh
node scripts/build-pdb-table.js            # writes data/pdb-cells.json (~40 MB, git-ignored)
node scripts/build-pdb-table.js --limit 5000 --out mini.json   # partial table (testing)
```

Then validate the space group determined from an HKL file:

```sh
node src/xrdspace.js --hklin data_XDS.HKL --valid
node src/xrdspace.js --hklin data_XDS.HKL --valid --pdb-table /path/to/pdb-cells.json
```

The determined space group is checked against the PDB entries whose
**Niggli-reduced cell** matches the query cell within `--tol` (length %) and
`--tol-angle` (degrees) — different settings of the same lattice are
recognised automatically, exactly like the cell search. The result is one of:

- **VERIFIED** — PDB assigns the same space group to matching cells.
- **MISMATCH** — matching PDB cells belong to a different space group; the
  assignment deserves a second look.
- **AMBIGUOUS** — PDB assigns the *enantiomorph* (e.g. `P 41 21 2` vs
  `P 43 21 2`): the two cannot be told apart from the diffraction pattern
  alone (only by anomalous scattering).
- **INDETERMINATE** — no PDB structure matches the cell (widen the tolerances
  or rebuild the table after the PDB grows).

Example output:

```
==============================================
  xrdspace  —  PDB space-group validation
==============================================
  Lookup table      : data/pdb-cells.json (2289 entries)
  Query cell        : 63.150 83.590 53.800  90.0 99.3 90.0
  Tolerances        : 1% lengths, 1.5 deg angles
  Matching PDB      : 2 entries
  PDB space groups  : P 1 21 1 (No. 4) x2
  Matches           : 1BAB (No. 4), 1A01 (No. 4)
----------------------------------------------
  Determined SG     : P 1 21/c 1 (No. 14)
  Result            : MISMATCH
    PDB structures with this cell are in space group 4;
    No. 14 was determined. Check the indexing / space-group assignment.
==============================================
```

### Chiral (Sohncke) space groups for macromolecular data

Protein and other macromolecular crystals are almost always in one of the
**65 chiral (Sohncke) space groups** — groups without inversion, mirrors,
glides or roto-inversions. To avoid reporting a non-chiral group (e.g.
`P 21/c`) for a protein dataset, xrdspace automatically restricts the
candidates to the Sohncke groups when the **unit-cell volume exceeds
64 000 Å³** (a 40 × 40 × 40 Å cube):

- The restriction is **not a hard veto**: if no Sohncke group is consistent
  with the systematic absences (i.e. every chiral candidate has violations),
  xrdspace falls back to the full candidate list, so a genuinely
  non-centrosymmetric-free large cell is still handled.
- The restriction is reported in the output:
  `Chiral restriction : on (Sohncke space groups only)`.
- Override it explicitly with `--chiral` (force on) or `--no-chiral`
  (force off).

### Output files

| File | Format | Use |
|---|---|---|
| `<input>_merged.hkl` | SHELX five-column `H K L I SIG(I)` | Feed directly to **SHELXD / SHELXT / SHELXS** |
| `<input>_xds.hkl` | Merged XDS_ASCII (`MERGE=TRUE`, `FRIEDEL'S_LAW=TRUE`) with cell, space group, wavelength and resolution-range header | Re-integration / further processing |
| `<input>_unmerged.hkl` | **Unmerged** XDS_ASCII (`MERGE=FALSE`, labelled `UNMERGED (all observations)`) with cell, space group, wavelength and resolution-range header. XDS_ASCII input records are copied verbatim (auxiliary `XD`, `YD`, `ZD`, `RLP`, `PEAK`, `CORR`, `PSI` columns preserved); other formats are written as 12-column records | Programs that prefer redundancy — e.g. **Phenix** (`phenix.refine` / `phenix.xtriage`) and anomalous-data work. **MOLREP** usually takes the merged file |
| `<input>_merged.ins` | SHELX instruction file: `TITL`, `CELL`, `LATT` (sign encodes centrosymmetry), `SYMM` (generating operations, one per inversion pair for centric groups), `SFAC`, `UNIT`, `HKLF 4`, `TREF 50` | Structure solution with SHELXT |
| `xrdspace.log` | Consolidated aligned text report: input/cell, space-group determination, Laue table, candidates, merging statistics, optional PDB validation, and the output-file list (no duplicated blocks) | Archive the analysis; default name, change with `--log <file>` |

### Model transform between space groups (`sg-model.js`)

Once a structure exists, the space group determined from the HKL file may differ
from the one the model was solved/refined in. `src/sg-model.js` rewrites a SHELX
`.res/.ins` model into a target space group: the model's asymmetric unit is
expanded under its **current** symmetry (`LATT` + `SYMM`) to reconstruct the full
cell content, then reduced to an asymmetric unit under the **target** group's
operators. Raising the symmetry therefore removes redundant (symmetry-related)
molecules; lowering it adds the symmetry partners — coordinates stay in the same
fractional frame, scattering-factor types, occupancies and ADP columns are
preserved, and atoms are renumbered.

```js
import { transformModelToSpaceGroup, parseShelxModel, opsFromLattSymm } from './src/sg-model.js';

const r = transformModelToSpaceGroup(resText, 'P -1');   // or 2, or 'P 21/c'
// r = { ok, hm, sgId, nOldAsu, nFull, nNewAsu, removed, added, report, res }
```

Example: a P-1 asymmetric unit forced to P 1 grows from 2 to 4 atoms (the
inversion partners become independent); transforming back to P -1 removes them
again (`removed = 2`).

---

## Library API

```js
import { analyzeHkl, loadSpaceGroups, resolveSpaceGroup, verdict } from './src/index.js';
```

### `analyzeHkl(text, options)`

Runs the full analysis on HKL file **text** and returns a result object.

**Options**

| Option | Type | Description |
|---|---|---|
| `cell` | `{a,b,c,alpha,beta,gamma}` | Unit cell, required when the file has none |
| `spaceGroup` | `number \| string` | Force a space group (number or Hermann–Mauguin / Hall symbol) |
| `laue` | `string` | Force a Laue class for merging (e.g. `'2/m'`, `'mmm'`) |
| `resolution` | `{dmin, dmax}` | Restrict analysis to a resolution range (Å) |
| `sigThreshold` | `number` | I/σ threshold for systematic absences (default `5`) |
| `xdsOutput` | `string` | `OUTPUT_FILE` name written into the merged XDS_ASCII header |
| `sfac` | `string[]` | Element symbols for the `.ins` `SFAC` line |
| `unit` | `number[]` | Counts per element for the `.ins` `UNIT` line |
| `chiral` | `boolean` | Restrict candidates to the 65 chiral (Sohncke) space groups. Default: `true` for cells with volume > 64 000 Å³ (≈ 40×40×40), `false` otherwise |
| `quality` | `boolean` | Also compute the per-resolution-shell table and the artifact flags (default `false`; the CLI always enables it) |

**Return value**

```js
{
  ok: true,
  cell: { a, b, c, alpha, beta, gamma },
  summary: {
    format, title, wavelength, nReflections,
    crystalSystem, metricSystem, uniqueAxis,
    laueClass, laueRSym, centering,
    centricity,          // 'centric' | 'acentric' | 'indeterminate'
    centricityScore,     // <|E^2-1|>
    chiral,              // true when the Sohncke (chiral) restriction was applied
    forced,              // true when a space group was forced
    bestSpaceGroup, bestSpaceGroupNumber,
    merged: { nUnique, nObs, completeness, rMerge, rPim, meanIsig, meanMultiplicity }
  },
  laueTable:        [{ name, order, rsym, nOrbits, chosen }],
  centeringResults: { P: {...}, A: {...}, B: {...}, C: {...}, I: {...}, F: {...}, R: {...} },
  candidates:       [{ id, hm, hs, laue, centric, violations, confirmedOps, confirmedAbsences }],  // top 30
  best:             { id, hm, hs },   // space group used (forced or determined)
  determined:       { id, hm, hs },   // space group determined from the data
  forced:           { id, hm, hs } | null,
  merge: {
    nUnique, nObs,
    shelxHkl,            // merged SHELX five-column text
    xdsAscii,            // merged XDS_ASCII text
    unmergedXdsAscii,    // UNMERGED XDS_ASCII text (all observations, MERGE=FALSE)
    inputWasMerged,      // true when the input declared MERGE=TRUE
    shelxIns,            // SHELX .ins text
    statistics: { dmin, dmax, dIsig1, nUniqueIsig1, completenessIsig1,
                  dCC30, nUniqueCC30, completenessCC30,
                  nObs, nUnique, meanMultiplicity, completeness,
                  rMerge, rMeas, rPim, meanIsig, meanI },
    report,              // human-readable merging report
    consistency: { violations, confirmedOps, confirmedAbsences },
    shells,              // (quality) [{ dLo, dHi, nObs, nUnique, completeness,
                         //   multiplicity, rMerge, rMeas, rPim, meanIsig, ccHalf, negFrac }]
    artifacts,           // (quality) { outliers, iceRings, anisotropy, anisoRatio, anisotropic }
  }
}
```

On failure the result is `{ ok: false, error }` where `error` is a message or
the special code `'NO_CELL'` (the file has no unit cell — supply `options.cell`).

### Other exports

| Export | Description |
|---|---|
| `loadSpaceGroups()` | The full 230-space-group dictionary (all settings): `{id, hm, hs, o, s[]}` |
| `getLaueGroups()` | The 11 Laue classes with reciprocal-space operation matrices |
| `resolveSpaceGroup(sgData, spec)` | Resolve a space group by number or symbol |
| `writeShelxIns(sg, cell, options)` | Generate a SHELX `.ins` file for a given space group and cell |
| `verdict(result)` | One-line verdict string, e.g. `"P 21/c (No. 14)"` |
| `isSohncke(sg)` | `true` when a space group is chiral (no op with negative rotation determinant) |
| `cellVolume(cell)` | Unit-cell volume in Å³ |
| `niggliReduce(cell)` | Niggli-reduced cell (Krivý–Gruber 1976): `{a,b,c,alpha,beta,gamma,A,B,C,D,E,F}` |
| `cellSettings(cell)` | Distinct standard settings of a cell (axis permutations / sign conventions) |
| `transformCell(cell, M)` | Cell after an integer change of basis M |
| `cellSimilarity(cellA, cellB)` | Reduced-cell match: `{match, dLmaxPct, dAmaxDeg, reducedQuery, reducedCandidate}` |
| `cellToleranceWindows(settings, opts)` | Per-setting search windows `{aMin,…}` from a tolerance |
| `searchCodByCell(cell, opts)` | Search the COD by unit cell (returns ranked entries incl. esds) |
| `searchPdbByCell(cell, opts)` | Search the RCSB PDB by unit cell (returns ranked entries) |
| `searchByCell(cell, opts)` | Search both databases: `{settings, results, total, errors}` |
| `buildPdbLookupTable(rows)` | Build the offline PDB lookup entries from raw `{id, cell, sg, hm}` rows (computes Niggli-reduced cells, sorts by reduced a) |
| `loadPdbLookup(file)` | Load the PDB lookup table JSON (cached) |
| `searchPdbLookup(table, cell, opts)` | Entries of the lookup table whose Niggli-reduced cell matches `cell` within `tolLen`/`tolAng` |
| `validateSpaceGroupAgainstPdb(table, cell, sgId, opts)` | Validate SG number `sgId` against the table: `{verdict, matches, total, sgNumbers, sgCounts, ...}` with `verdict` = `verified \| mismatch \| enantiomorph \| none` |
| `transformModelToSpaceGroup(text, sg)` | Rewrite a SHELX `.res/.ins` model into space group `sg` (number/symbol/object): `{ok, hm, added, removed, report, res}` — adds symmetry partners when lowering symmetry, removes redundant molecules when raising it |
| `parseShelxModel(text)` | Parse a SHELX model: `{title, cell, latt, symm, sfac, atoms, ...}` |
| `opsFromLattSymm(latt, symm)` | Full general-position operator set (closed) from `LATT` + `SYMM` lines |

---

## Supported input formats

| Format | Detection | Notes |
|---|---|---|
| **XDS_ASCII** | `!` header lines | Cell from `!UNIT_CELL_CONSTANTS=`, plus `SPACE_GROUP_NUMBER/NAME`, `X-RAY_WAVELENGTH`, `MERGE`, `FRIEDELS_LAW` |
| **SHELX five-column** | 5+ numeric columns `H K L I SIG(I)` | No cell in the file — provide `--cell` |
| **COD `.hkl`** | CIF `loop_` with `_refln_` keys | Reads `F²_meas` (+σ), `I_meas` (+σ), `F_meas` (+σ) or `f_obs` (+σ); cell must be supplied |

---

## How it works

1. **Crystal system from the metric** — length equality (0.5 %) and angle
   (1°) tolerances classify the cell as cubic / tetragonal / hexagonal /
   trigonal (rhombohedral) / orthorhombic / monoclinic (unique axis detected) /
   triclinic.
2. **Laue class by R(sym)** — for each of the 11 Laue classes (all settings of
   2/m are tried) the reflections are merged under the class's reciprocal
   matrices and R(sym) is computed (strongest 30 000 reflections for very large
   data). The chosen class is the highest-symmetry one compatible with the
   metric whose R(sym) is within `max(7 %, 2 × R(sym) of −1)`; otherwise the
   crystal system is downgraded to match the data (pseudo-symmetry handling).
3. **Centering by parity** — each Bravais condition (e.g. C: h+k = 2n, F:
   h,k,l unmixed, R: −h+k+l = 3n) is checked; reflections above the I/σ
   threshold that are forbidden count as violations. The most restrictive
   centering with zero violations (and no significant weak-forbidden presence)
   wins.
4. **Systematic absences** — for every candidate space group (crystal system +
   centering, all settings of a given number scored and the best kept) each
   screw/glide operation with a non-lattice translation is tested: reflections
   invariant under the operation must have phase `h·t ≈ 0 mod 1`. Strong
   violations reject the group; confirmed absences (weak or entirely missing
   forbidden reflections along the invariant axes/planes) rank it. For
   macromolecular cells (volume > 64 000 Å³) the candidate pool is first
   restricted to the 65 Sohncke (chiral) groups — a group is chiral when none
   of its operations has a negative rotation determinant — with automatic
   fallback to the full pool if no chiral group is consistent with the data.
5. **Centricity** — the Wilson statistic <|E²−1|> (E² = I/⟨I⟩) distinguishes
   centric (≈ 0.968) from acentric (≈ 0.736) data and breaks remaining ties.
   Intensities are normalized per resolution shell (a Wilson correction) and
   negative measurements are ignored; the intermediate band (0.80–0.90) is
   reported as **indeterminate**.
6. **Merging** — reflections are grouped into Laue orbits (canonical
   representative = lexicographically smallest image), merged with 1/σ²
   weights; the merged σ adds the standard error of the mean so inconsistent
   observations inflate the error. Completeness is computed against the exact
   count of unique reflections in the resolution shell (analytic per-(h,k)
   l-bounds, no full-cube scan), **corrected for the Bravais lattice centering**
   (P/A/B/C/I/F/R). The report also gives two **meaningful resolution limits** —
   **I/σ = 1** and **CC(1/2) = 0.30** — and the **completeness within each**, so
   the noisy outer shells do not drag the headline completeness down (a 2I/σ
   limit is almost trivially satisfied and is not reported). **Mean I/σ** is
   likewise reported only within the CC(1/2) = 0.30 limit.
7. **Diagnostics** — a per-resolution-shell quality table (#obs, #unique,
   completeness, multiplicity, R(merge)/R(meas)/R(pim), <I/σ>, CC(1/2), %
   negative) and artifact flags: discordant ("alien") measurements relative to
   the orbit median, a scan for hexagonal-ice rings (3.90/3.67/3.44/2.67/2.25 Å),
   and an anisotropy check comparing the resolution at I/σ = 2 along a*, b* and
   c*.

---

## Validation against the COD

`tests/xrdspace-cod.js` downloads reflection files for **2000 COD entries**
(one per space group where possible, cached in `HKLs/cod/`), runs the
determination with each entry's published unit cell, and compares the result
with the published space group:

- **PASS** — exact space-group number match
- **NEAR** — the published group is among the zero-violation candidates
  (symmetry/setting ambiguities that absences alone cannot always resolve)
- **FAIL** — the published group was not recovered
- **SKIP** — unusable data (powder pattern, no single-crystal reflections)

```sh
node tests/xrdspace-cod.js              # all 2000 entries
node tests/xrdspace-cod.js 1100908      # a specific COD entry
node tests/xrdspace-cod.js --limit 50   # first N entries
node tests/xrdspace-cod.js --sg 14      # only space group 14
```

Latest full run (see `tests/xrdspace-report.json` and the chart
`tests/xrdspace-report.svg`):

| Crystal system | Total | PASS | NEAR | FAIL | SKIP |
|---|---:|---:|---:|---:|---:|
| Triclinic | 40 | 26 | 14 | 0 | 0 |
| Monoclinic | 202 | 116 | 82 | 2 | 2 |
| Orthorhombic | 625 | 425 | 181 | 4 | 15 |
| Tetragonal | 451 | 210 | 213 | 15 | 13 |
| Trigonal | 296 | 112 | 181 | 1 | 2 |
| Hexagonal | 186 | 73 | 101 | 1 | 11 |
| Cubic | 200 | 118 | 66 | 14 | 2 |
| **Total** | **2000** | **1080** | **838** | **37** | **45** |

**98.1 %** of the 1955 assessed entries have the published space group either
determined exactly or present among the zero-violation candidates (the test
exits non-zero only if the rate drops below 90 %). The SKIP count is 45
because powder-pattern CIFs (a `_pd_*` loop with no single-crystal `_refln_`
list) are correctly rejected instead of being misread as SHELX five-column
data.

### Systematic-absence corrections

Three bugs in how reflection conditions were derived from the space-group
operations were hurting the trigonal (and, to a lesser extent, orthorhombic)
space groups. All three are fixed in `src/analyze.js`:

1. **Centering translations were treated as screw/glide conditions.** For
   R-centred groups in the hexagonal setting the operation list contains pure
   centering vectors such as `(2/3, 1/3, 1/3)`. The old code only skipped
   *integer* translations, so every R-centred group collected spurious
   conditions and `R3`/`R-3` were out-ranked by their `R3m`/`R-3m`
   supergroups. The translation is now reduced modulo the centering lattice.
2. **Equivalent screw operations were double-counted.** The 6-fold components
   of a `6₂` screw impose the *same* `00l: l = 3n` condition as the 3-fold
   components of a `3₁`/`3₂` screw. Counting each operation separately gave
   `P62 2 2` four confirmed conditions against `P32`'s two, so `P31`, `P32`,
   `P3121` and `P3221` were repeatedly promoted to `P6222`. Axis conditions
   with the same invariant axis and order are now de-duplicated.
3. **Settings with the wrong centering could rescue a candidate.** The
   per-number settings loop scanned *all* settings, so `C 1 2/c 1` could be
   scored through its `I 1 2/a 1`/`A 1 2/a 1` descriptions, whose conditions
   happened to fit the data (`C2/m` → `C2/c`). Settings are now restricted to
   the detected Bravais centering.

Together these raised the exact-match count from **928 → 1024** on the same
2000 entries, including **trigonal 56 → 117** (now ahead of XPREP's 107 on the
same entries) and orthorhombic 352 → 414. The recovered rate is unchanged at
98.0 % — the few remaining tetragonal/hexagonal dips are genuine ties where
several space groups have zero violations and identical evidence, so the
published group stays in the candidate list (NEAR) but is not the top pick.

### Wide-set (10 000-entry) validation and further fixes

`tests/xrdspace-compare-10k.js` validates against a **fresh, stratified set of
10 000 COD single-crystal entries** (`tests/cod-picks-10k.json`), disjoint from
the 2000 above and drawn to cover all 186 space groups in the pool. It runs the
full determination on every entry (space-group accuracy + output data-quality
metrics: R(merge), completeness, d(I/σ=1), d(CC½=0.30), mean I/σ, multiplicity)
and, on a stratified 500-entry subset, a head-to-head against **Bruker XPREP**
(the XPREP pass drives the same interactive binary via `scripts/xprep-run.py`
and records the chosen space group, its R(sym) and CFOM).

```sh
node tests/xrdspace-compare-10k.js --xrd            # full 10k xrdspace pass
node tests/xrdspace-compare-10k.js --xprep --picks tests/cod-picks-10k-500.json
node tests/xrdspace-compare-10k.js --report         # merge both + summary
```

Latest run (see `tests/compare-10k-report.json`):

| Crystal system | Assessed | PASS | NEAR | FAIL |
|---|---:|---:|---:|---:|
| Triclinic | 2146 | 1283 | 861 | 2 |
| Monoclinic | 5177 | 4936 | 205 | 36 |
| Orthorhombic | 1550 | 1257 | 272 | 21 |
| Tetragonal | 483 | 287 | 188 | 8 |
| Trigonal | 423 | 232 | 185 | 6 |
| Hexagonal | 74 | 36 | 38 | 0 |
| Cubic | 135 | 61 | 70 | 4 |
| **Total** | **9988** | **8092** | **1819** | **77** |

**99.2 %** of the 9988 assessed entries have the published space group
determined exactly or present among the zero-violation candidates; **81.0 %**
are exact matches. Head-to-head on the 500-entry subset (both programs
answered): xrdspace recovers the published group for **99.0 %** (exact 82.3 %)
versus XPREP's **85.9 %** exact — xrdspace's *recovery* rate (published group
determined or listed as a zero-violation candidate) is far higher because XPREP
offers only a short candidate list. The two agree on 396 entries; the remaining
gaps are concentrated in triclinic (the P 1 vs P −1 centricity call, where the
Wilson |E²−1| distributions of centric and acentric data genuinely overlap).

This wider set exposed two further bugs, both fixed in `src/analyze.js`:

4. **Centering tie-break crashed on B/C/I-centred cells.** For primitive data
   every centering has zero violations, and the old tie-break picked the most
   *restrictive* centering (B/C/I) instead of P. A B/C/I-centred triclinic cell
   has no stored settings in the dictionary (only `P 1` / `P −1`), so the
   candidate list came back empty and the tool returned **no space group at
   all**. The tie is now broken by the *mean intensity of the forbidden
   reflections* (the correct centering's forbidden reflections sit at background
   level), and a final fallback to the primitive setting guarantees a candidate
   always exists.
5. **Confirmed-absences outranked the Laue class.** A higher-symmetry group
   (e.g. `I 41/a m d`, Laue 4/mmm) beat the correct lower-symmetry one
   (`I 41/a`, Laue 4/m) purely because its extra mirrors add trivially-
   "confirmed" absences, even when the R-merge showed 4/mmm at 53 % R(sym)
   against 4/m at 0 %. A candidate whose Laue class was *rejected* by the R-merge
   (R(sym) above the Laue-selection cap) is now ranked below every candidate
   whose Laue class is supported, before the confirmed-absence counts are
   compared.

Fixes 4–5 raised the 10 000-entry exact-match rate **79.2 % → 81.0 %**
(tetragonal PASS 110 → 287) and the 2000-entry exact match **1024 → 1080**,
with no regressions.

![xrdspace space-group determination vs COD](tests/xrdspace-report.png?v=3)

---

## Comparison with Bruker XPREP

`tests/xrdspace-xprep.js` runs the same 2000 cached COD datasets through both
**xrdspace** and **Bruker XPREP** and compares their space-group determination
head to head. For each entry the harness:

1. reads the cached `HKLs/cod/<id>.hkl`,
2. writes a minimal SHELX `.fcf` (HKLF list code 3) in `tests/xprep-work/`,
3. drives XPREP headlessly with `scripts/xprep-run.py` (pexpect), navigating
   the interactive menus to the space-group determination and reading the
   option XPREP chooses from its `.prp` log,
4. runs `analyzeHkl()` on exactly the same reflections and unit cell,
5. compares both with each other and with the published COD space group.

Requires the XPREP binary (`XPREP_BIN`, default
`../xdsgo/executables/xprep`) and `python3` + `pexpect`.

```sh
npm run test:xprep                       # all 2000 entries
node tests/xrdspace-xprep.js --limit 50  # first N
node tests/xrdspace-xprep.js --id 1501632
node tests/xrdspace-xprep.js --sg 14
```

Full 2000-entry result (`tests/xrdspace-xprep-report.json`):

| Crystal system | Assessed | xrdspace PASS | XPREP PASS | agree |
|---|---:|---:|---:|---:|
| Triclinic | 40 | 26 | 27 | 39 |
| Monoclinic | 201 | 114 | 124 | 103 |
| Orthorhombic | 614 | 352 | 251 | 227 |
| Tetragonal | 439 | 169 | 135 | 32 |
| Trigonal | 295 | 56 | 107 | 52 |
| Hexagonal | 175 | 100 | 39 | 20 |
| Cubic | 198 | 111 | 69 | 51 |
| **Total** | **1962** | **928 (47.3 %)** | **752 (38.3 %)** | **524** |

> **Note:** this head-to-head was measured *before* the corrections described
> above. With them, xrdspace's exact-match total on the same 2000 entries
> rises to **1080** (including trigonal 56 → 112 and tetragonal 169 → 210),
> so the xrdspace column here is a conservative lower bound. Re-run
> `npm run test:xprep` to refresh it — or use the wider 10 000-entry head-to-
> head above (`tests/xrdspace-compare-10k.js`), which is the current
> benchmark.

- xrdspace recovers the published group (exact, or present among the
  zero-violation candidates) for **97.7 %** of assessed entries; XPREP for
  **94.0 %** of entries whose space group it was able to determine.
- Head-to-head on exact matches: **xrdspace 542**, **XPREP 366**,
  both correct 386; the two programs agree on the space group for 524 entries.
- XPREP's default choice is the lowest-CFOM option, which on sparse or
  pseudo-symmetric data often lands on a subgroup or supergroup (e.g. `P2221`
  vs `Pbam`, `I4/mmm` vs `I-42m`). When the tools disagree, the published
  group is usually among the candidates of the "wrong" program.
- 38 entries produced no XPREP answer through this harness (28 driver
  timeouts on very large or slow inputs, 10 with no space-group table, e.g.
  XPREP `F-superlattice` warnings).

---

## Validation against macromolecular (MX) data

`tests/xrdspace-mx.js` validates against **real protein diffraction data** —
XDS_ASCII.HKL files from a synchrotron campaign, each with a POINTLESS run
(`pointless.xml`) whose `<BestSolution Type="spacegroup">` is the ground truth.
These are large macromolecular cells (V ≈ 5.6 × 10⁵ – 1.1 × 10⁶ Å³), so the
**chiral (Sohncke) restriction is active** on every dataset.

```sh
# point XRDSPACE_MX_DIR at a directory containing hkl/ and xml/ subdirs
XRDSPACE_MX_DIR=/path/to/mx node tests/xrdspace-mx.js
# or: npm run test:mx
```

On a recent run of 18 macromolecular datasets:

| Result | Count |
|---|---:|
| PASS (exact space-group match) | 16 |
| NEAR (in zero-violation candidates) | 0 |
| FAIL | 2 |
| **Correct (PASS + NEAR)** | **88.9 %** |
| **Non-chiral determinations** | **0** |

The two non-PASS cases are weak-data situations where POINTLESS itself
returned a low-confidence, low-symmetry "safe" answer; xrdspace returned a
metrically-consistent higher-symmetry group. In every case the determined
space group was **chiral** — the original motivation for the Sohncke
restriction (never reporting a non-chiral group for a protein) is fully met.

The per-dataset report (`tests/xrdspace-mx-report.json`) is written locally
and is **not committed** (it identifies the samples).

---

## Project structure

```
xrdspace/
├── src/
│   ├── xrdspace.js        # command-line interface (POINTLESS-style arguments)
│   ├── index.js           # main entry point / public API (analyzeHkl, ...)
│   ├── hkl-parser.js      # XDS_ASCII / SHELX / COD .hkl parsers
│   ├── analyze.js         # crystal system, R(sym) Laue selection, centering,
│   │                      #   systematic-absence scoring, centricity, ranking
│   ├── laue.js            # the 11 Laue classes (built from the dictionary)
│   ├── merge.js           # reflection merging, statistics, SHELX/XDS/.ins writers
│   ├── op-math.js         # symmetry-operation parsing and direct↔reciprocal math
│   ├── cell-search.js     # unit-cell database search: Niggli reduction, cell
│   │                      #   similarity, standard settings, COD + PDB clients
│   ├── pdb-lookup.js      # offline PDB lookup: build/load table, match by
│   │                      #   Niggli-reduced cell, validate space group (--valid)
│   └── space-groups.js    # dictionary of all 230 space groups (all settings)
├── scripts/
│   ├── build-pdb-table.js # download PDB cell+space-group data from RCSB and
│   │                      #   write data/pdb-cells.json (for --valid)
│   └── xprep-run.py       # pexpect driver: run XPREP headlessly, read its .prp
├── tests/
│   ├── xrdspace-cod.js    # COD validation harness (2000 entries)
│   ├── cod-picks.json     # the 2000 COD entries (id, cell, published SG)
│   ├── xrdspace-compare-10k.js  # wide 10k COD validation + XPREP head-to-head
│   ├── cod-picks-10k.json # the 10000 stratified COD entries
│   ├── cod-picks-10k-500.json   # stratified 500-entry XPREP head-to-head subset
│   ├── xrdspace-xprep.js  # head-to-head xrdspace vs Bruker XPREP comparison
│   ├── xrdspace-mx.js     # macromolecular (MX) validation harness
│   ├── xrdspace-cellsearch.js  # offline tests of Niggli reduction / similarity
│   ├── xrdspace-pdbvalid.js    # offline tests of the PDB lookup (--valid)
│   ├── xrdspace-report.json  # latest COD validation results
│   ├── xrdspace-report.svg   # PASS/NEAR/FAIL chart per crystal system
│   └── xrdspace-report.png   # PNG render of the COD chart
├── data/                  # generated PDB lookup table (git-ignored, ~40 MB)
├── package.json
├── LICENSE                # MIT
└── README.md
```

---

## License

MIT — see [LICENSE](LICENSE). Copyright (c) 2026 Denis Spasyuk.
