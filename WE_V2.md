# wE version 2: what it computes, and how to use it in your own software

`wE_v2.py` (this directory) is a Python port of `molprobify_v2.com`, the version-2 weighted energy of the Untangle Challenge
(1aho / has5003). It reproduces the tcsh scorer's wE to the printed digit (validation below). This file gives the complete
formula, so you can also compute wE v2 from your own geometry code without phenix.

- **Source of truth:** octamus1 `/petabyte/home/jamesh/projects/build_challenge/1aho/claude_numbers/wE_v2/molprobify_v2.com`
  (md5 e247bc5f, 24 Nov 2026 chain + Rama-modulo fix; fe6e46ab / c461987c before; copy: `erudition/molprobify_v2.com`). It is
  a copy of the Challenge scorer `/petabyte/home/jamesh/Develop/molprobify_runme.com` (v1, md5 9397c563 since the 24 Nov chain
  fix; da46cf17 before) with the revisions below.
- **Evidence and history:** `NUMBERS.md` `[N:wE-v2]` (tests on the Challenge models, calibration of every constant).
- **Status:** the official Challenge number is still v1 (`molprobify_runme.com`, "weighted energy (wE)"). v2 is an
  exploration scorer; the counter-fit session uses it as its factory scorer since 15 Nov 2026. Say which one you report.

## Quick start

```bash
# scores in-process via mmtbx under phenix.python (~5.5 s/model; falls back to phenix CLI subprocesses if keep=1 or workdir=DIR)
phenix.python wE_v2.py model.pdb [model2.pdb ...] [ligand.cif ...] [key=value ...]
# score outputs you already have (no phenix needed; ~0.38 s on 1aho; the pdb is read only for occupancies)
python3 wE_v2.py model.pdb geo=model.geo rama=ramalyze.log rota=rotalyze.log cbdev=cbetadev.log coot=molprobity_coot.py
# extras: items=items.txt (every weighted item), json=out.json, workdir=DIR keep=1 (run CLI phenix and keep outputs), backend=mmtbx|cli, jobs=3
```

```python
import wE_v2
r = wE_v2.score_model('model.pdb')                                  # in-process mmtbx under phenix.python
rs = wE_v2.score_models(['m1.pdb', 'm2.pdb'])                       # batch scoring in one warm process
ms = wE_v2.ModelScorer('template.pdb')                              # build topology/restraints once (~1.7-2.7 s)
r = ms.score_sites(new_sites_cart)                                  # ~3.7 s with Probe, ~1.6 s with include_clashes=False
r = wE_v2.score_files(pdb, geo, rama, rota, cbdev, coot)            # existing outputs (~0.38 s)
r['wE']; r['terms']['RAMA']; r['items'][:10]; r['nb_phys']
r = wE_v2.score_items(items, averages, clashes)                     # your own deviates (see "The formula")
```

The last line printed is `weighted energy (wE): X`, as in the tcsh scorer, so existing log parsers work (in multi-PDB batch mode
each model prints `=== <pdb> ===` followed by its table and `weighted energy (wE): X`, and `json=` writes a list of summaries).
Needs only the Python 3.6+ standard library for `score_files`, `score_items`, and `backend=cli`. On octamus1 nodes there is no
`python3`: use `/programs/phenix-2.1rc2-6037/phenix_bin/phenix.python` (also the phenix on the default PATH there).

Options (defaults are the tcsh defaults): `nbmode=occ|bonds|stat|physical|both`, `nbN=occ|vdw120|vdw100`, `nboutw=0|1`,
`nbphysw=0.0`, `cbsigma=0.032`, `cbdof=3|2`, `cbmode=lovell|cdl`, `clashmode=probe|smooth|none`, `occw=0|1`, `nscale=1.0`,
`stable=0|1`, `pnn=exact|soft`, `csda=0.6`, `cdl=True|False`, `cdl_interp=0|1`, `ignore=ARGH[,HOHH,HBON]`,
`nbP0=1.045e-3`, `nbk=0.667`, `sigomega=4.01398721805631`, `round=1|0` (`round_numbers=True|False`), `backend=auto|mmtbx|cli`,
`jobs=1..5` (concurrent subprocesses when `backend=cli`). Invalid option names, invalid option values, unknown positional
arguments, incomplete `score_files` CLI inputs (`geo=`, `rama=`, `rota=`, `cbdev=`), or single-model output options
(`workdir=`, `items=`, `fullgeo=`, `clashes=`, `cbdevlog=`, `outprefix=`) combined with multiple PDB files raise a `ValueError`
(API) or exit non-zero (CLI). `overridefile=FILE` (and any `override*=` key) and counter-fit's files are described under
"Counter-fit files and overrides" below; the low-scatter/smooth options (`stable=1`, `nscale=`, `occw=`, `clashmode=`,
`nbphysw=`, `cbmode=`) are described under "Stabilization options (factory / trajectory mode)". Not ported: the refmac/phenix
fudge-restraint files (`writefudge`), Amber rst7 input, and the R-factor columns of `_score.txt`; none of them enters wE.

### Python API, return schema, and runtime behaviour

- `score_model(pdb, cifs=(), workdir=None, keep=False, fullgeo=None, clashes=None, cbdevlog=None, outprefix=None, **opts)` —
  with `backend='auto'` (default), runs everything in-process through `mmtbx` when imported under `phenix.python` (unless
  `workdir` or `keep` is set), including `overridefile=`, `cdl_interp=0|1`, and writing counter-fit's `fullgeo`/`clashes`/`cbdevlog`
  files, eliminating 5 subprocess startups and duplicate `pdb_interpretation`/MolProbity/CaBLAM passes while matching
  `molprobify_v2.com` and `code_v2/molprobify_runme.com` to the printed digit and byte. With `backend='cli'` (or `workdir`/`keep`),
  runs the five Phenix CLI programs (`jobs=N` runs up to `N` concurrently, heaviest first); on Linux when neither `workdir` nor
  `TMPDIR`/`TEMP`/`TMP` is set, temporary runs go under `/dev/shm/$USER/wE_v2_*` (avoiding node-local `/tmp`). Reusing an
  existing `workdir` unlinks prior logs/`.geo`/`molprobity_coot.py` before running. Non-zero exit codes from `ramalyze`,
  `rotalyze`, `cbetadev`, or `geometry_minimization` raise `RuntimeError`; if `phenix.molprobity` / `clashscore` fails, a
  `WARNING:` is written to `stderr` and `CLASH` scores with 0 probe clashes.
- `score_models(pdbs, cifs=(), **opts)` — scores multiple PDB files in one warm Python process, importing `mmtbx`/rotamer/Rama
  tables and parsing ligand `.cif` files once (honors `overridefile=` / `override*=` and `cdl_interp=0|1`).
- `ModelScorer(pdb, cifs=(), **opts)` — builds `mmtbx.model` + geometry restraints + CDL proxies once for a fixed topology;
  `score_sites(sites_cart=None, include_clashes=True, fullgeo=None, clashes=None, cbdevlog=None, outprefix=None, **override_opts)`
  updates coordinates, riding H (`idealize_h_riding`), and CDL targets in memory (`~3.7 s` with Probe clashes, `~1.6 s` with
  `include_clashes=False`, `~5.0 s` when also writing `outprefix=` files).
- `score_files(pdb, geo, rama, rota, cbdev, coot=None, **opts)` — scores precomputed Phenix files with zero external calls
  (`read_occupancies(pdb)` also skips `XP` records directly).
- `score_items(items, averages, clashes=(), nb_phys=0.0, opts=None)` — computes the v2 sum from `{category: [E_1, ...]}` (or
  weighted `[(E_1, w_1), ...]`, and optional `(seq, N_override)` tuple) and `{name: (mean, n)}` (average keys are
  case-insensitive, so both `'clash'` and `'CLASH'` apply the `+ n` clash offset).
- Returned dictionary (`r`) and `json=` schema:
  - `r['wE']`: `float` (`g6`-rounded when `round=1`, full `float64` when `round=0`).
  - `r['terms']`: `{category: {'outliers': float|None, 'n_gt_001': int|None, 'N': int|float|None, 'average': float}}`.
  - `r['items']`: `[(weighted_val, category, k, N, Pnn, energy), ...]` sorted descending (`k` = 1-based energy rank within
    `category`, matching column 3 of `<prefix>_v2items.txt`; omitted from `json=` for size, written by `items=FILE`).
  - `r['weighted']`: `[(weighted_val, 'avg_<name>' | 'clash' | 'nb_occ_phys'), ...]`.
  - `r['nb_phys']`: `(sum_soft_repulsive_lj, n_repulsive_contacts)`; `r['nb_occ_phys']`: occupancy-weighted sum
    `Σ min(occ_a, occ_b)·soft(LJ⁺)`.
  - `r['nonbond']`: `{'N': ..., 'total_occupancy': ..., 'nb_occ_phys': ...}` (in `occ`/`bonds` modes), `r['n_clashes']`: `int`,
    `r['options']`: normalized options dict.

### Counter-fit files and overrides

Ports of counter-fit's factory scorer `code_v2/molprobify_runme.com` (md5 `457d87ea` = `molprobify_v2.com` + `cdl=`,
`cdl_interp=`, and `symop=` in fullgeo), for its `fudgeify.py`. All of these work in-process under `backend=mmtbx` (the default
under `phenix.python`, **~6.7 s warm / ~8.0 s cold on 1 CPU**, vs 24 s on 5 CPUs via CLI) as well as under `backend=cli`.

- `overridefile=FILE` (also any `override*=` key on both the CLI and Python API, as in the tcsh): lines
  `KEY OVERRIDE ideal sigma | name_resnum ...`, e.g. `TORSION OVERRIDE -100 15 | N_12 CA_12 CB_12 SG_12`. A restraint whose key
  and atoms (name_resnum, so every altloc) match gets the lowest energy its overrides allow; this changes wE (including in
  batch `score_models` and `ModelScorer`). A TORSION override also drops that residue's RAMA (backbone atoms only), ROTA (any
  side-chain atom) or OMEGA (CA…CA) items. A missing file is ignored, as in the tcsh.
- `cdl_interp=0|1`: when `1`, interpolates CDL backbone bond/angle ideals between `phi`/`psi` bins (`interpolate=True` in
  `mmtbx.conformation_dependent_library.update_restraints`, matching `cdl_interp.py`); `0` (default) uses standard CDL bins.
- `fullgeo=FILE`: `<prefix>_fullgeo.txt` as the tcsh writes it with `keepgeo=1`: `key energy delta obs ideal sigma | atom alt
  resname chain resnum - ...`, H kept, ` symop=<op>` on symmetry-mate contacts and bonds, after overrides, ordered in pure
  ASCII byte order without mutating process locale (matches the tcsh `sort -k2gr` run in the C locale).
- `clashes=FILE`: `<prefix>_clashes.txt` (`CLASH ljenergy -overlap |atom16|atom16`, from a `flip_symmetric_amino_acids=True`
  hierarchy). `cbdevlog=FILE`: `<prefix>_cbetadev.log` in exact `phenix.cbetadev` format (model basename in column 1, no stderr
  banner in either `mmtbx` or `cli` mode). `outprefix=P`: all three, as `P_fullgeo.txt`, `P_clashes.txt`, `P_cbetadev.log`.
- Python: the same keywords on `score_model`, `ModelScorer.score_sites`, and `score_files`, or
  `write_counterfit_files(geo, coot, cbdev, overridefile, fullgeo=, clashes=, cbdevlog=, outprefix=)` on its own.
- Factory call equivalent to `molprobify_runme.com MODEL.pdb keepgeo=1 cdl=True cdl_interp=0 overridefile=overrides.txt`
  (runs in-process on 1 CPU in ~6.7–8.0 s, or pass `cdl_interp=1` when interpolating CDL bins):
  `phenix.python wE_v2.py MODEL.pdb cdl=True overridefile=overrides.txt outprefix=MODEL`.

### Stabilization options (factory / trajectory mode)

In default `wE_v2`, the largest of `N` noise deviates lands where `Pnn(z, N) ≈ 0.5` (the steepest slope of the outlier weight,
`NUMBERS.md` `[N:v2-stability]`), producing ~9× larger seed scatter (`sd ≈ 0.59–0.70`) and 3-decimal PDB round-off jitter
(`sd ≈ 0.55`, span `1.81`) than v1, while unweighted conformer sums inflate `wE` linearly with altloc count `K`, Probe's
`-0.400 Å` boundary jumps `0 → 2.93`, and mild nonbond bumps (`LJ < 3`) get `Pnn ≈ 0`. Six opt-in keywords address each source
while keeping all defaults identical to `molprobify_v2.com`:

- `nscale=1.0` (`float > 0`): multiplies each category's effective `N` inside `Pnn(z, nscale·N)` (leaving `chisq` average degrees
  of freedom unchanged). At `nscale=10`, noise maxima sit in the flat tail of `Pnn`, cutting `1aho` seed `sd` from `0.59–0.70` to
  `0.07–0.11` and removing the `Thr57` false-better (`-1.41 → +0.01`).
- `occw=0|1`: weights every category's outlier contributions *outside* `Pnn` (`w · Pnn(√E, N_eff) · soft(E)` with
  `N_eff = Σ w`, avoiding `nboutw=1`'s partial-occupancy clash hiding) and weights category means by `w = min_atom(occ)`.
  Splitting a conformer into `K` identical copies of occupancy `1/K` is then strictly score-invariant across all 10
  `mmtbx`/MolProbity categories (`41.345952` for every `K ∈ {1, 2, 4, 5, 10}`, vs `41.35 → 27.16` at `occw=0`).
- `clashmode=probe|smooth|none`: `probe` (default) uses v1/v2's `1 + n·worst` formula; `smooth` subtracts Probe's `-0.400 Å`
  threshold energy `lj(2.6, 3)` (`C^0`-continuous onset `0.000000` across `-0.400 Å` instead of `0 → 2.9308`) and scores each
  clash individually with `Pnn(z, nscale·N_nb)` (suppressing `[N:balanced-ens]`'s spurious `N = 6..12` Probe/Reduce switching);
  `none` omits the Probe `CLASH` term (and skips running `clashscore` in `mmtbx` unless `clashes=`/`outprefix=` is requested,
  halving `score_model` time to `~2.7 s`).
- `nbphysw=0.0` (`float >= 0`): adds `nbphysw · nb_occ_phys` (`nb_occ_phys = Σ min(occ_a, occ_b)·soft(LJ⁺)`) so mild packing
  bumps (`LJ < 3`: `s15_rw`, `blk_rw`, `RoPE`, `Accepted2`, `island`) separate cleanly without inflating conformer splits.
- `cbmode=lovell|cdl`: `lovell` (default) uses `phenix.cbetadev` fixed targets; `cdl` constructs `CB_ideal` analytically from
  each conformer's own `.geo`/`mmtbx` `CA–CB` bond and `N–CA–C`, `C–CA–CB`, `N–CA–CB` angle ideals (works in `mmtbx`, `cli`, and
  stdlib `score_files`, and writes CDL `_cbetadev.log` when `cbdevlog=`/`outprefix=` is set). Note: keep `cbmode=lovell` when
  comparing unrefined over-parameterized ensembles such as 10-conformer `cacbtb.pdb`, whose CDL residuals are artificially small.
- `stable=0|1`: convenience preset setting `round=0, cdl_interp=1, nscale=10, clashmode=smooth, occw=1, nbphysw=0.1` (any
  explicit `key=value` overrides individual preset keys). `stable=1` places `wE` on its own scale (`wE0 best` `10.96` vs default
  v2 `54.29`); never mix `stable=1` and default v2 numbers.

| metric / model test (`1aho`, 364 models of `[N:v2-stability]`) | default `wE_v2` | `stable=1` |
|---|---|---|
| `best_wE0` (9 seeds, mean ± sd) | 53.22 ± 0.592 | **10.85 ± 0.067** (8.8× lower sd) |
| `bestH1` (16 seeds, mean ± sd) | 53.01 ± 0.701 | **10.88 ± 0.105** (6.7× lower sd) |
| Sub-mÅ `0.001 Å` PDB round-off jitter (12 draws, sd / span) | 0.5501 / 1.8065 | **0.0721 / 0.2757** (7.6× lower sd) |
| 21 `wE0_dist` rows with mean below `best` (`[N:v2-stability]`) | 3 (`cac_otw39`, `cac_best`, `cac_otw1`) | **none (`0 / 20`)** |
| 16 relabel / model arms with negative mean diff vs control | 4 (`Thr57 −5.1`, `cac_best −6.8`, `s14 −1.8`, `s15 −1.1`) | **none (`0 / 16` negative)** |
| `Thr57` / `Ser40 OG` / `Gly4 O` / `Cys48` / `Ala45 O` (`t`) | −5.1 / +1.1 / +0.2 / +1.1 / −0.1 | **+0.2 / +0.1 / +0.9 / +0.7 / +1.3** (ties) |
| `Cys46` plane / `s14_rw` / `s15_rw` / `blk_rw` / `b14_rw` (`t`) | +1.7 / −1.8 / −1.1 / +1.7 / +2.7 | **+3.1 / +8.6 / +4.3 / +5.3 / +10.7** |
| `RoPE` / `Accepted2` / `island` / `terwilliger6` (`t`) | +4.7 / +6.3 / +9.1 / +19.4 | **+5.4 / +10.9 / +10.2 / +24.3** |
| `cf9` (`best` / `rope` / `cacbtb` / `fused` / `otw39` / `lotswrong`) | 54.29 / 54.78 / 61.58 / 110.96 / 219.60 / 836.41 | **10.96 / 11.22 / 13.23 / 27.98 / 72.48 / 290.45** |
| `b5k1B` ensembles `N = 1 → 10 → 20` (`[N:balanced-ens]`) | 226.2 → 1441.2 (`cl=25`) → 2848.5 (`cl=0`) | **88.74 → 149.77 → 176.85** (`clashmode=none`: `144.14` at `N=10`) |

## The formula

Every term is an energy on the scale of a squared Gaussian deviate, E = z². Two weights turn energies into wE:

- **Outlier weight** (probability not noise): Pnn(z, N) = erf(|z|/√2)^N, the probability that a deviate this large is not the
  largest of N noise draws. It is ≈ 0 for an ordinary item and → 1 for a real outlier.
- **Average weight** (chi-square): chisq(x, k) = P(χ²_k < x) = igamma(k/2, x/2) for k < 3200, else erf(√x − √k)/2 + 0.5.
  For a category with mean energy m over n items, x = n·m: ≈ 0 when m ≈ 1 or less, → 1 when the whole category is too wide.
- **Softening:** soft(E) = E for E ≤ 10, else 10 + ln(E/10). Applied to every outlier and every average.

```
wE = Σ_categories Σ_items(E_i > 0)  Pnn(√min(E_i,400), N_c) · soft(E_i)        (outliers: v2 revision 1)
   + Σ_averages  chisq(n·clip(m,0,400), n) · soft(max(m,0))                   (as v1)
   + clash term                                                              (as v1)
   [+ NB_PHYS if nbmode = physical or both]
```

N_c is the number of items in the category (zeros included). v1 applied Pnn to the single worst item of each category; v2 applies
it to every item, so each real outlier pays and a category of pure noise pays about what v1's worst term did.

### Categories and their energies

| category | items (source) | energy E |
|---|---|---|
| BOND | every bond restraint of `phenix.geometry_minimization macro_cycles=0` .geo, no H atoms | .geo residual, (Δ/σ)² |
| ANGLE | every angle restraint, no H | .geo residual |
| TORSION | every dihedral restraint, no H | .geo residual |
| CHIR | every chirality restraint, no H | .geo residual |
| PLANE | every atom of every plane restraint (one item per atom), no H | (δ/σ)² of that atom |
| OMEGA | every CA–C–N–CA dihedral of the .geo | [(sin ω / sin 4.014°)² + (1 + cos ω)^10] / (2·nPro + 1), nPro = prolines at residue number ±1 |
| RAMA | every residue line of `phenix.ramalyze`, any chain (chain A only before 24 Nov 2026; patch note 10) | z² with p = score%/100 = erfc(z/√2), i.e. 2·inverf(1 − p)²; 99 if p rounds to 0 |
| ROTA | every line of `phenix.rotalyze` with occupancy > 0 | same as RAMA |
| CBDEV | every line of `phenix.cbetadev` | x = dev/0.032 Å is a chi(3) deviate; E = Gaussian z² with the same two-sided tail: P = erfc(x/√2) + √(2/π)·x·e^(−x²/2), E = invnorm(P/2)² (E = x² for x > 37) |
| NONBOND | see below | LJ in well depths |

"No H" means: an item that involves any atom whose name starts with H is dropped (non-bonds keep their H). The .geo comes from
`phenix.geometry_minimization model.pdb macro_cycles=0 stop_for_unknowns=false allow_polymer_cross_special_position=True
const_shrink_donor_acceptor=0.6`.

### Non-bonds (v2 revision 3; default `nbmode=occ`)

- Each .geo non-bonded pair with model distance r and phenix vdW distance r0 gets
  LJ = lj0(r, r0) − lj0(6, r0), lj0(r, r0) = 4[(s/r)^12 − (s/r)^6], s = r0·2^(−1/6).
  This puts the LJ minimum (−1) at r0 and the zero crossing at 0.89·r0. Repulsive contacts are those with LJ > 0.
- H–H contacts between two ARG residues are dropped (`ignore=ARGH`).
- **Outlier items:** every repulsive contact's LJ (unweighted), padded with zeros to N = round(total atom occupancy of the
  model, waters and H included). N no longer counts repulsive contacts, which was v1's NONBOND-n flaw (an extra bump lowered wE).
- **Average term:** m = Σ w·LJ⁺ / N, w = min(occ_a, occ_b), with the usual chi-square weight (≈ 0 for every model tested:
  m ≈ 0.01–0.06).
- **NB_PHYS** = Σ soft(LJ⁺) over all contacts is printed, but enters wE only with `nbmode=physical` or `both`.
- Other modes: `bonds` (N = number of non-H bonds; average = Σ LJ⁺ / N_bonds); `stat` (every contact an item, N = all
  contacts, E = invnorm(P0·e^(−k·LJ)/2)² for LJ > 0 with P0 = 1.045e-3, k = 0.667); `nbN=vdw120|vdw100` (N = occupancy-
  weighted number of contacts within 1.2 / 1.0 × r0); `nboutw=1` (outlier items w·LJ).

### Averages and the clash term (unchanged from v1)

- Averages: BOND, ANGLE, TORSION, CHIR, PLANE, OMEGA, RAMA, ROTA, CBDEV (mean of the same item energies, n = their number),
  NONBOND (above), CLASH.
- Clashes: probe clashes from `phenix.molprobity ... output.coot=True` (`molprobity_coot.py`); each with overlap d gets
  E = LJ(3 + d, 3). With n clashes and worst E_w:
  clash term = softPnna(√min(1 + n·E_w, 400), n) · (soft(1 + n·E_w) + n), and the CLASH average adds
  chisq(n·clip(m,0,400), n) · (soft(m) + n). softPnna is v1's approximation: 1 − 2^(−(|z|/a(N))^e(z,N)), a = −0.0192266·L² + 0.751694·L
  + 1.12482, e = 0.736621·L + 0.21805·z + 1, L = log10 N. Clash-free models have no clash term.

### Constants and their origin (NUMBERS `[N:wE-v2]`)

- cbsigma 0.032 Å per component: rms/√3 of 57,000 Cβ deviations in 500 truth re-refinements.
- Rama/rota: z² with the two-sided relation (v1 used inverf(1 − p)² = z²/2).
- Omega σ 4.014°, softening at 10, cap 400 for the weights: as v1.

### Rounding

The tcsh scorer passes numbers between awk and gnuplot as text, so most intermediate values are rounded to 6 significant digits
(gnuplot prints 15). `wE_v2.py` copies this by default (`round=1`) so its wE matches the tcsh output digit for digit; pass
`round=0` (`round_numbers=False`) to disable intermediate rounding and the `%g == 1` Rama/Rota 99-cap. If you implement the
formula in your own code without that rounding, expect differences in the last printed digit only.

## Behaviour you should know about

- **Near the truth, v2 is mostly terms that correct models pay too.** wE0 best.pdb scores 54.29: RAMA 22.24, CBDEV 15.26
  (the truth itself scores Cβ 0.68), ROTA 12.41, TORSION 2.35, OMEGA 2.03, NONBOND 0. v1 gives ≈ 18 for the same models, so
  never compare a v2 number with a v1 number.
- **Seed scatter near the truth is large: sd 0.36–0.70** for correct models (v1: 0.05–0.11; `stable=1`: 0.07–0.11): 0.36–0.60
  across best, c-a-c best and RoPE (`[N:wE-v2]`, 27 Oct; best 0.59 over 9 seeds) and 0.70 for best_H1 over 16 seeds
  (`[N:v2-stability]`, 17 Nov). Almost all of it is the Cβ term: refined correct models'
  largest Cβ deviations (0.115–0.122 Å) sit on the steep part of the weight curve at σ 0.032 Å. Compare models with several
  refinement seeds, not one.
- **The LJ distance is v1's:** the well minimum sits at phenix's r0 (the onset of phenix's repulsion), so the zero crossing is
  0.3–0.4 Å short and overlaps up to ≈ 0.3 Å score as attractive (`[N:ser40-ljfix]`). With the minimum moved to 1.12·r0
  (σ = r0), the geometry-fixed Glu32 C···HG2 pair (2.36 Å vs r0 2.92) scores LJ 36 in every model and saturates NONBOND; such
  pairs would need excluding first. There is no attraction, electrostatics or hydrogen-bond term.
- **Mild extra bumps are invisible in the default mode:** contacts with LJ < 3 get Pnn ≈ 0. Only NB_PHYS (an aggregate) sees
  them. For example, the RoPE layout costs +3.25 ± 0.03 in the physical non-bond sum against best_H1's layout (worse in 32/32
  seeds), but wE v2 only +0.25 ± 0.16, which is not significant (`[N:wE-v2]`).
- **The phenix version is part of the score.** The geometry energies are phenix's restraint residuals (its library, CDL, and
  vdW distances), so hold one phenix version fixed within a comparison. Both the tcsh scorer and `wE_v2.py` use whatever phenix
  is on PATH (on octamus1 that is 2.1rc2-6037), and both pass `const_shrink_donor_acceptor=0.6`.
- **Quirks kept from the scorer** (reproduced, not fixed): omega's proline neighbour test uses residue numbers and ignores
  the chain; the CLASH terms add n to the energy before weighting; the probe CLASH term double-counts contacts that NONBOND
  also charges (JH's decision, 22 Oct). RAMA read chain A only until 24 Nov 2026; it now reads every chain (patch note 10).

## Validation

15 Nov 2026, octamus1 (phenix 2.1rc2-6037 on PATH, gnuplot 4.6, phenix.python 3.9; SLURM 624678). Each model was scored by
`molprobify_v2.com ... debug` (which keeps its .geo and `molprobity_coot.py`), by `wE_v2.py` on those same files, and by
`wE_v2.py` end to end (its own phenix runs):

| model | tcsh | wE_v2.py, same files | wE_v2.py, end to end |
|---|---|---|---|
| wE0 best | 54.2861 | 54.2861 | 54.2861 |
| best_H1 seed 2001 | 51.4465 | 51.4465 | 51.4465 |
| wE0 rope | 54.776 | 54.776 | 54.776 |
| cacbtb.pdb (as provided) | 61.5809 | 61.5809 | 61.5809 |
| wE0 fused | 110.959 | 110.959 | 110.959 |
| otw39 as distributed | 219.603 | 219.603 | 219.603 |
| wE0 qfit_best (clashes) | 300.745 | 300.745 | 300.745 |
| wE0 af2_combo24 (clashes) | 599.508 | 599.508 | 599.508 |
| wE0 lotswrong (clashes) | 836.406 | 836.406 | 836.406 |

Every option, on lotswrong (tcsh = wE_v2.py on the same files): default 836.406; `nbmode=bonds` 836.266; `nbmode=stat` 879.298;
`nbmode=physical` 916.068; `nbmode=both` 1001.92; `nbN=vdw120` 831.663; `nboutw=1` 825.152; `pnn=soft` 1178.44; `cbdof=2`
890.951.

In all 17 runs, the item list (`items=`, the analogue of the tcsh `<prefix>_v2items.txt`) matched the tcsh list across all six
columns (`weighted_val`, `category`, within-category rank `k`, `N`, `Pnn`, `energy` to 6 digits; in `pnn=soft`, 12 bottom-of-list
items around `10^-15`–`10^-20` differ by 1 ULP in `1 − 2^v` cancellation between C `pow` and Python `**`, with zero effect on
any category sum or `wE`), and so did every per-category term. Additional checks on octamus1 (SLURM 624798, 624958):
- `7KR0` under Phenix 2.2-6143 (`mac077`: `34.4368`, including all-underflow `NONBOND - - - 0.0000`; `mac095`: `183.011`) matched
  `molprobify_v2.com` across `wE`, `_v2terms.txt`, and all 6 columns of `_v2items.txt`.
- `wE0 lotswrong` + 200 Amber `XP` extra-point lines scored `836.406` (`NONBOND` `N = 1016`) end to end, matching
  `molprobify_v2.com` and creating its temporary workspace in `/dev/shm/$USER` when `TMPDIR` was unset.
- `wE0 best` with `round=0` (`round_numbers=False`) gave `54.2860604`; `cdl=False` gave `54.2861` (`NB_PHYS` `20.6244` vs
  `21.0973` with CDL on).

Matching the list needs one more gnuplot quirk: an item whose Pnn underflows to exactly 0 is dropped from the list (gnuplot
prints an error instead of a value); its value would be 0 anyway, so wE is unaffected.

## Files

- `wE_v2.py`: the scorer (this directory; initial unpatched copy in `history/erudition_wE_v2/wE_v2.py`, md5 `43cc9b65`).
- `README.md`: index of this directory; `studies/`: scripts and small results of the tests behind these notes (24 Nov).
- `../erudition/molprobify_v2.com`: copy of the tcsh source of truth (md5 `e247bc5f`; `fe6e46ab` / `c461987c` before 24 Nov).
- Tests on octamus1:
  - `/scratch/jamesh/erudition_wEv2py_2026-11-15/` (initial 17-run validation, SLURM 624678).
  - `/scratch/jamesh/agy_wEv2_check_2026-11-15/` (`7KR0` Phenix 2.2-6143, Amber `XP`, `/dev/shm`, `round=0`, `cdl=False`, and
    19-run 6-column `items=` regression check; SLURM 624798, 624958).
  - `/scratch/jamesh/erudition_wEv2py_2026-11-15/cf/` (counter-fit files and overrides vs code_v2; `cfrun.sh`, `cmp_cf.sh`,
    `check2.sh`; SLURM 625274, 625293).
  - `/scratch/jamesh/agy_wEv2_check_2026-11-15/cf_inproc/` (18 end-to-end in-process `backend=mmtbx` counter-fit reference runs,
    `cdl_interp=1` reference + `mmtbx`/`cli` checks, and `score_models`/`ModelScorer` override checks; SLURM 627890).
  - `/scratch/jamesh/agy_wEv2_stability_2026-11-17/` (364-model `[N:v2-stability]` cross-check + 164-model seed-scatter & relabel
    tables, `K`-split conformer invariance, `-0.400 Å` Probe continuity, `b5k1B` `N=1..20` ensembles, and sub-mÅ coordinate
    round-off jitter; SLURM 634794, 635015, 636319, 636681).
  - `/scratch/jamesh/wEv2_git_2026-11-24/` (blank-chain / 4-digit resnum `RAMA_LINE` + `RAMA_OVERRIDE_MODULO = 9999` Mac1
    parity check; SLURM 31115245).

## Patch notes (15–24 Nov 2026)

Changes applied to `wE_v2/` over the initial port (`43cc9b65`):
1. **Linux `/dev/shm/$USER` default temp dir (`_mkdtemp`):** avoids writing `.geo` and Phenix logs to node-local `/tmp` on
   `octamus1` when `workdir` and `TMPDIR` are unset.
2. **Amber `XP` extra-point handling:** `score_model()` no longer overwrites `f['pdb']` (`this.pdb`, `XP`-stripped) with the
   unstripped input path, and `read_occupancies()` also skips `XP` records directly for `score_files()`.
3. **Within-category rank `k` in `score_items()`:** category energies are sorted descending before `enumerate(..., 1)` so `k` in
   `r['items']` and `items=` matches `<prefix>_v2items.txt` (`1` = worst item in category).
4. **Residue `0` / negative residue numbers in `OMEGA`:** replaced Python floor modulo `% 9999` with `math.fmod(..., 9999)` to
   match `awk` remainder for N-terminal `<= 0` residue numbers adjacent to proline.
5. **Validation, safety, and precision options:** added `cdl=True|False` and `round=1|0` (`round_numbers=True|False`), Maclaurin
   `inverf(y)` series for `|y| < 10^-5`, strict option/CLI validation, `shutil.which` executable resolution (supporting Windows
   `.bat` wrappers), stale-output cleanup on `workdir=` reuse, non-zero exit checks (`RuntimeError` on required steps, `stderr`
   warning if `phenix.molprobity` omits `molprobity_coot.py`), and case-insensitive `'clash'`/`'CLASH'` handling in
   `score_items()`.
6. **Performance (`backend=mmtbx`, `ModelScorer`, batch CLI, fast `.geo` parser; SLURM 625267, 625273, 625283):**
   - **In-process `mmtbx` (`backend=auto|mmtbx`):** **6.8 s** cold / **5.5 s** warm on `1aho` (`3.9 s` on `fused`; all 9 `1aho`
     models in **~50 s** vs **~322 s** via subprocesses, **6.5–9.2× faster**), reproducing `molprobify_v2.com` to the printed
     digit across all 9 `1aho` models, all 8 `lotswrong` options, Amber `XP`, and `7KR0` (`mac077` / `mac095` under Phenix 2.2).
     With `round=0` under `mmtbx`, unquantized C++ float64 residuals/distances/scores are used directly (`53.99650357` on
     `wE0 best` vs `54.2860604` from 3-digit `.geo` text).
   - **`ModelScorer.score_sites()`:** **3.7 s**/call with Probe clashes (**9.6×**), **1.6 s**/call with `include_clashes=False`
     (**22×**) after one `1.7–2.7 s` topology build.
   - **Parallel CLI (`backend=cli jobs=3`):** **18.3 s** vs **35.8 s** serial (**1.95×**).
   - **Fast `parse_geo()` / `score_files()`:** **0.38 s** vs **0.95 s** on `wE0 best` (**2.5×**; `1.04 s` on 4,298-atom `mac077`).
7. **Counter-fit files and overrides (erudition, 15 Nov; merged onto `c7b323c9`):** `overridefile=`, `fullgeo=`, `clashes=`,
   `cbdevlog=`, `outprefix=`. Against code_v2 (`keepgeo=1 debug=1 cdl=True cdl_interp=0`, phenix 2.2-6143) on the 9 `1aho`
   models, without and with counter-fit's 1aho overrides (Cys12 chi1, Gln32/Gln37 chi3), all 18 runs gave the same wE (e.g.
   `wE0 best` 54.2861 / 42.8183) and byte-identical fullgeo, clashes, cbetadev log and items files; each fullgeo carried
   4,007–4,031 `symop=` lines and 5 overridden lines.
8. **In-process counter-fit files, `cdl_interp=1`, and merge fixes (16 Nov; SLURM 627890):**
   - Fixed `score_models()` / `ModelScorer` / multi-PDB batch CLI ignoring `overridefile=` (and rejected single-model file
     options `fullgeo=`, `clashes=`, `cbdevlog=`, `outprefix=` in multi-PDB CLI batch mode); accepted `override*=` keys in the
     Python API as well as the CLI; changed `_sort_k2gr()` to pure C-locale ASCII byte comparison without mutating global
     `LC_COLLATE`.
   - Enabled in-process `backend=mmtbx` (and `ModelScorer.score_sites`) for `overridefile=`, `fullgeo=` (`model.restraints_as_geo()`),
     `clashes=` (`cl_res.as_coot_data()`), `cbdevlog=` (`cb_res.show_old_output(..., prefix=stem)`), and `outprefix=`, and
     cleaned up `backend=cli` `cbetadev.log` (real PDB stem in column 1, `stderr` banner kept separate) so both backends write
     byte-identical files.
   - Ported `cdl_interp=0|1` in both `backend=mmtbx` (`_cdl_interp_context` injecting `interpolate=True` into
     `mmtbx.conformation_dependent_library.update_restraints`) and `backend=cli`.
   - Verified end-to-end in-process (`backend=mmtbx`, Phenix 2.2-6143) across all **18** `code_v2` reference runs (`9` `1aho`
     models × `{noovr, ovr}`): exact `wE` match and **byte-identical** `_fullgeo.txt`, `_clashes.txt`, `_cbetadev.log`, and
     `_v2items.txt` (`0` diffs across all 18 runs; **~6.7 s** warm / **~8.0 s** cold on **1 CPU**, `5.0 s` via `ModelScorer`
     with `outprefix=`, vs `24 s` on 5 CPUs via CLI). Verified `cdl_interp=1` against two `code_v2` + `cdl_interp.py` reference
     runs (`wE0_best` + overrides `42.8183`, `wE0_lotswrong` `857.206` vs `836.406` at `cdl_interp=0`), matching `wE` and all
     four output files byte-for-byte in both `backend=mmtbx` and `backend=cli`.
9. **Stabilization options (`nscale=`, `occw=`, `clashmode=`, `nbphysw=`, `cbmode=`, `stable=1`; 17–18 Nov, SLURM 634794, 635015, 636319, 636681):**
   - Added opt-in keywords `nscale=1.0`, `occw=0|1`, `clashmode=probe|smooth|none`, `nbphysw=0.0`, `cbmode=lovell|cdl`, and
     `stable=0|1` (`STABLE_PRESET`: `round=0, cdl_interp=1, nscale=10, clashmode=smooth, occw=1, nbphysw=0.1`) across
     `score_model`, `score_models`, `ModelScorer`, `score_files`, and `score_items`, keeping all default behavior unchanged
     (`364 / 364` exact match to `[N:v2-stability]` SLURM 633923).
   - Verified `mmtbx` vs `cli` `score_files` parity, exact `K`-conformer splitting invariance (`41.345952` across
     `K ∈ {1, 2, 4, 5, 10}`), `C^0` onset across Probe's `-0.400 Å` boundary, monotonic `b5k1B_n01..n20` ensemble scaling
     (`88.74 → 176.85` vs default v2 `226.17 → 2848.46`), **7.6×** lower sub-mÅ PDB round-off jitter (`sd = 0.072` vs `0.550`),
     and **8.8×** lower seed scatter (`best_wE0` `10.85 ± 0.067` vs `53.22 ± 0.592`) across all 364 models of `[N:v2-stability]`
     with **0 / 20** `wE0_dist` rows below `best` and **0 / 16** negative relabel arms.
10. **RAMA from every chain (24 Nov; JH "fix the chain-A-only bug"):** v1 and v2 kept only `phenix.ramalyze` lines starting
    with chain ` A` (`/^ A/`), so protein in any other chain had no Ramachandran items. Now every residue line counts:
    `RAMA_LINE = /^[ A-Za-z0-9][A-Za-z0-9] +-?[0-9]/` in `parse_rama` and the in-process path of `wE_v2.py` (md5 `69f3f4fa`;
    `241f26b5` before), the same regex in the tcsh `molprobify_v2.com` (`fe6e46ab`) and v1 `molprobify_runme.com`
    (`9397c563`). Models with all protein in chain A score exactly as before: best, `cacbtb.pdb`, lotswrong and counter-fit's
    2-conformer truth gave the same `wE` to every printed digit, tcsh and Python, default and `stable=1` (SLURM 31094541,
    31094548; octamus1 `/scratch/jamesh/erudition_chainA_2026-11-24/`). The Challenge truth `fmodelme.pdb` (conformer B as
    chain B) went from `24.6464` to `39.5085` (`stable=1`: `8.4895` → `10.0485`): RAMA `6.09` → `20.79`, N `62` → `124`,
    from chain B Asp8 (0.28 %), Gly17 (0.38 %) and Ser33 (0.45 %). Use `onechain.pdb` (the same coordinates as one chain) as
    the truth's reference (JH 24 Nov).
11. **Blank-chain / 4-digit resnum `RAMA_LINE` and `RAMA_OVERRIDE_MODULO = 9999` (24 Nov; SLURM 31115245):** broadened
    `RAMA_LINE` to `r'^[ A-Za-z0-9][ A-Za-z0-9] *-?[0-9]'` (scoring blank-chain PDBs and 4-digit residue numbers) and merged
    counter-fit's Mac1 fix (`RAMA_OVERRIDE_MODULO = 9999` instead of `64`, matching `$modulo` in `molprobify_v2.com` md5
    `e247bc5f` and `wE_v2.py` md5 `fdf62d77`; `1aho` `wE0_best` `54.2861`, `fmodelme.pdb` `39.5085`, and two-chain Mac1 +
    `His45`/`His86` overrides `70.9299` matched across tcsh and Python).
