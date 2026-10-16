# Plan — Timing-gate benchmark study on run 24386

**Status:** draft (iteration 1, pre-approval — awaiting user answers to §2
Open Questions before promoting to approved).
**Scope:** Sweep `benchmark_hycal_timing.py --pre-window` over run 24386
with the validated iter1 calibration; aggregate per-window metrics
(acceptance, μ, σ_E/E Gaussian + Crystal-Ball, χ²/dof proxy) into a
single TSV and trend plots; identify the optimal pre-window width.
**Hard constraint:** no new C++ or Python modules; the sweep is a thin
Python/bash driver that invokes the existing stable benchmark script
and parses its existing TSV output. No edits to
`benchmark_hycal_timing.py` itself.
**Deliverables:** aggregate TSV, three trend PNGs, per-window benchmark
artifacts under a dedicated directory, and a short markdown summary
with the optimal-window recommendation.

## Research inputs used

- `analysis/pyscripts/benchmark_hycal_timing.py` lines 237–305 (CLI:
  `--pre-window LO HI`, `--window W`, `--calibrate`, `--cal-method`,
  `--exclude-rowcol`, `--signal-min`, `--e-range`, `--bin-width`,
  `--max-events`, `--target-peak`, `--inner-radius`,
  `--min-seed-energy`, `--min-per-module`), 383–415 (pre-window applied
  as `pk.time <= pre_lo or pk.time >= pre_hi` — strict exclusive at
  both ends), 605–657 (TSV schema: calibration block has `path`,
  `n_modules`, `n_inner_events`, `mu_pre_MeV`, `sigma_pre_MeV`,
  `mu_post_MeV`, `sigma_post_MeV`, `sigma_over_E_post`, `mu_cb_MeV`,
  `sigma_cb_MeV`, `alpha_cb`, `n_cb`, `sigma_over_E_cb`; headline block
  has `path`, `n_events`, `gauss_mu_MeV`, `gauss_sigma_MeV`,
  `n_in_fit`).
- `analysis/pyscripts/_common.py` lines 81–152 (`discover_split_files`
  supports glob, directory, single-split), 204–254
  (`setup_pipeline`), 282–283 (default HC time window `[100, 200]` ns
  in library code — distinct from runinfo).
- `analysis/pyscripts/output/baseline_024386_iter1_test.tsv` lines 1–9
  (prior run: 50 000 events on split .00050, `pre_window=[100,200]`,
  `seed_time_window=4`, 33.4 s wall; headline σ=150 saturates fit
  bounds because `--signal-min 500` keeps the Compton tail — the
  physics 3.49 GeV peak is only recovered under `--calibrate`).
- `docs/technical_notes/hycal_clustering/hycal_clustering.md` lines
  328–443 (benchmark methodology, 24386 runinfo pre-window `[130, 200]`
  ns, `W = 16` ns production default, `--cal-method cb-peak` +
  `--exclude-rowcol 15 20 15 20` as the reference calibration recipe,
  σ_E/E = 1.90–1.98 % under those settings, calibration-equalised
  comparison is the methodologically correct way to isolate intrinsic
  resolution — raw μ shifts across pre-windows reflect pulse-selection
  changes under a fixed calibration, not detector physics).
- `/Users/L00338853/work/PRad/data/evio/` directory listing — 11
  splits present (`prad_024386.evio.000{10,20,30,40,50,60,70,80,90,100,110}`).
- `/Users/L00338853/work/PRad/data/physics_calib_out/Physics_calib/run24386/calib_factor_iter1.json`
  confirmed present.

---

## 1. Objective

Validate the optimal pre-window `[LO, HI]` (ns) for run-24386-style
clean PRad-II data by sweeping the lower edge (`LO`) and recording for
each configuration:

- **acceptance** — `n_events` (headline) and `n_inner_events`
  (calibration-ring sample) for the legacy path (gated path reported in
  parallel but not the primary observable, since the point of this
  study is pre-window width, not gate behaviour).
- **energy response** — raw peak μ, calibrated peak μ (post per-module
  gain refinement), both from Gaussian and Crystal-Ball fits.
- **energy resolution** — σ_E/E from Gaussian and Crystal-Ball on the
  calibration-equalised spectrum (the only comparable observable across
  pre-windows under a frozen per-module calibration).
- **fit quality** — fit-window event count `n_in_fit`, CB tail
  parameters (α, n) as a χ²-proxy on fit quality (the current script
  does not emit true χ²/dof — see §2 Open Questions).

Deliver an evidence-based recommendation for the production pre-window
for 24386 (and runs with similar beam/trigger conditions).

---

## 2. Open questions (block implementation)

The user prompt asked these explicitly. Resolve before step 4.

1. **Sweep grid.** Prompt suggests `LO ∈ {50, 100, 150, 200, 250}` ns;
   runinfo for 24386 is `[130, 200]`; the prior iter1 smoke test used
   `[100, 200]`. Confirm the grid. Proposed default (used in §4 below
   unless user overrides):

   - **Fixed `HI = 200` ns** (matches runinfo upper edge and prior test).
   - **Sweep `LO ∈ {50, 75, 100, 125, 130, 150, 175}` ns** — seven points
     bracketing the runinfo lower edge (130) with ±80 ns coverage.
     Window widths therefore are {150, 125, 100, 75, 70, 50, 25} ns.
   - Includes the runinfo exact baseline (130) and the previous smoke
     test (100) so the aggregate curve is directly comparable to the
     technical note's headline.

2. **Vary `HI` as well?** Prompt asks "Should post-window be fixed or
   also varied?". Proposed default: **fix `HI = 200` ns** for this
   study (uses the production upper edge and keeps the sweep
   1-dimensional for clean trend plots). A 2-D (`LO × HI`) scan would
   quadruple runtime and dilute the main physics question (does
   narrowing the gate from the left side degrade acceptance /
   resolution?). If the user wants a 2-D scan it is a trivial
   generalisation — see §9 "Extensions".

3. **EVIO file selection.** Prompt options: (a) representative middle
   file, (b) all 11 files separately, (c) all 11 files combined.
   Proposed default: **single representative split
   `prad_024386.evio.00050` with `--max-events 50000`**, matching the
   prior iter1 test (same ~30k physics events / split confirmed in the
   iter1 TSV). Rationale: 7 pre-window points × ~35 s wall each
   ≈ 4 minutes total; a per-split or all-files run would multiply
   runtime by 11× or more without changing the trend shape, since
   intra-run beam conditions are stable (confirmed by the technical
   note's 100 000-event-run headline on the full run matching the
   50 000-event split-50 smoke test to the σ-decimal). If the user
   wants a cross-split check we add it as an optional §9 extension.

4. **Reference baseline.** Prompt asks for a timing baseline from
   `hycal_clustering.md`. Confirmed: **runinfo pre-window for 24386 is
   `[130, 200]` ns**; the `--window` default is `4` ns in the script
   (vs. `W = 16` ns recommended in the technical note). Proposed
   default: use **`--window 16`** throughout the sweep (production
   default per the technical note). This is a fixed knob of the gated
   clusterer and does not affect the legacy headline, which is the
   primary observable of this study.

5. **"Timing bias" definition.** Prompt lists "timing bias (if
   available from benchmark output)" as a metric. The current
   `benchmark_hycal_timing.py` TSV **does not emit** per-event cluster
   timing — it only emits energy. Options:
   - (a) **Interpret "bias" as the pre-cal μ shift (`mu_pre_MeV`)
     relative to the target 3488.5 MeV**: already in the TSV, zero
     extra work. Report this in the aggregate.
   - (b) **Extend the benchmark script** to emit mean cluster time per
     path (requires editing `benchmark_hycal_timing.py` — out of this
     plan's scope; see §9).
   - (c) **Use `study_hycal_timing.py` separately** for a one-shot
     timing distribution at the chosen optimal window (post-hoc).
   Proposed default: **(a) for the sweep aggregate, (c) as a follow-up
   validation at the chosen optimum**.

6. **Calibration file wiring.** The script accepts
   `--hc-calib-file`. Proposed default: pass
   `--hc-calib-file /Users/L00338853/work/PRad/data/physics_calib_out/Physics_calib/run24386/calib_factor_iter1.json`
   explicitly on every invocation to lock the per-module gains across
   the sweep. (Without this the pipeline would fall back to runinfo,
   which may not point at iter1.)

---

## 3. Requirements summary

### Confirmed from user prompt

- Dataset: **run 24386**, 11 EVIO splits in `/Users/L00338853/work/PRad/data/evio/`.
- Calibration: **iter1** at
  `/Users/L00338853/work/PRad/data/physics_calib_out/Physics_calib/run24386/calib_factor_iter1.json`.
- Driver: existing `analysis/pyscripts/benchmark_hycal_timing.py`.
- Target σ_E/E observable: **Crystal-Ball + Gaussian under
  `--calibrate`** (both already emitted).
- Data cuts: single-cluster, `E > 500` MeV headline (default
  `--signal-min`), inner-ring (`r < 200` mm) exclusion for the
  per-module calibration, PbWO₄ beam-hole block
  `--exclude-rowcol 15 20 15 20`.
- Deliverables and output paths as listed in prompt §6.
- Success criterion: **optimal pre-window maximises σ_E/E** (i.e.
  minimises it) while maintaining **>90 % event acceptance** relative
  to the widest window in the sweep.

### Inferred (confirm as part of §2 open questions)

- Fixed `HI = 200` ns, swept `LO`, 7-point grid.
- `--window 16` ns (production gated-cluster default).
- `--cal-method cb-peak`, `--exclude-rowcol 15 20 15 20` (technical
  note's reference calibration).
- Single split `prad_024386.evio.00050`, `--max-events 50000`.
- Headline reference metric: **legacy σ_E/E (CB)**; gated reported in
  parallel but not the knob being swept.

### Acceptance criteria

- Aggregate TSV contains one row per pre-window configuration with all
  nine metrics (LO, HI, width, n_events, n_inner, μ_pre, σ_pre,
  μ_post, σ_post, σ/E_G, μ_CB, σ_CB, σ/E_CB, α_CB, n_CB) for both
  legacy and gated.
- Three PNG trend plots are generated with legacy + gated overlays and
  the runinfo baseline marked at `LO = 130`.
- Markdown summary identifies the recommended pre-window and quantifies
  the degradation slope per 25 ns of LO narrowing.
- All intermediate per-window benchmark outputs (one `.tsv`, one
  `.png`, one `_calibrated.png` per LO) survive in
  `output/timing_gate_benchmark_24386/`.

---

## 4. Experimental design

### 4.1 Window sweep configurations

| idx | `--pre-window LO HI` | width (ns) | notes |
|---|---|---:|---|
| 0 | `50 200`  | 150 | widest; acceptance baseline |
| 1 | `75 200`  | 125 | |
| 2 | `100 200` | 100 | prior iter1 smoke-test point |
| 3 | `125 200` |  75 | |
| 4 | `130 200` |  70 | **runinfo baseline** — reference row |
| 5 | `150 200` |  50 | moderate narrowing |
| 6 | `175 200` |  25 | aggressive narrowing — expect degradation |

### 4.2 Fixed parameters (same across all 7 runs)

- `--window 16` (gated clusterer seed_time_window; constant per §2.4)
- `--signal-min 500` (default; keeps the full physics spectrum for the
  headline Gaussian — the calibrated 3488.5 MeV peak is extracted
  separately via `--calibrate`)
- `--calibrate` (always on — all metrics of interest live in the
  calibration block)
- `--cal-method cb-peak`
- `--target-peak 3488.5`
- `--inner-radius 200`
- `--min-seed-energy 200`
- `--min-per-module 30`
- `--exclude-rowcol 15 20 15 20`
- `--bin-width 20`
- `--e-range 2500 4000`
- `--max-events 50000`
- `--hc-calib-file /Users/L00338853/work/PRad/data/physics_calib_out/Physics_calib/run24386/calib_factor_iter1.json`
- `evio_path = /Users/L00338853/work/PRad/data/evio/prad_024386.evio.00050`

### 4.3 Per-window invocation (template)

```bash
python analysis/pyscripts/benchmark_hycal_timing.py \
    /Users/L00338853/work/PRad/data/evio/prad_024386.evio.00050 \
    analysis/pyscripts/output/timing_gate_benchmark_24386/pre_${LO}_${HI} \
    --max-events 50000 \
    --pre-window ${LO} ${HI} \
    --window 16 \
    --calibrate --cal-method cb-peak \
    --exclude-rowcol 15 20 15 20 \
    --target-peak 3488.5 \
    --hc-calib-file /Users/L00338853/work/PRad/data/physics_calib_out/Physics_calib/run24386/calib_factor_iter1.json \
    2>&1 | tee analysis/pyscripts/output/timing_gate_benchmark_24386/pre_${LO}_${HI}.log
```

Each invocation emits:
- `pre_${LO}_${HI}.tsv` — the structured TSV with headline + calibration blocks
- `pre_${LO}_${HI}.png` — headline spectrum, legacy + gated Gaussian overlays
- `pre_${LO}_${HI}_calibrated.png` — inner-ring pre/post-calibration spectrum
  with Crystal-Ball overlay
- `pre_${LO}_${HI}.log` — stdout/stderr capture (progress + headline table)

### 4.4 Aggregation

A short Python driver (`sweep_pre_window_24386.py`, see §6) will:

1. Walk the list of (LO, HI) pairs.
2. Shell out to `benchmark_hycal_timing.py` with the fixed args.
3. After each run completes, parse the emitted `.tsv`, extract the
   headline and calibration rows for both `legacy` and `new`, and
   append one row per path to the aggregate TSV.
4. After all runs complete, read the aggregate back and produce the
   three trend PNGs.

The driver is a new standalone script — not an edit to
`benchmark_hycal_timing.py`. Hand off implementation to
`script-implementer` after plan approval.

### 4.5 TSV parsing contract (brittleness note)

`benchmark_hycal_timing.py` emits the calibration block only when
`--calibrate` is set. The driver must assert the calibration section is
present and fail loudly if it is not, to avoid silent NaN rows in the
aggregate. Section delimiters are the comment-prefixed lines
(`# --calibrate inner ring ...`); the driver can key off those as
anchors.

---

## 5. Aggregate TSV schema

Path: `analysis/pyscripts/output/timing_gate_benchmark_24386_summary.tsv`

Columns (one row per (path, pre-window) combination; 7 LOs × 2 paths =
14 rows):

```
pre_lo_ns   pre_hi_ns   width_ns   path            \
n_events    n_inner_events   n_modules_used    \
mu_pre_MeV   sigma_pre_MeV   mu_post_MeV   sigma_post_MeV   sigma_over_E_G   \
mu_cb_MeV    sigma_cb_MeV    alpha_cb    n_cb    sigma_over_E_CB   \
acceptance_vs_widest_pct   gauss_mu_headline_MeV   gauss_sigma_headline_MeV   n_in_fit
```

Where:
- `acceptance_vs_widest_pct = 100 × n_events / n_events_at_LO=50`
  (per path).
- `gauss_mu_headline_MeV`, `gauss_sigma_headline_MeV`, `n_in_fit` come
  from the TSV's headline block (bulk single-cluster spectrum, not
  the calibrated inner-ring sample).
- A `# gate_efficiency_pct` comment column may be added once per
  pre-window as `100 × n_events(gated) / n_events(legacy)` for the
  user-requested "gate efficiency" line item.

---

## 6. Trend plots

Path: `analysis/pyscripts/output/timing_gate_benchmark_24386_curves.png`

Single PNG, 3 vertically-stacked panels (one `plt.subplots(3, 1,
figsize=(8, 11), sharex=True)`), x-axis = pre-window LO (ns).

### Panel A — σ_E/E vs. LO
- Two curves: `legacy` (grey), `gated W=16` (blue).
- Both the Gaussian (`σ_over_E_G`, dashed) and Crystal-Ball
  (`σ_over_E_CB`, solid) resolutions overlaid.
- Horizontal line at the technical-note headline (1.90 % CB) for
  reference.
- Vertical line at `LO = 130` labelled "runinfo baseline".
- y-axis as percent.

### Panel B — Event acceptance vs. LO
- Two curves: `legacy` (grey), `gated` (blue), both plotted as
  `n_events / n_events_at_LO=50` × 100 %.
- Horizontal dashed line at 90 % (success-criterion threshold).
- Vertical line at `LO = 130`.

### Panel C — Energy peak position vs. LO
- Four curves: `legacy μ_pre`, `legacy μ_post`, `gated μ_pre`,
  `gated μ_post`.
- Horizontal dashed line at `--target-peak = 3488.5 MeV`.
- Vertical line at `LO = 130`.
- y-axis in MeV, zoomed to ±200 MeV around target.

Figure title encodes the fixed-knob context so the plot stands alone:
`"Run 24386 — pre-window LO sweep (HI=200 ns, W=16 ns, iter1 calib, split 00050, 50k events)"`.

---

## 7. Markdown summary deliverable

Path: `analysis/pyscripts/output/timing_gate_benchmark_24386_findings.md`

Sections (~1 page):
1. **Setup** — one-paragraph recap of §4 fixed knobs + the swept grid.
2. **Headline result** — recommended pre-window, with the σ_E/E and
   acceptance at that choice vs. the runinfo baseline `[130, 200]`.
3. **Degradation slope** — "each 25 ns reduction in width costs X pp
   in σ_E/E" (linear fit over the monotonic degradation region, i.e.
   excluding the widest end if it has saturated).
4. **Gate behaviour** — one-row comparison of legacy vs. gated σ_E/E
   at the recommended pre-window, confirming the two-decimal agreement
   expected from the technical note on this clean run.
5. **Energy peak stability** — comment on whether `μ_post` is stable
   across pre-windows (expected: yes, since calibration refits
   per-window) and whether `μ_pre` tracks the pre-window change
   (expected: yes, this is the knob's direct effect).
6. **Caveats** — single split, 50k events; no explicit cluster-time
   observable (see §2.5); no cross-run generalisation claim.
7. **Attachments** — relative links to the aggregate TSV, curves PNG,
   and the per-window output directory.

---

## 8. Resource requirements & runtime estimate

### Compute

- Python environment with `prad2py` installed and
  `PRAD2_DATABASE_DIR` set (`analysis/pyscripts/setup_env.sh`).
- `scipy` for the fits (already a benchmark-script requirement).
- One CPU core; the benchmark is sequential and I/O-bound on the EVIO
  reader. No GPU, no distributed compute.
- ~500 MB free RAM (event lists for the calibration branch — scales
  with 50 000 events × 2 paths × ~56 B).

### Disk

- Per invocation: ~500 kB (`.tsv` ≈ 10 kB, `.png` ≈ 150 kB,
  `_calibrated.png` ≈ 250 kB, `.log` ≈ 50 kB).
- 7 invocations × 500 kB ≈ 3.5 MB in
  `output/timing_gate_benchmark_24386/`.
- Aggregate TSV + curves PNG: < 500 kB.

### Runtime

- Prior iter1 test: 50 000 events in 33.4 s wall on split .00050.
- Each sweep run runs the identical pipeline; wall time per run is
  dominated by waveform analysis and is **independent of pre-window**
  (pre-window is applied after `wa.analyze`, see
  `benchmark_hycal_timing.py` lines 399–415).
- Expected total wall: **7 × 35 s ≈ 4 min** for the sweep + ~5 s for
  the aggregate step + plotting = **under 5 minutes total**.
- If the user opts for the per-split-×-11 extension (§9.1), budget
  **~45 min**. If they opt for the full `[LO × HI]` 2-D scan (§9.2),
  budget scales with the grid.

---

## 9. Extensions (out of initial scope; propose only if user requests)

1. **Per-split cross-check.** Run the sweep on 3 representative splits
   (e.g. .00010, .00050, .00100) and overlay all three on the trend
   plots to validate intra-run stability. +~10 min.
2. **2-D `LO × HI` scan.** Add a second axis sweeping `HI ∈ {175, 200,
   225, 250}`. Replace trend lines with 2-D heatmaps. +~15 min
   compute, +~1 day plot-design work.
3. **Cluster-time observable.** Patch `benchmark_hycal_timing.py` to
   emit mean cluster time per path per event into the TSV; re-run the
   sweep. Adds a true "timing bias vs. pre-window" curve. Requires a
   separate plan for the script edit.
4. **Full-run baseline.** Repeat the recommended optimum on the full
   11-file set with `--max-events 0` and reconfirm σ_E/E — gives the
   publication-ready number vs. the technical note's existing full-run
   result. +~6 min.

---

## 10. Dependencies & risks

- **Dependency: `benchmark_hycal_timing.py` TSV schema is stable.** If
  the script is edited to add/remove calibration columns, the driver's
  parser breaks silently. Mitigation: driver asserts on exact column
  headers and fails loudly on mismatch.
- **Dependency: iter1 calibration file is final.** Prompt states "no
  further iteration needed"; this is user-asserted and taken on
  faith. If a new iter is published mid-study, the sweep must be
  re-run end-to-end.
- **Risk: narrowest window (LO=175) may have too few per-module
  entries** to converge the inner-ring calibration (threshold:
  `--min-per-module 30`). If the calibration block shows
  `n_modules_used ≪ 141`, the σ_E/E for that row is apples-to-oranges
  vs. the wider windows. Mitigation: driver flags any row whose
  `n_modules_used` drops >10 % below the widest-window value and the
  summary document notes the caveat.
- **Risk: Crystal-Ball fit degenerates on low-statistics rows.** The
  script already reports NaN when the fit fails; the driver must
  propagate NaN into the aggregate (not treat as 0). Trend plot skips
  NaN points.
- **Risk: `discover_split_files` ambiguity on single-split input.**
  `_common.py` line 96–97 confirms a non-glob, non-directory path is
  returned verbatim; the single-split invocation is safe.

---

## 11. Deliverables inventory

| artifact | path |
|---|---|
| aggregate metrics TSV | `analysis/pyscripts/output/timing_gate_benchmark_24386_summary.tsv` |
| trend plots PNG | `analysis/pyscripts/output/timing_gate_benchmark_24386_curves.png` |
| per-window outputs dir | `analysis/pyscripts/output/timing_gate_benchmark_24386/` |
| ↳ per-window TSVs | `pre_{LO}_{HI}.tsv` × 7 |
| ↳ per-window headline PNGs | `pre_{LO}_{HI}.png` × 7 |
| ↳ per-window calibration PNGs | `pre_{LO}_{HI}_calibrated.png` × 7 |
| ↳ per-window logs | `pre_{LO}_{HI}.log` × 7 |
| findings markdown | `analysis/pyscripts/output/timing_gate_benchmark_24386_findings.md` |
| sweep driver script | `analysis/pyscripts/sweep_pre_window_24386.py` (new, standalone) |

---

## 12. Handoff after approval

- Hand the sweep-driver implementation to **`script-implementer`** —
  it is a standalone Python script with no library edits.
- Hand the findings-markdown writeup to **`doc-writer`** once the
  driver has produced the aggregate TSV and curves, since the headline
  numbers are not knowable until the sweep completes.
- Delegate TodoWrite management and overall sequencing to
  **`supervisor`**.

---

## Todo

### Phase A — Pre-flight (manual verification; do NOT skip)

- [ ] A1. Confirm `prad2py` is importable:
  `python -c "from prad2py import dec, det"` from
  `analysis/pyscripts/` with `setup_env.sh` sourced.
- [ ] A2. Confirm the iter1 calibration file is readable and
  non-empty:
  `test -s /Users/L00338853/work/PRad/data/physics_calib_out/Physics_calib/run24386/calib_factor_iter1.json`.
- [ ] A3. Confirm the target EVIO split exists and is readable:
  `test -r /Users/L00338853/work/PRad/data/evio/prad_024386.evio.00050`.
- [ ] A4. Re-run the prior iter1 smoke test and confirm numbers match
  `baseline_024386_iter1_test.tsv` (`n_phys ≈ 29923`, `pre_window
  [100,200]`, legacy `n_events ≈ 20383`). Catches environment drift
  before the sweep burns 5 min.
- [ ] A5. User resolves §2 Open Questions (sweep grid, HI-fixed vs.
  swept, split selection, `--window` value, timing-bias definition,
  calibration-file wiring).

### Phase B — Driver implementation (hand to `script-implementer`)

- [ ] B1. Create `analysis/pyscripts/sweep_pre_window_24386.py`:
  - `argparse` with: `--pre-window-los` (nargs='+', default
    `[50,75,100,125,130,150,175]`), `--pre-hi` (default `200`),
    `--window` (default `16`), `--evio-path` (default
    `/Users/L00338853/work/PRad/data/evio/prad_024386.evio.00050`),
    `--max-events` (default `50000`), `--hc-calib-file` (default
    iter1 path), `--out-dir` (default
    `analysis/pyscripts/output/timing_gate_benchmark_24386`),
    `--summary-tsv` (default
    `analysis/pyscripts/output/timing_gate_benchmark_24386_summary.tsv`),
    `--curves-png` (default
    `analysis/pyscripts/output/timing_gate_benchmark_24386_curves.png`),
    `--skip-existing` (bool, re-use prior per-window outputs if
    `.tsv` is present).
  - Build the per-window command list using §4.3 as the template;
    pin `--signal-min 500`, `--calibrate`, `--cal-method cb-peak`,
    `--exclude-rowcol 15 20 15 20`, `--target-peak 3488.5`,
    `--bin-width 20`, `--e-range 2500 4000`,
    `--inner-radius 200`, `--min-seed-energy 200`,
    `--min-per-module 30` as hardcoded constants with a comment
    anchoring them to the hycal_clustering.md reference recipe.
- [ ] B2. Implement `run_one(lo, hi, cfg) -> Path` that `subprocess.run`s
  the benchmark with `check=True`, redirects stdout+stderr to the
  per-window `.log`, honours `--skip-existing`.
- [ ] B3. Implement `parse_tsv(tsv_path) -> dict` that extracts the
  headline block (per path: `n_events`, `gauss_mu_MeV`,
  `gauss_sigma_MeV`, `n_in_fit`) and the calibration block (per
  path: `n_modules`, `n_inner_events`, `mu_pre_MeV`, `sigma_pre_MeV`,
  `mu_post_MeV`, `sigma_post_MeV`, `sigma_over_E_post`, `mu_cb_MeV`,
  `sigma_cb_MeV`, `alpha_cb`, `n_cb`, `sigma_over_E_cb`).
  - Assert both blocks are present; raise `RuntimeError` with the
    path if either is missing (fails loud on `--calibrate` omission).
  - Handle NaN values (string `"nan"`) gracefully; coerce to
    `float('nan')`.
- [ ] B4. Implement `write_summary(rows, out_tsv)` that emits the §5
  schema, 2 rows per LO (legacy + gated).
- [ ] B5. Implement `plot_curves(summary_tsv, out_png)` that reads
  the aggregate and emits the §6 three-panel figure. Use
  `matplotlib.pyplot` only (no seaborn; matches existing scripts).
  Mark `LO = 130` with a dotted vertical line labelled "runinfo
  baseline"; mark the 90 % acceptance line on panel B; mark the
  `--target-peak = 3488.5 MeV` line on panel C.
- [ ] B6. `__main__` orchestration: pre-flight (assert
  `--hc-calib-file` and `--evio-path` exist), create `--out-dir`,
  loop `run_one` for each LO (log progress), then
  `write_summary`, then `plot_curves`.
- [ ] B7. Add a module docstring summarising inputs, outputs, and
  usage example mirroring the §4.3 invocation (minus the
  per-pre-window loop; refer to this plan for the sweep grid).

### Phase C — Execution (hand to `supervisor`; monitor progress)

- [ ] C1. Create output directory:
  `mkdir -p analysis/pyscripts/output/timing_gate_benchmark_24386`.
- [ ] C2. Dry-run the driver with `--pre-window-los 100`
  (single-point) to catch wiring bugs cheaply (~35 s wall). Confirm
  the per-window TSV + both PNGs land in the output dir and the
  aggregate contains 2 rows.
- [ ] C3. Full sweep: run the driver with the approved grid (default
  7 LOs). ETA ~5 min. Capture overall log to
  `output/timing_gate_benchmark_24386_sweep.log`.
- [ ] C4. Smoke-check the aggregate TSV: `wc -l` ≈ 15 (14 data rows
  + header), no `nan` in `sigma_over_E_CB` for the four widest
  windows (LO ≤ 125).
- [ ] C5. Visually inspect `timing_gate_benchmark_24386_curves.png`:
  - Panel A: σ_E/E CB monotonically rises as LO narrows past
    runinfo baseline.
  - Panel B: acceptance plateau on the wide end, modest fall-off on
    the narrow end.
  - Panel C: `μ_post` stable (within ±10 MeV of 3488.5),
    `μ_pre` trends.
- [ ] C6. Flag any row where `n_modules_used` drops >10 % below the
  LO=50 reference (narrow-window calibration under-converged) and
  record in the findings doc caveats.

### Phase D — Findings writeup (hand to `doc-writer` after Phase C)

- [ ] D1. Author `timing_gate_benchmark_24386_findings.md` per §7
  outline. Use actual numbers from the sweep aggregate.
- [ ] D2. State the recommended pre-window and justify against the
  success criteria (σ_E/E minimised ∧ acceptance > 90 %).
- [ ] D3. Fit a linear trend to the degradation region of σ_E/E vs.
  width; report slope in `pp per 25 ns`.
- [ ] D4. Note caveats from §10 that materialised in execution.
- [ ] D5. Link all §11 deliverables with relative paths.

### Phase E — Optional follow-ups (user-directed)

- [ ] E1. If user requests the §9.1 per-split cross-check, re-run
  the driver with `--evio-path` set to each of .00010, .00050,
  .00100 into separate output subdirs; overlay all three on an
  updated curves plot.
- [ ] E2. If user requests the §9.3 cluster-time observable, open a
  new plan to extend `benchmark_hycal_timing.py` with a mean-time
  TSV field and re-run.
- [ ] E3. If user requests the §9.4 full-run baseline at the chosen
  optimum, invoke the benchmark directly (bypassing the sweep
  driver) with `--max-events 0` on the EVIO directory.
