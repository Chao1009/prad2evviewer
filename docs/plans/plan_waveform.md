# Plan — waveform-analysis arc (living summary)

**Status:** living document, revised as each phase lands.
**Scope:** the overall PRad-II waveform-analysis workstream — template
extraction, matched filter, deconvolution, CFD timing, pile-up studies,
and the clustering timing-gate validation that depends on all of the
above. This file is the single roll-up; per-piece plans live in their
own `plan-*.md` files.

---

## 1. Purpose

Track the state of the PRad-II HyCal waveform-analysis arc across its
components (template extraction, matched filter, per-pulse LM template
fit "deconvolution", CFD timing, synthetic pile-up studies, clustering
timing-gate validation on real data) and sequence the open work. The
document is updated whenever a component lands, a measurement result
comes in, or scope shifts between phases.

Three questions this arc exists to answer:

1. Does the per-pulse LM template fit (`WaveAnalyzer::Deconvolve`),
   turned on in production, actually recover better amplitudes than
   the raw peak integral when real pile-up is present?
2. Does CFD timing at `cfd_fraction = 0.5` eliminate the amplitude-
   dependent timing walk well enough that a seed-anchored
   coincidence cut improves the single-cluster energy resolution,
   rather than making it worse?
3. Where is the right pulse-detection trigger — the current
   `findPeaks` local-maxima path, or a matched-filter-based detector —
   and does the choice matter for the two questions above?

## 2. Current state (as of today)

**In place and working.**

- **Pulse-template extraction** from physics data via
  `analysis/pyscripts/fit_pulse_template.py` plus the trigger-scan
  grid driver `scripts/run_trigger_scan.sh`. Per-type (PbGlass /
  PbWO4 / Veto) `(τ_r, τ_f, t0)` + per-channel entries in the JSON
  schema consumed by `fdec::PulseTemplateStore::LoadFromFile`.
- **Template comparison across runs / triggers / cuts** via
  `analysis/pyscripts/compare_templates.py`.
- **Per-type pulse-template JSON files** (several, under
  `analysis/pyscripts/output/`); current best-match for Phase 1 is
  `pulse_templates_025308_h50_RawSum_max20modules.json`.
- **Synthetic pile-up generator** with `(ΔT, ratio)` grid in
  `analysis/pyscripts/pileup_generator.py` (plan:
  `docs/plans/plan-pileup-generator.md`).
- **Matched-filter prototype in Python** at
  `analysis/pyscripts/benchmark_matched_filter.py:205, 246, 276`.
  Operates on pedsub waveforms using `np.correlate` + local-maxima +
  prominence filter + greedy NMS. **Standalone benchmark only — NOT
  wired into the C++ production pipeline.**
- **Per-pulse LM template fit ("deconvolver")** implemented at
  `prad2dec/src/WaveAnalyzer.cpp::Deconvolve` and bound as
  `wa.deconvolve(samples, wres, template)` in
  `python/bind_dec.cpp:559–582`. The explicit `wa.deconvolve` call
  bypasses the production `nnls_deconv.enabled` config gate
  (`WaveAnalyzer.cpp:832–835`). **Currently disabled in production**
  (`database/daq_config.json:117` has `"enabled": false`).
- **Benchmark of deconv on synthetic pile-up** via
  `analysis/pyscripts/benchmark_deconv.py`.
- **Standard CFD timing** at `cfd_fraction = 0.5f` (hardcoded,
  `WaveAnalyzer.cpp:523`) baked into the private
  `WaveAnalyzer::findPeaks` at `WaveAnalyzer.cpp:518–562`. Applied
  unconditionally to every peak that `findPeaks` discovers; result
  stored as `peak.time` on the returned `WaveResult.peaks[]`. Not
  callable in isolation from Python.
- **Clustering timing-gate validation on real data** via
  `analysis/pyscripts/benchmark_hycal_timing.py` — legacy vs. gated
  `HyCalCluster` with per-event σ_E/E Gaussian + CB fits.
  **Current baseline (run 25320, 2.2 GeV e-p, no template fit):**
  gated `W = 4 ns` is **~14% worse** than legacy in σ_E/E.
  Hypothesised cause: amplitude-dependent timing walk on the raw
  `findPeaks` CFD time is clipping the timing-cut tails of low-
  amplitude neighbour pulses out of valid clusters, leaving a less
  energetic — and worse-measured — cluster for the one-cluster
  spectrum.
- **Per-module Δt diagnostic** via
  `analysis/pyscripts/study_hycal_timing.py` (TSV producer) plus
  `analysis/pyscripts/module_dt_diagnostic.py` (fit + flag;
  plan: `docs/plans/plan-module-dt-diagnostic.md`). Finding:
  primary ΔT peak at μ ≈ 0 with σ ≈ 1 ns across most HyCal
  modules; secondary ΔT peak at ΔT ≈ 15 ns on a fraction of
  modules diagnosed as PMT afterpulsing. Documented in
  `docs/analysis_notes/hycal_time_calib.md` §10.

## 3. Open question / motivation for Phase 1

Does swapping the raw-peak amplitude path for the template-fit
amplitude path (`wa.deconvolve`, explicit call bypassing the
production gate), while keeping the standard CFD timing
(`wres.peaks[k].time`) for the coincidence cut and introducing a
smallest-|Δt| neighbour-selection policy on the Python side,
restore the expected `σ_E/E(gated) ≤ σ_E/E(legacy)` ordering that
the baseline currently inverts?

Phase 1 answers this empirically on the same dataset
(`prad_025320.evio.*`) with a 2×2 config grid.

## 4. Phase 1 — immediate

**Deliverable.** Build
`analysis/pyscripts/benchmark_hycal_timing_deconv.py` per the
companion plan `docs/plans/plan-benchmark-hycal-timing-deconv.md`.

**Compares.** Four configurations on the same event stream:

- `legacy_largestE`, `legacy_smallestDt`,
- `gated_largestE`, `gated_smallestDt`.

**Scope boundary.** Python-only. No C++ changes, no binding changes.
Matched filter is NOT in the Phase 1 loop (candidate detection stays
on `WaveAnalyzer::findPeaks` via `wa.analyze_result`). The scope
reduction is deliberate and documented — the literal matched-filter
→ template-fit → CFD pipeline needs the three binding/C++ changes
listed under §5, none of which are in place today, and the
measurement in §3 is cleanly answerable without them.

**Expected outcome branches.**

Let Δ ≡ `σ_E/E(gated-template-fit) − σ_E/E(legacy-template-fit)`
(signed, in percentage points). Let σ(Δ) be the quadrature sum of
the two fits' statistical uncertainties on σ_E/E, read from the
Phase 1 Gaussian fit covariance.

- **Branch A — gated template-fit matches or beats legacy
  template-fit in σ_E/E** (Δ ≤ 0, or |Δ| < 0.3 pp regardless of
  sign). Template fit was the missing (or at worst neutral)
  ingredient. Phase 3 recommendation: enable
  `nnls_deconv.enabled = true` in production `daq_config.json`,
  keep `seed_time_window = 4 ns`. Consider `apply_to_all_peaks =
  true` for best timing-walk correction (approximate 10×
  slowdown — acceptable for offline, likely not for online).
  The 0.3 pp band is the "indistinguishable from a tie"
  threshold; a marginal result inside the band is interpreted as
  "template fit was sufficient" even if the sign nominally
  favours legacy.
- **Branch B — gated template-fit still loses to legacy
  template-fit** (Δ ≥ 0.3 pp AND Δ > 3·σ(Δ)). Template fit alone
  is insufficient. Investigate the gated clusterer's pulse-
  selection policy using the 2×2 grid outcome; consider widening
  `seed_time_window` to 6–8 ns; move to Phase 2 (matched-filter in
  the detection loop). Branch B requires BOTH the practical-
  significance threshold (0.3 pp) AND the statistical-
  significance threshold (3σ); either alone is insufficient to
  justify Phase 2's binding + C++ work.
- **Branch C — smallest-Δt consistently beats largest-E** (either
  in the legacy or the gated row of the 2×2 grid, or both; same
  0.3 pp / 3σ joint threshold applied to the within-row
  difference). Strong argument for changing the C++ cluster BFS
  pulse-selection policy (Phase 3 item at
  `prad2det/src/HyCalCluster.cpp:197–209`).

Branches are not exclusive: a run can land in A + C simultaneously
(both the template-fit fixed the gated regression AND smallest-Δt
beats largest-E), in which case Phase 3 opens two items.

## 5. Phase 2 — deferred, conditional on Phase 1 Branch B

Enable the literal matched-filter → template-fit → CFD pipeline:
use `benchmark_matched_filter.py:find_mf_peaks(...)` to produce the
candidate peak set, feed those candidates into
`WaveAnalyzer::Deconvolve` for amplitude, and recompute CFD timing
via a newly-exposed public C++ helper. Required work:

The missing pieces for Phase 2 are three — the `Peak` constructor is
already available (`python/bind_dec.cpp:204–217` binds `fdec::Peak`
with `.def(py::init<>())`, so a Phase 2 script can construct Peak
objects directly):

- **Binding changes** (`python/bind_dec.cpp`).
  - **(a)** Make `WaveResult.npeaks` writable (currently
    `.def_readonly("npeaks", ...)` at `bind_dec.cpp:366`). Needed
    so Python can set the active peak count on a `WaveResult` it
    has populated from matched-filter candidates.
  - **(b)** Add an `assign_peaks(list[Peak])` method (or
    equivalent) to `WaveResult` so a Python-constructed peak list
    can be safely injected into the fixed-size internal
    `peaks[MAX_PEAKS]` buffer without touching other `WaveResult`
    state (`ped`, baseline fit, etc.). A `set_peak(i, Peak)`
    helper is a leaner alternative but requires the caller to
    manage `npeaks` separately.
- **C++ changes.**
  - **(c)** `prad2dec/include/WaveAnalyzer.h` +
    `prad2dec/src/WaveAnalyzer.cpp`: extract the CFD block (currently
    at `WaveAnalyzer.cpp:518–562`, inside the private `findPeaks`)
    into a public method `compute_cfd_time(samples, nsamples,
    raw_pos, left, right, ped_mean, ped_rms, clk_mhz, cfd_fraction =
    0.5)` returning the interpolated peak time in ns.
  - Bind the new helper (`python/bind_dec.cpp`) as
    `WaveAnalyzer.compute_cfd_time(...)`.
- **Script changes.**
  - In `benchmark_hycal_timing_deconv.py` (or a Phase-2 sibling), swap
    the Phase 1 `wa.analyze_result` candidate detection for
    `find_mf_peaks(...)` from `benchmark_matched_filter.py`.
  - Construct a `WaveResult` from MF candidates (using the new
    binding helpers), feed to `wa.deconvolve`, compute CFD time per
    peak via `wa.compute_cfd_time(...)`.
  - Add a `--detector {findpeaks,matched-filter}` CLI switch so the
    same script covers both Phase 1 and Phase 2 modes.

**Estimated effort.** 1–2 days of C++ + binding work; 1 day for
the script refactor; 1 day for review and benchmarking.

**Precondition.** Phase 1 must first demonstrate that template-fit
alone is insufficient (Branch B).

## 6. Phase 3 — further deferred, contingent on Phase 1 / 2 outcomes

Each item here is one planner artifact (`plan-*.md` file) to be
drafted when the triggering condition is met.

- **Enable `nnls_deconv.enabled = true` (and optionally
  `apply_to_all_peaks = true`) in production
  `database/daq_config.json`.** Trigger: Phase 1 Branch A. Requires
  a short validation plan plus the config edit itself and a note in
  the technical-notes directory.
- **Change the C++ cluster BFS to smallest-|Δt| pulse selection at
  `prad2det/src/HyCalCluster.cpp:197–209`.** Trigger: Phase 1 Branch
  C sustained across more than one run / trigger mix. Needs a
  dedicated plan and a C++ change; probably introduce a runtime
  switch rather than a hard replacement so the existing largest-E
  policy stays reachable.
- **V1–V4 veto counter exclusion in
  `benchmark_hycal_timing_deconv.py`.** Not needed in Phase 1
  because `mod.is_hycal()` already drops Veto modules
  (`HyCalSystem.h:174`). If a downstream track-matching stage needs
  V-counter inhibit windows, that is a separate script. Trigger:
  explicit request from the track-matching workstream.
- **Technical note** on the amplitude-dependent timing walk
  observation and the template-fit improvement, likely under
  `docs/technical_notes/waveform_analysis/`. Trigger: completion of
  Phase 1 with a result (either branch) worth writing up.

## 7. References

- `docs/plans/plan-benchmark-hycal-timing-deconv.md` — Phase 1
  implementation plan.
- `docs/plans/plan-module-dt-diagnostic.md` — per-module ΔT
  diagnostic plan.
- `docs/plans/plan-pileup-generator.md` — synthetic pile-up
  generator plan.
- `docs/research-study-hycal-timing-time-offset.md` — end-to-end
  trace of the raw `pk.time` path through `study_hycal_timing.py`
  (confirming no `Module::time_offset` subtraction on the Python
  side).
- `docs/analysis_notes/hycal_time_calib.md` §10 — secondary ΔT
  coincidence peak diagnosed as PMT afterpulsing.
- `docs/technical_notes/waveform_analysis/wave_analysis.md` —
  waveform-algorithm technical notes.
- `analysis/pyscripts/Commands.md` — running log of commands.

## 8. Change log

| Date       | Phase     | Status                              | Notes                                                                                                                                   |
|------------|-----------|-------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------|
| 2026-10-13 | Phase 1   | planned, awaiting implementation    | `plan-benchmark-hycal-timing-deconv.md` drafted. 2×2 grid on run 25320, template file `pulse_templates_025308_h50_RawSum_max20modules.json`. |
| _future_   | Phase 1   | implemented / in validation         | record the four `sigma_over_E_pct` values and which outcome branch (A/B/C) fired.                                                        |
| _future_   | Phase 2/3 | opened / closed                     | one row per plan artifact opened or closed as a direct consequence of the Phase 1 result.                                                |
