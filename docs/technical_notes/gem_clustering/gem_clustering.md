# GEM Clustering in `prad2det`

**Author:** Chao Peng (Argonne National Laboratory)

`gem::GemCluster` (in [`prad2det/include/GemCluster.h`](../../../prad2det/include/GemCluster.h))
takes per-plane strip charges from the SSP/MPD readout and produces 2-D
GEM hits. It runs in two stages: a **per-plane 1-D clustering**
pipeline (group consecutive strips → recursively split at charge
valleys → charge-weighted position → cross-talk filter), followed by
an **X/Y matching** pass that turns paired X- and Y-clusters into
3-vector `GEMHit`s.

The algorithm is ported from `mpd_gem_view_ssp::GEMCluster` (Bai /
Gnanvo / Peng) and is invoked once per detector per event from
`gem::GemSystem::Reconstruct()`.

## GEM detector layout

PRad-II uses **four GEM detectors** in **two paired layers** for
redundant 2-D tracking and ghost-hit rejection. Each detector reads
out two orthogonal strip planes through APV25 chips at 0.4 mm pitch:

![layout](plots/gem_fig1_layout.png)

| field | value |
|---|---|
| detectors | GEM0, GEM1 (layer 1, z ≈ 5407 mm); GEM2, GEM3 (layer 2, z ≈ 5807 mm) |
| inter-detector spacing | 39.7 mm (within a layer) |
| X plane | 12 APVs × 128 ch, pitch 0.4 mm — 1 408 strips ≈ 563.2 mm. APV positions 10/11 share strips around the hole (`shared_pos: 10`, `pin_rotate: 16`) |
| Y plane | 24 APVs × 128 ch, pitch 0.4 mm — 3 072 strips ≈ 1 228.8 mm. Strips that cross the hole y-band are split into top + bottom segments |
| beam hole | 52 × 52 mm², centred vertically (y = 614.4 mm) and offset along x (x = 534.4 mm) — sits inside the APV pos 10/11 strip range |

Strip indices come from `GemSystem::ProcessEvent()` after pedestal
subtraction, common-mode correction, and zero-suppression. Each
`StripHit` carries the plane-wise strip number, the maximum
time-sample charge (`charge`), the time bin where that maximum lives
(`max_timebin`), the physical position in mm, and the full
time-sample ADC vector (used downstream for time-difference cuts).

## Algorithm

### Step 1 — Group consecutive strips

`groupHits()` sorts the input `StripHit`s by strip number, then walks
the sorted list and starts a new group whenever the gap to the next
strip exceeds `consecutive_thres` (default 1, i.e. only strictly
adjacent strips group). Each group is a candidate 1-D cluster but may
still need to be split if multiple showers share neighbouring strips.
Optional SBS-style strip cuts (off by default) can end a group at a
failing strip — see
[SBS-style quality variables and cuts](#sbs-style-quality-variables-and-cuts).

### Step 2 — Recursive valley split

For each group `splitCluster()` looks for an internal local minimum
preceded by a sufficiently steep descent and followed by a
sufficiently steep ascent — both gated by `split_thres` (default 14
ADC counts):

1. Walk the group left-to-right. Set `descending = true` as soon as
   `charge[i] − charge[i+1] > split_thres`.
2. Once descending, track the running minimum.
3. The first ascent step where `charge[i+1] − charge[i] > split_thres`
   confirms the valley.
4. Halve the charge of the valley strip (it is shared between the
   left and right sub-clusters), emit the left sub-cluster, and
   recurse on the right.

Groups smaller than 3 strips never split. Groups with multiple
valleys split recursively, producing one cluster per "peak".

### Step 3 — Charge-weighted position

Each cluster's position comes from `reconstructCluster()`:

```
position    = Σ (strip_pos_i · charge_i) / Σ charge_i      [mm]
peak_charge = max_i charge_i
total_charge = Σ charge_i
max_timebin  = time bin of the seed (highest-charge) strip
```

This is a plain centroid, not log-weighted — at 0.4 mm pitch the
strip density is high enough that linear weighting tracks the shower
position well, and unlike HyCal there's no module-scale grid to
worry about.

### Step 4 — Cross-talk identification

APV25 capacitive coupling and electronic cross-talk produce **ghost
clusters** at characteristic distances from a true cluster (a known
function of the chip's pin-out and sampling pattern). `setCrossTalk()`
sorts clusters by ascending `peak_charge` and walks the list: any
cluster all of whose strips are flagged `cross_talk` (set upstream by
`GemSystem` based on inter-strip ADC ratios) and whose distance to a
larger cluster matches one of the `charac_dists` entries within
±`cross_talk_width` mm gets its own `cross_talk` flag set.

The default characteristic distances mirror the original
`mpd_gem_view_ssp` values:

```cpp
charac_dists = { 6.4, 17.6, 24.4, 24.8, 25.2, 25.6,
                 26.0, 26.4, 26.8, 33.6, 44.8 };   // mm
```

### Step 5 — Filter

`filterClusters()` drops:

- clusters with fewer than `min_cluster_hits` (default 1) strips,
- clusters with more than `max_cluster_hits` (default 20) strips
  (these are typically noise bursts, not real showers),
- clusters flagged as cross-talk,
- clusters failing one of the optional SBS-style quality cuts (all off
  by default, see
  [below](#sbs-style-quality-variables-and-cuts)).

### Step 6 — Cartesian X/Y matching

`CartesianReconstruct()` turns a list of accepted X-clusters and a
list of accepted Y-clusters into 2-D `GEMHit`s. Two modes are
supported:

**Mode 0 — ADC-sorted 1:1.** Sort both lists by descending
`peak_charge` and pair index-wise. Always produces
`min(N_x, N_y)` hits with no rejected combinations. Useful when the
upstream strip clustering is clean enough that the brightest X
genuinely matches the brightest Y.

**Mode 1 — full Cartesian + cuts (default).** Form every X×Y pair,
then drop any pair that fails either of:

- **ADC asymmetry** — `|Q_X_peak − Q_Y_peak| / (Q_X_peak + Q_Y_peak) ≤ match_adc_asymmetry` (default 0.8). Real GEM hits deposit similar charge on both planes; large asymmetries flag ghosts from accidental coincidence.
- **Time difference** — the seed-strip ADC-weighted mean times must satisfy `|⟨t⟩_X − ⟨t⟩_Y| ≤ match_time_diff` (default 50 ns). The seed mean time is `Σ(adc_i · t_i) / Σ adc_i` over the time samples with `adc_i > 0`, where `t_i = (i + 1) · ts_period` ns (`StripMeanTime`, see [SBS-style quality variables](#sbs-style-quality-variables-and-cuts)). A pair is not rejected by this cut when either seed has no positive sample (undefined time).

Any pair passing both cuts becomes a `GEMHit` with `det_id` and
per-plane charge / size / max-timebin recorded for downstream use.

## Parameters

All settings live in `gem::ClusterConfig` (in
[`prad2det/include/GemSystem.h`](../../../prad2det/include/GemSystem.h)).
The system stores one `ClusterConfig` per detector via
`SetReconConfigs()`, so different detectors can have different cuts —
useful when the four GEMs have different APV gains or noise levels.

| field | default | unit | role |
|---|---:|---|---|
| `min_cluster_hits` | 1 | strips | Lower bound on cluster size — drops single-noisy-strip "clusters" if set ≥ 2. |
| `max_cluster_hits` | 20 | strips | Upper bound on cluster size — kills runaway clusters caused by noise bursts or beam halo. |
| `consecutive_thres` | 1 | strips | Maximum gap between adjacent strips that still keep them in the same group. `1` = strictly adjacent only. |
| `split_thres` | 14 | ADC counts | Charge-difference threshold that gates both descent detection and valley-confirmation in `splitCluster()`. Lower values split more aggressively. |
| `cross_talk_width` | 2 | mm | Tolerance on the characteristic-distance match in `setCrossTalk()`. |
| `charac_dists` | `{6.4, 17.6, 24.4..26.8, 33.6, 44.8}` | mm | APV25 cross-talk characteristic distances. |
| `match_mode` | 1 | — | `0` = ADC-sorted 1:1 matching, `1` = Cartesian product with cuts. |
| `match_adc_asymmetry` | 0.8 | fraction | Cap on `|Q_X − Q_Y|/(Q_X + Q_Y)` (mode 1). Set negative to disable. |
| `match_time_diff` | 50 | ns | Cap on `|⟨t⟩_X − ⟨t⟩_Y|` (mode 1). Set negative to disable. |
| `ts_period` | 25 | ns | Time-sample period (default = 1 / 40 MHz APV clock). Also used for the strip mean times of the quality variables. |
| `strip_time_min/max`, `strip_unimodal`, `seed_min_peak_adc`, `seed_min_sum_adc`, `strip_time_agreement`, `strip_ts_corr_min` | off | | SBS-style quality cuts, see [below](#sbs-style-quality-variables-and-cuts). |

## Worked example — strip clustering

A 75-strip window with three real showers and ~4 ADC noise:

![strip](plots/gem_fig2_strip_clustering.png)

The left panel shows the full above-threshold strip distribution
(threshold 30 ADC). DFS-style grouping with `consecutive_thres = 1`
produces three groups (110-114, 130-139, 158-160). The middle group
has two local maxima — `splitCluster()` walks from the 132 peak,
detects a descent into the 134/135 valley, then sees a 198 ADC upturn
into the 137 peak (»`split_thres = 14`) and partitions the group at
strip 135 (right panel). The valley strip's charge is halved and
contributes to both sub-clusters, so neither side double-counts the
shared edge.

| cluster | strips | position (mm) | Σ ADC | peak ADC |
|---|:---:|---:|---:|---:|
| 1 | 110-114 | 44.87 | 2143 | 778 |
| 2 (left of split) | 130-134 | 52.81 | 1750 | 599 |
| 2 (right of split) | 135-139 | 54.82 | 1182 | 480 |
| 3 | 158-160 | 63.60 |  566 | 297 |

Note that the position (charge-weighted centroid) lands between
strips at sub-pitch resolution — for cluster 1, x = 44.87 mm sits
between strips 112 (44.8 mm) and 113 (45.2 mm), reflecting the actual
shower offset within the strip pitch.

## Worked example — X/Y matching

Three X-plane clusters and three Y-plane clusters per detector, with
two prompt big showers and one late, small "out-of-time" pair (e.g.
backsplash or a delayed accidental):

![matching](plots/gem_fig3_xy_matching.png)

**Left panel — Mode 1 (default).** All 9 X×Y candidates are listed.
The big-prompt × big-prompt pairings (X0/X1 with Y0/Y1) easily pass
both cuts. The big × late pairings fail the time cut (`Δt > 50 ns`).
Note that the small × small (X2 ↔ Y2) pair *does* pass — both have
similar (small) charge and similar (late) timing, so it gets
reconstructed as a 2-D hit even though it is most likely noise. If
that's a problem in production, tighter `match_adc_asymmetry` plus a
minimum on `peak_charge` upstream filters it out.

**Right panel — Mode 0.** No physical cuts; just sort both lists by
peak ADC and pair X[rank] ↔ Y[rank]. Three hits, by construction, but
with no defence against ghost-pair formation when accidental
coincidence is significant. Mode 0 is mostly useful for debugging or
for very low-occupancy runs where the cuts of mode 1 would just throw
away good hits.

## Parameter sensitivity

![params](plots/gem_fig4_params.png)

**Left — `split_thres`.** Same multi-peak group, three different
thresholds:

- `split_thres = 200`: too coarse — the two real peaks aren't
  resolved (1 cluster).
- `split_thres = 50`: catches the deep valley between the two main
  peaks (2 clusters), correct for this trace.
- `split_thres = 5`: also splits a small secondary fluctuation (3
  clusters), some of which are spurious.

The default of 14 ADC is calibrated for the typical noise RMS (~5
ADC) plus a margin to avoid splitting on noise; lower it on quieter
detectors, raise it on noisy ones.

**Right — cross-talk.** A bright primary cluster at 50 mm with a
small ghost cluster 24.4 mm away. The dotted lines show all 11
characteristic distances — the densely-clustered 24.4–26.8 mm group
covers the most common APV25 cross-talk pattern. `setCrossTalk()`
identifies the small cluster as a cross-talk match (it sits within
`cross_talk_width = 2 mm` of one of the characteristic distances and
has only `cross_talk`-flagged strips), and `filterClusters()` drops
it. The primary survives.

## Output — `GEMHit`

Each X/Y match produces one `GEMHit`:

| field | type | meaning |
|---|---|---|
| `x`, `y`, `z` | `float` | Hit position (mm). `z` is set by the application from per-detector geometry; `GemCluster` itself sets `z = 0`. |
| `det_id` | `int` | 0..3 (GEM0..GEM3). |
| `x_charge`, `y_charge` | `float` | Total ADC of the X / Y cluster. |
| `x_peak`, `y_peak` | `float` | Max-strip ADC of the X / Y cluster. |
| `x_max_timebin`, `y_max_timebin` | `short` | Time-sample bin of the max-ADC strip on each plane. |
| `x_size`, `y_size` | `int` | Number of strips in the X / Y cluster. |
| `x_time`, `y_time`, `time_diff`, `adc_asym`, `x/y_max_strip_dt`, `x/y_min_ts_corr` | `float` | SBS-style X/Y quality, see [below](#sbs-style-quality-variables-and-cuts). NaN = undefined. |

These map directly onto the `gem_*` branches of the recon tree (see
[`docs/REPLAYED_DATA.md`](../../REPLAYED_DATA.md)).

## SBS-style quality variables and cuts

The SBS GEM code (`mpd_gem_view_ssp`, `gem/src/Cuts.cpp` +
`GEMCluster.cpp`) judges strips, clusters and X/Y pairs with a set of
time-sample quality variables. `GemCluster` computes the same
variables for **every** cluster and 2-D hit, and offers the matching
cuts as `ClusterConfig` knobs. **All of these cuts are off by
default**, both in the library and in the shipped
`reconstruction_config.json` / `reconstruction_config_x17.json`, so
existing reconstruction output is unchanged (bit-for-bit) until a cut
is switched on. The variables are filled whether or not a cut is
enabled; `NaN` always means "undefined".

### Definitions

The three helpers are free functions in `GemCluster.h` (namespace
`gem`), so offline code can reuse them on raw `ts_adc` vectors:

- **Strip mean time** `StripMeanTime(ts_adc, ts_period = 25)` [ns]:

  ```
  t = Σ_{i: a_i > 0} a_i · (i + 1) · ts_period  /  Σ_{i: a_i > 0} a_i
  ```

  Only positive samples contribute. Sample 0 maps to `1 · ts_period`
  (SBS convention), so a 6-sample strip lies in [25, 150] ns. `NaN`
  if no sample is positive. This is exactly the seed mean time the
  existing X/Y time cut (`match_time_diff`) has always used.
- **Time-sample correlation** `TimeSampleCorrelation(a, b)`: Pearson
  `r = Σ(a−ā)(b−b̄) / sqrt(Σ(a−ā)² · Σ(b−b̄)²)`, accumulated in double.
  `NaN` if the sizes differ, size < 2, or either vector is flat.
- **Unimodal pulse** `IsUnimodalPulse(ts_adc)` (SBS
  `Cuts::is_concave_shape`): samples strictly rise up to the first
  maximum and strictly fall after it. Any flat pair or tail bump
  fails; a peak in the first or last sample passes (SBS behaviour
  since `cdb0e93`; edge peaks are the job of `reject_first/last_timebin`).
  Empty `ts_adc` fails.

**Seed** = the first strip with the maximum `charge` (strict `>`,
from the lowest strip up), as in SBS `__get_seed_strip_index` and the
X/Y time cut. Per cluster (`StripCluster`, filled in
`reconstructCluster()`):

| field | definition |
|---|---|
| `seed_time` | `StripMeanTime` of the seed strip (ns) |
| `seed_peak_adc` | max of the seed's `ts_adc` (SBS "seed strip peak ADC"). Equals `peak_charge` whenever strip `charge` = max sample, i.e. on every PRad path (a halved valley strip can never be the seed); kept for the SBS naming |
| `seed_sum_adc` | sum of all the seed's `ts_adc` samples, negatives included (SBS "seed strip sum ADC") |
| `max_strip_dt` | max over **all** non-seed strips of `\|t_i − seed_time\|` (ns) |
| `min_ts_corr` | min over all non-seed strips of `TimeSampleCorrelation(seed, strip_i)` |

Non-finite per-strip values are skipped when forming the max / min;
`max_strip_dt` and `min_ts_corr` are `NaN` for single-strip clusters.
Per 2-D hit (`GEMHit`, both match modes):

| field | definition |
|---|---|
| `x_time`, `y_time` | `seed_time` of the X / Y cluster (ns) |
| `time_diff` | `x_time − y_time` (signed, ns) |
| `adc_asym` | `(x_peak − y_peak) / (x_peak + y_peak)`, signed; `NaN` if the sum ≤ 0. `\|adc_asym\|` is exactly what `match_adc_asymmetry` cuts |
| `x/y_max_strip_dt`, `x/y_min_ts_corr` | copied from the X / Y cluster |

In match mode 1, `GEMHit`s exist only for pairs that **passed** the
X/Y cuts, so per-hit `adc_asym` / `time_diff` distributions are
post-cut. For the full distributions, pair X/Y clusters by peak rank
from the per-cluster values (as SBS does), or use a QA config with
`match_adc_asymmetry < 0` and `match_time_diff < 0`. With
`replay_recon -gem_hit` the per-hit values and a per-cluster block are
written to the recon tree (see
[`docs/REPLAYED_DATA.md`](../../REPLAYED_DATA.md)).

### Cuts and config keys

All keys live in `reconstruction_config.json` under `gem.default`
(and can be overridden per detector in `gem."0".."3"`, like every
other `ClusterConfig` knob):

| JSON key | `ClusterConfig` field | default (off) | SBS value | stage | rejects when |
|---|---|---|---|---|---|
| `strip_mean_time_range` | `strip_time_min`, `strip_time_max` | `[]` (±∞) | `[-99999, 99999]` (open); `[25, 150]` historically | strip | `t` outside `[min, max]` (inclusive), or `t` = `NaN`. With positive samples only, a 6-sample `t` always lies in [25, 150] ns, so a useful window must be narrower |
| `strip_unimodal_shape` | `strip_unimodal` | `false` | `true`, but compiled out | strip | `!IsUnimodalPulse(ts_adc)` |
| `seed_min_peak_adc` | same | `0` (≤ 0 = off) | 30 | cluster | `seed_peak_adc < value` |
| `seed_min_sum_adc` | same | `0` (≤ 0 = off) | 60 | cluster | `seed_sum_adc < value` |
| `strip_time_agreement` | same | `-1` (< 0 = off) | 50 ns | cluster | `max_strip_dt > value` |
| `strip_ts_corr_min` | same | `-1` (≤ −1 = off) | 0.7, but never called | cluster | `min_ts_corr < value` |

- **Strip-level** cuts run while grouping (`groupHits()`), with SBS
  `IsGoodStrip` semantics: a failing strip is left out and **ends the
  current run** of consecutive strips; the next run starts after it.
  The strip is not erased from the plane hits, so `GetPlaneHits()`
  and the raw tree still show it. A `[]` or `null`
  `strip_mean_time_range` disables the window (useful in a
  per-detector override); any other malformed value is ignored. The
  four scalar keys must be numbers, like every other numeric
  `ClusterConfig` key (a `null` or string aborts the config load); use
  the off values from the table to disable them.
- **Cluster-level** cuts run in `filterClusters()` after the size and
  cross-talk checks. Cross-talk flagging runs first, so a strong
  cluster that later fails a quality cut can still flag its
  cross-talk partners. `NaN` quality values (single-strip cluster,
  empty `ts_adc`) **pass** the cluster-level cuts.
- Because these `ClusterConfig` cuts live in `GemCluster`, they apply
  on every path that clusters (unlike the `GemSystem` strip keys, see
  *Raw → recon replay* below): `GemSystem::Reconstruct()` (server, replay from EVIO,
  Python), the replay raw → recon path and
  `hycal_shower_profile`.
- `PipelineBuilder` prints the values on the `[GEMCFG]` line
  (`strip_t=[min,max] unimodal= seed_peak= seed_sum= strip_dt= ts_corr=`).

### Mapping to SBS

| SBS cut (`gem_tracking*.conf`; `gem.conf` currently loads `gem_tracking_moller.conf`, same cut values) | active in SBS? | PRad-II equivalent |
|---|---|---|
| `max time bin = 1, 2, 3, 4` (strip) | **yes** | `reject_first_timebin` + `reject_last_timebin` (`GemSystem` strip cuts; production false/false, X17 true/false) |
| `strip mean time range` (strip) | wired, but the range is open | `strip_mean_time_range` |
| `use concave shape cut for strip` | **no** (`USE_STRIP_SHAPE_CUT` commented out) | `strip_unimodal_shape` |
| `seed strip min peak ADC` | **yes** (30) | `seed_min_peak_adc` |
| `seed strip min sum ADC` | **yes** (60) | `seed_min_sum_adc` |
| `strip mean time agreement` (seed vs strips) | **yes** (50 ns), buggy — see below | `strip_time_agreement` |
| `time sample correlation coefficient` | **no** (never called) | `strip_ts_corr_min` |
| `min/max cluster size` | yes (1 / 20) | `min_cluster_hits` / `max_cluster_hits` |
| cross-talk cluster removal | **no** (`setCrossTalk` commented out) | `charac_dists` etc. (on in PRad production) |
| `xy cluster matching mode` | mode 0 (rank pairing, no cuts) | `match_mode` (PRad production: 1) |
| `2d cluster adc assymetry` (0.8) | **no** (mode 1 only, and a no-op there) | `match_adc_asymmetry`, value on `GEMHit::adc_asym` |
| X/Y seed time agreement (50 ns) | **no** (mode 1 only) | `match_time_diff`, value on `GEMHit::time_diff` |

**SBS bugs deliberately not ported:**

- `Cuts::cluster_strip_time_agreement` loops
  `for (i = 0; i < size && i != seed; i++)`, which stops at the seed:
  only strips *before* the seed are checked (none if the seed is the
  first strip). PRad checks every `i ≠ seed`, so it is stricter for
  clusters with a late strip after the seed.
- Unqualified `abs(float)` in `Cuts.cpp` resolves to `int abs(int)`:
  the time difference is truncated (effective cut `|Δt| < 51` ns) and
  the X/Y ADC asymmetry becomes an integer division that is 0 for any
  two positive peaks. PRad uses float throughout.
- `Cuts::__get_sum_adc` accumulates into an `int` (the running sum is
  truncated after each sample). PRad sums in float.
- The SBS correlation returns 0 (with a warning on stdout) for
  mismatched lengths; PRad returns `NaN`.

**Mean-time definition.** SBS `Cuts::__get_mean_time` uses *all*
samples, negatives included, divides by the int-truncated sum and
returns 0 when that sum is 0; the SBS replay histogram helper
(`generate_gem_histos.h`) uses a float sum and returns 0 when it is
≤ 0. PRad keeps the positive-samples-only definition (it is the one
the existing X/Y time cut has always used) and returns `NaN` instead
of 0. For a clean pulse without negative samples the definitions
agree up to the SBS int truncation; baseline undershoot pulls the SBS
value down.

**Strip-level vs seed-level minimum ADC.** PRad's `min_peak_adc` /
`min_sum_adc` (production 30 / 60) are `GemSystem` strip cuts applied
to **every** strip in `GemSystem::collectHits()`, before clustering:
a weak strip is dropped entirely, which trims cluster edges and, with
`consecutive_thres = 1`, splits clusters. SBS applies the same numbers
to the **seed strip only**. To emulate SBS, set `min_peak_adc` and
`min_sum_adc` to 0 and use `seed_min_peak_adc = 30`,
`seed_min_sum_adc = 60`. (The `GemSystem` strip cuts, like
`reject_*_timebin`, are read only from `gem.default`; the new keys
are per-detector `ClusterConfig` knobs.) An SBS-like QA configuration
in `gem.default` is therefore:

```json
"reject_first_timebin": true,  "reject_last_timebin": true,
"min_peak_adc": 0.0,           "min_sum_adc": 0.0,
"seed_min_peak_adc": 30.0,     "seed_min_sum_adc": 60.0,
"strip_time_agreement": 50.0,
"charac_dists": [],            "match_mode": 0
```

It reproduces the SBS cut *logic*, not SBS numbers exactly: the mean
time, the full seed-vs-strip loop and the float arithmetic differ as
described above.

**Raw → recon replay.** The `GemSystem` strip keys in this block
(`reject_*_timebin`, `min_peak_adc`, `min_sum_adc`, like the
zero-suppression threshold) act only where strips are built from
EVIO, in `GemSystem::collectHits()`: the server, Python,
`replay_rawdata`, and `replay_recon` on EVIO input. `replay_recon` on a
`_raw.root` file (and `hycal_shower_profile`) rebuilds the strips from
the raw-tree `gem.*` arrays, which already carry the strip cuts of the
config `replay_rawdata` ran with (production: 30 / 60, no time-bin
rejection). On that path only the `ClusterConfig` keys (`seed_min_*`,
`strip_*`, `match_*`, `charac_dists`, …) take effect, even though the
`[GEMSYS]` line still prints the `-r` config's strip keys. Run the full
recipe on EVIO input (`replay_recon -r <qa_config>.json <evio>`);
`replay_rawdata` has no `-r` option and reads
`reconstruction_config.json` from the database directory.

### What the cuts do on real data

Measured on run 24246 (split 0, first 20k events, September 2026;
noise-burst events with ≥ 400 clusters excluded from the
distributions; efficiencies from `gem_eff_audit.py`, leave-one-out):

- **Negligible or no effect at the SBS values:** `strip_time_agreement = 50` (the
  seed-vs-strip spread ends at ≈ 52 ns; 13 of 80k multi-strip clusters
  exceed 50), `seed_min_peak_adc = 30` / `seed_min_sum_adc = 60` (the
  5σ zero suppression already guarantees them, even with
  `min_peak_adc = min_sum_adc = 0`) and `strip_mean_time_range =
  [25, 150]` (always satisfied, see the table). Seed times sit at
  ≈ 66–115 ns (p5–p95).
- **Edge time bins** (`reject_first/last_timebin`, SBS "max time bin
  = 1..4"): removes 4–18 % of clusters, including an out-of-time
  population (peak in the last sample, seed time > 120 ns, coincident
  in X and Y). It moves the efficiency by at most ≈ 3 points.
- **`strip_unimodal_shape`**: costs 3–8 efficiency points (GEM0 worst).
  SBS switched this cut off in May 2026.
- **`strip_ts_corr_min = 0.7`: do not use it as a cluster cut.** Weak
  edge strips are noisy, so `min_ts_corr` falls with cluster size and
  peak (≥ 79 % of 5+ strip clusters fail). The cut removes ≈ 40 % of
  hits, keeps single-strip noise clusters, and drops the efficiency
  from 62–87 % to 12–15 %. SBS never calls this cut.
- X/Y agreement is good on every detector: rank-paired peak
  correlation r ≈ 0.92, leading-pair `tx − ty` median within ±1.2 ns
  (p16–p84 half-width 2.5–5 ns), Y peaks ≈ 10–25 % larger than X
  (median asymmetry −0.05 to −0.11). The production `match_time_diff = 50` and
  `match_adc_asymmetry = 0.7` remove ≤ 0.1 % and ≤ 0.8 % of rank pairs.

## Reproducing the plots

The detector geometry is read from
[`database/gem_map.json`](../../../database/gem_map.json); the
strip-clustering and matching algorithms are re-implemented in pure
Python in [`plot_gem_clustering.py`](plot_gem_clustering.py) (NumPy +
Matplotlib only).

```bash
cd docs/technical_notes/gem_clustering
python plot_gem_clustering.py
```

Regenerates `plots/gem_fig1_layout.png`, `plots/gem_fig2_strip_clustering.png`,
`plots/gem_fig3_xy_matching.png`, `plots/gem_fig4_params.png` and prints the
reconstructed cluster table to stdout.

## See also

- [`prad2det/include/GemCluster.h`](../../../prad2det/include/GemCluster.h),
  [`GemCluster.cpp`](../../../prad2det/src/GemCluster.cpp) — algorithm source
- [`prad2det/include/GemSystem.h`](../../../prad2det/include/GemSystem.h) —
  hierarchy, pedestal/CM/zero-suppression, strip mapping, per-detector
  `ClusterConfig` storage
- [`database/gem_map.json`](../../../database/gem_map.json) —
  APV mapping + plane / pitch / hole geometry
- [`database/reconstruction_config.json`](../../../database/reconstruction_config.json) —
  per-run cluster-config defaults
- [`docs/REPLAYED_DATA.md`](../../REPLAYED_DATA.md) —
  branch layout for the recon tree (where `GEMHit`s land as `gem_*`)
- mpd_gem_view_ssp `GEMCluster` — original implementation lineage
