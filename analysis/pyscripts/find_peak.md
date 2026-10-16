# Visual walkthrough of `findPeaks`

I'll use ASCII diagrams because they render cleanly in this chat and preserve alignment. Each step shows a small toy waveform (~40 samples, one clean pulse + one small piled pulse) and marks what the algorithm is doing at that step.

The toy waveform for all diagrams:

```
sample index:  0    5   10   15   20   25   30   35   40
                                 ┌─┐
                                 │ │  ← peak A (large)
                                 │ │
                             ┌───┘ └──┐
                             │         └─┐        ┌─┐
                             │           └──┐     │ │ ← peak B (small, on tail of A)
                             │              └─┐   │ └─┐
                             │                └──┬┘   └───┐
   ped ~~~~~~~~~~~~~~~~~~~~~~┘                   ...      └~~~~~~~~
   noise:      ·  ··  · ·                                         · ·

                                 ↑ peak A raw max at sample ~17
                                                     ↑ peak B raw max at sample ~28
```

`buf[]` is the smoothed version of this; `raw[]` is the noisy original.

---

## Step 1 — Setup: two waveforms, one threshold

```
raw[]:  jagged 12-bit ADC samples (has noise wiggles)
buf[]:  smoothed by triangular kernel (noise reduced)

thr  =  max(peak_nsigma·rms, min_peak_height)
     =  max(5·0.5, 10)  =  10 ADC   (horizontal cut line)

trend_tol = max(0.1, 0.5·rms) = 0.25 ADC   (noise band for slope tests)
```

```
    ADC
     │           ┌─┐
     │           │ │  ← buf (smooth)  ─── heavy line
     │           │ │
     │       ┌──┘ └──┐
     │       │        └─┐        ┌─┐
     │       │          └──┐     │ │
─────┼───────┴──────────────┴────┴─┴─── thr = 10 ADC
     │
─────┴──────────────────────────────── ped_mean
       sample →
```

The scan runs on `buf`; heights and timing get pulled from `raw`.

---

## Step 2 — Local-max scan with `trend()`

For each interior sample `i`, compute two trends:

```
       buf[i-1]   buf[i]   buf[i+1]
          ↑         ↑         ↑
          └── tr1 ──┘         │
                    └── tr2 ──┘

     tr1 = trend(buf[i], buf[i-1])   +1 if buf[i] > left neighbor
     tr2 = trend(buf[i], buf[i+1])   +1 if buf[i] > right neighbor
```

The four cases:

```
   RISING SLOPE          FALLING SLOPE         LOCAL MAX            LOCAL MIN
       ●                       ●                   ●              
      /                         \                 / \              
     ●                           ●               ●   ●              ●   ●
    /                             \                                  \ /
   ●                               ●                                  ●

   tr1=+1, tr2=−1        tr1=−1, tr2=+1     tr1=+1, tr2=+1     tr1=−1, tr2=−1
   product = −1          product = −1        product = +1       product = +1
   REJECT                REJECT              CANDIDATE          candidate (fails
                                                                 height cut later)
```

Rule: reject if `tr1 * tr2 < 0` (monotone slope) or both flat.

---

## Step 3 — Plateau handling

Digitization can flatten the peak:

```
              ┌──┐     ← plateau: buf[17] == buf[18]
             /    \
            /      \
           /        \
          ●          ●
                                                     
     i=17: tr1=+1, tr2=0 (flat) — walk forward to find flat_end
     
     flat_end walks:  17 → 18  (still flat)
                      18 → 19  (falls, stop)
     
     flat_end = 18
     peak_pos = (17 + 18) / 2 = 17   ← center of plateau
```

If the plateau doesn't eventually fall (runs into buffer edge or rises again), reject — it's a step, not a peak.

---

## Step 4 — Peak-finding range: shoulders

Walk outward from the peak while the slope continues in the expected direction:

```
              ┌─┐                            
     left ← rising      falling → right      
             │  │  \                         
        ●   /│  │\  \                        
         \ / │  │ \  ●                       
          ●  │  │  \ │\                      
             │  │   \│ \                     
             │  │    ●  ●                    
     └───────┼──┼─────────┘                  
        left │  │        right               
             i  flat_end
```

- Left walk: `left--` while `trend(buf[left], buf[left-1]) > 0` (still rising into peak).
- Right walk: `right++` while `trend(buf[right], buf[right+1]) >= 0` (still falling or flat).

The asymmetry (strict `>` on left, `>=` on right) is deliberate — pulses have sharp rises and slow tails.

Result: **peak-finding range `[left, right]`** — cached for later overlap tests, but *not* used for integration.

---

## Step 5 — Local baseline and `smooth_height`

Interpolate a line between the two shoulders and measure the peak's rise above it:

```
                    ┌─┐  ← buf[peak_pos]
                    │ │
                    │ │  smooth_height = buf[peak_pos] − base
                    │ │      ↑
                    │ │      │
              base  │ │      ↓
     buf[left] ●····│·│············● buf[right]
                ····│·│············
                    │ │
                    └─┴─
                left   right
                    ↑
                 peak_pos

     base = (buf[left] · (right − peak_pos)   +   buf[right] · (peak_pos − left))
            ─────────────────────────────────────────────────────────────────
                              (right − left)
```

This is a linear interpolation of the baseline at `peak_pos`. On a clean isolated pulse, `base ≈ ped_mean`. On a peak riding a tail:

```
                    ┌─┐  peak B
                    │ │
                    │·│·······────  base (follows the tail of peak A)
              base  │·│         ● buf[right]
     ●·············│·│·······
        \···········│·│
         \··········│·│
          buf[left] └─┘
                                    smooth_height measures only 
                                    the incremental bump above tail
```

Reject if `smooth_height < thr`.

---

## Step 6 — Pedestal-relative height (second cut)

Also require `ped_height = buf[peak_pos] − ped_mean ≥ thr`. On a normal pulse this is redundant with step 5; on a drifted baseline it isn't.

```
                    ┌─┐
                    │ │  ped_height (rise above pedestal)
                    │ │       ↑
                    │ │       │
                    │ │       │
    ────────────────┴─┴──── ped_mean
                              ↓
```

---

## Step 7 — Tail integration

Walk outward from `peak_pos` summing pedsub samples. Terminate on `N_break = 2` **consecutive** below-threshold samples:

```
                    peak_pos
                       ↓
                      ┌─┐
                      │ │
                     ┌┘ └┐
                    ┌┘   └┐
                   ┌┘     └┐            ← above tail_cut, commit
                  ┌┘       └─┐          ← above, commit
                 ┌┘          └─┐        ← above, commit
                ┌┘             └─┐      ← above, commit
   ─ ── ── ── ── ── ── ── ── ── ── ── ─  tail_cut = int_tail_ratio · ped_height
                                └─┐     ← BELOW (below_run=1), hold in `pending`
                                  └─●   ← ABOVE! commit pending + this sample,
                                          reset below_run=0, pending=0
                                   └─┐  ← BELOW (below_run=1), hold in `pending`
                                     └─ ← BELOW (below_run=2) → STOP
                                     
     int_left ↑                       ↑ int_right
     (last committed sample on each side)
```

Key mechanics:

- Below-threshold samples go into `pending` (not lost, not committed yet).
- A single noise dip doesn't truncate — the walk continues.
- `N_break = 2` consecutive dips finally terminate the walk.
- On termination, `pending` samples in the trailing dip are **discarded** (they were below threshold, and there was no recovery).
- `int_left` and `int_right` only advance when a sample is *actually* committed — so they always point to a real above-threshold sample.

Result: **integration range `[int_left, int_right]`** — inclusive bounds of samples that actually contributed to `integral`. Distinct from `[left, right]`.

---

## Step 8 — Raw-peak refinement

Smoothing shifts the maximum by a sample or two. Search a small window in `raw[]` for the true max:

```
    smoothed peak (buf)     
         ┌─┐                
         │ │                
        ┌┘ └┐               
        │    │              
        ↑                    
     peak_pos (smooth)      

    raw samples in the same region:  (noisy!)
         ● ● ●               
        ● ↓ ↑ ●              
       ●   ● raw max is one sample right
       ↓             
     peak_pos                 raw_pos = peak_pos + 1
                              raw_height = raw[raw_pos] − ped_mean
```

Search window is `smooth_order + plateau_half_width` samples on each side. `raw_pos` and `raw_height` are what get stored in `Peak::pos` and `Peak::height`.

---

## Step 9 — Overlap rejection (`min_peak_ratio`)

Before accepting the new peak, check if its peak-finding range overlaps any previous peak's range. If yes, and if the new peak is too small compared to the neighbor, reject it as a shoulder:

```
    Case A: no overlap → accept
    
       peak 0 range         peak 1 range
       ┌────────┐           ┌────────┐
       │        │           │        │
       ↓        ↓           ↓        ↓
       ●●●●●●●●●             ●●●●●●●●●    → both accepted


    Case B: overlap, new peak is tall enough → accept
    
       peak 0 range     peak 1 range
       ┌────────┐───────┐
       │      ┌─┴──┐    │
       ↓      │    ↓    ↓
       ●●●●●●● ●●●●●●●●●●●●●●
       (both smooth_heights are large)  → both accepted


    Case C: overlap, new peak is a small shoulder → reject
    
       peak 0 range   peak 1 range
       ┌────────┐──┐
       │       ┌┴┐ │
       │       │·│ │      smooth_height(new) < 0.3 · height(neighbor)
       ↓       │·│ ↓
       ●●●●●●●●●●●●●●     → peak 1 REJECTED as a shoulder of peak 0
```

The condition is `smooth_height < result.peaks[k].height · cfg.min_peak_ratio` (default 0.3).

---

## Step 10 — Sub-sample time (quadratic vertex)

Fit a parabola through three raw samples around the peak:

```
                    h_zero
                     ●
                    /│\
                   / │ \
                  /  │  \
          h_minus●   │   ●h_plus
                    /│\
                   / │ \    ← parabola vertex sits at δ from center
                  /  │  \
                     δ
                     
     δ = 0.5 · (h_minus − h_plus) / (h_minus − 2·h_zero + h_plus)
     
     δ clamped to [−1, +1]
     
     denom < 0 required (concave down = real peak)
```

Lifts timing resolution from the 4 ns sample grid to well under 1 ns.

---

## Step 11 — Linear CFD (50% constant fraction)

For pulses with `raw_height > 5·rms`, scan the rising edge for the 50% crossing:

```
        ┌───    v_thr = 0.5 · raw_height
         raw_height  ●●●●●● ← raw_pos
           ↑        /
           │       /
    v_thr─┼─── ●···  ← v1 (first sample ≥ v_thr)
           │   /│
           │  / │  ← linear interpolate the crossing
           │ / ·│
           │/·  │
     ● ● ●·    ●
         ↑     ↑
         v0    j+1
         j
    
    frac = (v_thr − v0) / (v1 − v0)
    t_pickoff = j + frac
```

Extra guard: `dv = v1 − v0 > 1.5·rms` — the crossing must be a real rise, not a noise wiggle.

Timing hierarchy so far:
- always: quadratic vertex (`~1 ns`)
- `raw_height > 5·rms`: linear CFD (`~100 ps`)
- `raw_height > 10·rms`: log-normal CFD fit (`~10 ps`)

---

## Step 12 — Pile-up flag

Compare **integration ranges** (not peak-finding ranges) of every accepted peak against this new one. Overlap or near-overlap → flag both:

```
    Peak 0 integration range         Peak 1 integration range
    ┌────────────────────┐           ┌──────────────┐
    │       ●●●●         │           │    ●●●       │
    │      ●    ●        │           │   ●   ●      │
    ↓     ●      ●●●●●   ↓           ↓  ●     ●●    ↓
    int_left           int_right     int_left    int_right (peak 1)
                          │              │
                          └── gap ───────┘
                       
    gap = int_left(1) − int_right(0)
    
    If  gap ≤ peak_pileup_gap (default 2 samples):
        peak 0: quality |= Q_PEAK_PILED
        peak 1: quality |= Q_PEAK_PILED
    
    (Flag is symmetric — both peaks get marked)
```

Test in code: `int_left ≤ prev.right + gap  AND  int_right ≥ prev.left − gap`.

---

## Step 13 — Store the peak, skip past it

```
    Peak fields:
    ┌─────────────────────────────────────────┐
    │  pos       = raw_pos                     │
    │  left      = int_left                    │  ← integration bounds
    │  right     = int_right                   │
    │  height    = raw_height                  │
    │  integral  = Σ pedsub samples            │
    │  time      = t_pickoff · 1000 / clk_mhz  │  ← ns
    │  overflow  = (raw[raw_pos] ≥ 4095)       │
    │  quality   = Q_PEAK_GOOD | Q_PEAK_PILED  │
    │  time_algo = which method produced time  │
    └─────────────────────────────────────────┘

    pk_range[k] = [left, right]   ← peak-finding bounds, cached for next iter

    result.npeaks++
    
    i = right   ← skip past this peak's range so the falling edge isn't rescanned
```

---

## Big picture: two ranges per peak

The single most confusing thing about this function is that every peak carries **two** ranges. Here's the difference visually:

```
                        peak_pos
                          ↓
                         ┌─┐    
                        ┌┘ └┐   
                       ┌┘   └─┐ 
                      ┌┘      └─┐
                     ┌┘         └─┐
                    ┌┘            └─┐
                   ┌┘               └───┐
                  ┌┘                    └────────
                                                
   [────────────────────────────────────]        ← [left, right]      peak-finding range
                                                    (slope-driven; used for overlap test)
                                                
        [────────────────────────────]            ← [int_left, int_right]  integration range
                                                    (threshold-driven; used for integral,
                                                     pile-up test, downstream shading)
```

- **Peak-finding range** tracks the *entire hump*, following the slope out to where it flattens.
- **Integration range** tracks *above-threshold samples*, tighter because the tail falls below `tail_cut` before it flattens.

Mixing them up is a subtle bug: overlap rejection needs the wider range (shoulders that touch), pile-up flagging needs the tighter range (integrated energy that overlaps). The source is careful to keep them separate — and now, hopefully, so is your mental model.
