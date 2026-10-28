#!/usr/bin/env python3
"""
gem_xy_corr.py — X/Y cluster time-sample correlation of GEM 2-D hits.

rho_cl is the Pearson correlation, over the APV25 time samples s = 0..5, of
the cluster-summed samples U_s / V_s of the X and Y cluster of a 2-D hit
(the SBS-offline corrcoeff_clust):

    rho = sum_s (U_s - <U>)(V_s - <V>)
          / sqrt( sum_s (U_s - <U>)^2 * sum_s (V_s - <V>)^2 )

`prad2ana_replay_recon -gem_hit` writes it per hit (gem_xy_corr) and writes
the cluster-summed samples of every 1-D cluster (gem_cl_ts_adc[n_gem_cl][6],
StripCluster::ts_adc_sum), so rho can be formed for ANY X/Y cluster pair,
including pairs the match_mode 1 cuts rejected.  ts_corr() repeats the
arithmetic of gem::TimeSampleCorrelation (float32 in, sequential float64
sums, NaN for a flat waveform), so it reproduces gem_xy_corr exactly.

Library use (from analysis/pyscripts):
    from gem_xy_corr import ts_corr, xy_pairs
    rho = ts_corr(u, v)                    # (..., n_samples) arrays -> (...)
    ix, iy, rho = xy_pairs(ev_cl_det, ev_cl_plane, ev_cl_ts_adc, det)

CLI — per GEM, quantiles (to about +-0.001) and fractions above a few cuts of
gem_xy_corr (stored hits) and of rho over all X x Y cluster pairs of the
per-cluster block.  Events whose hit or cluster block is truncated at 400
(noise bursts) are skipped, so both rows cover the same events.  The tree is
read in chunks, so merged run files are fine.
    python analysis/pyscripts/gem_xy_corr.py prad_024246_recon.root [-n 20000]

Recon files replayed before 2026-10 have neither branch: rerun
prad2ana_replay_recon -gem_hit on the EVIO files, or on _raw.root files
written since 2026-08 (per-strip gem.ts_adc; older raw files use the legacy
GEM schema and must be replayed from EVIO).
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

MAX_GEM_HITS = 400            # prad2::kMaxGemHits
MAX_GEM_CLUSTERS = 400        # prad2::kMaxGemClusters
QUANTILES = (5, 16, 50, 84, 95)
CUTS = (0.0, 0.5, 0.8, 0.9)
EDGES = np.linspace(-1., 1., 4001)   # 0.0005-wide bins for the CLI quantiles


def ts_corr(u, v) -> np.ndarray:
    """Pearson r over the last axis of two broadcastable sample arrays.

    Same arithmetic as gem::TimeSampleCorrelation: samples taken as float32,
    accumulated in float64 in sample order, clamped to [-1, 1], returned as
    float32.  NaN where either waveform is flat (or has NaN samples), or
    when there are fewer than 2 samples or the sample counts differ.
    """
    u = np.asarray(u, dtype=np.float32).astype(np.float64)
    v = np.asarray(v, dtype=np.float32).astype(np.float64)
    shape = np.broadcast_shapes(u.shape[:-1], v.shape[:-1])
    n = u.shape[-1]
    if n < 2 or v.shape[-1] != n:
        return np.full(shape, np.nan, dtype=np.float32)

    mean_u = np.zeros(u.shape[:-1])
    mean_v = np.zeros(v.shape[:-1])
    for s in range(n):
        mean_u = mean_u + u[..., s]
        mean_v = mean_v + v[..., s]
    mean_u = mean_u / n
    mean_v = mean_v / n

    s_uv = np.zeros(shape)
    s_uu = np.zeros(u.shape[:-1])
    s_vv = np.zeros(v.shape[:-1])
    for s in range(n):
        du = u[..., s] - mean_u
        dv = v[..., s] - mean_v
        s_uv = s_uv + du * dv
        s_uu = s_uu + du * du
        s_vv = s_vv + dv * dv
    s_uu, s_vv = np.broadcast_to(s_uu, shape), np.broadcast_to(s_vv, shape)

    ok = (s_uu > 0.) & (s_vv > 0.)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.clip(s_uv / np.sqrt(s_uu * s_vv), -1., 1.)
    return np.where(ok, r, np.nan).astype(np.float32)


def xy_pairs(cl_det, cl_plane, cl_ts_adc, det):
    """rho for every X x Y cluster pair of one GEM in one event.

    cl_det, cl_plane, cl_ts_adc are one event's gem_cl_det, gem_cl_plane and
    gem_cl_ts_adc (shape (n_gem_cl, 6)).  Returns (ix, iy, rho): the
    gem_cl_* indices of the X and Y clusters and rho[len(ix), len(iy)].
    """
    cl_det, cl_plane = np.asarray(cl_det), np.asarray(cl_plane)
    ts = np.asarray(cl_ts_adc, dtype=np.float32)
    ix = np.flatnonzero((cl_det == det) & (cl_plane == 0))
    iy = np.flatnonzero((cl_det == det) & (cl_plane == 1))
    return ix, iy, ts_corr(ts[ix][:, None, :], ts[iy][None, :, :])


class _Summary:
    """Histogram of rho on [-1, 1] plus exact counts above the CUTS."""

    def __init__(self):
        self.counts = np.zeros(len(EDGES) - 1, dtype=np.int64)
        self.n_ge = np.zeros(len(CUTS), dtype=np.int64)
        self.n_nan = 0

    def fill(self, rho):
        rho = np.asarray(rho, dtype=np.float32).ravel()
        ok = np.isfinite(rho)
        self.n_nan += int(np.count_nonzero(~ok))
        rho = rho[ok]
        self.counts += np.histogram(rho, bins=EDGES)[0]
        self.n_ge += [np.count_nonzero(rho >= c) for c in CUTS]

    def line(self, name: str) -> str:
        n = int(self.counts.sum())
        if n == 0:
            return f"  {name:<14} n=0 (NaN {self.n_nan})"
        cum = np.cumsum(self.counts)
        qs = []
        for p in QUANTILES:
            t = p / 100. * n
            i = min(int(np.searchsorted(cum, t)), len(self.counts) - 1)
            frac = (t - (cum[i] - self.counts[i])) / self.counts[i] if self.counts[i] else 0.
            qs.append(f"p{p}={EDGES[i] + frac * (EDGES[i + 1] - EDGES[i]):+.3f}")
        fs = " ".join(f"f>={c:g}:{k / n:.3f}" for c, k in zip(CUTS, self.n_ge))
        return f"  {name:<14} n={n:<9d} {' '.join(qs)}  {fs}"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="recon ROOT file (replay_recon -gem_hit)")
    ap.add_argument("-n", "--max-events", type=int, default=None,
                    help="read at most this many events")
    ap.add_argument("--tree", default="recon")
    args = ap.parse_args()
    if args.max_events is not None and args.max_events < 0:
        ap.error("-n must be >= 0")

    try:
        import uproot
    except ImportError as e:
        raise SystemExit(f"[gem_xy_corr] uproot required (pip install uproot): {e}")

    branches = ["n_gem_hits", "det_id", "gem_xy_corr",
                "n_gem_cl", "gem_cl_det", "gem_cl_plane", "gem_cl_ts_adc"]
    try:
        f = uproot.open(args.input)
    except (OSError, ValueError) as e:
        print(f"[gem_xy_corr] cannot open {args.input}: {e}", file=sys.stderr)
        return 1
    with f:
        if args.tree not in f:
            print(f"[gem_xy_corr] {args.input}: no '{args.tree}' tree", file=sys.stderr)
            return 1
        tree = f[args.tree]
        missing = [b for b in branches if b not in tree]
        if missing:
            print(f"[gem_xy_corr] {args.input}: missing {', '.join(missing)} — "
                  "replay with prad2ana_replay_recon -gem_hit (2026-10 or later)",
                  file=sys.stderr)
            return 1

        hits = [_Summary() for _ in range(4)]
        pairs = [_Summary() for _ in range(4)]
        n_used = n_trunc = 0
        for data in tree.iterate(branches, step_size=50_000,
                                 entry_stop=args.max_events, library="np"):
            hit_rho = [[] for _ in range(4)]
            pair_rho = [[] for _ in range(4)]
            for k in range(len(data["n_gem_cl"])):
                if (data["n_gem_hits"][k] >= MAX_GEM_HITS or
                        data["n_gem_cl"][k] >= MAX_GEM_CLUSTERS):
                    n_trunc += 1
                    continue
                n_used += 1
                det_id, xy_corr = data["det_id"][k], data["gem_xy_corr"][k]
                cl_det, cl_plane = data["gem_cl_det"][k], data["gem_cl_plane"][k]
                cl_ts = np.asarray(data["gem_cl_ts_adc"][k]).reshape(-1, 6)
                for d in range(4):
                    hit_rho[d].append(xy_corr[det_id == d])
                    pair_rho[d].append(xy_pairs(cl_det, cl_plane, cl_ts, d)[2].ravel())
            for d in range(4):
                if hit_rho[d]:
                    hits[d].fill(np.concatenate(hit_rho[d]))
                    pairs[d].fill(np.concatenate(pair_rho[d]))

    print(f"{args.input}: {n_used} events used, {n_trunc} with a truncated hit or "
          f"cluster block (n_gem_hits or n_gem_cl == 400) skipped")
    for d in range(4):
        print(f"GEM{d}")
        print(hits[d].line("hits"))
        print(pairs[d].line("all XxY pairs"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
