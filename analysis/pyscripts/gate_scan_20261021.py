#!/usr/bin/env python3
"""gate_scan_20261021.py -- timing-gate scan harness (plan_20261021.md, Phase 1).

Scans the cluster timing gate of benchmark_hycal_timing_test_20261021.py over
(time reference x W x pre-window) on a fixed set of runs/splits, in stages
A -> B1 -> B2 (tune, half A) and C (evaluate, half B).  The driver
benchmark_hycal_timing_scan.py is imported read-only and is not modified.

Subcommands
-----------
  run        launch the benchmark jobs of one stage (no aggregation)
  summarise  parse every result.tsv, write tables, plots, determinism check,
             and the stage-C markdown report
  select     apply the plan 3.5 selection rule on half A and write
             selection_<stage>.json

Usage (see plan 3.7)
  G="python3 gate_scan_20261021.py"; P="--out $OUT --parallel 12 --timeout-seconds 1800"
  $G run --stage A --half A $P --dry-run | tail -3
  nohup $G run --stage A --half A $P > $OUT/logs/A.out 2>&1
  $G summarise --out $OUT && $G select --stage A --out $OUT

Layout:  <OUT>/<cfg>/scan_RRRRRR/split_SSSSS/{result.tsv,benchmark.log,job.cmd}
Config tag:  t{raw|cal}_w{W}_p{LO-HI|map}   (pmap = driver's per-run pre-window)
"""
from __future__ import annotations

import argparse
import collections
import datetime
import hashlib
import json
import math
import os
import re
import shlex
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")      # non-GUI backend: avoids plotting hangs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_hycal_timing_scan import (          # read-only use of the driver
    _HC_CALIB_MAP, _PRE_WINDOW_MAP, _DEFAULT_PRE_WINDOW, _DATABASE_DIR,
    _BENCHMARK_SCRIPT, _DEFAULT_EVIO_DIR, gain_corr_path, _EXCLUDED_RUNS,
    parse_tsv,
)

NAN = float("nan")
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
REF_TAG = "traw_w8_pmap"        # today's setting; legacy row = reference
GATE_OFF_W = 100.0              # W >= this: gate-off reference, never a candidate
EPS_MIN = 0.99                  # eps_yield and eps_peak floor (D3)
STAGES = ("A", "B1", "B2", "C")
SEL_STAGES = ("A", "B1", "B2")

# ---------------------------------------------------------------------------
# Run types and splits (plan 3.3).  (half A = tune, half B = evaluate)
# ---------------------------------------------------------------------------
_S10_A, _S10_B = [10, 30, 50, 70], [20, 40, 60, 80]
_S50_A, _S50_B = [50, 150, 250, 350], [100, 200, 300, 400]
RUN_TYPES = [
    {"name": "3.49 GeV H2", "pref": "cb",
     "runs": {24503: (_S10_A, _S10_B)}},
    {"name": "2.24 GeV H2", "pref": "gauss",
     "runs": {25320: (_S50_A, _S50_B)}},
    {"name": "2.24 GeV C", "pref": "cb",
     "runs": {25308: (_S50_A, _S50_B)}},
    {"name": "0.73 GeV H2", "pref": "gauss",
     "runs": {24590: (_S10_A, _S10_B)}},
    {"name": "2.24 GeV Ta X17", "pref": "cb",
     "runs": {25800: ([100, 300, 500, 700], [200, 400, 600, 800]),
              26138: ([99, 299], [199, 399]),
              26230: (list(range(100, 1301, 200)), list(range(200, 1201, 200)))}},
]
RT_NAMES = [rt["name"] for rt in RUN_TYPES]
RT_BY_NAME = {rt["name"]: rt for rt in RUN_TYPES}
HALF_SPLITS: dict[str, dict[str, list]] = {}    # rtname -> half -> [(run, split)]
RT_OF_RUN: dict[int, str] = {}
HALF_OF: dict[tuple, str] = {}
for _rt in RUN_TYPES:
    HALF_SPLITS[_rt["name"]] = {"A": [], "B": []}
    for _run, (_sa, _sb) in _rt["runs"].items():
        if _run in _EXCLUDED_RUNS:
            continue
        RT_OF_RUN[_run] = _rt["name"]
        for _h, _lst in (("A", _sa), ("B", _sb)):
            for _s in _lst:
                HALF_SPLITS[_rt["name"]][_h].append((_run, _s))
                HALF_OF[(_run, _s)] = _h

# ---------------------------------------------------------------------------
# Configs
# ---------------------------------------------------------------------------
def _n(x) -> str:
    """Integral floats print as integers (tag convention)."""
    x = float(x)
    return str(int(x)) if x.is_integer() else repr(x)


def make_cfg(time_calib, W, LO=None, HI=None, bench_args="", tag=None) -> dict:
    if time_calib not in ("off", "auto"):
        raise ValueError(f"time_calib must be 'off' or 'auto', got {time_calib!r}")
    if (LO is None) != (HI is None):
        raise ValueError("LO and HI must both be given or both be None (pmap)")
    W = float(W)
    LO = None if LO is None else float(LO)
    HI = None if HI is None else float(HI)
    if tag is None:
        win = "map" if LO is None else f"{_n(LO)}-{_n(HI)}"
        tag = f"t{'raw' if time_calib == 'off' else 'cal'}_w{_n(W)}_p{win}"
    return {"tag": tag, "time_calib": time_calib, "W": W, "LO": LO, "HI": HI,
            "bench_args": bench_args or ""}


_TAG_RX = re.compile(r"^t(raw|cal)_w([0-9.]+)_p(?:map|([0-9.]+)-([0-9.]+))$")


def parse_tag(tag: str):
    """-> (time_calib, W, LO, HI) from a standard tag, or None."""
    m = _TAG_RX.match(tag)
    if not m:
        return None
    lo = None if m.group(3) is None else float(m.group(3))
    hi = None if m.group(4) is None else float(m.group(4))
    return ("off" if m.group(1) == "raw" else "auto", float(m.group(2)), lo, hi)


def cfg_from_tag(tag: str) -> dict:
    p = parse_tag(tag)
    if p is None:
        raise SystemExit(f"tag {tag!r} does not follow t{{raw|cal}}_w<W>_p<LO-HI|map>")
    return make_cfg(*p, tag=tag)


def load_configs_file(path: str) -> list[dict]:
    raw = json.loads(Path(path).read_text())
    if not isinstance(raw, list):
        raise SystemExit(f"--configs-file {path}: expected a JSON list of config dicts")
    cfgs, seen = [], set()
    for d in raw:
        c = make_cfg(d.get("time_calib", "off"), d.get("W", 8.0), d.get("LO"),
                     d.get("HI"), d.get("bench_args", ""), d.get("tag"))
        if c["tag"] in seen:
            raise SystemExit(f"--configs-file: duplicate tag {c['tag']}")
        seen.add(c["tag"])
        cfgs.append(c)
    return cfgs


def dedup(cfgs: list[dict]) -> list[dict]:
    out, seen = [], set()
    for c in cfgs:
        if c["tag"] not in seen:
            seen.add(c["tag"])
            out.append(c)
    return out


# ---------------------------------------------------------------------------
# Selections (selection_<stage>.json) and stage grids
# ---------------------------------------------------------------------------
def parse_overrides(items) -> dict:
    out = {}
    for it in items or []:
        if "=" not in it:
            raise SystemExit(f"--select-override expects STAGE=TAG, got {it!r}")
        st, tag = it.split("=", 1)
        if st not in SEL_STAGES:
            raise SystemExit(f"--select-override stage must be one of {SEL_STAGES}")
        out[st] = tag
    return out


def sel_path(d, stage) -> Path:
    return Path(d) / f"selection_{stage}.json"


def get_selection(stage, out, baseline, overrides) -> dict:
    """Winner config of an earlier stage: override > --baseline > --out."""
    if stage in overrides:
        return cfg_from_tag(overrides[stage])
    tried = []
    for d in ([baseline] if baseline else []) + [out]:
        p = sel_path(d, stage)
        tried.append(str(p))
        if p.exists():
            w = json.loads(p.read_text()).get("winner")
            if not w:
                raise SystemExit(f"{p} has no winner (selection failed); fix and re-run select")
            return make_cfg(w["time_calib"], w["W"], w.get("LO"), w.get("HI"),
                            "", w["tag"])
    raise SystemExit(
        f"missing stage-{stage} selection (looked for: {', '.join(tried)}).\n"
        f"Run `select --stage {stage} --out DIR` first, or pass --baseline DIR0, "
        f"or --select-override {stage}=TAG.")


def stage_cfgs(stage, sel) -> list[dict]:
    """The grid of one stage; sel(stage) -> winner cfg of an earlier stage."""
    if stage == "A":
        return ([make_cfg("off", w) for w in (4, 8, 16, 100)]
                + [make_cfg("auto", w) for w in (2, 3, 4, 6, 8, 12, 16)])
    if stage == "B1":
        a = sel("A")
        return [make_cfg(a["time_calib"], a["W"], lo, 200) for lo in (100, 120, 130, 140)]
    if stage == "B2":
        b = sel("B1")
        return [make_cfg(b["time_calib"], b["W"], b["LO"], hi) for hi in (180, 190, 220)]
    if stage == "C":
        a, b2 = sel("A"), sel("B2")
        return dedup([
            make_cfg("off", 8),                                       # C1
            make_cfg(a["time_calib"], a["W"]),                        # C2 (pmap)
            make_cfg(b2["time_calib"], b2["W"], b2["LO"], b2["HI"],
                     tag=b2["tag"]),                                  # C3
            make_cfg("auto", 4, 130, 200),                            # C4
        ])
    raise ValueError(stage)


def select_candidates(stage, sel) -> list[dict]:
    cfgs = stage_cfgs(stage, sel)
    if stage == "B2":                      # HI=200 comes from the B1 winner
        b = sel("B1")
        cfgs = cfgs + [make_cfg(b["time_calib"], b["W"], b["LO"], 200)]
    return cfgs


# ---------------------------------------------------------------------------
# Jobs
# ---------------------------------------------------------------------------
def build_cmd(cfg, run, evio, stem, lo, hi, bench_args) -> list[str]:
    """Driver command (run_benchmark) + --window/--time-calib + bench args."""
    cmd = [sys.executable, _BENCHMARK_SCRIPT, str(evio), str(stem),
           "--calibrate", "--exclude-rowcol", "15", "20", "15", "20",
           "--pre-window", str(float(lo)), str(float(hi)),
           "--inner-radius", "200", "--min-per-module", "30",
           "--lms-corr", "auto" if gain_corr_path(run).exists() else "off"]
    if run in _HC_CALIB_MAP:
        cmd += ["--hc-calib-file",
                str(_DATABASE_DIR / "calibration" / _HC_CALIB_MAP[run])]
    cmd += ["--window", str(float(cfg["W"])), "--time-calib", cfg["time_calib"]]
    cmd += shlex.split(bench_args)
    return cmd


def plan_jobs(cfgs, half, out: Path, evio_dir: Path, extra_args: str) -> list[dict]:
    jobs = []
    for cfg in cfgs:
        args = " ".join(x for x in (cfg["bench_args"], extra_args) if x).strip()
        for rtname in RT_NAMES:
            for run, split in HALF_SPLITS[rtname][half]:
                if cfg["LO"] is None:
                    lo, hi = _PRE_WINDOW_MAP.get(run, _DEFAULT_PRE_WINDOW)
                else:
                    lo, hi = cfg["LO"], cfg["HI"]
                sdir = out / cfg["tag"] / f"scan_{run:06d}" / f"split_{split:05d}"
                evio = evio_dir / f"prad_{run:06d}.evio.{split:05d}"
                jobs.append({
                    "cfg": cfg, "run": run, "split": split, "lo": lo, "hi": hi,
                    "bench_args": args, "evio": evio, "dir": sdir,
                    "tsv": sdir / "result.tsv", "log": sdir / "benchmark.log",
                    "cmd": build_cmd(cfg, run, evio, sdir / "result", lo, hi, args),
                })
    return jobs


# ---------------------------------------------------------------------------
# result.tsv parsing (driver parse_tsv for the tables + our own header scan)
# ---------------------------------------------------------------------------
def _f(s) -> float:
    try:
        return float(s)
    except (TypeError, ValueError):
        return NAN


_RX_PRE = re.compile(r"pre_window=\[([^,\]]+),([^\]]+)\]")
_RX_SEEDW = re.compile(r"seed_time_window=(\S+)")
_RX_NPHYS = re.compile(r"n_phys=(\S+)\s+n_kept=(\S+)\s+elapsed_s=(\S+)")
_RX_TARGET = re.compile(r"cal_method=\S+\s+target_mean_MeV=(\S+)")
_RX_TC = re.compile(r"time_calib=(\S+)(?:\s+n_offsets=(\S+))?")
_RX_CC = re.compile(r"cluster_counts\s+legacy\s+n0=(\S+)\s+n1=(\S+)\s+n2plus=(\S+)"
                    r"\s+new\s+n0=(\S+)\s+n1=(\S+)\s+n2plus=(\S+)")
_RX_PAIRED = re.compile(r"paired\s+n_same_seed=(\S+)\s+dE_frac_mean=(\S+)")
_RX_PEAK = re.compile(r"peak_yield\s.*?s_ref=(\S+)\s+modules_common=(\S+)"
                      r"\s+legacy=(\S+)\s+new=(\S+)")


def parse_header(path: Path) -> dict:
    """Scan all '#' lines (first match wins); missing -> NaN / None."""
    h = {"pre": (NAN, NAN), "W": NAN, "n_phys": NAN, "n_kept": NAN,
         "elapsed_s": NAN, "target_mean": NAN, "time_calib": None,
         "n_offsets": NAN,
         "cc": {"legacy": (NAN,) * 3, "new": (NAN,) * 3},
         "n_same_seed": NAN, "dE_frac_mean": NAN,
         "peak": {"s_ref": NAN, "modules": NAN, "legacy": NAN, "new": NAN}}
    seen = set()
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("#"):
            continue
        for key, rx in (("pre", _RX_PRE), ("W", _RX_SEEDW), ("np", _RX_NPHYS),
                        ("tg", _RX_TARGET), ("tc", _RX_TC), ("cc", _RX_CC),
                        ("pd", _RX_PAIRED), ("pk", _RX_PEAK)):
            if key in seen:
                continue
            m = rx.search(line)
            if not m:
                continue
            seen.add(key)
            g = m.groups()
            if key == "pre":
                h["pre"] = (_f(g[0]), _f(g[1]))
            elif key == "W":
                h["W"] = _f(g[0])
            elif key == "np":
                h["n_phys"], h["n_kept"], h["elapsed_s"] = _f(g[0]), _f(g[1]), _f(g[2])
            elif key == "tg":
                h["target_mean"] = _f(g[0])
            elif key == "tc":
                h["time_calib"], h["n_offsets"] = g[0], _f(g[1])
            elif key == "cc":
                v = [_f(x) for x in g]
                h["cc"] = {"legacy": tuple(v[:3]), "new": tuple(v[3:])}
            elif key == "pd":
                h["n_same_seed"], h["dE_frac_mean"] = _f(g[0]), _f(g[1])
            elif key == "pk":
                h["peak"] = {"s_ref": _f(g[0]), "modules": _f(g[1]),
                             "legacy": _f(g[2]), "new": _f(g[3])}
    return h


def _complete(t: dict) -> bool:
    return (all(p in t["headline"] for p in ("legacy", "new"))
            and all(p in t["calibration"] for p in ("legacy", "new")))


def load_result(tsv: Path, tag: str, run: int, split: int):
    try:
        t = parse_tsv(tsv)
        hdr = parse_header(tsv)
    except Exception as exc:                      # noqa: BLE001
        print(f"WARNING: cannot parse {tsv}: {exc}", file=sys.stderr)
        return None
    if not _complete(t):
        print(f"WARNING: {tsv}: headline/calibration block incomplete; skipped",
              file=sys.stderr)
        return None
    rtname = RT_OF_RUN.get(run, "?")
    pref = RT_BY_NAME[rtname]["pref"] if rtname in RT_BY_NAME else "gauss"
    return {"tag": tag, "run": run, "split": split, "rt": rtname, "pref": pref,
            "half": HALF_OF.get((run, split), "?"), "hdr": hdr,
            "head": t["headline"], "cal": t["calibration"], "path": tsv}


def load_tree(root: Path) -> dict:
    """tag -> {(run, split): record}"""
    tree: dict = {}
    root = Path(root)
    if not root.is_dir():
        return tree
    for cfgdir in sorted(p for p in root.iterdir() if p.is_dir()):
        for tsv in sorted(cfgdir.glob("scan_*/split_*/result.tsv")):
            m1 = re.match(r"scan_0*(\d+)$", tsv.parent.parent.name)
            m2 = re.match(r"split_(\d+)$", tsv.parent.name)
            if not (m1 and m2):
                continue
            rec = load_result(tsv, cfgdir.name, int(m1.group(1)), int(m2.group(1)))
            if rec is not None:
                tree.setdefault(cfgdir.name, {})[(rec["run"], rec["split"])] = rec
    return tree


_SIGCOLS = {"G_pre": ("sigma_pre_MeV", "mu_pre_MeV"),
            "G_post": ("sigma_post_MeV", "mu_post_MeV"),
            "CB_pre": ("sigma_cb_pre_MeV", "mu_cb_pre_MeV"),
            "CB_post": ("sigma_cb_MeV", "mu_cb_MeV")}


def sigma(rec, path, kind=None) -> float:
    """sigma/E (fraction) from the calibrated-block MeV columns."""
    kind = kind or ("CB_post" if rec["pref"] == "cb" else "G_post")
    c = rec["cal"].get(path)
    if not c:
        return NAN
    s, m = c.get(_SIGCOLS[kind][0], NAN), c.get(_SIGCOLS[kind][1], NAN)
    return s / m if (np.isfinite(s) and np.isfinite(m) and m > 0) else NAN


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
def mean_se(x):
    a = np.asarray([v for v in x if np.isfinite(v)], dtype=float)
    if len(a) == 0:
        return NAN, NAN, 0
    se = a.std(ddof=1) / math.sqrt(len(a)) if len(a) >= 2 else NAN
    return float(a.mean()), float(se), len(a)


def ref_map(tree: dict, base: dict) -> dict:
    """(run, split) -> record of the reference config (out first, then DIR0)."""
    m = dict(base.get(REF_TAG, {}))
    m.update(tree.get(REF_TAG, {}))
    return m


def collect(tree, tag, rtname, half, refmap, splits=None) -> list[dict]:
    rows = []
    for key in (splits if splits is not None else HALF_SPLITS[rtname][half]):
        rec = tree.get(tag, {}).get(key)
        ref = refmap.get(key)
        if rec is None or ref is None:
            continue
        rows.append({
            "key": key,
            "s_new": sigma(rec, "new"), "s_leg": sigma(rec, "legacy"),
            "s_ref": sigma(ref, "legacy"),
            "n_new": rec["head"]["new"].get("n_events", NAN),
            "n_ref": ref["head"]["legacy"].get("n_events", NAN),
            "pk_new": rec["hdr"]["peak"]["new"], "pk_leg": rec["hdr"]["peak"]["legacy"],
            "dE": rec["hdr"]["dE_frac_mean"],
            "n2_leg": rec["hdr"]["cc"]["legacy"][2], "n2_new": rec["hdr"]["cc"]["new"][2],
        })
    return rows


def stats_of(rows) -> dict:
    """Aggregate per-split rows (sigma/E in fractions; *_pp in percentage points)."""
    d = [r["s_new"] - r["s_ref"] for r in rows]
    dg = [r["s_new"] - r["s_leg"] for r in rows]
    dm, dse, n = mean_se(d)
    gm, gse, _ = mean_se(dg)
    nn = sum(r["n_new"] for r in rows)
    nr = sum(r["n_ref"] for r in rows)
    pn = sum(r["pk_new"] for r in rows)
    pl = sum(r["pk_leg"] for r in rows)
    return {
        "n": len(rows), "n_valid": n,
        "d_pp": dm * 100, "d_se_pp": dse * 100,
        "dg_pp": gm * 100, "dg_se_pp": gse * 100,
        "s_new": mean_se([r["s_new"] for r in rows])[0],
        "s_leg": mean_se([r["s_leg"] for r in rows])[0],
        "s_ref": mean_se([r["s_ref"] for r in rows])[0],
        "s_ref_se": mean_se([r["s_ref"] for r in rows])[1],
        "eps_yield": nn / nr if (rows and nr > 0) else NAN,
        "eps_peak": pn / pl if (rows and np.isfinite(pl) and pl > 0) else NAN,
        "dE_mean": mean_se([r["dE"] for r in rows])[0],
        "n2_leg": sum(r["n2_leg"] for r in rows), "n2_new": sum(r["n2_new"] for r in rows),
    }


def pooled(per_rt: dict, key: str, se_key: str):
    """Equal-weight mean over run types; SE = sqrt(sum SE^2)/n_rt."""
    v = [per_rt[n][key] for n in RT_NAMES if n in per_rt]
    s = [per_rt[n][se_key] for n in RT_NAMES if n in per_rt]
    if len(v) != len(RT_NAMES) or not all(np.isfinite(v)):
        return NAN, NAN
    se = math.sqrt(sum(x * x for x in s)) / len(s) if all(np.isfinite(s)) else NAN
    return float(np.mean(v)), se


def fm(x, nd=3, pct=False) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x * 100:.{nd}f}" if pct else f"{x:.{nd}f}"


def pm(m, se, nd=3) -> str:
    return f"{fm(m, nd)} +- {fm(se, nd)}"


# ---------------------------------------------------------------------------
# Selection (plan 3.5), half A
# ---------------------------------------------------------------------------
def pre_width(cfg) -> float:
    """Pre-window width HI-LO; pmap counts as the mean per-run map width over
    the representative runs (only matters when pmap competes with explicit
    windows; between pmap configs it is a constant)."""
    if cfg["LO"] is not None:
        return cfg["HI"] - cfg["LO"]
    w = [(_PRE_WINDOW_MAP.get(r, _DEFAULT_PRE_WINDOW)[1]
          - _PRE_WINDOW_MAP.get(r, _DEFAULT_PRE_WINDOW)[0]) for r in RT_OF_RUN]
    return float(np.mean(w))


def evaluate(tree, cfgs, half, refmap, trace):
    usable, missing = [], []
    for c in cfgs:
        have = tree.get(c["tag"], {})
        if all(any(k in have for k in HALF_SPLITS[n][half]) for n in RT_NAMES):
            usable.append(c)
        else:
            missing.append(c["tag"])
            trace.append(f"EXCLUDED {c['tag']}: no result for at least one run type")
    if not usable:
        raise SystemExit("no usable candidate config has results for every run type")
    common = {}
    for n in RT_NAMES:
        s = set(HALF_SPLITS[n][half]) & set(refmap)
        for c in usable:
            s &= set(tree[c["tag"]])
        if not s:
            raise SystemExit(f"run type {n}: no split common to all candidates and the reference")
        common[n] = sorted(s)
        nominal = len(HALF_SPLITS[n][half])
        if len(s) < nominal:
            trace.append(f"{n}: {len(s)}/{nominal} splits common to all candidates + reference")
    evals = []
    for c in usable:
        per = {n: stats_of(collect(tree, c["tag"], n, half, refmap, common[n]))
               for n in RT_NAMES}
        md, se = pooled(per, "d_pp", "d_se_pp")
        fails = []
        for n in RT_NAMES:
            for k in ("eps_yield", "eps_peak"):
                v = per[n][k]
                if not np.isfinite(v):
                    fails.append(f"{n}: {k}=NaN")
                elif v < EPS_MIN:
                    fails.append(f"{n}: {k}={v:.4f}")
        evals.append({"cfg": c, "tag": c["tag"], "per_rt": per, "mean_d_pp": md,
                      "se_pp": se, "eps_ok": not fails, "fails": fails,
                      "eps_yield_min": min(per[n]["eps_yield"] for n in RT_NAMES),
                      "eps_peak_min": min(per[n]["eps_peak"] for n in RT_NAMES)})
    return evals, missing, common


def choose(evals, trace):
    """Plan 3.5 steps 1-4.  Returns (winner eval | None, in_band tag set)."""
    cands = [e for e in evals if e["cfg"]["W"] < GATE_OFF_W]
    trace.append(f"step 1: {len(cands)} candidates after dropping W>={_n(GATE_OFF_W)} reference")
    ok = [e for e in cands if e["eps_ok"]]
    for e in cands:
        if not e["eps_ok"]:
            trace.append(f"step 2: FAIL {e['tag']}: " + "; ".join(e["fails"]))
    trace.append(f"step 2: {len(ok)} candidates keep eps_yield,eps_peak >= {EPS_MIN} in every run type")
    ok = [e for e in ok if np.isfinite(e["mean_d_pp"])]
    if not ok:
        return None, set()
    best = min(ok, key=lambda e: e["mean_d_pp"])
    se = best["se_pp"] if np.isfinite(best["se_pp"]) else 0.0
    thr = best["mean_d_pp"] + se
    trace.append(f"step 3: min mean Delta = {best['mean_d_pp']:.4f} pp ({best['tag']}), "
                 f"pooled SE {best['se_pp']:.4f} pp -> band <= {thr:.4f} pp")
    band = [e for e in ok if e["mean_d_pp"] <= thr + 1e-12]
    win = sorted(band, key=lambda e: (-e["cfg"]["W"], -pre_width(e["cfg"]),
                                      e["mean_d_pp"]))[0]
    trace.append("step 4: in band: " + ", ".join(e["tag"] for e in band)
                 + f"; largest W, then widest pre-window, then lowest Delta -> {win['tag']}")
    return win, {e["tag"] for e in band}


def cmd_select(a) -> int:
    out, base = Path(a.out), (Path(a.baseline) if a.baseline else None)
    ovr = parse_overrides(a.select_override)
    sel = lambda st: get_selection(st, out, base, ovr)          # noqa: E731
    cfgs = (dedup(load_configs_file(a.configs_file)) if a.configs_file
            else select_candidates(a.stage, sel))
    tree = load_tree(out)
    if base:
        for tag, d in load_tree(base).items():
            tree.setdefault(tag, d)
    refmap = ref_map(tree, {})
    if not refmap:
        raise SystemExit(f"reference config {REF_TAG} has no results under {out}"
                         + (f" or {base}" if base else ""))
    trace: list[str] = []
    evals, missing, common = evaluate(tree, cfgs, "A", refmap, trace)
    win, band = choose(evals, trace)

    optima = {}
    for n in RT_NAMES:
        c = [e for e in evals if e["cfg"]["W"] < GATE_OFF_W
             and np.isfinite(e["per_rt"][n]["d_pp"])]
        if c:
            b = min(c, key=lambda e: e["per_rt"][n]["d_pp"])
            p = b["per_rt"][n]
            optima[n] = {"tag": b["tag"], "d_pp": p["d_pp"], "se_pp": p["d_se_pp"],
                         "eps_yield": p["eps_yield"], "eps_peak": p["eps_peak"]}

    # readable table
    print(f"\nSelection stage {a.stage}, half A, reference legacy of {REF_TAG}; "
          f"splits per run type: " + ", ".join(f"{n}={len(common[n])}" for n in RT_NAMES))
    print(f"{'config':<22}{'mean D(pp)':>12}{'SE':>8}{'epsY_min':>10}{'epsP_min':>10}  flags")
    for e in sorted(evals, key=lambda e: (np.nan_to_num(e["mean_d_pp"], nan=9e9))):
        flags = []
        if e["cfg"]["W"] >= GATE_OFF_W:
            flags.append("gate-off ref")
        if not e["eps_ok"]:
            flags.append("EPS-FAIL")
        if e["tag"] in band:
            flags.append("in-band")
        if win and e["tag"] == win["tag"]:
            flags.append("<== WINNER")
        print(f"{e['tag']:<22}{fm(e['mean_d_pp'], 4):>12}{fm(e['se_pp'], 4):>8}"
              f"{fm(e['eps_yield_min'], 4):>10}{fm(e['eps_peak_min'], 4):>10}  "
              + " ".join(flags))
    print("\nPer-run-type optimum (diagnostic; min Delta in pp):")
    for n, o in optima.items():
        print(f"  {n:<18}{o['tag']:<22}{o['d_pp']:+.4f} +- {fm(o['se_pp'], 4)}  "
              f"epsY={fm(o['eps_yield'], 4)} epsP={fm(o['eps_peak'], 4)}")
    print("\nRule trace:")
    for t in trace:
        print("  " + t)

    def clean(o):
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, float) and not np.isfinite(o):
            return None
        return o

    doc = {
        "stage": a.stage, "half": "A", "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "reference": REF_TAG,
        "rule": {"eps_min": EPS_MIN, "band": "min mean Delta + pooled SE of the minimum",
                 "tie_break": "largest W, widest pre-window (HI-LO; pmap = mean per-run width), lowest mean Delta"},
        "winner": ({k: win["cfg"][k] for k in ("tag", "time_calib", "W", "LO", "HI")}
                   if win else None),
        "n_splits": {n: len(common[n]) for n in RT_NAMES},
        "candidates": [{
            "tag": e["tag"], "W": e["cfg"]["W"], "LO": e["cfg"]["LO"], "HI": e["cfg"]["HI"],
            "time_calib": e["cfg"]["time_calib"], "gate_off_reference": e["cfg"]["W"] >= GATE_OFF_W,
            "mean_d_pp": e["mean_d_pp"], "se_pp": e["se_pp"], "eps_ok": e["eps_ok"],
            "eps_fails": e["fails"], "in_band": e["tag"] in band,
            "per_run_type": {n: {k: e["per_rt"][n][k] for k in
                                 ("n", "d_pp", "d_se_pp", "dg_pp", "dg_se_pp", "eps_yield", "eps_peak")}
                             for n in RT_NAMES}} for e in evals],
        "excluded_missing": missing, "run_type_optima": optima, "trace": trace,
    }
    out.mkdir(parents=True, exist_ok=True)
    p = sel_path(out, a.stage)
    p.write_text(json.dumps(clean(doc), indent=2))
    if not win:
        print(f"\nERROR: no candidate passes the efficiency gate; wrote {p} with winner=null",
              file=sys.stderr)
        return 1
    print(f"\nWinner: {win['tag']}  -> {p}")
    return 0


# ---------------------------------------------------------------------------
# run: pool, manifest, failures
# ---------------------------------------------------------------------------
def git_info() -> dict:
    try:
        head = subprocess.check_output(["git", "-C", str(_REPO_ROOT), "rev-parse", "HEAD"],
                                       text=True, stderr=subprocess.DEVNULL).strip()
        dirty = bool(subprocess.check_output(["git", "-C", str(_REPO_ROOT), "status", "--porcelain"],
                                             text=True, stderr=subprocess.DEVNULL).strip())
        return {"head": head, "dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"head": None, "dirty": None}


def load_manifest(out: Path) -> dict:
    p = out / "manifest.json"
    if p.exists():
        try:
            return json.loads(p.read_text())
        except json.JSONDecodeError:
            pass
    return {}


def save_manifest(out: Path, m: dict):
    (out / "manifest.json").write_text(json.dumps(m, indent=2))


def job_done(job, prev_args) -> tuple[bool, str]:
    """Is this job's result.tsv complete and consistent with the job?"""
    tsv = job["tsv"]
    if not tsv.exists():
        return False, "no result.tsv"
    try:
        t = parse_tsv(tsv)
        h = parse_header(tsv)
    except Exception as exc:                      # noqa: BLE001
        return False, f"unparseable ({exc})"
    if not _complete(t):
        return False, "headline/calibration block incomplete"
    if not (np.isclose(h["pre"][0], job["lo"]) and np.isclose(h["pre"][1], job["hi"])):
        return False, f"pre_window {h['pre']} != ({job['lo']}, {job['hi']})"
    if not np.isclose(h["W"], job["cfg"]["W"]):
        return False, f"seed_time_window {h['W']} != {job['cfg']['W']}"
    tc_off = h["time_calib"] in (None, "off")        # old benchmark: no line = off
    if tc_off != (job["cfg"]["time_calib"] == "off"):
        return False, f"time_calib {h['time_calib']} != {job['cfg']['time_calib']}"
    if prev_args is not None and prev_args != job["bench_args"]:
        return False, "bench-args changed vs manifest"
    side = job["dir"] / "job.cmd"
    if side.exists() and side.read_text().strip() != shlex.join(job["cmd"]):
        return False, "command changed (job.cmd)"
    return True, "ok"


def kill_proc(proc):
    for sig, wait in ((signal.SIGTERM, 10.0), (signal.SIGKILL, 5.0)):
        try:
            os.killpg(proc.pid, sig)
        except (ProcessLookupError, PermissionError):
            pass
        t0 = time.monotonic()
        while proc.poll() is None and time.monotonic() - t0 < wait:
            time.sleep(0.2)
        if proc.poll() is not None:
            return


def _raise_interrupt(signum, frame):
    raise KeyboardInterrupt


def run_pool(jobs, parallel, timeout, stage, out: Path) -> tuple[int, int]:
    """Plain Popen pool, polling loop.  Returns (n_ok, n_failed)."""
    fail_path = out / "failures.tsv"
    if not fail_path.exists():
        fail_path.write_text("timestamp\tstage\tcfg\trun\tsplit\treason\trc\telapsed_s\tlog\n")

    def record_failure(j, reason, rc, el):
        row = [datetime.datetime.now().isoformat(timespec="seconds"), stage, j["cfg"]["tag"],
               str(j["run"]), f"{j['split']:05d}", reason, str(rc), f"{el:.1f}", str(j["log"])]
        with open(fail_path, "a") as fh:
            fh.write("\t".join(row) + "\n")
        print(f"FAILED {j['cfg']['tag']} run {j['run']} split {j['split']:05d}: {reason}", flush=True)

    pending = collections.deque(jobs)
    running: list = []
    total, n_ok, n_fail, t_sum = len(jobs), 0, 0, 0.0
    t_start = time.monotonic()
    old_term = signal.signal(signal.SIGTERM, _raise_interrupt)
    try:
        while pending or running:
            while pending and len(running) < parallel:
                j = pending.popleft()
                if not j["evio"].exists():
                    n_fail += 1
                    record_failure(j, "evio file missing", -1, 0.0)
                    continue
                j["dir"].mkdir(parents=True, exist_ok=True)
                (j["dir"] / "job.cmd").write_text(shlex.join(j["cmd"]) + "\n")
                lf = open(j["log"], "w")
                try:
                    proc = subprocess.Popen(j["cmd"], stdout=lf, stderr=subprocess.STDOUT,
                                            start_new_session=True)
                except OSError as exc:
                    lf.close()
                    n_fail += 1
                    record_failure(j, f"launch error: {exc}", -1, 0.0)
                    continue
                running.append((proc, j, time.monotonic(), lf))
            time.sleep(0.5)
            still = []
            for proc, j, t0, lf in running:
                rc, el = proc.poll(), time.monotonic() - t0
                if rc is None and timeout and el > timeout:
                    kill_proc(proc)
                    lf.write(f"\n[HARNESS] TIMEOUT after {el:.1f}s (limit {timeout}s)\n")
                    lf.close()
                    n_fail += 1
                    record_failure(j, "timeout", "killed", el)
                elif rc is None:
                    still.append((proc, j, t0, lf))
                    continue
                else:
                    lf.close()
                    if rc != 0:
                        n_fail += 1
                        record_failure(j, "nonzero exit", rc, el)
                    else:
                        ok, why = job_done(j, None)
                        if ok:
                            n_ok += 1
                            t_sum += el
                        else:
                            n_fail += 1
                            record_failure(j, f"bad result: {why}", rc, el)
                done = n_ok + n_fail
                eta = ((total - done) * (t_sum / max(n_ok, 1)) / parallel) if n_ok else NAN
                print(f"[{done}/{total}] ok={n_ok} failed={n_fail} running={len(running) - 1}"
                      f" elapsed={(time.monotonic() - t_start) / 60:.1f}min"
                      f" ETA={'?' if not np.isfinite(eta) else f'{eta / 60:.1f}min'}", flush=True)
            running = still
    except KeyboardInterrupt:
        print("INTERRUPTED: terminating running jobs", flush=True)
        for proc, j, t0, lf in running:
            kill_proc(proc)
            lf.close()
        raise
    finally:
        signal.signal(signal.SIGTERM, old_term)
    return n_ok, n_fail


def cmd_run(a) -> int:
    out, base = Path(a.out), (Path(a.baseline) if a.baseline else None)
    ovr = parse_overrides(a.select_override)
    sel = lambda st: get_selection(st, out, base, ovr)          # noqa: E731
    cfgs = (dedup(load_configs_file(a.configs_file)) if a.configs_file
            else stage_cfgs(a.stage, sel))
    expect = "B" if a.stage == "C" else "A"
    if a.half != expect:
        print(f"WARNING: stage {a.stage} normally runs on half {expect}", file=sys.stderr)
    jobs = plan_jobs(cfgs, a.half, out, Path(a.evio_dir), a.bench_args or "")
    prev_args = {t: c.get("bench_args") for t, c in load_manifest(out).get("cfgs", {}).items()}

    todo, skipped = [], 0
    for j in jobs:
        if not a.force:
            ok, _ = job_done(j, prev_args.get(j["cfg"]["tag"]))
            if ok:
                skipped += 1
                continue
        todo.append(j)

    if a.dry_run:
        for j in jobs:
            print(shlex.join(j["cmd"]))
        print(f"dry-run: {len(jobs)} jobs (stage {a.stage}, half {a.half}, "
              f"{len(cfgs)} configs, {skipped} already complete)")
        return 0

    for j in todo:                                   # never leave a stale result
        if j["tsv"].exists():
            j["tsv"].unlink()

    out.mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(exist_ok=True)
    man = load_manifest(out)
    man["git"] = git_info()
    man["benchmark_sha1"] = hashlib.sha1(Path(_BENCHMARK_SCRIPT).read_bytes()).hexdigest()
    man["splits"] = {n: {"A": HALF_SPLITS[n]["A"], "B": HALF_SPLITS[n]["B"]} for n in RT_NAMES}
    man.setdefault("cfgs", {})
    for c in cfgs:
        # job.cmd/header checks carry the per-job truth; this records the cfg
        extra = " ".join(x for x in (c["bench_args"], a.bench_args or "") if x).strip()
        man["cfgs"][c["tag"]] = {**c, "bench_args": extra}
    rec = {"stage": a.stage, "half": a.half, "configs": [c["tag"] for c in cfgs],
           "bench_args": a.bench_args or "", "n_jobs": len(jobs), "n_skipped": skipped,
           "parallel": a.parallel, "git": man["git"], "benchmark_sha1": man["benchmark_sha1"],
           "start": datetime.datetime.now().isoformat(timespec="seconds"),
           "end": None, "status": "running"}
    man.setdefault("stages", []).append(rec)
    save_manifest(out, man)

    print(f"stage {a.stage} half {a.half}: {len(jobs)} jobs, {skipped} already complete, "
          f"{len(todo)} to run, parallel={a.parallel}", flush=True)
    n_ok = n_fail = 0
    try:
        n_ok, n_fail = run_pool(todo, a.parallel, a.timeout_seconds, a.stage, out)
        rec["status"] = "done" if n_fail == 0 else "done_with_failures"
    except KeyboardInterrupt:
        rec["status"] = "interrupted"
        return 130
    finally:
        rec.update(end=datetime.datetime.now().isoformat(timespec="seconds"),
                   n_ok=n_ok, n_failed=n_fail)
        save_manifest(out, man)
    print(f"finished: ok={n_ok} failed={n_fail} skipped={skipped}", flush=True)
    return 1 if n_fail else 0


# ---------------------------------------------------------------------------
# summarise
# ---------------------------------------------------------------------------
def cfg_params(tag, manifest):
    p = parse_tag(tag)
    if p:
        return {"time_calib": p[0], "W": p[1], "LO": p[2], "HI": p[3]}
    c = manifest.get("cfgs", {}).get(tag)
    if c:
        return {k: c.get(k) for k in ("time_calib", "W", "LO", "HI")}
    return None


def write_tsv(path: Path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join("" if v is None else (f"{v:.6g}" if isinstance(v, float) else str(v))
                               for v in r) + "\n")


def write_jobs_tsv(tree, path: Path):
    hdr = ["cfg", "time_calib", "W", "LO", "HI", "run", "split", "run_type", "half", "path",
           "n_events", "n0", "n1", "n2plus", "sigE_G_pre_pct", "sigE_G_post_pct",
           "sigE_CB_pre_pct", "sigE_CB_post_pct", "sigE_pref_pct", "n_modules",
           "n_inner_events", "peak_yield", "dE_frac_mean", "elapsed_s"]
    rows = []
    for tag in sorted(tree):
        for (run, split), r in sorted(tree[tag].items()):
            h = r["hdr"]
            for p in ("legacy", "new"):
                c = r["cal"][p]
                cc = h["cc"][p]
                rows.append([tag, "off" if h["time_calib"] in (None, "off") else "auto",
                             h["W"], h["pre"][0], h["pre"][1], run, split, r["rt"], r["half"], p,
                             r["head"][p].get("n_events", NAN), *cc,
                             *[sigma(r, p, k) * 100 for k in ("G_pre", "G_post", "CB_pre", "CB_post")],
                             sigma(r, p) * 100, c.get("n_modules", NAN), c.get("n_inner_events", NAN),
                             h["peak"][p], h["dE_frac_mean"] if p == "new" else NAN, h["elapsed_s"]])
    write_tsv(path, hdr, rows)


def baseline_point(base_dir, base_tree):
    """DIR0's final working point: selection_B2 (else B1, else A) winner tag."""
    for st in ("B2", "B1", "A"):
        p = sel_path(base_dir, st)
        if p.exists():
            w = json.loads(p.read_text()).get("winner")
            if w and w["tag"] in base_tree:
                return w["tag"]
    return None


def write_runtype_tsv(tree, base_tree, refmap, manifest, base_tag, path: Path):
    hdr = ["cfg", "time_calib", "W", "LO", "HI", "run_type", "half", "n_splits",
           "sigE_legacy_pct", "sigE_new_pct", "sigE_ref_pct", "delta_pp", "delta_se_pp",
           "dgate_pp", "dgate_se_pp", "eps_yield", "eps_peak", "dE_frac_mean_pct",
           "n2plus_legacy", "n2plus_new",
           "dnew_vs_base_pp", "dnew_vs_base_se_pp", "dleg_vs_base_pp", "n_vs_base"]
    rows = []
    for tag in sorted(tree):
        pr = cfg_params(tag, manifest) or {}
        for n in RT_NAMES:
            for half in ("A", "B"):
                rr = collect(tree, tag, n, half, refmap)
                if not rr:
                    continue
                s = stats_of(rr)
                first = tree[tag][rr[0]["key"]]["hdr"]
                lo = "map" if pr.get("LO") is None and parse_tag(tag) else first["pre"][0]
                hi = "map" if pr.get("HI") is None and parse_tag(tag) else first["pre"][1]
                bn, bse, bnn, bleg = NAN, NAN, 0, NAN
                if base_tag:
                    bt = base_tree[base_tag]
                    dn, dl = [], []
                    for r_ in rr:
                        b = bt.get(r_["key"])
                        if b is not None:
                            dn.append(r_["s_new"] - sigma(b, "new"))
                            dl.append(r_["s_leg"] - sigma(b, "legacy"))
                    bn, bse, bnn = mean_se(dn)
                    bleg = mean_se(dl)[0]
                rows.append([tag, "off" if first["time_calib"] in (None, "off") else "auto",
                             first["W"], lo, hi, n, half, s["n"],
                             s["s_leg"] * 100, s["s_new"] * 100, s["s_ref"] * 100,
                             s["d_pp"], s["d_se_pp"], s["dg_pp"], s["dg_se_pp"],
                             s["eps_yield"], s["eps_peak"], s["dE_mean"] * 100,
                             s["n2_leg"], s["n2_new"],
                             bn * 100, bse * 100, bleg * 100, bnn])
    write_tsv(path, hdr, rows)


def determinism_check(tree, manifest):
    """Plan 3.5: legacy path must not depend on W (same time ref, same window)."""
    groups = collections.defaultdict(list)
    for tag, d in tree.items():
        ba = manifest.get("cfgs", {}).get(tag, {}).get("bench_args", "")
        for key, r in d.items():
            h = r["hdr"]
            groups[(h["time_calib"] in (None, "off"), h["pre"], ba, key)].append(r)
    viol, ngroups, ncomp = [], 0, 0
    for (toff, pre, ba, key), recs in sorted(groups.items(), key=lambda kv: str(kv[0])):
        recs = [r for r in recs]
        if len({r["hdr"]["W"] for r in recs}) < 2:
            continue
        ngroups += 1
        a0 = recs[0]
        for r in recs[1:]:
            if r["hdr"]["W"] == a0["hdr"]["W"]:
                continue
            ncomp += 1
            where = (f"run {key[0]} split {key[1]:05d} tc={'off' if toff else 'auto'} "
                     f"{a0['tag']} vs {r['tag']}")
            for fld in ("n_events", "gauss_mu_MeV", "gauss_sigma_MeV", "n_in_fit"):
                x, y = a0["head"]["legacy"].get(fld), r["head"]["legacy"].get(fld)
                if not (x == y or (np.isnan(x) and np.isnan(y))):
                    viol.append(f"{where}: legacy headline {fld} {x} != {y}")
            x, y = a0["cal"]["legacy"].get("n_inner_events"), r["cal"]["legacy"].get("n_inner_events")
            if x != y:
                viol.append(f"{where}: legacy n_inner_events {x} != {y}")
            same_t = abs(a0["hdr"]["target_mean"] - r["hdr"]["target_mean"]) < 1e-9
            for kind, fld in (("G_post", "sigma_post_MeV"), ("CB_post", "sigma_cb_MeV")):
                if same_t:
                    x, y = a0["cal"]["legacy"].get(fld), r["cal"]["legacy"].get(fld)
                    if x != y:
                        viol.append(f"{where}: legacy {fld} {x} != {y} (same target_mean)")
                else:
                    dx = abs(sigma(a0, "legacy", kind) - sigma(r, "legacy", kind)) * 100
                    if not dx <= 0.005:
                        viol.append(f"{where}: legacy sigma/E ({kind}) differs by {dx:.4f} pp "
                                    f"> 0.005 (target_mean {a0['hdr']['target_mean']} vs "
                                    f"{r['hdr']['target_mean']})")
    ok = not viol
    lines = [f"Determinism check (plan 3.5): {ngroups} (split, time-ref, window) groups with "
             f"several W, {ncomp} pairwise comparisons -> {'PASS' if ok else 'FAIL'}"]
    if ngroups == 0:
        lines.append("  (nothing to compare: no group has more than one W)")
    lines += ["  " + v for v in viol[:100]]
    if len(viol) > 100:
        lines.append(f"  ... {len(viol) - 100} more")
    return ok, "\n".join(lines)


# --- plots -----------------------------------------------------------------
C_RAW, C_CAL, C_REF, C_INK = "#2b6cb0", "#d9731f", "#6b7280", "#1f2937"


def _common_splits(tree, tags, rtname, half, refmap):
    s = set(HALF_SPLITS[rtname][half]) & set(refmap)
    for t in tags:
        s &= set(tree.get(t, {}))
    return sorted(s)


def _ref_band(ax, refmap, rtname, half, splits):
    vals = [sigma(refmap[k], "legacy") * 100 for k in splits]
    m, se, _ = mean_se(vals)
    if np.isfinite(m):
        ax.axhline(m, color=C_REF, lw=1.2, ls="-", zorder=1)
        if np.isfinite(se):
            ax.axhspan(m - se, m + se, color=C_REF, alpha=0.18, lw=0, zorder=0)
    return m


def _empty(ax, msg="no data"):
    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes, color=C_REF)
    ax.set_xticks([])
    ax.set_yticks([])


def plot_sigma_vs_W(tree, refmap, manifest, path: Path):
    stageA = {}
    for tag in tree:
        p = cfg_params(tag, manifest)
        if p and p["LO"] is None and tag in tree:
            stageA[tag] = p
    # linear W axis: scanned W at its true position; the gate-off reference
    # (W >= GATE_OFF_W) is placed in a fixed slot past the last scanned W
    ws = sorted({p["W"] for p in stageA.values() if p["W"] < GATE_OFF_W}) or [0.0]
    span = (ws[-1] - ws[0]) if len(ws) > 1 else 1.0
    x_off = ws[-1] + max(0.2 * span, 1.0)
    x_lo, x_hi = ws[0] - 0.5, x_off + 0.8
    ticks = ws + [x_off]
    tick_lab = [_n(w) for w in ws] + [f"off\n(W={_n(GATE_OFF_W)})"]
    fig, axes = plt.subplots(2, len(RT_NAMES), figsize=(4 * len(RT_NAMES), 7.2),
                             gridspec_kw={"height_ratios": [3, 1.3]}, squeeze=False)
    for j, n in enumerate(RT_NAMES):
        ax, ax2 = axes[0][j], axes[1][j]
        tags = [t for t in stageA if any(k in tree[t] for k in HALF_SPLITS[n]["A"])]
        sp = _common_splits(tree, tags, n, "A", refmap)
        ax.set_title(f"{n}  (half A, n={len(sp)})", fontsize=10)
        if not tags or not sp:
            _empty(ax)
            _empty(ax2)
            continue
        _ref_band(ax, refmap, n, "A", sp)
        for tc, col, mk in (("off", C_RAW, "o"), ("auto", C_CAL, "s")):
            pts = []
            for t in tags:
                if stageA[t]["time_calib"] == tc:
                    s = stats_of(collect(tree, t, n, "A", refmap, sp))
                    pts.append((stageA[t]["W"], s))
            pts.sort(key=lambda x: x[0])
            lab = "raw times" if tc == "off" else "offset-corrected"
            cur = [(w, s) for w, s in pts if w < GATE_OFF_W]
            off = [(x_off, s) for w, s in pts if w >= GATE_OFF_W]
            if cur:
                ax.errorbar([w for w, _ in cur], [s["s_new"] * 100 for _, s in cur],
                            yerr=[np.nan_to_num(s["d_se_pp"]) for _, s in cur],
                            color=col, marker=mk, ms=5, lw=1.4, capsize=2, label=lab, zorder=3)
                ax2.plot([w for w, _ in cur], [s["eps_yield"] * 100 for _, s in cur],
                         color=col, marker=mk, ms=4, lw=1.2, label=f"eps_yield {lab}")
                ax2.plot([w for w, _ in cur], [s["eps_peak"] * 100 for _, s in cur],
                         color=col, marker=mk, ms=4, lw=1.0, ls="--", mfc="none",
                         label=f"eps_peak {lab}")
            if off:
                ax.errorbar([w for w, _ in off], [s["s_new"] * 100 for _, s in off],
                            yerr=[np.nan_to_num(s["d_se_pp"]) for _, s in off], color=col,
                            marker="X", ms=9, mfc="none", ls="none", capsize=2, zorder=3,
                            label=f"{lab}, W={_n(GATE_OFF_W)} (gate off)")
                ax2.plot([w for w, _ in off], [s["eps_yield"] * 100 for _, s in off],
                         color=col, marker="X", ms=7, mfc="none", ls="none")
                ax2.plot([w for w, _ in off], [s["eps_peak"] * 100 for _, s in off],
                         color=col, marker="x", ms=7, ls="none")
        ax2.axhline(EPS_MIN * 100, color=C_INK, lw=0.8, ls=":")
        ax.set_ylabel("preferred post sigma/E (%)" if j == 0 else "")
        ax2.set_ylabel("efficiency (%)" if j == 0 else "")
        ax2.set_xlabel("W (ns)")
        ax.tick_params(labelbottom=False)
        for a_ in (ax, ax2):
            a_.set_xlim(x_lo, x_hi)
            a_.set_xticks(ticks)
            a_.axvline((ws[-1] + x_off) / 2, color=C_REF, lw=0.6, ls=":", zorder=0)
            a_.grid(alpha=0.25, lw=0.5)
        ax2.set_xticklabels(tick_lab, fontsize=7)
        if j == 0:
            ax.legend(fontsize=7, loc="best", title="grey: legacy of traw_w8_pmap +-SE\n"
                      "bars: SE of paired Delta", title_fontsize=7)
            ax2.legend(fontsize=6, loc="lower right", ncol=1)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def _prewindow_group(tree, manifest, base_dir, overrides):
    """(time_calib, W, LO_winner) for the pre-window scan plot."""
    expl = {t: cfg_params(t, manifest) for t in tree}
    expl = {t: p for t, p in expl.items() if p and p["LO"] is not None}
    if not expl:
        return None
    win = None
    for d in ([base_dir] if base_dir else []) + [overrides.get("_out")]:
        if d and sel_path(d, "B1").exists():
            w = json.loads(sel_path(d, "B1").read_text()).get("winner")
            if w:
                win = (w["time_calib"], w["W"], w["LO"])
                break
    if win is None:
        cnt = collections.Counter((p["time_calib"], p["W"]) for p in expl.values())
        (tc, W), _ = cnt.most_common(1)[0]
        los = sorted(p["LO"] for p in expl.values() if p["HI"] != 200.0
                     and (p["time_calib"], p["W"]) == (tc, W))
        win = (tc, W, los[0] if los else None)
    return win


def plot_prewindow(tree, refmap, manifest, base_dir, out_dir, path: Path):
    grp = _prewindow_group(tree, manifest, base_dir, {"_out": out_dir})
    fig, axes = plt.subplots(len(RT_NAMES), 2, figsize=(10, 2.6 * len(RT_NAMES)), squeeze=False)
    for j, n in enumerate(RT_NAMES):
        for k, (title, xkey) in enumerate((("LO scan (HI=200)", "LO"), ("HI scan", "HI"))):
            ax = axes[j][k]
            ax.set_title(f"{n}: {title}" if grp else n, fontsize=9)
            if grp is None:
                _empty(ax, "no explicit-window configs")
                continue
            tc, W, lo_w = grp
            sel = []
            for t in tree:
                p = cfg_params(t, manifest)
                if not p or p["LO"] is None or (p["time_calib"], p["W"]) != (tc, W):
                    continue
                if xkey == "LO" and p["HI"] == 200.0:
                    sel.append((p["LO"], t))
                if xkey == "HI" and lo_w is not None and p["LO"] == lo_w:
                    sel.append((p["HI"], t))
            sel.sort()
            tags = [t for _, t in sel]
            sp = _common_splits(tree, tags, n, "A", refmap)
            if not sel or not sp:
                _empty(ax)
                continue
            _ref_band(ax, refmap, n, "A", sp)
            xs = [x for x, _ in sel]
            ss = [stats_of(collect(tree, t, n, "A", refmap, sp)) for t in tags]
            ax.errorbar(xs, [s["s_new"] * 100 for s in ss],
                        yerr=[np.nan_to_num(s["d_se_pp"]) for s in ss], color=C_CAL if tc == "auto" else C_RAW,
                        marker="o", ms=5, capsize=2, lw=1.4, label=f"gated ({tc}, W={_n(W)})")
            ax.plot(xs, [s["s_leg"] * 100 for s in ss], color=C_INK, marker="^", ms=5, lw=1.0,
                    ls="--", label="legacy (same pre-window)")
            ax.set_xlabel(xkey + " (ns)")
            ax.set_ylabel("post sigma/E (%)")
            ax.grid(alpha=0.25, lw=0.5)
            if j == 0:
                ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def c_configs(tree, out, base_dir, overrides):
    """[(label, tag)] for C1..C4, unresolvable ones omitted; dedup merges labels."""
    items = [("C1", REF_TAG)]
    for lab, st in (("C2", "A"), ("C3", "B2")):
        try:
            items.append((lab, get_selection(st, out, base_dir, overrides)["tag"]))
        except SystemExit:
            pass
    items.append(("C4", "tcal_w4_p130-200"))
    merged: dict = {}
    for lab, tag in items:
        merged.setdefault(tag, []).append(lab)
    return [("/".join(v), t) for t, v in merged.items()]


def plot_eval(tree, refmap, manifest, out, base_dir, overrides, path: Path):
    cl = [(lab, t) for lab, t in c_configs(tree, out, base_dir, overrides) if t in tree]
    fig, axes = plt.subplots(1, len(RT_NAMES), figsize=(4 * len(RT_NAMES), 4), squeeze=False)
    for j, n in enumerate(RT_NAMES):
        ax = axes[0][j]
        ax.set_title(f"{n} (half B)", fontsize=10)
        have = [(lab, t) for lab, t in cl if collect(tree, t, n, "B", refmap)]
        if not have:
            _empty(ax)
            continue
        spB = _common_splits(tree, [t for _, t in have], n, "B", refmap)
        if spB:
            _ref_band(ax, refmap, n, "B", spB)
        for x, (lab, t) in enumerate(have):
            sB = stats_of(collect(tree, t, n, "B", refmap, spB or None))
            ax.errorbar([x], [sB["s_new"] * 100], yerr=[np.nan_to_num(sB["d_se_pp"])],
                        color=C_CAL, marker="o", ms=6, capsize=3, ls="none",
                        label="half B (new)" if x == 0 else None)
            rA = collect(tree, t, n, "A", refmap)
            if rA:
                sA = stats_of(rA)
                ax.errorbar([x + 0.12], [sA["s_new"] * 100], yerr=[np.nan_to_num(sA["d_se_pp"])],
                            color=C_INK, marker="o", mfc="none", ms=6, capsize=3, ls="none",
                            label="half A (new)" if x == 0 else None)
        ax.set_xticks(range(len(have)))
        ax.set_xticklabels([f"{lab}\n{t}" for lab, t in have], fontsize=6, rotation=0)
        ax.set_ylabel("preferred post sigma/E (%)" if j == 0 else "")
        ax.grid(alpha=0.25, lw=0.5)
        if j == 0:
            ax.legend(fontsize=7, title="grey: legacy of traw_w8_pmap, half B", title_fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# --- stage-C report ----------------------------------------------------------
def stage_c_report(tree, base_tree, refmap, out, base_dir, overrides, base_tag) -> str:
    L = ["# Stage C evaluation (half B)", "",
         f"Reference: legacy row of `{REF_TAG}` on the same splits.  sigma/E in %, Delta in "
         "percentage points (pp).  Delta = new - legacy(ref); Dgate = new - legacy of the same job; "
         "SE = SD/sqrt(n) over splits.", ""]
    cl = c_configs(tree, out, base_dir, overrides)
    present = [(lab, t) for lab, t in cl if t in tree and any(
        collect(tree, t, n, "B", refmap) for n in RT_NAMES)]
    for lab, t in cl:
        if (lab, t) not in present:
            L.append(f"- {lab} `{t}`: no half-B results")
    L.append("")
    if not present:
        return "\n".join(L + ["No stage-C results found."])
    head = ["run type", "config", "ref sigma/E (legacy)", "new sigma/E", "Dgate +- SE (pp)",
            "Delta +- SE (pp)", "eps_yield", "eps_peak", "n"]
    if base_tag:
        head.append(f"Delta vs base `{base_tag}` (pp)")
    L.append("| " + " | ".join(head) + " |")
    L.append("|" + "---|" * len(head))
    byc: dict = {}
    for lab, t in present:
        per = {}
        for n in RT_NAMES:
            rr = collect(tree, t, n, "B", refmap)
            if not rr:
                continue
            s = stats_of(rr)
            per[n] = s
            row = [n, f"{lab} `{t}`", fm(s["s_ref"], 3, True), fm(s["s_new"], 3, True),
                   pm(s["dg_pp"], s["dg_se_pp"]), pm(s["d_pp"], s["d_se_pp"]),
                   fm(s["eps_yield"], 4), fm(s["eps_peak"], 4), str(s["n"])]
            if base_tag:
                dn = [r["s_new"] - sigma(base_tree[base_tag][r["key"]], "new")
                      for r in rr if r["key"] in base_tree[base_tag]]
                m, se, nb = mean_se(dn)
                row.append(f"{pm(m * 100, se * 100)} (n={nb})")
            L.append("| " + " | ".join(row) + " |")
        byc[(lab, t)] = per
        if len(per) == len(RT_NAMES):
            dg, dgse = pooled(per, "dg_pp", "dg_se_pp")
            d, dse = pooled(per, "d_pp", "d_se_pp")
            row = ["**mean of 5**", f"{lab} `{t}`", "", "", pm(dg, dgse), pm(d, dse),
                   fm(min(p["eps_yield"] for p in per.values()), 4) + " (min)",
                   fm(min(p["eps_peak"] for p in per.values()), 4) + " (min)", ""]
            if base_tag:
                row.append("")
            L.append("| " + " | ".join(row) + " |")
    L += ["", "## Verdict: does the gate improve sigma/E? (Dgate_B < 0 by more than 2 SE)", "",
          "| config | run type | Dgate_B +- SE (pp) | verdict |", "|---|---|---|---|"]
    for (lab, t), per in byc.items():
        for n, s in list(per.items()):
            L.append(f"| {lab} `{t}` | {n} | {pm(s['dg_pp'], s['dg_se_pp'])} | "
                     f"{_verdict(s['dg_pp'], s['dg_se_pp'])} |")
        if len(per) == len(RT_NAMES):
            m, se = pooled(per, "dg_pp", "dg_se_pp")
            L.append(f"| {lab} `{t}` | **mean of 5** | {pm(m, se)} | {_verdict(m, se)} |")

    # over-tuning: final (C3) winner, A vs B
    fin = [(lab, t) for lab, t in present if "C3" in lab.split("/")]
    L += ["", "## Over-tuning guard (final winner, half A vs half B)", ""]
    if not fin:
        L.append("C3 (final winner) unavailable: selection_B2 missing or no half-B data.")
    else:
        lab, t = fin[0]
        L += [f"Config {lab} `{t}`.  FLAG if |Delta_A - Delta_B| > 2 sqrt(SE_A^2 + SE_B^2).", "",
              "| run type | Delta_A +- SE | Delta_B +- SE | |A-B| | 2 x combined SE | flag |",
              "|---|---|---|---|---|---|"]
        perA, perB = {}, byc[fin[0]]
        for n in RT_NAMES:
            rA = collect(tree, t, n, "A", refmap)
            if rA:
                perA[n] = stats_of(rA)
        rows = [(n, perA[n]["d_pp"], perA[n]["d_se_pp"], perB[n]["d_pp"], perB[n]["d_se_pp"])
                for n in RT_NAMES if n in perA and n in perB]
        if len(perA) == len(RT_NAMES) and len(perB) == len(RT_NAMES):
            ma, sa = pooled(perA, "d_pp", "d_se_pp")
            mb, sb = pooled(perB, "d_pp", "d_se_pp")
            rows.append(("**mean of 5**", ma, sa, mb, sb))
        for n, ma, sa, mb, sb in rows:
            comb = math.sqrt(sa ** 2 + sb ** 2) if np.isfinite(sa) and np.isfinite(sb) else NAN
            diff = abs(ma - mb)
            flag = "n/a" if not np.isfinite(comb) else ("**FLAG**" if diff > 2 * comb else "ok")
            L.append(f"| {n} | {pm(ma, sa)} | {pm(mb, sb)} | {fm(diff)} | {fm(2 * comb)} | {flag} |")

    # gate-off reference, half A
    off = [t for t in tree if (cfg_params(t, {}) or {}).get("W", 0) >= GATE_OFF_W
           and (cfg_params(t, {}) or {}).get("LO") is None]
    if off:
        L += ["", "## Gate-off reference (W=100, half A) against legacy", "",
              "| config | run type | Dgate +- SE (pp) | Delta +- SE (pp) |", "|---|---|---|---|"]
        for t in sorted(off):
            for n in RT_NAMES:
                rr = collect(tree, t, n, "A", refmap)
                if rr:
                    s = stats_of(rr)
                    L.append(f"| `{t}` | {n} | {pm(s['dg_pp'], s['dg_se_pp'])} | "
                             f"{pm(s['d_pp'], s['d_se_pp'])} |")
    return "\n".join(L) + "\n"


def _verdict(m, se) -> str:
    if not (np.isfinite(m) and np.isfinite(se)):
        return "undetermined"
    if m < -2 * se:
        return "improves"
    if m > 2 * se:
        return "WORSE"
    return "no significant effect"


def cmd_summarise(a) -> int:
    out = Path(a.out)
    base_dir = Path(a.baseline) if a.baseline else None
    ovr = parse_overrides(a.select_override)
    tree = load_tree(out)
    base_tree = load_tree(base_dir) if base_dir else {}
    if not tree:
        raise SystemExit(f"no result.tsv found under {out}")
    manifest = load_manifest(out)
    refmap = ref_map(tree, base_tree)
    if not refmap:
        print(f"WARNING: reference {REF_TAG} has no results; Delta columns will be NA",
              file=sys.stderr)
    base_tag = baseline_point(base_dir, base_tree) if base_dir else None
    if base_dir:
        print(f"baseline working point: {base_tag}" if base_tag
              else f"WARNING: no usable selection json in {base_dir}", flush=True)
    n_rec = sum(len(v) for v in tree.values())
    print(f"loaded {n_rec} results in {len(tree)} configs from {out}", flush=True)

    write_jobs_tsv(tree, out / "gate_scan_jobs.tsv")
    write_runtype_tsv(tree, base_tree, refmap, manifest, base_tag, out / "gate_scan_runtype.tsv")
    ok, text = determinism_check(tree, manifest)
    (out / "determinism_check.txt").write_text(text + "\n")
    print(text, flush=True)
    plot_sigma_vs_W(tree, refmap, manifest, out / "gate_scan_sigma_vs_W.png")
    plot_prewindow(tree, refmap, manifest, base_dir, out, out / "gate_scan_sigma_vs_prewindow.png")
    plot_eval(tree, refmap, manifest, out, base_dir, ovr, out / "gate_scan_eval.png")
    md = stage_c_report(tree, base_tree, refmap, out, base_dir, ovr, base_tag)
    (out / "gate_scan_stageC.md").write_text(md)
    print("\n" + md)
    print(f"wrote gate_scan_jobs.tsv, gate_scan_runtype.tsv, determinism_check.txt, "
          f"gate_scan_sigma_vs_W.png, gate_scan_sigma_vs_prewindow.png, gate_scan_eval.png, "
          f"gate_scan_stageC.md in {out}")
    return 0


# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--out", required=True, help="scan output directory")
        sp.add_argument("--baseline", metavar="DIR0",
                        help="finished scan whose selection_*.json (and results) are reused")
        sp.add_argument("--select-override", action="append", metavar="STAGE=TAG",
                        help="use TAG as the winner of STAGE (A, B1, B2) instead of its "
                             "selection json; repeatable")

    r = sub.add_parser("run", help="launch the jobs of one stage (no aggregation)")
    common(r)
    r.add_argument("--stage", required=True, choices=STAGES)
    r.add_argument("--half", required=True, choices=("A", "B"))
    r.add_argument("--configs-file", help="JSON list of {tag,time_calib,W,LO,HI,bench_args}; "
                                          "replaces the stage grid")
    r.add_argument("--bench-args", default="", help="string appended to every benchmark command")
    r.add_argument("--parallel", type=int, default=12)
    r.add_argument("--timeout-seconds", type=int, default=1800)
    r.add_argument("--evio-dir", default=_DEFAULT_EVIO_DIR)
    r.add_argument("--dry-run", action="store_true", help="print commands and the job count")
    r.add_argument("--force", action="store_true", help="rerun jobs that already completed")
    r.set_defaults(fn=cmd_run)

    s = sub.add_parser("summarise", help="tables, plots, determinism check, stage-C report")
    common(s)
    s.set_defaults(fn=cmd_summarise)

    q = sub.add_parser("select", help="apply the plan 3.5 selection rule on half A")
    common(q)
    q.add_argument("--stage", required=True, choices=SEL_STAGES)
    q.add_argument("--configs-file", help="candidates from a JSON list instead of the stage grid")
    q.set_defaults(fn=cmd_select)
    return p


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    for i, tok in enumerate(argv[:-1]):      # allow --bench-args "--flag x"
        if tok == "--bench-args":
            argv[i:i + 2] = [f"--bench-args={argv[i + 1]}"]
            break
    a = build_parser().parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
