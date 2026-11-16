#!/usr/bin/env python3
"""
farm_merge.py - P14.7: merge the per-file farm products (farm_pass.py) into the SS/out-format tables.

    python farm_merge.py --in PRODUCTS_DIR --out OUT_DIR [--steps dwell,edge,lms,mip,transfer] [--pileup w1|central|none]
           [--lms-runs RUN ...] [--edge-runs RUN ...] [--min-n N] [--rt RT ...] [--lprod-dir DIR] [--gains PATH] [--tag TAG]

PRODUCTS_DIR holds <run>/<file>.npz.  Steps (all by default; each reads only what it needs, so `--steps transfer` can be run
later on another machine against OUT_DIR from an earlier merge):
  dwell     snake_scan.blocks_to_rows (unchanged) on the stored blocks/EPICS of every scan file -> dwell_<run>.tsv, dwell_scan.tsv
  edge      per dwell module (dwell rows of all files, runs without --edge-runs excluded) the stored events of the dwell
            blocks, pile-up setting --pileup, then edge_pass.make_hist / fit_edge / second_pass unchanged
            -> edge_scan.tsv (columns of edge_23701.tsv), edge_hist/edge_<module>.npz
  lms       lms_pass.reduce_group (unchanged) over the LMS/alpha series of ALL scan runs (or --lms-runs) pooled as one epoch
            -> lms_scan[_TAG].tsv, plus lms_scan_per_run.tsv (one row per run: step detection) and lms_scan_per_run.npz
  mip       cosmic_mip.group_table / write_mip_tsv / write_lms_npz (unchanged) over the 23744 products -> mip_23744.tsv, lms_23744.npz
  transfer  transfer.build (unchanged) against edge_scan.tsv, lms_scan[_TAG].tsv and the production tables
            -> transfer_<rt>[_TAG].tsv, prior_<rt>.tsv (without a tag)
Pile-up settings (per-event peak counts of the 3x3 members are stored, so no re-decode is needed to change it):
  w1       wave 1 (default): single event = the central module has a peak and no member of the 3x3 has >= 2 peaks
  central  only the central module must have exactly one peak (pile-up in the neighbours is tolerated)
  none     no rejection: every event with a central peak (all peaks of a member summed)
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import warnings
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path

import numpy as np

FARM = Path(__file__).resolve().parent
SS = FARM.parent
PYSCRIPTS = SS.parents[1]
for p in (str(SS), str(PYSCRIPTS)):
    if p not in sys.path:
        sys.path.insert(0, p)
import snake_scan as SSN        # noqa: E402
import edge_pass as EP          # noqa: E402
import lms_pass as LP           # noqa: E402
import cosmic_mip as CM         # noqa: E402

warnings.filterwarnings("ignore", message="All-NaN slice")
LG_RUNS = (23721, 23722, 23723)         # lead-glass scan runs: the PbWO4 argmax classifier is meaningless there
SETTINGS = ("w1", "central", "none")
STEPS = ("dwell", "edge", "lms", "mip", "transfer")
MIN_N_DEFAULT = 20000                    # exact single events of a module's dwell(s) below which no edge is fitted
COSMIC_RUN = 23744


# ---------------------------------------------------------------------------------------------------------------
# index of the products
# ---------------------------------------------------------------------------------------------------------------

class Prod:
    """One product npz; keys are read on demand and the file is closed between reads (hundreds of products stay open otherwise)."""

    def __init__(self, path):
        self.path = Path(path)
        with np.load(self.path, allow_pickle=False) as z:
            self.files = list(z.files)
            self.meta = json.loads(str(z["meta_json"]))
        self.name = self.meta["file"]
        self.run = int(self.meta["run"])
        self.mode = self.meta["mode"]

    def __getitem__(self, k):
        with np.load(self.path, allow_pickle=False) as z:
            return z[k]

    def load(self, keys=None) -> dict:
        with np.load(self.path, allow_pickle=False) as z:
            return {k: z[k] for k in (keys or self.files)}

    def blocks(self):
        d = self.load(("b_counts_hi_json", "b_evfirst", "b_evlast", "b_tsfirst", "b_tslast", "b_n", "b_nhi"))
        ch = json.loads(str(d["b_counts_hi_json"]))
        evf, evl, tsf, tsl = (d[k] for k in ("b_evfirst", "b_evlast", "b_tsfirst", "b_tslast"))
        n, nhi = d["b_n"], d["b_nhi"]
        return [dict(ev_first=int(evf[i]), ev_last=int(evl[i]), ts_first=int(tsf[i]), ts_last=int(tsl[i]), n=int(n[i]),
                     n_hi=int(nhi[i]), counts={}, counts_hi={k: int(v) for k, v in ch[i].items()}) for i in range(len(n))]

    def epics(self):
        return json.loads(str(self["epics_json"]))


def index(in_dir) -> list:
    out = []
    for p in sorted(Path(in_dir).glob("*/*.npz")):
        if p.name.startswith("."):
            continue
        try:
            out.append(Prod(p))
        except Exception as ex:                                        # truncated/corrupt product: report, skip
            print(f"[merge] WARNING: cannot read {p}: {ex} (rerun that file with farm_job.sh --force)", file=sys.stderr)
    return out


def compact(files) -> str:
    """prad_023701.evio.{00000,00003} style source list, one group per run."""
    by = {}
    for f in sorted(files):
        m = re.match(r"(prad_0*\d+\.evio)\.(\d+)", f)
        if m:
            by.setdefault(m.group(1), []).append(m.group(2))
        else:
            by.setdefault(f, [])
    return " ".join(f"{k}.{{{','.join(v)}}}" if v else k for k, v in by.items())


# ---------------------------------------------------------------------------------------------------------------
# dwell
# ---------------------------------------------------------------------------------------------------------------

def step_dwell(prods, out, map_path):
    scan = [p for p in prods if p.mode == "scan"]
    rows_run = {}
    allrows = {}
    for p in scan:
        rows = SSN.blocks_to_rows(p.name, p.blocks(), p.epics())
        rows_run.setdefault(p.run, []).extend(rows)
        allrows[p.name] = rows
    for run, rows in sorted(rows_run.items()):
        SSN.write_dwell_tsv(rows, out / f"dwell_{run}.tsv", [p.name for p in scan if p.run == run], map_path)
    flat = [r for run in sorted(rows_run) for r in rows_run[run]]
    SSN.write_dwell_tsv(flat, out / "dwell_scan.tsv", [p.name for p in scan], map_path)
    nd = sum(r["kind"] == "dwell" for r in flat)
    print(f"[merge] dwell: {len(scan)} scan files, {len(flat)} rows, {nd} dwells, {len(rows_run)} runs -> dwell_<run>.tsv, dwell_scan.tsv")
    return flat


# ---------------------------------------------------------------------------------------------------------------
# edge
# ---------------------------------------------------------------------------------------------------------------

class Cache:
    def __init__(self, n=3):
        self.d, self.n = OrderedDict(), n

    def get(self, p: Prod):
        if p.name in self.d:
            self.d.move_to_end(p.name)
            return self.d[p.name]
        v = dict(evn=p["ev_evn"].astype(np.int64), ctr=p["ev_ctr"], I=p["ev_I"], npk=p["ev_npk"], cnt=p["cnt"],
                 bf=p["b_evfirst"], bl=p["b_evlast"])
        self.d[p.name] = v
        while len(self.d) > self.n:
            self.d.popitem(last=False)
        return v


def gather(module, ranges, mid, members, cache, setting):
    """ranges: [(Prod, ev_first, ev_last)] -> dict(module, members, I, n_pileup, n_exact_single, n_cpk, keep_every)."""
    Is, n_cpk, n_w1, n_cen, ke = [], 0, 0, 0, set()
    c = members.index(module)
    for p, ef, el in ranges:
        d = cache.get(p)
        m = (d["evn"] >= ef) & (d["evn"] <= el) & (d["ctr"] == mid[module])
        I, npk = d["I"][m][:, :len(members)], d["npk"][m][:, :len(members)]
        if setting == "w1":
            I = I[~(npk >= 2).any(axis=1)]
        elif setting == "central":
            I = I[npk[:, c] == 1]
        Is.append(I)
        bsel = np.where((d["bf"] >= ef) & (d["bl"] <= el))[0]
        cn = d["cnt"]
        cn = cn[np.isin(cn[:, 0], bsel) & (cn[:, 1] == mid[module])]
        n_cpk += int(cn[:, 2].sum()); n_w1 += int(cn[:, 3].sum()); n_cen += int(cn[:, 4].sum())
        ke.add(int(p.meta.get("keep_every", 1)))
    n_single = {"w1": n_w1, "central": n_cen, "none": n_cpk}[setting]
    return dict(module=module, members=members, I=np.concatenate(Is) if Is else np.zeros((0, len(members)), np.float32),
                n_pileup=0, n_exact_single=n_single, n_exact_pileup=n_cpk - n_single if setting != "none" else 0,
                n_cpk=n_cpk, keep_every=max(ke) if ke else 1)


def safe(fn, *a):
    try:
        return fn(*a), ""
    except Exception as ex:
        return None, f"{type(ex).__name__}: {ex}"


def step_edge(prods, rows, out, map_path, setting, min_n, edge_runs):
    byname = {p.name: p for p in prods if p.mode == "scan"}
    nbm = EP.neighbours(map_path)
    first = next(iter(byname.values()))
    names = [str(s) for s in first["mod_names"]]
    mid = {n: i for i, n in enumerate(names)}
    calibs = {p.meta["calib"] for p in byname.values()}
    if len(calibs) != 1:
        print(f"[merge] WARNING: several calibration files in the products: {sorted(calibs)}; using the first", file=sys.stderr)
    calib = sorted(calibs)[0]
    W = EP.load_calib(calib)
    mods = {}
    for r in rows:
        if r["kind"] != "dwell" or r["module"] not in nbm or r["file"] not in byname:
            continue
        p = byname[r["file"]]
        if edge_runs is not None and p.run not in edge_runs:
            continue
        mods.setdefault(r["module"], []).append((p, r["ev_first"], r["ev_last"]))
    order = sorted(mods, key=lambda m: (mods[m][0][0].run, mods[m][0][0].name, m))
    cache = Cache(3)
    fail, small, hist1, s1, exact = {}, [], {}, {}, {}
    for m in order:
        a = gather(m, mods[m], mid, nbm[m], cache, setting)
        exact[m] = (a["n_exact_single"], a["n_exact_pileup"])
        if a["n_exact_single"] < min_n:
            small.append((m, a["n_exact_single"]))
            continue
        h, err = safe(EP.make_hist, a, W)
        f, err = safe(EP.fit_edge, h) if h is not None else (None, err)
        if f is None:
            fail[m] = err
            continue
        s1[m] = f[0] / EP.E_BEAM
    print(f"[merge] edge pass 1: {len(s1)} modules fitted, {len(small)} below --min-n {min_n}, {len(fail)} fit failures")
    rows_out, hist_dir = [], out / "edge_hist"
    for m in order:
        if m not in s1:
            continue
        a = gather(m, mods[m], mid, nbm[m], cache, setting)
        scale = {n: s1[n] for n in a["members"] if n != m and n in s1}
        h = EP.make_hist(a, W, scale)
        f, err = safe(EP.fit_edge, h)
        if f is None:
            fail[m] = err
            continue
        E, sig, pl, sl, e, spread = f
        sc, _ = safe(EP.fit_edge, h, "central")
        s = E / EP.E_BEAM
        EP.save_hist(h, hist_dir / f"edge_{m}.npz")
        n1, np1 = exact[m]
        rows_out.append(dict(module=m, s=s, s_err=float(np.hypot(e / EP.E_BEAM, 0.5 * spread * s)), window_spread=spread,
                             n_single=n1, n_pileup=np1, s_central=(sc[0] / EP.E_BEAM if sc else float("nan")),
                             s_pass1=s1[m]))
    rows_out.sort(key=lambda r: int(r["module"][1:]))
    files = sorted({p.name for m in s1 for p, _, _ in mods[m]})
    ke = sorted({int(p.meta.get("keep_every", 1)) for m in s1 for p, _, _ in mods[m]})
    note = (f"; farm_merge: pile-up setting '{setting}' ({dict(w1='wave 1: no member of the 3x3 with >= 2 peaks', central='central module exactly one peak', none='no rejection')[setting]}); "
            f"n_single/n_pileup are exact counts of the dwell blocks (events with a central peak), the histograms use every {max(ke) if ke else 1}-th event "
            f"(per-event storage); modules below {min_n} single events ({len(small)}) or with failed fits ({len(fail)}) are omitted; "
            f"edge runs: {'all scan runs' if edge_runs is None else sorted(edge_runs)}")
    EP.write_edge_tsv(out / "edge_scan.tsv", rows_out, [compact(files)], map_path, calib, note)
    for m, err in fail.items():
        print(f"[merge] edge failure {m}: {err}", file=sys.stderr)
    if small:
        print("[merge] edge omitted (too few single events): " + ", ".join(f"{m}:{n}" for m, n in small[:20]) + (" ..." if len(small) > 20 else ""))
    print(f"[merge] edge: {len(rows_out)} modules -> edge_scan.tsv, edge_hist/")
    return rows_out


# ---------------------------------------------------------------------------------------------------------------
# lms
# ---------------------------------------------------------------------------------------------------------------

def lms_dict(p: Prod) -> dict:
    return {k[2:]: v for k, v in p.load([k for k in p.files if k.startswith("L_") and k not in ("L_source", "L_map")]).items()}


def _robust_z(x, err=None):
    """(x - median)/sigma, sigma = 1.48 MAD of the series (in quadrature with the point's own error when given)."""
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    med = np.nanmedian(x) if ok.any() else np.nan
    sig = 1.4826 * np.nanmedian(np.abs(x[ok] - med)) if ok.sum() > 2 else np.nan
    if err is not None:
        sig = np.hypot(sig, np.asarray(err, float)) if np.isfinite(sig) else np.full(len(x), np.nan)
    with np.errstate(all="ignore"):
        z = (x - med) / sig
    return np.where(np.isfinite(z), z, np.nan) if np.ndim(z) else np.full(len(x), np.nan)


def step_lms(prods, out, lms_runs, tag, c23744=None):
    scan = [p for p in prods if p.mode == "scan"]
    runs = sorted({p.run for p in scan})
    sel = [r for r in runs if lms_runs is None or r in lms_runs]
    files = {r: [lms_dict(p) for p in scan if p.run == r] for r in sel}
    names, L, err, n, info = LP.reduce_group([f for r in sel for f in files[r]])
    sfx = f"_{tag}" if tag else ""
    hdr = [f"farm scan epoch: runs {sel[0]}-{sel[-1]} ({len(sel)} runs: {sel}) pooled: {info['n_lms']} LMS events in {info['n_files']} files, "
           f"{info['n_alpha']} alpha events; {info['n_files_group_alpha']} files used the pooled alpha median (< {LP.MIN_ALPHA} alpha events)",
           "ASSUMPTION: the PMT HV and the LMS light are the same in all pooled runs (one scan epoch); lms_scan_per_run.tsv tests it per run "
           "(a step in L_rel_med, R_LMS2/3 or ped_med means the runs must be split)",
           "primary wave-1 L_scan = 23701 alone (--lms-runs 23701); L = median of the alpha-normalised single-peak integrals, "
           "err = 1.25 x 1.48 MAD/sqrt(n)"]
    src = np.array([f"{len(sel)}runs"] * len(names), dtype=object)
    srcfiles = compact([p.name for p in scan if p.run in sel])
    LP.write_tsv(out / f"lms_scan{sfx}.tsv", names, L, err, n, srcfiles + "; map hycal_map_til_run023849.json", hdr, src, "source")
    # per run
    rows, Lr = [], []
    a_all = LP._alpha_median([f for r in sel for f in files[r]])
    for r in sel:
        nm, Lq, eq, nq, iq = LP.reduce_group(files[r])
        Lr.append(Lq)
        ok = np.isfinite(Lq) & np.isfinite(L)
        rel = Lq[ok] / L[ok]
        al = np.concatenate([f["alpha_ref"] for f in files[r] if len(f["alpha_ref"])] or [np.empty((0, 3))])
        lr = np.concatenate([f["lms_ref"] for f in files[r] if len(f["lms_ref"])] or [np.empty((0, 3))])
        with np.errstate(all="ignore"):
            R = np.nanmedian(lr, axis=0) / (np.nanmedian(al, axis=0) if len(al) else np.nan) if len(lr) else np.full(3, np.nan)
        ped = np.nanmedian([np.nanmedian(f["ped"]) for f in files[r]]) if files[r] else np.nan
        rows.append(dict(run=r, epoch="scan", n_files=iq["n_files"], n_lms=iq["n_lms"], n_alpha=iq["n_alpha"], R_LMS2=R[1], R_LMS3=R[2],
                         ped_med=ped, n_mod=int(ok.sum()),
                         L_rel_med=float(np.median(rel)) if ok.any() else np.nan,
                         L_rel_err=float(1.4826 * np.median(np.abs(rel - np.median(rel))) / np.sqrt(ok.sum())) if ok.sum() > 1 else np.nan))
    if c23744 is not None:
        nm, Lq, eq, nq, iq = c23744
        ok = np.isfinite(Lq) & np.isfinite(L)
        rel = Lq[ok] / L[ok]
        rows.append(dict(run=COSMIC_RUN, epoch="cosmic", n_files=iq["n_files"], n_lms=iq["n_lms"], n_alpha=iq["n_alpha"],
                         R_LMS2=np.nan, R_LMS3=np.nan, ped_med=np.nan, n_mod=int(ok.sum()),
                         L_rel_med=float(np.median(rel)) if ok.any() else np.nan,
                         L_rel_err=float(1.4826 * np.median(np.abs(rel - np.median(rel))) / np.sqrt(max(ok.sum(), 1))) if ok.sum() > 1 else np.nan))
    sc = [r for r in rows if r["epoch"] == "scan"]
    for k in ("L_rel_med", "R_LMS2", "R_LMS3", "ped_med"):
        z = _robust_z([r[k] for r in sc], [r["L_rel_err"] for r in sc] if k == "L_rel_med" else None)
        for r, v in zip(sc, z):
            r["z_" + k] = float(v)
    cols = ["run", "epoch", "n_files", "n_lms", "n_alpha", "R_LMS2", "R_LMS3", "ped_med", "n_mod", "L_rel_med", "L_rel_err",
            "z_L_rel_med", "z_R_LMS2", "z_R_LMS3", "z_ped_med"]
    steps = [r["run"] for r in sc if abs(r.get("z_L_rel_med", 0) or 0) > 3 and np.isfinite(r.get("z_L_rel_med", np.nan))]
    with open(out / "lms_scan_per_run.tsv", "w") as f:
        f.write(f"# source: {srcfiles}\n")
        f.write("# per run: L_rel_med = median over modules of L_run/L_scan(pooled) (alpha-normalised single-peak LMS integral); R_LMS2/3 = "
                "median(LMS ref)/median(alpha ref) of the run; ped_med = median of the per-file median PbWO4 pedestal; z_* = (x - median "
                "over scan runs)/(1.48 MAD; z_L_rel_med also in quadrature with L_rel_err); |z| > 3 flags a step. Run 23744 (cosmic epoch) is shown for comparison only (not pooled)\n")
        f.write(f"# runs with |z_L_rel_med| > 3: {steps if steps else 'none'}\n")
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) if c in ("run", "epoch", "n_files", "n_lms", "n_alpha", "n_mod") else
                              (f"{r[c]:.5g}" if np.isfinite(r.get(c, np.nan)) else "nan") for c in cols) + "\n")
    np.savez_compressed(out / "lms_scan_per_run.npz", runs=np.array(sel), mod_names=np.array(names), L=np.array(Lr, np.float32),
                        L_pooled=L.astype(np.float32))
    print(f"[merge] lms: {int(np.isfinite(L).sum())} modules, {info['n_lms']} LMS events, {len(sel)} runs, steps(|z|>3): {steps or 'none'}"
          f" -> lms_scan{sfx}.tsv, lms_scan_per_run.tsv")
    return names, L, err, n


# ---------------------------------------------------------------------------------------------------------------
# mip (23744)
# ---------------------------------------------------------------------------------------------------------------

def cosmic_result(p: Prod) -> dict:
    z = p.load()
    counts = json.loads(str(z["counts_json"]))
    hits = np.zeros(len(z["hits_module"]), CM.HIT_DTYPE)
    hits["module"], hits["peak_adc"], hits["integral_adc"], hits["column_len"] = (z["hits_module"], z["hits_peak_adc"],
                                                                                  z["hits_integral_adc"], z["hits_column_len"])
    return dict(hits=hits, counts=counts, trig_hi={int(k): v for k, v in json.loads(str(z["trig_hi_json"])).items()},
                lms=dict(ts=z["lms_ts"], ev=z["lms_ev"], I=z["lms_I"], ref=z["lms_ref"]),
                alpha=dict(ts=z["alpha_ts"], ev=z["alpha_ev"], ref=z["alpha_ref"]),
                names=[str(s) for s in z["names"]], mod_x=z["mod_x"], mod_y=z["mod_y"], W=z["W"])


def step_mip(prods, out, map_path, in_dir):
    cos = sorted([p for p in prods if p.mode == "cosmic"], key=lambda p: p.name)
    if not cos:
        print("[merge] mip: no cosmic (23744) products")
        return None
    res = [cosmic_result(p) for p in cos]
    files = [Path(in_dir) / str(COSMIC_RUN) / p.name for p in cos]
    groups = CM.region_groups(map_path)
    rows = CM.group_table([r["hits"] for r in res], res[0]["names"], res[0]["W"], groups, files)
    stats = []
    for p, r in zip(cos, res):
        c = r["counts"]
        stats.append(f"{c.get('n_all', 0)} events, {c.get('n_tracks', 0)} tracks, {len(r['hits'])} interior hits, {c.get('n_multi', 0)} "
                     f"multi-peak, {len(r['lms']['ts'])} LMS, {c.get('n_flag_phys', 0)} flagged nch>900 without bit 24/25, {len(r['alpha']['ts'])} alpha")
    CM.write_mip_tsv(out / "mip_23744.tsv", rows, files, map_path, stats)
    CM.write_lms_npz(out / "lms_23744.npz", res, files)
    for r in rows:
        print(f"[merge] mip {r['group']:5s} mpv {r['mpv']:.1f} +- {r['err']:.1f} ADC n {r['n_hits']} modules {r['n_modules']}")
    # LMS of 23744 reduced like a group (for the per-run comparison)
    d = dict(mod_names=np.array(res[0]["names"]), lms_I=np.concatenate([r["lms"]["I"].reshape(-1, len(res[0]["names"])) for r in res]),
             lms_ref=np.concatenate([r["lms"]["ref"].reshape(-1, 3) for r in res]),
             alpha_ref=np.concatenate([r["alpha"]["ref"].reshape(-1, 3) for r in res]),
             lms_ts=np.concatenate([np.asarray(r["lms"]["ts"], np.int64) for r in res]))
    try:
        return LP.reduce_group([d])
    except Exception as ex:
        print(f"[merge] mip: 23744 LMS reduction failed: {ex}", file=sys.stderr)
        return None


# ---------------------------------------------------------------------------------------------------------------
# transfer
# ---------------------------------------------------------------------------------------------------------------

def step_transfer(out, rts, lprod_dir, gains, tag):
    import transfer as TR
    sfx = f"_{tag}" if tag else ""
    edge, lscan = out / "edge_scan.tsv", out / f"lms_scan{sfx}.tsv"
    for f in (edge, lscan):
        if not f.exists():
            print(f"[merge] transfer: {f} missing (run the edge and lms steps first)", file=sys.stderr)
            return
    for rt in rts:
        lp = Path(lprod_dir) / f"lms_prod_{rt}.tsv"
        g = Path(gains) if gains else TR.V3 / rt / "ring" / "gains.tsv"
        if not lp.exists() or not g.exists():
            print(f"[merge] transfer {rt}: skipped, missing {lp if not lp.exists() else g} "
                  f"(run `farm_merge.py --steps transfer` on the machine that has them, README step 6)")
            continue
        stage = Path(tempfile.mkdtemp(prefix="farm_merge_stage_"))
        try:
            shutil.copy(edge, stage / "edge_23701.tsv")           # transfer.build reads this fixed name
            shutil.copy(lp, stage / lp.name)
            res = TR.build(rt, g, False, lscan, stage, tag, prior=not tag, quiet=True)
            for f in stage.glob("transfer_*.tsv"):
                lines = f.read_text().splitlines()
                lines.insert(1, f"# farm_merge: 'edge edge_23701.tsv' above is the staged copy of {edge.name} (all scan runs, "
                                f"pile-up/dwell selection of that table); L_scan = {lscan.name}")
                (out / f.name).write_text("\n".join(lines) + "\n")
            for f in stage.glob("prior_*.tsv"):
                shutil.copy(f, out / f.name)
            print(f"[merge] transfer {rt}: N={res['N']:.5f} (n={res['n_N']}), pool spread {res['spread']:.4f} over {res['n_pool']} modules, "
                  f"{len(res['rows'])} modules -> transfer_{rt}{sfx}.tsv")
        finally:
            shutil.rmtree(stage, ignore_errors=True)


# ---------------------------------------------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", default="", help="products dir (<run>/<file>.npz); not needed for --steps transfer")
    ap.add_argument("--out", required=True)
    ap.add_argument("--steps", default=",".join(STEPS))
    ap.add_argument("--map", default=str(FARM / "hycal_map_til_run023849.json"))
    ap.add_argument("--pileup", choices=SETTINGS, default="w1")
    ap.add_argument("--lms-runs", type=int, nargs="*", default=None, help="runs pooled into L_scan (default: all scan runs); 23701 = wave-1 primary")
    ap.add_argument("--edge-runs", type=int, nargs="*", default=None,
                    help=f"runs whose dwells give edges (default: all scan runs except the lead-glass runs {LG_RUNS})")
    ap.add_argument("--min-n", type=int, default=MIN_N_DEFAULT)
    ap.add_argument("--rt", nargs="*", default=["2p24GeV_H2", "3p49GeV_H2"])
    ap.add_argument("--lprod-dir", default=str(SS / "out"), help="dir with lms_prod_<rt>.tsv (transfer step)")
    ap.add_argument("--gains", default=None, help="ring/gains.tsv (default transfer.V3/<rt>/ring/gains.tsv; with several --rt leave unset)")
    ap.add_argument("--tag", default="", help="suffix of lms_scan_TAG.tsv and transfer_<rt>_TAG.tsv (sensitivities)")
    a = ap.parse_args(argv)
    steps = [s for s in a.steps.split(",") if s]
    bad = [s for s in steps if s not in STEPS]
    if bad:
        ap.error(f"unknown step {bad}")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    map_path = str(Path(a.map).resolve())
    prods = index(a.inp) if a.inp and set(steps) & {"dwell", "edge", "lms", "mip"} else []
    if prods:
        part = sorted(p.name for p in prods if p.meta.get("max_events", 0))
        if part:
            print(f"[merge] WARNING: {len(part)} products were made with max_events != 0 (partial files), e.g. {part[:3]}: "
                  "their counters and arrays are incomplete; delete them and rerun with farm_job.sh --force", file=sys.stderr)
        vers = {p.meta.get("version") for p in prods}
        if len(vers) > 1:
            print(f"[merge] WARNING: products of different VERSION {sorted(vers, key=str)}: rerun the older ones (farm_job.sh --force)", file=sys.stderr)
        nscan = sum(p.mode == "scan" for p in prods)
        print(f"[merge] {len(prods)} products ({nscan} scan, {len(prods) - nscan} cosmic) in {a.inp}")
    edge_runs = None if a.edge_runs is None else set(a.edge_runs)
    if a.edge_runs is None:
        edge_runs = {p.run for p in prods if p.mode == "scan"} - set(LG_RUNS) if prods else None
    rows = []
    if "dwell" in steps or "edge" in steps:
        rows = step_dwell(prods, out, map_path)
    if "edge" in steps:
        step_edge(prods, rows, out, map_path, a.pileup, a.min_n, edge_runs)
    c44 = step_mip(prods, out, map_path, a.inp) if "mip" in steps else None
    if "lms" in steps:
        step_lms(prods, out, None if a.lms_runs is None else set(a.lms_runs), a.tag, c44)
    if "transfer" in steps:
        step_transfer(out, a.rt, a.lprod_dir, a.gains, a.tag)
    return 0


if __name__ == "__main__":
    sys.exit(main())
