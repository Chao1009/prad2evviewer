"""P16 shared helpers: the thr_clusters.npz dump of a report (THR_DUMP=True), module thresholds, window geometry, per-split fits.
Imported by p16_slice1.py, p16_table.py, p16_d3.py (source setup_env.sh first)."""
import math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import energy_resolution_benchmark as bm

ROOT = Path(__file__).resolve().parents[4]                       # analysis/pyscripts
P12RAW = ROOT / "output/energy_resolution/variants/P12/raw"
P16 = ROOT / "output/energy_resolution/variants/P16"
EDGE_RT = {"3p49GeV_H2_old": "3p49GeV_H2_old"}                 # edge table per run type (P12 made the three H2 energies; P16 adds 3.49 old)
TET_ADC = 10.0                                                  # daq_config.json TET above pedestal
_GEO = {}


def geometry(rt):
    """Module geometry and beam of the run type (first scan file of output/energy_resolution/scan/<rt>/)."""
    if rt not in _GEO:
        f = sorted((ROOT / "output/energy_resolution/scan" / rt).glob("*.npz"))[0]
        with np.load(f) as z:
            _GEO[rt] = bm.clusterers(int(z["run"]), str(z["calib"]), [])[1]
    return _GEO[rt]


def thresholds(rt):
    """[n modules] threshold in MeV = TET x k (median MeV per ADC of the module's stored hits, P12 edge table); NaN = no edge (dead or < 50 hits)."""
    f = P12RAW / f"edge_{rt}.npz"
    if not f.exists():
        f = P16 / "raw" / f"edge_{rt}.npz"
    if not f.exists():                                          # control v3 (2026-11-14): Ta edge table
        f = ROOT / "output/energy_resolution/control_phase2_v3_20261114/prep" / f"edge_{rt}.npz"
    return TET_ADC * np.load(f)["k_mev_per_adc"]


def load(path):
    """Dump as a dict; members as object arrays mi / me (raw energies)."""
    z = dict(np.load(path, allow_pickle=False))
    ptr = np.r_[0, np.cumsum(z["nmem"])]
    z["mi"] = bm.objects([z["mi"][a:b] for a, b in zip(ptr[:-1], ptr[1:])])
    z["me"] = bm.objects([z["me"][a:b] for a, b in zip(ptr[:-1], ptr[1:])])
    z["clip"] = np.zeros(z["e"].size, bool)
    return z


def window_d(z, geo):
    """[n, 24] distance (module sizes) of each window module centre from the cluster position (NaN where no module)."""
    m = z["win_m"].astype(np.int64)
    mm = np.where(m >= 0, m, 0)
    sz = float(np.median(geo["sx"][geo["pwo4"]]))
    d = np.hypot(geo["x"][mm] - z["x"][:, None], geo["y"][mm] - z["y"][:, None]) / sz
    return np.where(m >= 0, d, np.nan), sz


def consts(geo):
    kE = round(geo["ebeam"] / bm.E_REF, 6)
    return kE, round(math.sqrt(geo["ebeam"] / bm.E_REF), 6), bm.pool_mass("H2")


def fit(z, sel, geo, rt, e_key="e", oos=False):
    """fit_path of the clusters sel (bool) of one split, energies z[e_key], module set = the dump's used (used_oos) seeds; fields as fit_path."""
    kE, kS, mass = consts(geo)
    ep_mod = bm.e_elastic(bm.theta_deg(geo["x"], geo["y"], geo["hz"]), geo["ebeam"], mass)
    u = sel & (z["used_oos"] if oos else z["used"])
    seeds = np.unique(z["seed"][u].astype(np.int64))
    used = {int(s): float(ep_mod[s]) for s in seeds}
    if not used:
        return None
    pool = dict(seed=z["seed"][sel], e=z[e_key][sel], x=z["x"][sel], y=z["y"][sel], pu=z["pu"][sel], clip=z["clip"][sel])
    ref = float(z["ref_oos" if oos else "ref"][sel][0])
    return bm.fit_path(pool, used, float(np.mean(list(used.values()))), ref, geo, kE, kS, False, bm.SEP_CUT_MEV[rt])[0]


def mse(v):
    v = np.asarray([x for x in v if x is not None and math.isfinite(x)], float)
    return (float(v.mean()), float(v.std(ddof=1) / math.sqrt(v.size)), v.size) if v.size > 1 else (float("nan"), float("nan"), v.size)


def gains(root, rt, n_mod, names):
    """{run: [n modules] gain} of the default arm (moller 0, nonlinear 0) from <root>/<rt>/gains.tsv; modules not listed get the run's g0
    (the median of the listed prior-dominated rows)."""
    import csv
    rows = [r for r in csv.DictReader((l for l in open(f"{root}/{rt}/gains.tsv") if not l.startswith("#")), delimiter="\t")
            if r["moller"] == "0" and r["nonlinear"] == "0"]
    idx = {n: i for i, n in enumerate(names)}
    out = {}
    for run in sorted({r["run"] for r in rows}):
        rr = [r for r in rows if r["run"] == run]
        pri = [float(r["g"]) for r in rr if float(r["prior_frac"]) > 0.99]
        g = np.full(n_mod, float(np.median(pri)) if pri else 1.0)
        for r in rr:
            g[idx[r["module"]]] = float(r["g"])
        out[int(run)] = g
    return out


def egain(z, G):
    """Per cluster sum of g x me with the gains of its run (no position correction, no shape factor)."""
    out = np.empty(z["e"].size)
    for run, g in G.items():
        s = np.flatnonzero(z["run"] == run)
        out[s] = [g[m] @ w for m, w in zip(z["mi"][s], z["me"][s])]
    return out
