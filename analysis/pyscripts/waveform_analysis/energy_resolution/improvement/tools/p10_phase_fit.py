"""P10: fit the residual onset-phase factor from a PHASE_STUDY dump (phase_clusters.npz) and write a PHASE_RESID table.

Usage: p10_phase_fit.py <phase_clusters.npz> <scan dir> <form> <out.tsv> [--lam L] [--cv]
Model on the clusters with 0.95 < E / E' < 1.05 (the pool's control gate path, corrections applied):
  E / E' / median - 1 = sum over the three largest members j of f_j * delta_m(x_j)   (f_j: energy fraction, x_j: onset phase)
  delta = global harmonics (form G1: k = 1; G2: k = 1, 2) [+ per-module deviations with ridge weight lam: M1 = G2 + module k = 1, M2 = G2 + module k = 1, 2].
Parity: split ordinal within its run (sorted splits of the scan file) % 2.  par k rows are fitted on parity k clusters; the report
applies them to the other parity (PHASE_RESID), so the evaluation is cross-fitted.  --cv: print the held-out residual std per lam instead.
"""
import argparse, glob, sys
import numpy as np
import scipy.sparse as sp

ap = argparse.ArgumentParser()
ap.add_argument("dump"); ap.add_argument("scans"); ap.add_argument("form", choices=("G1", "G2", "M1", "M2")); ap.add_argument("out")
ap.add_argument("--lam", type=float, default=100.0); ap.add_argument("--cv", action="store_true")
ap.add_argument("--window", type=float, nargs=2, default=(0.95, 1.05))
a = ap.parse_args()
z = np.load(a.dump)
run, split = z["run"], z["split"]
par = np.zeros(run.size, int)
for r in np.unique(run):
    sps = sorted(x for f in glob.glob(f"{a.scans}/run_{int(r)}*.npz") for x in np.load(f)["splits"].tolist())   # 3.49 GeV old: several files per run
    m = run == r
    par[m] = np.array([sps.index(int(s)) % 2 for s in split[m]])
rr = z["r"]
ok = (rr > a.window[0]) & (rr < a.window[1]) & np.isfinite(z["x0"])
names = np.load(sorted(glob.glob(f"{a.scans}/run_{int(run[0])}*.npz"))[0])["mod_name"]
n_mod = len(names)
K = 1 if a.form == "G1" else 2
ng = 2 * K
MOD = {"M1": 1, "M2": 2}.get(a.form, 0)           # per-module harmonics
nd = 2 * MOD


def design(sel, use_mod):
    y = rr[sel] / np.median(rr[sel]) - 1.0
    n = sel.sum()
    rows, cols, vals = [], [], []
    G = np.zeros((n, ng))
    for j in range(3):
        xj, fj, mj = z[f"x{j}"][sel], z[f"f{j}"][sel], z[f"m{j}"][sel]
        v = np.isfinite(xj)
        th = 2 * np.pi * np.nan_to_num(xj) / 256
        for k in range(1, K + 1):
            G[:, 2 * (k - 1)] += np.where(v, fj * np.cos(k * th), 0)
            G[:, 2 * (k - 1) + 1] += np.where(v, fj * np.sin(k * th), 0)
        if use_mod:
            i = np.flatnonzero(v & (mj >= 0))
            for k in range(1, MOD + 1):
                rows += [i, i]
                cols += [ng + nd * mj[i] + 2 * (k - 1), ng + nd * mj[i] + 2 * (k - 1) + 1]
                vals += [fj[i] * np.cos(k * th[i]), fj[i] * np.sin(k * th[i])]
    A = sp.csr_matrix(G)
    if use_mod:
        M = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, ng + nd * n_mod))
        A = sp.csr_matrix(sp.hstack([sp.csr_matrix(G), sp.csr_matrix((n, nd * n_mod))])) + M
    return A, y


def solve(sel, lam):
    A, y = design(sel, MOD > 0)
    A = sp.hstack([A, sp.csr_matrix(np.ones((A.shape[0], 1)))], format="csr")       # free offset
    pen = np.zeros(A.shape[1])
    if MOD > 0:
        pen[ng:ng + nd * n_mod] = lam
    c = np.linalg.solve((A.T @ A).toarray() + np.diag(pen) + 1e-12 * np.eye(A.shape[1]), A.T @ y)
    return c[:-1], A, y


def pred(c, sel):
    A, y = design(sel, MOD > 0)
    return y - A @ c, y


if a.cv:
    for lam in ([0.0] if not MOD else [3.0, 10.0, 30.0, 100.0, 300.0, 1000.0]):
        out = []
        for k in (0, 1):
            c, _, _ = solve(ok & (par == k), lam)
            e, y = pred(c, ok & (par != k))
            out += [e.std() * 1000, y.std() * 1000]
        print(f"{a.form} lam {lam:g}: held-out resid std {np.mean(out[0::2]):.4f} (raw {np.mean(out[1::2]):.4f}) permille")
    sys.exit()
rows = []
for k in (0, 1):
    c, _, _ = solve(ok & (par == k), a.lam)
    g = np.zeros(4)
    g[:ng] = c[:ng]
    rows.append((k, "*", *g))
    if MOD > 0:
        dev = c[ng:].reshape(-1, nd)
        for m in np.flatnonzero(np.abs(dev).sum(1) > 0):
            tot = g.copy()
            tot[:nd] += dev[m]
            rows.append((k, str(names[m]), *tot))
with open(a.out, "w") as fh:
    fh.write(f"# P10 residual phase table: form {a.form}, lam {a.lam:g}, dump {a.dump}\n")
    fh.write("par\tname\tc1\ts1\tc2\ts2\n")
    for r in rows:
        fh.write("\t".join(str(x) if isinstance(x, (str, int)) else f"{x:.6g}" for x in r) + "\n")
print(a.form, "rows", len(rows), "global par0", np.round(np.array(rows[0][2:]) * 1000, 2), "permille")
