#!/usr/bin/env python3
"""p1b_corr_table.py <out npz> <column> <scan npz...>: hit-level correction tables of a peak column against the template amplitude h_a (AMP P1b).

Single-pulse hits (h_k == 1, not extra, h_a >= 30 ADC): r = column / reference (h_a, or h_e with `column@h_e`: the table then makes the column agree with the integral energy, whose own phase correction ph(x) is inside h_e); module scale = median r over h_a > 200; the module's sampling-phase curve
(16 bins of the onset phase x = ((h_t - t50) mod 4 ns) / 4 ns, median of r / scale over h_a > 200; the pooled curve for modules with < 2000 hits);
the amplitude curve g(A) = median over all modules of r / scale / curve in 16 log bins of A.  Written: curve [n_mod, 16], g_a (bin centres, median column value, ADC),
g [bins], names.  Used by p1b_report.py (P1B_CORR) as column / (curve[module, phase] * g(A)), A = column / (median r).
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import energy_resolution_benchmark as bm        # noqa: E402

NB = 16
EDGES = np.array([30, 40, 50, 65, 80, 100, 130, 170, 220, 300, 400, 550, 750, 1000, 1400, 2000, 4000.0])


def main(out, col, *files):
    col, _, ref = col.partition("@")             # column[@reference], default reference h_a
    ref = ref or "h_a"
    zs = [np.load(f) for f in files]
    names = zs[0]["mod_name"]
    t50 = bm.template_columns(bm.PYDIR / str(zs[0]["bank"]), names)[0]
    d = {q: np.concatenate([z[q] for z in zs]) for q in ("h_mod", "h_t", "h_a", "h_e", "h_k", "h_extra", col)}
    m, ha, c = d["h_mod"].astype(int), d["h_a"].astype(np.float64), d[col].astype(np.float64)
    ref_v = d[ref].astype(np.float64)
    one = (d["h_k"] == 1) & ~d["h_extra"] & (ha >= 30) & np.isfinite(c) & np.isfinite(ha) & (ref_v > 0)
    m, ha, c, ref_v = m[one], ha[one], c[one], ref_v[one]
    x = ((d["h_t"][one].astype(np.float64) - t50[m]) % bm.DT) / bm.DT
    b = np.minimum((x * NB).astype(int), NB - 1)
    r = c / ref_v
    n_mod = names.size
    o = np.argsort(m, kind="stable")
    off = np.r_[0, np.cumsum(np.bincount(m, minlength=n_mod))]
    scale = np.full(n_mod, np.nan)
    for k in range(n_mod):
        g = o[off[k]:off[k + 1]]
        s = g[ha[g] > 200]
        if s.size > 200:
            scale[k] = np.median(r[s])
    rn = r / scale[m]
    hi = (ha > 200) & np.isfinite(rn)
    pooled = np.array([np.median(rn[hi & (b == q)]) for q in range(NB)])
    curve = np.tile(pooled, (n_mod, 1))
    for k in range(n_mod):
        g = o[off[k]:off[k + 1]]
        if g.size >= 2000:
            for q in range(NB):
                t = g[(b[g] == q) & hi[g]]
                if t.size > 20:
                    curve[k, q] = np.median(rn[t])
    rr = rn / curve[m, b]
    ok = np.isfinite(rr)
    cen, gv = [], []
    for lo, up in zip(EDGES[:-1], EDGES[1:]):
        s = ok & (ha >= lo) & (ha < up)
        cen.append(np.median(c[s]))
        gv.append(np.median(rr[s]))
    np.savez(out, curve=curve, g_a=np.array(cen), g=np.array(gv), scale=scale, names=names)
    print(f"{col}: {len(m)} hits; pooled phase curve p2p {pooled.max() - pooled.min():.4f}; g(A) {np.round(gv, 4)}")


if __name__ == "__main__":
    main(*sys.argv[1:])
