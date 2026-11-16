"""P11: T clusters by seed layer (0 outermost PbWO4 layer, 1 second) and lead-glass coincidence: share, median and IQR/1.349 of E/m, from t_clusters.npz.
Run from analysis/pyscripts: python3 <this file> -> improvement/P11_layers.tsv (control gate path)"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import energy_resolution_benchmark as B
O = Path("output/energy_resolution")
out = []
for rt in B.RUN_TYPES:
    z = np.load(next((O / "scan" / rt).glob("*.npz")))
    _, geo = B.clusterers(int(z["run"]), str(z["calib"]), ())
    pw = geo["pwo4"]; row, col = geo["row"], geo["col"]
    layer = np.minimum.reduce([row - row[pw].min(), row[pw].max() - row, col - col[pw].min(), col[pw].max() - col])
    d = np.load(O / "variants/P11" / rt / "t_clusters.npz"); p = f"g{B.CONTROL_GATE[rt]:g}_"
    e, lg, seed = d[p + "e"] / d[p + "m"], d[p + "lg"].astype(bool), d[p + "seed"].astype(int)
    L = layer[seed]; win = (e > 0.5) & (e < 1.3)
    for l in (0, 1):
        for name, s in (("all", L == l), ("lg", (L == l) & lg), ("nolg", (L == l) & ~lg)):
            v = e[s & win]
            q = np.percentile(v, (25, 50, 75)) if v.size > 20 else (np.nan,) * 3
            out.append((rt, l, name, int(s.sum()), 100.0 * s.sum() / max((L == l).sum(), 1), q[1], 100.0 * (q[2] - q[0]) / 1.349 / q[1]))
with open(Path("waveform_analysis/energy_resolution/improvement/P11_layers.tsv"), "w") as f:
    f.write("run_type\tlayer\tsubset\tn\tpct_of_layer\tmedian_E_over_m\tiqr_res_pct\n")
    for r in out:
        f.write("\t".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in r) + "\n")
        print(*(f"{x:.3f}" if isinstance(x, float) else x for x in r))
