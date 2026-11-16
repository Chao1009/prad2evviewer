#!/usr/bin/env python3
"""p1b_scale.py <out tsv> <column> <rawpk amp_scale.tsv> <scan npz>: AMP_SCALE table of a peak column on the module calibration of the raw peak (AMP P1b).

The 8-split pass-1 gain solve of the raw peak (P1 report arm rawpk_cfd, `amp_scale.tsv`: run, module, global factor x solved gain) is the module
calibration; a column x gets scale_x(m) = scale_raw(m) * median_m(h_apk_raw / column) over the module's single-pulse hits with h_a > 200 ADC
(the median over modules for fewer than 100 hits), so that every estimator starts the report's gain solve from the same module constants.
For column h_apk_raw the table is copied unchanged."""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import energy_resolution_benchmark as bm        # noqa: E402


def main(out, col, raw_tab, scan):
    rows = bm.tsv_rows(Path(raw_tab))
    z = np.load(scan)
    names = [str(n) for n in z["mod_name"]]
    m = z["h_mod"].astype(int)
    ha, raw = z["h_a"].astype(np.float64), z["h_apk_raw"].astype(np.float64)
    c = z[col].astype(np.float64)
    ok = (z["h_k"] == 1) & ~z["h_extra"] & (ha > 200) & np.isfinite(raw) & np.isfinite(c) & (c > 0)
    ratio = {}
    for k in np.unique(m[ok]):
        s = ok & (m == k)
        if s.sum() >= 100:
            ratio[names[k]] = float(np.median(raw[s] / c[s]))
    fall = float(np.median(list(ratio.values())))
    bm.write_tsv(out, ("run", "module", "scale"),
                 [[r["run"], r["module"], float(r["scale"]) * (1.0 if col == "h_apk_raw" else ratio.get(r["module"], fall))] for r in rows])
    print(f"{col}: {len(rows)} modules, ratio raw/column median {fall:.4f}, {len(ratio)} modules with their own ratio")


if __name__ == "__main__":
    main(*sys.argv[1:])
