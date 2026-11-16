#!/usr/bin/env python3
"""P14.3 slice-3 tests (pytest or `python test_p14_3.py`).  Per-file results come from output/snake_scan/lms_npz (made by
`snake_scan.py lms FILE`; computed here when absent)."""
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lms_pass as LP  # noqa: E402
import snake_scan as SSN  # noqa: E402

EVIO = Path("/Users/L00338853/work/PRad/data/evio")
F23701 = EVIO / "prad_023701.evio.00010"
F25320 = EVIO / "prad_025320.evio.00100"


def _npz(path):
    return LP.load(LP.lms_pass(path, LP.default_map(path)))


def test_lms_count_23701():
    f = _npz(F23701)
    n = len(f["lms_ts"])
    nch = int(np.median(f["lms_nch"]))
    pk = float(np.nanmedian(f["lms_pk_med"]))
    print(f"lms_count_23701: {n} LMS events, channels per event median {nch} (range {f['lms_nch'].min()}-"
          f"{f['lms_nch'].max()}), PbWO4 median peak {pk:.0f} ADC, {np.isfinite(f['lms_I']).sum(1).min()}+ PbWO4 pulses/event")
    assert 8 <= n <= 9, n
    assert 1560 <= nch <= 1575, nch              # plan: 1569
    assert 300 <= pk <= 450, pk


def test_lms_ref_by_name():
    a = LP.ref_names_by_daq(SSN.DERIVED_MAP, str(F23701))
    b = LP.ref_names_by_daq("", str(F25320))
    print("lms_ref_by_name: 23701", a, " 25320", b)
    assert a[(3, 3, 6)] == "LMS2", a
    assert b[(3, 13, 11)] == "LMS2", b


def test_lms_prod_25320():
    f = _npz(F25320)
    n = len(f["lms_ts"])
    nmod = int(np.isfinite(f["lms_I"]).any(axis=0).sum())
    names, L, err, nn, info = LP.reduce_group([f])
    nL = int((np.isfinite(L)).sum())
    print(f"lms_prod_25320: {n} LMS events, {nmod} PbWO4 modules with a pulse, L_prod (n >= {LP.MIN_N}) for {nL}, "
          f"alpha events {len(f['alpha_ref'])}")
    assert 61 <= n <= 72, n
    assert nL >= 1140, nL


def test_segments_none():
    with tempfile.TemporaryDirectory() as td:
        steps, rows, zs = LP.reduce_segments(23701, out=Path(td) / "segments_23701.tsv")
    assert len(rows) == 9, len(rows)
    print("segments_none: steps", steps or "none")
    assert not steps, steps


if __name__ == "__main__":
    ok = True
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            try:
                fn()
                print("PASS", name)
            except Exception as e:  # report honestly
                ok = False
                print("FAIL", name, repr(e)[:400])
    sys.exit(0 if ok else 1)
