#!/usr/bin/env python3
"""P14.1 slice-1 tests (pytest or `python test_p14_1.py`).  dwell_00/dwell_30/epics_module read
out/dwell_23701.tsv (made by `snake_scan.py path ... --out`); set P14_RERUN=1 to recompute files 00 and 30."""
import json
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import snake_scan as SSN  # noqa: E402

EVIO = Path("/Users/L00338853/work/PRad/data/evio")
TSV = HERE / "out" / "dwell_23701.tsv"


def _rows(fname):
    if os.environ.get("P14_RERUN"):
        return SSN.dwell_table(EVIO / fname, SSN.DERIVED_MAP)
    return [r for r in SSN.read_dwell_tsv(TSV) if r["file"] == fname]


def test_map_diff():
    with tempfile.TemporaryDirectory() as d:
        out = Path(d) / "m.json"
        ents = SSN.derive_map(SSN.DB / "hycal_map.json", SSN.HISTORY, out)
        diff = SSN.map_diff(SSN.DB / "hycal_map.json", json.loads(out.read_text()))
    names = sorted(n for n, _, _ in diff)
    assert names == sorted(["LMS1", "LMS2", "LMS3", "LMSP", "W155", "W188"]), names
    # the committed derived map equals a fresh derivation
    assert json.loads(SSN.DERIVED_MAP.read_text()) == json.loads(json.dumps(ents))
    new = {n: d for n, _, d in diff}
    assert new["LMS2"] == {"crate": 3, "slot": 3, "channel": 6} and new["W155"]["channel"] == 14


def test_dwell_00():
    r = [x for x in _rows("prad_023701.evio.00000") if x["kind"] == "dwell"]
    assert len(r) == 1, r
    assert r[0]["module"] == "W599" and r[0]["share"] >= 0.95 and r[0]["n_events"] >= 600_000, r[0]


def test_dwell_30():
    r = _rows("prad_023701.evio.00030")
    i = next(k for k, x in enumerate(r) if x["kind"] == "dwell")
    assert r[i]["module"] == "W585", r
    assert any(x["kind"] == "ramp" and x["module"] == "W586" for x in r[:i]), r
    assert [x for x in r if x["kind"] == "dwell"] == [r[i]]
    # ramp blocks are excluded from the dwell: the dwell starts after the ramp's last event
    assert all(x["ev_last"] < r[i]["ev_first"] for x in r[:i])


def test_epics_module():
    rows = SSN.read_dwell_tsv(TSV)
    mods = SSN.load_modules()
    aff = SSN.fit_affine(rows, mods)
    d = [x for x in rows if x["kind"] == "dwell"]
    worst = 0.0
    for x in d:
        nm, dist = SSN.ptrans_to_module(x["ptrans_x"], x["ptrans_y"], aff, mods)
        mx = aff["x"]["a"] * x["ptrans_x"] + aff["x"]["b"] - mods[x["module"]][0]
        my = aff["y"]["a"] * x["ptrans_y"] + aff["y"]["b"] - mods[x["module"]][1]
        worst = max(worst, (mx * mx + my * my) ** 0.5)
    print(f"epics_module: worst |ptrans-mapped position - dwell module centre| = {worst:.2f} mm over {len(d)} dwells")
    assert worst < 5.0, worst


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
