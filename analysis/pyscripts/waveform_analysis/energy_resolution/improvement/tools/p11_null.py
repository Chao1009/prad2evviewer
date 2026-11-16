"""P11 pool-null check: pool rows (kind ep, every path) of variants/P11/<RT>/regions_splits.tsv against control v2, all columns, per (run, split, path).
Run from analysis/pyscripts: python3 <this file>"""
import csv
from pathlib import Path
O = Path("output/energy_resolution")
rd = lambda p: {(r["run"], r["split"], r["path"]): r for r in csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")
                if r["kind"] == "ep" and r["region"] == "pool"}
for rt in ("0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "3p49GeV_H2_old", "2p24GeV_Ta_X17"):
    c, v = rd(O / "control_phase2_v2_20261113" / rt / "regions_splits.tsv"), rd(O / "variants/P11" / rt / "regions_splits.tsv")
    cols = [k for k in next(iter(c.values())) if k not in ("run_type", "run", "split", "path", "kind", "region")]
    bad = [(k, q) for k in c if k in v for q in cols if c[k][q] != v[k][q]]
    print(rt, "rows", len(c), len(v), "differences", len(bad), bad[:3])
