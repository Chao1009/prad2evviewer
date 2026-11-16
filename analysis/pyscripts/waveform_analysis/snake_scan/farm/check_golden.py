#!/usr/bin/env python3
"""
check_golden.py EVIO [--write]  -- README step 1b: the bindings built on ifarm give the same numbers as the Mac build.
Runs farm_pass.scan_pass on the first 10000 physics events of the given file (prad_023701.evio.00000) with all events kept and
compares counters and sums with golden_23701_00000.json (made on the Mac with --write).  Takes about 10-20 s.
"""
import json
import sys
from pathlib import Path

import numpy as np

FARM = Path(__file__).resolve().parent
sys.path.insert(0, str(FARM))
import farm_pass as FP   # noqa: E402
import edge_pass as EP   # noqa: E402

evio = Path(sys.argv[1])
r = FP.scan_pass(evio, FP.DEFAULT_MAP.resolve(), EP.resolve_calib(FP.run_of(evio)), max_events=10000, keep_every=1)["arrays"]
cnt = r["cnt"]
g = dict(n_blocks=int(len(r["b_n"])), n_hi=int(r["b_nhi"].sum()), cnt_sum=[int(x) for x in cnt[:, 2:].sum(axis=0)],
         n_ev=int(len(r["ev_evn"])), I_sum=float(r["ev_I"].astype(np.float64).sum()), npk_sum=int(r["ev_npk"].astype(np.int64).sum()),
         ev_first=int(r["b_evfirst"][0]), ev_last=int(r["b_evlast"][-1]))
ref_path = FARM / "golden_23701_00000.json"
if "--write" in sys.argv:
    ref_path.write_text(json.dumps(g, indent=1) + "\n")
    print("wrote", ref_path, g)
    sys.exit(0)
ref = json.loads(ref_path.read_text())
bad = [k for k in ref if (abs(g[k] - ref[k]) > 1e-6 * max(1.0, abs(ref[k])) if not isinstance(ref[k], list) else g[k] != ref[k])]
print("golden:", ref)
print("here:  ", g)
print("OK: identical" if not bad else f"MISMATCH in {bad}: the bindings (bind_dec.cpp), database or map differ from the Mac build")
sys.exit(1 if bad else 0)
