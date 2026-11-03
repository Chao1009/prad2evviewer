#!/usr/bin/env python3
"""Compare rates pickles against gb_cache (good/bad/mid/nofit per module)."""
import glob, pickle, os, sys
O = "output/waveform_rates_20261029/pkl/"
C = "output/waveform_examples_20261029/gb_cache/"
n = bad = 0
for f in sorted(glob.glob(C + "*.pkl")):
    g = O + os.path.basename(f)
    if not os.path.exists(g):
        continue
    a, b = pickle.load(open(f, "rb")), pickle.load(open(g, "rb"))
    ok = a["tot"] == b["tot"] and a["cnt"] == b["cnt"] and a["mods"].keys() == b["mods"].keys() and all(
        all(a["mods"][m][k] == b["mods"][m][k] for k in ("good", "bad", "mid", "nofit")) for m in a["mods"])
    n += 1; bad += (not ok)
    print(os.path.basename(f), "OK" if ok else "MISMATCH", a["tot"])
print(f"compared {n}, mismatches {bad}")
sys.exit(bad > 0)
