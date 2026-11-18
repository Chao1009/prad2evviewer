#!/bin/bash
# Status of run_rawsum_full.sh: finished runs (npz present), running processes, last progress line per log.
PYDIR=$(cd "$(dirname "$0")/../.." && pwd)
D=${OUT_DIR:-$PYDIR/output/pileup_rates}/rawsum_full
echo "finished: $(ls "$D"/run_*.npz 2>/dev/null | wc -l) / 13"; ls "$D"/run_*.npz 2>/dev/null | xargs -n1 basename
echo "running: $(pgrep -f 'pileup_rates.py scan .*rawsum_full' | wc -l) scan processes"
for f in "$D"/logs/run_*.log; do echo "$(basename "$f"): $(grep -E '\[progress\]|\[done\]|ERROR|Traceback' "$f" | tail -1)"; done
tail -2 "$D/job.log" 2>/dev/null
