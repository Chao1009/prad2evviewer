#!/bin/bash
# Status of run_exotic_full.sh: finished files (npz present), running processes, last progress line of the running logs, free disk.
PYDIR=$(cd "$(dirname "$0")/../../.." && pwd)
D=${OUT_DIR:-$PYDIR/output/pileup_rates}/exotic_full
echo "finished: $(ls "$D"/run_*_s*.npz "$D"/run_*_s*.p4.pkl 2>/dev/null | grep -v part | sed 's/\.p4\.pkl$//; s/\.npz$//' | sort -u | wc -l) / $(wc -l < "$D/files.txt" 2>/dev/null)"
echo "running: $(pgrep -f 'pileup_rates.py scan .*exotic_full' | wc -l) scan processes"
for f in "$D"/logs/run_*.log; do l=$(grep -E '\[progress\]|\[done\]|ERROR|Traceback' "$f" | tail -1); case "$l" in *"[done]"*) ;; *) echo "$(basename "$f"): $l";; esac; done
grep -l -E 'ERROR|Traceback' "$D"/logs/*.log 2>/dev/null | sed 's/^/FAILED: /'
du -sh "$D" | cut -f1; df -h "$D" | tail -1 | awk '{print "free: "$4}'
