#!/bin/bash
# P9 variant runner: p9_run.sh <variant> <run type>   (variant "control" = no extra factor); EVAL_REGIONS on; control gate; <=2 at a time by the caller
cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh >/dev/null 2>&1 && export MPLBACKEND=Agg
V=output/energy_resolution/variants/P9; mkdir -p $V/logs
X=(--set 'PC_EXTRA=""'); [ "$1" != control ] && X=(--set "PC_EXTRA=\"$1\"")   # control = no factor (the benchmark default is now per run type, "auto")
BM=${BM:-waveform_analysis/energy_resolution_benchmark.py}
t0=$(date +%s)
python $BM report --out $V/$1 --set EVAL_REGIONS=True "${X[@]}" output/energy_resolution/scan/$2 >| $V/logs/$1_$2.log 2>&1
echo "$1 $2 rc=$? wall=$(( $(date +%s)-t0 ))s" >> $V/logs/walltimes.txt
