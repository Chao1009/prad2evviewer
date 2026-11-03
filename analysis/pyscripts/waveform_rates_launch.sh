#!/bin/bash
# Resumable launcher: waveform_good_bad.py rates over every evio file.
# usage: source setup_env.sh && export MPLBACKEND=Agg && bash waveform_rates_launch.sh [P=6]
# Existing outputs are skipped (by the script itself).  Aggregate afterwards with
#   python3 waveform_rates_aggregate.py
cd "$(dirname "$0")"
P=${1:-6}
OUT=output/waveform_rates_20261029
T=output/templates_spline_best50_20261026
bank() { case $1 in
 024386|024503|024512) echo $T/3p49GeV_H2/templates_3p49GeV_H2.npz;;
 024590|024603|024578) echo $T/0p73GeV_H2/templates_0p73GeV_H2.npz;;
 025228|025257|025260|025299|025320) echo $T/2p24GeV_H2/templates_2p24GeV_H2.npz;;
 025308) echo $T/2p24GeV_C/templates_2p24GeV_C.npz;;
 025800|026138|026230) echo output/templates_spline_ta_fix_20261028/n200/templates_2p24GeV_Ta_X17.npz;;
esac; }
JOBS=$OUT/jobs.txt; : > $JOBS
# R5 runs first, 024578 (excluded from R5 for non-standard trigger) last
for r in 024386 024503 024512 024590 024603 025228 025257 025260 025299 025308 025320 025800 026138 026230 024578; do
  for f in /Users/L00338853/work/PRad/data/evio/prad_$r.evio.*; do
    n=${f##*.}; echo "$f $OUT/pkl/${r#0}_$n.pkl $(bank $r) $OUT/logs/${r#0}_$n.log" >> $JOBS
  done
done
cat $JOBS | xargs -P "$P" -L 1 bash -c 'python3 waveform_good_bad.py rates "$0" "$1" --mf-templates "$2" > "$3" 2>&1'
