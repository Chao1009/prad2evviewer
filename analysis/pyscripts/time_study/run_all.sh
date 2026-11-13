#!/bin/bash
# usage: SCRATCH=<dir holding ct_<type>_e100/> bash run_all.sh   (run from time_study/)
D=../output/module_timing_20261030/dump
mkdir -p figures tables
{ echo "| run type, run_split | veto | N | mu (ns) | sigma (ns) | background (events/bin) |"
  echo "|---|---|---|---|---|---|"
  for a in "0p73GeV_H2 024590_00050" "2p24GeV_H2 025228_00100" "3p49GeV_H2 024503_00050"; do
    set -- $a
    python3 tcluster_minus_tveto.py $SCRATCH/ct_${1}_e100/$2/result_clusters.tsv $D/$2.npz \
        --label "$1 $2" --out figures/tcluster_minus_tveto_$1.png
  done; } | tee tables/tcluster_minus_tveto_fits.md
# pooled V1-V4 (BEST recipe, per-counter alignment, 0.2 ns bins): writes tables/pooled_vetos.md and figures/*_pooled.png
python3 tcluster_minus_tveto_pooled.py --scratch $SCRATCH
# all splits of the run (same recipe, mu_k per run or per split; needs ct_0p73GeV_H2_e100, ct_3p49GeV_H2_e100, ts_025320_e100 in SCRATCH): tables/pooled_vetos_allsplits.md, figures/*_pooled_allsplits.png
python3 tcluster_minus_tveto_allsplits.py --scratch $SCRATCH
