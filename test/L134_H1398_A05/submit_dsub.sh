#!/bin/bash
#DSUB -q arm
#DSUB -aa
#DSUB -n L134_H1398_A05_threshold_test
#DSUB -e logs/%J.err
#DSUB -o logs/%J.out
#DSUB -R "cpu=10"

set -euo pipefail

ROOT="/share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff/test/L134_H1398_A05"
cd "/share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff"

mkdir -p "$ROOT/logs" "$ROOT/results"

python3 "$ROOT/scripts/threshold_overlap_sweep.py" \
  --chrom "HC04_A05" \
  --comparison "H1398_vs_L134" \
  --enriched "results/knn_empirical_pvalue/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv" \
  --global-enriched "results/knn_empirical_pvalue/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv" \
  --reference "test/L134_H1398_A05/reference/reference_windows.tsv" \
  --out-dir "test/L134_H1398_A05/results" \
  --log-file "test/L134_H1398_A05/logs/threshold_overlap_test.log"
