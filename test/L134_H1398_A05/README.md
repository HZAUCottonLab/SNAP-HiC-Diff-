# L134 vs H1398 — A05 threshold overlap test

Goal: sweep metric-threshold combinations on `HC04_A05`, evaluate overlap with reference windows (exact match by `chrom/start/end`).

## Inputs

- Candidate windows: `results/knn_empirical_pvalue/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv`
- Reference windows: `reference/reference_windows.tsv` (parsed from `reference/*.png` filenames)
- Threshold grid: `config/threshold_grid.yaml` (186 combinations)

## Outputs

Written under `results/`:

- `threshold_overlap_summary.tsv` — 1 row per `threshold_id`
- `threshold_overlap_detail.tsv` — 1 row per (threshold_id, pass window)

Logs under `logs/`:

- `threshold_overlap_test.log`

## Run (must use dsub)

```bash
cd /share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff/test/L134_H1398_A05

dsub -s submit_dsub.sh
```

Or direct form:

```bash
dsub -q arm -aa \
  -e logs/%J.err \
  -o logs/%J.out \
  -n L134_H1398_A05_threshold_test \
  -R "cpu=10" \
  "python3 scripts/threshold_overlap_sweep.py"
```

## Notes

- Filter modes in grid: `pvalue`, `fdr_chrom` (BH-FDR within A05), `fdr_global` (BH-FDR across all chromosomes).
- Overlap definition: exact window equality (same `chrom/start/end`).
