# Plan 03-01 Summary: Plotting Script + Snakemake Rules

**Phase:** 03-knn-significant-heatmap
**Plan:** 01
**Status:** COMPLETE
**Wave:** 1

## Tasks Completed

### Task 1: Write scripts/plot_knn_significant_heatmap.py

**Files:** `scripts/plot_knn_significant_heatmap.py` (395 lines)

**Key decisions made during execution:**
- `--chrom` must be string (e.g., `HC04_A01`), not int
- `cooler.Cooler().matrix().fetch()` uses region string syntax: `f"{chrom}:{start}-{end}"`, not separate args
- BH-FDR implementation had a critical bug: sorted index mapping was reversed; fixed
- `knn_enriched.tsv` has 57,872 windows; 395 pass all 5 metrics FDR < 0.05 (0.68%)
- numpy thread warning suppressed with `NUMEXPR_MAX_THREADS=64`

**Verification:**
```bash
python scripts/plot_knn_significant_heatmap.py --scan-only \
    --enriched results/knn_empirical_pvalue/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv \
    --manifest results/logs/knn_significant_heatmap/H1398_vs_L134_significant_windows.tsv
# Output: H1398_vs_L134: 57872 total windows, 395 significant (FDR < 0.05)

snakemake -n results/knn_significant_heatmap/H1398_vs_L134/HC04_A01_161.png
# 395 PNG jobs + 1 checkpoint = 396 total

snakemake results/knn_significant_heatmap/H1398_vs_L134/HC04_A01_161.png
# PNG: 194K, PDF: 68K
```

### Task 2: Write rules/knn_significant_heatmap.smk

**Files:** `rules/knn_significant_heatmap.smk` (76 lines)

**Architecture decisions:**
- `checkpoint` writes manifest TSV (not .done file) — simpler and avoids path confusion
- `rule` depends on `checkpoint_manifest` as regular input (not lambda)
- `--scan-only` + `--manifest` flags for checkpoint mode
- Standard `checkpoint_manifest = OUTDIR + "/logs/knn_significant_heatmap/{comparison}_significant_windows.tsv"` input
- Snakemake 8.30.0 does not auto-expand checkpoint wildcards in `rule all`; `use rule` cannot override shell block; used standard approach with `checkpoint_manifest` input

**Verification:**
```bash
snakemake -n 2>&1 | grep "total"
# total 396  (1 checkpoint + 395 PNG jobs)
```

## Deviations from Plan

1. **checkpoint output changed** — plan used `.done` marker; actual uses manifest TSV directly (simpler)
2. **`use rule` syntax** — not usable because shell block cannot be overridden in snakemake 8.30.0; used standard approach
3. **`rule all` target** — plan referenced `.done`; actual uses `_all_knn_heatmap_pngs()` Python helper to expand all PNGs
4. **BH-FDR bug** — fixed during execution; implementation was reversed
5. **`--window-id` format** — plan specified `:06d` padding; actual uses plain integer (`{window_id}` without padding)
