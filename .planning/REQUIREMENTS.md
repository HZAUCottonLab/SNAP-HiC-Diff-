# Requirements

## Phase 01: KNN Empirical P-Value

| ID | Requirement | Source |
|----|-------------|--------|
| REQ-01 | Each comparison (e.g. H1398_vs_L134) runs as an independent Snakemake job | User (B) |
| REQ-02 | Input: `{comparison}_window_pair_metrics.tsv` from `rules/window_compare.smk` merge output | User (A) |
| REQ-03 | Output: `{comparison}_knn_enriched.tsv` in same directory | User (D) |
| REQ-04 | KNN uses EC, IB, SP as quality-matching features (no P-value calc on these) | User |
| REQ-05 | K=10000, metric=euclidean, scaler=robust, all configurable via config.yaml | User (C) |
| REQ-06 | Score columns: MSE, d_SSIM, d_Spearman, tw_score, d_NCC | User |
| REQ-07 | d_SSIM = 1 - SSIM; d_Spearman = 1 - Spearman; tw_score = -log10(tw_pvalue); d_NCC = 1 - NCC | User |
| REQ-08 | All P-values: right-tail empirical test against KNN background | User |
| REQ-09 | No FDR correction | User |
| REQ-10 | Background file: `/share/home/xhhuang24/project/PanTE/background_pool/results/window_compare/HiC_window_background.tsv` | User |
