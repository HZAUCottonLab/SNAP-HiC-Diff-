# Phase 02: KNN 结果可视化 — Summary

**Completed:** 2026-05-30
**Plans:** 1 (02-01)

---

## What Was Built

### 1. `jupyter/knn_viz.py` — KNNVizPlotter 工具类

8 个 KDE/density 绘图方法 + 数据加载 + auto-discover CLI：

| Method | Description |
|--------|-------------|
| `kde_per_score()` | Univariate KDE per score (5 panels) |
| `chrom_kde_heatmap(col)` | Chromosome × -log10(p) density heatmap |
| `hexbin_pairs()` | 5×5 pairwise hexbin density matrix (subsampled 5000) |
| `knn_vs_score()` | knn_dist_median vs -log10(p) per score (hexbin, subsampled 10000) |
| `pca_contour()` | PCA scatter + KDE density contour (subsampled 15000) |
| `corr_heatmap()` | Pearson correlation heatmap |
| `bg_dev_kde()` | Background deviation KDE per score |
| `summary_table()` | Stats table figure (% p<0.05/0.01/0.001, median) |

CLI: `python3 jupyter/knn_viz.py --comparison H1398_vs_L134`
Auto-discover: `python3 jupyter/knn_viz.py`

### 2. `jupyter/knn_eda.ipynb` — 21-cell Notebook

| Section | Cell Type |
|---------|-----------|
| Header + imports | markdown + code |
| Data loading + plotter init | code |
| Sections 1–8 | markdown (header) + code (figure call) |
| Saved paths | code |

### 3. `results/figures/H1398_vs_L134/` — 24 figures

12 PNG (dpi=150) + 12 SVG for H1398_vs_L134 comparison.

## Memory Notes

- `hexbin_pairs`: 5×5 panels, subsampled to 5000 rows to avoid OOM
- `knn_vs_score`: 5 hexbin panels, subsampled to 10000 rows
- `pca_contour`: subsampled to 15000 rows
- `chrom_kde_heatmap`: uses `np.histogram` (not seaborn) — memory-safe
- `plt.close(fig)` + `plt.close("all")` after each `_save()` critical for preventing accumulation

## Key Decisions

1. **KDE density over scatter/histogram**: seaborn `kdeplot` + matplotlib `hexbin` avoid OOM that scatter plots cause with 57K rows
2. **`jupyter/` not `notebooks/`**: aligned with existing project directory (`jupyter/window_comparison_visualization.ipynb`)
3. **Subsampling thresholds**: hexbin 5000, knn_vs_score 10000, PCA 15000 — all pass end-to-end

## Artifacts

| File | Purpose |
|------|---------|
| `jupyter/knn_viz.py` | Plotter class (8 methods) |
| `jupyter/knn_eda.ipynb` | EDA Notebook (21 cells) |
| `jupyter/__init__.py` | Empty marker — resolves `from jupyter.knn_viz import` via sys.path |
| `results/figures/H1398_vs_L134/*.png` | 12 PNG figures (dpi=150) |
| `results/figures/H1398_vs_L134/*.svg` | 12 SVG figures (vector) |

---

*Phase: 02-knn-viz*
