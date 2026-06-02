# SNAP-HiC-Diff Roadmap

## Overview

SNAP-HiC-Diff is a pipeline for differential chromatin interaction analysis between biological conditions using Hi-C data.

## Phases

| # | Phase | Status | Plans |
|---|-------|--------|-------|
| 01 | KNN Empirical P-Value | Completed | 5 |
| 02 | KNN 结果可视化 | Completed | 1 |
| 03 | KNN 显著窗口热图绘制 | Planned | 2 |

**Plans:**
- [ ] 03-01-PLAN.md — Plotting script + Snakemake rules (2 tasks, wave 1)
- [ ] 03-02-PLAN.md — Snakefile update (1 task, wave 1)

---

## Phase 01: KNN Empirical P-Value

**Goal:** Add a Snakemake rule that computes KNN-based empirical P-values for window pair metrics, using a background pool of ~2.15M windows.

**Depends on:** None

**Deliverables:**
- `config/config.yaml` — `knn_empirical_pvalue` section with K=10000, euclidean, robust scaler
- `scripts/compute_knn_empirical_pvalue.py` — KNN + empirical P-value computation
- `rules/window_compare.smk` — new `rule knn_empirical_pvalue`
- `Snakefile` `rule all` — include `{comparison}_knn_enriched.tsv` outputs

**Criteria:**
- Each comparison (H1398_vs_L134) gets an independent job
- P-values computed from MSE, d_SSIM, d_Spearman, tw_score, d_NCC (all right-tail)
- Output: `{comparison}_knn_enriched.tsv` with empirical P-values + KNN quality metrics

---

## Phase 02: KNN 结果可视化

**Goal:** 构建 Jupyter Notebook 对 KNN 经验 P 值结果进行全面的探索性数据分析，包括 P 值分布、Manhattan 图、散点矩阵、热力图、KNN 质量诊断、PCA 等 10+ 张图。

**Depends on:** Phase 01 (knn_enriched.tsv 数据)

**Deliverables:**
- `jupyter/knn_viz.py` — KNNVizPlotter 工具类，11 个绘图方法
- `jupyter/knn_eda.ipynb` — 主 EDA Notebook，调用 Plotter

**Plans:**
- [x] 02-PLAN.md — KNNVizPlotter + 主 Notebook（2 tasks，wave 1）

**Criteria:**
- jupyter/knn_viz.py: KNNVizPlotter class with 11 plotting methods + auto-discover + CLI
- jupyter/knn_eda.ipynb: 25+ cells covering setup through all 10+ figure types
- results/figures/H1398_vs_L134/: 11+ PNG files (dpi=150) and 11+ SVG files
- All 5 score metrics (MSE, d_SSIM, d_Spearman, tw_score, d_NCC) visualized

---

## Phase 03: KNN 显著窗口热图绘制

**Goal:** 新增 Snakemake rule，筛选 5 个指标 FDR < 0.05 的显著窗口，每个窗口绘制一张 2x2 热图：上面为 control/treatment iced（O/E），下面为 control/treatment gauss（Gaussian smoothed）。输出 PNG + PDF。

**Depends on:** Phase 01 (knn_enriched.tsv 数据，FDR checkpoint 内 inline 计算)

**Deliverables:**
- `rules/knn_significant_heatmap.smk` — Snakemake rule
- `scripts/plot_knn_significant_heatmap.py` — 热图绘制脚本
- `Snakefile` — include rule + rule all target

**Criteria:**
- 筛选条件：5 个 `*_empirical_pvalue` 列经 BH 校正 FDR 全部 < 0.05（checkpoint 内 inline 计算，不依赖文件预存 FDR 列）
- 布局：2x2（上行 iced / 下行 gauss；左列 control / 右列 treatment）
- 颜色条：iced 共享一个，gauss 共享一个
- 坐标轴：基因组坐标（bp），每轴最多 5 个刻度
- 区域范围：窗口自身范围
- 数据源：iced → `results/balanced/{sample}/{sample}__{chrom}_iced.cool`（per-chrom 文件，用 `cooler.fetch(chrom:start-end)` 提取）；gauss → `results/windows/{sample}/{sample}__{chrom}_{window_id}_gaussian.h5`
- 输出格式：PNG + PDF
- 输出路径：`results/knn_significant_heatmap/{comparison}/{chrom}_{window_id}.png/pdf`
