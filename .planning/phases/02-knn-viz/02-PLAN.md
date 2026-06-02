# Phase 02: KNN 结果可视化 — Plan 01

---
phase: 02-knn-viz
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - jupyter/knn_viz.py
  - jupyter/knn_eda.ipynb
autonomous: true
requirements: []
must_haves:
  truths:
    - "jupyter/knn_viz.py 包含 KNNVizPlotter 类，实现全部 10+ 绘图方法"
    - "jupyter/knn_eda.ipynb 是合法的 Jupyter Notebook (JSON)，包含 10+ 个 figure cells"
    - "执行 notebook 后 results/figures/H1398_vs_L134/ 目录包含 10+ PNG 文件"
    - "全部 5 个 score 指标 (MSE, d_SSIM, d_Spearman, tw_score, d_NCC) 均被可视化"
    - "每张图同时输出 PNG (dpi=150) 和 SVG"
  artifacts:
    - path: "jupyter/knn_viz.py"
      provides: "KNNVizPlotter 工具类，封装数据加载和 11 个绘图方法"
      exports: ["KNNVizPlotter", "auto_discover_comparisons", "load_data"]
    - path: "jupyter/knn_eda.ipynb"
      provides: "主 EDA Notebook，调用 Plotter 方法生成 10+ 张图"
      cells: "15+ cells (setup, load, plotter init, 11 figure cells, summary)"
    - path: "results/figures/H1398_vs_L134/"
      provides: "所有输出图片 PNG + SVG"
  key_links:
    - from: "jupyter/knn_eda.ipynb"
      to: "jupyter/knn_viz.py"
      via: "import knn_viz"
      pattern: "from knn_viz import"
    - from: "jupyter/knn_viz.py"
      to: "results/window_compare/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv"
      via: "load_data()"
      pattern: "knn_enriched.tsv"
---

<objective>
构建 KNN 经验 P 值 EDA 可视化工具：
1. `jupyter/knn_viz.py` — KNNVizPlotter 工具类（含 11 个绘图方法 + 数据加载 + auto-discover）
2. `jupyter/knn_eda.ipynb` — 主 Notebook，调用 Plotter 生成 10+ 张图

两个文件可以并行开发（ipynb 调用 py，故 py 先完成），但同一 wave 可并行。
</objective>

<context>
@.planning/phases/02-knn-viz/02-CONTEXT.md
@.planning/phases/02-knn-viz/02-DISCUSSION-LOG.md
@config/config.yaml
</context>

<read_first>
- `jupyter/knn_viz.py` — (new file, no read needed)
- `jupyter/window_comparison_visualization.ipynb` — existing notebook pattern reference (熟悉 jupyter/ 目录结构和 notebook JSON 格式)
- `results/window_compare/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv` — data schema reference (18 columns, 57,872 rows)
</read_first>

<tasks>

<task type="auto">
  <name>Task 1: 创建 jupyter/knn_viz.py — KNNVizPlotter 工具类</name>
  <files>jupyter/knn_viz.py</files>
  <action>
Create `jupyter/knn_viz.py` as a standalone executable script + importable module.

**Data loading:**
- `auto_discover_comparisons(base_dir)` — glob scan `results/window_compare/*/knn_enriched.tsv`, return list of comparison names
- `load_data(comparison)` — load `{comparison}_knn_enriched.tsv`, drop rows with NA/inf in SCORE_COLS or knn_dist_median

**KNNVizPlotter class:**
```python
class KNNVizPlotter:
    SCORE_COLS = [
        "MSE_empirical_pvalue", "d_SSIM_empirical_pvalue",
        "d_Spearman_empirical_pvalue", "tw_score_empirical_pvalue",
        "d_NCC_empirical_pvalue"
    ]
    SCORE_LABELS = {
        "MSE_empirical_pvalue": "MSE",
        "d_SSIM_empirical_pvalue": "d_SSIM",
        "d_Spearman_empirical_pvalue": "d_Spearman",
        "tw_score_empirical_pvalue": "tw_score",
        "d_NCC_empirical_pvalue": "d_NCC"
    }
    BG_COLS = {
        "MSE_empirical_pvalue": "bg_MSE_median",
        "d_SSIM_empirical_pvalue": "bg_d_SSIM_median",
        "d_Spearman_empirical_pvalue": "bg_d_Spearman_median",
        "tw_score_empirical_pvalue": "bg_tw_score_median",
        "d_NCC_empirical_pvalue": "bg_d_NCC_median"
    }

    def __init__(self, df, comparison, outdir):
        self.df = df
        self.comparison = comparison
        self.outdir = os.path.join(outdir, comparison)
        self.saved_paths = []
        os.makedirs(self.outdir, exist_ok=True)
        plt.style.use("seaborn-v0_8-whitegrid")
        plt.rcParams.update({"font.size": 12, "figure.dpi": 150})

    def _save(self, name, fig):
        png = os.path.join(self.outdir, f"{name}.png")
        svg = os.path.join(self.outdir, f"{name}.svg")
        fig.savefig(png, dpi=150, bbox_inches="tight")
        fig.savefig(svg, bbox_inches="tight")
        self.saved_paths.append((name, png, svg))
        return fig

    def _neg_log10(self, pvals, floor=1e-30):
        return -np.log10(np.maximum(pvals, floor))
```

**11 绘图方法 (all return (fig, ax_or_axes)):**

1. `hist_pvalues(self)` — 2×3 grid of histograms (5 scores + knn_dist), x-axis in log scale, each subplot has median annotation

2. `manhattan_style(self, score_col)` — scatter: x=window_start (or window_id), y=-log10(pval), colored by chrom, one plot per score via subplot panels or separate calls

3. `scatter_matrix(self)` — 5×5 seaborn pairplot / matplotlib subplots showing pairwise scatter of the 5 p-values, sample to 5000 rows for speed

4. `heatmap_chrom_window(self, score_col, n_chroms=20)` — heatmap of -log10(pval) with chroms as rows, windows as columns (or pivot on chrom+window_id), use imshow with diverging colormap centered at 1 (p=0.1)

5. `knn_dist_vs_pvalue(self)` — 5 subplots: knn_dist_median vs each p-value, colored by knn_valid_k, add Spearman R annotation

6. `bg_deviation(self)` — for each score: (actual - bg_median) / bg_median, plot as scatter or hexbin, 5 subplots

7. `pca_plot(self, n_components=2)` — PCA on 5 p-values (StandardScaler → PCA), scatter colored by -log10(MSE_pvalue) quantile, add explained variance bar

8. `chrom_boxplot(self)` — boxplot of -log10(pval) per chrom, one subplot per score, rotate x-labels

9. `qq_plot(self)` — Q-Q plot per score (observed p-values vs uniform [0,1]), add diagonal reference line

10. `correlation_heatmap(self)` — Pearson correlation matrix of the 5 p-values, annotate with r values, triangular heatmap

11. `summary_stats(self)` — text/table figure showing: row count, chrom count, % significant (p<0.05, p<0.01, p<0.001) per metric, knn_dist range

**main() function:**
- argparse: `--comparison` (optional), `--data-dir` (default "results/window_compare"), `--outdir` (default "results/figures")
- If `--comparison` not provided, auto-discover and process all comparisons
- For each comparison: load data → create plotter → call all 11 methods → print saved paths

**Style rules:**
- rcParams: `font.size=12`, `figure.dpi=150`, `axes.labelsize=11`, `axes.titlesize=13`
- Color palette: seaborn default or "muted"
- Chrom colors: use seaborn color palette cycling through at least 15 distinct colors
- All figures save both PNG (dpi=150, bbox_inches="tight") and SVG
- All methods return (fig, ax_or_axes) for notebook inline display
  </action>
  <verify>
    <automated>cd /share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff && python jupyter/knn_viz.py --comparison H1398_vs_L134 2>&1 | tail -5 && ls results/figures/H1398_vs_L134/*.png 2>/dev/null | wc -l</automated>
  </verify>
  <done>
    jupyter/knn_viz.py exists, KNNVizPlotter class has all 11 methods,
    python jupyter/knn_viz.py --help prints usage,
    python jupyter/knn_viz.py --comparison H1398_vs_L134 generates 11+ PNG files in results/figures/H1398_vs_L134/
  </done>
</task>

<task type="auto">
  <name>Task 2: 创建 jupyter/knn_eda.ipynb — 主 Notebook</name>
  <files>jupyter/knn_eda.ipynb</files>
  <action>
Create `jupyter/knn_eda.ipynb` as a valid Jupyter Notebook (JSON format, no nbformat dependency needed — hand-craft the JSON).

**Notebook structure (16 cells):**

Cell 1 — Markdown: "# KNN Empirical P-value EDA\n\nAnalysis of KNN-based empirical P-values for window comparison results.\n\n**Comparison:** `{comparison}`\n**Data:** `results/window_compare/{comparison}/{comparison}_knn_enriched.tsv`"

Cell 2 — Code: Setup imports
```python
import sys, os
# Project root is cwd — jupyter/ is a subdir, add cwd to resolve "jupyter.knn_viz"
PROJECT_ROOT = os.getcwd()
sys.path.insert(0, PROJECT_ROOT)
os.makedirs("results/figures", exist_ok=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from jupyter.knn_viz import KNNVizPlotter, auto_discover_comparisons, load_data
plt.style.use("seaborn-v0_8-whitegrid")
print("Available comparisons:", auto_discover_comparisons("results/window_compare"))
```

Cell 3 — Code: Select comparison (edit this cell to switch)
```python
COMPARISON = "H1398_vs_L134"  # <-- change this to switch comparison
df = load_data(COMPARISON)
print(f"Loaded {len(df):,} windows for {COMPARISON}")
print("Chromosomes:", sorted(df["chrom"].unique()))
print(df[KNNVizPlotter.SCORE_COLS].describe())
```

Cell 4 — Code: Initialize plotter
```python
plotter = KNNVizPlotter(df, COMPARISON, "results/figures")
```

Cell 5 — Markdown: "## 1. P-value Distribution Histograms"
Cell 6 — Code: `fig, ax = plotter.hist_pvalues(); plt.show()`

Cell 7 — Markdown: "## 2. Manhattan-Style Plots"
Cell 8 — Code: loop over SCORE_COLS and call `plotter.manhattan_style(col)` per score

Cell 9 — Markdown: "## 3. Scatter Matrix (5 Metrics)"
Cell 10 — Code: `fig, axes = plotter.scatter_matrix(); plt.show()` (may need `plt.tight_layout()`)

Cell 11 — Markdown: "## 4. KNN Distance vs P-values"
Cell 12 — Code: `fig, axes = plotter.knn_dist_vs_pvalue(); plt.show()`

Cell 13 — Markdown: "## 5. Background Deviation Analysis"
Cell 14 — Code: `fig, axes = plotter.bg_deviation(); plt.show()`

Cell 15 — Markdown: "## 6. PCA Projection"
Cell 16 — Code: `fig, ax = plotter.pca_plot(); plt.show()`

Cell 17 — Markdown: "## 7. Chromosome Boxplots"
Cell 18 — Code: `fig, axes = plotter.chrom_boxplot(); plt.show()`

Cell 19 — Markdown: "## 8. Q-Q Plots"
Cell 20 — Code: `fig, axes = plotter.qq_plot(); plt.show()`

Cell 21 — Markdown: "## 9. P-value Correlation Heatmap"
Cell 22 — Code: `fig, ax = plotter.correlation_heatmap(); plt.show()`

Cell 23 — Markdown: "## 10. Summary Statistics"
Cell 24 — Code: `fig, ax = plotter.summary_stats(); plt.show()`

Cell 25 — Code: Print saved paths
```python
print(f"\nSaved figures for {COMPARISON}:")
for name, png, svg in plotter.saved_paths:
    print(f"  {name}: {png}")
```

**JSON notebook format requirements:**
- `"nbformat": 4`, `"nbformat_minor": 5`
- Each cell has `cell_type` ("markdown" or "code"), `metadata`, `source` (list of strings)
- Code cells have `execution_count: null` (or integer), `outputs: []`
- File saved as UTF-8 JSON
  </action>
  <verify>
    <automated>cd /share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff && python -c "import json; json.load(open('jupyter/knn_eda.ipynb')); print('Valid JSON')" && grep -c '"cell_type"' jupyter/knn_eda.ipynb</automated>
  </verify>
  <done>
    jupyter/knn_eda.ipynb is valid JSON (parses without error),
    contains 15+ cells (mix of markdown and code),
    has COMPARISON variable, plotter instantiation, and calls to all figure methods
  </done>
</task>

</tasks>

<verification>
- `python jupyter/knn_viz.py --help` exits 0 and shows usage
- `python jupyter/knn_viz.py --comparison H1398_vs_L134` generates 11 PNG + 11 SVG in `results/figures/H1398_vs_L134/`
- `python jupyter/knn_viz.py` (auto-discover) processes H1398_vs_L134
- `jupyter/knn_eda.ipynb` is valid JSON notebook
- `jupyter nbconvert --to notebook --execute jupyter/knn_eda.ipynb --output jupyter/knn_eda_executed.ipynb` exits 0 (or at minimum, notebook runs without import errors)
</verification>

<success_criteria>
- jupyter/knn_viz.py: KNNVizPlotter class with 11 plotting methods + auto-discover + CLI
- jupyter/knn_eda.ipynb: 25+ cells covering setup through all 10+ figure types
- results/figures/H1398_vs_L134/: 11+ PNG files (dpi=150) and 11+ SVG files
- All 5 score metrics visualized across multiple figure types
- saved_paths records all generated figures with PNG and SVG paths
</success_criteria>

<output>
After completion, create `.planning/phases/02-knn-viz/02-SUMMARY.md`
</output>
