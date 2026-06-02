# Phase 02: KNN 结果可视化 - Context

**Gathered:** 2026-05-30
**Status:** Ready for planning

<domain>
## Phase Boundary

为 KNN 经验 P 值结果构建 Jupyter Notebook，实现全面的交互式探索性数据分析（EDA）。每个 comparison（目前有 H1398_vs_L134）生成一个独立可视化 notebook。核心数据为 `results/window_compare/{comparison}/{comparison}_knn_enriched.tsv`。

</domain>

<decisions>
## Implementation Decisions

### 可视化类型组合
- **全面探索型（10+ 张图）**，覆盖：
  - P 值分布直方图（每个指标一张，log 刻度）
  - Manhattan-style 基因组位置图（按基因组位置 coloring by -log10(p)）
  - 指标两两散点图矩阵（5x5 pairwise scatter matrix）
  - 热力图（染色体 × 窗口，P 值或 -log10(p) 值）
  - knn_dist_median vs P-value 散点（KNN 质量诊断）
  - bg_median vs 实际值偏离度散点
  - 5 指标 PCA / UMAP 可视化
  - 额外探索图（待 planner 细化）

### 交互性深度
- **纯静态图**，matplotlib/seaborn
- 不使用 Plotly / ipywidgets，保持简单可靠

### Notebook 结构
- **工具类封装**，Plotter 类封装所有绘图方法
- Notebook 只调用 `plotter.hist_pvalues()`, `plotter.manhattan()`, etc.
- 数据加载逻辑也封装到 Plotter 或独立的 loader 模块

### 多 comparison 支持
- **自动发现**，使用 `os.listdir()` 或 `glob` 扫描 `results/window_compare/`
- Notebook 自动检测所有 available comparisons
- 生成 comparison 汇总对比表（跨 comparison 的指标分布）

### 输出形式
- **inline + SVG + PNG + 路径记录**
- 保存到 `results/figures/{comparison}/`
- 同时输出矢量图（出版质量）和位图（dpi=150）
- 记录所有输出路径到 DataFrame 便于后续引用

### Claude's Discretion
- 具体的 matplotlib 配色方案和风格（rcParams 配置）
- 各图的子图布局（1 列 vs 2 列 vs grid）
- 字体大小、标题、轴标签的规范
- 异常值处理策略（是否剔除极端 P 值用于可视化）
- SVG vs PNG 的具体 DPI 设置

</decisions>

<specifics>
## Specific Ideas

- 数据文件路径：`results/window_compare/{comparison}/{comparison}_knn_enriched.tsv`
- 57,872 窗口 × 18 列（当前 H1398_vs_L134）
- 5 个 P 值指标：MSE, d_SSIM, d_Spearman, tw_score, d_NCC
- KNN 质量指标：knn_valid_k, knn_dist_median, knn_dist_max
- 背景基线：bg_{score}_median
- 窗口 ID 格式：`chrom`（如 HC04_A01）+ window_id（整数）
- P 值范围：~1e-5 到接近 1，log 刻度必要

</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase 01 Context
- `.planning/phases/01-knn-empirical-pvalue/01-CONTEXT.md` — 数据 schema、5 指标含义、背景基线逻辑

### Notebook 规范
- `docs/snakemake/wildcards.md` — 无直接关联，但项目规范参考

### 无外部依赖
- 不需要额外库，matplotlib/seaborn/scikit-learn 为 pipeline 环境已有

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- Phase 01 产生的 `knn_enriched.tsv` 即为 notebook 的输入数据
- Pipeline 的 conda 环境（`env/SNAP-HiC-Diff.yaml`）已包含 matplotlib/seaborn/numpy/pandas/scikit-learn

### Established Patterns
- Snakemake 的 `results/` 目录结构规范：`results/window_compare/{comparison}/`
- Output 保存路径参考：`results/figures/{comparison}/`

### Integration Points
- Notebook 读取已生成的 `.tsv` 文件，不修改 pipeline
- 可作为 pipeline 最终 step（Snakemake rule 触发 notebook 执行，或独立运行）

</code_context>

<deferred>
## Deferred Ideas

- Snakemake rule 调用 notebook 执行（nbconvert）— Phase 02 只做 notebook 本身
- 自动化报告生成（PDF/HTML export from notebook）— future phase
- 与 UCSC Genome Browser 或 HiGlass 的集成可视化 — future phase

</deferred>

---

*Phase: 02-knn-viz*
*Context gathered: 2026-05-30*
