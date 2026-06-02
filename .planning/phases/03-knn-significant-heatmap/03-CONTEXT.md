# Phase 03: KNN 显著窗口热图绘制 - Context

**Gathered:** 2026-05-30
**Status:** Ready for planning

<domain>
## Phase Boundary

新增一个 Snakemake rule，筛选 `knn_enriched.tsv` 中 5 个指标 FDR < 0.05 的显著窗口，针对每个窗口绘制一张 2x2 热图：上面两张为 control/treatment 的 iced（O/E matrix），下面两张为 control/treatment 的 gaussian smoothed heatmap。输出 PNG + PDF。

</domain>

<decisions>
## Implementation Decisions

### 数据源
- **iced 热图**: `results/balanced/{sample}/{sample}__{chrom}_iced.cool`
  - **per-chrom 文件，非 per-window**！用 `cooler.Cooler().matrix().fetch(chrom, window_start, window_end)` 提取窗口区域
  - 扫描 `results/balanced/` 目录匹配 `{sample}__{chrom}_iced.cool`
- **gauss 热图**: `results/windows/{sample}/{sample}__{chrom}_{window_id}_gaussian.h5`
  - 每个窗口一个 h5 文件，存储该窗口的高斯平滑矩阵
- **显著窗口筛选**: 来自 `results/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv`
  - **注意**：`knn_enriched.tsv` 无 FDR 列（`compute_fdr: false`），checkpoint 脚本内 inline BH 计算
  - 筛选条件：5 个 `*_empirical_pvalue` 列经 BH 校正后 FDR 全部 < 0.05

### 热图坐标轴
- X/Y 轴均为基因组坐标轴（bp），非窗口 index
- 每个轴最多 5 个刻度

### 布局
- 2 行 2 列：
  - 上: control iced（左）/ treatment iced（右）
  - 下: control gauss（左）/ treatment gauss（右）
- 颜色条：
  - iced（control + treatment）共享同一颜色条
  - gauss（control + treatment）共享同一颜色条

### 热图区域范围
- 固定窗口大小，即窗口自身的 `window_start` ~ `window_end` 范围

### 运行方式
- Snakemake rule 驱动
- 输入：knn_enriched.tsv（筛选） + cool/h5 文件
- 输出：`results/knn_significant_heatmap/{comparison}/{chrom}_{window_id}.png/pdf`

### 输出格式
- PNG + PDF 双格式
- 布局：2x2 subplot

### Snakemake 集成
- rule 名：`plot_knn_significant_heatmap`
- wildcards: `{comparison}`, `{chrom}`, `{window_id}`
- checkpoint 架构：checkpoint 扫描 knn_enriched.tsv 并 inline BH 计算 FDR，输出 manifest；regular rule 每个窗口一个 job
- checkpoint 脚本内 inline BH 计算（`scipy.stats.false_discovery_control`），不改已有 pipeline

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase 01 / 02 Context
- `.planning/phases/01-knn-empirical-pvalue/01-CONTEXT.md` — 数据 schema、5 指标含义、FDR 列名（`*_empirical_FDR`）
- `.planning/phases/02-knn-viz/02-CONTEXT.md` — KNNVizPlotter 封装模式、plt.rcParams 规范

### Snakemake 规范
- `rules/window_compare.smk` — 现有 rule 写法参考（conda、params、log、benchmark 字段）
- `Snakefile` — include 规范、OUTDIR 变量、rule all target 格式
- `config/config.yaml` — knn_empirical_pvalue section

### 数据文件路径规范
- Cool 文件: `results/balanced/{sample}/{sample}__{chrom}_{window_id}_iced.cool`
- Gaussian h5: `results/windows/{sample}/{sample}__{chrom}_{window_id}_gaussian.h5`
- Knenriched: `results/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv`

### Pipeline 环境
- `env/SNAP-HiC-Diff.yaml` — 必须包含 cooltools（用于读写 .cool 文件）和 h5py（用于读写 .h5 文件）

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `jupyter/knn_viz.py` 的 `KNNVizPlotter` 封装模式：`_save(name, fig)` 同时保存 PNG + SVG
- `plt.rcParams` 配置：font.size=12, figure.dpi=150
- Snakemake rule 模板：参考 `rules/knn_empirical_pvalue.smk`

### Established Patterns
- `results/` 目录结构：`results/{module}/{comparison}/`
- Snakemake wildcards: `{comparison}`, `{chrom}`, `{window_id}`
- 输出文件格式：`{comparison}_{chrom}_{window_id}_heatmap.png/pdf`

### Integration Points
- 从 `knn_enriched.tsv` 读取窗口列表和筛选条件
- 从 `results/balanced/` 扫描 cool 文件
- 从 `results/windows/` 读取 gaussian h5 文件

</code_context>

<specifics>
## Specific Ideas

- 窗口路径规律: `{sample}__{chrom}_{window_id}_iced.cool`（如 `L128__HC04_D04_iced.cool`）
- Gaussian h5 文件每个窗口一个
- FDR 列名: `MSE_empirical_FDR`, `d_SSIM_empirical_FDR`, `d_Spearman_empirical_FDR`, `tw_score_empirical_FDR`, `d_NCC_empirical_FDR`
- 5 个 FDR 全部 < 0.05 才算显著窗口

</specifics>

<deferred>
## Deferred Ideas

- 每个显著窗口单独一个图 → 显著窗口数量可能很多（如 500+），未来考虑聚合视图（一张图多个窗口）
- HiGlass 集成可视化 — future phase

---

*Phase: 03-knn-significant-heatmap*
*Context gathered: 2026-05-30*
