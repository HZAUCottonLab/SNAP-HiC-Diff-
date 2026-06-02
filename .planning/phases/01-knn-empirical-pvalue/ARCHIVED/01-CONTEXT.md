# Phase 01: KNN Empirical P-Value - Context

**Gathered:** 2026-05-30
**Status:** Ready for planning (E deferred)

<domain>
## Phase Boundary

新增一个 Snakemake rule `knn_empirical_pvalue`，对每个比较窗口（comparison）根据其 `EC`、`IB`、`SP` 质量特征，从背景窗口池中找到数据质量最相近的 10000 个背景窗口，用这 10000 个窗口的统计量构建局部背景分布，计算该比较窗口的经验 P 值。

`EC`、`IB`、`SP` 仅用于 KNN 质量匹配，不参与显著性检验。
P 值计算方向统一：值越大表示差异越强，右尾检验。

**Scope 外：**
- 计算 MSE/SSIM/Spearman/NCC/TW p-value（已有）
- 计算 EC/IB/SP（已有）
- FDR 校正（用户明确不需要）
- E 待讨论（match_quality_flag 阈值）
</domain>

<decisions>
## Implementation Decisions

### A. 输入绑定
- **D-01:** 输入固定为 `{OUTDIR}/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv`，即 `rules/window_compare.smk` 中 `rule merge_window_pair_metrics` 的输出文件
- **D-02:** 该文件已包含所有必需列：`chrom, start, end, sample_a, sample_b, EC, IB, SP, MSE, SSIM, Spearman, tw_pvalue, NCC`；列名在 merge 脚本中确认（`NCC` 列存在于 merge 输出）

### B. 规则粒度
- **D-03:** 每个 comparison 跑一个独立的 Snakemake job
- **D-04:** Wildcard 为 `{comparison}`，与现有 `merge_window_pair_metrics` 一致
- **D-05:** 36 个 comparison → 36 个 job，可并行调度

### C. Config 暴露程度
- **D-06:** 以下参数暴露给 `config/config.yaml` 的 `knn_empirical_pvalue` 小节：

```yaml
knn_empirical_pvalue:
  background_file: "/share/home/xhhuang24/project/PanTE/background_pool/results/window_compare/HiC_window_background.tsv"
  knn_k: 10000
  knn_metric: "euclidean"
  scaler: "robust"
  score_columns:
    - MSE
    - d_SSIM
    - d_Spearman
    - tw_score
    - d_NCC
  derive_scores:
    d_SSIM: "1 - SSIM"
    d_Spearman: "1 - Spearman"
    tw_score: "-log10(tw_pvalue)"
    d_NCC: "1 - NCC"
  pvalue_direction:
    MSE: "greater"
    d_SSIM: "greater"
    d_Spearman: "greater"
    tw_score: "greater"
    d_NCC: "greater"
  compute_fdr: false
```

- **D-07:** `knn_features` 固定为 `["EC", "IB", "SP"]`，不暴露（用户未要求可配置）
- **D-08:** `tw_pvalue` 安全下限固定为 `1e-30`，不暴露

### D. 输出结构
- **D-09:** 输出路径：`{OUTDIR}/window_compare/{comparison}/{comparison}_knn_enriched.tsv`
- **D-10:** 新建 rule `knn_empirical_pvalue`，不覆盖原始 `{comparison}_window_pair_metrics.tsv`
- **D-11:** 输出保留原始列，新增以下字段：

| 字段 | 说明 |
|------|------|
| `knn_valid_k` | 实际用于计算 P 值的背景窗口数量（≤ knn_k） |
| `knn_dist_median` | 候选窗口到 10000 个匹配背景窗口的距离中位数 |
| `knn_dist_max` | 第 10000 个背景窗口的距离（最远匹配距离） |
| `bg_MSE_median` | 匹配背景窗口的 MSE 中位数 |
| `bg_d_SSIM_median` | 匹配背景窗口的 d_SSIM 中位数 |
| `bg_d_Spearman_median` | 匹配背景窗口的 d_Spearman 中位数 |
| `bg_tw_score_median` | 匹配背景窗口的 tw_score 中位数 |
| `bg_d_NCC_median` | 匹配背景窗口的 d_NCC 中位数 |
| `MSE_empirical_pvalue` | 右尾经验 P 值 |
| `d_SSIM_empirical_pvalue` | 右尾经验 P 值 |
| `d_Spearman_empirical_pvalue` | 右尾经验 P 值 |
| `tw_score_empirical_pvalue` | 右尾经验 P 值 |
| `d_NCC_empirical_pvalue` | 右尾经验 P 值 |
| `match_quality_flag` | 匹配质量标记（见 D-13） |

- **D-12:** 日志文件：`{OUTDIR}/logs/window_compare/{comparison}/knn_empirical_pvalue.log`，记录：背景窗口有效数量、EC/IB/SP 的 median 和 IQR、BallTree 参数、KNN 的 K 值、各统计量 P 值计算方向、输出文件路径
- **D-13:** `match_quality_flag` 标记规则：**待讨论（E）**——暂定至少支持 `PASS`、`LOW_BG_COUNT`、`HIGH_KNN_DISTANCE`

### 算法核心
- **D-14:** 背景窗口数量约 2,152,388，使用 `sklearn.neighbors.BallTree`（内存高效）
- **D-15:** 标准化参数（median / IQR）**仅从背景窗口拟合**，不混入比较窗口
- **D-16:** 标准化后背景窗口和比较窗口使用同一套参数：`x_scaled = (x - median_bg) / IQR_bg`
- **D-17:** BallTree 用标准化后的 `EC_scaled, IB_scaled, SP_scaled` 构建
- **D-18:** 对每个比较窗口查询 K=10000 个最近邻
- **D-19:** 右尾经验 P 值公式：`p = (count(background_value >= query_value) + 1) / (K + 1)`
- **D-20:** 去除 `EC, IB, SP, MSE, SSIM, Spearman, tw_pvalue, NCC` 中包含 `NA`、`inf`、`-inf` 的记录

### Snakemake Rule 接入点
- **D-21:** 新 rule 放在 `rules/window_compare.smk` 底部，与 `merge_window_pair_metrics` 相邻
- **D-22:** 新 rule 依赖 `merge_window_pair_metrics` 的输出（`{comparison}_window_pair_metrics.tsv`）
- **D-23:** 新增 `rule all` 的 input 声明，或在 `merge_all_window_pair_metrics` 后串联

</decisions>

<specifics>
## Specific Ideas

- 用户明确采用 `d_NCC = 1 - NCC` 统一所有差异指标方向（右尾检验）
- 用户明确不需要 FDR 校正
- `tw_pvalue_safe = max(tw_pvalue, 1e-30)` 防止 `-log10(0)`
- 背景文件路径：`/share/home/xhhuang24/project/PanTE/background_pool/results/window_compare/HiC_window_background.tsv`
- 不对 `EC`、`IB`、`SP` 做 log transform（理由：EC 分布集中；IB 已是 abs(log2 ratio)；SP 是 0-1 比例）
</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### 本项目代码
- `rules/window_compare.smk` — 新 rule 的挂载位置和 wildcard 模式
- `scripts/merge_window_pair_metrics.py` — 输入文件列名参考（merge_keys = chrom/window_id/window_start/window_end；NCC 列存在于输出）
- `scripts/compute_window_pair_covariates.py` — EC/IB/SP 计算公式参考
- `config/config.yaml` — 新增 `knn_empirical_pvalue` 配置小节的位置
- `Snakefile` — `rule all` 接入点确认

### 项目规范
- `.cursor/rules/011-snakemake-wildcards.mdc` — wildcard 提取规范（shell block 开头 `SAMPLE="{wildcards.sample}"`）
- `.cursor/rules/010-snakemake-rule-format.mdc` — 标准 rule 字段顺序、benchmark/log 规范
- `docs/snakemake/wildcards.md` — `${{VAR}}` 转义规范

### 背景数据
- `/share/home/xhhuang24/project/PanTE/background_pool/results/window_compare/HiC_window_background.tsv` — 背景文件格式（含列名：chrom, window_id, window_start, window_end, EC, IB, SP, raw_sum_cpm_a, raw_sum_cpm_b, nz_frac_a, nz_frac_b, MSE, SSIM, Spearman, tw_statistic, tw_pvalue, comparison, sample_a, sample_b, NCC）

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `scripts/merge_window_pair_metrics.py` 的 merge_keys 模式（`chrom, window_id, window_start, window_end`）可直接复用
- `rules/window_compare.smk` 的 rule 结构（字段顺序、shell 块模式、conda 环境）可直接参照

### Established Patterns
- 每个 `{comparison}` wildcard job 的目录结构：`{OUTDIR}/window_compare/{comparison}/`
- benchmark/log 输出路径约定：`benchmark/window_compare/{comparison}/` 和 `logs/window_compare/{comparison}/`
- conda 环境固定为 `"../env/SNAP-HiC-Diff.yaml"`（BallTree 依赖 sklearn）

### Integration Points
- 输入：`rules.merge_window_pair_metrics.output.merged`
- 输出：新文件 `{comparison}_knn_enriched.tsv`
- 配置：`config["knn_empirical_pvalue"]`
- 线程：`config["parallel"]["threads_per_job"]`（KNN 脚本内可能需要 joblib 并行）
</code_context>

<deferred>
## Deferred Ideas

### E. match_quality_flag 阈值（待讨论）
- `HIGH_KNN_DISTANCE` 的判断标准未确定：
  - 基于 KNN 距离分布的百分位（如 > 95th percentile）
  - 基于绝对阈值（需确定数值）
  - 基于标准化空间的最大距离（如 > 3 IQR）
- `LOW_BG_COUNT` 触发条件未确定（当有效背景窗口 < knn_k 时）
- Flag 的实际使用场景未明确（下游如何消费这个 flag？）
- **行动：** 下次讨论时先确认 downstream 用途，再定阈值策略

</deferred>

---

*Phase: 01-knn-empirical-pvalue*
*Context gathered: 2026-05-30*
