# Phase 04: Meltron 融解/凝结分析 - Context

**Gathered:** 2026-06-03
**Status:** Ready for planning

<domain>
## Phase Boundary

在 KNN 经验 P 值筛选出的显著差异窗口（FDR < 0.05）上，通过 Meltron 框架计算 chromatin 融解/凝结 score。Meltron 通过 KS 检验比较两组 insulation score（IS）分布，判断基因组区间在条件间是否发生了 chromatin 拓扑变化。

**Phase 04 交付物：**
- `rules/compute_IS.smk` — IS 计算 Snakemake rule（cooltools CLI + Python 格式转换）
- `rules/meltron.smk` — Meltron 分析 Snakemake rule（R subprocess 调用）
- `scripts/compute_IS_batch.py` — 多档 IS 计算 + 格式转换脚本
- `scripts/run_meltron.py` — Meltron CLI 包装脚本
- `config/config.yaml` — IS 距离和 threshold 配置
- `Snakefile` — 包含新 rule 和 all target

**不在本 phase：**
- Meltron 结果可视化（Phase 05）

</domain>

<decisions>
## Implementation Decisions

### ROI 来源
- **D-01:** ROI = KNN 显著窗口。筛选条件：`{comparison}_knn_enriched.tsv` 中 5 个 empirical_pvalue 列经 BH 校正 FDR 全部 < 0.05 的窗口

### ROI ID 格式
- **D-02:** 双 ID 格式。Meltron genomic_intervals 文件包含四列 `(chrom, start, end, ID)`
- **D-03:** `ID = "{chrom}_{window_idx}"`（窗口在染色体上的序号）
- **D-04:** 基因组坐标作为额外列保留：`{chrom}:{start}-{end}`，便于结果回溯

### IS 计算架构
- **D-05:** 新建独立 Snakemake rule（`compute_IS`），不在现有 extract_expected 或其他 rule 内扩展
- **D-06:** IS 结果按 sample 输出，所有 comparison 共用（Snakemake DAG 自然去重，不重复计算）

### IS 计算工具
- **D-07:** 使用 `cooltools insulation` CLI，Python 脚本包装（**实测修正**：命令是 `insulation` 不是 `cis-insulation`）
- **D-08:** WINDOW 是**位置参数**，一次可传多个：`cooltools insulation sample.cool 100000 200000 500000 1000000 -o out.tsv -p 8`
- **D-08b:** cooltools 输出多列：`chrom, start, end, log2_insulation_score_{WINDOW}, n_valid_pixels_{WINDOW}, boundary_strength_{WINDOW}, is_boundary_{WINDOW}, ...`，需 Python 转换脚本转成 Meltron 需要的 4 列 `(chrom, viewpoint, IS_distance, IS)`
- **D-08c:** 25 个染色体 × 1 次命令 = 25 次 cooltools 调用（不是 per-distance 循环的 250 次）

### IS 距离尺度
- **D-09:** 多档 IS distances，通过 `config["meltron"]["is_distances"]` 配置
- **D-10:** 默认值：100kb 到 1Mb，10 档（步长 100kb），与 Meltron 原论文一致
- **D-11:** IS square 大小 = IS distance（即用 IS 距离本身作为 square 大小，Meltron 原版逻辑）

### IS 计算粒度
- **D-12:** Per-chromosome IS 计算。每个染色体内独立跑 cooltools cis-insulation
- **D-13:** Meltron 对每个窗口做 KS test 时，只用该窗口所在染色体的 IS 数据（合理，TAD insulation 主要是染色体内现象）

### Comparison pairs
- **D-14:** 与 KNN 保持一致，从 `window_compare_pairs.yaml` 读取
- **D-15:** 每个 comparison（A_vs_B）：A 的 IS 作为 reference，B 的 IS 作为 input，跑 Meltron

### Meltron CLI 包装
- **D-16:** R subprocess 包装。Python 脚本调用 `Rscript command_line_apps/Meltron.R`，通过 Snakemake conda 环境管理 R 依赖（tidyverse, DescTools, data.table）
- **D-17:** Meltron 参数：`--referenceIS`（A 的 IS）、`--inputIS`（B 的 IS）、`--genomicIntervals`（KNN 显著窗口）、`--outfile`、`--comparison`（decondensation 或 condensation）、`--cutoff_NA`（默认 0.5）
- **D-18:** `--comparison` 方向：A_vs_B 表示 A 相对 B 融解（decondensation）还是凝结（condensation），通过 config 配置

### Score threshold
- **D-19:** melting_score threshold 通过 `config["meltron"]["score_threshold"]` 配置，不写死。**默认 5.0 是 advisory，不是 genome-wide significance 阈值**
- **D-20:** 输出列包含原始 p-value、Bonferroni 校正后 p-value、melting_score（-log10(corrected_pval)）。未来应用 matched null windows 做 empirical calibration

### 输出文件
- **D-21:** `meltron_merged.tsv` — 最终输出（decondensation + condensation 合并）
- **D-22:** 输出路径：`results/meltron/{comparison}/{comparison}_meltron_merged.tsv`
- **D-22b:** 中间文件：`_intervals.tsv`（4 列）、`_interval_map.tsv`（窗口映射）、`_decondensation.tsv`、`_condensation.tsv`
- **D-23:** 输出列：`chrom, start, end, ID, meltron_melting_score, meltron_compaction_score, meltron_best_score, meltron_best_direction, meltron_is_significant, window_id, n_significant_metrics`

### Snakemake Rule DAG
```
# Phase 04 DAG
compute_IS (checkpoint, per sample)
  └── compute_IS_h5 (rule, declares per-chrom IS .npz/.tsv files)
        ↓
  └── meltron (rule, per comparison)
        input: knn_enriched.tsv (significant windows) + IS tables for both samples
        output: {comparison}_meltron.tsv
```

### Meltron IS 格式要求
- **D-24:** IS table 格式：`(chrom, viewpoint, IS_distance, IS)`，4 列 tab 分隔（Meltron 源码验证）
- **D-25:** IS distance 以 bp 为单位（如 100000 表示 100kb）
- **D-26:** 每个 sample 输出一个 IS table 文件（所有 chromosome 合并）

### Meltron 运行模式（实测修正）
- **D-27:** `cooltools insulation` 是正确命令（不是 `cis-insulation`）
- **D-28:** WINDOW 是位置参数，一次可传多个：`cooltools insulation sample.cool 100000 200000 500000 1000000 -o out.tsv -p 8`
- **D-29:** 每个样本 × 每个染色体 = 1 次 cooltools 调用（25 次，不是 per-distance 循环的 250 次）
- **D-30:** cooltools 输出格式：`chrom, start, end, log2_IS_{WINDOW}, n_valid_pixels_{WINDOW}, boundary_strength_{WINDOW}, is_boundary_{WINDOW}`（每档一列组）

### Meltron interval 格式（严格 4 列）
- **D-31:** Meltron genomic intervals 必须严格 4 列：`(chrom, start, end, ID)`，无额外列
- **D-32:** 窗口映射通过独立文件 `*_interval_map.tsv` 保存（包含 window_id, n_significant_metrics）

### 运行方向
- **D-33:** 每个 comparison 同时跑 decondensation 和 condensation 两个方向，不只跑一个
- **D-34:** decondensation → 融解事件（sample B 相对 A 失去 TAD boundary），condensation → 凝结事件

### 候选窗口策略（post-selection bias 处理）
- **D-35:** min_significant=1（任一指标 FDR < 0.05 即入 Meltron），Meltron score 作为 support evidence，不做 genome-wide significance 宣称
- **D-36:** score_threshold=5.0 advisory，后续应通过 matched null windows 做 empirical FDR 校准

</decisions>

<specifics>
## Specific Ideas

- 用户要求先完成 Phase 03（KNN significant heatmap），再开始 Phase 04
- 用户提到"meltron 可以深入分析"——意味着后续可能有 Phase 05 的可视化
- Meltron 原文 DOI: 10.1038/s41586-021-04081-2
- 参考：pombo-lab/Meltron (GitHub)
- KNN 结果表（knn_enriched.tsv）已有窗口位置信息，复用为 Meltron 的 genomic_intervals

**Phase 分解（根据漏洞分析修正）：**
```
Phase 04: Meltron insulation-level validation（当前 phase）
Phase 05: boundary calling + boundary matching（未来）
Phase 06: TAD rearrangement type classification（gain/loss/shift/split/merge）（未来）
```

**当前 Phase 04 scope（修正后）：**
- 不是独立的 TAD 重排 caller
- 是 KNN 候选窗口的 insulation-level 二次验证模块
- 输出 meltron support score，标注每个候选窗口
- 后续 phase 用 boundary calling + DI + classifier 完成 TAD 重排分类

</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

- `.planning/ROADMAP.md` — Phase 03 未完成，Phase 04 依赖 Phase 01 的 knn_enriched 数据
- `.planning/phases/01-knn-empirical-pvalue/01-CONTEXT.md` — KNN pipeline 上下文
- `.planning/phases/03-knn-significant-heatmap/03-CONTEXT.md` — significant window 筛选逻辑
- `config/window_compare_pairs.yaml` — comparison pairs 定义
- `rules/compute_gap_qc.smk` — 现有 checkpoint 模式参考
- `rules/compute_oe_matrix.smk` — 现有 per-sample rule 模式参考
- `https://github.com/pombo-lab/Meltron` — Meltron 源码和 CLI 文档

### Meltron 源码
- `command_line_apps/Meltron.R` — CLI 入口，参数解析，KS test + Bonferroni 校正
- `command_line_apps/long_matrix_to_IS.R` — 多档 IS 计算参考（但用 cooltools 替代）
- `command_line_apps/matrix_wide_to_long.R` — 格式转换参考

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `compute_gap_qc.py` — checkpoint shell block 循环处理多染色体的模式，可参考
- `compute_oe_matrix_batch.py` — 多染色体批量处理 Python 脚本，可参考其参数传递模式
- `convert_cool_to_mcool.py` — mcool 输入 + shell block 循环模式，可参考

### Established Patterns
- **checkpoint 架构**：checkpoint 声明真实 output 文件，DAG 依赖 checkpoint 而非 rule，DAG 自然等待
- **conda 环境**：所有 rule 用 `../env/SNAP-HiC-Diff.yaml`，R 依赖通过 conda R 包管理
- **shell block 规范**：`set -euo pipefail`，变量用 `${{VAR}}` 转义，checkpoint manifest 作为完成信号

### Integration Points
- 上游输入：`results/knn_significant/{comparison}/{comparison}_knn_significant.tsv`（FDR < 0.05 窗口）
- IS 输入：`results/IS/{sample}/{sample}_IS.tsv`（新 rule 输出）
- 输出：`results/meltron/{comparison}/{comparison}_meltron.tsv`

</code_context>

<deferred>
## Deferred Ideas

- Phase 05: boundary calling + boundary matching — 基于 Meltron intervals 和 cooltools insulation boundary calling
- Phase 06: TAD rearrangement type classification (gain/loss/shift/split/merge/complex) — 基于 boundary matching + DI
- Per-chromosome IS table 的 chunked 输出（如果内存压力大）— 后续优化

</deferred>

---

*Phase: 04-meltron*
*Context gathered: 2026-06-03*
