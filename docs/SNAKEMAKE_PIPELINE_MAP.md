# SNAKEMAKE_PIPELINE_MAP

## 1. Snakefile 组织结构

### include 顺序
```
Snakefile
  ├─ config/config.yaml         # 主配置
  ├─ config/chroms.yaml         # 染色体列表
  ├─ rules/validate_inputs.smk
  ├─ rules/build_windows.smk
  ├─ rules/extract_chrom.smk    # checkpoint extract_and_balance_sample
  ├─ rules/extract_expected.smk # extract_expected_tsv + extract_expected_npz
  ├─ rules/compute_oe_matrix.smk
  ├─ rules/compute_gap_qc.smk  # checkpoint compute_gap_qc
  ├─ rules/sliding_window.smk  # checkpoint sliding_window + rule sliding_window_h5
  └─ rules/window_compare.smk   # 5 rules + 1 rule merge_all
```

### Helper Functions（Snakefile 层）
- `get_chroms()` → 返回 CHROMS
- `get_sample_mcool(sample)` → 从 manifest 读取 cool_uri
- `get_sample_pairs()` → 返回 SAMPLES 两两组合
- `get_comparison_pairs()` → 读取 `window_compare_pairs.yaml`
- `get_window_compare_pairs()` → 返回 comparison 名称列表
- `get_pair_samples(comparison)` → 返回 (sample_a, sample_b)

### Global Variables
- `CHROMS`：从 `config/chroms.yaml` 加载，26 条染色体
  - **⚠️ 命名必须和 mcool 文件一致**（`hictk dump --table chroms` 的输出）。项目历史上用过两套：`Ghir_*`（旧）和 `HC04_*`（当前）。不一致会导致 `extract_and_balance_sample` 全量失败。详见 `docs/BUG_FIX_GUIDE.md`
- `SAMPLES`：从 `manifest.tsv` 动态加载，74 个样本
- `OUTDIR`：来自 config["outdir"]，默认 "results"
- `RESOLUTION`：来自 config["resolution"]，默认 4000
- `WINDOW_SIZE`：来自 config["window_size"]，默认 300000
- `STEP_SIZE`：来自 config["step_size"]，默认等于 window_size
- `TMP`：来自 config["tmpdir"]，默认 "tmp"

## 2. Rule 执行顺序 (DAG)

```
rule all (final target)
    │
    ├─► validate_inputs
    │
    ├─► build_windows
    │
    ├─► [extract_and_balance_sample] × 74 samples (checkpoint)
    │       │
    │       ├─► [extract_expected_tsv] × 74 samples
    │       │       │
    │       │       └─► [extract_expected_npz] × 74 samples
    │       │
    │       ├─► [compute_oe_matrix] × 74 samples
    │       │
    │       └─► [compute_gap_qc] × 74 samples (checkpoint)
    │               │
    │               └─► [sliding_window] × 74 samples (checkpoint)
    │                       │
    │                       └─► [sliding_window_h5] × 74 samples (rule)
    │
    └─► [merge_all_window_pair_metrics] (最终汇总)
            │
            └─► [merge_window_pair_metrics] × 36 comparisons
                    ├─► [compute_window_pair_covariates] × 36
                    ├─► [compare_gaussian_smoothed_windows] × 36
                    └─► [compare_iced_windows_tw] × 36
```

### 并行度
- `extract_and_balance_sample`：每个 sample 一个 job × 74 = 74 jobs（每个内部循环 26 条染色体）
- `extract_expected_tsv / npz`：每个 sample 一个 job × 74
- `compute_oe_matrix`：每个 sample 一个 job × 74
- `compute_gap_qc`：每个 sample 一个 job × 74（checkpoint）
- `sliding_window`：每个 sample 一个 job × 74（checkpoint）
- `window_compare`：每个 comparison 一个 job × 36（3 个并行的 sub-rule）

## 3. Rule 详细说明

### rule validate_inputs
- **文件**：`rules/validate_inputs.smk`
- **Input**：config["manifest"] + `scripts/validate_inputs.py`
- **Output**：`{OUTDIR}/validated/manifest.validated.tsv`
- **Threads**：`config["parallel"]["threads_per_job"]`
- **Conda**：`../env/SNAP-HiC-Diff.yaml`
- **功能**：`hictk validate` 检查所有 mcool 文件有效性

### rule build_windows
- **文件**：`rules/build_windows.smk`
- **Input**：`validate_inputs` output + `scripts/build_windows.py`
- **Output**：`{OUTDIR}/windows/windows_{RESOLUTION}.parquet`
- **功能**：生成基因组窗口坐标（parquet）

### checkpoint extract_and_balance_sample
- **文件**：`rules/extract_chrom.smk`
- **Wildcard**：`{sample}`（74 个）
- **Input**：从 manifest 动态读取 mcool URI
- **Output**：`{OUTDIR}/logs/extract_balance/{sample}.done`
- **Threads**：`config["hictk"]["threads"]`（10）
- **功能**：每个 sample 循环 26 条染色体，`hictk dump` + `hictk load` + `hictk balance ice`
- **⚠️ 依赖 `config/chroms.yaml` 的值**：染色体名必须和 mcool 文件一致。详见 `docs/BUG_FIX_GUIDE.md` B-01

### rule extract_expected_tsv
- **文件**：`rules/extract_expected.smk`
- **Wildcard**：`{sample}`（74 个）
- **Input**：依赖 `extract_and_balance_sample` checkpoint
- **Output**：`{OUTDIR}/logs/expected/{sample}_cooltools.done`
- **功能**：`cooltools expected-cis` 生成距离分层期望值

### rule extract_expected_npz
- **文件**：`rules/extract_expected.smk`
- **Wildcard**：`{sample}`（74 个）
- **Input**：依赖 `extract_expected_tsv` + `scripts/extract_expected_batch.py`
- **Output**：`{OUTDIR}/logs/expected/{sample}_to_npz.done`
- **功能**：TSV → NPZ 格式转换

### rule compute_oe_matrix
- **文件**：`rules/compute_oe_matrix.smk`
- **Wildcard**：`{sample}`（74 个）
- **Input**：依赖 `extract_expected_npz` + `scripts/compute_oe_matrix_batch.py`
- **Output**：`{OUTDIR}/logs/oe/{sample}.done` + NPZ 文件
- **功能**：O/E 归一化矩阵计算

### checkpoint compute_gap_qc
- **文件**：`rules/compute_gap_qc.smk`
- **Wildcard**：`{sample}`（74 个）
- **Input**：依赖 `extract_and_balance_sample` checkpoint + `scripts/compute_gap_qc.py`
- **Output**：`{OUTDIR}/logs/gap_qc/{sample}.done`
- **功能**：gap 检测 + 动态 P20 QC

### checkpoint sliding_window
- **文件**：`rules/sliding_window.smk`
- **Wildcard**：`{sample}`（74 个）
- **Input**：依赖 `compute_gap_qc` checkpoint（gaps.tsv）+ O/E NPZ 文件
- **Output**：`{OUTDIR}/logs/windows/{sample}.done`
- **功能**：gap fill → z-score → Gaussian 平滑

### rule sliding_window_h5
- **文件**：`rules/sliding_window.smk`
- **Wildcard**：`{sample}`（74 个）
- **Input**：依赖 `sliding_window` checkpoint
- **Output**：所有 H5 文件（`ln_oe`, `gap_filled`, `zscore`, `gaussian`）
- **功能**：Snakemake 依赖追踪（mtime 追踪）

### rule compute_window_pair_covariates
- **文件**：`rules/window_compare.smk`
- **Wildcard**：`{comparison}`（36 个）
- **Input**：`scripts/compute_window_pair_covariates.py` + window_qc.tsv（双方样本）
- **Output**：`{OUTDIR}/window_compare/{comparison}/{comparison}_covariates.tsv`
- **功能**：EC、IB、SP 协变量计算

### rule compare_gaussian_smoothed_windows
- **文件**：`rules/window_compare.smk`
- **Wildcard**：`{comparison}`（36 个）
- **Input**：`covariates` + `sliding_window` checkpoint（双方）
- **Output**：`{OUTDIR}/window_compare/{comparison}/{comparison}_gaussian_metrics.tsv`
- **Threads**：`config["window_compare"]["gaussian_metrics"]["n_jobs"]`（80）
- **功能**：MSE、SSIM、Spearman、NCC

### rule compare_iced_windows_tw
- **文件**：`rules/window_compare.smk`
- **Wildcard**：`{comparison}`（36 个）
- **Input**：`covariates` + `extract_and_balance_sample` checkpoint（双方）
- **Output**：`{OUTDIR}/window_compare/{comparison}/{comparison}_tw.tsv`
- **Threads**：`config["window_compare"]["tw"]["n_jobs"]`（80）
- **功能**：Tracy-Widom 统计量

### rule merge_window_pair_metrics
- **文件**：`rules/window_compare.smk`
- **Wildcard**：`{comparison}`（36 个）
- **Input**：`covariates` + `gaussian_metrics` + `tw`
- **Output**：`{OUTDIR}/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv`
- **功能**：合并 3 路指标 + 加 comparison 元数据

### rule merge_all_window_pair_metrics
- **文件**：`rules/window_compare.smk`
- **Input**：所有 comparison 的 `window_pair_metrics.tsv`
- **Output**：`{OUTDIR}/window_compare/merged_window_pair_metrics.tsv`
- **功能**：合并所有 comparison 结果（awk 实现）

## 4. Checkpoint 使用说明

Pipeline 使用了 3 个 checkpoint 动态追踪产出：

| Checkpoint | 触发条件 | 动态产出 |
|-----------|---------|---------|
| `extract_and_balance_sample` | mcool URI 存在 | `.done` marker |
| `compute_gap_qc` | balanced cool 文件存在 | `gaps.tsv` + `window_qc.tsv` |
| `sliding_window` | `compute_gap_qc` checkpoint 完成 | `.done` marker |

下游规则通过 `checkpoints.<name>.get(...)` 引用 checkpoint，确保 Snakemake 在 checkpoint 真正产生文件后再执行。

**注意**：`compare_iced_windows_tw` 依赖的是 `extract_and_balance_sample` checkpoint（读取 iced.cool），而不是 `sliding_window` checkpoint。

## 5. Wildcard 使用模式

| Wildcard | 来源 | 取值范围 | 用于 |
|---------|------|---------|------|
| `{sample}` | checkpoint `extract_and_balance_sample` | 74 个样本（从 manifest.tsv） | extract_chrom, extract_expected, compute_oe_matrix, compute_gap_qc, sliding_window |
| `{comparison}` | `window_compare_pairs.yaml` | 36 个 comparison | window_compare 下的 5 个 rule |
| `{chrom}` | `config/chroms.yaml` | 26 条染色体 | 通过 expand() 批量生成输入 |
| `{stage}` | 硬编码 | `ln_oe`, `gap_filled`, `zscore`, `gaussian` | `sliding_window_h5` 输出 |

## 6. 关键设计决策

### checkpoint 双层结构（extract_chrom）
`extract_and_balance_sample` 既是 checkpoint（等待动态 mcool URI）又是执行 rule（shell 循环染色体）。.done marker 作为下游 checkpoint 的依赖追踪。

### checkpoint 双层结构（sliding_window）
`sliding_window` checkpoint 触发 Python 脚本执行；`sliding_window_h5` rule 仅声明 H5 文件路径（mtime 追踪），不重复执行计算。

### subprocess 子脚本模式
为避免 Snakemake NamedList 的 pickle 问题，跨文件调用统一使用 `subprocess.run([sys.executable, script_path, ...])`：
- `extract_expected_batch.py` 调用 `extract_expected.py`
- `compute_oe_matrix_batch.py` 调用 `compute_oe_matrix.py`
