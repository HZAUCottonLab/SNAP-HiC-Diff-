# RULE_SCRIPT_CONFIG_MAP

## 1. Rule → Script 调用关系

| Rule | Script | 调用方式 | 功能 |
|------|--------|---------|------|
| `validate_inputs` | `scripts/validate_inputs.py` | `python {script} --input ... --output ...` | mcool 文件验证 |
| `build_windows` | `scripts/build_windows.py` | `python {script} --input ... --output ...` | 生成基因组窗口 parquet |
| `extract_and_balance_sample` | — | 直接 shell 命令 `hictk dump/load/balance` | 提取并 ICE 平衡染色体 |
| `extract_expected_tsv` | — | 直接 shell 命令 `cooltools expected-cis` | 生成期望值 TSV |
| `extract_expected_npz` | `scripts/extract_expected_batch.py` | `python {script} --chroms ... --inputs ... --outputs ...` | TSV → NPZ 批量转换 |
| `compute_oe_matrix` | `scripts/compute_oe_matrix_batch.py` | `python {script} --chroms ... --inputs ... --expected ...` | O/E 归一化批量计算 |
| `compute_gap_qc` | `scripts/compute_gap_qc.py` | `python {script} --chroms ... --balanced-cool-files ... --outdir ...` | Gap 检测 + 窗口 QC |
| `sliding_window` | `scripts/sliding_window.py` | `python {script} --chroms ... --gap-tsv-files ... --ln-oe-files ... --outdir ...` | Gap fill → zscore → Gaussian |
| `sliding_window_h5` | — | 无（仅声明输出文件路径）| Snakemake mtime 追踪 |
| `compute_window_pair_covariates` | `scripts/compute_window_pair_covariates.py` | `python {script} --comparison ... --sample-a ... --sample-b ...` | EC、IB、SP 协变量 |
| `compare_gaussian_smoothed_windows` | `scripts/compare_gaussian_smoothed_windows.py` | `python {script} --comparison ... --covariates ...` | MSE、SSIM、Spearman、NCC |
| `compare_iced_windows_tw` | `scripts/compare_iced_windows_tw.py` | `python {script} --comparison ... --covariates ...` | Tracy-Widom 检验 |
| `merge_window_pair_metrics` | `scripts/merge_window_pair_metrics.py` | `python {script} --comparison ... --covariates ... --gaussian ... --tw ...` | 合并 3 路指标 |
| `merge_all_window_pair_metrics` | — | 直接 shell 命令 `awk` | 合并所有 comparison |

## 2. Script → Rule 反向调用关系

| Script | 被调用的 Rule(s) |
|--------|-----------------|
| `scripts/validate_inputs.py` | `rule validate_inputs` |
| `scripts/build_windows.py` | `rule build_windows` |
| `scripts/extract_expected.py` | `scripts/extract_expected_batch.py`（subprocess）→ `rule extract_expected_npz` |
| `scripts/extract_expected_batch.py` | `rule extract_expected_npz` |
| `scripts/compute_oe_matrix.py` | `scripts/compute_oe_matrix_batch.py`（subprocess）→ `rule compute_oe_matrix` |
| `scripts/compute_oe_matrix_batch.py` | `rule compute_oe_matrix` |
| `scripts/compute_gap_qc.py` | `checkpoint compute_gap_qc` |
| `scripts/sliding_window.py` | `checkpoint sliding_window` |
| `scripts/compute_window_pair_covariates.py` | `rule compute_window_pair_covariates` |
| `scripts/compare_gaussian_smoothed_windows.py` | `rule compare_gaussian_smoothed_windows` |
| `scripts/compare_iced_windows_tw.py` | `rule compare_iced_windows_tw` |
| `scripts/merge_window_pair_metrics.py` | `rule merge_window_pair_metrics` |
| `scripts/detect_breakpoint_by_nan_density.py` | `scripts/compute_gap_qc.py`（import）|
| `scripts/chunk_windows.py` | **未被任何 rule 调用**（standalone utility）|

## 3. Config → Rule → Script 参数流向

### resolution / window_size / step_size
| Config Key | Value | Rule | Script Param |
|-----------|-------|------|-------------|
| `resolution` | `4000` | `build_windows`, `compute_gap_qc`, `sliding_window` | `--resolution` |
| `window_size` | `300000` | `build_windows`, `compute_gap_qc`, `sliding_window` | `--window-size` |
| `step_size` | `32000` | `build_windows`, `compute_gap_qc`, `sliding_window` | `--step-size` |

### smoothing
| Config Key | Value | Rule | Script Param |
|-----------|-------|------|-------------|
| `smoothing.sigma` | `1.0` | `sliding_window` | `--sigma` |

### gap_detection
| Config Key | Value | Rule | Script Param |
|-----------|-------|------|-------------|
| `gap.gradient_threshold` | `10.0` | `compute_gap_qc` | `--gradient-threshold` |
| `gap.min_low_ratio` | `0.001` | `compute_gap_qc` | `--min-low-ratio` |

### qc
| Config Key | Value | Rule | Script Param |
|-----------|-------|------|-------------|
| `qc.dynamic_raw_sum_cpm_pct` | `10` | `compute_gap_qc` | `--dynamic-raw-sum-cpm-pct` |
| `qc.max_bad_bin_frac` | `0.5` | `compute_gap_qc` | `--max-bad-bin-frac` |
| `nan_density_breakpoint.abnormal_threshold` | `0.60` | `compute_gap_qc` | `--nan-abnormal-threshold` |

### tracy_widom
| Config Key | Value | Rule | Script Param |
|-----------|-------|------|-------------|
| `tracy_widom.min_nbin` | `10` | `compare_iced_windows_tw` | `--min-nbin` |
| `tracy_widom.na_fraction` | `0.5` | `compare_iced_windows_tw` | `--filter-threshold` |

### window_compare.covariates
| Config Key | Value | Rule | Script Param |
|-----------|-------|------|-------------|
| `window_compare.covariates.cpm_column` | `raw_sum_cpm` | `compute_window_pair_covariates` | `--cpm-column` |
| `window_compare.covariates.nz_column` | `nz_frac` | `compute_window_pair_covariates` | `--nz-column` |
| `window_compare.covariates.pseudocount_cpm` | `0.01` | `compute_window_pair_covariates` | `--pseudocount` |
| `window_compare.covariates.weights.EC` | `0.33` | `compute_window_pair_covariates` | `--weight-ec` |
| `window_compare.covariates.weights.IB` | `0.33` | `compute_window_pair_covariates` | `--weight-ib` |
| `window_compare.covariates.weights.SP` | `0.33` | `compute_window_pair_covariates` | `--weight-sp` |

### parallelization
| Config Key | Value | Rule | Script Param |
|-----------|-------|------|-------------|
| `hictk.threads` | `10` | `extract_and_balance_sample`, `extract_expected_tsv` | `--threads` (hictk/cooltools CLI) |
| `parallel.threads_per_job` | `80` | `validate_inputs`, `build_windows`, `compute_oe_matrix`, `compute_gap_qc`, `sliding_window`, `compute_window_pair_covariates` | `--n-jobs` (joblib) |
| `window_compare.gaussian_metrics.n_jobs` | `80` | `compare_gaussian_smoothed_windows` | `--n-jobs` |
| `window_compare.tw.n_jobs` | `80` | `compare_iced_windows_tw` | `--n-jobs` |

## 4. 输入输出关系 (I/O Chain)

### Stage 1: Validate
```
manifest.tsv → validate_inputs.py → manifest.validated.tsv
```

### Stage 2: Build Windows
```
manifest.validated.tsv + hictk dump chroms
    → build_windows.py → windows_{RESOLUTION}.parquet
```

### Stage 3: Extract & Balance
```
mcool::/resolutions/{RESOLUTION}
    → hictk dump/load/balance
    → balanced/{sample}/{sample}__{chrom}_iced.cool
    → raw/{sample}/{sample}__{chrom}.raw.cool
```

### Stage 4: Expected
```
iced.cool → cooltools expected-cis
    → extract_expected_tsv: {chrom}_cooltools.tsv
    → extract_expected_batch.py: {chrom}_expected.npz
```

### Stage 5: O/E
```
iced.cool + expected.npz
    → compute_oe_matrix_batch.py
    → {sample}__{chrom}_oe.npz (sparse, symmetric)
```

### Stage 6: Gap QC + Sliding Window
```
iced.cool
    → compute_gap_qc.py
        → {sample}__{chrom}_gaps.tsv
        → {sample}__{chrom}_window_qc.tsv (with 8 breakpoint columns)

iced.cool + gaps.tsv + oe.npz
    → sliding_window.py
        → {sample}__{chrom}_ln_oe.h5
        → {sample}__{chrom}_gap_filled.h5
        → {sample}__{chrom}_zscore.h5
        → {sample}__{chrom}_gaussian.h5
```

### Stage 7: Window Compare
```
window_qc.tsv (sample_a) + window_qc.tsv (sample_b)
    → compute_window_pair_covariates.py
        → {comparison}_covariates.tsv (EC, IB, SP)

gaussian.h5 (sample_a) + gaussian.h5 (sample_b)
    → compare_gaussian_smoothed_windows.py
        → {comparison}_gaussian_metrics.tsv (MSE, SSIM, Spearman, NCC)

iced.cool (sample_a) + iced.cool (sample_b)
    → compare_iced_windows_tw.py
        → {comparison}_tw.tsv (tw_statistic, tw_pvalue)

covariates + gaussian + tw
    → merge_window_pair_metrics.py
        → {comparison}_window_pair_metrics.tsv (final per-comparison)
```

## 5. 影响范围分析

### 修改某 Rule/Script/Config 的影响

| 修改位置 | 影响范围 | 风险 |
|---------|---------|------|
| `config/config.yaml` — `resolution` | `build_windows`, `compute_gap_qc`, `sliding_window`, `extract_chrom`, `compare_iced_windows_tw` | 高：需要重新运行整个 pipeline |
| `config/config.yaml` — `window_size` / `step_size` | `build_windows`, `compute_gap_qc`, `sliding_window` | 高：窗口定义全改 |
| `config/config.yaml` — `smoothing.sigma` | `sliding_window` | 中：仅需重新跑 sliding window stage |
| `config/config.yaml` — `window_compare.covariates` | `compute_window_pair_covariates` | 中：仅需重新跑 comparison stage |
| `rules/window_compare.smk` | 所有 window_compare sub-rules | 中：影响 DAG 结构 |
| `scripts/sliding_window.py` | `checkpoint sliding_window` | 中：gap fill / zscore / gaussian 算法变更 |
| `scripts/compute_gap_qc.py` | `checkpoint compute_gap_qc` | 高：QC 阈值变化影响所有 PASS/FAIL 判断 |
| `scripts/detect_breakpoint_by_nan_density.py` | `compute_gap_qc.py` | 中：断点检测逻辑变更 |

### 新增指标的影响评估

如果新增指标（如 Kolmogorov-Smirnov 检验、DIN 归一化等），需要在以下位置修改：

| 新增内容 | 需要修改的文件 |
|---------|--------------|
| 新指标脚本 | `scripts/compare_NEW.py` |
| 新 rule | `rules/window_compare.smk` |
| 新 rule 调用 | `Snakefile` — `rule all` inputs |
| 新 merge 逻辑 | `scripts/merge_window_pair_metrics.py` |
| 新 config 参数 | `config/config.yaml` — `window_compare` section |

## 6. 脚本间依赖关系

```
scripts/detect_breakpoint_by_nan_density.py
    ↑
    │ (import)
    │
scripts/compute_gap_qc.py
    │ gap.tsv + window_qc.tsv
    ↓
scripts/sliding_window.py
    │ gaussian.h5
    ↓
scripts/compare_gaussian_smoothed_windows.py
    │ (reads gaussian.h5)
    ↑
scripts/compute_window_pair_covariates.py
    │ (reads window_qc.tsv)
    ↓
scripts/merge_window_pair_metrics.py

scripts/compare_iced_windows_tw.py
    │ (reads iced.cool directly)
    ↑
scripts/merge_window_pair_metrics.py

scripts/extract_expected_batch.py
    │ (subprocess call)
    ↓
scripts/extract_expected.py
    │
scripts/compute_oe_matrix_batch.py
    │ (subprocess call)
    ↓
scripts/compute_oe_matrix.py
```

### 关键：无跨脚本 import
除了 `compute_gap_qc.py → detect_breakpoint_by_nan_density.py` 的 import 关系外，所有脚本均为独立执行单元，通过文件路径传递数据。
