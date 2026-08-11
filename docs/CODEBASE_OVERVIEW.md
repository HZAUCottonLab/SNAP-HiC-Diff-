# CODEBASE_OVERVIEW

## 1. 项目核心目标

SNAP-HiC-Diff 是一个基于 Snakemake 的 Hi-C 差异分析流程，从 `.mcool`（多分辨率接触图）出发，通过 O/E 归一化、滑动窗口分解、gap 检测和统计比较，检测基因组中的差异染色质相互作用（Differential Chromatin Interactions）。

- 输入：ICE 平衡后的 Hi-C 接触图（`.mcool`）
- 核心产物：窗口级别的差异指标（`window_pair_metrics.tsv`）
- 指标：MSE、SSIM、Spearman、NCC（高斯平滑矩阵）；Tracy-Widom 统计量（ICE 矩阵）

## 2. 目录结构

| 路径 | 作用 | 关键输入/输出 | 备注 |
|------|------|--------------|------|
| `Snakefile` | 入口编排（含 `compute_IS`、`meltron`） | config + manifest | `rule all` 定义最终产物 |
| `run_snakemake.sh` | 本地/HPC 启动脚本 | — | 支持 `--dsub` 走 HPC |
| `config/config.yaml` | 主配置：分辨率、窗口、QC 阈值 | — | 74 个参数 |
| `config/chroms.yaml` | 染色体列表 | — | 26 条（HC04_A01-A13, D01-D13）|
| `data/manifest.tsv` | 74 个样本元数据 | `.mcool` URI | 样本 = rep1 + rep2 配对 |
| `data/window_compare_pairs.yaml` | 36 对样本比较配置 | — | 36 个 comparison |
| `env/SNAP-HiC-Diff.yaml` | Conda 环境 | — | snakemake>=8.0, hictk, cooltools |
| `rules/*.smk` | Snakemake 规则定义 | checkpoint/rule | 含 per-sample CWorld IS 与 Meltron 分支 |
| `scripts/*.py` | Python 执行脚本（13 个）| argparse / snakemake 双模式 | subprocess 调子脚本避免 pickle |
| `jupyter/` | 探索性 notebook | — | 不参与主流程 |

## 3. 数据流

```
manifest.tsv
    │
    ▼
validate_inputs.smk → manifest.validated.tsv
    │
    ▼
build_windows.smk → windows_{RESOLUTION}.parquet
    │
    ▼
extract_chrom.smk (checkpoint × 74 samples)
    │ 每个 sample: mcool → balanced/{sample}/{sample}__{chrom}_iced.cool
    │
    ├─► extract_expected.smk
    │       Step1: cooltools expected-cis → {chrom}_cooltools.tsv
    │       Step2: extract_expected_batch.py → {chrom}_expected.npz
    │
    ├─► compute_oe_matrix.smk
    │       compute_oe_matrix_batch.py → {chrom}_oe.npz
    │
    └─► compute_gap_qc.smk (checkpoint × 74 samples)
            compute_gap_qc.py → gaps.tsv + window_qc.tsv
                │
                ▼
            sliding_window.smk (checkpoint × 74 samples)
                sliding_window.py → ln_oe.h5, gap_filled.h5, zscore.h5, gaussian.h5
                    │
                    ▼
                window_compare.smk
                    ├─ compute_window_pair_covariates → {cmp}_covariates.tsv
                    ├─ compare_gaussian_smoothed_windows → {cmp}_gaussian_metrics.tsv
                    ├─ compare_iced_windows_tw → {cmp}_tw.tsv
                    └─ merge_window_pair_metrics → {cmp}_window_pair_metrics.tsv

extract_chrom.smk (configured-resolution balanced `iced.cool` per sample × chromosome; default 4 kb)
    │
    ▼
compute_IS.smk: `compute_IS_chrom` (26 four-CPU jobs per sample)
    │ each job: one chromosome `.cool` → `tmp/cworld_matrix/{sample}/` matrix
    │          → ten sequential CWorld scales (32000, 60000, 92000, 120000,
    │            152000, 180000, 212000, 240000, 272000, 300000 bp)
    │          → ten `.insulation.bedGraph` files + one manifest in
    │             `results/cworld_insulation/{sample}/`
    ▼
compute_IS.smk: `merge_IS` (one four-CPU job per sample)
    │ consumes all 26 manifests through convert_cworld_IS_to_meltron_IS.py
    ▼
results/IS/{sample}/{sample}_IS.tsv
    │ `chrom<TAB>viewpoint<TAB>IS_distance<TAB>IS`
    ▼
meltron.smk (per comparison; reads the unchanged TSV contract for both samples)
```

`compute_IS.smk` 由两个普通 rule 组成，而非旧的 per-sample `compute_IS` checkpoint：`compute_IS_chrom` 按 sample × chromosome 运行，使用 `meltron.cworld_is.source_resolution_bp` 配置的输入分辨率（整数且不低于 1 kb；默认 4 kb）并在一个 job 中顺序计算十个 scale；`merge_IS` 仅在该 sample 的 26 个 manifest 均存在时生成最终 TSV。任一 chromosome 失败只阻塞所属 sample 的 merge 及其下游比较。

## 4. 关键文件说明

### Snakemake 相关

| 文件 | 作用 |
|------|------|
| `Snakefile` | include 顺序：validate_inputs → build_windows → extract_chrom → extract_expected → compute_oe_matrix → compute_gap_qc → sliding_window → window_compare → compute_IS → meltron |
| `rules/extract_chrom.smk` | **checkpoint** `extract_and_balance_sample`；每个 sample 一个 job，shell 循环 26 条染色体；`hictk load + balance` |
| `rules/extract_expected.smk` | 2 个 rule：`extract_expected_tsv`（cooltools）和 `extract_expected_npz`（子脚本）|
| `rules/compute_gap_qc.smk` | **checkpoint** `compute_gap_qc`；调用 `compute_gap_qc.py` |
| `rules/sliding_window.smk` | **checkpoint** `sliding_window` + `rule sliding_window_h5`（声明 H5 输出）|
| `rules/compute_IS.smk` | 两个 rule：`compute_IS_chrom` 为每个 sample × chromosome 使用 `meltron.cworld_is.source_resolution_bp` 配置的 `iced.cool`（整数且不低于 1 kb；默认 4 kb），顺序生成十个 scale（32000, 60000, 92000, 120000, 152000, 180000, 212000, 240000, 272000, 300000 bp）的 bedGraph 和 manifest；`merge_IS` 汇总 26 个 manifest，写出 `results/IS/{sample}/{sample}_IS.tsv` |
| `rules/meltron.smk` | `meltron_per_comparison`；读取双方不变的四列 IS TSV 契约和 KNN-enriched windows |

### Python 脚本

| 脚本 | 核心功能 |
|------|---------|
| `validate_inputs.py` | `hictk validate` 检查 mcool，提取染色体列表 |
| `build_windows.py` | `hictk dump chroms`，生成 parquet 窗口坐标 |
| `compute_gap_qc.py` | 梯度断点检测 + 动态 P20 QC 阈值；import `detect_breakpoint_by_nan_density` |
| `sliding_window.py` | 4 stage：`ln_oe` → `gap_filled`（2D 插值）→ `zscore`（global robust MAD）→ `gaussian`（高斯平滑）|
| `compare_gaussian_smoothed_windows.py` | MSE、SSIM、Spearman、NCC（从 gaussian.h5）|
| `compare_iced_windows_tw.py` | Tracy-Widom 检验（从 iced.cool）；`normdiff_by_mean_sd` 归一化 |
| `detect_breakpoint_by_nan_density.py` | 5×5 滑动窗口 NaN 密度断点检测；被 `compute_gap_qc.py` 调用 |
| `compute_window_pair_covariates.py` | EC（调和均值 CPM）、IB（|log2 比率|）、SP（nz_frac 均值）|
| `merge_window_pair_metrics.py` | 合并 3 路指标 + 加 comparison 元数据 |
| `convert_cworld_IS_to_meltron_IS.py` | 由 `merge_IS` 读取每个 chromosome 的 manifest 并汇总 native CWorld insulation bedGraph，写出 Meltron 兼容的 `chrom, viewpoint, IS_distance, IS` TSV |

## 5. Config 关键参数

| 参数路径 | 默认值 | 作用 |
|---------|--------|------|
| `resolution` | `4000` | Hi-C bin 大小（bp）|
| `meltron.cworld_is.source_resolution_bp` | `4000` | CWorld/Meltron 输入 bin 大小（bp；整数且不低于 1000，必须与 `.cool` 一致）|
| `window_size` | `300000` | 滑动窗口宽度（bp）|
| `step_size` | `32000` | 滑动步长（bp）|
| `smoothing.sigma` | `1.0` | 高斯平滑 sigma |
| `gap.gradient_threshold` | `10.0` | 断点梯度比阈值 |
| `qc.dynamic_raw_sum_cpm_pct` | `10` | 动态 P20 阈值百分位（0=禁用）|
| `tracy_widom.min_nbin` | `10` | TW 最小有效 bin 数 |
| `tracy_widom.na_fraction` | `0.5` | TW 行 NaN 过滤阈值 |
| `window_compare.gaussian_metrics.n_jobs` | `80` | 高斯比较并行 job 数 |
| `window_compare.tw.n_jobs` | `80` | TW 检验并行 job 数 |

## 6. 后续开发最需要理解的核心模块

1. **`sliding_window.py`** — gap fill 算法（`_fill_gaps_2d_nearest`）、`_build_2d_gap_masks` 的 1D→2D 转换语义；`global_robust_zscore` 是全局 MAD z-score（**不是** distance-stratified）
2. **`compute_gap_qc.py`** — 动态 P20 阈值逻辑（per-chrom 应用）、`detect_breakpoint_by_nan_density` 的 5×5 扫描断点算法
3. **`compare_iced_windows_tw.py`** — `normdiff_by_mean_sd`（按距离分层归一化）和 Tracy-Widom 统计量的计算公式
4. **`compute_window_pair_covariates.py`** — sample 名称解析的 `_rep1`/`_rep2` 字符串 hack（`split("_rep1")[0]`）

## 7. 扩展建议

### MSE/SSIM/Spearman
已实现于 `compare_gaussian_smoothed_windows.py`，**无需新增**。

扩展方向：
- 增加 SSIM 的窗口大小参数（当前硬编码）
- 添加 `window_compare.metrics_whitelist` 配置项控制输出列

### Tracy-Widom
已实现于 `compare_iced_windows_tw.py`，**无需新增**。

扩展方向：
- 支持 beta=2（复数矩阵）；当前硬编码 `beta=1`
- 添加 FDR 校正的 p-value 输出

### 新增指标
最佳挂载点：`rules/window_compare.smk` + `scripts/compare_NEW_metrics.py`，参考 `compare_gaussian_smoothed_windows.smk` 的 DAG 模式：
1. 依赖 `compute_window_pair_covariates.output.covariates`
2. 读取对应 stage 的 H5 文件
3. 在 `merge_window_pair_metrics.py` 中追加 merge 逻辑

## 5a. Chromosome Naming Convention — CRITICAL

**Two naming schemes exist in this project:**

| Scheme | Prefix | Example | Used by |
|--------|--------|---------|---------|
| Old | `Ghir_` | `Ghir_A01` | `Ghir_Pyb` mcool files |
| Current | `HC04_` | `HC04_A01` | `compare_mtx_v2` mcool files |

**`config/chroms.yaml` must match the actual mcool files.** Mismatch causes all 36 `extract_and_balance_sample` jobs to fail with `exit 1`.

**How to verify:**

```bash
# Check what names the mcool files actually use
# (extract mcool path from data/manifest_mcool.tsv first column)
hictk dump /path/to/sample.mcool::/resolutions/4000 --table chroms | head -3

# Compare with chroms.yaml
cat config/chroms.yaml
```

**`dsub` runs from git HEAD, NOT the working tree.** Working-tree-only changes to `config/chroms.yaml` are invisible to dsub. Always commit before submitting:

```bash
git add config/chroms.yaml && git commit -m "fix: align chroms.yaml"
dsub -s run_snakemake.sh
```

See `docs/BUG_FIX_GUIDE.md` for the full fix procedure.

## 6. 需要人工确认

1. `extract_expected.py` 支持 `snakemake` 对象（与 `chunk_windows.py` 相同的双模式设计），但从未在 rule 中以 snakemake 方式调用——是否应删除 snakemake 分支？
2. `compute_gap_qc.py` 的动态 P20 阈值设计为 per-chromosome 级别应用，但 QC PASS/Fail 状态在每个样本内部不一致——是否有全局统一阈值的需求？
3. 36 对 `window_compare_pairs.yaml` 目前是 rep1 vs rep2 的同样本重复，是否需要跨样本比较（如 H1398_rep1 vs L128_rep1）？
