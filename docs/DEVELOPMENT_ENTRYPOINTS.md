# DEVELOPMENT_ENTRYPOINTS

## 1. 新增模块接入位置概览

| 分析类型 | 推荐挂载位置 | 说明 |
|---------|------------|------|
| 新归一化方法 | `scripts/compute_oe_matrix.py` + `compute_oe_matrix_batch.py` | 在 O/E 计算前/后添加新归一化 stage |
| 新 gap fill 算法 | `scripts/sliding_window.py` — `_fill_gaps_2d_nearest` | 替换 `_build_2d_gap_masks` 或新增 stage |
| 新 QC 指标 | `scripts/compute_gap_qc.py` — `process_window()` | 添加新列到 window_qc.tsv |
| 新比较指标 | `scripts/compare_gaussian_smoothed_windows.py` — `compute_window_metrics()` | 添加新 metric dict key |
| 新 TW 指标 | `scripts/compare_iced_windows_tw.py` — `compute_tracy_widom_stat()` | 添加新统计量 |
| 新协变量 | `scripts/compute_window_pair_covariates.py` — vectorized 计算块 | 添加新列到 covariates.tsv |
| 新 merge 逻辑 | `scripts/merge_window_pair_metrics.py` | 在 merge_keys 后追加新列 |
| 全新 pipeline stage | 新建 `rules/new_stage.smk` + `scripts/new_stage.py`，在 Snakefile 添加 include | 最高侵入量 |

## 2. 不推荐修改的位置

| 位置 | 原因 |
|------|------|
| `Snakefile` — `rule all` inputs | 删改会破坏 pipeline 目标定义；Snakemake 不会检测到缺失的 rule |
| `Snakefile` — `CHROMS` / `SAMPLES` 加载逻辑 | 硬编码依赖关系，任何错误会导致整个 DAG 失效 |
| `rules/extract_chrom.smk` — `hictk` shell 命令 | 直接调 CLI，修改可能导致数据不一致 |
| `rules/compute_gap_qc.smk` — checkpoint 依赖链 | checkpoint 间依赖顺序敏感，错误会导致上游文件不可达 |
| `rules/sliding_window.smk` — checkpoint `sliding_window` | 改 checkpoint 名称会影响下游所有 rule 的 `checkpoints.sliding_window.get()` 调用 |
| `scripts/compute_window_pair_covariates.py` — sample 名称解析 | 依赖 `_rep1`/`_rep2` 字符串 hack，修改会破坏所有 comparison |
| `env/SNAP-HiC-Diff.yaml` | conda 环境锁定，版本变更需要重建环境 |

## 3. 新增 MSE/SSIM/Spearman 模块

**当前状态：已实现**，无需新增。

- **MSE**：`scripts/compare_gaussian_smoothed_windows.py` — `compute_window_metrics()` — `np.mean((vals_a - vals_b) ** 2)`
- **SSIM**：`scripts/compare_gaussian_smoothed_windows.py` — `compute_window_metrics()` — `skimage.metrics.structural_similarity`
- **Spearman**：`scripts/compare_gaussian_smoothed_windows.py` — `compute_window_metrics()` — `scipy.stats.spearmanr`

**扩展方向：**
- 配置化输出列：在 `config.yaml` 添加 `window_compare.metrics_whitelist: [MSE, SSIM, Spearman, NCC]`，在 script 中过滤
- SSIM 窗口参数化：当前 `data_range=1.0` 硬编码，改为 `config["smoothing"]["ssim_data_range"]`
- 添加 KS 检验：在同一脚本中添加 `scipy.stats.ks_2samp` 指标

## 4. 新增 Tracy-Widom 相关模块

**当前状态：已实现（beta=1, real symmetric matrix）**。

`scripts/compare_iced_windows_tw.py` 使用 `TracyWidom(beta=1)`。

**扩展方向：**
- 支持 `beta=2`（complex Hermitian 矩阵，对应某些 Hi-C 归一化方法）
- 添加 eigenvalue 谱统计量（最大/最小特征值、谱距）
- FDR 校正的 p-value 多重检验
- 将 Diffmat 计算改为 log-scale：在 `compare_iced_windows_tw.py` 的 `process_single_window()` 中 `np.log(mat_a / mat_b)` 替换当前 `mat_a / mat_b`

## 5. Config 扩展指南

### 新增顶级配置项

在 `config/config.yaml` 末尾添加：

```yaml
new_feature:
  param1: 123
  param2: "value"
```

在 `Snakefile` 中通过 `config["new_feature"]["param1"]` 访问。

### 在 Rule 中使用新配置

在 rule 的 `params` 块中引用：

```snakemake
rule new_rule:
    params:
        my_param = config["new_feature"]["param1"],
    shell:
        "python {input.script} --param {params.my_param}"
```

### 配置验证

当前 pipeline 无显式配置验证。建议新增 `scripts/validate_config.py` 在 `validate_inputs` rule 中调用：

```python
def validate_config(cfg):
    required_keys = ["resolution", "window_size", "step_size"]
    for key in required_keys:
        if key not in cfg:
            raise ValueError(f"Missing required config key: {key}")
    if cfg["window_size"] < cfg["resolution"] * 5:
        raise ValueError("window_size must be >= 5 * resolution")
```

## 6. 文件修改规范

### Rule 文件规范（Snakemake v8.30.0）

```snakemake
rule <name>:
    input:
        # ...
    output:
        # ...
    log:
        # ...
    benchmark:
        "{OUTDIR}/benchmark/<name>/<name>.tsv"
    threads:
        config["..."]
    params:
        # ...
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        # use ${{VAR}} for shell variables, {{wildcard}} for escaping
        """
```

关键规范：
- **字段顺序**：input → output → log → benchmark → threads → params → conda → shell
- **benchmark**：每个 rule 必须有
- **shell 模板**：用 `r"""..."""` raw string，`{{wildcard}}` 转义
- **bash 变量**：`${{VAR}}` 格式（双大括号外层，单大括号变量名）
- **f-string 禁止**：Snakemake v8 不支持在 rule 块内 f-string

### Script 文件规范

所有脚本应同时支持 **snakemake 对象**（用于 rule 内联）和 **argparse**（用于 subprocess 调用）：

```python
try:
    INPUT = snakemake.input[0]
    OUTPUT = snakemake.output[0]
    LOG = snakemake.log[0]
except NameError:
    args = parse_args()
    INPUT = args.input
    OUTPUT = args.output
    LOG = args.log
```

建议：`extract_expected.py`、`chunk_windows.py` 已经示范此模式。

### Config 文件规范

- `config/config.yaml`：顶级参数（resolution、window_size、parallel 等）
- `config/chroms.yaml`：染色体列表（独立文件）
- `data/window_compare_pairs.yaml`：样本比较对配置（独立文件）
- 新增 config 遵循相同分组原则

## 7. 风险点和注意事项

| 修改类型 | 风险等级 | 缓解措施 |
|---------|---------|---------|
| 修改 `config/resolution` | **高** | 整个 pipeline 需重跑；输出目录包含 resolution 维度 |
| 修改 checkpoint 名称 | **高** | 下游所有 `checkpoints.<name>.get()` 引用需同步更新 |
| 修改 `smoothing.sigma` | 中 | 仅需重跑 `sliding_window` 和 `window_compare` stages |
| 修改 `qc` 阈值 | 中 | 仅需重跑 `compute_gap_qc` 和后续 stages |
| 修改 `window_compare` 配置 | 中 | 仅需重跑 `window_compare` rule |
| 添加新 script | 低 | 需要同步更新 rule 的 `input.script` 路径 |
| 修改 sliding_window 4 stage 输出 | 低 | 需同步更新 `sliding_window_h5` 的 `stage` wildcard |

## 8. 测试建议

### 本地测试

```bash
# Dry-run（不执行 job）
snakemake -n -p

# 单 rule 测试
snakemake --until validate_inputs --forcerun validate_inputs

# 指定 comparison 测试
snakemake --until window_compare/H1398_rep1_vs_H1398_rep2/window_pair_metrics.tsv \
    --forcerun window_compare
```

### 单元测试

建议为关键函数添加 pytest：

```
tests/
├── test_global_robust_zscore.py
├── test_detect_breakpoint_by_nan_density.py
├── test_compute_window_metrics.py
├── test_tracy_widom.py
└── test_merge_metrics.py
```

### 集成测试

```bash
# 1. 用小样本子集测试
# 2. 对比已知 ground truth 的输出
# 3. 验证 checkpoint 恢复（删除中间文件后重跑）
snakemake --until compute_oe_matrix --forcerun
snakemake --until window_compare
```

## 9. 文档更新建议

### 必须更新

| 文件 | 何时更新 |
|------|---------|
| `docs/CODEBASE_OVERVIEW.md` | 新增指标、新配置参数、新算法 |
| `docs/SNAKEMAKE_PIPELINE_MAP.md` | 新增 rule、checkpoint、wildcard |
| `docs/RULE_SCRIPT_CONFIG_MAP.md` | 新增 script、修改参数流向 |
| `docs/DEVELOPMENT_ENTRYPOINTS.md` | 新增接入点、风险点 |
| `config/config.yaml` | 添加注释说明新参数 |

### 可选更新

| 文件 | 何时更新 |
|------|---------|
| `README.md` | 新增功能说明、新的使用示例 |
| `docs/project_design.md` | 新增设计决策记录 |
| `docs/rule_index.md` | 新增 rule 文档 |
| `docs/snakemake_v8_30_0_standard.md` | Snakemake 规范更新 |
