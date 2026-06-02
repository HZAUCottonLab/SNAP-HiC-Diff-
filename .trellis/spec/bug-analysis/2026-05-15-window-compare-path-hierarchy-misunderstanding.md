# Bug: Script Path Hierarchy Misunderstanding — window_compare Pipeline

> 日期: 2026-05-15
> 标签: snakemake, python, path, window_compare, integration
> 状态: 已解决

## 问题描述

Snakemake pipeline 运行成功（exit code 0），但所有输出文件只有表头、无数据行。日志显示 Rule 1 的 `compute_window_pair_covariates` 对所有 26 个染色体都输出 `"Skipping: QC file not found"`。

## 报错信息

```
# 日志 (results/logs/window_compare/L128_rep1_vs_L128_rep2/compute_covariates.log):
[L128_rep1_vs_L128_rep2] Starting covariate computation across 26 chromosomes
  Processing HC04_A01
    Skipping HC04_A01: QC file not found   ← 26 个染色体全部跳过
  Processing HC04_A02
    Skipping HC04_A02: QC file not found
  ...
  No windows across any chromosome, writing empty output
```

## 根因分析

### 直接原因

脚本 `compute_window_pair_covariates.py` 把 **输出文件路径** 当作 **父目录** 来拼接 QC 文件路径：

```python
# 脚本拼接的路径
qc_a_path = f"{outdir}/windows/{sample_a}/{sample_a}__{chrom}_window_qc.tsv"

# 实际传入的 outdir = covariates.tsv 文件路径
# results/window_compare/L128_rep1_vs_L128_rep2/L128_rep1_vs_L128_rep2_covariates.tsv

# 拼接结果 = 错误路径
# results/window_compare/L128_rep1_vs_L128_rep2/L128_rep1_vs_L128_rep2_covariates.tsv/windows/L128_rep1/...
#                                       ↑ 文件名被当作目录名 ↑
```

### 正确的路径层级

```
outdir = results/window_compare/{comparison}/{comparison}_covariates.tsv
parent      = results/window_compare/{comparison}/          (1 级)
parent.parent    = results/window_compare/                      (2 级)
parent.parent.parent = results/                                 (3 级) ← 正确

results_dir = results/
qc_path = results/windows/{sample}/{sample}__{chrom}_window_qc.tsv  ✓ 存在
```

### 为什么前三轮都没发现

此 bug 隐藏在三层代码问题之下：

| 轮次 | 报错 | 根因 | 脚本是否执行到路径拼接 |
|------|------|------|----------------------|
| 1 | `NameError: name '_os' is not defined` | 重构时 `import os` 移出了函数 | ❌ 第 1 行就 crash |
| 2 | `NameError: name 'chrom' is not defined` | merge 脚本 copy-paste 残留变量 | ❌ 第 1 行就 crash |
| 3 | 0 行输出，`QC file not found` | `parent.parent` 只向上 2 级而非 3 级 | ✓ 路径错误但脚本继续执行 |

直到层 1-2 的 crash 全部修复后，层 3 的静默失败才暴露出来。

## 解决方案

```python
# outdir = results/window_compare/{comparison}/{comparison}_covariates.tsv
# parent = results/window_compare/{comparison}/
# parent.parent = results/window_compare/
# parent.parent.parent = results/
results_dir = Path(outdir).parent.parent.parent
qc_a_path = f"{results_dir}/windows/{sample_a}/{sample_a}__{chrom}_window_qc.tsv"
```

修复后结果：
```
H1398: 61,390 行
L128:  59,441 行
L134:  59,293 行
```

## 相关文件

- `scripts/compute_window_pair_covariates.py` — 主 bug 所在
- `scripts/compare_gaussian_smoothed_windows.py` — 类似路径构造，已正确（但注释错误：写 `.., ".."` 实为 2 级）
- `scripts/compare_iced_windows_tw.py` — 类似路径构造，已正确
- `scripts/merge_window_pair_metrics.py` — 无路径拼接问题（纯合并）

## 经验教训

1. **Snakemake 成功 ≠ 脚本逻辑正确**：静默失败（0 行输出）更难发现
2. **多层 bug 叠加时先修 crash 再查逻辑**：层 1-2 的 crash 掩盖了层 3 的静默失败
3. **路径构造必须有注释标明层级数**：`# parent.parent.parent = results/` 可以防止此类 bug
4. **integration test 必要性**：dry-run 和单元测试都不执行脚本实际路径拼接逻辑
