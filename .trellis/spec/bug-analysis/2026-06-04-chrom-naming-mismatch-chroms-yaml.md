# Bug Analysis: extract_and_balance_sample 染色体名不匹配导致全量失败

**Date**: 2026-06-04
**Severity**: CRITICAL
**Status**: FIXED
**Dsub Job**: 36930160 → 36931564

## Summary

36 个样本的 `extract_and_balance_sample` 全部失败，同一错误：

```
[ERROR] Chromosome Ghir_A01 not found in /path/to/L282.mcool::/resolutions/4000
exit 1
```

根因：`config/chroms.yaml` 的 committed 版本使用 `Ghir_*` 染色体名，working tree 已改为 `HC04_*` 但未 commit。Snakemake 从 committed 版本读取配置，`hictk dump` 无法在 mcool 文件中找到对应染色体。

## Root Cause

### 直接原因

`config/chroms.yaml` 的 committed 版本（git HEAD）中染色体名为 `Ghir_*`：

```yaml
chroms:
  - Ghir_A01
  - Ghir_A02
  # ...
```

但 working tree 中已修改为 `HC04_*`：

```yaml
chroms:
  - HC04_A01
  - HC04_A02
  # ...
```

Snakemake 从 committed 版本读取配置，通过 `params.chroms` 传入 shell block：

```snakemake
# extract_chrom.smk (正确模式，已参数化)
params:
    chroms = " ".join(CHROMS)   # CHROMS = config["chroms"]

shell:
    for CHR in {params.chroms}; do
        hictk dump "$URI" --table chroms \
            | awk -v chr="$CHR" '$1 == chr' > "$CHROMSIZES"
        if [ ! -s "$CHROMSIZES" ]; then
            echo "[ERROR] Chromosome $CHR not found in $URI"
            exit 1
        fi
```

`hictk dump` 输出的实际染色体名是 `HC04_A01`，awk 匹配 `Ghir_A01` 失败 → CHROMSIZES 为空 → `exit 1`。

### 间接原因

**染色体命名约定未在文档中显式声明。**

本项目在不同阶段使用过两套染色体命名体系：

| 数据来源 | 命名 | 备注 |
|---------|------|------|
| 旧 mcool 文件（Ghir_Pyb 项目） | `Ghir_A01` ~ `Ghir_D13` | 早期 convert_to_mcool 阶段 |
| 新 mcool 文件（compare_mtx_v2） | `HC04_A01` ~ `HC04_D13` | 当前 pipeline 数据 |
| `config/chroms.yaml` (committed) | `Ghir_*` | commit 时未同步更新 |
| `config/chroms.yaml` (working tree) | `HC04_*` | 用户手动修改但未 commit |

`rules/extract_chrom.smk` 本身代码正确（使用 `{params.chroms}` 而非硬编码），但依赖 `chroms.yaml` 的值。染色体命名约定只在 working tree 中修改，未同步到 git commit。

## Symptom

```
TASK_INFO:
    JOB_NAME: SNAPHiCDiff_background_build
    EXIT_CODE: 10001
    EXIT_MESSAGE: Job command execution failed, errCode: 1

[ERROR] Chromosome Ghir_A01 not found in /share/.../L282.mcool::/resolutions/4000
exit 1

# 所有 36 个样本均以相同错误失败
Error in rule extract_and_balance_sample:
    message: For further error details see the cluster/cloud log and log files of the involved rule(s).
```

## Fix Applied

```bash
git checkout -b fix/extract-chrom-chrom-name-20260604_090736
git add config/chroms.yaml
git commit -m "fix: update chroms.yaml to HC04_* chromosome naming"

# 清理残留文件
rm -f tmp/*__Ghir_*
rm -f results/logs/extract_balance/*.log results/logs/extract_balance/*.done

# 重新提交
dsub -s run_snakemake.sh   # job 36931564
```

Pipeline 运行正常，36 个样本全部完成：

```
36931564.err: "Config file config/chroms.yaml is extended by additional config specified via the command line."
# output paths 全部使用 HC04_A01 ~ HC04_D13
# 111/286 steps (39%) done, 无报错
```

## Why the Bug Was Not Caught Earlier

1. **Working tree 未 commit 的变更无法被 dsub 运行捕获**：dsub 从 git HEAD 构建执行环境，不读取 working tree
2. **无染色体名验证**：Snakemake DAG 构建阶段不验证 `hictk dump` 能否匹配染色体名
3. **两次修改未原子 commit**：commit `c83f43b` 同时修改了 `chroms.yaml`（→ `Ghir_*`）和 `extract_chrom.smk`（→ 参数化），但两个文件都基于旧命名体系

## Prevention

| Priority | Mechanism | Specific Action | Status |
|----------|-----------|-----------------|--------|
| P0 | Pre-submission validation | 提交前 `grep -r "Ghir_" config/ rules/` 验证无残留旧命名 | TODO |
| P0 | CI/CD gate | dsub job 启动后立即验证染色体名匹配：`hictk dump {URI} --table chroms \| awk '{print $1}'` vs `chroms.yaml` | TODO |
| P1 | Configfile atomicity | 同一 commit 中修改 config 和对应 rule，确保两者同步 | TODO |
| P1 | 数据源清单 | 在 `docs/` 记录每批数据的染色体命名来源 | TODO |
| P2 | Manifest schema | 在 manifest TSV 中增加染色体命名列，供 pipeline 自动推断 | Future |

## Related Bugs

- 2026-05-20-sliding-window-checkpoint-missing-oe-dependency.md (checkpoint 依赖完整性)
- 2026-05-15-window-compare-path-hierarchy-misunderstanding.md (路径命名约定)

## Trellis

- **Symptom**: 36 个 extract_and_balance_sample 样本全部失败，相同的 "Chromosome not found" 错误
- **Root Cause**: committed `chroms.yaml` 使用 `Ghir_*`，working tree 使用 `HC04_*`；dsub 从 git HEAD 运行，读不到 working tree 变更
- **Evidence**: `.err` 中 `for CHR in Ghir_A01 Ghir_A02...`，`hictk dump` 实际输出 `HC04_A01`
- **Fix**: commit working tree `chroms.yaml`
- **DSub Validation**: job 36931564 正常运行，111/286 steps，无阻断性报错
- **Git Branch**: `fix/extract-chrom-chrom-name-20260604_090736` (尚未合并)
- **Reusable Pattern**: 任何影响所有样本的配置变更必须 commit 后才能通过 dsub 生效；working tree 变更对 dsub 不可见
