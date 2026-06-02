# Bug Analysis: sliding_window checkpoint 竞态条件

**Date**: 2026-05-20
**Severity**: HIGH
**Status**: FIXED

## Summary

`sidling_window` checkpoint 因缺少 `oe_files` 依赖声明，导致并行执行时出现 FileNotFoundError。

## Root Cause

checkpoint 只声明 `gaps` 依赖，未声明 `oe_files`：

```python
# 修复前
checkpoint sliding_window:
    input:
        gaps = expand(..., chrom=CHROMS),
        script = SCRIPT_DIR + "/scripts/sliding_window.py"  # 缺少 oe_files!
```

当 `compute_oe_matrix` 和 `sliding_window` 并行运行时，OE 文件可能在 sliding_window 启动后才创建。

## Symptom

```
FileNotFoundError: [Errno 2] No such file or directory: 'results/oe/H1398_rep2__HC04_A01_oe.npz'
```

时间戳分析：
- `sliding_window` 启动：16:32:48
- OE 文件创建：16:33:40
- **竞态**：~52 秒延迟

## Fix

```python
# 修复后
checkpoint sliding_window:
    input:
        gaps = expand(..., chrom=CHROMS),
        oe_files = expand(OUTDIR + "/oe/{{sample}}__{chrom}_oe.npz", chrom=CHROMS),
        script = SCRIPT_DIR + "/scripts/sliding_window.py"
```

## Prevention

1. **checkpoint/input 必须显式声明所有上游依赖**
2. **修改任何规则后，检查下游 checkpoint 是否需要更新**
3. **并行执行前验证 DAG 依赖完整性**

## Related Bugs

- 2026-05-20-checkpoint-real-file-mtime-tracking.md (相关 checkpoint 架构问题)
