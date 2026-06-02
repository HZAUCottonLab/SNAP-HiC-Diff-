# Bug Analysis: 2D Gap Interpolation — 1D/2D Mask 语义混淆

**Date**: 2026-05-11
**Severity**: MEDIUM (wrong biological values in output)
**Status**: FIXED

## Problem

用 2D 最近邻插值替换 Stage 2 的 naive 0-fill 时，单元测试在测试框架内通过，但在真实数据上失败。

表现：
- 无 gap 时：部分矩阵位置被错误地设为 0（预期保留原值）
- 有 gap 时：插值源来自 NaN 区域（预期来自有效值）

## Root Cause

**类别**: E — 隐式假设（`gap_bin_mask` 维度语义未明确文档化）

### 核心语义混淆

`gap_bin_mask` 是 **1D 布尔数组**，但 2D 矩阵中 gap 的含义是：

```
gap_bin_mask[k] = True
        ↓
row k 全为 gap 区域（整个第 k 行）
AND
col k 全为 gap 区域（整个第 k 列）
```

即：2D 矩阵中 `(i, j)` 是 gap 当且仅当 `gap_bin_mask[i] == True OR gap_bin_mask[j] == True`。

**错误假设**：将 `gap_bin_mask` 直接当作 2D mask 使用——`gap_bin_mask[i, j]` 不存在（1D），但代码逻辑隐式地将它当作 2D 使用。

### 具体错误

第一次实现 `_classify_bins` 时：

```python
# 错误：gap_bin_mask 是 1D，却在循环中当 2D 用
for i in range(n):
    row_gaps = gap_bin_mask.copy()   # ← 整个 1D 数组复制到每一行
    j = 0
    while j < n:
        if row_gaps[j]:              # ← row_gaps[j] 只检查列索引 j 是否 gap
            ...
```

实际 `gap_bin_mask[j]` 只描述"列 j 本身是否是 gap 列"，但对角线元素 `(i, j)` 的 gap 状态应该是 `gap_bin_mask[i] or gap_bin_mask[j]`。

### 第二个错误：RGI 的边界行为

第一次用 `scipy.interpolate.RegularGridInterpolator` 时，gap row 的值是 NaN → 被设为 0 作为占位符 → RGI 的 `method='nearest'` 将 0 识别为最近邻（最近的有效像素本身就是 0）。

## Why Fixes Failed

1. **Surface Fix #1**（替换 `gap_bin_mask` 为 2D mask）：仍然用了错误的循环逻辑，`row_gaps = gap_bin_mask.copy()` 仍然每次都是整个 1D 数组。仍未理解 gap_bin_mask 的 1D 语义。

2. **Surface Fix #2**（跳过 gap row 不 zeroing）：gap row 中仍有非 gap 列的元素需要插值，但跳过 zeroing 后，gap row 的 NaN 值被 RGI 识别为最近邻（因为 RGI nearest 不区分 NaN 和 0）。

3. **Correct Fix**：彻底重构为两步：
   - 第一步：`_build_2d_gap_masks(gap_bin_mask, n)` — 用 run-length 分析从 1D gap_bin_mask 构造 2D masks
   - 第二步：`_interpolate_2d_nearest` — 显式最近邻，搜索 `fill_input != 0.0` 的位置作为源，不依赖 RGI 的插值行为

## Prevention

| 优先级 | 机制 | 具体动作 | 状态 |
|--------|------|----------|------|
| P0 | 文档化数据结构 | 在 `sliding_window.py` 注释中明确标注 `gap_bin_mask` 的维度语义 | DONE |
| P0 | 类型注解 | 给 `gap_bin_mask` 添加 numpy 类型注解 `np.ndarray[n]` vs `np.ndarray[n,n]` | TODO |
| P1 | 集成测试 | 补充真实数据上的冒烟测试（不仅仅是单元测试） | DONE |
| P1 | 维度检查断言 | 在 helper 函数入口添加 shape assertion | TODO |

## Systematic Expansion

### 类似问题

- 1D/2D 混淆在 Hi-C 矩阵处理中非常普遍（bin index vs row/col index）
- 其他可能存在同样问题的地方：
  - `distance_stratified_zscore` 中的对角线处理（`|i-j|` 是 1D 距离 vs 2D 矩阵位置）
  - `compute_oe_matrix.smk` 中的 expected value 索引（bin distance vs 矩阵索引）

### 架构改进

- 所有处理 Hi-C 矩阵的函数，应在 docstring 中明确标注：
  - 输入 `gap_bin_mask` 的维度（1D bin-level）
  - 在 2D 矩阵中 `gap_bin_mask[i]` 的确切含义（row i AND col i 都是 gap）

## Knowledge Capture

- [x] 写入 bug-analysis 文档（本文）
- [x] 更新 `sliding_window.py` 中 `gap_bin_mask` 的 docstring，标注维度语义
- [x] 在 `cross-layer-thinking-guide.md` 中增加 Hi-C 矩阵 1D/2D 语义 case study
- [ ] 在 `.trellis/spec/guides/` 中补充 Hi-C 矩阵维度语义独立 checklist
