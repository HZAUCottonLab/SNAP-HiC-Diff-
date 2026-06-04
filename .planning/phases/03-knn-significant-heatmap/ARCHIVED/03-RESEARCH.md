# Phase 03: KNN 显著窗口热图绘制 - Research

**Researched:** 2026-05-30
**Domain:** Hi-C visualization + Snakemake dynamic wildcards
**Confidence:** HIGH (based on verified codebase patterns + established library APIs)

## Summary

Phase 03 adds a Snakemake rule to plot 2×2 heatmaps for windows passing 5-metric FDR < 0.05 significance. The top row shows iced O/E matrices (control left, treatment right); the bottom row shows Gaussian-smoothed matrices. Each significant window outputs one PNG + one PDF.

Three technical challenges dominate: (1) **FDR column availability** — the current `knn_enriched.tsv` has raw p-values but no pre-computed FDR columns, so the plotting script must either add FDR computation inline or the knn_empirical_pvalue script must be updated to pre-compute them; (2) **dynamic wildcards** — Snakemake needs to know the full set of significant windows (which varies per comparison) before constructing the DAG, requiring a checkpoint-based two-tier architecture; (3) **shared colorbars** — iced and gauss sub-pairs must each share a colorbar but with independent normalization.

The codebase already has all required libraries in `env/SNAP-HiC-Diff.yaml`: `cooler` for .cool reading, `h5py` for .h5 reading, `scipy` for BH-FDR, and `matplotlib` for plotting.

**Primary recommendation:** Use a checkpoint to scan `knn_enriched.tsv` and emit a manifest, then a regular rule with `{comparison}/{chrom}/{window_id}` wildcards plots one window per job. Compute FDR inline in the plotting script using `scipy.stats.false_discovery_control`.

---

## User Constraints (from CONTEXT.md)

### Locked Decisions

- **数据源 iced**: `results/balanced/{sample}/{sample}__{chrom}_{window_id}_iced.cool`
- **数据源 gauss**: `results/windows/{sample}/{sample}__{chrom}_{window_id}_gaussian.h5`
- **筛选来源**: `results/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv`
- **筛选条件**: 5 个 FDR 列（`*_empirical_FDR`）全部 < 0.05
- **布局**: 2 行 2 列（上行 iced / 下行 gauss；左列 control / 右列 treatment）
- **颜色条**: iced 共享一个，gauss 共享一个
- **坐标轴**: 基因组坐标（bp），每轴最多 5 个刻度
- **热图区域**: 固定窗口自身范围（`window_start` ~ `window_end`）
- **输出格式**: PNG + PDF
- **Snakemake rule**: `plot_knn_significant_heatmap`，wildcards: `{comparison}`, `{chrom}`, `{window_id}`
- **输出路径**: `results/knn_significant_heatmap/{comparison}/{chrom}_{window_id}.png/pdf`

### Deferred Ideas (OUT OF SCOPE)

- 聚合视图（多窗口一张图）— future phase
- HiGlass 集成 — future phase

---

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| FDR 筛选逻辑 | Python script | — | 纯计算逻辑，在 Snakemake job 内完成 |
| .cool 文件读取 | Python script (cooler) | — | 读取指定基因组区域的矩阵 |
| .h5 文件读取 | Python script (h5py) | — | 读取单个窗口的 gauss 矩阵 |
| 2×2 热图绘制 | Python script (matplotlib) | — | 纯可视化逻辑 |
| Snakemake 调度 | Snakemake engine | — | 管理动态 wildcards 和 job 并行化 |
| checkpoint manifest | Snakemake checkpoint | — | 动态发现显著窗口数量 |

---

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| cooler | ≥0.9.0 | 读取 .cool 文件的基因组区域矩阵 | 已在 `compute_oe_matrix.py` 验证 |
| h5py | ≥3.10 | 读取 `*_gaussian.h5` 窗口矩阵 | 已在 `sliding_window.py` 验证 |
| scipy | ≥1.10 | BH-FDR 校正 (`false_discovery_control`) | 已在 `analyze_fdr_pvalue.py` 验证 |
| matplotlib | ≥3.7 | 2×2 热图 + 共享颜色条 | 已在 `knn_viz.py` 验证 |
| numpy | ≥1.24 | 矩阵操作、NaN 处理 | 已在多个脚本验证 |
| pandas | ≥2.0 | TSV 读取、窗口筛选 | 已在 `knn_viz.py` 验证 |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `cooler.Cooler().matrix().fetch()` | `cooltools.api.expected` fetch | `cooler` 是更底层、更稳定的 API，项目已用 |
| `plt.colorbar(im, ax=axs)` | `make_axes_locatable` (axes_grid1) | 内置 API 更简洁，不需要额外 import |
| inline FDR 计算 | 改 `compute_knn_empirical_pvalue.py` 输出 FDR 列 | inline 更轻量，不改已有 pipeline 输出格式 |

**Installation:** 无需额外安装 — 所有库已存在于 `env/SNAP-HiC-Diff.yaml`。

---

## Architecture Patterns

### System Architecture Diagram

```
knn_enriched.tsv (input)
         │
         ▼
┌─────────────────────────┐
│  checkpoint:            │
│  scan knn_enriched.tsv  │  ──→ 筛选 FDR < 0.05 的窗口
│  emit manifest           │  ──→ manifest = {comparison}_significant_windows.tsv
└─────────────────────────┘
         │
         ▼  (Snakemake DAG expands based on manifest)
┌──────────────────────────────────────────────────────────┐
│  rule plot_knn_significant_heatmap (per window job):   │
│                                                          │
│  input:                                                 │
│    - knn_enriched.tsv  (提供窗口坐标)                   │
│    - {sample_a}__{chrom}_{window_id}_iced.cool         │
│    - {sample_b}__{chrom}_{window_id}_gaussian.h5       │
│                                                          │
│  process:                                               │
│    1. 读取 iced .cool → 矩阵（cooler）                 │
│    2. 读取 gauss .h5  → 矩阵（h5py）                   │
│    3. 计算两行各自 vmin/vmax                           │
│    4. 绘制 2×2 subplot + 2 个共享颜色条                │
│    5. 保存 PNG + PDF                                   │
│                                                          │
│  output:                                                │
│    - results/knn_significant_heatmap/{comparison}/      │
│      {chrom}_{window_id}.png                           │
│    - results/knn_significant_heatmap/{comparison}/      │
│      {chrom}_{window_id}.pdf                           │
└──────────────────────────────────────────────────────────┘
```

### Recommended Project Structure

```
scripts/
└── plot_knn_significant_heatmap.py   # 主绘图脚本

rules/
└── knn_significant_heatmap.smk       # Snakemake rule

results/
└── knn_significant_heatmap/
    └── {comparison}/
        ├── {chrom}_{window_id}.png
        └── {chrom}_{window_id}.pdf
```

### Pattern 1: Two-Tier Snakemake Architecture for Dynamic Wildcards

**What:** Checkpoint discovers dynamic file list; regular rule declares files with wildcards.
**When to use:** Number of output files is data-dependent and unknown at Snakefile parse time.
**Reference:** `rules/sliding_window.smk` — `checkpoint sliding_window` + `rule sliding_window_h5`

```python
# Tier 1: checkpoint — scans knn_enriched.tsv, writes manifest
checkpoint plot_knn_significant_heatmap:
    input:
        script = SCRIPT_DIR + "/scripts/plot_knn_significant_heatmap.py",
        enriched = OUTDIR + "/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
    output:
        manifest = OUTDIR + "/logs/knn_significant_heatmap/{comparison}.done"
    shell:
        r"""
        python {input.script} --scan-only \
            --enriched "{input.enriched}" \
            --out "{output.manifest}" \
            > "{log}" 2>&1
        """

# Tier 2: regular rule — one job per significant window
rule plot_knn_significant_window:
    input:
        lambda w: checkpoints.plot_knn_significant_heatmap.get(
            comparison=w.comparison
        ).output.manifest,
        # ...cool and h5 files...
    output:
        png = OUTDIR + "/knn_significant_heatmap/{comparison}/{chrom}_{window_id}.png",
        pdf = OUTDIR + "/knn_significant_heatmap/{comparison}/{chrom}_{window_id}.pdf",
    shell:
        r"""
        python {input.script} \
            --comparison "$COMPARISON" \
            --chrom "$CHR" \
            --window-id "$WID" \
            --out-dir "$OUTDIR" \
            > "{log}" 2>&1
        """
```

**Anti-pattern:** Using `expand()` with dynamic lists in a regular rule's `output:` — the expand must happen at DAG construction time, but the list is unknown until checkpoint runs.

### Pattern 2: Shared Colorbar Per Row

**What:** 2×2 subplot where top row shares one colorbar, bottom row shares another.
**Reference:** `matplotlib.org` multi_image example — `colors.Normalize(vmin, vmax)` per row.

```python
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

def plot_2x2_heatmap(m_iced_ctrl, m_iced_treat, m_gauss_ctrl, m_gauss_treat,
                      win_start, win_end, chrom, out_png, out_pdf):
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))

    # ── Row 1: iced (shared colorbar) ──────────────────────────────────
    vmin_iced = min(np.nanmin(m_iced_ctrl), np.nanmin(m_iced_treat))
    vmax_iced = max(np.nanmax(m_iced_ctrl), np.nanmax(m_iced_treat))
    norm_iced = mcolors.Normalize(vmin=vmin_iced, vmax=vmax_iced)

    im_iced_ctrl = axes[0, 0].imshow(
        m_iced_ctrl, norm=norm_iced, cmap="RdBu_r",
        extent=[win_start, win_end, win_end, win_start], origin="upper"
    )
    axes[0, 1].imshow(
        m_iced_treat, norm=norm_iced, cmap="RdBu_r",
        extent=[win_start, win_end, win_end, win_start], origin="upper"
    )

    # ── Row 2: gauss (shared colorbar) ──────────────────────────────────
    vmin_gauss = min(np.nanmin(m_gauss_ctrl), np.nanmin(m_gauss_treat))
    vmax_gauss = max(np.nanmax(m_gauss_ctrl), np.nanmax(m_gauss_treat))
    norm_gauss = mcolors.Normalize(vmin=vmin_gauss, vmax=vmax_gauss)

    im_gauss_ctrl = axes[1, 0].imshow(
        m_gauss_ctrl, norm=norm_gauss, cmap="RdBu_r",
        extent=[win_start, win_end, win_end, win_start], origin="upper"
    )
    axes[1, 1].imshow(
        m_gauss_treat, norm=norm_gauss, cmap="RdBu_r",
        extent=[win_start, win_end, win_end, win_start], origin="upper"
    )

    # ── Axis formatting ──────────────────────────────────────────────────
    for ax in axes.flat:
        ax.set_xlabel(f"{chrom} position (bp)")
        ax.set_ylabel(f"{chrom} position (bp)")
        # 每轴最多 5 个刻度
        ax.locator_params(axis='both', nbins=5)

    # ── Colorbars ────────────────────────────────────────────────────────
    fig.colorbar(im_iced_ctrl, ax=axes[0, :], orientation='vertical',
                 fraction=0.02, pad=0.04, label='O/E (iced)')
    fig.colorbar(im_gauss_ctrl, ax=axes[1, :], orientation='vertical',
                 fraction=0.02, pad=0.04, label='z-score (gauss)')

    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
```

### Pattern 3: Reading .cool Region Matrix

**Reference:** `scripts/compute_oe_matrix.py` lines 73-76.

```python
import cooler

def load_iced_matrix(cool_path: str, chrom: str, window_start: int, window_end: int) -> np.ndarray:
    """Load ICE O/E matrix for a genomic region from a .cool file."""
    clr = cooler.Cooler(cool_path)
    # 验证染色体名
    chroms_in_file = list(clr.chromnames)
    actual_chrom = chroms_in_file[0] if chrom in chroms_in_file else chrom
    # fetch 返回密集 numpy 数组
    mat = clr.matrix(balance=True, sparse=True).fetch(actual_chrom).toarray()
    return mat
```

### Pattern 4: Reading Gaussian h5 Matrix

**Reference:** `scripts/compare_gaussian_smoothed_windows.py` lines 38-43 and `sliding_window.py` lines 335-347.

```python
import h5py

def load_gaussian_matrix(h5_path: str, window_id: int) -> np.ndarray:
    """Load Gaussian-smoothed matrix from per-window HDF5 file."""
    with h5py.File(h5_path, 'r') as f:
        key = f"w{window_id:06d}"
        matrix = f[key]['matrix'][:]
    return matrix
```

H5 file structure (per `sliding_window.py`):
```
{window_id:06d}/
  matrix      # shape [n_bins, n_bins], float32
  window_start # scalar, genomic start (bp)
  window_end   # scalar, genomic end (bp)
```

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Hi-C matrix region fetch | 自己解析 .cool 二进制格式 | `cooler.Cooler().matrix().fetch()` | .cool 是 HDF5 格式族，有成熟 API |
| Gaussian smooth | 自己写高斯卷积 | `scipy.ndimage.gaussian_filter` | 已有 `sliding_window.py` 验证 |
| FDR 校正 | 手写 BH 算法 | `scipy.stats.false_discovery_control` | 已在 `analyze_fdr_pvalue.py` 验证 |
| 共享颜色条 | 手动计算 colorbar 位置 | `colors.Normalize(vmin, vmax)` + `fig.colorbar(im, ax=...)` | matplotlib 内置支持 |

**Key insight:** Hi-C 数据格式（HDF5 变体）和统计校正是高度标准化的领域，主流工具有成熟实现。

---

## Common Pitfalls

### Pitfall 1: FDR 列不存在导致筛选失败
**What goes wrong:** 脚本期望 `*_empirical_FDR` 列但实际 `knn_enriched.tsv` 只有原始 p 值列。
**Why it happens:** `compute_knn_empirical_pvalue.py` 默认 `compute_fdr: false`（config.yaml line 94），且 03-CONTEXT.md 声称有 FDR 列但实际文件没有。
**How to avoid:** 在绘图脚本内用 `scipy.stats.false_discovery_control` 实时计算 FDR，不依赖文件预存列。
**Warning signs:** `KeyError: '*_empirical_FDR'` 或 0 个显著窗口。

### Pitfall 2: Snakemake checkpoint 输出文件作为下游 input
**What goes wrong:** checkpoint 的 `output.manifest` 被标记为 input 后，manifest 本身（而非真实输出文件）的 mtime 触发下游重跑。
**Why it happens:** Snakemake 用 manifest 的 mtime 决定下游 job 是否需要重跑，但 manifest 只是 touch 文件。
**How to avoid:** 使用标准的 checkpoint → regular rule 二级架构（见 Pattern 1），不把 manifest 声明为最终 output。
**Warning signs:** 删除 PNG 后 Snakemake 仍认为 job 最新（因为 manifest 没变）。

### Pitfall 3: 颜色条数值范围不一致
**What goes wrong:** iced 和 gauss 行各自颜色条数值范围差异巨大（iced 是 log O/E ~[-3, 3]，gauss 是 z-score ~[-5, 5]），共用 vmin/vmax 导致对比困难。
**Why it happens:** 两个矩阵的数值分布完全不同，但代码可能用全局 vmin/vmax。
**How to avoid:** 分别对 iced 行和 gauss 行计算独立的 `Normalize(vmin, vmax)`（见 Pattern 2 代码示例）。
**Warning signs:** 热图全是单一颜色，或颜色条显示 [0, 1] 默认范围。

### Pitfall 4: 窗口矩阵与染色体范围不匹配
**What goes wrong:** 从 .cool fetch 的区域包含额外的 bin（来自 cooler 的 bin 对齐），与 Gaussian h5 的矩阵 size 不一致。
**Why it happens:** `.cool` 文件以 bin 为单位存储，`fetch(chrom:start-end)` 可能返回包含额外 bin 的矩阵。
**How to avoid:** 使用 `window_start`/`window_end` 计算目标 bin 范围，然后用 `matrix[bin_start:bin_end, bin_start:bin_end]` 切片匹配 h5 矩阵大小。
**Warning signs:** `ValueError: shapes not aligned` 或图像变形。

### Pitfall 5: Hi-C 矩阵不对称（NaN 对角线）
**What goes wrong:** 矩阵在可视化时显示为一半有数据一半为空。
**Why it happens:** Hi-C 是对称矩阵，但 `sliding_window.py` 在 ln_oe 阶段将 d ≤ 1 对角线设为 NaN（line 238-240）。
**How to avoid:** 绘制时用 `np.triu` + `np.tril` 填充对称，或使用 `mat = (mat + mat.T) / 2` 强制对称化。
**Warning signs:** 图像左下角或右上角全是白色/空白。

---

## Code Examples

### FDR Computation (Inline)

```python
# Source: scripts/analyze_fdr_pvalue.py (verified)
from scipy.stats import false_discovery_control
import numpy as np

PVALUE_COLS = [
    "MSE_empirical_pvalue",
    "d_SSIM_empirical_pvalue",
    "d_Spearman_empirical_pvalue",
    "tw_score_empirical_pvalue",
    "d_NCC_empirical_pvalue",
]

def filter_significant_windows(df, fdr_threshold=0.05):
    """Return windows where ALL 5 metrics pass FDR < threshold."""
    pval_matrix = df[PVALUE_COLS].values.astype(float)
    fdr_matrix = false_discovery_control(pval_matrix, method='bh')
    # 所有 5 个指标 FDR < threshold
    all_pass = np.all(fdr_matrix < fdr_threshold, axis=1)
    return df[all_pass].copy()
```

### Snakemake Wildcard Extraction

```python
# Source: rules/window_compare.smk (verified pattern)
shell:
    r"""
    set -euo pipefail

    COMPARISON="{wildcards.comparison}"
    CHR="{wildcards.chrom}"
    WID="{wildcards.window_id}"

    # shell 中使用 ${{VAR}} 避免与 wildcard 语法冲突
    BASE="${{CHR}}_${{WID}}"
    """
```

### matplotlib 2×2 GridSpec with Two Colorbars

```python
# Source: matplotlib.org multi_image example + StackOverflow #46040631
import matplotlib.pyplot as plt
from matplotlib import colors

fig, axes = plt.subplots(2, 2, figsize=(10, 10))

# Row 1: iced — shared norm
vmin1 = min(np.nanmin(ax1_data), np.nanmin(ax2_data))
vmax1 = max(np.nanmax(ax1_data), np.nanmax(ax2_data))
norm1 = colors.Normalize(vmin=vmin1, vmax=vmax1)
axes[0, 0].imshow(ax1_data, norm=norm1, cmap="RdBu_r")
axes[0, 1].imshow(ax2_data, norm=norm1, cmap="RdBu_r")

# Row 2: gauss — shared norm
vmin2 = min(np.nanmin(ax3_data), np.nanmin(ax4_data))
vmax2 = max(np.nanmax(ax3_data), np.nanmax(ax4_data))
norm2 = colors.Normalize(vmin=vmin2, vmax=vmax2)
axes[1, 0].imshow(ax3_data, norm=norm2, cmap="RdBu_r")
axes[1, 1].imshow(ax4_data, norm=norm2, cmap="RdBu_r")

# Two separate colorbars
fig.colorbar(axes[0, 0].images[0], ax=axes[0, :], location='right',
             fraction=0.02, pad=0.04)
fig.colorbar(axes[1, 0].images[0], ax=axes[1, :], location='right',
             fraction=0.02, pad=0.04)
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| 手动写 .cool 解析 | `cooler.Cooler().matrix().fetch()` | 2018+ | 标准化，数据格式变化时只需升级 cooler |
| 固定 vmin/vmax | `colors.Normalize(vmin, vmax)` per row | matplotlib 3.x | 每个子图自动缩放，颜色对比更准确 |
| f-string 路径 | `.format()` + `{{}}` 转义 | Snakemake 8.x | 避免与 wildcard 语法冲突 |
| checkpoint 直接输出真实文件 | 二级架构（manifest + regular rule） | 已知 pattern | mtime 追踪更可靠 |

**Deprecated/outdated:**
- `plt.colormaps["RdBu_r"]` → 直接用 `cmap="RdBu_r"` 字符串（无变化，just cleaner）
- Snakemake 7.x 动态 wildcards → Snakemake 8.x 改进但 checkpoint 模式兼容

---

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | `knn_enriched.tsv` 没有预存 FDR 列，需 inline 计算 | Common Pitfalls | 中 — 不影响实现，只影响是否需要改 knn_empirical_pvalue 脚本 |
| A2 | iced .cool 文件中染色体名与 chrom 列匹配 | Pattern 3 | 低 — 可用 `clr.chromnames[0]` 兜底 |
| A3 | Gaussian h5 矩阵与窗口坐标一致 | Pattern 4 | 低 — sliding_window.py 保证 |
| A4 | 每个显著窗口的 cool/h5 文件都存在 | Architecture | 低 — 如果文件缺失会报错，不产生假阳性 |

---

## Open Questions

1. **FDR 列来源** ✅ (RESOLVED)
   - `knn_enriched.tsv` 无 FDR 列，`compute_fdr: false` 已在 config.yaml 中确认
   - 绘图脚本内 inline BH 计算：见 Pattern 1（`false_discovery_control`）
   - 或者改 `compute_knn_empirical_pvalue.py` 将 `compute_fdr: false` 改为 `true`
   - **已选方案**：checkpoint 脚本内 inline BH，不改已有 pipeline 输出格式

2. **显著窗口数量估计** ✅ (RESOLVED)
   - `H1398_vs_L134` 有 57,872 个总窗口
   - 假设 ~100-500 个通过 5 指标 FDR < 0.05（最终以实际执行为准）
   - Snakemake checkpoint 在运行时发现，不在规划阶段预设

3. **iced 矩阵来源** ✅ (RESOLVED)
   - **per-chrom cool 文件，不存在 per-window cool 文件**
   - 路径：`results/balanced/{sample}/{sample}__{chrom}_iced.cool`（如 `L134__HC04_A01_iced.cool`）
   - 用 `cooler.Cooler().matrix().fetch(chrom, window_start, window_end)` 提取窗口区域
   - ROADMAP.md criterion 已更新为此路径

---

## Environment Availability

Step 2.6: SKIPPED — no external dependencies beyond project conda environment.

All required tools are verified present in `env/SNAP-HiC-Diff.yaml`:
- `cooler>=0.9.0` ✓
- `cooltools>=0.5.4` ✓
- `h5py>=3.10` ✓
- `scipy>=1.10` ✓
- `matplotlib>=3.7` ✓
- `numpy>=1.24` ✓
- `pandas>=2.0` ✓

No additional installation needed.

---

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest + inline Python assertions |
| Config file | none — tests run via dsub |
| Quick run command | `dsub -q arm -aa -n test -e %J.err -o %J.out -R "cpu=10" python -c "from scripts.plot_knn_significant_heatmap import *; print('import ok')"` |
| Full suite command | `dsub -s scripts/test_knn_significant_heatmap.sh` |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| REQ-03-1 | FDR < 0.05 筛选正确（所有 5 指标） | unit | `pytest tests/test_fdr_filter.py::test_all_5_fdr_pass -x` | ❌ Wave 0 |
| REQ-03-2 | 从 .cool 正确提取窗口区域矩阵 | unit | `pytest tests/test_cool_fetch.py::test_window_region -x` | ❌ Wave 0 |
| REQ-03-3 | 从 .h5 正确读取 gaussian 矩阵 | unit | `pytest tests/test_h5_read.py::test_window_matrix -x` | ❌ Wave 0 |
| REQ-03-4 | 2×2 热图输出 PNG + PDF | integration | `dsub ... python scripts/plot_knn_significant_heatmap.py --comparison H1398_vs_L134 --chrom HC04_A01 --window-id 26` | ❌ Wave 0 |
| REQ-03-5 | Snakemake DAG 正确展开（checkpoint + jobs） | integration | `snakemake -n results/knn_significant_heatmap/H1398_vs_L134/HC04_A01_26.png` | ❌ Wave 0 |

### Sampling Rate
- **Per task commit:** Quick unit test via dsub
- **Per wave merge:** Full suite via dsub script
- **Phase gate:** Full suite green before `/gsd-verify-work`

### Wave 0 Gaps
- [ ] `tests/test_fdr_filter.py` — covers REQ-03-1
- [ ] `tests/test_cool_fetch.py` — covers REQ-03-2
- [ ] `tests/test_h5_read.py` — covers REQ-03-3
- [ ] `scripts/test_knn_significant_heatmap.sh` — dsub submit script for REQ-03-4, REQ-03-5
- [ ] Framework install: already in `env/SNAP-HiC-Diff.yaml` — no action needed

*(Project has no existing `tests/` directory — Wave 0 creates it.)*

---

## Security Domain

N/A — this phase is visualization-only with no user input, network access, or sensitive data paths. No authentication, authorization, or cryptographic operations are involved.

---

## Sources

### Primary (HIGH confidence)
- `scripts/compute_oe_matrix.py` — cooler.Cooler usage verified
- `scripts/sliding_window.py` — h5 structure verified (lines 335-347)
- `scripts/compare_gaussian_smoothed_windows.py` — h5 reading verified (lines 38-43)
- `scripts/analyze_fdr_pvalue.py` — scipy.stats.false_discovery_control usage verified
- `jupyter/knn_viz.py` — matplotlib plotting patterns verified
- `rules/sliding_window.smk` — checkpoint architecture verified
- `rules/window_compare.smk` — wildcard + shell pattern verified
- `env/SNAP-HiC-Diff.yaml` — library availability verified

### Secondary (MEDIUM confidence)
- `matplotlib.org/stable/gallery/images_contours_and_fields/multi_image.html` — shared colorbar pattern
- `stackoverflow.com/questions/46040631` — two colorbars per grid pattern
- `compute_oe_matrix.smk` — iced cool file path pattern verified

### Tertiary (LOW confidence)
- FDR column naming assumption — based on 03-CONTEXT.md, not verified against actual file

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all libraries verified in project codebase
- Architecture: HIGH — checkpoint + regular rule pattern verified in `sliding_window.smk`
- Pitfalls: MEDIUM — FDR column discrepancy identified from file inspection

**Research date:** 2026-05-30
**Valid until:** 2026-06-30 (stable domain — Snakemake + matplotlib patterns well-established)
