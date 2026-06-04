# Phase 03: KNN 显著窗口热图绘制 — Validation Strategy

**Phase:** 03-knn-significant-heatmap
**Date:** 2026-05-30
**Status:** Ready for planning

---

## Validation Dimensions

### Dimension 1 — Functional Correctness

| Criterion | Method | Pass Condition |
|-----------|--------|----------------|
| FDR < 0.05 筛选逻辑 | pytest unit test | 5 个指标 FDR 全部 < 0.05 的窗口被选出；至少 1 个指标 FDR ≥ 0.05 的窗口被过滤 |
| cool 文件读取 | pytest unit test | 返回非空 numpy 数组，shape 与窗口 bin 数一致 |
| h5 文件读取 | pytest unit test | 返回非空 numpy 数组，shape 为 [n_bins, n_bins] |
| 2x2 热图绘制 | 集成测试 | 生成的 PNG 和 PDF 文件存在且非空（> 10KB） |

### Dimension 2 — Integration

| Criterion | Method | Pass Condition |
|-----------|--------|----------------|
| Snakemake DAG | `snakemake -n` dry-run | 无语法错误，checkpoint 展开后 rule 数量 = 显著窗口数量 |
| PNG+PDF 输出 | glob 检查 | `results/knn_significant_heatmap/{comparison}/` 下每个显著窗口有 .png 和 .pdf |

### Dimension 3 — Performance

| Criterion | Method | Pass Condition |
|-----------|--------|----------------|
| 单窗口绘图时间 | benchmark log | < 60s per window |
| 并行 job 数 | config | threads_per_job=80 |

### Dimension 4 — Reproducibility

| Criterion | Method | Pass Condition |
|-----------|--------|----------------|
| 环境一致性 | conda env export diff | `env/SNAP-HiC-Diff.yaml` 无新增依赖 |
| 随机种子 | 无随机操作 | N/A |

### Dimension 5 — Security

| Criterion | Method | Pass Condition |
|-----------|--------|----------------|
| 路径遍历 | 代码审查 | 无用户可控路径拼接 |
| 输入验证 | 代码审查 | 窗口坐标在合理范围内 |

### Dimension 6 — Observability

| Criterion | Method | Pass Condition |
|-----------|--------|----------------|
| 日志输出 | 检查 log 文件 | 包含 significant window count、cool/h5 路径、FDR threshold |
| benchmark | benchmark TSV | 每 job 有 walltime、peak_rss 记录 |

### Dimension 7 — Fault Tolerance

| Criterion | Method | Pass Condition |
|-----------|--------|----------------|
| 缺失文件处理 | 代码审查 + 模拟测试 | cool/h5 文件不存在时报错，不静默失败 |
| NaN 处理 | 代码审查 | 热图矩阵含 NaN 时用 `np.nanmin/nanmax` 计算颜色范围 |

### Dimension 8 — Automated Verification

| Criterion | Method | Pass Condition |
|-----------|--------|----------------|
| 单元测试 | pytest | `tests/` 下所有测试通过 |
| Snakemake dry-run | `snakemake -n` | 0 个预期外规则 |
| 输出文件存在 | glob + file size | 每个显著窗口有 PNG + PDF，> 5KB |

### Dimension 9 — Usability

| Criterion | Method | Pass Condition |
|-----------|--------|----------------|
| 错误信息 | 手动测试缺失文件 | 报错信息包含文件名和缺失原因 |
| 使用文档 | 代码审查 | 脚本 `--help` 输出完整参数说明 |

---

## Test Plan

### Wave 0 (Pre-implementation — optional)

```bash
# FDR 筛选逻辑
pytest tests/test_fdr_filter.py -x

# cool 文件读取
pytest tests/test_cool_fetch.py -x

# h5 文件读取
pytest tests/test_h5_read.py -x
```

### Wave 1 (Post-implementation)

```bash
# Snakemake dry-run
snakemake -n results/knn_significant_heatmap/H1398_vs_L134.done

# 实际运行（1 个窗口）
snakemake results/knn_significant_heatmap/H1398_vs_L134/HC04_A01_26.png

# 检查输出
python -c "
import glob, os
pngs = glob.glob('results/knn_significant_heatmap/H1398_vs_L134/*.png')
pdfs = glob.glob('results/knn_significant_heatmap/H1398_vs_L134/*.pdf')
print(f'PNG: {len(pngs)}, PDF: {len(pdfs)}')
sizes = [os.path.getsize(f) for f in pngs]
print(f'Min PNG size: {min(sizes)/1024:.1f} KB')
"
```

---

## Coverage Summary

| Dimension | Covered By |
|-----------|------------|
| Functional correctness | pytest unit tests + integration |
| Integration | snakemake -n + dry-run |
| Performance | benchmark TSV |
| Security | 代码审查 |
| Observability | log 文件检查 |
| Fault tolerance | 代码审查 |
| Automated verification | pytest + snakemake -n |
| Usability | --help + 错误测试 |
