#!/usr/bin/env python3
"""分析 KNN enriched 文件中各 P 值指标 FDR 校正后通过 0.05 阈值的窗口数量."""

import sys
import pandas as pd
import numpy as np
from scipy.stats import false_discovery_control

INFILE = "/share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff/results/knn_empirical_pvalue/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv"
THRESHOLD = 0.05

PVALUE_COLS = [
    "MSE_empirical_pvalue",
    "d_SSIM_empirical_pvalue",
    "d_Spearman_empirical_pvalue",
    "tw_score_empirical_pvalue",
    "d_NCC_empirical_pvalue",
]

def main():
    print(f"Loading: {INFILE}")
    df = pd.read_csv(INFILE, sep="\t")
    total_windows = len(df)
    print(f"Total windows: {total_windows}\n")

    # 1. 原始 P 值 < 0.05 的统计
    print("=" * 60)
    print("一、原始 P 值 < 0.05 的统计（未校正）")
    print("=" * 60)
    raw_pass = {}
    for col in PVALUE_COLS:
        vals = df[col].values
        n_pass = np.sum(vals < THRESHOLD)
        raw_pass[col] = n_pass
        pct = n_pass / total_windows * 100
        print(f"  {col:35s}: {n_pass:6d} / {total_windows} ({pct:5.2f}%)")

    # 2. Benjamini-Hochberg FDR 校正
    print("\n" + "=" * 60)
    print("二、BH-FDR 校正后的统计")
    print("=" * 60)

    fdr_pass = {}
    for col in PVALUE_COLS:
        raw_pvals = df[col].values.astype(float)

        # scipy.stats.false_discovery_control with method='bh'
        fdr_pvals = false_discovery_control(raw_pvals, method='bh')

        n_pass = np.sum(fdr_pvals < THRESHOLD)
        fdr_pass[col] = (n_pass, fdr_pvals)
        pct = n_pass / total_windows * 100
        print(f"  {col:35s}: {n_pass:6d} / {total_windows} ({pct:5.2f}%)")

    # 3. 任意一个指标 FDR < 0.05 的窗口数
    print("\n" + "=" * 60)
    print("三、任意一个指标 FDR < 0.05 的窗口数（去重）")
    print("=" * 60)

    # 重新计算一次，用矩阵方式
    pval_matrix = df[PVALUE_COLS].values.astype(float)
    fdr_matrix = false_discovery_control(pval_matrix, method='bh')

    any_pass_mask = np.any(fdr_matrix < THRESHOLD, axis=1)
    n_any_pass = np.sum(any_pass_mask)
    pct_any = n_any_pass / total_windows * 100
    print(f"  任意一个指标 FDR < {THRESHOLD} 的窗口: {n_any_pass:6d} / {total_windows} ({pct_any:5.2f}%)")

    # 各指标分别通过的交集情况
    print("\n" + "=" * 60)
    print("四、各指标通过情况的交集分析")
    print("=" * 60)

    # 统计只通过1个、2个、...、5个指标的窗口数
    n_pass_per_window = np.sum(fdr_matrix < THRESHOLD, axis=1)
    for k in range(1, len(PVALUE_COLS) + 1):
        n_windows = np.sum(n_pass_per_window == k)
        if n_windows > 0:
            pct = n_windows / total_windows * 100
            print(f"  通过 {k} 个指标 FDR < {THRESHOLD}: {n_windows:6d} / {total_windows} ({pct:5.2f}%)")

    n_zero = np.sum(n_pass_per_window == 0)
    print(f"  通过 0 个指标 FDR < {THRESHOLD}: {n_zero:6d} / {total_windows} ({n_zero/total_windows*100:5.2f}%)")

    # 5. 每个指标 FDR < 0.05 的详细分布
    print("\n" + "=" * 60)
    print("五、各指标 FDR 校正 P 值的分位数分布")
    print("=" * 60)
    print(f"{'指标':35s} {'Min':>8} {'25%':>8} {'50%':>8} {'75%':>8} {'Max':>8}")
    for i, col in enumerate(PVALUE_COLS):
        vals = fdr_matrix[:, i]
        print(f"  {col:35s} {vals.min():8.4f} {np.percentile(vals,25):8.4f} {np.percentile(vals,50):8.4f} {np.percentile(vals,75):8.4f} {vals.max():8.4f}")

    print("\nDone.")

if __name__ == "__main__":
    main()
