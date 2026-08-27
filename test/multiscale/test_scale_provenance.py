from pathlib import Path

import pandas as pd

from scripts.compute_knn_empirical_pvalue import (
    SCALE_COLUMNS,
    attach_scale_metadata,
    empty_knn_output,
)
from scripts.plot_knn_significant_heatmap import METRICS, run_scan, scale_h5_path


def test_attach_scale_metadata_adds_the_full_constant_identity():
    frame = pd.DataFrame([{
        "chrom": "HC04_A01",
        "window_id": 7,
        "window_start": 20000,
        "window_end": 120000,
        "comparison": "a_vs_b",
    }])
    result = attach_scale_metadata(frame, 4000, 25, 5)
    assert result.loc[0, SCALE_COLUMNS].to_dict() == {
        "resolution_bp": 4000,
        "window_bins": 25,
        "window_size_bp": 100000,
        "step_bins": 5,
        "step_size_bp": 20000,
    }


def test_empty_knn_output_keeps_the_final_schema():
    query = pd.DataFrame(columns=[
        "chrom", "window_id", "window_start", "window_end", "comparison",
    ])
    query = attach_scale_metadata(query, 4000, 25, 5)
    result = empty_knn_output(query, ["MSE", "d_SSIM"])
    assert result.empty
    assert result.columns.tolist() == [
        *SCALE_COLUMNS,
        "chrom", "window_id", "window_start", "window_end", "comparison",
        "knn_valid_k", "knn_dist_median", "knn_dist_max",
        "MSE_empirical_pvalue", "bg_MSE_median",
        "d_SSIM_empirical_pvalue", "bg_d_SSIM_median",
    ]


def test_scale_h5_path_never_falls_back_to_legacy_results_windows():
    root = Path("results/window_scales/r4000/w25")
    assert scale_h5_path(root, "H1394", "HC04_A01", "gaussian") == (
        root / "windows/H1394/H1394__HC04_A01_gaussian.h5"
    )


def test_empty_knn_scan_writes_a_valid_empty_summary(tmp_path):
    query = pd.DataFrame(columns=[
        "chrom", "window_id", "window_start", "window_end", "comparison",
    ])
    empty_knn = empty_knn_output(
        attach_scale_metadata(query, 4000, 25, 5), METRICS
    )
    enriched_path = tmp_path / "empty_knn.tsv"
    empty_knn.to_csv(enriched_path, sep="\t", index=False)

    out_dir = tmp_path / "scan"
    run_scan(enriched_path, out_dir)

    assert {
        "significant_windows.tsv", "pval_table.tsv", "pval_summary.tsv",
    } == {path.name for path in out_dir.iterdir()}

    summary = pd.read_csv(out_dir / "pval_summary.tsv", sep="\t")
    assert summary.columns.tolist() == [
        "metric", "n_total", "n_pass", "pct_pass", "pass_median_p",
        "fail_median_p", "pass_min_p", "pass_max_p", "n_pass_0_metrics",
        "n_pass_1_metric", "n_pass_2_metrics", "n_pass_3_metrics",
        "n_pass_4_metrics", "n_pass_5_metrics", "n_pass_any",
        "n_pass_all5", "pvalue_threshold",
    ]
    assert summary["metric"].tolist() == METRICS
    assert (summary[["n_total", "n_pass"]] == 0).all().all()
    assert summary[[
        "n_pass_0_metrics", "n_pass_1_metric", "n_pass_2_metrics",
        "n_pass_3_metrics", "n_pass_4_metrics", "n_pass_5_metrics",
        "n_pass_any", "n_pass_all5",
    ]].eq(0).all().all()
    assert summary[[
        "pct_pass", "pass_median_p", "fail_median_p", "pass_min_p",
        "pass_max_p",
    ]].isna().all().all()
    assert summary["pvalue_threshold"].eq(0.05).all()
