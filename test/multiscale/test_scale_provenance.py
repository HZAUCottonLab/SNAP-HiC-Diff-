from pathlib import Path

import pandas as pd

from scripts.compute_knn_empirical_pvalue import (
    SCALE_COLUMNS,
    attach_scale_metadata,
    empty_knn_output,
)
from scripts.plot_knn_significant_heatmap import scale_h5_path


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
