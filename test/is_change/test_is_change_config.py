from pathlib import Path
import sys

import pytest


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from cworld_is_config import parse_cworld_is_settings
from is_change_config import ConfigError, parse_is_change_settings


@pytest.fixture
def valid_config():
    return {
        "meltron": {
            "cworld_is": {
                "source_resolution_bp": 4000,
                "is_distances_bp": [
                    32000, 60000, 92000, 120000, 152000,
                    180000, 212000, 240000, 272000, 300000,
                ],
                "ids_divisor": 2,
                "chrom_job_threads": 4,
                "merge_job_threads": 4,
            }
        },
        "is_change": {
            "raw_knn_pvalue_threshold": 0.05,
            "meltron_score_threshold": 5.0,
            "merge_gap_bp": 0,
            "context_bp": 300000,
            "local_radius_bins": 1,
            "min_local_valid_fraction": 0.6666666667,
            "flank_fraction_of_scale": 0.25,
            "min_flank_bins": 2,
            "min_flank_valid_fraction": 0.80,
            "max_missing_run_bins": 1,
            "high_percentile": 0.95,
            "broad_percentile": 0.90,
            "spatial_cluster_tolerance_bp": 12000,
            "min_primary_scale_bands": 2,
            "scale_bands": {
                "small": [32000, 60000, 92000],
                "medium": [120000, 152000, 180000],
                "large": [212000, 240000, 272000, 300000],
            },
            "diagnostic_plot_top_n": 50,
        },
    }


def test_parses_scales_and_bands(valid_config):
    cworld = parse_cworld_is_settings(valid_config)
    settings = parse_is_change_settings(valid_config, cworld)
    assert settings.resolution_bp == 4000
    assert settings.merge_gap_bp == 0
    assert settings.min_primary_scale_bands == 2
    assert settings.scale_bands["small"] == (32000, 60000, 92000)


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda c: c["is_change"].update(merge_gap_bp=1), "merge_gap_bp"),
        (lambda c: c["is_change"].update(high_percentile=1.0), "high_percentile"),
        (
            lambda c: c["is_change"].update(min_flank_valid_fraction=0.79),
            "min_flank_valid_fraction",
        ),
        (lambda c: c["is_change"].update(scale_bands={"small": [32000]}), "partition"),
    ],
)
def test_rejects_invalid_values(valid_config, mutate, message):
    mutate(valid_config)
    with pytest.raises(ConfigError, match=message):
        parse_is_change_settings(valid_config, parse_cworld_is_settings(valid_config))
