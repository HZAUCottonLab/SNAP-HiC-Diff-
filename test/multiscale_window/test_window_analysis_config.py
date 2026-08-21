from pathlib import Path
import sys

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from window_analysis_config import WindowAnalysisConfigError, parse_window_analysis_settings


def make_config(tmp_path, resolution=4000, bins=(25, 50, 75, 100, 125), fraction=0.20):
    presets = []
    for window_bins in bins:
        path = tmp_path / f"r{resolution}_w{window_bins}.tsv"
        path.write_text("comparison\tchrom\twindow_start\twindow_end\n")
        step_bins = __import__("math").ceil(window_bins * fraction)
        presets.append({
            "resolution_bp": resolution,
            "window_bins": window_bins,
            "step_bins": step_bins,
            "file": str(path),
        })
    return {
        "resolution": resolution,
        "window_analysis": {"window_bins": list(bins), "step_fraction": fraction},
        "knn_empirical_pvalue": {"background_presets": presets},
    }


def test_default_4kb_window_math(tmp_path):
    settings = parse_window_analysis_settings(make_config(tmp_path))
    assert [s.window_size_bp for s in settings.scales] == [100000, 200000, 300000, 400000, 500000]
    assert [s.step_bins for s in settings.scales] == [5, 10, 15, 20, 25]
    assert [s.step_size_bp for s in settings.scales] == [20000, 40000, 60000, 80000, 100000]


def test_fractional_step_rounds_up(tmp_path):
    settings = parse_window_analysis_settings(make_config(tmp_path, bins=(17,)))
    assert settings.scales[0].step_bins == 4


def test_10kb_changes_only_physical_sizes(tmp_path):
    settings = parse_window_analysis_settings(make_config(tmp_path, resolution=10000))
    assert [s.window_size_bp for s in settings.scales] == [250000, 500000, 750000, 1000000, 1250000]
    assert [s.step_bins for s in settings.scales] == [5, 10, 15, 20, 25]


def test_missing_exact_preset_names_expected_tuple(tmp_path):
    config = make_config(tmp_path)
    config["knn_empirical_pvalue"]["background_presets"].pop()
    with pytest.raises(WindowAnalysisConfigError, match=r"r4000/w125.*step_bins=25"):
        parse_window_analysis_settings(config)


def test_legacy_only_window_config_is_rejected():
    with pytest.raises(WindowAnalysisConfigError, match="window_analysis"):
        parse_window_analysis_settings({"resolution": 4000, "window_size": 300000, "step_size": 32000})
