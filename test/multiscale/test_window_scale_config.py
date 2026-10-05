from copy import deepcopy

import pytest

from scripts.window_scale_config import build_window_scales


def valid_config():
    root = "/share/home/xhhuang24/project/PanTE/background_pool/results/multiscale_4k/window_scales"
    bins_and_labels = ((25, "100kb"), (50, "200kb"), (75, "300kb"), (100, "400kb"), (125, "500kb"))
    return {
        "resolution": 4000,
        "window_analysis": {
            "window_bins": [value for value, _ in bins_and_labels],
            "step_fraction": 0.2,
        },
        "knn_empirical_pvalue": {
            "background_presets": [
                {
                    "resolution_bp": 4000,
                    "window_bins": value,
                    "step_bins": value // 5,
                    "file": f"{root}/{label}/window_compare/HiC_window_background.tsv",
                }
                for value, label in bins_and_labels
            ]
        },
    }


def test_build_window_scales_returns_the_approved_order_and_geometry():
    scales = build_window_scales(valid_config())
    assert list(scales) == ["w25", "w50", "w75", "w100", "w125"]
    assert scales["w25"] == {
        "resolution_bp": 4000,
        "window_bins": 25,
        "window_size_bp": 100000,
        "step_bins": 5,
        "step_size_bp": 20000,
        "background_file": (
            "/share/home/xhhuang24/project/PanTE/background_pool/results/"
            "multiscale_4k/window_scales/100kb/window_compare/"
            "HiC_window_background.tsv"
        ),
    }
    assert scales["w125"]["step_size_bp"] == 100000


def test_build_window_scales_rounds_fractional_steps_up():
    config = valid_config()
    config["window_analysis"]["window_bins"] = [16]
    config["knn_empirical_pvalue"]["background_presets"] = [{
        "resolution_bp": 4000,
        "window_bins": 16,
        "step_bins": 4,
        "file": "/background/r4000_w16.tsv",
    }]
    assert build_window_scales(config)["w16"]["step_bins"] == 4


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda c: c.update(resolution=0), "resolution must be a positive integer"),
        (lambda c: c.update(window_size=300000), "remove legacy config key: window_size"),
        (lambda c: c.update(step_size=32000), "remove legacy config key: step_size"),
        (lambda c: c["window_analysis"].update(window_bins=[]), "window_bins must be a non-empty list"),
        (lambda c: c["window_analysis"].update(window_bins=[25, 25]), "window_bins must not contain duplicates"),
        (lambda c: c["window_analysis"].update(step_fraction=float("inf")), "step_fraction must be finite"),
        (lambda c: c["window_analysis"].update(step_fraction=0), "step_fraction must be greater than 0 and at most 1"),
        (lambda c: c["knn_empirical_pvalue"].update(background_presets=[]), "background_presets must be a non-empty list"),
        (
            lambda c: c["knn_empirical_pvalue"]["background_presets"].pop(),
            "missing background preset for resolution_bp=4000, window_bins=125, step_bins=25",
        ),
        (
            lambda c: c["knn_empirical_pvalue"]["background_presets"].append(
                deepcopy(c["knn_empirical_pvalue"]["background_presets"][0])
            ),
            "duplicate background preset",
        ),
        (
            lambda c: c["knn_empirical_pvalue"]["background_presets"][0].update(step_bins=6),
            "unused background preset for active resolution",
        ),
    ],
)
def test_build_window_scales_rejects_invalid_configuration(mutate, message):
    config = valid_config()
    mutate(config)
    with pytest.raises(ValueError, match=message):
        build_window_scales(config)
