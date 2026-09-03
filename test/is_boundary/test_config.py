from copy import deepcopy
from pathlib import Path
import sys

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from cworld_is_config import parse_cworld_is_settings
from is_boundary_config import parse_is_boundary_settings
from tad_config import parse_tad_settings


def repository_config():
    return yaml.safe_load((REPO / "config" / "config.yaml").read_text())


def test_repository_boundary_annotation_contract_is_exact():
    config = repository_config()
    settings = parse_is_boundary_settings(
        config,
        parse_cworld_is_settings(config),
        parse_tad_settings(config),
    )
    assert settings.base_bin_bp == 4000
    assert settings.jitter_quantile == 0.95
    assert settings.jitter_min_pairs == 100
    assert settings.jitter_fallback_bp == 20000
    assert settings.max_shift_bp == 100000
    assert settings.shift_sensitivity_bp == (60000, 80000, 100000)
    assert settings.match_tie_tolerance_bp == 4000
    assert settings.strength.robust_z_threshold == 2.0
    assert settings.strength.stable_abs_z_threshold == 1.0
    assert settings.permutation.enabled is False
    assert settings.permutation.n_permutations == 10000


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("base_bin_bp", 8000, "base_bin_bp"),
        ("jitter_quantile", 1.0, "jitter_quantile"),
        ("jitter_min_pairs", 0, "jitter_min_pairs"),
        ("max_shift_bp", 90000, "shift_sensitivity_bp"),
        ("match_tie_tolerance_bp", 2000, "4 kb"),
    ],
)
def test_rejects_invalid_boundary_annotation_values(field, value, message):
    config = deepcopy(repository_config())
    config["is_boundary_annotation"][field] = value
    with pytest.raises(ValueError, match=message):
        parse_is_boundary_settings(
            config,
            parse_cworld_is_settings(config),
            parse_tad_settings(config),
        )


def test_permutation_cannot_be_enabled_when_annotation_is_disabled():
    config = deepcopy(repository_config())
    config["is_boundary_annotation"]["enabled"] = False
    config["is_boundary_annotation"]["permutation"]["enabled"] = True
    with pytest.raises(ValueError, match="permutation.enabled"):
        parse_is_boundary_settings(
            config,
            parse_cworld_is_settings(config),
            parse_tad_settings(config),
        )


@pytest.mark.parametrize(
    "path",
    [
        ("base_bin_bp",),
        ("strength", "robust_z_threshold"),
        ("permutation", "seed"),
    ],
)
def test_rejects_missing_boundary_annotation_keys_with_stable_error(path):
    config = deepcopy(repository_config())
    value = config["is_boundary_annotation"]
    for key in path[:-1]:
        value = value[key]
    del value[path[-1]]
    with pytest.raises(ValueError, match=r"ERROR\[IS_BOUNDARY_CONFIG_INVALID\]"):
        parse_is_boundary_settings(
            config,
            parse_cworld_is_settings(config),
            parse_tad_settings(config),
        )


@pytest.mark.parametrize(
    "sensitivity",
    [
        [],
        100000,
        [-4000, 60000, 80000, 100000],
    ],
)
def test_rejects_empty_wrong_or_negative_shift_sensitivity(sensitivity):
    config = deepcopy(repository_config())
    config["is_boundary_annotation"]["shift_sensitivity_bp"] = sensitivity
    with pytest.raises(
        ValueError, match=r"ERROR\[IS_BOUNDARY_CONFIG_INVALID\].*shift_sensitivity_bp"
    ):
        parse_is_boundary_settings(
            config,
            parse_cworld_is_settings(config),
            parse_tad_settings(config),
        )


@pytest.mark.parametrize(
    ("field", "value", "sensitivity_thresholds"),
    [
        ("robust_z_threshold", float("inf"), [1.5, 2.0, float("inf")]),
        ("stable_abs_z_threshold", -1.0, None),
        ("sensitivity_thresholds", [1.5, 2.0, float("inf")], None),
        ("sensitivity_thresholds", [-1.0, 1.5, 2.0], None),
    ],
)
def test_rejects_nonfinite_or_negative_strength_thresholds(
    field, value, sensitivity_thresholds
):
    config = deepcopy(repository_config())
    strength = config["is_boundary_annotation"]["strength"]
    strength[field] = value
    if sensitivity_thresholds is not None:
        strength["sensitivity_thresholds"] = sensitivity_thresholds
    with pytest.raises(ValueError, match=r"ERROR\[IS_BOUNDARY_CONFIG_INVALID\]"):
        parse_is_boundary_settings(
            config,
            parse_cworld_is_settings(config),
            parse_tad_settings(config),
        )


def test_rejects_strength_threshold_that_overflows_float_conversion():
    config = deepcopy(repository_config())
    huge_integer = 10**1000
    strength = config["is_boundary_annotation"]["strength"]
    strength["robust_z_threshold"] = huge_integer
    strength["sensitivity_thresholds"] = [1.5, 2.0, huge_integer]
    with pytest.raises(
        ValueError,
        match=r"ERROR\[IS_BOUNDARY_CONFIG_INVALID\].*must be finite and positive",
    ):
        parse_is_boundary_settings(
            config,
            parse_cworld_is_settings(config),
            parse_tad_settings(config),
        )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("shift_sensitivity_bp", [60000, 100000]),
        ("shift_sensitivity_bp", [60000, 80000, 120000]),
    ],
)
def test_rejects_nonpublication_shift_sensitivity(field, value):
    config = deepcopy(repository_config())
    config["is_boundary_annotation"][field] = value

    with pytest.raises(ValueError, match="fixed publication contract"):
        parse_is_boundary_settings(
            config, parse_cworld_is_settings(config), parse_tad_settings(config)
        )


def test_rejects_nonpublication_strength_sensitivity():
    config = deepcopy(repository_config())
    config["is_boundary_annotation"]["strength"]["sensitivity_thresholds"] = [1.0, 2.0, 3.0]

    with pytest.raises(ValueError, match="fixed publication contract"):
        parse_is_boundary_settings(
            config, parse_cworld_is_settings(config), parse_tad_settings(config)
        )
