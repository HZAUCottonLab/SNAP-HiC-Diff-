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
