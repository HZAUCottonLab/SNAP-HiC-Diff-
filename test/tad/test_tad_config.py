from copy import deepcopy
from pathlib import Path
import sys

import pytest
import yaml


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from tad_config import TadConfigError, parse_tad_settings


def repository_config():
    return yaml.safe_load((REPO / "config" / "config.yaml").read_text())


def test_repository_tad_contract_is_exact():
    settings = parse_tad_settings(repository_config())
    assert settings.resolutions == (4000, 8000, 20000)
    assert settings.signals == ("raw", "balanced")
    assert settings.threads == 20
    assert settings.excluded_chromosomes == ()
    assert settings.badregionfilter is False
    assert settings.parameter_set == "D"
    assert settings.d["resolution_overrides"][4000]["size_threshold_bins"] == 10


@pytest.mark.parametrize(
    "mutate,match",
    [
        (lambda config: config["tad"].update(resolutions=[4000, 20000]), "resolutions"),
        (lambda config: config["tad"].update(signals=["raw"]), "signals"),
        (lambda config: config["tad"].update(threads=0), "threads"),
        (lambda config: config["tad"].update(parameter_set="A"), "parameter_set"),
        (lambda config: config["tad"]["d"].pop("size_bins"), "size_bins"),
    ],
)
def test_invalid_tad_configuration_has_stable_parse_error(mutate, match):
    config = deepcopy(repository_config())
    mutate(config)
    with pytest.raises(TadConfigError, match=match) as error:
        parse_tad_settings(config)
    assert str(error.value).startswith("ERROR[TAD_CONFIG_INVALID]")


@pytest.mark.parametrize(
    "mutate",
    [
        lambda config: config["tad"].pop("run_root"),
        lambda config: config["tad"].pop("code_dir"),
        lambda config: config.pop("convert"),
        lambda config: config["convert"].pop("resolutions"),
        lambda config: config["tad"].update(d=None),
        lambda config: config["tad"]["d"].update(size_bins=None),
    ],
    ids=[
        "missing_run_root",
        "missing_code_dir",
        "missing_convert",
        "missing_convert_resolutions",
        "non_mapping_d",
        "malformed_size_bins",
    ],
)
def test_malformed_tad_structures_have_stable_parse_errors(mutate):
    config = deepcopy(repository_config())
    mutate(config)
    with pytest.raises(TadConfigError) as error:
        parse_tad_settings(config)


@pytest.mark.parametrize(
    ("field", "invalid"),
    [
        ("region_size_bins", 0),
        ("overlap_bins", 150),
        ("logged", 1),
        ("scale", ""),
        ("plateau", 0),
        ("chaosfilter", "true"),
        ("chaos_pct", 1.1),
        ("max_gamma", 0),
        ("num_part", 1.5),
        ("plots", 0),
        ("pctile_threshold", -0.1),
        ("pct_value", 1.1),
        ("final_consistent_domains", "true"),
        ("size_threshold_bins", 0),
        ("boundary_buffer_bins", -1),
        ("size_bins", [100, 200, 400, 750]),
        ("variance_thresholds", [45, 45, 60, 45, -1]),
    ],
)
def test_fixed_d_values_are_validated_before_compute(field, invalid):
    config = deepcopy(repository_config())
    config["tad"]["d"][field] = invalid
    with pytest.raises(TadConfigError, match="TAD_CONFIG_INVALID"):
        parse_tad_settings(config)

@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: d.update(resolution_overrides=[]),
        lambda d: d["resolution_overrides"][4000].update(size_threshold_bins=0),
        lambda d: d["resolution_overrides"][4000].update(variance_thresholds=[20] * 4),
    ],
)
def test_4kb_d_overrides_are_validated_before_compute(mutate):
    config = deepcopy(repository_config())
    mutate(config["tad"]["d"])
    with pytest.raises(TadConfigError, match="TAD_CONFIG_INVALID"):
        parse_tad_settings(config)
