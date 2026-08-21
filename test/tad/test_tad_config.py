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
