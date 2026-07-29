from pathlib import Path
import sys

import pytest


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from cworld_is_config import ConfigError, parse_cworld_is_settings


DISTANCES = (
    32000, 60000, 92000, 120000, 152000,
    180000, 212000, 240000, 272000, 300000,
)
IDS = (
    16000, 30000, 46000, 60000, 76000,
    90000, 106000, 120000, 136000, 150000,
)


def valid_config():
    return {
        "meltron": {
            "cworld_is": {
                "source_resolution_bp": 4000,
                "is_distances_bp": list(DISTANCES),
                "ids_divisor": 2,
                "chrom_job_threads": 4,
                "merge_job_threads": 4,
            }
        }
    }


def test_parses_exact_scale_series_and_ids():
    settings = parse_cworld_is_settings(valid_config())

    assert settings.source_resolution_bp == 4000
    assert settings.is_distances_bp == DISTANCES
    assert settings.ids_distances_bp == IDS
    assert settings.chrom_job_threads == settings.merge_job_threads == 4


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda c: c["meltron"]["cworld_is"].update(source_resolution_bp=20000), "4000"),
        (lambda c: c["meltron"]["cworld_is"].update(is_distances_bp=list(DISTANCES[:-1])), "exactly 10"),
        (lambda c: c["meltron"]["cworld_is"].update(is_distances_bp=[*DISTANCES[:-1], 30001]), "divisible"),
        (lambda c: c["meltron"]["cworld_is"].update(is_distances_bp=[*DISTANCES[:-1], 272000]), "unique"),
        (lambda c: c["meltron"]["cworld_is"].update(ids_divisor=3), "ids_divisor"),
        (lambda c: c["meltron"]["cworld_is"].update(chrom_job_threads=0), "positive"),
    ],
)
def test_rejects_invalid_constraints(mutate, message):
    config = valid_config()
    mutate(config)
    with pytest.raises(ConfigError, match=message):
        parse_cworld_is_settings(config)
