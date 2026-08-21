import importlib.util
from pathlib import Path
import sys

import cooler
import numpy as np
import pandas as pd
import pytest


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from tad_errors import TadError


def load_adapter():
    path = REPO / "scripts" / "tad_from_mcool.py"
    spec = importlib.util.spec_from_file_location("tad_from_mcool", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def tiny_mcool(tmp_path):
    path = tmp_path / "tiny.mcool"
    bins = pd.DataFrame({
        "chrom": ["A02", "A02", "A03"],
        "start": [0, 4000, 0],
        "end": [4000, 8000, 4000],
        "weight": [0.25, 0.5, np.nan],
    })
    pixels = pd.DataFrame({
        "bin1_id": [0, 0, 0, 1],
        "bin2_id": [0, 1, 2, 1],
        "count": [8.0, 4.0, 4.0, 0.0],
    })
    cooler.create_cooler(
        f"{path}::/resolutions/4000", bins, pixels,
        dtypes={"count": "float64"},
    )
    return path


def test_raw_values_and_balanced_values_with_nonfinite_weight_filter(tiny_mcool):
    adapter = load_adapter()
    clr = cooler.Cooler(f"{tiny_mcool}::/resolutions/4000")
    assert list(adapter.iter_pixels_as_triples(clr, "raw")) == [
        (1, 1, 8.0),
        (1, 2, 4.0),
        (1, 3, 4.0),
    ]
    assert list(adapter.iter_pixels_as_triples(clr, "balanced")) == [
        (1, 1, 0.5),
        (1, 2, 0.5),
    ]


def test_4kb_override_is_identical_for_both_signal_runs():
    adapter = load_adapter()
    d_config = {
        "size_threshold_bins": 4,
        "variance_thresholds": [45, 45, 60, 45, 45],
        "resolution_overrides": {
            "4000": {
                "size_threshold_bins": 10,
                "variance_thresholds": [20, 20, 30, 20, 20],
            }
        },
    }
    raw_effective = adapter.effective_d_config(d_config, 4000)
    balanced_effective = adapter.effective_d_config(d_config, 4000)
    assert raw_effective == balanced_effective
    assert raw_effective["size_threshold_bins"] == 10
    assert raw_effective["variance_thresholds"] == [20, 20, 30, 20, 20]


def test_balanced_mode_reports_dynamic_missing_weight_context(tmp_path):
    adapter = load_adapter()
    path = tmp_path / "unbalanced.mcool"
    bins = pd.DataFrame({
        "chrom": ["A02"], "start": [0], "end": [4000]
    })
    pixels = pd.DataFrame({"bin1_id": [0], "bin2_id": [0], "count": [1]})
    cooler.create_cooler(f"{path}::/resolutions/4000", bins, pixels)

    with pytest.raises(TadError) as error:
        adapter.validate_mcool_resolution(path, 4000, "sample-X", "balanced")
    assert str(error.value).startswith("ERROR[TAD_MISSING_WEIGHT]")
    assert "sample=sample-X" in str(error.value)
    assert "signal=balanced" in str(error.value)
    assert "resolution=4000" in str(error.value)


def test_missing_resolution_and_empty_contacts_have_stable_codes(tiny_mcool, tmp_path):
    adapter = load_adapter()
    with pytest.raises(TadError, match=r"ERROR\[TAD_MISSING_RESOLUTION\]"):
        adapter.validate_mcool_resolution(tiny_mcool, 8000, "S", "raw")

    clr = cooler.Cooler(f"{tiny_mcool}::/resolutions/4000")
    clr.pixels()[:]
    with pytest.raises(TadError, match=r"ERROR\[TAD_EMPTY_INTERACTIONS\]"):
        adapter.write_interactions(
            tmp_path / "empty.interactions",
            [],
            sample="S",
            signal="raw",
            resolution=4000,
        )


def test_global_bed_ids_are_one_based(tiny_mcool, tmp_path):
    adapter = load_adapter()
    clr = cooler.Cooler(f"{tiny_mcool}::/resolutions/4000")
    output = tmp_path / "S.bed"
    assert adapter.write_bed_from_bins(output, clr) == 3
    assert output.read_text() == (
        "A02\t0\t4000\t1\n"
        "A02\t4000\t8000\t2\n"
        "A03\t0\t4000\t3\n"
    )
