import importlib.util
import json
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


@pytest.fixture
def d_config():
    return {
        "size_bins": [2, 4],
        "variance_thresholds": [20, 30],
        "region_size_bins": 10,
        "overlap_bins": 2,
        "logged": True,
        "scale": False,
        "plateau": 0.5,
        "chaosfilter": False,
        "chaos_pct": 0.1,
        "max_gamma": 1.5,
        "num_part": 2,
        "plots": False,
        "pctile_threshold": 0.95,
        "pct_value": 0.5,
        "final_consistent_domains": True,
        "size_threshold_bins": 3,
        "boundary_buffer_bins": 1,
    }


def cli_args(mcool_path, run_dir, ready_flag, d_config, **overrides):
    values = {
        "mcool": str(mcool_path),
        "sample": "sample-X",
        "signal": "raw",
        "resolution": "4000",
        "run_dir": str(run_dir),
        "ready_flag": str(ready_flag),
        "bad_region_bed": str(run_dir / "bad-regions.bed"),
        "badregionfilter": "true",
        "excluded_chromosomes": "",
        "d_config_json": json.dumps(d_config),
        "processors": "2",
    }
    values.update(overrides)
    arguments = []
    for key, value in values.items():
        if key == "excluded_chromosomes" and not value:
            continue
        arguments.extend([f"--{key.replace('_', '-')}", value])
    return arguments


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


def test_balanced_values_reject_negative_and_zero_endpoint_weights(tmp_path):
    adapter = load_adapter()
    path = tmp_path / "signed-weights.mcool"
    bins = pd.DataFrame({
        "chrom": ["A02"] * 4,
        "start": [0, 4000, 8000, 12000],
        "end": [4000, 8000, 12000, 16000],
        "weight": [-0.5, -0.25, 0.0, 0.5],
    })
    pixels = pd.DataFrame({
        "bin1_id": [0, 0, 3],
        "bin2_id": [1, 2, 3],
        "count": [8.0, 8.0, 8.0],
    })
    cooler.create_cooler(f"{path}::/resolutions/4000", bins, pixels)

    clr = cooler.Cooler(f"{path}::/resolutions/4000")
    assert list(adapter.iter_pixels_as_triples(clr, "balanced")) == [(4, 4, 2.0)]


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


def test_cli_creates_signal_artifacts_and_invalidates_stale_ready_before_bad_rerun(
    tiny_mcool, tmp_path, d_config, capsys
):
    adapter = load_adapter()
    run_dir = tmp_path / "run"
    ready_flag = run_dir / "ready.txt"

    raw_args = cli_args(tiny_mcool, run_dir, ready_flag, d_config, signal="raw")
    assert adapter.cli(raw_args) == 0
    assert capsys.readouterr().err == ""

    input_dir = run_dir / "input"
    raw_interactions = input_dir / "sample-X.raw.interactions"
    raw_settings = run_dir / "settings_sample-X.raw.txt"
    assert raw_interactions.read_text() == "1\t1\t8.0\n1\t2\t4.0\n1\t3\t4.0\n"
    assert "counts_file_1\tsample-X.raw.interactions\n" in raw_settings.read_text()
    assert "signal=raw\n" in ready_flag.read_text()
    assert (input_dir / "sample-X.bed").read_text().endswith("A03\t0\t4000\t3\n")

    balanced_args = cli_args(
        tiny_mcool, run_dir, ready_flag, d_config, signal="balanced"
    )
    assert adapter.cli(balanced_args) == 0
    assert capsys.readouterr().err == ""
    assert (input_dir / "sample-X.balanced.interactions").read_text() == (
        "1\t1\t0.5\n1\t2\t0.5\n"
    )
    assert "counts_file_1\tsample-X.balanced.interactions\n" in (
        run_dir / "settings_sample-X.balanced.txt"
    ).read_text()
    assert "signal=balanced\n" in ready_flag.read_text()

    raw_interactions.write_text("previous raw interactions\n")
    raw_settings.write_text("previous raw settings\n")
    ready_flag.write_text("stale ready manifest\n")
    bad_json_args = cli_args(
        tiny_mcool,
        run_dir,
        ready_flag,
        d_config,
        d_config_json="{not-valid-json",
    )
    assert adapter.cli(bad_json_args) == 1
    error = capsys.readouterr().err
    assert error.startswith("ERROR[TAD_CONFIG_INVALID]")
    assert "sample=sample-X" in error
    assert "signal=raw" in error
    assert "resolution=4000" in error
    assert f"path={run_dir.resolve()}" in error
    assert not ready_flag.exists()
    assert raw_interactions.read_text() == "previous raw interactions\n"
    assert raw_settings.read_text() == "previous raw settings\n"


def test_changed_d_config_archives_old_stage_artifacts_before_new_prepare(
    tiny_mcool, tmp_path, d_config, capsys
):
    adapter = load_adapter()
    run_dir = tmp_path / "run"
    ready_flag = tmp_path / "prepare" / "ready.txt"

    assert adapter.cli(cli_args(tiny_mcool, run_dir, ready_flag, d_config)) == 0
    assert capsys.readouterr().err == ""

    previous_calls = run_dir / "FINAL_DOMAIN_CALLS" / "old_master_consistent_domains.txt"
    previous_calls.parent.mkdir()
    previous_calls.write_text("A02\t0\t4000\n")
    previous_log = run_dir / "GPS" / "stage.log"
    previous_log.parent.mkdir()
    previous_log.write_text("failed-stage-diagnostic\n")

    changed_d_config = dict(d_config, size_threshold_bins=99)
    assert adapter.cli(
        cli_args(tiny_mcool, run_dir, ready_flag, changed_d_config)
    ) == 0
    assert capsys.readouterr().err == ""

    archives = sorted((run_dir / ".previous_attempts").iterdir())
    assert len(archives) == 1
    assert (archives[0] / "FINAL_DOMAIN_CALLS" / previous_calls.name).read_text() == (
        "A02\t0\t4000\n"
    )
    assert (archives[0] / "GPS" / "stage.log").read_text() == (
        "failed-stage-diagnostic\n"
    )
    assert not (run_dir / "FINAL_DOMAIN_CALLS").exists()
    assert not (run_dir / "GPS").exists()


@pytest.mark.parametrize(
    ("overrides", "expected_context"),
    [
        ({"signal": "unsupported"}, "signal=unsupported"),
        ({"resolution": "not-an-integer"}, "sample=sample-X"),
        ({"badregionfilter": "not-a-boolean"}, "signal=raw"),
    ],
)
def test_cli_configuration_errors_use_stable_diagnostic(
    tiny_mcool, tmp_path, d_config, capsys, overrides, expected_context
):
    adapter = load_adapter()
    run_dir = tmp_path / "run"
    result = adapter.cli(
        cli_args(tiny_mcool, run_dir, run_dir / "ready.txt", d_config, **overrides)
    )

    error = capsys.readouterr().err
    assert result == 1
    assert error.startswith("ERROR[TAD_CONFIG_INVALID]")
    assert expected_context in error
    assert f"path={run_dir.resolve()}" in error
