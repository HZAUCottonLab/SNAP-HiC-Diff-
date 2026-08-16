from pathlib import Path
import subprocess
import sys

import cooler
import numpy as np
import pandas as pd
import pytest
import yaml


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "call_is_change_loci.py"
sys.path.insert(0, str(REPO / "scripts"))
from call_is_change_loci import read_balanced_bins

SCALES = (32000, 60000, 92000, 120000, 152000, 180000, 212000, 240000, 272000, 300000)
PRIMARY_COLUMNS = [
    "comparison", "region_id", "chrom", "event_start", "event_end", "summit",
    "is_change_type", "confidence_tier", "support_scale_bands", "support_scales",
    "max_amplitude_percentile", "max_transition_percentile", "max_pulse_percentile",
    "broad_shift_length_bp", "mask_status", "coverage_state_summary", "rank_within_region",
]


def write_config(root):
    config = {
        "meltron": {"cworld_is": {
            "source_resolution_bp": 4000, "is_distances_bp": list(SCALES),
            "ids_divisor": 2, "chrom_job_threads": 1, "merge_job_threads": 1,
        }},
        "is_change": {
            "raw_knn_pvalue_threshold": 0.05, "meltron_score_threshold": 5.0,
            "merge_gap_bp": 0, "context_bp": 300000, "local_radius_bins": 1,
            "min_local_valid_fraction": 2 / 3, "flank_fraction_of_scale": 0.25,
            "min_flank_bins": 2, "min_flank_valid_fraction": 0.80,
            "max_missing_run_bins": 1, "high_percentile": 0.90,
            "broad_percentile": 0.80, "spatial_cluster_tolerance_bp": 12000,
            "min_primary_scale_bands": 2,
            "scale_bands": {
                "small": [32000, 60000, 92000],
                "medium": [120000, 152000, 180000],
                "large": [212000, 240000, 272000, 300000],
            },
            "diagnostic_plot_top_n": 2,
        },
    }
    path = root / "config.yaml"
    path.write_text(yaml.safe_dump(config))
    return path


def make_balanced_coolers(root, weights_a, weights_b):
    starts = np.arange(len(weights_a)) * 4000
    bins = pd.DataFrame({"chrom": "chr1", "start": starts, "end": starts + 4000})
    pixels = pd.DataFrame({"bin1_id": np.arange(len(bins)), "bin2_id": np.arange(len(bins)), "count": 1})
    paths = []
    for name, weights in (("a.cool", weights_a), ("b.cool", weights_b)):
        path = root / name
        table = bins.assign(weight=np.asarray(weights, dtype=float))
        cooler.create_cooler(str(path), table, pixels, ordered=True)
        paths.append(path)
    return paths

def make_chromosome_cooler(root, name, chrom, weights):
    starts = np.arange(len(weights)) * 4000
    bins = pd.DataFrame({"chrom": chrom, "start": starts, "end": starts + 4000, "weight": weights})
    pixels = pd.DataFrame({"bin1_id": np.arange(len(bins)), "bin2_id": np.arange(len(bins)), "count": 1})
    path = root / name
    cooler.create_cooler(str(path), bins, pixels, ordered=True)
    return path

def test_read_balanced_bins_concatenates_chromosome_shards_and_rejects_duplicates(tmp_path):
    chr1 = make_chromosome_cooler(tmp_path, "chr1.cool", "chr1", [1.0, np.nan])
    chr2 = make_chromosome_cooler(tmp_path, "chr2.cool", "chr2", [1.0, 1.0])

    bins = read_balanced_bins([chr1, chr2], 4000)

    assert bins.to_dict("records") == [
        {"chrom": "chr1", "start": 0, "available": True},
        {"chrom": "chr1", "start": 4000, "available": False},
        {"chrom": "chr2", "start": 0, "available": True},
        {"chrom": "chr2", "start": 4000, "available": True},
    ]
    duplicate = make_chromosome_cooler(tmp_path, "chr1_duplicate.cool", "chr1", [1.0, 1.0])
    with pytest.raises(ValueError, match="duplicate chrom/start"):
        read_balanced_bins([chr1, duplicate], 4000)


def test_read_balanced_bins_accepts_terminal_partial_bin(tmp_path):
    bins = pd.DataFrame({
        "chrom": ["chr1", "chr1"],
        "start": [0, 4000],
        "end": [4000, 6500],
        "weight": [1.0, np.nan],
    })
    pixels = pd.DataFrame({"bin1_id": [0, 1], "bin2_id": [0, 1], "count": [1, 1]})
    path = tmp_path / "chr1_terminal_partial.cool"
    cooler.create_cooler(str(path), bins, pixels, ordered=True)

    observed = read_balanced_bins([path], 4000)

    assert observed.to_dict("records") == [
        {"chrom": "chr1", "start": 0, "available": True},
        {"chrom": "chr1", "start": 4000, "available": False},
    ]


def write_all_scale_is_tracks(root, n_bins=160):
    starts = np.arange(n_bins) * 4000
    rows_a, rows_b = [], []
    baseline = np.tile(np.array([-0.2, 0.2]), n_bins // 2)
    for scale in SCALES:
        delta = baseline.copy()
        if scale in (32000, 120000):  # two independent bands: transition at 320 kb
            delta[72:80] = -2.0
            delta[80:89] = 2.0
        elif scale == 60000:  # one-band pulse
            delta[48:50] = 2.0
        elif scale == 92000:  # broad-shift candidate, separate from the transition
            delta[100:114] = 2.0
        rows_a.extend(("chr1", start, scale, 0.0) for start in starts)
        rows_b.extend(("chr1", start, scale, value) for start, value in zip(starts, delta))
    columns = ["chrom", "viewpoint", "IS_distance", "IS"]
    a_path, b_path = root / "is_a.tsv", root / "is_b.tsv"
    pd.DataFrame(rows_a, columns=columns).to_csv(a_path, sep="\t", index=False)
    pd.DataFrame(rows_b, columns=columns).to_csv(b_path, sep="\t", index=False)
    return a_path, b_path


def write_candidate_regions(root, rows):
    path = root / "candidate_regions.tsv"
    pd.DataFrame(rows, columns=["region_id", "chrom", "region_start", "region_end"]).to_csv(
        path, sep="\t", index=False
    )
    return path


def run_caller(root, config_path, cool_a, cool_b, regions):
    return subprocess.run(
        [
            sys.executable, str(SCRIPT), "--candidate-regions", str(regions),
            "--is-a", str(root / "is_a.tsv"), "--is-b", str(root / "is_b.tsv"),
            "--cool-a", str(cool_a), "--cool-b", str(cool_b), "--config", str(config_path),
            "--out-primary", str(root / "primary.tsv"),
            "--out-preliminary", str(root / "preliminary.tsv"),
            "--out-masked", str(root / "masked.bed"), "--log", str(root / "caller.log"),
        ], text=True, capture_output=True, check=False,
    )


def test_caller_keeps_zero_masks_invalid_weight_and_calls_two_bands(tmp_path):
    config = write_config(tmp_path)
    cool_a, cool_b = make_balanced_coolers(
        tmp_path, [1.0, np.nan] + [1.0] * 158, [1.0] * 160
    )
    write_all_scale_is_tracks(tmp_path)
    regions = write_candidate_regions(
        tmp_path, [(1, "chr1", 0, 32000), (2, "chr1", 320000, 324000)]
    )

    result = run_caller(tmp_path, config, cool_a, cool_b, regions)

    assert result.returncode == 0, result.stderr
    masked = pd.read_csv(tmp_path / "masked.bed", sep="\t", header=None)
    assert (masked[1] == 4000).any()
    assert "A_EMPTY" in set(masked[3])
    primary = pd.read_csv(tmp_path / "primary.tsv", sep="\t")
    assert not primary.empty
    assert (primary["mask_status"] == "clean").all()
    assert set(primary["support_scale_bands"].str.split(",").explode()) >= {"small", "medium"}


def test_caller_writes_empty_schemas_and_not_evaluable_masks(tmp_path):
    config = write_config(tmp_path)
    cool_a, cool_b = make_balanced_coolers(tmp_path, [np.nan] * 160, [np.nan] * 160)
    write_all_scale_is_tracks(tmp_path)
    regions = write_candidate_regions(tmp_path, [(1, "chr1", 0, 32000)])

    result = run_caller(tmp_path, config, cool_a, cool_b, regions)

    assert result.returncode == 0, result.stderr
    assert pd.read_csv(tmp_path / "primary.tsv", sep="\t").columns.tolist() == PRIMARY_COLUMNS
    assert pd.read_csv(tmp_path / "preliminary.tsv", sep="\t").columns.tolist() == PRIMARY_COLUMNS
    reasons = set(pd.read_csv(tmp_path / "masked.bed", sep="\t", header=None)[3])
    assert "BOTH_EMPTY" in reasons


def test_caller_aligns_non_grid_context_to_track_coordinates(tmp_path):
    config = write_config(tmp_path)
    raw = yaml.safe_load(config.read_text())
    raw["is_change"]["context_bp"] = 300001
    config.write_text(yaml.safe_dump(raw))
    cool_a, cool_b = make_balanced_coolers(tmp_path, [1.0] * 160, [1.0] * 160)
    write_all_scale_is_tracks(tmp_path)
    regions = write_candidate_regions(tmp_path, [(1, "chr1", 320000, 324000)])

    result = run_caller(tmp_path, config, cool_a, cool_b, regions)

    assert result.returncode == 0, result.stderr
    assert not pd.read_csv(tmp_path / "primary.tsv", sep="\t").empty


def test_caller_removes_temp_file_after_serialization_failure(tmp_path, monkeypatch):
    sys.path.insert(0, str(REPO / "scripts"))
    from call_is_change_loci import atomic_write

    def fail(*args, **kwargs):
        raise OSError("write failed")

    monkeypatch.setattr(pd.DataFrame, "to_csv", fail)
    with pytest.raises(OSError, match="write failed"):
        atomic_write(pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), (
            tmp_path / "primary.tsv", tmp_path / "preliminary.tsv", tmp_path / "masked.bed",
        ))
    assert not list(tmp_path.glob(".*.tmp"))


def test_caller_ranks_region_events_by_evidence_before_coordinates(monkeypatch):
    sys.path.insert(0, str(REPO / "scripts"))
    import call_is_change_loci as caller

    clusters = pd.DataFrame([
        {
            "chrom": "chr1", "event_start": 4000, "event_end": 8000, "summit": 4000,
            "event_type": "transition", "confidence_tier": "primary",
            "support_scale_bands": "small", "support_scales": "32000",
            "max_amplitude_percentile": 0.99, "max_transition_percentile": 0.99,
            "max_pulse_percentile": 0.99,
        },
        {
            "chrom": "chr1", "event_start": 8000, "event_end": 12000, "summit": 8000,
            "event_type": "transition", "confidence_tier": "primary",
            "support_scale_bands": "small,medium", "support_scales": "32000,120000",
            "max_amplitude_percentile": 0.50, "max_transition_percentile": 0.50,
            "max_pulse_percentile": 0.50,
        },
    ])
    monkeypatch.setattr(caller, "_event_rows_for_region", lambda *args: (clusters, []))
    regions = pd.DataFrame([{"region_id": 1, "chrom": "chr1", "region_start": 0, "region_end": 12000}])

    primary, _, _ = caller.output_rows(regions, None, None, None, None, None, "A__vs__B")

    assert primary.loc[primary["event_start"] == 8000, "rank_within_region"].item() == 1
