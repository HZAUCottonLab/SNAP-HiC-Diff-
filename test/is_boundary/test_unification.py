from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import cooler
import numpy as np
import pandas as pd
import pytest


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from is_boundary_io import (
    BOUNDARY_VIEW_COLUMNS,
    MULTIRES_INPUT_COLUMNS,
    TAD_UNEVALUABLE_COLUMNS,
    atomic_write_tsv,
    read_multires_boundaries,
    read_resolution_boundaries,
)
from is_boundary_unification import (
    build_jitter_offsets,
    build_sensitivity_views,
    maximum_cardinality_pairs,
    unify_boundaries,
)
from unify_tad_boundaries import _tad_unevaluable


SCRIPT = REPO / "scripts" / "unify_tad_boundaries.py"


def write_bed7(path, rows):
    path.write_text("".join("\t".join(map(str, row)) + "\n" for row in rows))
    return path


def write_bed6(path, rows):
    path.write_text("".join("\t".join(map(str, row)) + "\n" for row in rows))
    return path


def boundary_frame(signal, rows):
    return pd.DataFrame([
        {
            "chrom": chrom, "start": start, "end": end, "boundary_id": boundary_id,
            "support_resolutions": support_resolutions,
            "n_supporting_resolutions": len(support_resolutions.split(",")),
            "n_source_boundaries": count, "signal": signal,
        }
        for chrom, start, end, boundary_id, support_resolutions, count in rows
    ])


def test_reads_bed7_and_unifies_overlap_with_provenance(tmp_path):
    raw_path = write_bed7(tmp_path / "S.multires.raw.boundaries.bed", [
        ("A02", 0, 12000, "raw_1", "4000,8000", 2, 3),
        ("A02", 40000, 48000, "raw_2", "4000", 1, 1),
    ])
    balanced_path = write_bed7(tmp_path / "S.multires.balanced.boundaries.bed", [
        ("A02", 4000, 16000, "bal_1", "4000,20000", 2, 2),
    ])
    raw = read_multires_boundaries(raw_path, "raw")
    balanced = read_multires_boundaries(balanced_path, "balanced")

    observed = unify_boundaries(raw, balanced, tie_tolerance_bp=4000)

    assert observed[[
        "chrom", "start", "end", "signal_support", "support_resolutions",
        "n_source_boundaries", "component_state",
    ]].to_dict("records") == [
        {
            "chrom": "A02", "start": 0, "end": 16000,
            "signal_support": "both", "support_resolutions": "4000,8000,20000",
            "n_source_boundaries": 5, "component_state": "simple",
        },
        {
            "chrom": "A02", "start": 40000, "end": 48000,
            "signal_support": "raw_only", "support_resolutions": "4000",
            "n_source_boundaries": 1, "component_state": "simple",
        },
    ]


def test_adjacency_does_not_unify_and_many_to_many_is_complex(tmp_path):
    raw = pd.DataFrame([
        {"chrom": "A02", "start": 0, "end": 8000, "boundary_id": "r1", "support_resolutions": "4000", "n_supporting_resolutions": 1, "n_source_boundaries": 1, "signal": "raw"},
        {"chrom": "A02", "start": 6000, "end": 14000, "boundary_id": "r2", "support_resolutions": "4000", "n_supporting_resolutions": 1, "n_source_boundaries": 1, "signal": "raw"},
        {"chrom": "A02", "start": 30000, "end": 38000, "boundary_id": "r3", "support_resolutions": "4000", "n_supporting_resolutions": 1, "n_source_boundaries": 1, "signal": "raw"},
    ])
    balanced = pd.DataFrame([
        {"chrom": "A02", "start": 4000, "end": 10000, "boundary_id": "b1", "support_resolutions": "4000", "n_supporting_resolutions": 1, "n_source_boundaries": 1, "signal": "balanced"},
        {"chrom": "A02", "start": 8000, "end": 16000, "boundary_id": "b2", "support_resolutions": "4000", "n_supporting_resolutions": 1, "n_source_boundaries": 1, "signal": "balanced"},
        {"chrom": "A02", "start": 38000, "end": 46000, "boundary_id": "b3", "support_resolutions": "4000", "n_supporting_resolutions": 1, "n_source_boundaries": 1, "signal": "balanced"},
    ])

    observed = unify_boundaries(raw, balanced, tie_tolerance_bp=4000)

    assert (observed["component_state"] == "complex_within_sample").sum() == 1
    complex_row = observed.loc[
        observed["component_state"] == "complex_within_sample"
    ].iloc[0]
    assert complex_row[["raw_boundary_ids", "balanced_boundary_ids"]].to_dict() == {
        "raw_boundary_ids": "r1,r2",
        "balanced_boundary_ids": "b1,b2",
    }
    assert set(observed.loc[observed.start >= 30000, "signal_support"]) == {"raw_only", "balanced_only"}


def test_readers_reject_bad_bed_schema_and_atomic_writer_removes_temporary_file(tmp_path, monkeypatch):
    malformed = write_bed7(tmp_path / "malformed.bed", [
        ("A02", 0, 8000, "b1", "4000", 1, 1, "extra"),
    ])
    with pytest.raises(RuntimeError, match="IS_BOUNDARY_INPUT_SCHEMA_INVALID"):
        read_multires_boundaries(malformed, "raw")

    bad_resolution = write_bed6(tmp_path / "bad_resolution.bed", [
        ("A02", 0, 8000, "b1", 8000, 1),
    ])
    with pytest.raises(RuntimeError, match="IS_BOUNDARY_INPUT_SCHEMA_INVALID"):
        read_resolution_boundaries([(4000, bad_resolution)], "raw")

    output = tmp_path / "output.tsv"

    def fail(*args, **kwargs):
        raise OSError("write failed")

    monkeypatch.setattr(pd.DataFrame, "to_csv", fail)
    with pytest.raises(OSError, match="write failed"):
        atomic_write_tsv(pd.DataFrame(columns=MULTIRES_INPUT_COLUMNS), output, MULTIRES_INPUT_COLUMNS)
    assert not list(tmp_path.glob(".output.tsv.*.tmp"))


def test_sensitivity_views_keep_source_inputs_and_drop_leave_one_resolution_out():
    raw = boundary_frame("raw", [
        ("A02", 0, 8000, "raw_4", "4000", 1),
        ("A02", 16000, 24000, "raw_8", "8000", 1),
        ("A02", 32000, 40000, "raw_20", "20000", 1),
    ])
    balanced = boundary_frame("balanced", [
        ("A02", 0, 8000, "balanced_4", "4000", 1),
        ("A02", 16000, 24000, "balanced_8", "8000", 1),
        ("A02", 32000, 40000, "balanced_20", "20000", 1),
    ])
    per_resolution = pd.DataFrame([
        {"chrom": chrom, "start": start, "end": end, "boundary_id": boundary_id,
         "resolution": resolution, "n_source_boundaries": 1, "signal": signal}
        for signal, prefix in (("raw", "raw"), ("balanced", "balanced"))
        for chrom, start, end, boundary_id, resolution in (
            ("A02", 0, 8000, f"{prefix}_r4", 4000),
            ("A02", 16000, 24000, f"{prefix}_r8", 8000),
            ("A02", 32000, 40000, f"{prefix}_r20", 20000),
        )
    ])

    observed = build_sensitivity_views(
        "S", raw, balanced, per_resolution,
        SimpleNamespace(match_tie_tolerance_bp=4000),
    )

    assert observed.columns.tolist() == list(BOUNDARY_VIEW_COLUMNS)
    assert set(observed["view"]) == {
        "unified", "raw_only_input", "balanced_only_input", "leave_r4000_out",
        "leave_r8000_out", "leave_r20000_out",
    }
    raw_input = observed[observed["view"] == "raw_only_input"]
    assert raw_input[["start", "end", "raw_boundary_ids", "balanced_boundary_ids"]].to_dict("records") == [
        {"start": 0, "end": 8000, "raw_boundary_ids": "raw_4", "balanced_boundary_ids": ""},
        {"start": 16000, "end": 24000, "raw_boundary_ids": "raw_8", "balanced_boundary_ids": ""},
        {"start": 32000, "end": 40000, "raw_boundary_ids": "raw_20", "balanced_boundary_ids": ""},
    ]
    assert not ((observed["view"] == "leave_r4000_out") & (observed["start"] == 0)).any()
    assert ((observed["view"] == "leave_r8000_out") & (observed["start"] == 0)).any()
    assert ((observed["view"] == "leave_r20000_out") & (observed["start"] == 0)).any()


def test_jitter_offsets_include_raw_balanced_pair_and_propagate_exclusions():
    raw = boundary_frame("raw", [("A02", 0, 8000, "raw_1", "4000", 1)])
    balanced = boundary_frame("balanced", [("A02", 4000, 12000, "balanced_1", "4000", 1)])
    per_resolution = pd.DataFrame(columns=[
        "chrom", "start", "end", "boundary_id", "resolution",
        "n_source_boundaries", "signal",
    ])

    observed = build_jitter_offsets(raw, balanced, per_resolution, base_bin_bp=4000)

    assert observed.columns.tolist() == [
        "sample", "chrom", "source", "comparison", "resolution_a", "resolution_b",
        "boundary_a_id", "boundary_b_id", "center_a", "center_b", "offset_bp", "eligible",
    ]
    signal_pair = observed[observed["comparison"] == "raw_vs_balanced"].iloc[0]
    assert signal_pair[["boundary_a_id", "boundary_b_id", "offset_bp", "eligible"]].to_dict() == {
        "boundary_a_id": "raw_1", "boundary_b_id": "balanced_1",
        "offset_bp": 4000, "eligible": True,
    }

    excluded = build_jitter_offsets(
        raw.assign(touches_chromosome_end=True), balanced, per_resolution, base_bin_bp=4000
    )
    assert not excluded.loc[excluded["comparison"] == "raw_vs_balanced", "eligible"].iloc[0]

    complex_offsets = build_jitter_offsets(
        boundary_frame("raw", [
            ("A02", 0, 8000, "raw_a", "4000", 1),
            ("A02", 4000, 12000, "raw_b", "4000", 1),
        ]),
        boundary_frame("balanced", [
            ("A02", 2000, 10000, "balanced_a", "4000", 1),
            ("A02", 6000, 14000, "balanced_b", "4000", 1),
        ]),
        per_resolution,
        base_bin_bp=4000,
    )
    assert not complex_offsets.loc[
        complex_offsets["comparison"] == "raw_vs_balanced", "eligible"
    ].any()

    near_tie_raw = boundary_frame("raw", [
        ("A02", 0, 8000, "near_tie_raw", "4000", 1),
    ])
    near_tie_balanced = boundary_frame("balanced", [
        ("A02", 0, 8000, "near_tie_balanced_1", "4000", 1),
        ("A02", 2000, 10000, "near_tie_balanced_2", "4000", 1),
    ])
    assert (
        unify_boundaries(near_tie_raw, near_tie_balanced, tie_tolerance_bp=4000)
        ["component_state"].eq("complex_within_sample").all()
    )
    near_tie_offsets = build_jitter_offsets(
        near_tie_raw, near_tie_balanced, per_resolution, base_bin_bp=4000
    )
    assert not near_tie_offsets.loc[
        near_tie_offsets["comparison"] == "raw_vs_balanced", "eligible"
    ].iloc[0]

    resolution_metadata = pd.DataFrame([
        {"chrom": "A02", "start": 0, "end": 8000, "boundary_id": "raw_r4",
         "resolution": 4000, "n_source_boundaries": 1, "signal": "raw",
         "has_nonfinite_base_bins": False},
        {"chrom": "A02", "start": 0, "end": 8000, "boundary_id": "raw_r8",
         "resolution": 8000, "n_source_boundaries": 1, "signal": "raw",
         "has_nonfinite_base_bins": True},
        {"chrom": "A02", "start": 0, "end": 8000, "boundary_id": "raw_r20",
         "resolution": 20000, "n_source_boundaries": 1, "signal": "raw",
         "has_nonfinite_base_bins": False},
    ])
    leave_one_out = build_jitter_offsets(raw, balanced, resolution_metadata, base_bin_bp=4000)
    assert not leave_one_out.loc[
        (leave_one_out["source"] == "raw")
        & (leave_one_out["comparison"] == "full_vs_leave_r4000_out"),
        "eligible",
    ].iloc[0]


def test_unified_component_retains_terminal_and_bad_coverage_state():
    raw = boundary_frame("raw", [("A02", 0, 8000, "raw", "4000", 1)])
    balanced = boundary_frame("balanced", [("A02", 0, 8000, "balanced", "4000", 1)])

    terminal = unify_boundaries(
        raw.assign(terminal_partial_coordinate=True), balanced, tie_tolerance_bp=4000
    )
    unavailable = unify_boundaries(
        raw.assign(overlaps_bad_region=True), balanced, tie_tolerance_bp=4000
    )

    assert terminal.loc[0, "component_state"] == "simple_terminal_partial_raw"
    assert unavailable.loc[0, "component_state"] == "unevaluable_coverage_raw_bad_region"


def test_tad_unevaluable_mask_publishes_terminal_nonfinite_and_bad_regions():
    bins = pd.DataFrame({
        "chrom": ["A02", "A02", "A02"],
        "start": [0, 4000, 8000],
        "end": [4000, 8000, 12000],
        "weight": [np.nan, 1.0, 1.0],
    })

    observed = _tad_unevaluable(
        "B", {"A02": 12000}, bins, [("A02", 4000, 8000)], 4000
    )

    assert observed.to_dict("records") == [
        {"sample": "B", "chrom": "A02", "start": 0, "end": 4000, "reason": "nonfinite_base_bins"},
        {"sample": "B", "chrom": "A02", "start": 0, "end": 4000, "reason": "terminal_insufficient_tad_coverage"},
        {"sample": "B", "chrom": "A02", "start": 4000, "end": 8000, "reason": "bad_region"},
        {"sample": "B", "chrom": "A02", "start": 8000, "end": 12000, "reason": "terminal_insufficient_tad_coverage"},
    ]


def test_sparse_matcher_scales_with_disjoint_local_overlaps():
    count = 1500
    raw = boundary_frame("raw", [
        ("A02", index * 20000, index * 20000 + 8000, f"raw_{index}", "4000", 1)
        for index in range(count)
    ]).assign(center_bp=lambda frame: (frame.start + frame.end) // 2)
    balanced = boundary_frame("balanced", [
        ("A02", index * 20000 + 4000, index * 20000 + 12000, f"balanced_{index}", "4000", 1)
        for index in range(count)
    ]).assign(center_bp=lambda frame: (frame.start + frame.end) // 2)

    assert maximum_cardinality_pairs(raw, balanced, lambda left, right: left.start < right.end and right.start < left.end) == [
        (index, index) for index in range(count)
    ]
    offsets = build_jitter_offsets(
        raw.drop(columns="center_bp"), balanced.drop(columns="center_bp"),
        pd.DataFrame(columns=["chrom", "start", "end", "boundary_id", "resolution", "n_source_boundaries", "signal"]),
        base_bin_bp=4000,
    )
    assert len(offsets) == count


def make_mcool(path):
    bins = pd.DataFrame({
        "chrom": ["A02"] * 12,
        "start": np.arange(12) * 4000,
        "end": np.arange(1, 13) * 4000,
        "weight": np.ones(12),
    })
    pixels = pd.DataFrame({"bin1_id": np.arange(12), "bin2_id": np.arange(12), "count": 1})
    cooler.create_cooler(f"{path}::/resolutions/4000", bins, pixels, ordered=True)
    return path


def test_cli_writes_stable_tables_without_temporary_files(tmp_path):
    raw_multi = write_bed7(tmp_path / "raw.multires.bed", [
        ("A02", 0, 8000, "raw_1", "4000", 1, 1),
    ])
    balanced_multi = write_bed7(tmp_path / "balanced.multires.bed", [
        ("A02", 4000, 12000, "balanced_1", "4000", 1, 1),
    ])
    raw_paths, balanced_paths = [], []
    for resolution, start in ((4000, 0), (8000, 16000), (20000, 32000)):
        raw_paths.append(write_bed6(tmp_path / f"raw.{resolution}.bed", [
            ("A02", start, start + 8000, f"raw_{resolution}", resolution, 1),
        ]))
        balanced_paths.append(write_bed6(tmp_path / f"balanced.{resolution}.bed", [
            ("A02", start + 4000, start + 12000, f"balanced_{resolution}", resolution, 1),
        ]))
    mcool = make_mcool(tmp_path / "sample.mcool")
    output_unified = tmp_path / "unified.tsv"
    output_offsets = tmp_path / "offsets.tsv"
    output_views = tmp_path / "views.tsv"
    output_tad_unevaluable = tmp_path / "tad_unevaluable.tsv"

    result = subprocess.run(
        [
            sys.executable, str(SCRIPT), "--sample", "S", "--raw-multires", str(raw_multi),
            "--balanced-multires", str(balanced_multi), "--raw-resolution-boundaries",
            *map(str, raw_paths), "--balanced-resolution-boundaries", *map(str, balanced_paths),
            "--mcool", str(mcool), "--config", str(REPO / "config" / "config.yaml"),
            "--output-unified", str(output_unified), "--output-offsets", str(output_offsets),
            "--output-views", str(output_views),
            "--output-tad-unevaluable", str(output_tad_unevaluable),
            "--log", str(tmp_path / "unify.log"),
        ], text=True, capture_output=True, check=False,
    )

    assert result.returncode == 0, result.stderr
    unified = pd.read_csv(output_unified, sep="\t")
    views = pd.read_csv(output_views, sep="\t")
    assert unified["unified_boundary_id"].tolist() == ["unified_boundary_000001"]
    assert set(views["view"]) == {
        "unified", "raw_only_input", "balanced_only_input", "leave_r4000_out",
        "leave_r8000_out", "leave_r20000_out",
    }
    assert pd.read_csv(output_offsets, sep="\t").columns.tolist()[0] == "sample"
    assert pd.read_csv(output_tad_unevaluable, sep="\t").columns.tolist() == list(
        TAD_UNEVALUABLE_COLUMNS
    )
    assert not list(tmp_path.glob(".*.tmp"))
