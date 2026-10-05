from pathlib import Path
import sys

import pytest


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from tad_boundaries import (
    Interval,
    boundaries_from_tads,
    merge_intervals,
    write_multires_boundaries,
    write_resolution_boundaries,
)
from tad_errors import TadError


def interval(start, end, resolution, count=1):
    return Interval("A02", start, end, {resolution}, count)


def test_endpoint_expansion_clips_at_chromosome_edges():
    boundaries = boundaries_from_tads(
        [("A02", 4000, 12000)],
        {"A02": 16000},
        4000,
        sample="S",
        signal="raw",
    )
    assert boundaries == [
        interval(0, 8000, 4000),
        interval(8000, 16000, 4000),
    ]


def test_same_resolution_and_transitive_overlap_take_union():
    merged = merge_intervals(
        [
            interval(0, 8000, 4000),
            interval(4000, 12000, 4000),
            interval(10000, 20000, 8000),
        ],
        {"A02": 0},
    )
    assert merged == [Interval("A02", 0, 20000, {4000, 8000}, 3)]


def test_endpoint_adjacency_does_not_merge():
    assert merge_intervals(
        [interval(0, 8000, 4000), interval(8000, 16000, 8000)],
        {"A02": 0},
    ) == [interval(0, 8000, 4000), interval(8000, 16000, 8000)]


def test_invalid_tad_coordinate_has_dynamic_stable_error():
    with pytest.raises(TadError) as error:
        boundaries_from_tads(
            [("A02", -1, 4000)],
            {"A02": 16000},
            4000,
            sample="actual-sample",
            signal="balanced",
        )
    assert str(error.value).startswith("ERROR[TAD_INVALID_COORDINATE]")
    assert "sample=actual-sample" in str(error.value)
    assert "signal=balanced" in str(error.value)


def test_provenance_counts_survive_second_merge_pass():
    merged = merge_intervals(
        [
            Interval("A02", 0, 12000, {4000}, 2),
            Interval("A02", 10000, 20000, {8000, 20000}, 4),
        ],
        {"A02": 0},
    )
    assert merged == [Interval("A02", 0, 20000, {4000, 8000, 20000}, 6)]


def test_bed6_and_bed7_contracts_have_no_header(tmp_path):
    interval_row = Interval("A02", 8000, 20000, {4000, 8000}, 3)
    bed6 = tmp_path / "resolution.bed"
    bed7 = tmp_path / "multires.bed"
    write_resolution_boundaries(bed6, [interval_row], 4000)
    write_multires_boundaries(bed7, [interval_row])
    assert bed6.read_text() == "A02\t8000\t20000\tboundary_000001\t4000\t3\n"
    assert bed7.read_text() == (
        "A02\t8000\t20000\tboundary_000001\t4000,8000\t2\t3\n"
    )
