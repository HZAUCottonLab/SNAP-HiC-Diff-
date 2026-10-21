import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from is_boundary_point_state import (
    BOTH, GAINED, LOST, NONE, SHIFTED, STRENGTHENED, WEAKENED,
    BoundaryIndex, classify_point, locus_point, mark_shifts, point_rows,
)


def _boundaries(rows):
    return pd.DataFrame([
        {
            "chrom": chrom, "start": start, "end": end, "boundary_id": boundary_id,
            "center_bp": (start + end) // 2, "component_state": "simple",
            "signal_support": "both", "support_resolutions": "4000,8000",
            "n_supporting_resolutions": 2, "n_source_boundaries": 2,
        }
        for chrom, start, end, boundary_id in rows
    ])


def _locus(summit=100000, event_start=96000, event_end=104000,
           region_id=1, discovery_mode="refined", chrom="A02"):
    return {
        "region_id": region_id, "discovery_mode": discovery_mode, "chrom": chrom,
        "summit": summit, "event_start": event_start, "event_end": event_end,
    }


def test_point_is_the_summit():
    assert locus_point(_locus(summit=123456), 4000) == (123456, "summit")


def test_missing_summit_falls_back_to_event_start_and_says_so():
    assert locus_point(_locus(summit=float("nan"), event_start=77000), 4000) == (77000, "event_start")


def test_boundary_interval_is_half_open():
    index = BoundaryIndex(_boundaries([("A02", 99000, 101000, "b1")]))
    assert index.query("A02", 99000).n_hits == 1
    assert index.query("A02", 100999).n_hits == 1
    assert index.query("A02", 101000).n_hits == 0


def test_point_outside_every_boundary_reports_edge_distance():
    index = BoundaryIndex(_boundaries([("A02", 90000, 98000, "b1")]))
    query = index.query("A02", 102000)
    assert query.n_hits == 0
    assert query.edge_distance == 4000


def test_locus_on_another_chromosome_never_hits():
    index = BoundaryIndex(_boundaries([("A03", 90000, 110000, "b1")]))
    assert index.query("A02", 100000).n_hits == 0


def test_two_overlapping_resolutions_both_count():
    index = BoundaryIndex(_boundaries([
        ("A02", 96000, 112000, "coarse"),
        ("A02", 99000, 105000, "fine"),
    ]))
    assert index.query("A02", 100000).n_hits == 2


def test_rows_keep_coordinates_for_every_containing_boundary():
    index = BoundaryIndex(_boundaries([
        ("A02", 96000, 112000, "coarse"),
        ("A02", 99000, 105000, "fine"),
    ]))
    rows, n_hits, edge, point, source = point_rows(
        _locus(), index, sample="H1398", comparison="A_vs_B", base_bin_bp=4000
    )
    assert n_hits == 2
    assert [row["boundary_id"] for row in rows] == ["coarse", "fine"]
    assert all(row["inside"] for row in rows)
    assert all(row["center_offset_bp"] == point - row["center_bp"] for row in rows)


def test_locus_off_boundary_still_publishes_the_nearest_one():
    index = BoundaryIndex(_boundaries([("A02", 60000, 68000, "far")]))
    rows, n_hits, edge, _, _ = point_rows(
        _locus(), index, sample="H1398", comparison="A_vs_B", base_bin_bp=4000
    )
    assert n_hits == 0
    assert len(rows) == 1
    assert rows[0]["inside"] is False
    assert rows[0]["edge_distance_bp"] == edge > 0


def test_locus_with_no_boundary_on_its_chromosome_publishes_nothing():
    index = BoundaryIndex(_boundaries([("A03", 60000, 68000, "far")]))
    rows, n_hits, _, _, _ = point_rows(
        _locus(), index, sample="H1398", comparison="A_vs_B", base_bin_bp=4000
    )
    assert rows == []
    assert n_hits == 0


def test_state_matrix():
    assert classify_point(True, False, "") == LOST
    assert classify_point(False, True, "") == GAINED
    assert classify_point(False, False, "") == NONE
    assert classify_point(True, True, "stronger_in_B") == STRENGTHENED
    assert classify_point(True, True, "weaker_in_B") == WEAKENED
    assert classify_point(True, True, "") == BOTH


def test_shift_pairs_only_valleys_within_tolerance():
    states = {
        1: {"state": LOST, "chrom": "A02", "point": 100000, "direction": "valley"},
        2: {"state": GAINED, "chrom": "A02", "point": 110000, "direction": "valley"},
        3: {"state": LOST, "chrom": "A02", "point": 300000, "direction": "valley"},
        4: {"state": GAINED, "chrom": "A02", "point": 340000, "direction": "valley"},
    }
    partners = mark_shifts(states, tolerance_bp=12000)
    assert partners == {1: 2, 2: 1}


def test_shift_ignores_peaks():
    states = {
        1: {"state": LOST, "chrom": "A02", "point": 100000, "direction": "peak"},
        2: {"state": GAINED, "chrom": "A02", "point": 105000, "direction": "peak"},
    }
    assert mark_shifts(states, tolerance_bp=12000) == {}


def test_shift_never_crosses_chromosomes():
    states = {
        1: {"state": LOST, "chrom": "A02", "point": 100000, "direction": "valley"},
        2: {"state": GAINED, "chrom": "A03", "point": 101000, "direction": "valley"},
    }
    assert mark_shifts(states, tolerance_bp=12000) == {}
