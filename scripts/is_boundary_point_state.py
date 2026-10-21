"""Point-in-boundary classification for IS-change loci.

An IS-change locus is a point. A boundary is an interval whose width follows
the resolution that called it, so asking which boundary centre best matches
which produces artificial ties. The question that matters is simpler: does the
locus point fall inside a boundary at all, in each sample? A locus is at a
boundary when any resolution, any signal, says so.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

POINT_BOUNDARY_COLUMNS = (
    "comparison", "locus_key", "region_id", "discovery_mode", "chrom", "locus_summit",
    "point_source", "sample", "boundary_id", "start", "end", "center_bp",
    "inside", "center_offset_bp", "edge_distance_bp",
    "signal_support", "support_resolutions", "n_supporting_resolutions",
    "component_state",
)
POINT_STATE_COLUMNS = (
    "point_bp", "point_source",
    "point_at_boundary_a", "point_at_boundary_b",
    "point_n_boundaries_a", "point_n_boundaries_b",
    "point_edge_distance_a", "point_edge_distance_b",
    "point_state", "point_shift_partner",
)
STRENGTHENED = "strengthened_in_B"
WEAKENED = "weakened_in_B"
LOST = "lost_in_B"
GAINED = "gained_in_B"
SHIFTED = "shifted_in_B"
BOTH = "boundary_in_both"
NONE = "not_at_boundary"


def locus_point(locus, base_bin_bp):
    """Return ``(point, source)`` for a locus.

    ``summit`` is the refined event centre and is the only point a locus should
    contribute. An unrefined locus has no summit, so its event start stands in
    and the source is recorded so those rows stay identifiable.
    """
    summit = locus.get("summit") if isinstance(locus, dict) else getattr(locus, "summit", None)
    if summit is not None and not pd.isna(summit):
        return int(summit), "summit"
    return int(locus["event_start"] if isinstance(locus, dict) else locus.event_start), "event_start"


@dataclass
class BoundaryQuery:
    """Result of asking one boundary table about one point."""

    rows: pd.DataFrame
    n_hits: int
    edge_distance: int


class BoundaryIndex:
    """Chromosome-partitioned boundary intervals for point lookups."""

    def __init__(self, frame):
        self._by_chrom = {}
        if frame.empty:
            return
        # The unified table names its identifier ``unified_boundary_id``; the
        # per-resolution tables use ``boundary_id``. Accept either.
        if "boundary_id" not in frame.columns and "unified_boundary_id" in frame.columns:
            frame = frame.rename(columns={"unified_boundary_id": "boundary_id"})
        for chrom, group in frame.groupby("chrom", sort=False, observed=True):
            ordered = group.sort_values(["start", "end"], kind="stable").reset_index(drop=True)
            self._by_chrom[chrom] = ordered

    def has_chromosome(self, chrom):
        return chrom in self._by_chrom

    def query(self, chrom, point):
        frame = self._by_chrom.get(chrom)
        if frame is None or frame.empty:
            return BoundaryQuery(pd.DataFrame(columns=frame.columns if frame is not None else []), 0, -1)
        starts = frame["start"].to_numpy(dtype="int64")
        ends = frame["end"].to_numpy(dtype="int64")
        # Half-open intervals throughout: [start, end).
        inside = (starts <= point) & (point < ends)
        if inside.any():
            return BoundaryQuery(frame.loc[inside].reset_index(drop=True), int(inside.sum()), 0)
        # Distance to the nearest interval edge. Adjacent coordinates are 0:
        # both loci and boundary edges sit on the same 4 kb grid, so an
        # abutting boundary is adjacent, not 1 bp away.
        before = starts - point
        after = point - ends
        distances = np.where(before > 0, before, np.where(after > 0, after, 0))
        nearest = int(distances.min())
        return BoundaryQuery(frame.iloc[0:0].reset_index(drop=True), 0, nearest)


def classify_point(at_a, at_b, insulation_direction):
    """States derived purely from boundary presence at the locus point.

    When both samples carry a boundary the locus is not gained or lost; what
    changed is its strength, and ``insulation_direction`` says which way. That
    column comes from the IS-change caller: a valley in B means B's insulation
    score dropped, which is a stronger boundary.
    """
    if at_a and at_b:
        if insulation_direction == "stronger_in_B":
            return STRENGTHENED
        if insulation_direction == "weaker_in_B":
            return WEAKENED
        return BOTH
    if at_a:
        return LOST
    if at_b:
        return GAINED
    return NONE


def point_rows(locus, boundary_index, sample, comparison, base_bin_bp):
    """Long-table rows for one locus and one sample.

    Every boundary that contains the point is emitted, so overlapping
    resolutions stay visible. When none contains it, the single nearest
    boundary is emitted instead with ``inside`` false, which keeps every locus
    joinable to a breakpoint table instead of silently dropping it.
    """
    point, source = locus_point(locus, base_bin_bp)
    result = boundary_index.query(locus["chrom"], point)
    frame = result.rows
    if frame.empty:
        nearest = _nearest_row(boundary_index, locus["chrom"], point)
        if nearest is None:
            return [], 0, result.edge_distance, point, source
        frame = nearest.to_frame().T
    records = []
    for item in frame.itertuples(index=False):
        start, end = int(item.start), int(item.end)
        center = int(item.center_bp) if not pd.isna(item.center_bp) else (start + end) // 2
        records.append({
            "comparison": comparison,
            "locus_key": _locus_key(locus),
            "region_id": locus["region_id"],
            "discovery_mode": locus["discovery_mode"],
            "chrom": locus["chrom"],
            "locus_summit": point,
            "point_source": source,
            "sample": sample,
            "boundary_id": item.boundary_id,
            "start": start,
            "end": end,
            "center_bp": center,
            "inside": bool(start <= point < end),
            "center_offset_bp": point - center,
            "edge_distance_bp": 0 if start <= point < end else min(point - start, point - end),
            "signal_support": item.signal_support,
            "support_resolutions": item.support_resolutions,
            "n_supporting_resolutions": int(item.n_supporting_resolutions),
            "component_state": item.component_state,
        })
    return records, result.n_hits, result.edge_distance, point, source


def _nearest_row(boundary_index, chrom, point):
    frame = boundary_index._by_chrom.get(chrom)
    if frame is None or frame.empty:
        return None
    starts = frame["start"].to_numpy(dtype="int64")
    ends = frame["end"].to_numpy(dtype="int64")
    before = starts - point
    after = point - ends
    distances = np.where(before > 0, before, np.where(after > 0, after, 0))
    return frame.iloc[int(distances.argmin())]


def _locus_key(locus):
    """Stable per-row key.

    ``region_id`` repeats across loci and is not unique on its own, so the key
    combines it with the event interval, which is unique within a comparison.
    """
    if isinstance(locus, dict):
        values = (locus["region_id"], locus["event_start"], locus["event_end"],
                  locus["discovery_mode"])
    else:
        values = (locus.region_id, locus.event_start, locus.event_end,
                  locus.discovery_mode)
    return "{}:{}:{}:{}".format(*values)


def _locus_keys(frame):
    return [_locus_key(row._asdict()) for row in frame.itertuples(index=False)]


def _blank_point():
    return {
        "point_bp": float("nan"),
        "point_source": "",
        "point_at_boundary_a": False,
        "point_at_boundary_b": False,
        "point_n_boundaries_a": 0,
        "point_n_boundaries_b": 0,
        "point_edge_distance_a": -1,
        "point_edge_distance_b": -1,
        "point_state": "",
        "point_shift_partner": "",
    }


def attach_point_state(loci, boundaries_a, boundaries_b, sample_a, sample_b,
                       comparison, base_bin_bp, shift_tolerance_bp):
    """Attach the point-in-boundary state to annotated loci.

    Returns ``(annotated, long_rows)``. ``annotated`` gains every column in
    ``POINT_STATE_COLUMNS``; ``long_rows`` carries one row per locus, sample,
    and contributing boundary so each call stays joinable to a breakpoint
    table.
    """
    index_a = boundaries_a if isinstance(boundaries_a, BoundaryIndex) else BoundaryIndex(boundaries_a)
    index_b = boundaries_b if isinstance(boundaries_b, BoundaryIndex) else BoundaryIndex(boundaries_b)
    frame = loci.copy()
    long_rows = []
    states = {}
    rows = []
    for position, locus in enumerate(frame.to_dict("records")):
        records_a, hits_a, edge_a, point, source = point_rows(
            locus, index_a, sample_a, comparison, base_bin_bp
        )
        records_b, hits_b, edge_b, _, _ = point_rows(
            locus, index_b, sample_b, comparison, base_bin_bp
        )
        long_rows.extend(records_a)
        long_rows.extend(records_b)
        state = classify_point(hits_a > 0, hits_b > 0, locus.get("insulation_direction", ""))
        states[position] = {
            "state": state, "chrom": locus["chrom"], "point": point,
            "direction": locus.get("refinement_direction", ""),
        }
        rows.append({
            "point_bp": point,
            "point_source": source,
            "point_at_boundary_a": hits_a > 0,
            "point_at_boundary_b": hits_b > 0,
            "point_n_boundaries_a": hits_a,
            "point_n_boundaries_b": hits_b,
            "point_edge_distance_a": edge_a,
            "point_edge_distance_b": edge_b,
            "point_state": state,
            "point_shift_partner": "",
        })
    for column in POINT_STATE_COLUMNS:
        frame[column] = [row[column] for row in rows]
    partners = mark_shifts(states, shift_tolerance_bp)
    for position, partner in partners.items():
        frame.at[position, "point_state"] = SHIFTED
        partner_point, _ = locus_point(frame.iloc[partner].to_dict(), base_bin_bp)
        frame.at[position, "point_shift_partner"] = str(partner_point)
    long_frame = pd.DataFrame(long_rows, columns=POINT_BOUNDARY_COLUMNS)
    return frame, long_frame


def mark_shifts(states, tolerance_bp):
    """Relabel a lost/gained pair that sits within ``tolerance_bp`` as a shift.

    Both members must be valleys: a boundary that disappeared here and
    reappeared there is one boundary that moved. The returned mapping is
    locus index to partner locus index, empty when nothing moved.
    """
    lost = [index for index, item in states.items() if item["state"] == LOST]
    gained = [index for index, item in states.items() if item["state"] == GAINED]
    partners = {}
    for left in lost:
        candidate = None
        for right in gained:
            if right in partners.values():
                continue
            if states[left]["chrom"] != states[right]["chrom"]:
                continue
            if states[left]["direction"] != "valley" or states[right]["direction"] != "valley":
                continue
            distance = abs(states[left]["point"] - states[right]["point"])
            if distance > tolerance_bp:
                continue
            if candidate is None or distance < candidate[0]:
                candidate = (distance, right)
        if candidate is not None:
            partners[left] = candidate[1]
            partners[candidate[1]] = left
    return partners
