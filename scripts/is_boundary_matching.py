"""Calibrate technical jitter and classify ordered cross-sample boundary events."""

from dataclasses import dataclass
import heapq
import math

import numpy as np
import pandas as pd

from is_boundary_io import BOUNDARY_EVENT_COLUMNS, BOUNDARY_STABILITY_COLUMNS
from is_boundary_unification import _overlap_pairs, maximum_cardinality_pairs, overlaps


@dataclass(frozen=True)
class JitterEstimate:
    threshold_bp: int
    source: str
    n_offsets: int
    raw_quantile_bp: float


def estimate_jitter(offsets_a, offsets_b, settings) -> JitterEstimate:
    values = pd.concat(
        [
            offsets_a.loc[offsets_a["eligible"], "offset_bp"],
            offsets_b.loc[offsets_b["eligible"], "offset_bp"],
        ],
        ignore_index=True,
    ).dropna().astype(float)
    if len(values) < settings.jitter_min_pairs:
        return JitterEstimate(
            settings.jitter_fallback_bp,
            "fallback_max_tad_resolution",
            len(values),
            float("nan"),
        )
    raw = float(np.quantile(values, settings.jitter_quantile))
    aligned = int(math.ceil(raw / settings.base_bin_bp) * settings.base_bin_bp)
    return JitterEstimate(
        max(settings.base_bin_bp, aligned), "empirical_q95", len(values), raw
    )


def align_ordered(boundaries_a, boundaries_b, max_distance_bp: int):
    """Return deterministic noncrossing index pairs within ``max_distance_bp``."""
    m = len(boundaries_b)
    nodes = [(0, 0, None, None)]
    tree = [0] * (m + 1)

    def better(left, right):
        left_node, right_node = nodes[left], nodes[right]
        return left if (left_node[0], left_node[1], left) > (
            right_node[0], right_node[1], right
        ) else right

    def query(index_b):
        best = 0
        while index_b:
            best = better(best, tree[index_b])
            index_b -= index_b & -index_b
        return best

    def update(index_b, node):
        while index_b <= m:
            tree[index_b] = better(tree[index_b], node)
            index_b += index_b & -index_b

    pairs = _local_pairs(boundaries_a, boundaries_b, max_distance_bp)
    start = 0
    while start < len(pairs):
        index_a = pairs[start][0]
        end = start
        pending = []
        while end < len(pairs) and pairs[end][0] == index_a:
            _, index_b, distance = pairs[end]
            previous = query(index_b)
            score, count, _, _ = nodes[previous]
            nodes.append((score + max_distance_bp - distance, count + 1, previous, (index_a, index_b)))
            pending.append((index_b + 1, len(nodes) - 1))
            end += 1
        for index_b, node in pending:
            update(index_b, node)
        start = end

    node = query(m)
    matches = []
    while node:
        _, _, node, pair = nodes[node]
        matches.append(pair)
    return list(reversed(matches))


def _sorted(frame):
    return frame.sort_values(
        ["center_bp", "start", "end", "unified_boundary_id"], kind="stable"
    ).reset_index(drop=True)


def _local_pairs(boundaries_a, boundaries_b, max_distance_bp):
    centers_a = boundaries_a["center_bp"].to_numpy(dtype=int)
    centers_b = boundaries_b["center_bp"].to_numpy(dtype=int)
    pairs = []
    for index_a, center_a in enumerate(centers_a):
        start = np.searchsorted(centers_b, center_a - max_distance_bp, side="left")
        end = np.searchsorted(centers_b, center_a + max_distance_bp, side="right")
        pairs.extend(
            (index_a, index_b, abs(int(center_a) - int(centers_b[index_b])))
            for index_b in range(start, end)
        )
    return pairs


def _candidate_data(boundaries_a, boundaries_b, max_distance_bp, *, include_overlaps=False):
    by_a = {index: [] for index in range(len(boundaries_a))}
    by_b = {index: [] for index in range(len(boundaries_b))}
    edges = _local_pairs(boundaries_a, boundaries_b, max_distance_bp)
    if include_overlaps:
        # Fixed-overlap matches are independent of the shift radius. Wide
        # boundary intervals can overlap while their centers exceed that
        # radius; keep their metadata available during sensitivity checks.
        existing = {(index_a, index_b) for index_a, index_b, _ in edges}
        centers_a = boundaries_a["center_bp"].to_numpy(dtype=int)
        centers_b = boundaries_b["center_bp"].to_numpy(dtype=int)
        for index_a, index_b in _overlap_pairs(boundaries_a, boundaries_b):
            if (index_a, index_b) not in existing:
                edges.append((index_a, index_b, abs(int(centers_a[index_a]) - int(centers_b[index_b]))))
    for index_a, index_b, distance in edges:
        by_a[index_a].append((distance, index_b))
        by_b[index_b].append((distance, index_a))
    for values in (*by_a.values(), *by_b.values()):
        values.sort()

    components = []
    seen_a, seen_b = set(), set()
    for start_a in by_a:
        if start_a in seen_a or not by_a[start_a]:
            continue
        indices_a, indices_b = set(), set()
        pending = [("a", start_a)]
        while pending:
            side, index = pending.pop()
            if side == "a":
                if index in seen_a:
                    continue
                seen_a.add(index)
                indices_a.add(index)
                pending.extend(("b", other) for _, other in by_a[index])
            else:
                if index in seen_b:
                    continue
                seen_b.add(index)
                indices_b.add(index)
                pending.extend(("a", other) for _, other in by_b[index])
        component_edges = [
            (index_a, index_b, distance)
            for index_a in indices_a
            for distance, index_b in by_a[index_a]
        ]
        components.append((indices_a, indices_b, component_edges))
    component_by_a = {}
    component_by_b = {}
    for component in components:
        for index in component[0]:
            component_by_a[index] = component
        for index in component[1]:
            component_by_b[index] = component
    return by_a, by_b, components, component_by_a, component_by_b


def _match_info(index_a, index_b, candidate_data, tie_tolerance_bp):
    by_a, by_b, _, component_by_a, _ = candidate_data
    distances_a = by_a[index_a]
    distances_b = by_b[index_b]
    first_a, first_b = distances_a[0][0], distances_b[0][0]
    reciprocal = distances_a[0][1] == index_b and distances_b[0][1] == index_a
    second_a = distances_a[1][0] if len(distances_a) > 1 else float("inf")
    second_b = distances_b[1][0] if len(distances_b) > 1 else float("inf")
    margin = min(second_a - first_a, second_b - first_b)
    component = component_by_a[index_a]
    return {
        "reciprocal": reciprocal,
        "margin": pd.NA if math.isinf(margin) else int(margin),
        "is_tied": margin < tie_tolerance_bp,
        "n_local_candidates": len(component[2]),
    }


def _reason(info):
    if info["is_tied"]:
        return "match_tie_within_tolerance"
    return "nonreciprocal_closest_match"


def _remaining_frame(boundaries, indices):
    frame = boundaries.iloc[sorted(indices)].copy()
    frame["_original_index"] = frame.index
    return frame.reset_index(drop=True)


def _missing_side_unevaluable_reason(mask, chrom, start, end, sample):
    """Return the audited reason when absence at an endpoint cannot be called."""
    if mask is None or mask.empty:
        return ""
    overlapping = mask.loc[
        mask["chrom"].eq(chrom)
        & (mask["start"] < end)
        & (mask["end"] > start),
        "reason",
    ]
    reasons = sorted({str(reason) for reason in overlapping if str(reason)})
    return f"{sample}:{','.join(reasons)}" if reasons else ""


def _support(rows, column, separator=","):
    values = [str(item[column]) for item in rows if pd.notna(item[column])]
    if column == "support_resolutions":
        return ",".join(map(str, sorted({
            int(value) for item in values for value in item.split(",") if value
        })))
    return separator.join(dict.fromkeys(values))


def _event_record(
    rows_a,
    rows_b,
    *,
    comparison,
    sample_a,
    sample_b,
    position_state,
    position_substate="",
    match_stage="",
    match_info=None,
    jitter=None,
    max_shift_bp=None,
    ambiguity_reason="",
):
    rows_a, rows_b = list(rows_a), list(rows_b)
    coordinates = [
        coordinate
        for row in (*rows_a, *rows_b)
        for coordinate in (int(row.start), int(row.end))
    ]
    record = {
        "comparison": comparison,
        "sample_a": sample_a,
        "sample_b": sample_b,
        "chrom": (rows_a or rows_b)[0].chrom,
        "boundary_event_id": "",
        "position_state": position_state,
        "position_substate": position_substate,
        "boundary_a_id": _support(rows_a, "unified_boundary_id"),
        "boundary_b_id": _support(rows_b, "unified_boundary_id"),
        "a_start": min(int(row.start) for row in rows_a) if rows_a else pd.NA,
        "a_end": max(int(row.end) for row in rows_a) if rows_a else pd.NA,
        "b_start": min(int(row.start) for row in rows_b) if rows_b else pd.NA,
        "b_end": max(int(row.end) for row in rows_b) if rows_b else pd.NA,
        "event_span_start": min(coordinates),
        "event_span_end": max(coordinates),
        "shift_bp": (
            int(rows_b[0].center_bp) - int(rows_a[0].center_bp)
            if len(rows_a) == len(rows_b) == 1 else pd.NA
        ),
        "match_stage": match_stage,
        "match_is_reciprocal": (
            match_info["reciprocal"] if match_info is not None else pd.NA
        ),
        "second_best_margin_bp": (
            match_info["margin"] if match_info is not None else pd.NA
        ),
        "n_local_candidates": (
            match_info["n_local_candidates"] if match_info is not None else 0
        ),
        "jitter_threshold_bp": jitter.threshold_bp,
        "jitter_source": jitter.source,
        "max_shift_bp": max_shift_bp,
        "ambiguity_reason": ambiguity_reason,
        "state_at_60000": "",
        "state_at_80000": "",
        "state_at_100000": "",
        "position_state_stable": False,
        "signal_support_a": _support(rows_a, "signal_support"),
        "signal_support_b": _support(rows_b, "signal_support"),
        "resolution_support_a": _support(rows_a, "support_resolutions"),
        "resolution_support_b": _support(rows_b, "support_resolutions"),
    }
    return record


def _interval_associations(reference, candidates, start_column, end_column, associations):
    """Add overlaps for one endpoint role with an ordered interval sweep."""
    valid_reference = reference[start_column].notna() & reference[end_column].notna()
    valid_candidates = candidates[start_column].notna() & candidates[end_column].notna()
    shared_chromosomes = sorted(
        set(reference.loc[valid_reference, "chrom"])
        & set(candidates.loc[valid_candidates, "chrom"])
    )
    for chrom in shared_chromosomes:
        reference_mask = valid_reference & reference["chrom"].eq(chrom)
        candidate_mask = valid_candidates & candidates["chrom"].eq(chrom)
        reference_positions = np.flatnonzero(reference_mask.to_numpy())
        candidate_positions = np.flatnonzero(candidate_mask.to_numpy())
        left = reference.iloc[reference_positions].assign(
            _position=reference_positions
        ).sort_values([start_column, end_column], kind="stable")
        right = candidates.iloc[candidate_positions].assign(
            _position=candidate_positions
        ).sort_values([start_column, end_column], kind="stable")
        right_rows = list(right[[start_column, end_column, "_position"]].itertuples(index=False))
        active, next_right = [], 0
        for item in left[[start_column, end_column, "_position"]].itertuples(index=False):
            left_start, left_end, left_position = map(int, item)
            while next_right < len(right_rows) and int(right_rows[next_right][0]) < left_end:
                right_start, right_end, right_position = right_rows[next_right]
                heapq.heappush(active, (int(right_end), int(right_position)))
                next_right += 1
            while active and active[0][0] <= left_start:
                heapq.heappop(active)
            associations[left_position].update(position for _, position in active)


def _association_indices(reference, candidates):
    associations = [set() for _ in range(len(reference))]
    if reference.empty or candidates.empty:
        return associations
    for start_column, end_column in (("a_start", "a_end"), ("b_start", "b_end")):
        _interval_associations(
            reference, candidates, start_column, end_column, associations
        )
    return associations


def _association_summaries(reference, candidates):
    """Summarize endpoint-associated events without per-reference table scans."""
    states = candidates["position_state"].astype(str).to_numpy()
    rows = []
    for indices in _association_indices(reference, candidates):
        if not indices:
            rows.append({"n_associated": 0, "position_state": "no_boundary_in_view"})
            continue
        associated_states = [states[index] for index in sorted(indices)]
        rows.append({
            "n_associated": len(indices),
            "position_state": (
                "complex"
                if len(indices) > 1 or "complex" in associated_states
                else associated_states[0]
            ),
        })
    return pd.DataFrame(rows, columns=("n_associated", "position_state"))


def _associated_events(reference, candidates):
    single_reference = pd.DataFrame([{
        "chrom": reference.chrom,
        "a_start": reference.a_start,
        "a_end": reference.a_end,
        "b_start": reference.b_start,
        "b_end": reference.b_end,
    }])
    indices = _association_indices(single_reference, candidates)[0]
    return candidates.iloc[sorted(indices)]


def _associated_state(reference, candidates):
    associated = _associated_events(reference, candidates)
    if associated.empty:
        return "no_boundary_in_view"
    if len(associated) > 1 or associated["position_state"].eq("complex").any():
        return "complex"
    return associated.iloc[0].position_state


def _events_at_distance(
    boundaries_a,
    boundaries_b,
    jitter,
    settings,
    comparison,
    sample_a,
    sample_b,
    max_shift_bp,
    tad_unevaluable_a=None,
    tad_unevaluable_b=None,
):
    records = []
    chromosomes = tuple(sorted(set(boundaries_a["chrom"]) | set(boundaries_b["chrom"])))
    for chrom in chromosomes:
        all_a = _sorted(boundaries_a[boundaries_a["chrom"] == chrom])
        all_b = _sorted(boundaries_b[boundaries_b["chrom"] == chrom])
        simple_a = _sorted(all_a[all_a["component_state"].str.startswith("simple")])
        simple_b = _sorted(all_b[all_b["component_state"].str.startswith("simple")])
        for _, row in all_a[~all_a["component_state"].str.startswith("simple")].iterrows():
            unevaluable = row.component_state.startswith("unevaluable_coverage")
            records.append(_event_record(
                [row], [], comparison=comparison, sample_a=sample_a, sample_b=sample_b,
                position_state=("ambiguous" if unevaluable else "complex"),
                match_stage=("coverage_unevaluable" if unevaluable else "within_sample_complex"),
                jitter=jitter, max_shift_bp=max_shift_bp,
                ambiguity_reason=(
                    "unevaluable_boundary_coverage:" + row.component_state.removeprefix("unevaluable_coverage_")
                    if unevaluable else "complex_within_sample"
                ),
            ))
        for _, row in all_b[~all_b["component_state"].str.startswith("simple")].iterrows():
            unevaluable = row.component_state.startswith("unevaluable_coverage")
            records.append(_event_record(
                [], [row], comparison=comparison, sample_a=sample_a, sample_b=sample_b,
                position_state=("ambiguous" if unevaluable else "complex"),
                match_stage=("coverage_unevaluable" if unevaluable else "within_sample_complex"),
                jitter=jitter, max_shift_bp=max_shift_bp,
                ambiguity_reason=(
                    "unevaluable_boundary_coverage:" + row.component_state.removeprefix("unevaluable_coverage_")
                    if unevaluable else "complex_within_sample"
                ),
            ))

        candidates = _candidate_data(simple_a, simple_b, max_shift_bp, include_overlaps=True)
        consumed_a, consumed_b = set(), set()
        for component_a, component_b, edges in candidates[2]:
            if len(component_a) > 1 and len(component_b) > 1 and any(
                _match_info(index_a, index_b, candidates, settings.match_tie_tolerance_bp)["is_tied"]
                for index_a, index_b, _ in edges
            ):
                records.append(_event_record(
                    [simple_a.iloc[index] for index in sorted(component_a)],
                    [simple_b.iloc[index] for index in sorted(component_b)],
                    comparison=comparison, sample_a=sample_a, sample_b=sample_b,
                    position_state="complex", match_stage="local_tied_component",
                    jitter=jitter, max_shift_bp=max_shift_bp,
                    ambiguity_reason="multiple_tied_non_crossing_matches",
                    match_info={"reciprocal": False, "margin": 0, "n_local_candidates": len(edges)},
                ))
                consumed_a.update(component_a)
                consumed_b.update(component_b)

        remaining_a = set(range(len(simple_a))).difference(consumed_a)
        remaining_b = set(range(len(simple_b))).difference(consumed_b)
        exact_a = simple_a.iloc[sorted(remaining_a)]
        exact_b = simple_b.iloc[sorted(remaining_b)]
        for local_a, local_b in maximum_cardinality_pairs(exact_a, exact_b, overlaps):
            index_a, index_b = int(exact_a.index[local_a]), int(exact_b.index[local_b])
            info = _match_info(index_a, index_b, candidates, settings.match_tie_tolerance_bp)
            state = "conserved" if info["reciprocal"] and not info["is_tied"] else "ambiguous"
            records.append(_event_record(
                [simple_a.iloc[index_a]], [simple_b.iloc[index_b]],
                comparison=comparison, sample_a=sample_a, sample_b=sample_b,
                position_state=state,
                position_substate="conserved_exact" if state == "conserved" else "",
                match_stage="fixed_overlap", match_info=info, jitter=jitter,
                max_shift_bp=max_shift_bp,
                ambiguity_reason="" if state == "conserved" else _reason(info),
            ))
            remaining_a.discard(index_a)
            remaining_b.discard(index_b)

        jitter_a = _remaining_frame(simple_a, remaining_a)
        jitter_b = _remaining_frame(simple_b, remaining_b)
        jitter_candidates = _candidate_data(jitter_a, jitter_b, jitter.threshold_bp)
        jitter_pairs = []
        for index_a, values in jitter_candidates[0].items():
            if not values:
                continue
            _, index_b = values[0]
            info = _match_info(
                index_a, index_b, jitter_candidates, settings.match_tie_tolerance_bp
            )
            if info["reciprocal"] and not info["is_tied"] and not overlaps(
                jitter_a.iloc[index_a], jitter_b.iloc[index_b]
            ):
                jitter_pairs.append((index_a, index_b))
        for local_a, local_b in sorted(set(jitter_pairs)):
            index_a = int(jitter_a.iloc[local_a]._original_index)
            index_b = int(jitter_b.iloc[local_b]._original_index)
            info = _match_info(
                local_a, local_b, jitter_candidates, settings.match_tie_tolerance_bp
            )
            records.append(_event_record(
                [simple_a.iloc[index_a]], [simple_b.iloc[index_b]],
                comparison=comparison, sample_a=sample_a, sample_b=sample_b,
                position_state="conserved",
                position_substate="conserved_within_jitter",
                match_stage="fixed_jitter", match_info=info, jitter=jitter,
                max_shift_bp=max_shift_bp,
            ))
            remaining_a.discard(index_a)
            remaining_b.discard(index_b)

        shifted_a = _remaining_frame(simple_a, remaining_a)
        shifted_b = _remaining_frame(simple_b, remaining_b)
        residual_candidates = _candidate_data(shifted_a, shifted_b, max_shift_bp)
        residual_a = dict(zip(shifted_a["_original_index"], range(len(shifted_a))))
        residual_b = dict(zip(shifted_b["_original_index"], range(len(shifted_b))))
        for local_a, local_b in align_ordered(shifted_a, shifted_b, max_shift_bp):
            index_a = int(shifted_a.iloc[local_a]._original_index)
            index_b = int(shifted_b.iloc[local_b]._original_index)
            item_a, item_b = simple_a.iloc[index_a], simple_b.iloc[index_b]
            distance = abs(int(item_a.center_bp) - int(item_b.center_bp))
            info = _match_info(
                local_a, local_b, residual_candidates, settings.match_tie_tolerance_bp
            )
            if jitter.threshold_bp >= max_shift_bp:
                state, reason = "ambiguous", "jitter_not_smaller_than_shift_radius"
            elif distance <= jitter.threshold_bp:
                state, reason = "ambiguous", "within_jitter_not_unique"
            elif not info["reciprocal"] or info["is_tied"]:
                state, reason = "ambiguous", _reason(info)
            else:
                state, reason = "shifted", ""
            records.append(_event_record(
                [item_a], [item_b], comparison=comparison, sample_a=sample_a,
                sample_b=sample_b, position_state=state, match_stage="ordered_dp",
                match_info=info, jitter=jitter, max_shift_bp=max_shift_bp,
                ambiguity_reason=reason,
            ))
            remaining_a.discard(index_a)
            remaining_b.discard(index_b)

        for index_a in sorted(remaining_a):
            info = {
                "reciprocal": False,
                "margin": pd.NA,
                "n_local_candidates": len(residual_candidates[0][residual_a[index_a]]),
            }
            item_a = simple_a.iloc[index_a]
            unevaluable = _missing_side_unevaluable_reason(
                tad_unevaluable_b, chrom, int(item_a.start), int(item_a.end), sample_b
            )
            records.append(_event_record(
                [item_a], [], comparison=comparison, sample_a=sample_a,
                sample_b=sample_b,
                position_state="ambiguous" if unevaluable else "lost_in_B",
                match_stage="missing_side_unevaluable" if unevaluable else "unmatched",
                match_info=info, jitter=jitter, max_shift_bp=max_shift_bp,
                ambiguity_reason=(f"missing_side_unevaluable:{unevaluable}" if unevaluable else ""),
            ))
        for index_b in sorted(remaining_b):
            info = {
                "reciprocal": False,
                "margin": pd.NA,
                "n_local_candidates": len(residual_candidates[1][residual_b[index_b]]),
            }
            item_b = simple_b.iloc[index_b]
            unevaluable = _missing_side_unevaluable_reason(
                tad_unevaluable_a, chrom, int(item_b.start), int(item_b.end), sample_a
            )
            records.append(_event_record(
                [], [item_b], comparison=comparison, sample_a=sample_a,
                sample_b=sample_b,
                position_state="ambiguous" if unevaluable else "gained_in_B",
                match_stage="missing_side_unevaluable" if unevaluable else "unmatched",
                match_info=info, jitter=jitter, max_shift_bp=max_shift_bp,
                ambiguity_reason=(f"missing_side_unevaluable:{unevaluable}" if unevaluable else ""),
            ))

    events = pd.DataFrame(records, columns=BOUNDARY_EVENT_COLUMNS)
    if events.empty:
        return events
    ranks = {chrom: index for index, chrom in enumerate(chromosomes)}
    events = events.assign(_chrom_rank=events["chrom"].map(ranks)).sort_values(
        ["_chrom_rank", "event_span_start", "event_span_end", "position_state"], kind="stable"
    ).drop(columns="_chrom_rank").reset_index(drop=True)
    events["boundary_event_id"] = [
        f"boundary_event_{index:06d}" for index in range(1, len(events) + 1)
    ]
    return events.loc[:, BOUNDARY_EVENT_COLUMNS]


def call_boundary_events(
    boundaries_a, boundaries_b, jitter, settings, comparison, sample_a, sample_b,
    tad_unevaluable_a=None, tad_unevaluable_b=None,
) -> pd.DataFrame:
    """Call conserved, shifted, gained, lost, complex, and ambiguous events."""
    events = _events_at_distance(
        boundaries_a, boundaries_b, jitter, settings, comparison, sample_a, sample_b,
        settings.max_shift_bp, tad_unevaluable_a, tad_unevaluable_b,
    )
    for radius in settings.shift_sensitivity_bp:
        column = f"state_at_{radius}"
        sensitivity_events = _events_at_distance(
            boundaries_a, boundaries_b, jitter, settings, comparison, sample_a, sample_b,
            radius, tad_unevaluable_a, tad_unevaluable_b,
        )
        events[column] = _association_summaries(
            events, sensitivity_events
        )["position_state"].to_numpy()
    state_columns = [f"state_at_{radius}" for radius in settings.shift_sensitivity_bp]
    events["position_state_stable"] = events.apply(
        lambda row: all(row[column] == row.position_state for column in state_columns),
        axis=1,
    )
    return events.loc[:, BOUNDARY_EVENT_COLUMNS]


def compare_event_views(
    reference_events, views_a, views_b, jitter, settings, comparison, sample_a, sample_b,
    tad_unevaluable_a=None, tad_unevaluable_b=None,
) -> pd.DataFrame:
    """Compare reference endpoint states against shared source/LOO boundary views."""
    shared_views = sorted(set(views_a["view"]) & set(views_b["view"]) - {"unified"})
    records = []
    for view in shared_views:
        events = _events_at_distance(
            views_a[views_a["view"] == view], views_b[views_b["view"] == view],
            jitter, settings, comparison, sample_a, sample_b, settings.max_shift_bp,
            tad_unevaluable_a, tad_unevaluable_b,
        )
        summaries = _association_summaries(reference_events, events)
        for reference, summary in zip(reference_events.itertuples(index=False), summaries.itertuples(index=False)):
            if summary.n_associated == 0:
                reason = "no_endpoint_overlap"
            elif summary.n_associated > 1:
                reason = "multiple_endpoint_events"
            elif summary.position_state == "complex":
                reason = "complex_component"
            else:
                reason = "single_endpoint_association"
            records.append({
                "comparison": comparison,
                "boundary_event_id": reference.boundary_event_id,
                "view": view,
                "reference_position_state": reference.position_state,
                "sensitivity_position_state": summary.position_state,
                "concordant": summary.position_state == reference.position_state,
                "reason": reason,
            })
    stability = pd.DataFrame(records, columns=BOUNDARY_STABILITY_COLUMNS)
    stability = stability.sort_values(
        ["comparison", "boundary_event_id", "view"], kind="stable"
    ).reset_index(drop=True)
    if not reference_events.empty:
        agreement = stability.groupby("boundary_event_id", sort=False)["concordant"].all()
        reference_events["position_state_stable"] = reference_events.apply(
            lambda row: bool(row.position_state_stable)
            and bool(agreement.get(row.boundary_event_id, True)),
            axis=1,
        )
    return stability
