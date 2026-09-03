"""Calibrate technical jitter and classify ordered cross-sample boundary events."""

from dataclasses import dataclass
import math

import numpy as np
import pandas as pd

from is_boundary_io import BOUNDARY_EVENT_COLUMNS, BOUNDARY_STABILITY_COLUMNS
from is_boundary_unification import maximum_cardinality_pairs, overlaps


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
    n, m = len(boundaries_a), len(boundaries_b)
    cost = np.full((n + 1, m + 1), np.inf)
    trace = np.full((n + 1, m + 1), "", dtype=object)
    cost[:, 0] = np.arange(n + 1) * 0.5
    cost[0, :] = np.arange(m + 1) * 0.5
    trace[1:, 0] = "delete"
    trace[0, 1:] = "insert"
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            distance = abs(
                int(boundaries_a.iloc[i - 1].center_bp)
                - int(boundaries_b.iloc[j - 1].center_bp)
            )
            candidates = [
                (cost[i - 1, j] + 0.5, 1, "delete"),
                (cost[i, j - 1] + 0.5, 2, "insert"),
            ]
            if distance <= max_distance_bp:
                candidates.append(
                    (cost[i - 1, j - 1] + distance / max_distance_bp, 0, "match")
                )
            best_cost, _, action = min(candidates)
            cost[i, j] = best_cost
            trace[i, j] = action
    pairs = []
    i, j = n, m
    while i or j:
        action = trace[i, j]
        if action == "match":
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif action == "delete":
            i -= 1
        else:
            j -= 1
    return list(reversed(pairs))


def _sorted(frame):
    return frame.sort_values(
        ["start", "end", "unified_boundary_id"], kind="stable"
    ).reset_index(drop=True)


def _candidate_data(boundaries_a, boundaries_b, max_distance_bp):
    by_a = {index: [] for index in range(len(boundaries_a))}
    by_b = {index: [] for index in range(len(boundaries_b))}
    edges = []
    for index_a, item_a in boundaries_a.iterrows():
        for index_b, item_b in boundaries_b.iterrows():
            distance = abs(int(item_a.center_bp) - int(item_b.center_bp))
            if distance <= max_distance_bp:
                edge = (int(index_a), int(index_b), distance)
                edges.append(edge)
                by_a[int(index_a)].append((distance, int(index_b)))
                by_b[int(index_b)].append((distance, int(index_a)))
    for values in (*by_a.values(), *by_b.values()):
        values.sort()

    unseen = {(index_a, index_b) for index_a, index_b, _ in edges}
    components = []
    while unseen:
        index_a, index_b = unseen.pop()
        indices_a, indices_b = {index_a}, {index_b}
        changed = True
        while changed:
            changed = False
            for edge_a, edge_b in list(unseen):
                if edge_a in indices_a or edge_b in indices_b:
                    indices_a.add(edge_a)
                    indices_b.add(edge_b)
                    unseen.remove((edge_a, edge_b))
                    changed = True
        component_edges = [
            edge for edge in edges if edge[0] in indices_a and edge[1] in indices_b
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


def _associated_events(reference, candidates):
    def overlap(left_start, left_end, right_start, right_end):
        return (
            pd.notna(left_start) and pd.notna(left_end)
            and pd.notna(right_start) and pd.notna(right_end)
            and max(int(left_start), int(right_start)) < min(int(left_end), int(right_end))
        )

    return candidates.loc[
        (candidates["chrom"] == reference.chrom)
        & (
            candidates.apply(
                lambda candidate: overlap(
                    reference.a_start, reference.a_end,
                    candidate.a_start, candidate.a_end,
                )
                or overlap(
                    reference.b_start, reference.b_end,
                    candidate.b_start, candidate.b_end,
                ),
                axis=1,
            )
        )
    ]


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
):
    records = []
    chromosomes = tuple(sorted(set(boundaries_a["chrom"]) | set(boundaries_b["chrom"])))
    for chrom in chromosomes:
        all_a = _sorted(boundaries_a[boundaries_a["chrom"] == chrom])
        all_b = _sorted(boundaries_b[boundaries_b["chrom"] == chrom])
        simple_a = _sorted(all_a[all_a["component_state"] == "simple"])
        simple_b = _sorted(all_b[all_b["component_state"] == "simple"])
        for _, row in all_a[all_a["component_state"] != "simple"].iterrows():
            records.append(_event_record(
                [row], [], comparison=comparison, sample_a=sample_a, sample_b=sample_b,
                position_state="complex", match_stage="within_sample_complex", jitter=jitter,
                max_shift_bp=max_shift_bp, ambiguity_reason="complex_within_sample",
            ))
        for _, row in all_b[all_b["component_state"] != "simple"].iterrows():
            records.append(_event_record(
                [], [row], comparison=comparison, sample_a=sample_a, sample_b=sample_b,
                position_state="complex", match_stage="within_sample_complex", jitter=jitter,
                max_shift_bp=max_shift_bp, ambiguity_reason="complex_within_sample",
            ))

        candidates = _candidate_data(simple_a, simple_b, max_shift_bp)
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

        jitter_a = simple_a.iloc[sorted(remaining_a)].copy()
        jitter_b = simple_b.iloc[sorted(remaining_b)].copy()
        jitter_a["_original_index"] = jitter_a.index
        jitter_b["_original_index"] = jitter_b.index
        jitter_a = jitter_a.reset_index(drop=True)
        jitter_b = jitter_b.reset_index(drop=True)
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

        shifted_a = simple_a.iloc[sorted(remaining_a)].copy()
        shifted_b = simple_b.iloc[sorted(remaining_b)].copy()
        shifted_a["_original_index"] = shifted_a.index
        shifted_b["_original_index"] = shifted_b.index
        shifted_a = shifted_a.reset_index(drop=True)
        shifted_b = shifted_b.reset_index(drop=True)
        for local_a, local_b in align_ordered(shifted_a, shifted_b, max_shift_bp):
            index_a = int(shifted_a.iloc[local_a]._original_index)
            index_b = int(shifted_b.iloc[local_b]._original_index)
            item_a, item_b = simple_a.iloc[index_a], simple_b.iloc[index_b]
            distance = abs(int(item_a.center_bp) - int(item_b.center_bp))
            info = _match_info(index_a, index_b, candidates, settings.match_tie_tolerance_bp)
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
            info = {"reciprocal": False, "margin": pd.NA, "n_local_candidates": len(candidates[0][index_a])}
            records.append(_event_record(
                [simple_a.iloc[index_a]], [], comparison=comparison, sample_a=sample_a,
                sample_b=sample_b, position_state="lost_in_B", match_stage="unmatched",
                match_info=info, jitter=jitter, max_shift_bp=max_shift_bp,
            ))
        for index_b in sorted(remaining_b):
            info = {"reciprocal": False, "margin": pd.NA, "n_local_candidates": len(candidates[1][index_b])}
            records.append(_event_record(
                [], [simple_b.iloc[index_b]], comparison=comparison, sample_a=sample_a,
                sample_b=sample_b, position_state="gained_in_B", match_stage="unmatched",
                match_info=info, jitter=jitter, max_shift_bp=max_shift_bp,
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
    boundaries_a, boundaries_b, jitter, settings, comparison, sample_a, sample_b
) -> pd.DataFrame:
    """Call conserved, shifted, gained, lost, complex, and ambiguous events."""
    events = _events_at_distance(
        boundaries_a, boundaries_b, jitter, settings, comparison, sample_a, sample_b,
        settings.max_shift_bp,
    )
    for radius in settings.shift_sensitivity_bp:
        column = f"state_at_{radius}"
        sensitivity_events = _events_at_distance(
            boundaries_a, boundaries_b, jitter, settings, comparison, sample_a, sample_b,
            radius,
        )
        events[column] = [
            _associated_state(row, sensitivity_events) for row in events.itertuples(index=False)
        ]
    state_columns = [f"state_at_{radius}" for radius in settings.shift_sensitivity_bp]
    events["position_state_stable"] = events.apply(
        lambda row: all(row[column] == row.position_state for column in state_columns),
        axis=1,
    )
    return events.loc[:, BOUNDARY_EVENT_COLUMNS]


def compare_event_views(
    reference_events, views_a, views_b, jitter, settings, comparison, sample_a, sample_b
) -> pd.DataFrame:
    """Compare reference endpoint states against shared source/LOO boundary views."""
    view_order = list(dict.fromkeys(views_a["view"].tolist()))
    shared_views = [
        view for view in view_order
        if view != "unified" and view in set(views_b["view"])
    ]
    records = []
    for view in shared_views:
        events = _events_at_distance(
            views_a[views_a["view"] == view], views_b[views_b["view"] == view],
            jitter, settings, comparison, sample_a, sample_b, settings.max_shift_bp,
        )
        for reference in reference_events.itertuples(index=False):
            associated = _associated_events(reference, events)
            state = _associated_state(reference, events)
            if associated.empty:
                reason = "no_endpoint_overlap"
            elif len(associated) > 1:
                reason = "multiple_endpoint_events"
            elif state == "complex":
                reason = "complex_component"
            else:
                reason = "single_endpoint_association"
            records.append({
                "comparison": comparison,
                "boundary_event_id": reference.boundary_event_id,
                "view": view,
                "reference_position_state": reference.position_state,
                "sensitivity_position_state": state,
                "concordant": state == reference.position_state,
                "reason": reason,
            })
    stability = pd.DataFrame(records, columns=BOUNDARY_STABILITY_COLUMNS)
    if not reference_events.empty:
        agreement = stability.groupby("boundary_event_id", sort=False)["concordant"].all()
        reference_events["position_state_stable"] = reference_events.apply(
            lambda row: bool(row.position_state_stable)
            and bool(agreement.get(row.boundary_event_id, True)),
            axis=1,
        )
    return stability
