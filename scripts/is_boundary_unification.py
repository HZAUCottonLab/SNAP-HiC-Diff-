import itertools

import pandas as pd

from is_boundary_io import (
    BOUNDARY_VIEW_COLUMNS,
    MULTIRES_INPUT_COLUMNS,
    UNIFIED_BOUNDARY_COLUMNS,
)
from tad_boundaries import Interval, merge_intervals


JITTER_OFFSET_COLUMNS = (
    "sample", "chrom", "source", "comparison", "resolution_a", "resolution_b",
    "boundary_a_id", "boundary_b_id", "center_a", "center_b", "offset_bp", "eligible",
)


def overlaps(left, right):
    return max(int(left.start), int(right.start)) < min(int(left.end), int(right.end))


def maximum_cardinality_pairs(left, right, eligible):
    n, m = len(left), len(right)
    score = [[(0, 0) for _ in range(m + 1)] for _ in range(n + 1)]
    trace = [[None for _ in range(m + 1)] for _ in range(n + 1)]
    for i in range(n + 1):
        for j in range(m + 1):
            if i < n and score[i][j] > score[i + 1][j]:
                score[i + 1][j] = score[i][j]
                trace[i + 1][j] = (i, j, "skip_left")
            if j < m and score[i][j] > score[i][j + 1]:
                score[i][j + 1] = score[i][j]
                trace[i][j + 1] = (i, j, "skip_right")
            if i < n and j < m and eligible(left.iloc[i], right.iloc[j]):
                distance = abs(int(left.iloc[i].center_bp) - int(right.iloc[j].center_bp))
                candidate = (score[i][j][0] + 1, score[i][j][1] - distance)
                if candidate > score[i + 1][j + 1]:
                    score[i + 1][j + 1] = candidate
                    trace[i + 1][j + 1] = (i, j, "match")
    pairs = []
    i, j = n, m
    while i or j:
        previous = trace[i][j]
        if previous is None:
            if i:
                i -= 1
            else:
                j -= 1
            continue
        previous_i, previous_j, action = previous
        if action == "match":
            pairs.append((previous_i, previous_j))
        i, j = previous_i, previous_j
    return list(reversed(pairs))


def _resolution_union(rows):
    return ",".join(map(str, sorted({
        int(value) for item in rows for value in str(item.support_resolutions).split(",")
    })))


def _row(raw_rows, balanced_rows, component_state="simple"):
    rows = list(raw_rows) + list(balanced_rows)
    resolutions = _resolution_union(rows)
    return {
        "chrom": rows[0].chrom,
        "start": min(int(item.start) for item in rows),
        "end": max(int(item.end) for item in rows),
        "unified_boundary_id": "",
        "center_bp": 0,
        "component_state": component_state,
        "signal_support": "both" if raw_rows and balanced_rows else ("raw_only" if raw_rows else "balanced_only"),
        "raw_boundary_ids": ",".join(str(item.boundary_id) for item in raw_rows),
        "balanced_boundary_ids": ",".join(str(item.boundary_id) for item in balanced_rows),
        "support_resolutions": resolutions,
        "n_supporting_resolutions": len(resolutions.split(",")),
        "n_source_boundaries": sum(int(item.n_source_boundaries) for item in rows),
    }


def _components(raw, balanced):
    edges = {
        (raw_index, balanced_index)
        for raw_index, raw_row in raw.iterrows()
        for balanced_index, balanced_row in balanced.iterrows()
        if overlaps(raw_row, balanced_row)
    }
    components = []
    unseen = set(edges)
    while unseen:
        seed = unseen.pop()
        raw_indices, balanced_indices = {seed[0]}, {seed[1]}
        changed = True
        while changed:
            changed = False
            for edge in list(unseen):
                if edge[0] in raw_indices or edge[1] in balanced_indices:
                    raw_indices.add(edge[0])
                    balanced_indices.add(edge[1])
                    unseen.remove(edge)
                    changed = True
        components.append((sorted(raw_indices), sorted(balanced_indices)))
    return components


def _has_tie(raw, balanced, tie_tolerance_bp):
    for left, right in ((raw, balanced), (balanced, raw)):
        for _, item in left.iterrows():
            distances = sorted(
                abs(int(item.center_bp) - int(other.center_bp))
                for _, other in right.iterrows()
                if overlaps(item, other)
            )
            if len(distances) > 1 and distances[1] - distances[0] < tie_tolerance_bp:
                return True
    return False


def _chrom_rank(chromosomes, chrom_order):
    if chrom_order is None:
        return {chrom: index for index, chrom in enumerate(chromosomes)}
    if isinstance(chrom_order, dict):
        ranks = dict(chrom_order)
    else:
        ranks = {chrom: index for index, chrom in enumerate(chrom_order)}
    missing = set(chromosomes).difference(ranks)
    if missing:
        raise ValueError(f"chrom_order missing chromosomes: {','.join(sorted(missing))}")
    return ranks


def unify_boundaries(raw, balanced, tie_tolerance_bp, chrom_order=None) -> pd.DataFrame:
    raw = raw.copy().assign(center_bp=lambda frame: (frame["start"] + frame["end"]) // 2)
    balanced = balanced.copy().assign(center_bp=lambda frame: (frame["start"] + frame["end"]) // 2)
    chromosomes = list(dict.fromkeys([*raw.get("chrom", []), *balanced.get("chrom", [])]))
    ranks = _chrom_rank(chromosomes, chrom_order)
    output = []
    for chrom in sorted(chromosomes, key=ranks.__getitem__):
        left = raw[raw["chrom"] == chrom].sort_values(["start", "end", "boundary_id"], kind="stable").reset_index(drop=True)
        right = balanced[balanced["chrom"] == chrom].sort_values(["start", "end", "boundary_id"], kind="stable").reset_index(drop=True)
        used_left, used_right = set(), set()
        for left_indices, right_indices in _components(left, right):
            component_left = left.iloc[left_indices].reset_index(drop=True)
            component_right = right.iloc[right_indices].reset_index(drop=True)
            used_left.update(left_indices)
            used_right.update(right_indices)
            if (
                len(component_left) > 1 and len(component_right) > 1
            ) or _has_tie(component_left, component_right, tie_tolerance_bp):
                output.append(_row(component_left.itertuples(index=False), component_right.itertuples(index=False), "complex_within_sample"))
                continue
            pairs = maximum_cardinality_pairs(component_left, component_right, overlaps)
            matched_left = {index for index, _ in pairs}
            matched_right = {index for _, index in pairs}
            for left_index, right_index in pairs:
                output.append(_row([component_left.iloc[left_index]], [component_right.iloc[right_index]]))
            for index, item in component_left.iterrows():
                if index not in matched_left:
                    output.append(_row([item], []))
            for index, item in component_right.iterrows():
                if index not in matched_right:
                    output.append(_row([], [item]))
        for index, item in left.iterrows():
            if index not in used_left:
                output.append(_row([item], []))
        for index, item in right.iterrows():
            if index not in used_right:
                output.append(_row([], [item]))
    frame = pd.DataFrame(output, columns=UNIFIED_BOUNDARY_COLUMNS)
    if frame.empty:
        return frame
    frame["center_bp"] = (frame["start"] + frame["end"]) // 2
    frame = frame.assign(_chrom_rank=frame["chrom"].map(ranks)).sort_values(
        ["_chrom_rank", "start", "end", "raw_boundary_ids", "balanced_boundary_ids"], kind="stable"
    ).drop(columns="_chrom_rank").reset_index(drop=True)
    frame["unified_boundary_id"] = [
        f"unified_boundary_{index:06d}" for index in range(1, len(frame) + 1)
    ]
    return frame.loc[:, UNIFIED_BOUNDARY_COLUMNS]


def _empty_multires():
    return pd.DataFrame(columns=(*MULTIRES_INPUT_COLUMNS, "signal"))


def _merge_resolution_boundaries(per_resolution, signal, excluded_resolution=None):
    frame = per_resolution[per_resolution["signal"] == signal]
    if excluded_resolution is not None:
        frame = frame[frame["resolution"] != excluded_resolution]
    if frame.empty:
        return _empty_multires()
    chromosomes = list(dict.fromkeys(frame["chrom"]))
    chrom_order = {chrom: index for index, chrom in enumerate(chromosomes)}
    intervals = [
        Interval(
            row.chrom, int(row.start), int(row.end), {int(row.resolution)},
            int(row.n_source_boundaries),
        )
        for row in frame.itertuples(index=False)
    ]
    eligibility_columns = [
        column for column in (
            "touches_chromosome_end", "has_nonfinite_base_bins",
            "overlaps_bad_region", "complex_component",
        ) if column in frame
    ]
    rows = []
    for index, item in enumerate(merge_intervals(intervals, chrom_order), 1):
        resolutions = ",".join(map(str, sorted(item.resolutions)))
        source_rows = frame[
            (frame["chrom"] == item.chrom)
            & (frame["start"] < item.end)
            & (frame["end"] > item.start)
        ]
        row = {
            "chrom": item.chrom,
            "start": item.start,
            "end": item.end,
            "boundary_id": f"{signal}_boundary_{index:06d}",
            "support_resolutions": resolutions,
            "n_supporting_resolutions": len(item.resolutions),
            "n_source_boundaries": item.n_source_boundaries,
            "signal": signal,
        }
        row.update({column: bool(source_rows[column].any()) for column in eligibility_columns})
        if "sample" in frame:
            row["sample"] = frame["sample"].dropna().astype(str).iloc[0]
        rows.append(row)
    columns = (*MULTIRES_INPUT_COLUMNS, "signal", *eligibility_columns)
    if "sample" in frame:
        columns += ("sample",)
    return pd.DataFrame(rows, columns=columns)


def _candidate_distances(left, right, item):
    return sorted(
        (abs(int(item.center_bp) - int(other.center_bp)), index)
        for index, other in right.iterrows()
        if overlaps(item, other)
    )


def _reciprocal_unique(left, right, left_index, right_index):
    left_item = left.iloc[left_index]
    right_item = right.iloc[right_index]
    left_distances = _candidate_distances(left, right, left_item)
    right_distances = _candidate_distances(right, left, right_item)
    if not left_distances or not right_distances:
        return False
    return (
        len(left_distances) == 1 or left_distances[0][0] < left_distances[1][0]
    ) and (
        len(right_distances) == 1 or right_distances[0][0] < right_distances[1][0]
    ) and left_distances[0][1] == right_index and right_distances[0][1] == left_index


def _is_excluded(item):
    return (
        item.get("component_state") == "complex_within_sample"
        or any(
            bool(item.get(column, False))
            for column in (
                "touches_chromosome_end", "has_nonfinite_base_bins",
                "overlaps_bad_region", "complex_component",
            )
        )
    )


def _in_complex_component(left, right, left_index, right_index):
    left_indices, right_indices = {left_index}, {right_index}
    changed = True
    while changed:
        changed = False
        for index, item in left.iterrows():
            if index not in left_indices and any(
                overlaps(item, right.iloc[other]) for other in right_indices
            ):
                left_indices.add(index)
                changed = True
        for index, item in right.iterrows():
            if index not in right_indices and any(
                overlaps(item, left.iloc[other]) for other in left_indices
            ):
                right_indices.add(index)
                changed = True
    return len(left_indices) > 1 and len(right_indices) > 1


def _sample(left, right):
    values = []
    for frame in (left, right):
        if "sample" in frame and not frame.empty:
            values.extend(frame["sample"].dropna().astype(str).unique())
    return values[0] if values else ""


def _offset_rows(left, right, source, comparison, resolution_a, resolution_b):
    left = left.copy().assign(center_bp=lambda frame: (frame["start"] + frame["end"]) // 2)
    right = right.copy().assign(center_bp=lambda frame: (frame["start"] + frame["end"]) // 2)
    rows = []
    for chrom in sorted(set(left["chrom"]) & set(right["chrom"])):
        left_chrom = left[left["chrom"] == chrom].sort_values(
            ["start", "end", "boundary_id"], kind="stable"
        ).reset_index(drop=True)
        right_chrom = right[right["chrom"] == chrom].sort_values(
            ["start", "end", "boundary_id"], kind="stable"
        ).reset_index(drop=True)
        for left_index, right_index in maximum_cardinality_pairs(
            left_chrom, right_chrom, overlaps
        ):
            item_a = left_chrom.iloc[left_index]
            item_b = right_chrom.iloc[right_index]
            rows.append({
                "sample": _sample(left_chrom, right_chrom),
                "chrom": chrom,
                "source": source,
                "comparison": comparison,
                "resolution_a": resolution_a,
                "resolution_b": resolution_b,
                "boundary_a_id": item_a.boundary_id,
                "boundary_b_id": item_b.boundary_id,
                "center_a": int(item_a.center_bp),
                "center_b": int(item_b.center_bp),
                "offset_bp": abs(int(item_a.center_bp) - int(item_b.center_bp)),
                "eligible": bool(
                    _reciprocal_unique(left_chrom, right_chrom, left_index, right_index)
                    and not _in_complex_component(
                        left_chrom, right_chrom, left_index, right_index
                    )
                    and not _is_excluded(item_a)
                    and not _is_excluded(item_b)
                ),
            })
    return rows


def build_jitter_offsets(raw, balanced, per_resolution, base_bin_bp: int) -> pd.DataFrame:
    if isinstance(base_bin_bp, bool) or not isinstance(base_bin_bp, int) or base_bin_bp <= 0:
        raise ValueError("base_bin_bp must be a positive integer")
    rows = _offset_rows(
        raw, balanced, "raw_balanced", "raw_vs_balanced", "multires", "multires"
    )
    resolutions = sorted(set(per_resolution.get("resolution", pd.Series(dtype=int))))
    for signal, full in (("raw", raw), ("balanced", balanced)):
        for resolution in resolutions:
            leave_one_out = _merge_resolution_boundaries(
                per_resolution, signal, excluded_resolution=resolution
            )
            rows.extend(_offset_rows(
                full, leave_one_out, signal, f"full_vs_leave_r{resolution}_out",
                "full", f"leave_r{resolution}_out",
            ))
    for signal in ("raw", "balanced"):
        signal_rows = per_resolution[per_resolution["signal"] == signal]
        for resolution_a, resolution_b in itertools.combinations(resolutions, 2):
            rows.extend(_offset_rows(
                signal_rows[signal_rows["resolution"] == resolution_a],
                signal_rows[signal_rows["resolution"] == resolution_b],
                signal, f"r{resolution_a}_vs_r{resolution_b}", resolution_a, resolution_b,
            ))
    return pd.DataFrame(rows, columns=JITTER_OFFSET_COLUMNS)


def _view_rows(sample, view, frame):
    if frame.empty:
        return pd.DataFrame(columns=BOUNDARY_VIEW_COLUMNS)
    return frame.assign(sample=sample, view=view).loc[:, BOUNDARY_VIEW_COLUMNS]


def build_sensitivity_views(sample, raw, balanced, per_resolution, settings) -> pd.DataFrame:
    tie_tolerance_bp = settings.match_tie_tolerance_bp
    empty = _empty_multires()
    views = [
        _view_rows(sample, "unified", unify_boundaries(raw, balanced, tie_tolerance_bp)),
        _view_rows(sample, "raw_only_input", unify_boundaries(raw, empty, tie_tolerance_bp)),
        _view_rows(sample, "balanced_only_input", unify_boundaries(empty, balanced, tie_tolerance_bp)),
    ]
    for resolution in (4000, 8000, 20000):
        views.append(_view_rows(
            sample,
            f"leave_r{resolution}_out",
            unify_boundaries(
                _merge_resolution_boundaries(per_resolution, "raw", resolution),
                _merge_resolution_boundaries(per_resolution, "balanced", resolution),
                tie_tolerance_bp,
            ),
        ))
    return pd.concat(views, ignore_index=True).loc[:, BOUNDARY_VIEW_COLUMNS]
