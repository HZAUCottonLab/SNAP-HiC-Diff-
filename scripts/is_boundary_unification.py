import itertools

import pandas as pd

from is_boundary_io import (
    BOUNDARY_VIEW_COLUMNS,
    MULTIRES_INPUT_COLUMNS,
    UNIFIED_BOUNDARY_COLUMNS,
)
JITTER_OFFSET_COLUMNS = (
    "sample", "chrom", "source", "comparison", "resolution_a", "resolution_b",
    "boundary_a_id", "boundary_b_id", "center_a", "center_b", "offset_bp", "eligible",
)


def overlaps(left, right):
    return max(int(left.start), int(right.start)) < min(int(left.end), int(right.end))


def _overlap_pairs(left, right):
    """Enumerate cross-signal interval overlaps without a chromosome-wide grid."""
    events = []
    for side, frame in enumerate((left, right)):
        for index, row in enumerate(frame.itertuples(index=False)):
            events.extend(((int(row.start), 1, side, index), (int(row.end), 0, side, index)))
    active = (set(), set())
    pairs = []
    for _, action, side, index in sorted(events):
        if action == 0:
            active[side].discard(index)
        else:
            other = 1 - side
            pairs.extend(
                (index, other_index) if side == 0 else (other_index, index)
                for other_index in sorted(active[other])
            )
            active[side].add(index)
    return pairs


def maximum_cardinality_pairs(left, right, eligible):
    """Maximize noncrossing overlap matches, then minimize their total shift."""
    m = len(right)
    nodes = [(0, 0, 0, None)]
    tree = [0] * (m + 1)

    def better(left_node, right_node):
        left_score, right_score = nodes[left_node][:2], nodes[right_node][:2]
        return left_node if left_score > right_score or (left_score == right_score and left_node > right_node) else right_node

    def query(index):
        best = 0
        while index:
            best = better(best, tree[index])
            index -= index & -index
        return best

    def update(index, node):
        while index <= m:
            tree[index] = better(tree[index], node)
            index += index & -index

    candidates = []
    for index_left, index_right in _overlap_pairs(left, right):
        row_left, row_right = left.iloc[index_left], right.iloc[index_right]
        if eligible(row_left, row_right):
            candidates.append((index_left, index_right, abs(int(row_left.center_bp) - int(row_right.center_bp))))
    start = 0
    candidates.sort()
    while start < len(candidates):
        index_left = candidates[start][0]
        pending = []
        end = start
        while end < len(candidates) and candidates[end][0] == index_left:
            _, index_right, distance = candidates[end]
            previous = query(index_right)
            count, total_distance = nodes[previous][:2]
            nodes.append((count + 1, total_distance - distance, previous, (index_left, index_right)))
            pending.append((index_right + 1, len(nodes) - 1))
            end += 1
        for index_right, node in pending:
            update(index_right, node)
        start = end

    pairs = []
    node = query(m)
    while node:
        _, _, node, pair = nodes[node]
        pairs.append(pair)
    return list(reversed(pairs))


def _resolution_union(rows):
    return ",".join(map(str, sorted({
        int(value) for item in rows for value in str(item.support_resolutions).split(",")
    })))


def _row(raw_rows, balanced_rows, component_state=None):
    raw_rows = list(raw_rows)
    balanced_rows = list(balanced_rows)
    rows = raw_rows + balanced_rows
    resolutions = _resolution_union(rows)
    if component_state is None:
        coverage = [
            f"{signal}_{reason}"
            for signal, signal_rows in (("raw", raw_rows), ("balanced", balanced_rows))
            for reason, column in (("bad_region", "overlaps_bad_region"), ("nonfinite_base_bins", "has_nonfinite_base_bins"))
            if any(bool(getattr(row, column, False)) for row in signal_rows)
        ]
        terminal = [
            signal for signal, signal_rows in (("raw", raw_rows), ("balanced", balanced_rows))
            if any(bool(getattr(row, "terminal_partial_coordinate", False)) for row in signal_rows)
        ]
        if coverage:
            component_state = "unevaluable_coverage_" + "_and_".join(coverage)
        elif terminal:
            component_state = "simple_terminal_partial_" + "_and_".join(terminal)
        else:
            component_state = "simple"
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
    edges = _overlap_pairs(raw, balanced)
    adjacent_raw = {index: [] for index in range(len(raw))}
    adjacent_balanced = {index: [] for index in range(len(balanced))}
    for raw_index, balanced_index in edges:
        adjacent_raw[raw_index].append(balanced_index)
        adjacent_balanced[balanced_index].append(raw_index)
    components = []
    seen_raw, seen_balanced = set(), set()
    for seed in sorted(index for index, values in adjacent_raw.items() if values):
        if seed in seen_raw:
            continue
        raw_indices, balanced_indices = set(), set()
        pending = [(0, seed)]
        while pending:
            side, index = pending.pop()
            if side == 0:
                if index in seen_raw:
                    continue
                seen_raw.add(index)
                raw_indices.add(index)
                pending.extend((1, other) for other in adjacent_raw[index])
            elif index not in seen_balanced:
                seen_balanced.add(index)
                balanced_indices.add(index)
                pending.extend((0, other) for other in adjacent_balanced[index])
        components.append((sorted(raw_indices), sorted(balanced_indices)))
    return components


def _has_tie(raw, balanced, tie_tolerance_bp):
    return any(
        len(distances) > 1 and distances[1][0] - distances[0][0] < tie_tolerance_bp
        for candidates in _overlap_distances(raw, balanced)
        for distances in candidates.values()
    )


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
    eligibility_columns = [
        column for column in (
            "touches_chromosome_end", "has_nonfinite_base_bins",
            "overlaps_bad_region", "complex_component", "terminal_partial_coordinate",
        ) if column in frame
    ]
    sample = None
    if "sample" in frame:
        sample_values = frame["sample"].dropna().astype(str)
        sample = sample_values.iloc[0] if not sample_values.empty else None

    def start_component(item):
        return {
            "chrom": item.chrom,
            "start": int(item.start),
            "end": int(item.end),
            "resolutions": {int(item.resolution)},
            "n_source_boundaries": int(item.n_source_boundaries),
            "eligibility": {
                column: False if pd.isna(getattr(item, column)) else bool(getattr(item, column))
                for column in eligibility_columns
            },
        }

    def append_component(item, index):
        resolutions = ",".join(map(str, sorted(item["resolutions"])))
        row = {
            "chrom": item["chrom"],
            "start": item["start"],
            "end": item["end"],
            "boundary_id": f"{signal}_boundary_{index:06d}",
            "support_resolutions": resolutions,
            "n_supporting_resolutions": len(item["resolutions"]),
            "n_source_boundaries": item["n_source_boundaries"],
            "signal": signal,
        }
        row.update(item["eligibility"])
        if sample is not None:
            row["sample"] = sample
        rows.append(row)

    ordered = frame.assign(_chrom_rank=frame["chrom"].map(chrom_order)).sort_values(
        ["_chrom_rank", "start", "end"], kind="stable"
    )
    rows = []
    component = None
    for item in ordered.itertuples(index=False):
        if (
            component is None
            or component["chrom"] != item.chrom
            or int(item.start) >= component["end"]
        ):
            if component is not None:
                append_component(component, len(rows) + 1)
            component = start_component(item)
            continue
        component["end"] = max(component["end"], int(item.end))
        component["resolutions"].add(int(item.resolution))
        component["n_source_boundaries"] += int(item.n_source_boundaries)
        for column in eligibility_columns:
            value = getattr(item, column)
            component["eligibility"][column] |= (
                False if pd.isna(value) else bool(value)
            )
    if component is not None:
        append_component(component, len(rows) + 1)
    columns = (*MULTIRES_INPUT_COLUMNS, "signal", *eligibility_columns)
    if sample is not None:
        columns += ("sample",)
    return pd.DataFrame(rows, columns=columns)


def _overlap_distances(left, right):
    by_left = {index: [] for index in range(len(left))}
    by_right = {index: [] for index in range(len(right))}
    for left_index, right_index in _overlap_pairs(left, right):
        distance = abs(int(left.iloc[left_index].center_bp) - int(right.iloc[right_index].center_bp))
        by_left[left_index].append((distance, right_index))
        by_right[right_index].append((distance, left_index))
    for values in (*by_left.values(), *by_right.values()):
        values.sort()
    return by_left, by_right


def _reciprocal_unique(left_index, right_index, distances, tie_tolerance_bp):
    left_distances, right_distances = distances[0][left_index], distances[1][right_index]
    if not left_distances or not right_distances:
        return False
    return (
        (len(left_distances) == 1 or left_distances[1][0] - left_distances[0][0] >= tie_tolerance_bp)
        and (len(right_distances) == 1 or right_distances[1][0] - right_distances[0][0] >= tie_tolerance_bp)
        and left_distances[0][1] == right_index
        and right_distances[0][1] == left_index
    )


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


def _sample(left, right):
    values = []
    for frame in (left, right):
        if "sample" in frame and not frame.empty:
            values.extend(frame["sample"].dropna().astype(str).unique())
    return values[0] if values else ""


def _offset_rows(
    left, right, source, comparison, resolution_a, resolution_b, tie_tolerance_bp
):
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
        distances = _overlap_distances(left_chrom, right_chrom)
        complex_left = set()
        complex_right = set()
        for left_indices, right_indices in _components(left_chrom, right_chrom):
            if len(left_indices) > 1 and len(right_indices) > 1:
                complex_left.update(left_indices)
                complex_right.update(right_indices)
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
                    _reciprocal_unique(left_index, right_index, distances, tie_tolerance_bp)
                    and not (left_index in complex_left and right_index in complex_right)
                    and not _is_excluded(item_a)
                    and not _is_excluded(item_b)
                ),
            })
    return rows


def build_jitter_offsets(
    raw, balanced, per_resolution, base_bin_bp: int, tie_tolerance_bp=None
) -> pd.DataFrame:
    if isinstance(base_bin_bp, bool) or not isinstance(base_bin_bp, int) or base_bin_bp <= 0:
        raise ValueError("base_bin_bp must be a positive integer")
    if tie_tolerance_bp is None:
        tie_tolerance_bp = base_bin_bp
    if (
        isinstance(tie_tolerance_bp, bool)
        or not isinstance(tie_tolerance_bp, int)
        or tie_tolerance_bp <= 0
    ):
        raise ValueError("tie_tolerance_bp must be a positive integer")
    rows = _offset_rows(
        raw, balanced, "raw_balanced", "raw_vs_balanced", "multires", "multires",
        tie_tolerance_bp,
    )
    resolutions = sorted(set(per_resolution.get("resolution", pd.Series(dtype=int))))
    for signal, full in (("raw", raw), ("balanced", balanced)):
        for resolution in resolutions:
            leave_one_out = _merge_resolution_boundaries(
                per_resolution, signal, excluded_resolution=resolution
            )
            rows.extend(_offset_rows(
                full, leave_one_out, signal, f"full_vs_leave_r{resolution}_out",
                "full", f"leave_r{resolution}_out", tie_tolerance_bp,
            ))
    for signal in ("raw", "balanced"):
        signal_rows = per_resolution[per_resolution["signal"] == signal]
        for resolution_a, resolution_b in itertools.combinations(resolutions, 2):
            rows.extend(_offset_rows(
                signal_rows[signal_rows["resolution"] == resolution_a],
                signal_rows[signal_rows["resolution"] == resolution_b],
                signal, f"r{resolution_a}_vs_r{resolution_b}", resolution_a, resolution_b,
                tie_tolerance_bp,
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
