import os
from pathlib import Path
import tempfile

import pandas as pd

from is_boundary_errors import ISBoundaryError


MULTIRES_INPUT_COLUMNS = (
    "chrom", "start", "end", "boundary_id", "support_resolutions",
    "n_supporting_resolutions", "n_source_boundaries",
)
UNIFIED_BOUNDARY_COLUMNS = (
    "chrom", "start", "end", "unified_boundary_id", "center_bp",
    "component_state", "signal_support", "raw_boundary_ids",
    "balanced_boundary_ids", "support_resolutions",
    "n_supporting_resolutions", "n_source_boundaries",
)
BOUNDARY_VIEW_COLUMNS = (
    "sample", "view", "chrom", "start", "end", "unified_boundary_id",
    "center_bp", "component_state", "signal_support", "raw_boundary_ids",
    "balanced_boundary_ids", "support_resolutions",
    "n_supporting_resolutions", "n_source_boundaries",
)
TAD_UNEVALUABLE_COLUMNS = (
    "sample", "chrom", "start", "end", "reason",
)
BOUNDARY_EVENT_COLUMNS = (
    "comparison", "sample_a", "sample_b", "chrom", "boundary_event_id",
    "position_state", "position_substate", "boundary_a_id", "boundary_b_id",
    "a_start", "a_end", "b_start", "b_end", "event_span_start",
    "event_span_end", "shift_bp", "match_stage", "match_is_reciprocal",
    "second_best_margin_bp", "n_local_candidates", "jitter_threshold_bp",
    "jitter_source", "max_shift_bp", "ambiguity_reason", "state_at_60000",
    "state_at_80000", "state_at_100000", "position_state_stable",
    "signal_support_a", "signal_support_b", "resolution_support_a",
    "resolution_support_b",
)
BOUNDARY_STABILITY_COLUMNS = (
    "comparison", "boundary_event_id", "view", "reference_position_state",
    "sensitivity_position_state", "concordant", "reason",
)
ANNOTATION_COLUMNS = (
    "boundary_event_id", "position_state", "position_substate", "strength_state",
    "strength_support_bands", "is_consistency_state", "evidence_tier", "boundary_a_id",
    "boundary_b_id", "shift_bp", "signal_support_a", "signal_support_b",
    "resolution_support_a", "resolution_support_b", "n_boundary_candidates",
    "jitter_threshold_bp", "max_shift_bp", "annotation_reason",
)
STRENGTH_LONG_COLUMNS = (
    "comparison", "boundary_event_id", "chrom", "position_state",
    "coordinate_role", "IS_distance", "scale_band", "a_start", "a_end",
    "b_start", "b_end", "strength_a", "strength_b", "delta_strength",
    "robust_z", "valley_a_start", "valley_b_start", "detection_a",
    "detection_b", "valid_fraction_a", "valid_fraction_b", "reference_n",
    "reference_median", "reference_mad",
)
STRENGTH_SUMMARY_COLUMNS = (
    "comparison", "boundary_event_id", "position_state", "strength_state",
    "strength_support_bands", "is_consistency_state", "small_band_score",
    "medium_band_score", "large_band_score", "n_evaluable_bands",
    "state_at_z1_5", "state_at_z2_0", "state_at_z2_5",
    "strength_state_stable",
)
RESOLUTION_INPUT_COLUMNS = (
    "chrom", "start", "end", "boundary_id", "resolution", "n_source_boundaries",
)


def _resolution_string(value):
    values = tuple(sorted({int(item) for item in str(value).split(",") if item}))
    if not values or any(item <= 0 for item in values):
        raise ValueError("resolution list must contain positive integers")
    return ",".join(map(str, values))


def _invalid(path, **context):
    raise ISBoundaryError("IS_BOUNDARY_INPUT_SCHEMA_INVALID", path=path, **context)


def _read_bed(path, columns):
    path = Path(path)
    try:
        frame = pd.read_csv(path, sep="\t", header=None, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        # A TAD caller may legitimately yield no boundaries for one signal.
        # Keep that state distinct from a malformed or unreadable BED.
        return pd.DataFrame(columns=columns)
    except Exception as error:
        raise ISBoundaryError(
            "IS_BOUNDARY_INPUT_SCHEMA_INVALID", path=path, reason=type(error).__name__
        ) from error
    if frame.empty:
        return pd.DataFrame(columns=columns)
    if frame.shape[1] != len(columns):
        _invalid(path, reason="column_count")
    frame.columns = columns
    return frame


def _integer_column(frame, column, path):
    try:
        values = pd.to_numeric(frame[column], errors="raise")
    except Exception as error:
        raise ISBoundaryError(
            "IS_BOUNDARY_INPUT_SCHEMA_INVALID", path=path, reason=f"invalid_{column}"
        ) from error
    if not values.astype(float).map(float.is_integer).all():
        _invalid(path, reason=f"invalid_{column}")
    return values.astype(int)


def _validate_coordinates(frame, path):
    frame["start"] = _integer_column(frame, "start", path)
    frame["end"] = _integer_column(frame, "end", path)
    if (frame["start"] < 0).any() or (frame["end"] <= frame["start"]).any():
        _invalid(path, reason="coordinates")


def _validate_order(frame, path):
    chrom_order = {chrom: index for index, chrom in enumerate(dict.fromkeys(frame["chrom"]))}
    ordered = (
        frame.assign(_chrom_rank=frame["chrom"].map(chrom_order))
        .sort_values(["_chrom_rank", "start", "end", "boundary_id"], kind="stable")
        .drop(columns="_chrom_rank")
    )
    if ordered.index.tolist() != frame.index.tolist():
        _invalid(path, reason="unsorted_rows")


def read_multires_boundaries(path, signal: str) -> pd.DataFrame:
    if signal not in {"raw", "balanced"}:
        raise ValueError("signal must be raw or balanced")
    frame = _read_bed(path, MULTIRES_INPUT_COLUMNS)
    if frame.empty:
        return frame.assign(signal=pd.Series(dtype=str))
    _validate_coordinates(frame, path)
    if not frame["chrom"].str.len().all() or not frame["boundary_id"].str.len().all():
        _invalid(path, reason="blank_identifier")
    if frame["boundary_id"].duplicated().any():
        _invalid(path, reason="duplicate_boundary_id")
    try:
        frame["support_resolutions"] = frame["support_resolutions"].map(_resolution_string)
    except (TypeError, ValueError) as error:
        raise ISBoundaryError(
            "IS_BOUNDARY_INPUT_SCHEMA_INVALID", path=path, reason="support_resolutions"
        ) from error
    frame["n_supporting_resolutions"] = _integer_column(
        frame, "n_supporting_resolutions", path
    )
    frame["n_source_boundaries"] = _integer_column(frame, "n_source_boundaries", path)
    expected_support = frame["support_resolutions"].str.count(",") + 1
    if (
        (frame["n_supporting_resolutions"] != expected_support).any()
        or (frame["n_source_boundaries"] <= 0).any()
    ):
        _invalid(path, reason="support_counts")
    _validate_order(frame, path)
    frame["signal"] = signal
    return frame.reset_index(drop=True)


def read_resolution_boundaries(paths, signal: str) -> pd.DataFrame:
    if signal not in {"raw", "balanced"}:
        raise ValueError("signal must be raw or balanced")
    frames = []
    for resolution, path in paths:
        if isinstance(resolution, bool) or not isinstance(resolution, int) or resolution <= 0:
            raise ValueError("resolution must be a positive integer")
        frame = _read_bed(path, RESOLUTION_INPUT_COLUMNS)
        if frame.empty:
            continue
        _validate_coordinates(frame, path)
        if not frame["chrom"].str.len().all() or not frame["boundary_id"].str.len().all():
            _invalid(path, reason="blank_identifier")
        frame["resolution"] = _integer_column(frame, "resolution", path)
        frame["n_source_boundaries"] = _integer_column(frame, "n_source_boundaries", path)
        if (frame["resolution"] != resolution).any() or (frame["n_source_boundaries"] <= 0).any():
            _invalid(path, reason="resolution_or_source_count")
        if frame["boundary_id"].duplicated().any():
            _invalid(path, reason="duplicate_boundary_id")
        _validate_order(frame, path)
        frames.append(frame.assign(signal=signal))
    if not frames:
        return pd.DataFrame(columns=(*RESOLUTION_INPUT_COLUMNS, "signal"))
    return pd.concat(frames, ignore_index=True)


def _read_tsv(path, columns):
    path = Path(path)
    try:
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except Exception as error:
        raise ISBoundaryError(
            "IS_BOUNDARY_INPUT_SCHEMA_INVALID", path=path, reason=type(error).__name__
        ) from error
    if tuple(frame.columns) != tuple(columns):
        _invalid(path, reason="header")
    return frame


def _typed_boundary_frame(path, columns, unique_columns=("unified_boundary_id",)):
    frame = _read_tsv(path, columns)
    for column in ("start", "end", "center_bp"):
        frame[column] = _integer_column(frame, column, path)
    _validate_coordinates(frame, path)
    if not frame["chrom"].str.len().all() or not frame["unified_boundary_id"].str.len().all():
        _invalid(path, reason="blank_identifier")
    if unique_columns and frame.duplicated(list(unique_columns)).any():
        _invalid(path, reason="duplicate_boundary_id")
    if not frame["center_bp"].eq((frame["start"] + frame["end"]) // 2).all():
        _invalid(path, reason="center_bp")
    return frame


def read_unified_boundaries(path) -> pd.DataFrame:
    return _typed_boundary_frame(path, UNIFIED_BOUNDARY_COLUMNS)


def read_boundary_views(path) -> pd.DataFrame:
    frame = _typed_boundary_frame(
        path, BOUNDARY_VIEW_COLUMNS, ("sample", "view", "unified_boundary_id")
    )
    if not frame["sample"].str.len().all() or not frame["view"].str.len().all():
        _invalid(path, reason="blank_view")
    return frame


def read_tad_unevaluable(path) -> pd.DataFrame:
    """Read the per-sample intervals where a TAD absence call is invalid."""
    frame = _read_tsv(path, TAD_UNEVALUABLE_COLUMNS)
    if frame.empty:
        return frame
    _validate_coordinates(frame, path)
    if not (
        frame["sample"].str.len().all()
        and frame["chrom"].str.len().all()
        and frame["reason"].str.len().all()
    ):
        _invalid(path, reason="blank_tad_unevaluable_field")
    return frame.reset_index(drop=True)


def read_jitter_offsets(path) -> pd.DataFrame:
    columns = (
        "sample", "chrom", "source", "comparison", "resolution_a", "resolution_b",
        "boundary_a_id", "boundary_b_id", "center_a", "center_b", "offset_bp", "eligible",
    )
    frame = _read_tsv(path, columns)
    for column in ("center_a", "center_b", "offset_bp"):
        frame[column] = _integer_column(frame, column, path)
    if (frame["offset_bp"] < 0).any() or not frame["sample"].str.len().all():
        _invalid(path, reason="offsets")
    if not frame["eligible"].isin(["True", "False", "true", "false", "1", "0"]).all():
        _invalid(path, reason="eligible")
    frame["eligible"] = frame["eligible"].isin(["True", "true", "1"])
    return frame


def atomic_write_tsv(frame: pd.DataFrame, path, columns) -> None:
    if tuple(frame.columns) != tuple(columns):
        raise ISBoundaryError(
            "IS_BOUNDARY_INPUT_SCHEMA_INVALID", path=path, reason="output_schema"
        )
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=target.parent, prefix=f".{target.name}.", suffix=".tmp"
    )
    os.close(descriptor)
    try:
        frame.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, target)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def atomic_write_tables(tables) -> None:
    """Stage schema-checked TSVs before replacing their published paths."""
    temporary = []
    try:
        for frame, path, columns in tables:
            if tuple(frame.columns) != tuple(columns):
                raise ISBoundaryError(
                    "IS_BOUNDARY_INPUT_SCHEMA_INVALID", path=path, reason="output_schema"
                )
            target = Path(path)
            target.parent.mkdir(parents=True, exist_ok=True)
            descriptor, name = tempfile.mkstemp(
                dir=target.parent, prefix=f".{target.name}.", suffix=".tmp"
            )
            os.close(descriptor)
            temporary.append((Path(name), target))
            frame.to_csv(name, sep="\t", index=False)
        for name, target in temporary:
            os.replace(name, target)
    except BaseException:
        for name, _ in temporary:
            name.unlink(missing_ok=True)
        raise
