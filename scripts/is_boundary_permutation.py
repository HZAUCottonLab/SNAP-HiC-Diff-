"""Constrained set-level permutation statistics for IS-boundary association."""

import numpy as np
import pandas as pd

from is_boundary_annotation import associated_event_ids, locus_anchor
from is_boundary_errors import ISBoundaryError


PERMUTATION_COLUMNS = (
    "comparison", "confidence_tier", "statistic", "observed", "null_mean",
    "null_q025", "null_q975", "fold_enrichment", "empirical_p", "bh_qvalue",
    "n_permutations", "seed", "detail",
)
_STATISTICS = (
    "any_boundary_fraction", "rearranged_given_boundary_fraction",
    "gained_in_B_fraction", "lost_in_B_fraction", "shifted_fraction",
    "median_nearest_boundary_distance_bp",
)
_REARRANGED = {"gained_in_B", "lost_in_B", "shifted"}


def _invalid(reason, **context):
    raise ISBoundaryError("IS_BOUNDARY_INPUT_SCHEMA_INVALID", reason=reason, **context)


def _require_columns(frame, columns, label):
    missing = set(columns).difference(frame.columns)
    if missing:
        _invalid(f"{label}_missing_columns", columns=",".join(sorted(missing)))


def _locus_id(locus):
    return str(getattr(locus, "locus_id", getattr(locus, "region_id", "unknown")))


def _mask_frame(masked):
    columns = ("chrom", "start", "end", "reason")
    if masked.empty:
        return pd.DataFrame(columns=columns)
    _require_columns(masked, columns, "masked")
    return masked.loc[:, columns]


def randomize_loci_once(loci, candidate_regions, masked, settings, rng):
    """Relocate each clean locus inside its original candidate region."""
    _require_columns(
        loci,
        ("region_id", "chrom", "event_start", "event_end", "summit", "is_change_type", "mask_status"),
        "loci",
    )
    _require_columns(candidate_regions, ("region_id", "chrom", "region_start", "region_end"), "regions")
    if candidate_regions["region_id"].duplicated().any():
        _invalid("duplicate_candidate_region_id")
    regions = candidate_regions.set_index("region_id")
    masked = _mask_frame(masked)
    randomized = []
    for locus in loci.itertuples(index=False):
        locus_id = _locus_id(locus)
        if locus.mask_status != "clean":
            _invalid("permutation_requires_clean_locus", locus_id=locus_id)
        try:
            region = regions.loc[locus.region_id]
        except KeyError as error:
            _invalid("unknown_candidate_region", locus_id=locus_id)
        if str(region.chrom) != str(locus.chrom):
            _invalid("region_chromosome_mismatch", locus_id=locus_id)
        try:
            start, end = int(locus.event_start), int(locus.event_end)
            region_start, region_end = int(region.region_start), int(region.region_end)
        except (TypeError, ValueError, OverflowError) as error:
            _invalid("permutation_coordinates", locus_id=locus_id)
        length, grid = end - start, settings.base_bin_bp
        if (
            length <= 0 or start % grid or end % grid or region_start < 0
            or region_end <= region_start or region_start % grid or region_end % grid
        ):
            _invalid("permutation_coordinates", locus_id=locus_id)
        first = ((region_start + grid - 1) // grid) * grid
        last = ((region_end - length) // grid) * grid
        candidates = np.arange(first, last + grid, grid, dtype=int)
        local_mask = masked.loc[masked["chrom"].astype(str).eq(str(locus.chrom))]
        eligible = [
            int(candidate) for candidate in candidates
            if not ((local_mask["start"] < candidate + length) & (local_mask["end"] > candidate)).any()
        ]
        if not eligible:
            raise ISBoundaryError(
                "IS_BOUNDARY_ANNOTATION_FAILED", stage="permutation", locus_id=locus_id
            )
        record = locus._asdict()
        record["event_start"] = int(rng.choice(eligible))
        record["event_end"] = record["event_start"] + length
        if record["is_change_type"] != "broad_shift":
            record["summit"] = record["event_start"]
        randomized.append(record)
    return pd.DataFrame(randomized, columns=loci.columns)


def empirical_upper_p(observed, null):
    return float((1 + np.count_nonzero(null >= observed)) / (len(null) + 1))


def empirical_lower_p(observed, null):
    return float((1 + np.count_nonzero(null <= observed)) / (len(null) + 1))


def bh_adjust(pvalues):
    values = np.asarray(pvalues, dtype=float)
    order = np.argsort(values, kind="stable")
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    output = np.empty_like(ranked)
    output[order] = np.clip(ranked, 0.0, 1.0)
    return output


def _nearest_distance(locus, events, base_bin_bp):
    anchor_start, anchor_end = locus_anchor(locus, base_bin_bp)
    distances = []
    for event in events.itertuples(index=False):
        for start, end in ((event.a_start, event.a_end), (event.b_start, event.b_end)):
            if pd.notna(start) and pd.notna(end):
                distances.append(max(int(start) - anchor_end, anchor_start - int(end), 0))
    return min(distances) if distances else float("nan")


def association_statistics(loci, events, base_bin_bp):
    """Return the six pre-specified statistics using endpoint-only association."""
    if loci.empty:
        return {statistic: float("nan") for statistic in _STATISTICS}
    _require_columns(loci, ("chrom", "event_start", "event_end", "summit", "is_change_type"), "loci")
    _require_columns(events, ("chrom", "boundary_event_id", "position_state", "a_start", "a_end", "b_start", "b_end"), "events")
    if events["boundary_event_id"].duplicated().any():
        _invalid("duplicate_boundary_event_id")
    state_by_id = events.set_index("boundary_event_id")["position_state"].to_dict()
    associated, rearranged, gained, lost, shifted, distances = [], [], [], [], [], []
    for locus in loci.itertuples(index=False):
        local_events = events.loc[events["chrom"].astype(str).eq(str(locus.chrom))]
        identifiers = associated_event_ids(locus, local_events, base_bin_bp)
        state = state_by_id[identifiers[0]] if len(identifiers) == 1 else ""
        associated.append(bool(identifiers))
        rearranged.append(state in _REARRANGED)
        gained.append(state == "gained_in_B")
        lost.append(state == "lost_in_B")
        shifted.append(state == "shifted")
        distances.append(_nearest_distance(locus, local_events, base_bin_bp))
    associated_count = sum(associated)
    return {
        "any_boundary_fraction": float(np.mean(associated)),
        "rearranged_given_boundary_fraction": (
            float(sum(rearranged) / associated_count) if associated_count else float("nan")
        ),
        "gained_in_B_fraction": float(np.mean(gained)),
        "lost_in_B_fraction": float(np.mean(lost)),
        "shifted_fraction": float(np.mean(shifted)),
        "median_nearest_boundary_distance_bp": float(np.nanmedian(distances)),
    }


def _stratum_value(loci, column):
    _require_columns(loci, (column,), "loci")
    values = [str(value) for value in loci[column].dropna().unique()]
    if len(values) != 1:
        _invalid(f"permutation_requires_one_{column}")
    return values[0]


def _result_row(comparison, tier, statistic, observed, null, settings):
    valid = np.asarray(null, dtype=float)
    valid = valid[np.isfinite(valid)]
    enough_null = len(valid) >= 0.95 * settings.permutation.n_permutations
    row = {
        "comparison": comparison,
        "confidence_tier": tier,
        "statistic": statistic,
        "observed": observed,
        "null_mean": float("nan"),
        "null_q025": float("nan"),
        "null_q975": float("nan"),
        "fold_enrichment": float("nan"),
        "empirical_p": float("nan"),
        "bh_qvalue": float("nan"),
        "n_permutations": settings.permutation.n_permutations,
        "seed": settings.permutation.seed,
        "detail": "",
    }
    if not enough_null:
        row["detail"] = "insufficient_valid_null_draws"
        return row
    null_mean = float(np.mean(valid))
    row.update(
        null_mean=null_mean,
        null_q025=float(np.quantile(valid, 0.025)),
        null_q975=float(np.quantile(valid, 0.975)),
        fold_enrichment=(float(observed / null_mean) if null_mean else float("nan")),
        empirical_p=(
            empirical_lower_p(observed, valid)
            if statistic == "median_nearest_boundary_distance_bp"
            else empirical_upper_p(observed, valid)
        ),
    )
    return row


def run_permutations(loci, candidate_regions, events, masked, settings):
    """Compare one locus tier with matched candidate-region randomizations."""
    comparison = _stratum_value(loci, "comparison")
    tier = _stratum_value(loci, "confidence_tier")
    _require_columns(events, ("comparison",), "events")
    event_comparisons = {str(value) for value in events["comparison"].dropna()}
    if event_comparisons and event_comparisons != {comparison}:
        _invalid("event_comparison_mismatch")
    observed = association_statistics(loci, events, settings.base_bin_bp)
    null = {statistic: [] for statistic in _STATISTICS}
    rng = np.random.default_rng(settings.permutation.seed)
    for _ in range(settings.permutation.n_permutations):
        randomized = randomize_loci_once(loci, candidate_regions, masked, settings, rng)
        statistics = association_statistics(randomized, events, settings.base_bin_bp)
        for statistic, value in statistics.items():
            null[statistic].append(value)
    rows = [
        _result_row(comparison, tier, statistic, observed[statistic], null[statistic], settings)
        for statistic in _STATISTICS
    ]
    result = pd.DataFrame(rows, columns=PERMUTATION_COLUMNS)
    valid = result["empirical_p"].notna()
    if valid.any():
        result.loc[valid, "bh_qvalue"] = bh_adjust(result.loc[valid, "empirical_p"].to_numpy())
    return result
