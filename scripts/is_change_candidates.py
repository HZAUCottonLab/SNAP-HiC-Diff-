"""Select scale-qualified IS-change windows and union their genomic evidence."""

import math

import pandas as pd


SOURCE_IDENTITY_KEYS = [
    "resolution_bp",
    "window_bins",
    "chrom",
    "window_start",
    "window_end",
]
SORT_COLUMNS = ["chrom", "window_start", "window_end", "window_bins"]


def _require_columns(frame, required, label):
    missing = set(required).difference(frame.columns)
    if missing:
        raise ValueError(f"{label} missing columns: {', '.join(sorted(missing))}")


def _require_unique_identity(frame, label):
    if frame[SOURCE_IDENTITY_KEYS].isna().any().any():
        raise ValueError(f"{label} has missing scale-qualified source identity keys")
    if frame.duplicated(SOURCE_IDENTITY_KEYS).any():
        raise ValueError(f"{label} has duplicate scale-qualified source identities")


def _empirical_pvalue_columns(frame):
    columns = sorted(
        column for column in frame.columns if column.endswith("_empirical_pvalue")
    )
    if len(columns) != 5:
        raise ValueError("KNN input must contain exactly five empirical P-value columns")
    return columns


def _meltron_interval_columns(meltron):
    if {"start", "end"}.issubset(meltron.columns):
        return "start", "end"
    if {"interval_start", "interval_end"}.issubset(meltron.columns):
        return "interval_start", "interval_end"
    raise ValueError("Meltron missing columns: start/end or interval_start/interval_end")


def _member_window_id(window):
    window_id = window["window_id"]
    if isinstance(window_id, float) and math.isfinite(window_id) and window_id.is_integer():
        window_id = int(window_id)
    return (
        f"r{int(window['resolution_bp'])}_w{int(window['window_bins'])}_{window_id}"
    )


def _minimum(current, next_value):
    if pd.isna(current):
        return next_value
    if pd.isna(next_value):
        return current
    return min(current, next_value)


def select_final_windows(knn, meltron, settings, return_audit=False):
    """Return jointly passing windows joined on the complete source identity."""
    _require_columns(
        knn, SOURCE_IDENTITY_KEYS + ["window_id"], "KNN"
    )
    _require_columns(
        meltron,
        SOURCE_IDENTITY_KEYS + ["ID", "meltron_best_score"],
        "Meltron",
    )
    _require_unique_identity(knn, "KNN")
    _require_unique_identity(meltron, "Meltron")
    interval_start, interval_end = _meltron_interval_columns(meltron)
    p_columns = _empirical_pvalue_columns(knn)

    knn = knn.copy()
    knn["knn_any_raw_pvalue_pass"] = (
        knn[p_columns] < settings.raw_knn_pvalue_threshold
    ).any(axis=1)

    source_identities = knn[SOURCE_IDENTITY_KEYS].merge(
        meltron[SOURCE_IDENTITY_KEYS],
        on=SOURCE_IDENTITY_KEYS,
        how="outer",
        indicator=True,
        validate="one_to_one",
    )
    if source_identities["_merge"].eq("right_only").any():
        raise ValueError("Meltron windows missing KNN scale-qualified source identities")
    if source_identities["_merge"].eq("left_only").any():
        raise ValueError("KNN windows missing Meltron scale-qualified source identities")

    merged = knn.merge(
        meltron,
        on=SOURCE_IDENTITY_KEYS,
        how="inner",
        validate="one_to_one",
        suffixes=("", "_meltron"),
    )
    if not merged["window_start"].eq(merged[interval_start]).all():
        raise ValueError("KNN and Meltron start coordinates disagree for matched windows")

    audit = {
        "matched_windows": len(merged),
        "meltron_end_delta_bp": {
            int(delta): int(count)
            for delta, count in (
                merged[interval_end].sub(merged["window_end"])
                .value_counts()
                .sort_index()
                .items()
            )
        },
    }
    selected = merged.loc[
        merged["knn_any_raw_pvalue_pass"]
        & (merged["meltron_best_score"] >= settings.meltron_score_threshold)
    ].sort_values(SORT_COLUMNS, kind="stable").reset_index(drop=True)
    audit["selected_windows"] = len(selected)
    if return_audit:
        return selected, audit
    return selected


def _new_region(window, p_columns):
    return {
        "chrom": window["chrom"],
        "region_start": window["window_start"],
        "region_end": window["window_end"],
        "support_window_bins": {int(window["window_bins"])},
        "member_window_ids": [_member_window_id(window)],
        "member_ids": [str(window["ID"])],
        "source_window_count": 1,
        "best_meltron_score": window["meltron_best_score"],
        **{f"min_{column}": window[column] for column in p_columns},
    }


def _add_window(region, window, p_columns):
    region["region_end"] = max(region["region_end"], window["window_end"])
    region["support_window_bins"].add(int(window["window_bins"]))
    region["member_window_ids"].append(_member_window_id(window))
    region["member_ids"].append(str(window["ID"]))
    region["source_window_count"] += 1
    region["best_meltron_score"] = max(
        region["best_meltron_score"], window["meltron_best_score"]
    )
    for column in p_columns:
        region[f"min_{column}"] = _minimum(
            region[f"min_{column}"], window[column]
        )


def _finalize_region(region):
    region = region.copy()
    bins = sorted(region.pop("support_window_bins"))
    region["support_window_bins"] = ",".join(str(value) for value in bins)
    region["support_scale_count"] = len(bins)
    region["member_window_ids"] = ",".join(region["member_window_ids"])
    region["member_ids"] = ",".join(region["member_ids"])
    return region


def merge_candidate_windows(selected, settings):
    """Merge overlapping or BED-adjacent selected windows by chromosome."""
    del settings
    _require_columns(
        selected,
        SOURCE_IDENTITY_KEYS + ["window_id", "ID", "meltron_best_score"],
        "selected candidates",
    )
    _require_unique_identity(selected, "selected candidates")
    p_columns = _empirical_pvalue_columns(selected)
    columns = [
        "region_id",
        "chrom",
        "region_start",
        "region_end",
        "support_window_bins",
        "support_scale_count",
        "source_window_count",
        "member_window_ids",
        "member_ids",
        "best_meltron_score",
        *[f"min_{column}" for column in p_columns],
    ]
    if selected.empty:
        return pd.DataFrame(columns=columns)

    regions = []
    current = None
    for _, window in selected.sort_values(SORT_COLUMNS, kind="stable").iterrows():
        if (
            current is None
            or window["chrom"] != current["chrom"]
            or window["window_start"] > current["region_end"]
        ):
            if current is not None:
                regions.append(_finalize_region(current))
            current = _new_region(window, p_columns)
        else:
            _add_window(current, window, p_columns)
    regions.append(_finalize_region(current))

    for index, region in enumerate(regions, start=1):
        region["region_id"] = index
    return pd.DataFrame(regions, columns=columns)
