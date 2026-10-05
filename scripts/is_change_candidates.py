"""Select and merge final IS-change candidate windows."""

import pandas as pd


SCALE_COLUMNS = [
    "resolution_bp", "window_bins", "window_size_bp", "step_bins", "step_size_bp",
]
IDENTITY_KEYS = [
    "resolution_bp", "window_bins", "chrom", "window_start", "window_end",
]
KNN_COORDS = ["chrom", "window_start", "window_end", "window_bins"]


def _require_columns(frame, required, label):
    missing = set(required).difference(frame.columns)
    if missing:
        raise ValueError(f"{label} missing columns: {', '.join(sorted(missing))}")


def _require_unique_identity(frame, label):
    if frame[IDENTITY_KEYS].isna().any().any():
        raise ValueError(f"{label} has missing full source identity keys")
    if frame.duplicated(IDENTITY_KEYS).any():
        raise ValueError(f"{label} has duplicate full source identity keys")


def select_final_windows(knn, meltron, settings, return_audit=False):
    """Return threshold-passing windows, optionally with matching audit data."""
    _require_columns(
        knn,
        [*SCALE_COLUMNS, "chrom", "window_id", "window_start", "window_end"],
        "KNN",
    )
    _require_columns(
        meltron,
        [
            *SCALE_COLUMNS,
            "chrom",
            "window_id",
            "window_start",
            "window_end",
            "start",
            "end",
            "ID",
            "meltron_best_score",
        ],
        "Meltron",
    )
    _require_unique_identity(knn, "KNN")
    _require_unique_identity(meltron, "Meltron")

    p_columns = sorted(c for c in knn.columns if c.endswith("_empirical_pvalue"))
    if len(p_columns) != 5:
        raise ValueError("KNN input must contain exactly five empirical P-value columns")
    knn = knn.copy()
    knn["knn_any_raw_pvalue_pass"] = (
        knn[p_columns] < settings.raw_knn_pvalue_threshold
    ).any(axis=1)

    missing_keys = meltron[IDENTITY_KEYS].merge(
        knn[IDENTITY_KEYS],
        on=IDENTITY_KEYS,
        how="left",
        indicator=True,
        validate="one_to_one",
    )
    if missing_keys["_merge"].eq("left_only").any():
        raise ValueError("Meltron windows missing KNN full source identity keys")

    merged = knn.merge(
        meltron,
        on=IDENTITY_KEYS,
        how="inner",
        validate="one_to_one",
        suffixes=("", "_meltron"),
    )
    if not merged["window_start"].eq(merged["start"]).all():
        raise ValueError("KNN and Meltron start coordinates disagree for matched windows")

    audit = {
        "matched_windows": len(merged),
        "meltron_end_delta_bp": {
            int(delta): int(count)
            for delta, count in (
                merged["end"].sub(merged["window_end"])
                .value_counts()
                .sort_index()
                .items()
            )
        },
    }
    selected = merged.loc[
        merged["knn_any_raw_pvalue_pass"]
        & (merged["meltron_best_score"] >= settings.meltron_score_threshold)
    ].sort_values(KNN_COORDS, kind="stable").reset_index(drop=True)
    audit["selected_windows"] = len(selected)
    if return_audit:
        return selected, audit
    return selected


def merge_candidate_windows(selected, settings):
    """Merge overlapping or BED-adjacent selected windows within chromosomes."""
    p_columns = sorted(c for c in selected.columns if c.endswith("_empirical_pvalue"))
    columns = [
        "region_id", "chrom", "region_start", "region_end",
        "support_window_bins", "support_scale_count", "source_window_count",
        "member_window_ids", "member_ids", "best_meltron_score",
        *[f"min_{column}" for column in p_columns],
    ]
    regions = []
    region_id = 1

    for chrom, windows in selected.sort_values(KNN_COORDS, kind="stable").groupby(
        "chrom", sort=False
    ):
        current = None
        for _, window in windows.iterrows():
            if current is None or window["window_start"] > current["region_end"]:
                if current is not None:
                    regions.append(current)
                    region_id += 1
                current = {
                    "region_id": region_id,
                    "chrom": chrom,
                    "region_start": window["window_start"],
                    "region_end": window["window_end"],
                    "_window_bins": set(),
                    "member_window_ids": [],
                    "member_ids": [],
                    "source_window_count": 0,
                    "best_meltron_score": window["meltron_best_score"],
                    **{f"min_{column}": window[column] for column in p_columns},
                }
            else:
                current["region_end"] = max(current["region_end"], window["window_end"])
                current["best_meltron_score"] = max(
                    current["best_meltron_score"], window["meltron_best_score"]
                )
                for column in p_columns:
                    current[f"min_{column}"] = min(
                        current[f"min_{column}"], window[column]
                    )
            member_window_id = (
                f"r{int(window['resolution_bp'])}_w{int(window['window_bins'])}:"
                f"{window['chrom']}:{window['window_id']}"
            )
            current["_window_bins"].add(int(window["window_bins"]))
            current["member_window_ids"].append(member_window_id)
            current["member_ids"].append(str(window["ID"]))
            current["source_window_count"] += 1
        if current is not None:
            regions.append(current)
            region_id += 1

    for region in regions:
        region["support_window_bins"] = ",".join(
            str(value) for value in sorted(region.pop("_window_bins"))
        )
        region["support_scale_count"] = len(region["support_window_bins"].split(","))
        region["member_window_ids"] = ",".join(region["member_window_ids"])
        region["member_ids"] = ",".join(region["member_ids"])
    return pd.DataFrame(regions, columns=columns)
