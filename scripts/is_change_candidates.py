"""Select and merge final IS-change candidate windows."""

import pandas as pd


KNN_COORDS = ["chrom", "window_start", "window_end"]


def select_final_windows(knn, meltron, settings):
    """Return coordinate-matched windows passing raw KNN and Meltron thresholds."""
    p_columns = sorted(c for c in knn.columns if c.endswith("_empirical_pvalue"))
    if len(p_columns) != 5:
        raise ValueError("KNN input must contain exactly five empirical P-value columns")
    knn = knn.copy()
    knn["knn_any_raw_pvalue_pass"] = (
        knn[p_columns] < settings.raw_knn_pvalue_threshold
    ).any(axis=1)
    meltron = meltron.rename(columns={"start": "window_start", "end": "window_end"})
    merged = knn.merge(meltron, on=KNN_COORDS, how="inner", validate="one_to_one")
    return merged.loc[
        merged["knn_any_raw_pvalue_pass"]
        & (merged["meltron_best_score"] >= settings.meltron_score_threshold)
    ].sort_values(KNN_COORDS, kind="stable").reset_index(drop=True)


def merge_candidate_windows(selected, settings):
    """Merge overlapping or BED-adjacent selected windows within chromosomes."""
    p_columns = sorted(c for c in selected.columns if c.endswith("_empirical_pvalue"))
    columns = [
        "region_id", "chrom", "region_start", "region_end",
        "member_window_ids", "member_ids", "n_candidate_windows",
        "max_meltron_best_score", *[f"min_{column}" for column in p_columns],
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
                    "member_window_ids": [],
                    "member_ids": [],
                    "n_candidate_windows": 0,
                    "max_meltron_best_score": window["meltron_best_score"],
                    **{f"min_{column}": window[column] for column in p_columns},
                }
            else:
                current["region_end"] = max(current["region_end"], window["window_end"])
                current["max_meltron_best_score"] = max(
                    current["max_meltron_best_score"], window["meltron_best_score"]
                )
                for column in p_columns:
                    current[f"min_{column}"] = min(
                        current[f"min_{column}"], window[column]
                    )
            current["member_window_ids"].append(str(window["window_id"]))
            current["member_ids"].append(str(window["ID"]))
            current["n_candidate_windows"] += 1
        if current is not None:
            regions.append(current)
            region_id += 1

    for region in regions:
        region["member_window_ids"] = ",".join(region["member_window_ids"])
        region["member_ids"] = ",".join(region["member_ids"])
    return pd.DataFrame(regions, columns=columns)
