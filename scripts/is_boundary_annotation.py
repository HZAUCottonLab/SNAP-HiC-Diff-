"""Associate IS-change loci with endpoint-supported TAD-boundary events."""

import pandas as pd

from is_boundary_io import ANNOTATION_COLUMNS


SIMPLE_POSITION_STATES = {"conserved", "gained_in_B", "lost_in_B", "shifted"}
_SUMMARY_COLUMNS = (
    "confidence_tier", "position_state", "strength_state", "is_consistency_state",
    "evidence_tier", "n_loci", "denominator", "fraction",
)
_STATE_ORDERS = {
    "confidence_tier": ("primary", "single_band_preliminary"),
    "position_state": (
        "conserved", "gained_in_B", "lost_in_B", "shifted", "complex", "ambiguous",
        "no_boundary_association",
    ),
    "strength_state": (
        "strengthened_in_B", "weakened_in_B", "stable_strength",
        "discordant_or_unevaluable", "not_applicable",
    ),
    "is_consistency_state": (
        "direction_concordant", "bidirectional_shift_support", "partial_shift_support",
        "discordant", "unevaluable", "not_applicable",
    ),
    "evidence_tier": ("high", "moderate", "low", "not_applicable"),
}


def locus_anchor(locus, base_bin_bp):
    if locus.is_change_type == "broad_shift":
        return int(locus.event_start), int(locus.event_end)
    return int(locus.summit), int(locus.summit) + base_bin_bp


def _interval_overlap(start_a, end_a, start_b, end_b):
    return max(start_a, start_b) < min(end_a, end_b)


def associated_event_ids(locus, events, base_bin_bp):
    """Return event IDs whose A or B endpoint overlaps the locus anchor."""
    anchor_start, anchor_end = locus_anchor(locus, base_bin_bp)
    identifiers = []
    for event in events.itertuples(index=False):
        supports = []
        if pd.notna(event.a_start) and pd.notna(event.a_end):
            supports.append((int(event.a_start), int(event.a_end)))
        if pd.notna(event.b_start) and pd.notna(event.b_end):
            supports.append((int(event.b_start), int(event.b_end)))
        if any(_interval_overlap(anchor_start, anchor_end, start, end) for start, end in supports):
            identifiers.append(event.boundary_event_id)
    return identifiers


def resolution_count(value):
    if pd.isna(value) or str(value).strip() == "":
        return 0
    return len({int(item) for item in str(value).split(",")})


def direction_support_band_count(value):
    bands = set()
    for role in str(value).split(";"):
        parts = role.split(":", 1)
        for band in parts[-1].split(","):
            if band and band != "nan":
                bands.add(band)
    return len(bands)


def evidence_tier(row):
    if row.position_state in {"complex", "ambiguous"}:
        return "not_applicable"
    if row.confidence_tier == "single_band_preliminary":
        return "low"
    if row.jitter_source != "empirical_q95":
        return "low"
    support_a = row.signal_support_a == "both" or resolution_count(row.resolution_support_a) >= 2
    support_b = row.signal_support_b == "both" or resolution_count(row.resolution_support_b) >= 2
    present_support = {
        "conserved": support_a and support_b,
        "shifted": support_a and support_b,
        "lost_in_B": support_a,
        "gained_in_B": support_b,
        "no_boundary_association": True,
    }[row.position_state]
    if not present_support:
        return "low"
    direction_ok = row.is_consistency_state in {
        "direction_concordant", "bidirectional_shift_support", "not_applicable"
    }
    if row.is_consistency_state != "not_applicable" and direction_support_band_count(row.strength_support_bands) < 2:
        return "low"
    strength_stable = (
        row.position_state == "no_boundary_association"
        or (pd.notna(row.strength_state_stable) and bool(row.strength_state_stable))
    )
    if row.position_state_stable and strength_stable and direction_ok:
        return "high"
    return "moderate"


def _event_order(events):
    return events.sort_values(
        ["chrom", "event_span_start", "event_span_end", "boundary_event_id"],
        kind="stable", na_position="last",
    )


def _mask_reasons(locus, masked, base_bin_bp):
    if masked.empty:
        return ()
    start, end = locus_anchor(locus, base_bin_bp)
    rows = masked[masked["chrom"].eq(locus.chrom)]
    return tuple(sorted({
        str(row.reason) for row in rows.itertuples(index=False)
        if _interval_overlap(start, end, int(row.start), int(row.end))
    }))


def _tad_coverage_reasons(locus, tad_unevaluable, base_bin_bp):
    if tad_unevaluable is None or tad_unevaluable.empty:
        return ()
    start, end = locus_anchor(locus, base_bin_bp)
    rows = tad_unevaluable[tad_unevaluable["chrom"].eq(locus.chrom)]
    return tuple(sorted({
        f"{row.sample}:{row.reason}" for row in rows.itertuples(index=False)
        if _interval_overlap(start, end, int(row.start), int(row.end))
    }))


def _blank_event():
    return {
        "boundary_event_id": "", "position_substate": "", "boundary_a_id": "",
        "boundary_b_id": "", "shift_bp": float("nan"), "signal_support_a": "",
        "signal_support_b": "", "resolution_support_a": "", "resolution_support_b": "",
        "jitter_threshold_bp": float("nan"), "max_shift_bp": float("nan"),
        "jitter_source": "empirical_q95", "position_state_stable": True,
    }


def _blank_strength(position_state):
    return {
        "strength_state": "discordant_or_unevaluable" if position_state == "conserved" else "not_applicable",
        "strength_support_bands": "",
        "is_consistency_state": "not_applicable" if position_state == "conserved" else "unevaluable",
        "strength_state_stable": False,
    }


def _event_values(event):
    return {
        key: getattr(event, key)
        for key in _blank_event()
        if key in event._fields
    }


def annotate_loci(loci, events, strength_summary, masked, settings, tad_unevaluable=None):
    """Annotate each locus once; shifted corridors are never association targets."""
    result_columns = [*loci.columns, *ANNOTATION_COLUMNS]
    if loci.empty:
        return pd.DataFrame(columns=result_columns)

    ordered = _event_order(events)
    event_by_id = {
        event.boundary_event_id: event for event in ordered.itertuples(index=False)
    }
    strengths = {
        row.boundary_event_id: row
        for row in strength_summary.itertuples(index=False)
    }
    rows = []
    for locus in loci.itertuples(index=False):
        local_events = ordered[ordered["chrom"].eq(locus.chrom)]
        ids = associated_event_ids(locus, local_events, settings.base_bin_bp)
        mask_reasons = _mask_reasons(locus, masked, settings.base_bin_bp)
        tad_coverage_reasons = _tad_coverage_reasons(
            locus, tad_unevaluable, settings.base_bin_bp
        )
        event = event_by_id[ids[0]] if len(ids) == 1 else None
        summary = strengths.get(event.boundary_event_id) if event is not None else None
        values = _blank_event()
        if event is not None:
            values.update(_event_values(event))

        if mask_reasons:
            position_state, reason = "ambiguous", f"mask_overlap:{','.join(mask_reasons)}"
        elif not ids and tad_coverage_reasons:
            position_state, reason = "ambiguous", f"tad_coverage:{','.join(tad_coverage_reasons)}"
        elif not ids:
            position_state, reason = "no_boundary_association", "no_overlapping_boundary_event"
        elif len(ids) > 1:
            position_state, reason = "complex", "multiple_boundary_events"
        elif event.position_state == "complex":
            position_state, reason = "complex", "complex_boundary_event"
        elif event.position_state == "ambiguous":
            position_state, reason = "ambiguous", (
                event.ambiguity_reason or "ambiguous_boundary_event"
            )
        else:
            position_state, reason = event.position_state, "unique_boundary_event"
        if event is None or position_state in {"complex", "ambiguous"}:
            strength = _blank_strength(position_state)
        else:
            strength = _blank_strength(position_state) if summary is None else {
                key: getattr(summary, key)
                for key in (
                    "strength_state", "strength_support_bands", "is_consistency_state",
                    "strength_state_stable",
                )
            }
        annotation = {
            **values,
            **strength,
            "boundary_event_id": ",".join(ids),
            "position_state": position_state,
            "evidence_tier": "",
            "n_boundary_candidates": len(ids),
            "annotation_reason": reason,
        }
        tier_row = pd.Series({
            **locus._asdict(), **annotation,
            "jitter_source": values["jitter_source"],
            "position_state_stable": values["position_state_stable"],
        })
        annotation["evidence_tier"] = evidence_tier(tier_row)
        rows.append({**locus._asdict(), **{key: annotation[key] for key in ANNOTATION_COLUMNS}})
    return pd.DataFrame(rows, columns=result_columns)


def summarize_annotations(primary, preliminary):
    combined = pd.concat((primary, preliminary), ignore_index=True)
    if combined.empty:
        return pd.DataFrame(columns=_SUMMARY_COLUMNS)
    groups = [
        "confidence_tier", "position_state", "strength_state", "is_consistency_state",
        "evidence_tier",
    ]
    counts = combined.groupby(groups, dropna=False, sort=False).size().rename("n_loci").reset_index()
    denominators = combined.groupby("confidence_tier", dropna=False).size().rename("denominator")
    counts = counts.join(denominators, on="confidence_tier")
    counts["fraction"] = counts["n_loci"] / counts["denominator"]
    for column, order in _STATE_ORDERS.items():
        counts[column] = pd.Categorical(counts[column], categories=order, ordered=True)
    return counts.sort_values(groups, kind="stable").loc[:, _SUMMARY_COLUMNS].reset_index(drop=True)


def annotation_qc(primary, preliminary, events, strength_long, strength_summary):
    """Return concise evaluability and sensitivity diagnostics for annotation outputs."""
    loci = pd.concat((primary, preliminary), ignore_index=True)
    rows = [
        {"metric": "masked_loci", "value": int(loci["annotation_reason"].str.startswith("mask_overlap:", na=False).sum()), "detail": ""},
        {"metric": "unevaluable_strength_scales", "value": int(strength_long["robust_z"].isna().sum()), "detail": ""},
        {"metric": "fallback_jitter_events", "value": int(events["jitter_source"].eq("fallback_max_tad_resolution").sum()), "detail": ""},
        {"metric": "position_label_changes_across_shift_thresholds", "value": int((~events["position_state_stable"].astype(bool)).sum()), "detail": ""},
        {"metric": "strength_label_changes_across_z_thresholds", "value": int((~strength_summary["strength_state_stable"].astype(bool)).sum()), "detail": ""},
    ]
    mask_counts = {}
    for reason in loci.loc[
        loci["annotation_reason"].str.startswith("mask_overlap:", na=False), "annotation_reason"
    ]:
        for item in str(reason).split(":", 1)[1].split(","):
            mask_counts[item] = mask_counts.get(item, 0) + 1
    for reason, count in sorted(mask_counts.items()):
        rows.append({"metric": "masked_loci_reason", "value": count, "detail": reason})
    for count, total in loci["n_boundary_candidates"].value_counts(sort=False).sort_index().items():
        rows.append({"metric": "candidate_multiplicity", "value": int(total), "detail": str(int(count))})
    return pd.DataFrame(rows, columns=("metric", "value", "detail"))
