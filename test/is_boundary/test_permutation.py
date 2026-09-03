from dataclasses import replace
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
import yaml


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from cworld_is_config import parse_cworld_is_settings
from is_boundary_config import parse_is_boundary_settings
from is_boundary_errors import ISBoundaryError
from is_boundary_io import BOUNDARY_EVENT_COLUMNS
from is_boundary_permutation import (
    PERMUTATION_COLUMNS,
    bh_adjust,
    empirical_lower_p,
    empirical_upper_p,
    randomize_loci_once,
    run_permutations,
)
from tad_config import parse_tad_settings


@pytest.fixture
def settings():
    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    return parse_is_boundary_settings(
        config, parse_cworld_is_settings(config), parse_tad_settings(config)
    )


def loci_frame():
    return pd.DataFrame([
        {
            "locus_id": "l1", "comparison": "A_vs_B", "confidence_tier": "primary",
            "region_id": "r1", "chrom": "A02", "event_start": 0, "event_end": 4000,
            "summit": 0, "is_change_type": "transition", "mask_status": "clean",
        },
        {
            "locus_id": "l2", "comparison": "A_vs_B", "confidence_tier": "primary",
            "region_id": "r2", "chrom": "A03", "event_start": 100000, "event_end": 108000,
            "summit": np.nan, "is_change_type": "broad_shift", "mask_status": "clean",
        },
    ])


def regions_frame():
    return pd.DataFrame([
        {"region_id": "r1", "chrom": "A02", "region_start": 0, "region_end": 40000},
        {"region_id": "r2", "chrom": "A03", "region_start": 100000, "region_end": 160000},
    ])


def masked_frame():
    return pd.DataFrame([
        {"chrom": "A02", "start": 12000, "end": 16000, "reason": "IS_NA_A"},
        {"chrom": "A02", "start": 16000, "end": 20000, "reason": "IS_NA_B"},
    ])


def events_frame():
    rows = []
    for event_id, chrom, state, a_start, b_start in (
        ("conserved", "A02", "conserved", 0, 0),
        ("gain", "A03", "gained_in_B", np.nan, 100000),
    ):
        rows.append({
            "comparison": "A_vs_B", "sample_a": "A", "sample_b": "B", "chrom": chrom,
            "boundary_event_id": event_id, "position_state": state,
            "position_substate": state, "boundary_a_id": "a" if pd.notna(a_start) else "",
            "boundary_b_id": "b" if pd.notna(b_start) else "",
            "a_start": a_start, "a_end": a_start + 4000 if pd.notna(a_start) else np.nan,
            "b_start": b_start, "b_end": b_start + 4000 if pd.notna(b_start) else np.nan,
            "event_span_start": min(value for value in (a_start, b_start) if pd.notna(value)),
            "event_span_end": max(value + 4000 for value in (a_start, b_start) if pd.notna(value)),
            "shift_bp": 0, "match_stage": "fixed", "match_is_reciprocal": True,
            "second_best_margin_bp": 4000, "n_local_candidates": 1,
            "jitter_threshold_bp": 4000, "jitter_source": "empirical_q95",
            "max_shift_bp": 100000, "ambiguity_reason": "", "state_at_60000": state,
            "state_at_80000": state, "state_at_100000": state,
            "position_state_stable": True, "signal_support_a": "both",
            "signal_support_b": "both", "resolution_support_a": "4000,8000",
            "resolution_support_b": "4000,8000",
        })
    return pd.DataFrame(rows, columns=BOUNDARY_EVENT_COLUMNS)


def test_randomized_loci_stay_in_candidate_region_and_avoid_mask(settings):
    randomized = randomize_loci_once(
        loci_frame(), regions_frame(), masked_frame(), settings, np.random.default_rng(7)
    )
    joined = randomized.merge(regions_frame(), on="region_id", suffixes=("", "_region"))
    assert (joined.chrom == joined.chrom_region).all()
    assert (joined.event_start >= joined.region_start).all()
    assert (joined.event_end <= joined.region_end).all()
    assert ((joined.event_start % 4000) == 0).all()
    assert (randomized.event_end - randomized.event_start).tolist() == [4000, 8000]
    overlaps = randomized.merge(masked_frame(), on="chrom")
    assert not ((overlaps.event_start < overlaps.end) & (overlaps.event_end > overlaps.start)).any()


def test_randomization_rejects_loci_without_an_eligible_start(settings):
    loci = loci_frame().iloc[[0]].copy()
    regions = regions_frame().iloc[[0]].copy()
    masked = pd.DataFrame([{"chrom": "A02", "start": 0, "end": 40000, "reason": "IS_NA_A"}])

    with pytest.raises(ISBoundaryError, match=r"ERROR\[IS_BOUNDARY_ANNOTATION_FAILED\].*stage=permutation.*l1"):
        randomize_loci_once(loci, regions, masked, settings, np.random.default_rng(1))


def test_empirical_p_has_plus_one_correction_and_bh_is_monotone():
    null = np.array([0.1, 0.2, 0.3, 0.4])
    assert empirical_upper_p(0.35, null) == 0.4
    assert empirical_lower_p(0.15, null) == 0.4
    adjusted = bh_adjust(np.array([0.01, 0.04, 0.03]))
    assert np.all((0 <= adjusted) & (adjusted <= 1))
    assert adjusted[0] <= adjusted[2]


def test_permutations_return_only_set_level_statistics(settings):
    settings = replace(settings, permutation=replace(settings.permutation, n_permutations=4, seed=9))

    output = run_permutations(loci_frame(), regions_frame(), events_frame(), masked_frame(), settings)

    assert output.columns.tolist() == list(PERMUTATION_COLUMNS)
    assert output.statistic.tolist() == [
        "any_boundary_fraction", "rearranged_given_boundary_fraction",
        "gained_in_B_fraction", "lost_in_B_fraction", "shifted_fraction",
        "median_nearest_boundary_distance_bp",
    ]
    assert output.n_permutations.eq(4).all()
    assert output.seed.eq(9).all()
    assert set(output.columns).isdisjoint({"locus_id", "locus_pvalue", "p_value"})
    assert output.loc[output.statistic.eq("gained_in_B_fraction"), "observed"].iloc[0] == 0.5


def test_conditional_rearrangement_statistic_reports_insufficient_null_draws(settings):
    settings = replace(settings, permutation=replace(settings.permutation, n_permutations=4))
    loci = loci_frame().iloc[[0]].assign(event_start=200000, event_end=204000, summit=200000)
    regions = regions_frame().iloc[[0]].assign(region_start=200000, region_end=240000)
    output = run_permutations(loci, regions, events_frame(), masked_frame().iloc[0:0], settings)

    conditional = output.loc[output.statistic.eq("rearranged_given_boundary_fraction")].iloc[0]
    assert conditional.detail == "insufficient_valid_null_draws"
    assert np.isnan(conditional.empirical_p)
    assert np.isnan(conditional.bh_qvalue)


def test_cli_seed_is_byte_deterministic(tmp_path):
    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    config["is_boundary_annotation"]["permutation"]["n_permutations"] = 3
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))
    loci_frame().to_csv(tmp_path / "loci.tsv", sep="\t", index=False)
    regions_frame().to_csv(tmp_path / "regions.tsv", sep="\t", index=False)
    events_frame().to_csv(tmp_path / "events.tsv", sep="\t", index=False)
    masked_frame().to_csv(tmp_path / "masked.tsv", sep="\t", index=False)

    base = [
        sys.executable, str(REPO / "scripts" / "permute_is_boundary_association.py"),
        "--loci", str(tmp_path / "loci.tsv"), "--candidate-regions", str(tmp_path / "regions.tsv"),
        "--events", str(tmp_path / "events.tsv"), "--masked", str(tmp_path / "masked.tsv"),
        "--config", str(config_path), "--log", str(tmp_path / "permutation.log"),
    ]
    outputs = [tmp_path / "one.tsv", tmp_path / "two.tsv"]
    for output in outputs:
        result = subprocess.run(base + ["--output", str(output)], text=True, capture_output=True)
        assert result.returncode == 0, result.stderr
    assert outputs[0].read_bytes() == outputs[1].read_bytes()
