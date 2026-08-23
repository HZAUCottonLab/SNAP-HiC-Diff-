from pathlib import Path
from types import SimpleNamespace
import sys

import pandas as pd
import pytest


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

import build_is_change_candidates as candidate_builder
from is_change_candidates import select_final_windows


P_COLUMNS = [
    "MSE_empirical_pvalue",
    "d_SSIM_empirical_pvalue",
    "d_Spearman_empirical_pvalue",
    "tw_score_empirical_pvalue",
    "d_NCC_empirical_pvalue",
]


@pytest.fixture
def settings():
    return SimpleNamespace(
        raw_knn_pvalue_threshold=0.05,
        meltron_score_threshold=5.0,
    )


def knn_row(window_id, start, end, window_bins=75, chrom="chr1", **values):
    return {
        "resolution_bp": 4000,
        "window_bins": window_bins,
        "window_size_bp": window_bins * 4000,
        "step_bins": window_bins // 5,
        "step_size_bp": window_bins * 800,
        "chrom": chrom,
        "window_id": window_id,
        "window_start": start,
        "window_end": end,
        **{column: values.get(column, 0.8) for column in P_COLUMNS},
    }


def meltron_row(window_id, start, end, score, window_bins=75, chrom="chr1"):
    return {
        "resolution_bp": 4000,
        "window_bins": window_bins,
        "window_size_bp": window_bins * 4000,
        "step_bins": window_bins // 5,
        "step_size_bp": window_bins * 800,
        "chrom": chrom,
        "window_id": window_id,
        "window_start": start,
        "window_end": end - 4000,
        "interval_start": start,
        "interval_end": end,
        "ID": f"r4000_w{window_bins}_{chrom}_{window_id}",
        "meltron_best_score": score,
    }


def test_selects_by_complete_scale_qualified_identity_and_audits_endpoint_offset(settings):
    knn = pd.DataFrame(
        [
            knn_row(1, 0, 300000, MSE_empirical_pvalue=0.01),
            knn_row(2, 300000, 600000),
        ]
    )
    meltron = pd.DataFrame(
        [
            meltron_row(1, 0, 304000, 5.0),
            meltron_row(2, 300000, 604000, 6.0),
        ]
    )

    selected, audit = select_final_windows(knn, meltron, settings, return_audit=True)

    assert selected[["window_bins", "window_id"]].to_dict("records") == [
        {"window_bins": 75, "window_id": 1}
    ]
    assert audit == {
        "matched_windows": 2,
        "selected_windows": 1,
        "meltron_end_delta_bp": {4000: 2},
    }


@pytest.mark.parametrize("n_pvalues", [4, 6])
def test_rejects_non_five_empirical_pvalue_schema(settings, n_pvalues):
    knn = pd.DataFrame([knn_row(1, 0, 300000)])
    if n_pvalues == 4:
        knn = knn.drop(columns=[P_COLUMNS[-1]])
    else:
        knn["extra_empirical_pvalue"] = 0.01
    meltron = pd.DataFrame([meltron_row(1, 0, 304000, 6.0)])

    with pytest.raises(ValueError, match="exactly five empirical P-value columns"):
        select_final_windows(knn, meltron, settings)


def test_rejects_duplicate_source_identity_even_when_window_ids_differ(settings):
    knn = pd.DataFrame(
        [
            knn_row(1, 0, 300000, MSE_empirical_pvalue=0.01),
            knn_row(2, 0, 300000, MSE_empirical_pvalue=0.01),
        ]
    )
    meltron = pd.DataFrame([meltron_row(1, 0, 304000, 6.0)])

    with pytest.raises(ValueError, match="KNN has duplicate scale-qualified source identities"):
        select_final_windows(knn, meltron, settings)


def test_rejects_meltron_source_identity_absent_from_knn(settings):
    knn = pd.DataFrame([knn_row(1, 0, 300000)])
    meltron = pd.DataFrame([meltron_row(1, 300000, 604000, 6.0)])

    with pytest.raises(ValueError, match="Meltron windows missing KNN"):
        select_final_windows(knn, meltron, settings)


def test_rejects_knn_source_identity_absent_from_meltron(settings):
    knn = pd.DataFrame(
        [
            knn_row(1, 0, 300000),
            knn_row(2, 300000, 600000),
        ]
    )
    meltron = pd.DataFrame([meltron_row(1, 0, 304000, 6.0)])

    with pytest.raises(ValueError, match="KNN windows missing Meltron"):
        select_final_windows(knn, meltron, settings)


def test_rejects_shifted_meltron_algorithm_interval_start(settings):
    knn = pd.DataFrame([knn_row(1, 0, 300000)])
    meltron = pd.DataFrame([meltron_row(1, 4000, 304000, 6.0)])
    meltron["window_start"] = 0
    meltron["window_end"] = 300000

    with pytest.raises(ValueError, match="start coordinates disagree"):
        select_final_windows(knn, meltron, settings)


def test_output_publication_restores_every_existing_file_after_third_promotion_fails(
    monkeypatch, tmp_path
):
    outputs = [
        tmp_path / "all.tsv",
        tmp_path / "regions.tsv",
        tmp_path / "regions.bed",
    ]
    original_bytes = [
        b"old all candidates\n",
        b"old union regions\n",
        b"old bed intervals\n",
    ]
    for output, content in zip(outputs, original_bytes):
        output.write_bytes(content)

    real_replace = candidate_builder.os.replace
    promotions = []

    def fail_third_promotion(source, destination):
        source = Path(source)
        destination = Path(destination)
        if destination in outputs and source.suffix == ".tmp":
            promotions.append((source, destination))
            if len(promotions) == 3:
                raise OSError("injected third promotion failure")
        return real_replace(source, destination)

    monkeypatch.setattr(candidate_builder.os, "replace", fail_third_promotion)
    all_candidates = pd.DataFrame({"candidate": ["new"]})
    regions = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "region_start": [0],
            "region_end": [4000],
            "region_id": [1],
        }
    )

    with pytest.raises(OSError, match="injected third promotion failure"):
        candidate_builder.write_outputs_atomically(
            all_candidates, regions, *outputs
        )

    assert len(promotions) == 3
    assert [output.read_bytes() for output in outputs] == original_bytes
    assert not list(tmp_path.glob(".*.tmp"))
    assert not list(tmp_path.glob(".*.bak"))
