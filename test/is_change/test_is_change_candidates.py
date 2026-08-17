from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys

import pandas as pd
import pytest


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from is_change_candidates import merge_candidate_windows, select_final_windows


SCRIPT = REPO / "scripts" / "build_is_change_candidates.py"
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
        merge_gap_bp=0,
    )


def knn_row(window_id, start, end, chrom="chr1", **values):
    return {
        "chrom": chrom,
        "window_id": window_id,
        "window_start": start,
        "window_end": end,
        **{column: values.get(column, 0.8) for column in P_COLUMNS},
    }


def meltron_row(window_id, start, end, score, chrom="chr1"):
    return {
        "chrom": chrom,
        "window_id": window_id,
        "start": start,
        "end": end,
        "ID": f"{chrom}_{window_id}",
        "meltron_best_score": score,
    }


def test_selects_window_identity_despite_meltron_endpoint_offset(settings):
    knn = pd.DataFrame(
        [
            knn_row(1, 0, 300000, MSE_empirical_pvalue=0.01),
            knn_row(2, 300000, 600000),
            knn_row(3, 600000, 900000, MSE_empirical_pvalue=0.05),
        ]
    )
    knn["MSE_empirical_FDR"] = 1.0
    meltron = pd.DataFrame(
        [
            meltron_row(1, 0, 304000, 5.0),
            meltron_row(2, 300000, 604000, 6.0),
            meltron_row(3, 600000, 904000, 6.0),
        ]
    )

    selected, audit = select_final_windows(knn, meltron, settings, return_audit=True)

    assert selected[["window_id", "window_start", "window_end"]].to_dict("records") == [
        {"window_id": 1, "window_start": 0, "window_end": 300000}
    ]
    assert audit == {
        "matched_windows": 3,
        "selected_windows": 1,
        "meltron_end_delta_bp": {4000: 3},
    }


@pytest.mark.parametrize("n_pvalues", [4, 6])
def test_rejects_non_five_empirical_pvalue_schema(settings, n_pvalues):
    columns = P_COLUMNS[:n_pvalues]
    if n_pvalues == 6:
        columns.append("extra_empirical_pvalue")
    knn = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_id": 1,
                "window_start": 0,
                "window_end": 300000,
                **{column: 0.01 for column in columns},
            }
        ]
    )
    meltron = pd.DataFrame([meltron_row(1, 0, 300000, 6.0)])

    with pytest.raises(ValueError, match="exactly five empirical P-value columns"):
        select_final_windows(knn, meltron, settings)


def test_selects_exact_coordinate_identity_match(settings):
    knn = pd.DataFrame([knn_row(1, 0, 300000, MSE_empirical_pvalue=0.01)])
    meltron = pd.DataFrame([meltron_row(1, 0, 300000, 6.0)])

    selected, audit = select_final_windows(knn, meltron, settings, return_audit=True)

    assert selected["window_id"].tolist() == [1]
    assert audit["meltron_end_delta_bp"] == {0: 1}


def test_rejects_duplicate_knn_identity(settings):
    knn = pd.DataFrame(
        [knn_row(1, 0, 300000, MSE_empirical_pvalue=0.01)] * 2
    )
    meltron = pd.DataFrame([meltron_row(1, 0, 304000, 6.0)])

    with pytest.raises(ValueError, match="KNN has duplicate chrom/window_id keys"):
        select_final_windows(knn, meltron, settings, return_audit=True)


def test_rejects_duplicate_meltron_identity(settings):
    knn = pd.DataFrame([knn_row(1, 0, 300000)])
    meltron = pd.DataFrame(
        [
            meltron_row(1, 0, 304000, 6.0),
            meltron_row(1, 0, 304000, 6.0),
        ]
    )

    with pytest.raises(ValueError, match="Meltron has duplicate chrom/window_id keys"):
        select_final_windows(knn, meltron, settings, return_audit=True)


def test_rejects_meltron_identity_absent_from_knn(settings):
    knn = pd.DataFrame([knn_row(1, 0, 300000)])
    meltron = pd.DataFrame([meltron_row(2, 0, 304000, 6.0)])

    with pytest.raises(ValueError, match="Meltron windows missing KNN"):
        select_final_windows(knn, meltron, settings, return_audit=True)


def test_rejects_shifted_meltron_start(settings):
    knn = pd.DataFrame([knn_row(1, 0, 300000)])
    meltron = pd.DataFrame([meltron_row(1, 4000, 304000, 6.0)])

    with pytest.raises(ValueError, match="start coordinates disagree"):
        select_final_windows(knn, meltron, settings, return_audit=True)


def test_merges_touching_but_not_positive_gap(settings):
    selected = pd.DataFrame(
        [
            {**knn_row(1, 0, 300000, MSE_empirical_pvalue=0.01), **meltron_row(1, 0, 300000, 5.2)},
            {**knn_row(2, 300000, 600000, MSE_empirical_pvalue=0.02), **meltron_row(2, 300000, 600000, 7.0)},
            {**knn_row(3, 604000, 904000, MSE_empirical_pvalue=0.03), **meltron_row(3, 604000, 904000, 8.0)},
        ]
    )

    regions = merge_candidate_windows(selected, settings)

    assert regions[["region_start", "region_end", "n_candidate_windows"]].to_dict("records") == [
        {"region_start": 0, "region_end": 600000, "n_candidate_windows": 2},
        {"region_start": 604000, "region_end": 904000, "n_candidate_windows": 1},
    ]
    assert regions.loc[0, "max_meltron_best_score"] == 7.0
    assert regions.loc[0, "member_window_ids"] == "1,2"
    assert regions.loc[0, "member_ids"] == "chr1_1,chr1_2"


def test_merges_overlaps_without_crossing_chromosome_boundaries(settings):
    selected = pd.DataFrame(
        [
            {**knn_row(1, 0, 400000), **meltron_row(1, 0, 400000, 6.0)},
            {**knn_row(2, 300000, 600000), **meltron_row(2, 300000, 604000, 6.0)},
            {**knn_row(1, 0, 300000, chrom="chr2"), **meltron_row(1, 0, 300000, 6.0, chrom="chr2")},
            {**knn_row(2, 300000, 600000, chrom="chr2"), **meltron_row(2, 300000, 600000, 6.0, chrom="chr2")},
        ]
    )

    regions = merge_candidate_windows(selected, settings)

    assert regions[["region_id", "chrom", "region_start", "region_end"]].to_dict("records") == [
        {"region_id": 1, "chrom": "chr1", "region_start": 0, "region_end": 600000},
        {"region_id": 2, "chrom": "chr2", "region_start": 0, "region_end": 600000},
    ]


def test_preserves_every_metric_regional_minimum(settings):
    first = knn_row(1, 0, 300000)
    second = knn_row(2, 300000, 600000)
    for index, column in enumerate(P_COLUMNS):
        first[column] = 0.5 - index / 100
        second[column] = 0.4 - index / 100
    selected = pd.DataFrame(
        [
            {**first, **meltron_row(1, 0, 300000, 6.0)},
            {**second, **meltron_row(2, 300000, 600000, 7.0)},
        ]
    )

    region = merge_candidate_windows(selected, settings).iloc[0]

    assert {
        f"min_{column}": region[f"min_{column}"] for column in P_COLUMNS
    } == {f"min_{column}": second[column] for column in P_COLUMNS}


def test_cli_writes_sorted_region_schema(tmp_path):
    knn = pd.DataFrame(
        [
            knn_row(2, 300000, 600000, chrom="chr1", MSE_empirical_pvalue=0.02),
            knn_row(1, 0, 300000, chrom="chr1", MSE_empirical_pvalue=0.01),
            knn_row(1, 0, 300000, chrom="chr2", MSE_empirical_pvalue=0.03),
        ]
    )
    meltron = pd.DataFrame(
        [
            meltron_row(2, 300000, 604000, 6.0, chrom="chr1"),
            meltron_row(1, 0, 304000, 6.0, chrom="chr1"),
            meltron_row(1, 0, 304000, 6.0, chrom="chr2"),
        ]
    )
    knn_path = tmp_path / "knn.tsv"
    meltron_path = tmp_path / "meltron.tsv"
    output = tmp_path / "candidate_regions.tsv"
    log = tmp_path / "build.log"
    knn.to_csv(knn_path, sep="\t", index=False)
    meltron.to_csv(meltron_path, sep="\t", index=False)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--knn-enriched", str(knn_path),
            "--meltron-merged", str(meltron_path),
            "--config", str(REPO / "config" / "config.yaml"),
            "--out", str(output),
            "--log", str(log),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    regions = pd.read_csv(output, sep="\t")
    assert regions.columns.tolist() == [
        "region_id", "chrom", "region_start", "region_end",
        "member_window_ids", "member_ids", "n_candidate_windows",
        "max_meltron_best_score", *[f"min_{column}" for column in sorted(P_COLUMNS)],
    ]
    assert regions[["chrom", "region_start", "region_end"]].to_dict("records") == [
        {"chrom": "chr1", "region_start": 0, "region_end": 600000},
        {"chrom": "chr2", "region_start": 0, "region_end": 300000},
    ]
    log_text = log.read_text()
    assert "Matched 3 KNN/Meltron windows by chrom/window_id" in log_text
    assert "Selected 3 windows after KNN and Meltron thresholds" in log_text
    assert "Meltron minus KNN end delta(s) in bp: {4000: 3}" in log_text
    assert "Wrote 2 candidate regions" in log_text


def test_cli_does_not_replace_previous_output_after_invalid_input(tmp_path):
    knn = pd.DataFrame([knn_row(1, 0, 300000, MSE_empirical_pvalue=0.01)]).drop(
        columns=P_COLUMNS[-1]
    )
    meltron = pd.DataFrame([meltron_row(1, 0, 300000, 6.0)])
    knn_path = tmp_path / "knn.tsv"
    meltron_path = tmp_path / "meltron.tsv"
    output = tmp_path / "candidate_regions.tsv"
    knn.to_csv(knn_path, sep="\t", index=False)
    meltron.to_csv(meltron_path, sep="\t", index=False)
    output.write_text("previous valid output\n")

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--knn-enriched", str(knn_path),
            "--meltron-merged", str(meltron_path),
            "--config", str(REPO / "config" / "config.yaml"),
            "--out", str(output),
            "--log", str(tmp_path / "build.log"),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode != 0
    assert output.read_text() == "previous valid output\n"
