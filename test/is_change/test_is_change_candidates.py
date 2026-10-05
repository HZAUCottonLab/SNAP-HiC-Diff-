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
SCALE_COLUMNS = [
    "resolution_bp", "window_bins", "window_size_bp", "step_bins", "step_size_bp",
]


@pytest.fixture
def settings():
    return SimpleNamespace(
        raw_knn_pvalue_threshold=0.05,
        meltron_score_threshold=5.0,
        merge_gap_bp=0,
    )


def scale_values(window_bins=75):
    step_bins = window_bins // 5
    return {
        "resolution_bp": 4000,
        "window_bins": window_bins,
        "window_size_bp": 4000 * window_bins,
        "step_bins": step_bins,
        "step_size_bp": 4000 * step_bins,
    }


def knn_row(window_id, start, end, chrom="chr1", window_bins=75, **values):
    return {
        **scale_values(window_bins),
        "chrom": chrom,
        "window_id": window_id,
        "window_start": start,
        "window_end": end,
        **{column: values.get(column, 0.8) for column in P_COLUMNS},
    }


def meltron_row(
    window_id,
    start,
    end,
    score,
    chrom="chr1",
    window_bins=75,
    source_start=None,
    source_end=None,
):
    return {
        **scale_values(window_bins),
        "chrom": chrom,
        "window_id": window_id,
        "window_start": start if source_start is None else source_start,
        "window_end": end if source_end is None else source_end,
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
            meltron_row(1, 0, 304000, 5.0, source_end=300000),
            meltron_row(2, 300000, 604000, 6.0, source_end=600000),
            meltron_row(3, 600000, 904000, 6.0, source_end=900000),
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
                **scale_values(),
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
    meltron = pd.DataFrame([meltron_row(1, 0, 304000, 6.0, source_end=300000)])

    with pytest.raises(ValueError, match="KNN has duplicate full source identity keys"):
        select_final_windows(knn, meltron, settings, return_audit=True)


def test_rejects_duplicate_meltron_identity(settings):
    knn = pd.DataFrame([knn_row(1, 0, 300000)])
    meltron = pd.DataFrame(
        [
            meltron_row(1, 0, 304000, 6.0, source_end=300000),
            meltron_row(1, 0, 304000, 6.0, source_end=300000),
        ]
    )

    with pytest.raises(ValueError, match="Meltron has duplicate full source identity keys"):
        select_final_windows(knn, meltron, settings, return_audit=True)


@pytest.mark.parametrize(
    ("source", "column", "label"),
    [
        ("knn", "resolution_bp", "KNN"),
        ("knn", "window_bins", "KNN"),
        ("knn", "chrom", "KNN"),
        ("knn", "window_start", "KNN"),
        ("knn", "window_end", "KNN"),
        ("meltron", "resolution_bp", "Meltron"),
        ("meltron", "window_bins", "Meltron"),
        ("meltron", "chrom", "Meltron"),
        ("meltron", "window_start", "Meltron"),
        ("meltron", "window_end", "Meltron"),
    ],
)
def test_rejects_missing_identity_values(settings, source, column, label):
    knn = pd.DataFrame([knn_row(1, 0, 300000)])
    meltron = pd.DataFrame([meltron_row(1, 0, 304000, 6.0, source_end=300000)])
    (knn if source == "knn" else meltron).loc[0, column] = None

    with pytest.raises(ValueError, match=f"{label} has missing full source identity keys"):
        select_final_windows(knn, meltron, settings, return_audit=True)


def test_rejects_meltron_identity_absent_from_knn(settings):
    knn = pd.DataFrame([knn_row(1, 0, 300000)])
    meltron = pd.DataFrame([meltron_row(2, 0, 304000, 6.0)])

    with pytest.raises(ValueError, match="Meltron windows missing KNN"):
        select_final_windows(knn, meltron, settings, return_audit=True)


def test_rejects_shifted_meltron_start(settings):
    knn = pd.DataFrame([knn_row(1, 0, 300000)])
    meltron = pd.DataFrame([
        meltron_row(
            1, 4000, 304000, 6.0, source_start=0, source_end=300000
        )
    ])

    with pytest.raises(ValueError, match="start coordinates disagree"):
        select_final_windows(knn, meltron, settings, return_audit=True)


def test_same_window_id_at_two_scales_is_not_a_collision(settings):
    knn = pd.DataFrame([
        knn_row(1, 0, 100000, window_bins=25, MSE_empirical_pvalue=0.01),
        knn_row(1, 0, 200000, window_bins=50, MSE_empirical_pvalue=0.02),
    ])
    meltron = pd.DataFrame([
        meltron_row(1, 0, 104000, 6.0, window_bins=25, source_end=100000),
        meltron_row(1, 0, 204000, 7.0, window_bins=50, source_end=200000),
    ])

    selected = select_final_windows(knn, meltron, settings)

    assert selected["window_bins"].tolist() == [25, 50]


def test_cross_scale_union_records_support_and_keeps_single_scale_regions(settings):
    selected = pd.DataFrame([
        {**knn_row(1, 0, 100000, window_bins=25), **meltron_row(
            1, 0, 104000, 6.0, window_bins=25, source_end=100000
        )},
        {**knn_row(1, 0, 200000, window_bins=50), **meltron_row(
            1, 0, 204000, 8.0, window_bins=50, source_end=200000
        )},
        {**knn_row(2, 400000, 500000, window_bins=25), **meltron_row(
            2, 400000, 504000, 5.0, window_bins=25, source_end=500000
        )},
    ])

    regions = merge_candidate_windows(selected, settings)

    assert regions[[
        "region_start", "region_end", "support_window_bins",
        "support_scale_count", "source_window_count", "best_meltron_score",
    ]].to_dict("records") == [
        {
            "region_start": 0,
            "region_end": 200000,
            "support_window_bins": "25,50",
            "support_scale_count": 2,
            "source_window_count": 2,
            "best_meltron_score": 8.0,
        },
        {
            "region_start": 400000,
            "region_end": 500000,
            "support_window_bins": "25",
            "support_scale_count": 1,
            "source_window_count": 1,
            "best_meltron_score": 5.0,
        },
    ]


def test_merges_touching_but_not_positive_gap(settings):
    selected = pd.DataFrame(
        [
            {**knn_row(1, 0, 300000, MSE_empirical_pvalue=0.01), **meltron_row(1, 0, 300000, 5.2)},
            {**knn_row(2, 300000, 600000, MSE_empirical_pvalue=0.02), **meltron_row(2, 300000, 600000, 7.0)},
            {**knn_row(3, 604000, 904000, MSE_empirical_pvalue=0.03), **meltron_row(3, 604000, 904000, 8.0)},
        ]
    )

    regions = merge_candidate_windows(selected, settings)

    assert regions[["region_start", "region_end", "source_window_count"]].to_dict("records") == [
        {"region_start": 0, "region_end": 600000, "source_window_count": 2},
        {"region_start": 604000, "region_end": 904000, "source_window_count": 1},
    ]
    assert regions.loc[0, "best_meltron_score"] == 7.0
    assert regions.loc[0, "support_window_bins"] == "75"
    assert regions.loc[0, "support_scale_count"] == 1
    assert regions.loc[0, "member_window_ids"] == "r4000_w75:chr1:1,r4000_w75:chr1:2"
    assert regions.loc[0, "member_ids"] == "chr1_1,chr1_2"


def test_merges_overlaps_without_crossing_chromosome_boundaries(settings):
    selected = pd.DataFrame(
        [
            {**knn_row(1, 0, 400000), **meltron_row(1, 0, 400000, 6.0)},
            {**knn_row(2, 300000, 600000), **meltron_row(2, 300000, 604000, 6.0, source_end=600000)},
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


def test_cli_writes_sorted_evidence_and_union_outputs(tmp_path):
    knn = pd.DataFrame(
        [
            knn_row(2, 300000, 600000, chrom="chr1", MSE_empirical_pvalue=0.02),
            knn_row(1, 0, 300000, chrom="chr1", MSE_empirical_pvalue=0.01),
            knn_row(1, 0, 300000, chrom="chr2", MSE_empirical_pvalue=0.03),
        ]
    )
    meltron = pd.DataFrame(
        [
            meltron_row(2, 300000, 604000, 6.0, chrom="chr1", source_end=600000),
            meltron_row(1, 0, 304000, 6.0, chrom="chr1", source_end=300000),
            meltron_row(1, 0, 304000, 6.0, chrom="chr2", source_end=300000),
        ]
    )
    knn_path = tmp_path / "knn.tsv"
    meltron_path = tmp_path / "meltron.tsv"
    candidates_output = tmp_path / "candidates.tsv"
    union_output = tmp_path / "candidate_union.tsv"
    bed_output = tmp_path / "candidate_union.bed"
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
            "--out-candidates", str(candidates_output),
            "--out-union-tsv", str(union_output),
            "--out-union-bed", str(bed_output),
            "--log", str(log),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    candidates = pd.read_csv(candidates_output, sep="\t")
    assert set(SCALE_COLUMNS).issubset(candidates.columns)
    assert candidates[["chrom", "window_start", "window_end"]].to_dict("records") == [
        {"chrom": "chr1", "window_start": 0, "window_end": 300000},
        {"chrom": "chr1", "window_start": 300000, "window_end": 600000},
        {"chrom": "chr2", "window_start": 0, "window_end": 300000},
    ]

    regions = pd.read_csv(union_output, sep="\t")
    assert regions.columns.tolist() == [
        "region_id", "chrom", "region_start", "region_end",
        "support_window_bins", "support_scale_count", "source_window_count",
        "member_window_ids", "member_ids", "best_meltron_score",
        *[f"min_{column}" for column in sorted(P_COLUMNS)],
    ]
    assert regions[["chrom", "region_start", "region_end"]].to_dict("records") == [
        {"chrom": "chr1", "region_start": 0, "region_end": 600000},
        {"chrom": "chr2", "region_start": 0, "region_end": 300000},
    ]
    assert [line.split("\t") for line in bed_output.read_text().splitlines()] == [
        ["chr1", "0", "600000", "1"],
        ["chr2", "0", "300000", "2"],
    ]
    log_text = log.read_text()
    assert f"{knn_path}: matched=3 selected=3" in log_text


def test_cli_does_not_replace_previous_outputs_after_invalid_input(tmp_path):
    knn = pd.DataFrame([knn_row(1, 0, 300000, MSE_empirical_pvalue=0.01)]).drop(
        columns=P_COLUMNS[-1]
    )
    meltron = pd.DataFrame([meltron_row(1, 0, 300000, 6.0)])
    knn_path = tmp_path / "knn.tsv"
    meltron_path = tmp_path / "meltron.tsv"
    candidates_output = tmp_path / "candidates.tsv"
    union_output = tmp_path / "candidate_union.tsv"
    bed_output = tmp_path / "candidate_union.bed"
    knn.to_csv(knn_path, sep="\t", index=False)
    meltron.to_csv(meltron_path, sep="\t", index=False)
    candidates_output.write_text("previous candidates\n")
    union_output.write_text("previous union\n")
    bed_output.write_text("previous bed\n")

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--knn-enriched", str(knn_path),
            "--meltron-merged", str(meltron_path),
            "--config", str(REPO / "config" / "config.yaml"),
            "--out-candidates", str(candidates_output),
            "--out-union-tsv", str(union_output),
            "--out-union-bed", str(bed_output),
            "--log", str(tmp_path / "build.log"),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode != 0
    assert candidates_output.read_text() == "previous candidates\n"
    assert union_output.read_text() == "previous union\n"
    assert bed_output.read_text() == "previous bed\n"


def test_cli_writes_header_only_outputs_when_no_windows_are_selected(tmp_path):
    knn = pd.DataFrame(columns=knn_row(1, 0, 300000).keys())
    meltron = pd.DataFrame(columns=meltron_row(1, 0, 304000, 6.0).keys())
    knn_path = tmp_path / "knn.tsv"
    meltron_path = tmp_path / "meltron.tsv"
    candidates_output = tmp_path / "candidates.tsv"
    union_output = tmp_path / "candidate_union.tsv"
    bed_output = tmp_path / "candidate_union.bed"
    knn.to_csv(knn_path, sep="\t", index=False)
    meltron.to_csv(meltron_path, sep="\t", index=False)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--knn-enriched", str(knn_path),
            "--meltron-merged", str(meltron_path),
            "--config", str(REPO / "config" / "config.yaml"),
            "--out-candidates", str(candidates_output),
            "--out-union-tsv", str(union_output),
            "--out-union-bed", str(bed_output),
            "--log", str(tmp_path / "build.log"),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert pd.read_csv(candidates_output, sep="\t").empty
    assert pd.read_csv(union_output, sep="\t").empty
    assert bed_output.read_text() == ""
