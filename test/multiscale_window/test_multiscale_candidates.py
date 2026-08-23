from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys

import pandas as pd


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from is_change_candidates import merge_candidate_windows, select_final_windows


P_COLUMNS = [
    "MSE_empirical_pvalue",
    "d_SSIM_empirical_pvalue",
    "d_Spearman_empirical_pvalue",
    "tw_score_empirical_pvalue",
    "d_NCC_empirical_pvalue",
]


def settings():
    return SimpleNamespace(
        raw_knn_pvalue_threshold=0.05,
        meltron_score_threshold=5.0,
        resolution_bp=4000,
    )


def knn_row(window_bins, window_id, start, end, chrom="chr1", **pvalues):
    return {
        "resolution_bp": 4000,
        "window_bins": window_bins,
        "window_size_bp": window_bins * 4000,
        "step_bins": max(1, window_bins // 5),
        "step_size_bp": max(1, window_bins // 5) * 4000,
        "chrom": chrom,
        "window_start": start,
        "window_end": end,
        "window_id": window_id,
        **{column: pvalues.get(column, 0.8) for column in P_COLUMNS},
    }


def meltron_row(window_bins, window_id, start, end, score, chrom="chr1"):
    return {
        "resolution_bp": 4000,
        "window_bins": window_bins,
        "window_size_bp": window_bins * 4000,
        "step_bins": max(1, window_bins // 5),
        "step_size_bp": max(1, window_bins // 5) * 4000,
        "chrom": chrom,
        "window_start": start,
        "window_end": end,
        "interval_start": start,
        "interval_end": end,
        "window_id": window_id,
        "ID": f"r4000_w{window_bins}_{chrom}_{window_id}",
        "meltron_best_score": score,
    }


def test_cross_scale_candidate_contract_is_exercised_through_public_module():
    """The implementation must expose a scale-qualified candidate API."""
    from is_change_candidates import SOURCE_IDENTITY_KEYS

    assert SOURCE_IDENTITY_KEYS == [
        "resolution_bp", "window_bins", "chrom", "window_start", "window_end",
    ]


def test_cross_scale_selection_and_bed_adjacency_union_preserve_scale_evidence():
    """Three jointly passing scales form one region; one-scale support is valid."""
    per_scale = [
        (knn_row(25, 1, 0, 100, MSE_empirical_pvalue=0.01), meltron_row(25, 1, 0, 100, 5.5)),
        (knn_row(50, 1, 50, 250, d_SSIM_empirical_pvalue=0.02), meltron_row(50, 1, 50, 250, 7.0)),
        (knn_row(75, 1, 250, 550, d_Spearman_empirical_pvalue=0.03), meltron_row(75, 1, 250, 550, 6.5)),
        (knn_row(100, 1, 600, 1000, tw_score_empirical_pvalue=0.01), meltron_row(100, 1, 600, 1000, 4.0)),
        (knn_row(125, 1, 1200, 1700), meltron_row(125, 1, 1200, 1700, 8.0)),
    ]

    selected = pd.concat(
        [
            select_final_windows(pd.DataFrame([knn]), pd.DataFrame([meltron]), settings())
            for knn, meltron in per_scale
        ],
        ignore_index=True,
    )
    regions = merge_candidate_windows(selected, settings())

    assert selected["window_bins"].tolist() == [25, 50, 75]
    assert regions.to_dict("records") == [
        {
            "region_id": 1,
            "chrom": "chr1",
            "region_start": 0,
            "region_end": 550,
            "support_window_bins": "25,50,75",
            "support_scale_count": 3,
            "source_window_count": 3,
            "member_window_ids": "r4000_w25_1,r4000_w50_1,r4000_w75_1",
            "member_ids": "r4000_w25_chr1_1,r4000_w50_chr1_1,r4000_w75_chr1_1",
            "best_meltron_score": 7.0,
            "min_MSE_empirical_pvalue": 0.01,
            "min_d_NCC_empirical_pvalue": 0.8,
            "min_d_SSIM_empirical_pvalue": 0.02,
            "min_d_Spearman_empirical_pvalue": 0.03,
            "min_tw_score_empirical_pvalue": 0.8,
        }
    ]
SCRIPT = REPO / "scripts" / "build_is_change_candidates.py"


def empty_knn_frame(window_bins):
    return pd.DataFrame(columns=knn_row(window_bins, 1, 0, window_bins * 4).keys())


def empty_meltron_frame(window_bins):
    return pd.DataFrame(
        columns=meltron_row(window_bins, 1, 0, window_bins * 4, 5.0).keys()
    )


def run_cli(knn_paths, meltron_paths, expected_bins, out_all, out_regions, out_bed, log):
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--knn-enriched",
            *map(str, knn_paths),
            "--meltron-merged",
            *map(str, meltron_paths),
            "--expected-window-bins",
            *map(str, expected_bins),
            "--out-all",
            str(out_all),
            "--out-regions",
            str(out_regions),
            "--out-bed",
            str(out_bed),
            "--config",
            str(REPO / "config" / "config.yaml"),
            "--log",
            str(log),
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def test_duplicate_scale_qualified_identity_fails_selection():
    duplicate = knn_row(25, 2, 0, 100, MSE_empirical_pvalue=0.01)
    knn = pd.DataFrame(
        [knn_row(25, 1, 0, 100, MSE_empirical_pvalue=0.01), duplicate]
    )
    meltron = pd.DataFrame([meltron_row(25, 1, 0, 100, 6.0)])

    import pytest

    with pytest.raises(ValueError, match="duplicate scale-qualified source identities"):
        select_final_windows(knn, meltron, settings())


def test_union_never_merges_adjacent_windows_across_chromosomes():
    selected = pd.concat(
        [
            select_final_windows(
                pd.DataFrame([knn_row(25, 1, 0, 100, chrom="chr1", MSE_empirical_pvalue=0.01)]),
                pd.DataFrame([meltron_row(25, 1, 0, 100, 6.0, chrom="chr1")]),
                settings(),
            ),
            select_final_windows(
                pd.DataFrame([knn_row(50, 1, 100, 300, chrom="chr2", MSE_empirical_pvalue=0.01)]),
                pd.DataFrame([meltron_row(50, 1, 100, 300, 6.0, chrom="chr2")]),
                settings(),
            ),
        ],
        ignore_index=True,
    )

    regions = merge_candidate_windows(selected, settings())

    assert regions[["chrom", "region_start", "region_end", "support_scale_count"]].to_dict("records") == [
        {"chrom": "chr1", "region_start": 0, "region_end": 100, "support_scale_count": 1},
        {"chrom": "chr2", "region_start": 100, "region_end": 300, "support_scale_count": 1},
    ]


def test_all_empty_scale_inputs_have_stable_ledger_and_region_schemas():
    selected = pd.concat(
        [
            select_final_windows(
                empty_knn_frame(window_bins), empty_meltron_frame(window_bins), settings()
            )
            for window_bins in (25, 50, 75, 100, 125)
        ],
        ignore_index=True,
    )
    regions = merge_candidate_windows(selected, settings())

    assert selected.empty
    assert regions.empty
    assert regions.columns.tolist() == [
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
        "min_MSE_empirical_pvalue",
        "min_d_NCC_empirical_pvalue",
        "min_d_SSIM_empirical_pvalue",
        "min_d_Spearman_empirical_pvalue",
        "min_tw_score_empirical_pvalue",
    ]


def test_three_output_cli_preserves_all_existing_outputs_when_validation_fails(tmp_path):
    bins = [25, 50, 75, 100, 125]
    knn_paths = []
    meltron_paths = []
    for window_bins in bins:
        knn_path = tmp_path / f"knn_w{window_bins}.tsv"
        meltron_path = tmp_path / f"meltron_w{window_bins}.tsv"
        knn = pd.DataFrame([knn_row(window_bins, 1, 0, window_bins * 4)])
        if window_bins == 25:
            knn = knn.drop(columns=P_COLUMNS[-1])
        knn.to_csv(knn_path, sep="\t", index=False)
        pd.DataFrame([meltron_row(window_bins, 1, 0, window_bins * 4, 6.0)]).to_csv(
            meltron_path, sep="\t", index=False
        )
        knn_paths.append(knn_path)
        meltron_paths.append(meltron_path)

    outputs = [tmp_path / name for name in ("all.tsv", "regions.tsv", "regions.bed")]
    for output in outputs:
        output.write_text(f"previous {output.name}\n")

    result = run_cli(
        knn_paths,
        meltron_paths,
        bins,
        outputs[0],
        outputs[1],
        outputs[2],
        tmp_path / "build.log",
    )

    assert result.returncode != 0
    assert [path.read_text() for path in outputs] == [
        "previous all.tsv\n",
        "previous regions.tsv\n",
        "previous regions.bed\n",
    ]


def test_three_output_cli_accepts_all_empty_scales_and_writes_bed_without_header(tmp_path):
    bins = [25, 50, 75, 100, 125]
    knn_paths = []
    meltron_paths = []
    for window_bins in bins:
        knn_path = tmp_path / f"knn_w{window_bins}.tsv"
        meltron_path = tmp_path / f"meltron_w{window_bins}.tsv"
        empty_knn_frame(window_bins).to_csv(knn_path, sep="\t", index=False)
        empty_meltron_frame(window_bins).to_csv(meltron_path, sep="\t", index=False)
        knn_paths.append(knn_path)
        meltron_paths.append(meltron_path)

    out_all = tmp_path / "all.tsv"
    out_regions = tmp_path / "regions.tsv"
    out_bed = tmp_path / "regions.bed"
    result = run_cli(
        knn_paths,
        meltron_paths,
        bins,
        out_all,
        out_regions,
        out_bed,
        tmp_path / "build.log",
    )

    assert result.returncode == 0, result.stderr
    assert pd.read_csv(out_all, sep="\t").empty
    assert pd.read_csv(out_regions, sep="\t").empty
    assert out_bed.read_text() == ""

def test_three_output_cli_rejects_an_observed_scale_mismatch(tmp_path):
    bins = [25, 50, 75, 100, 125]
    knn_paths = []
    meltron_paths = []
    for index, window_bins in enumerate(bins):
        knn_path = tmp_path / f"knn_w{window_bins}.tsv"
        meltron_path = tmp_path / f"meltron_w{window_bins}.tsv"
        observed_bins = 50 if index == 0 else window_bins
        pd.DataFrame(
            [knn_row(observed_bins, 1, 0, observed_bins * 4)]
        ).to_csv(knn_path, sep="\t", index=False)
        pd.DataFrame(
            [meltron_row(observed_bins, 1, 0, observed_bins * 4, 6.0)]
        ).to_csv(meltron_path, sep="\t", index=False)
        knn_paths.append(knn_path)
        meltron_paths.append(meltron_path)

    result = run_cli(
        knn_paths,
        meltron_paths,
        bins,
        tmp_path / "all.tsv",
        tmp_path / "regions.tsv",
        tmp_path / "regions.bed",
        tmp_path / "build.log",
    )

    assert result.returncode != 0
    assert "reports r4000/w50, expected r4000/w25" in result.stderr
