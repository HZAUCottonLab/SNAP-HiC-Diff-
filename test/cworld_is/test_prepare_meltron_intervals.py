from pathlib import Path
import subprocess
import sys

import pandas as pd


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "prepare_meltron_intervals.py"
P_VALUE_COLUMNS = (
    "MSE_empirical_pvalue",
    "d_SSIM_empirical_pvalue",
    "d_Spearman_empirical_pvalue",
    "tw_score_empirical_pvalue",
    "d_NCC_empirical_pvalue",
)

SCALE_COLUMNS = (
    "resolution_bp",
    "window_bins",
    "window_size_bp",
    "step_bins",
    "step_size_bp",
)

INTERVAL_MAP_COLUMNS = [
    *SCALE_COLUMNS,
    "chrom",
    "interval_start",
    "interval_end",
    "window_start",
    "window_end",
    "window_id",
    "ID",
    "n_significant_metrics",
]


def write_enriched(path, **coordinates):
    row = {
        "resolution_bp": 4000,
        "window_bins": 25,
        "window_size_bp": 100000,
        "step_bins": 5,
        "step_size_bp": 20000,
        "chrom": "HC04_A01",
        "window_id": 17,
        "comparison": "H1398_vs_L134",
        "knn_valid_k": 5,
        **{column: 0.001 for column in P_VALUE_COLUMNS},
        **coordinates,
    }
    pd.DataFrame([row]).to_csv(path, sep="\t", index=False)


def run_preparer(enriched, output):
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--enriched",
            str(enriched),
            "--out",
            str(output),
            "--min-significant",
            "1",
            "--raw-pvalue-threshold",
            "0.05",
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def test_preparer_maps_knn_window_coordinates_to_canonical_intervals(tmp_path):
    enriched = tmp_path / "knn_enriched.tsv"
    output = tmp_path / "meltron" / "H1398_vs_L134_intervals.tsv"
    write_enriched(
        enriched,
        window_start=4000,
        window_end=8000,
    )

    result = run_preparer(enriched, output)

    assert result.returncode == 0, result.stderr
    intervals = pd.read_csv(output, sep="\t")
    assert list(intervals.columns) == ["chrom", "start", "end", "ID"]
    assert intervals.to_dict("records") == [
        {"chrom": "HC04_A01", "start": 4000, "end": 8000, "ID": "r4000_w25_HC04_A01_17"}
    ]


def test_preparer_preserves_legacy_canonical_coordinates(tmp_path):
    enriched = tmp_path / "knn_enriched.tsv"
    output = tmp_path / "meltron" / "H1398_vs_L134_intervals.tsv"
    write_enriched(
        enriched,
        start=3000,
        end=9000,
        window_start=4000,
        window_end=8000,
    )

    result = run_preparer(enriched, output)

    assert result.returncode == 0, result.stderr
    intervals = pd.read_csv(output, sep="\t")
    assert intervals.loc[0, ["start", "end"]].tolist() == [3000, 9000]
    interval_map = pd.read_csv(
        output.with_name("H1398_vs_L134_interval_map.tsv"), sep="\t"
    )
    assert interval_map.loc[0, ["interval_start", "interval_end"]].tolist() == [
        3000,
        9000,
    ]


def test_preparer_reports_coordinate_schema_when_both_forms_are_missing(tmp_path):
    enriched = tmp_path / "knn_enriched.tsv"
    output = tmp_path / "meltron" / "H1398_vs_L134_intervals.tsv"
    write_enriched(enriched)

    result = run_preparer(enriched, output)

    assert result.returncode != 0
    assert (
        "input coordinate schema requires either start/end or "
        "window_start/window_end" in result.stderr
    )


def test_preparer_accepts_legacy_canonical_coordinates_without_window_columns(tmp_path):
    enriched = tmp_path / "knn_enriched.tsv"
    output = tmp_path / "meltron" / "H1398_vs_L134_intervals.tsv"
    write_enriched(
        enriched,
        start=4000,
        end=8000,
    )

    result = run_preparer(enriched, output)

    assert result.returncode == 0, result.stderr
    intervals = pd.read_csv(output, sep="\t")
    assert intervals.to_dict("records") == [
        {"chrom": "HC04_A01", "start": 4000, "end": 8000, "ID": "r4000_w25_HC04_A01_17"}
    ]


def test_preparer_rejects_partial_coordinate_schema(tmp_path):
    enriched = tmp_path / "knn_enriched.tsv"
    output = tmp_path / "meltron" / "H1398_vs_L134_intervals.tsv"
    write_enriched(enriched, start=4000)

    result = run_preparer(enriched, output)

    assert result.returncode != 0
    assert (
        "input coordinate schema requires either start/end or "
        "window_start/window_end" in result.stderr
    )


def test_preparer_records_full_scale_and_source_identity_in_interval_map(tmp_path):
    enriched = tmp_path / "knn_enriched.tsv"
    output = tmp_path / "meltron" / "H1398_vs_L134_intervals.tsv"
    write_enriched(
        enriched,
        start=3000,
        end=9000,
        window_start=4000,
        window_end=8000,
    )

    result = run_preparer(enriched, output)

    assert result.returncode == 0, result.stderr
    interval_map = pd.read_csv(
        output.with_name("H1398_vs_L134_interval_map.tsv"), sep="\t"
    )
    assert list(interval_map.columns) == INTERVAL_MAP_COLUMNS
    assert interval_map.to_dict("records") == [
        {
            "resolution_bp": 4000,
            "window_bins": 25,
            "window_size_bp": 100000,
            "step_bins": 5,
            "step_size_bp": 20000,
            "chrom": "HC04_A01",
            "interval_start": 3000,
            "interval_end": 9000,
            "window_start": 4000,
            "window_end": 8000,
            "window_id": 17,
            "ID": "r4000_w25_HC04_A01_17",
            "n_significant_metrics": 5,
        }
    ]


def test_preparer_writes_stable_empty_interval_map_at_raw_threshold(tmp_path):
    enriched = tmp_path / "knn_enriched.tsv"
    output = tmp_path / "meltron" / "H1398_vs_L134_intervals.tsv"
    pvalues = {column: 0.05 for column in P_VALUE_COLUMNS}
    write_enriched(
        enriched,
        window_start=4000,
        window_end=8000,
        **pvalues,
    )

    result = run_preparer(enriched, output)

    assert result.returncode == 0, result.stderr
    intervals = pd.read_csv(output, sep="\t")
    interval_map = pd.read_csv(output.with_name("H1398_vs_L134_interval_map.tsv"), sep="\t")
    assert list(intervals.columns) == ["chrom", "start", "end", "ID"]
    assert intervals.empty
    assert list(interval_map.columns) == INTERVAL_MAP_COLUMNS
    assert interval_map.empty


def test_preparer_requires_exactly_the_canonical_raw_pvalue_columns(tmp_path):
    enriched = tmp_path / "knn_enriched.tsv"
    output = tmp_path / "meltron" / "H1398_vs_L134_intervals.tsv"
    write_enriched(
        enriched,
        window_start=4000,
        window_end=8000,
        extra_empirical_pvalue=0.001,
    )

    result = run_preparer(enriched, output)

    assert result.returncode != 0


def test_preparer_selects_raw_pvalue_candidate_despite_stale_fdr_and_bh_columns(
    tmp_path,
):
    enriched = tmp_path / "knn_enriched.tsv"
    output = tmp_path / "meltron" / "H1398_vs_L134_intervals.tsv"
    raw_pvalues = {column: 0.5 for column in P_VALUE_COLUMNS}
    raw_pvalues["MSE_empirical_pvalue"] = 0.001
    write_enriched(
        enriched,
        window_start=4000,
        window_end=8000,
        MSE_pvalue_BH_FDR=1.0,
        d_SSIM_pvalue_BH_FDR=1.0,
        stale_fdr=1.0,
        stale_BH=1.0,
        **raw_pvalues,
    )

    result = run_preparer(enriched, output)

    assert result.returncode == 0, result.stderr
    intervals = pd.read_csv(output, sep="\t")
    interval_map = pd.read_csv(
        output.with_name("H1398_vs_L134_interval_map.tsv"), sep="\t"
    )
    assert intervals["ID"].tolist() == ["r4000_w25_HC04_A01_17"]
    assert interval_map["n_significant_metrics"].tolist() == [1]


def test_preparer_rejects_input_missing_a_canonical_raw_pvalue_column(tmp_path):
    enriched = tmp_path / "knn_enriched.tsv"
    output = tmp_path / "meltron" / "H1398_vs_L134_intervals.tsv"
    write_enriched(enriched, window_start=4000, window_end=8000)
    table = pd.read_csv(enriched, sep="\t").drop(columns=["d_NCC_empirical_pvalue"])
    table.to_csv(enriched, sep="\t", index=False)

    result = run_preparer(enriched, output)

    assert result.returncode != 0
    assert "input requires exactly the canonical empirical p-value columns" in result.stderr
