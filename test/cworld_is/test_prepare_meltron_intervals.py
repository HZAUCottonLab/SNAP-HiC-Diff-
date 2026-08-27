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


def write_enriched(path, **coordinates):
    row = {
        "chrom": "HC04_A01",
        "window_id": 17,
        "comparison": "H1398_vs_L134",
        "knn_valid_k": 5,
        "resolution_bp": 4000,
        "window_bins": 25,
        "window_size_bp": 100000,
        "step_bins": 5,
        "step_size_bp": 20000,
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
            "--fdr-threshold",
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
        {
            "chrom": "HC04_A01",
            "start": 4000,
            "end": 8000,
            "ID": "r4000_w25_HC04_A01_17",
        }
    ]
    interval_map = pd.read_csv(
        output.with_name("H1398_vs_L134_interval_map.tsv"), sep="\t"
    )
    assert interval_map.loc[0, [
        "resolution_bp", "window_bins", "window_size_bp",
        "step_bins", "step_size_bp", "window_start", "window_end",
    ]].tolist() == [4000, 25, 100000, 5, 20000, 4000, 8000]


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
        {
            "chrom": "HC04_A01",
            "start": 4000,
            "end": 8000,
            "ID": "r4000_w25_HC04_A01_17",
        }
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


def test_preparer_writes_header_only_outputs_when_no_window_passes(tmp_path):
    enriched = tmp_path / "knn_enriched.tsv"
    output = tmp_path / "meltron" / "H1398_vs_L134_intervals.tsv"
    write_enriched(
        enriched,
        window_start=4000,
        window_end=8000,
        **{column: 1.0 for column in P_VALUE_COLUMNS},
    )

    result = run_preparer(enriched, output)

    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == ["chrom\tstart\tend\tID"]
    interval_map = output.with_name("H1398_vs_L134_interval_map.tsv")
    assert interval_map.read_text().splitlines() == [
        "resolution_bp\twindow_bins\twindow_size_bp\tstep_bins\tstep_size_bp"
        "\tchrom\tinterval_start\tinterval_end\twindow_start\twindow_end"
        "\twindow_id\tID\tn_significant"
    ]
