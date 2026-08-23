from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "merge_meltron_scores.py"
KEY_COLUMNS = ["chrom", "start", "end", "ID"]

SCALE_COLUMNS = [
    "resolution_bp",
    "window_bins",
    "window_size_bp",
    "step_bins",
    "step_size_bp",
]

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

MERGED_COLUMNS = [
    *INTERVAL_MAP_COLUMNS,
    "meltron_melting_score",
    "meltron_compaction_score",
    "meltron_best_score",
    "meltron_best_direction",
    "meltron_is_significant",
    "perc_na_or_missing_ref_decondensation",
    "perc_na_or_missing_inp_decondensation",
    "perc_na_or_missing_ref_condensation",
    "perc_na_or_missing_inp_condensation",
]


def write_tsv(path, records, columns=None):
    pd.DataFrame(records, columns=columns).to_csv(path, sep="\t", index=False)


def run_merge(decondensation, condensation, interval_map, output):
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--comparison",
            "H1398_vs_L134",
            "--decondensation",
            str(decondensation),
            "--condensation",
            str(condensation),
            "--interval-map",
            str(interval_map),
            "--out",
            str(output),
            "--score-threshold",
            "7",
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def write_interval_map(path):
    write_tsv(
        path,
        [
            {
                "resolution_bp": 4000,
                "window_bins": 25,
                "window_size_bp": 100000,
                "step_bins": 5,
                "step_size_bp": 20000,
                "chrom": "HC04_A01",
                "interval_start": 4000,
                "interval_end": 8000,
                "window_start": 4000,
                "window_end": 8000,
                "window_id": 17,
                "ID": "HC04_A01_17",
                "n_significant_metrics": 3,
            }
        ],
    )


def test_merge_normalizes_direction_specific_meltron_scores(tmp_path):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_tsv(
        decondensation,
        [
            {
                "chrom": "HC04_A01",
                "start": 4000,
                "end": 8000,
                "ID": "HC04_A01_17",
                "melting_score": 6.0,
                "perc_na_or_missing_ref": 0.1,
                "perc_na_or_missing_inp": 0.2,
            }
        ],
    )
    write_tsv(
        condensation,
        [
            {
                "chrom": "HC04_A01",
                "start": 4000,
                "end": 8000,
                "ID": "HC04_A01_17",
                "compaction_score": 8.0,
                "perc_na_or_missing_ref": 0.3,
                "perc_na_or_missing_inp": 0.4,
            }
        ],
    )
    write_interval_map(interval_map)

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode == 0, result.stderr
    merged = pd.read_csv(output, sep="\t")
    row = merged.iloc[0]
    assert row["meltron_melting_score"] == 6.0
    assert row["meltron_compaction_score"] == 8.0
    assert row["meltron_best_score"] == 8.0
    assert row["meltron_best_direction"] == "condensation"
    assert row["meltron_is_significant"] == 1
    assert row[
        [
            "perc_na_or_missing_ref_decondensation",
            "perc_na_or_missing_inp_decondensation",
            "perc_na_or_missing_ref_condensation",
            "perc_na_or_missing_inp_condensation",
        ]
    ].tolist() == [0.1, 0.2, 0.3, 0.4]
    assert row["window_id"] == 17
    assert row["n_significant_metrics"] == 3


def test_merge_uses_compaction_score_when_only_condensation_is_available(tmp_path):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_tsv(decondensation, [], columns=[*KEY_COLUMNS, "melting_score"])
    write_tsv(
        condensation,
        [
            {
                "chrom": "HC04_A01",
                "start": 4000,
                "end": 8000,
                "ID": "HC04_A01_17",
                "compaction_score": 4.0,
                "perc_na_or_missing_ref": 0.3,
                "perc_na_or_missing_inp": 0.4,
            }
        ],
    )
    write_interval_map(interval_map)

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode == 0, result.stderr
    row = pd.read_csv(output, sep="\t").iloc[0]
    assert pd.isna(row["meltron_melting_score"])
    assert row["meltron_compaction_score"] == 4.0
    assert row["meltron_best_score"] == 4.0
    assert row["meltron_best_direction"] == "condensation"
    assert row["meltron_is_significant"] == 0
    assert row["perc_na_or_missing_ref_condensation"] == 0.3
    assert row["perc_na_or_missing_inp_condensation"] == 0.4


def test_merge_reports_missing_directional_score_schema(tmp_path):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_tsv(
        decondensation,
        [{"chrom": "HC04_A01", "start": 4000, "end": 8000, "ID": "HC04_A01_17"}],
    )
    write_tsv(condensation, [], columns=[*KEY_COLUMNS, "compaction_score"])
    write_interval_map(interval_map)

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode != 0
    assert "decondensation result requires melting_score" in result.stderr


def test_merge_assigns_outer_join_rows_to_the_direction_with_a_score(tmp_path):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_tsv(
        decondensation,
        [
            {
                "chrom": "HC04_A01",
                "start": 0,
                "end": 4000,
                "ID": "HC04_A01_A",
                "melting_score": 8.0,
            },
            {
                "chrom": "HC04_A01",
                "start": 8000,
                "end": 12000,
                "ID": "HC04_A01_C",
                "melting_score": 6.0,
            },
        ],
    )
    write_tsv(
        condensation,
        [
            {
                "chrom": "HC04_A01",
                "start": 4000,
                "end": 8000,
                "ID": "HC04_A01_B",
                "compaction_score": 5.0,
            },
            {
                "chrom": "HC04_A01",
                "start": 8000,
                "end": 12000,
                "ID": "HC04_A01_C",
                "compaction_score": 9.0,
            },
        ],
    )
    write_tsv(
        interval_map,
        [
            {"ID": "HC04_A01_A", "window_id": 1, "n_significant_metrics": 1},
            {"ID": "HC04_A01_B", "window_id": 2, "n_significant_metrics": 1},
            {"ID": "HC04_A01_C", "window_id": 3, "n_significant_metrics": 1},
        ],
    )

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode == 0, result.stderr
    rows = pd.read_csv(output, sep="\t").set_index("ID")
    assert rows.loc["HC04_A01_A", "meltron_melting_score"] == 8.0
    assert pd.isna(rows.loc["HC04_A01_A", "meltron_compaction_score"])
    assert rows.loc["HC04_A01_A", "meltron_best_score"] == 8.0
    assert rows.loc["HC04_A01_A", "meltron_best_direction"] == "decondensation"
    assert pd.isna(rows.loc["HC04_A01_B", "meltron_melting_score"])
    assert rows.loc["HC04_A01_B", "meltron_compaction_score"] == 5.0
    assert rows.loc["HC04_A01_B", "meltron_best_score"] == 5.0
    assert rows.loc["HC04_A01_B", "meltron_best_direction"] == "condensation"
    assert rows.loc["HC04_A01_C", "meltron_best_score"] == 9.0
    assert rows.loc["HC04_A01_C", "meltron_best_direction"] == "condensation"


def test_merge_retains_na_metrics_when_only_decondensation_is_available(tmp_path):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_tsv(
        decondensation,
        [
            {
                "chrom": "HC04_A01",
                "start": 4000,
                "end": 8000,
                "ID": "HC04_A01_17",
                "melting_score": 4.0,
                "perc_na_or_missing_ref": 0.1,
                "perc_na_or_missing_inp": 0.2,
            }
        ],
    )
    write_tsv(condensation, [], columns=[*KEY_COLUMNS, "compaction_score"])
    write_interval_map(interval_map)

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode == 0, result.stderr
    row = pd.read_csv(output, sep="\t").iloc[0]
    assert row["meltron_melting_score"] == 4.0
    assert pd.isna(row["meltron_compaction_score"])
    assert row["meltron_best_score"] == 4.0
    assert row["meltron_best_direction"] == "decondensation"
    assert row["perc_na_or_missing_ref_decondensation"] == 0.1


@pytest.mark.parametrize(
    "metric_column", ["n_significant", "n_significant_metrics"]
)
def test_merge_preserves_preparer_significant_metric_count(tmp_path, metric_column):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_tsv(decondensation, [], columns=[*KEY_COLUMNS, "melting_score"])
    write_tsv(
        condensation,
        [
            {
                "chrom": "HC04_A01",
                "start": 4000,
                "end": 8000,
                "ID": "HC04_A01_17",
                "compaction_score": 8.0,
            }
        ],
    )
    write_tsv(
        interval_map,
        [
            {
                "chrom": "HC04_A01",
                "interval_start": 4000,
                "interval_end": 8000,
                "window_start": 4000,
                "window_end": 8000,
                "window_id": 17,
                "ID": "HC04_A01_17",
                metric_column: 3,
            }
        ],
    )

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode == 0, result.stderr
    row = pd.read_csv(output, sep="\t").iloc[0]
    assert row["window_id"] == 17
    assert row["n_significant_metrics"] == 3


def test_merge_reports_missing_interval_metric_count_schema(tmp_path):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_tsv(decondensation, [], columns=[*KEY_COLUMNS, "melting_score"])
    write_tsv(
        condensation,
        [
            {
                "chrom": "HC04_A01",
                "start": 4000,
                "end": 8000,
                "ID": "HC04_A01_17",
                "compaction_score": 8.0,
            }
        ],
    )
    write_tsv(
        interval_map,
        [{"ID": "HC04_A01_17", "window_id": 17}],
    )

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode != 0
    assert (
        "interval map requires n_significant or n_significant_metrics"
        in result.stderr

    )

def test_merge_carries_full_scale_and_source_schema_from_interval_map(tmp_path):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_tsv(decondensation, [{"chrom": "HC04_A01", "start": 4000, "end": 8000, "ID": "HC04_A01_17", "melting_score": 6.0}])
    write_tsv(condensation, [{"chrom": "HC04_A01", "start": 4000, "end": 8000, "ID": "HC04_A01_17", "compaction_score": 8.0}])
    write_interval_map(interval_map)

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode == 0, result.stderr
    merged = pd.read_csv(output, sep="\t")
    assert list(merged.columns) == MERGED_COLUMNS
    row = merged.iloc[0]
    assert row[SCALE_COLUMNS].tolist() == [4000, 25, 100000, 5, 20000]
    assert row[["interval_start", "interval_end", "window_start", "window_end"]].tolist() == [4000, 8000, 4000, 8000]


def test_merge_writes_full_stable_schema_when_both_directions_are_empty(tmp_path):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_tsv(decondensation, [], columns=[*KEY_COLUMNS, "melting_score"])
    write_tsv(condensation, [], columns=[*KEY_COLUMNS, "compaction_score"])
    write_interval_map(interval_map)

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode == 0, result.stderr
    merged = pd.read_csv(output, sep="\t")
    assert list(merged.columns) == MERGED_COLUMNS
    assert merged.empty


def test_merge_rejects_duplicate_interval_map_ids_with_one_to_one_validation(tmp_path):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_tsv(decondensation, [], columns=[*KEY_COLUMNS, "melting_score"])
    write_tsv(
        condensation,
        [
            {
                "chrom": "HC04_A01",
                "start": 4000,
                "end": 8000,
                "ID": "HC04_A01_17",
                "compaction_score": 8.0,
            }
        ],
    )
    write_interval_map(interval_map)
    duplicated_map = pd.read_csv(interval_map, sep="\t")
    pd.concat([duplicated_map, duplicated_map], ignore_index=True).to_csv(
        interval_map, sep="\t", index=False
    )

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode != 0
    assert "not a one-to-one merge" in result.stderr
