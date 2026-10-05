from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "merge_meltron_scores.py"
KEY_COLUMNS = ["chrom", "start", "end", "ID"]
MELTRON_ID = "r4000_w25_HC04_A01_17"
SOURCE_COLUMNS = [
    "resolution_bp",
    "window_bins",
    "window_size_bp",
    "step_bins",
    "step_size_bp",
    "window_start",
    "window_end",
    "window_id",
]
EMPTY_OUTPUT_COLUMNS = [
    *KEY_COLUMNS,
    *SOURCE_COLUMNS,
    "meltron_melting_score",
    "meltron_compaction_score",
    "meltron_best_score",
    "meltron_best_direction",
    "meltron_is_significant",
    "n_significant_metrics",
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


def interval_map_record(
    *,
    window_id=17,
    window_start=4000,
    window_end=8000,
    metric_column="n_significant_metrics",
    n_significant_metrics=3,
):
    return {
        "ID": f"r4000_w25_HC04_A01_{window_id}",
        "resolution_bp": 4000,
        "window_bins": 25,
        "window_size_bp": 100000,
        "step_bins": 5,
        "step_size_bp": 20000,
        "window_start": window_start,
        "window_end": window_end,
        "window_id": window_id,
        metric_column: n_significant_metrics,
    }


def write_interval_map(path):
    write_tsv(path, [interval_map_record()])


def write_empty_directional_results(decondensation, condensation):
    write_tsv(decondensation, [], columns=[*KEY_COLUMNS, "melting_score"])
    write_tsv(condensation, [], columns=[*KEY_COLUMNS, "compaction_score"])


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
                "ID": MELTRON_ID,
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
                "ID": MELTRON_ID,
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


def test_merge_restores_full_source_window_identity(tmp_path):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_tsv(
        decondensation,
        [{
            "chrom": "HC04_A01",
            "start": 4000,
            "end": 8000,
            "ID": MELTRON_ID,
            "melting_score": 6.0,
        }],
    )
    write_tsv(
        condensation,
        [{
            "chrom": "HC04_A01",
            "start": 4000,
            "end": 8000,
            "ID": MELTRON_ID,
            "compaction_score": 8.0,
        }],
    )
    write_interval_map(interval_map)

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode == 0, result.stderr
    row = pd.read_csv(output, sep="\t").iloc[0]
    assert row[[
        "resolution_bp", "window_bins", "window_size_bp",
        "step_bins", "step_size_bp", "window_start", "window_end", "window_id",
    ]].tolist() == [4000, 25, 100000, 5, 20000, 4000, 8000, 17]


def test_merge_writes_header_only_output_for_empty_directional_results(tmp_path):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_empty_directional_results(decondensation, condensation)
    write_tsv(interval_map, [], columns=["ID", *SOURCE_COLUMNS, "n_significant"])

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == ["\t".join(EMPTY_OUTPUT_COLUMNS)]


@pytest.mark.parametrize(
    ("direction", "score_column"),
    [("decondensation", "melting_score"), ("condensation", "compaction_score")],
)
def test_merge_rejects_directional_ids_missing_from_interval_map(
    tmp_path, direction, score_column
):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    unmatched_id = "r4000_w25_HC04_A01_18"
    directional_record = {
        "chrom": "HC04_A01",
        "start": 4000,
        "end": 8000,
        "ID": unmatched_id,
        score_column: 8.0,
    }
    if direction == "decondensation":
        write_tsv(decondensation, [directional_record])
        write_tsv(condensation, [], columns=[*KEY_COLUMNS, "compaction_score"])
    else:
        write_tsv(decondensation, [], columns=[*KEY_COLUMNS, "melting_score"])
        write_tsv(condensation, [directional_record])
    write_interval_map(interval_map)

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode != 0
    assert (
        f"interval map missing source identity for Meltron IDs: {unmatched_id}"
        in result.stderr
    )


def test_merge_rejects_interval_map_with_missing_id(tmp_path):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_empty_directional_results(decondensation, condensation)
    record = interval_map_record()
    record["ID"] = None
    write_tsv(interval_map, [record])

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode != 0
    assert "interval map ID must be non-missing and unique" in result.stderr


def test_merge_rejects_interval_map_with_duplicate_id(tmp_path):
    decondensation = tmp_path / "decondensation.tsv"
    condensation = tmp_path / "condensation.tsv"
    interval_map = tmp_path / "interval_map.tsv"
    output = tmp_path / "merged.tsv"
    write_empty_directional_results(decondensation, condensation)
    record = interval_map_record()
    write_tsv(interval_map, [record, record.copy()])

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode != 0
    assert "interval map ID must be non-missing and unique" in result.stderr


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
                "ID": MELTRON_ID,
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
        [{"chrom": "HC04_A01", "start": 4000, "end": 8000, "ID": MELTRON_ID}],
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
                "ID": "r4000_w25_HC04_A01_1",
                "melting_score": 8.0,
            },
            {
                "chrom": "HC04_A01",
                "start": 8000,
                "end": 12000,
                "ID": "r4000_w25_HC04_A01_3",
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
                "ID": "r4000_w25_HC04_A01_2",
                "compaction_score": 5.0,
            },
            {
                "chrom": "HC04_A01",
                "start": 8000,
                "end": 12000,
                "ID": "r4000_w25_HC04_A01_3",
                "compaction_score": 9.0,
            },
        ],
    )
    write_tsv(
        interval_map,
        [
            interval_map_record(
                window_id=1,
                window_start=0,
                window_end=4000,
                n_significant_metrics=1,
            ),
            interval_map_record(
                window_id=2,
                window_start=4000,
                window_end=8000,
                n_significant_metrics=1,
            ),
            interval_map_record(
                window_id=3,
                window_start=8000,
                window_end=12000,
                n_significant_metrics=1,
            ),
        ],
    )

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode == 0, result.stderr
    rows = pd.read_csv(output, sep="\t").set_index("ID")
    assert rows.loc["r4000_w25_HC04_A01_1", "meltron_melting_score"] == 8.0
    assert pd.isna(rows.loc["r4000_w25_HC04_A01_1", "meltron_compaction_score"])
    assert rows.loc["r4000_w25_HC04_A01_1", "meltron_best_score"] == 8.0
    assert rows.loc["r4000_w25_HC04_A01_1", "meltron_best_direction"] == "decondensation"
    assert pd.isna(rows.loc["r4000_w25_HC04_A01_2", "meltron_melting_score"])
    assert rows.loc["r4000_w25_HC04_A01_2", "meltron_compaction_score"] == 5.0
    assert rows.loc["r4000_w25_HC04_A01_2", "meltron_best_score"] == 5.0
    assert rows.loc["r4000_w25_HC04_A01_2", "meltron_best_direction"] == "condensation"
    assert rows.loc["r4000_w25_HC04_A01_3", "meltron_best_score"] == 9.0
    assert rows.loc["r4000_w25_HC04_A01_3", "meltron_best_direction"] == "condensation"


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
                "ID": MELTRON_ID,
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
                "ID": MELTRON_ID,
                "compaction_score": 8.0,
            }
        ],
    )
    write_tsv(
        interval_map,
        [interval_map_record(metric_column=metric_column)],
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
                "ID": MELTRON_ID,
                "compaction_score": 8.0,
            }
        ],
    )
    record = interval_map_record()
    record.pop("n_significant_metrics")
    write_tsv(interval_map, [record])

    result = run_merge(decondensation, condensation, interval_map, output)

    assert result.returncode != 0
    assert (
        "interval map requires n_significant or n_significant_metrics"
        in result.stderr
    )
