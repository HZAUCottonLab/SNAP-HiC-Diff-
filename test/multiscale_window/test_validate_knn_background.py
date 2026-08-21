import json
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "validate_knn_background.py"
sys.path.insert(0, str(REPO / "scripts"))

from validate_knn_background import validate_background_table


REQUIRED = [
    "comparison", "chrom", "window_start", "window_end",
    "EC", "IB", "SP", "MSE", "SSIM", "Spearman", "tw_pvalue", "NCC",
]


def background_row(comparison, start, end):
    return {
        "comparison": comparison, "chrom": "chr1",
        "window_start": start, "window_end": end,
        "EC": 1.0, "IB": 2.0, "SP": 3.0, "MSE": 4.0,
        "SSIM": 0.8, "Spearman": 0.7, "tw_pvalue": 0.1, "NCC": 0.6,
    }


def write_background(path, rows, columns=None):
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def validate(path, knn_k=1):
    return validate_background_table(
        path, resolution_bp=4000, window_bins=25, step_bins=5, knn_k=knn_k
    )


def test_allows_matching_coordinates_for_distinct_comparisons(tmp_path):
    background = tmp_path / "background.tsv"
    write_background(
        background,
        [
            background_row("rep1", 0, 100000),
            background_row("rep2", 0, 100000),
        ],
    )

    audit = validate(background, knn_k=2)

    assert audit["total_rows"] == 2
    assert audit["comparisons"] == 2


def test_rejects_duplicate_comparison_chrom_window_coordinates(tmp_path):
    background = tmp_path / "background.tsv"
    row = background_row("rep1", 0, 100000)
    write_background(background, [row, row])

    with pytest.raises(ValueError, match="duplicate comparison/chrom/window coordinates"):
        validate(background, knn_k=2)


@pytest.mark.parametrize("column", REQUIRED)
def test_reports_each_missing_required_column(tmp_path, column):
    background = tmp_path / "background.tsv"
    row = background_row("rep1", 0, 100000)
    row.pop(column)
    write_background(background, [row])

    with pytest.raises(ValueError, match=f"background missing columns: {column}"):
        validate(background)


def test_rejects_background_with_fewer_finite_rows_than_knn_k(tmp_path):
    background = tmp_path / "background.tsv"
    invalid = background_row("rep1", 20000, 120000)
    invalid["NCC"] = "not-a-number"
    write_background(
        background,
        [
            background_row("rep1", 0, 100000),
            invalid,
            background_row("rep1", 40000, 140000),
        ],
    )

    with pytest.raises(ValueError, match="2 valid rows, fewer than knn_k=3"):
        validate(background, knn_k=3)


def test_accepts_complete_r4000_w25_windows_and_reports_expected_audit(tmp_path):
    background = tmp_path / "background.tsv"
    write_background(
        background,
        [
            background_row("rep1", 0, 100000),
            background_row("rep1", 20000, 120000),
        ],
    )

    audit = validate(background, knn_k=2)

    assert audit == {
        "path": str(background.resolve()),
        "resolution_bp": 4000,
        "window_bins": 25,
        "window_size_bp": 100000,
        "step_bins": 5,
        "step_size_bp": 20000,
        "total_rows": 2,
        "valid_rows": 2,
        "complete_rows": 2,
        "terminal_partial_rows": 0,
        "comparisons": 1,
        "chromosomes": 1,
    }


def test_rejects_nonterminal_96kb_span(tmp_path):
    background = tmp_path / "background.tsv"
    write_background(
        background,
        [
            background_row("rep1", 0, 96000),
            background_row("rep1", 20000, 120000),
        ],
    )

    with pytest.raises(ValueError, match="partial background window must be terminal"):
        validate(background, knn_k=2)


def test_allows_final_52kb_partial_window_only_at_end_of_group(tmp_path):
    valid_background = tmp_path / "valid-background.tsv"
    write_background(
        valid_background,
        [
            background_row("rep1", 0, 100000),
            background_row("rep1", 20000, 120000),
            background_row("rep1", 40000, 92000),
        ],
    )

    audit = validate(valid_background, knn_k=3)

    assert audit["complete_rows"] == 2
    assert audit["terminal_partial_rows"] == 1

    invalid_background = tmp_path / "invalid-background.tsv"
    write_background(
        invalid_background,
        [
            background_row("rep1", 0, 52000),
            background_row("rep1", 20000, 120000),
        ],
    )

    with pytest.raises(ValueError, match="partial background window must be terminal"):
        validate(invalid_background, knn_k=2)


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        ([background_row("rep1", 2000, 102000)], "starts must align to 4000 bp"),
        (
            [
                background_row("rep1", 0, 100000),
                background_row("rep1", 24000, 124000),
            ],
            "start grid mismatch",
        ),
    ],
)
def test_rejects_off_grid_starts_and_steps(tmp_path, rows, message):
    background = tmp_path / "background.tsv"
    write_background(background, rows)

    with pytest.raises(ValueError, match=message):
        validate(background, knn_k=len(rows))


def test_cli_failure_preserves_existing_audit(tmp_path):
    background = tmp_path / "background.tsv"
    row = background_row("rep1", 0, 100000)
    row.pop("EC")
    write_background(background, [row])
    audit = tmp_path / "audit.json"
    original_audit = '{"previous": true}\n'
    audit.write_text(original_audit)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--background", str(background),
            "--resolution-bp", "4000",
            "--window-bins", "25",
            "--step-bins", "5",
            "--knn-k", "1",
            "--audit", str(audit),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode != 0
    assert audit.read_text() == original_audit
    assert json.loads(audit.read_text()) == {"previous": True}
