import csv
from pathlib import Path

import pytest

from scripts.validate_knn_background import REQUIRED_COLUMNS, validate_background


SCALE = {
    "resolution_bp": 4000,
    "window_bins": 25,
    "window_size_bp": 100000,
    "step_bins": 5,
    "step_size_bp": 20000,
}


def row(start=0, end=100000, comparison="a_vs_b"):
    values = {
        "comparison": comparison,
        "chrom": "HC04_A01",
        "window_start": str(start),
        "window_end": str(end),
    }
    for column in REQUIRED_COLUMNS:
        values.setdefault(column, "0.5")
    return values


def write_table(path: Path, rows, columns=REQUIRED_COLUMNS):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def test_validate_background_returns_scale_audit(tmp_path):
    path = tmp_path / "background.tsv"
    write_table(path, [row(), row(20000, 120000, "c_vs_d")])
    audit = validate_background(path, SCALE, knn_k=2)
    assert audit["row_count"] == 2
    assert audit["valid_knn_rows"] == 2
    assert audit["resolution_bp"] == 4000
    assert audit["window_bins"] == 25
    assert audit["step_bins"] == 5


def test_validate_background_accepts_the_existing_partial_tail_rule(tmp_path):
    path = tmp_path / "background.tsv"
    write_table(path, [row(20000, 67999)])
    assert validate_background(path, SCALE, knn_k=1)["row_count"] == 1


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        ([row(1, 100001)], "not on the 20000 bp step grid"),
        ([row(0, 40000)], "window length 40000 is outside"),
        ([row(), row()], "duplicate logical background row"),
        ([{**row(), "MSE": "nan"}], "valid KNN rows 0 is smaller than knn_k 1"),
    ],
)
def test_validate_background_rejects_invalid_rows(tmp_path, rows, message):
    path = tmp_path / "background.tsv"
    write_table(path, rows)
    with pytest.raises(ValueError, match=message):
        validate_background(path, SCALE, knn_k=1)


def test_validate_background_rejects_missing_columns(tmp_path):
    path = tmp_path / "background.tsv"
    write_table(path, [row()], columns=["chrom", "window_start", "window_end"])
    with pytest.raises(ValueError, match="background missing columns"):
        validate_background(path, SCALE, knn_k=1)
