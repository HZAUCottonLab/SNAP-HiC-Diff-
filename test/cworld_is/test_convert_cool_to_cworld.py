from pathlib import Path
import gzip
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
import convert_cool_to_cworld as converter


def test_rejects_wrong_bin_size_before_loading_matrix(monkeypatch, tmp_path):
    fake = SimpleNamespace(chromnames=["HC04_A01"], info={"bin-size": 20000})
    monkeypatch.setattr(converter.cooler, "Cooler", lambda _: fake)

    with pytest.raises(ValueError, match="expected 4000 bp, found 20000 bp"):
        converter.convert_cool_to_cworld(
            "wrong.cool",
            "HC04_A01",
            tmp_path / "wrong.matrix.gz",
            expected_bin_size=4000,
        )


def test_writes_cworld_header_and_labels(monkeypatch, tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["HC04_A01", "HC04_A01"],
            "start": [0, 4000],
            "end": [4000, 8000],
        }
    )

    class FakeMatrix:
        def fetch(self, chrom):
            assert chrom == "HC04_A01"
            return np.array([[1.0, np.nan], [2.0, 3.0]])

    fake = SimpleNamespace(
        chromnames=["HC04_A01"],
        info={"bin-size": 4000},
        bins=lambda: bins,
        matrix=lambda balance: FakeMatrix(),
    )
    monkeypatch.setattr(converter.cooler, "Cooler", lambda _: fake)
    output = tmp_path / "nested" / "matrix.gz"

    converter.convert_cool_to_cworld(
        "ok.cool", "HC04_A01", output, expected_bin_size=4000
    )

    with gzip.open(output, "rt") as handle:
        lines = handle.read().splitlines()
    assert lines[0] == "\tHiC_1|aaa|HC04_A01:0-4000\tHiC_2|aaa|HC04_A01:4000-8000"
    assert lines[1] == "HiC_1|aaa|HC04_A01:0-4000\t1\t0"
    assert lines[2] == "HiC_2|aaa|HC04_A01:4000-8000\t2\t3"
    assert output.parent.is_dir()


def test_cli_reports_bin_size_mismatch(monkeypatch, capsys, tmp_path):
    fake = SimpleNamespace(chromnames=["HC04_A01"], info={"bin-size": 20000})
    monkeypatch.setattr(converter.cooler, "Cooler", lambda _: fake)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "convert_cool_to_cworld.py",
            "--cool", "wrong.cool",
            "--chrom", "HC04_A01",
            "--out", str(tmp_path / "out.matrix.gz"),
            "--expected-bin-size", "4000",
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        converter.run()

    assert exc_info.value.code == 2
    assert "expected 4000 bp, found 20000 bp" in capsys.readouterr().err
