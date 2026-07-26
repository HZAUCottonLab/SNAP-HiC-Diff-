from pathlib import Path
import subprocess
import sys

import pandas as pd


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "convert_cworld_IS_to_meltron_IS.py"


def write_bedgraph(root, sample, chrom, distance):
    path = (
        root
        / "cworld_insulation"
        / sample
        / f"{sample}__{chrom}--is{distance}--nt0--ids16000--ss0--immean.insulation.bedGraph"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "track type=bedGraph name=IS\n"
        f"{chrom}\t0\t16000\t0.25\n"
        f"{chrom}\t16000\t32000\t0.75\n"
    )


def run_converter(root, chroms):
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--sample", "S1",
            "--chroms", *chroms,
            "--outdir", str(root),
            "--is-distances", "32000",
            "--log", str(root / "converter.log"),
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def test_converter_writes_complete_meltron_tsv(tmp_path):
    write_bedgraph(tmp_path, "S1", "HC04_A01", 32000)

    result = run_converter(tmp_path, ["HC04_A01"])

    assert result.returncode == 0, result.stderr
    result_tsv = tmp_path / "IS" / "S1" / "S1_IS.tsv"
    frame = pd.read_csv(result_tsv, sep="\t")
    assert list(frame.columns) == ["chrom", "viewpoint", "IS_distance", "IS"]
    assert frame.to_dict("records") == [
        {"chrom": "HC04_A01", "viewpoint": 8000, "IS_distance": 32000, "IS": 0.25},
        {"chrom": "HC04_A01", "viewpoint": 24000, "IS_distance": 32000, "IS": 0.75},
    ]


def test_converter_rejects_missing_chromosome_distance(tmp_path):
    write_bedgraph(tmp_path, "S1", "HC04_A01", 32000)

    result = run_converter(tmp_path, ["HC04_A01", "HC04_A02"])

    assert result.returncode != 0
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()
