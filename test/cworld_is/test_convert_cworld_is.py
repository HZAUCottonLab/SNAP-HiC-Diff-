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
        / f"{sample}__{chrom}--is{distance + 1}--nt0--ids{distance // 2 + 1}--ss1--immean.insulation.bedGraph"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "track type=bedGraph name=IS\n"
        f"{chrom}\t0\t16000\t0.25\n"
        f"{chrom}\t16000\t32000\t0.75\n"
    )
    return path


def write_manifest(root, sample, chrom, entries):
    path = root / "cworld_insulation" / sample / f"{sample}__{chrom}.is_manifest.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{distance}\t{filename}\n" for distance, filename in entries))
    return path


def run_converter(root, chroms, manifests, distances="32000"):
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--sample", "S1",
            "--chroms", *chroms,
            "--manifests", *(str(path) for path in manifests),
            "--outdir", str(root),
            "--is-distances", distances,
            "--log", str(root / "converter.log"),
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def test_converter_writes_complete_meltron_tsv(tmp_path):
    bedgraph = write_bedgraph(tmp_path, "S1", "HC04_A01", 32000)
    manifest = write_manifest(tmp_path, "S1", "HC04_A01", [(32000, bedgraph.name)])

    result = run_converter(tmp_path, ["HC04_A01"], [manifest])

    assert result.returncode == 0, result.stderr
    result_tsv = tmp_path / "IS" / "S1" / "S1_IS.tsv"
    frame = pd.read_csv(result_tsv, sep="\t")
    assert list(frame.columns) == ["chrom", "viewpoint", "IS_distance", "IS"]
    assert frame.to_dict("records") == [
        {"chrom": "HC04_A01", "viewpoint": 8000, "IS_distance": 32000, "IS": 0.25},
        {"chrom": "HC04_A01", "viewpoint": 24000, "IS_distance": 32000, "IS": 0.75},
    ]


def test_converter_preserves_manifest_distance_despite_cworld_suffix(tmp_path):
    bedgraph = write_bedgraph(tmp_path, "S1", "HC04_A01", 64000)
    manifest = write_manifest(tmp_path, "S1", "HC04_A01", [(32000, bedgraph.name)])

    result = run_converter(tmp_path, ["HC04_A01"], [manifest])

    assert result.returncode == 0, result.stderr
    frame = pd.read_csv(tmp_path / "IS" / "S1" / "S1_IS.tsv", sep="\t")
    assert set(frame["IS_distance"]) == {32000}


def test_converter_rejects_duplicate_distance_in_manifest(tmp_path):
    bedgraph = write_bedgraph(tmp_path, "S1", "HC04_A01", 32000)
    manifest = write_manifest(
        tmp_path,
        "S1",
        "HC04_A01",
        [(32000, bedgraph.name), (32000, bedgraph.name)],
    )

    result = run_converter(tmp_path, ["HC04_A01"], [manifest])

    assert result.returncode != 0
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()


def test_converter_rejects_missing_manifest_target(tmp_path):
    manifest = write_manifest(
        tmp_path,
        "S1",
        "HC04_A01",
        [(32000, "missing.bedGraph")],
    )

    result = run_converter(tmp_path, ["HC04_A01"], [manifest])

    assert result.returncode != 0
    assert "manifest target missing" in result.stderr
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()


def test_converter_rejects_unequal_chromosome_and_manifest_lists(tmp_path):
    bedgraph = write_bedgraph(tmp_path, "S1", "HC04_A01", 32000)
    manifest = write_manifest(tmp_path, "S1", "HC04_A01", [(32000, bedgraph.name)])

    result = run_converter(tmp_path, ["HC04_A01", "HC04_A02"], [manifest])

    assert result.returncode != 0
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()


def test_converter_never_falls_back_to_legacy_bedgraph_glob(tmp_path):
    write_bedgraph(tmp_path, "S1", "HC04_A01", 32000)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--sample", "S1",
            "--chroms", "HC04_A01",
            "--outdir", str(tmp_path),
            "--is-distances", "32000",
            "--log", str(tmp_path / "converter.log"),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode != 0
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()


def test_converter_rejects_duplicate_manifest_paths(tmp_path):
    bedgraph = write_bedgraph(tmp_path, "S1", "HC04_A01", 32000)
    manifest = write_manifest(tmp_path, "S1", "HC04_A01", [(32000, bedgraph.name)])

    result = run_converter(
        tmp_path,
        ["HC04_A01", "HC04_A02"],
        [manifest, manifest],
    )

    assert result.returncode != 0
    assert "duplicate manifest path" in result.stderr
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()


def test_converter_rejects_empty_manifest_frame(tmp_path):
    bedgraph = write_bedgraph(tmp_path, "S1", "HC04_A01", 32000)
    bedgraph.write_text("track type=bedGraph name=IS\n")
    manifest = write_manifest(tmp_path, "S1", "HC04_A01", [(32000, bedgraph.name)])

    result = run_converter(tmp_path, ["HC04_A01"], [manifest])

    assert result.returncode != 0
    assert "empty frame" in result.stderr
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()
