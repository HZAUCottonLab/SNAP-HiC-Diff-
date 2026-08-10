from pathlib import Path
import math
import subprocess
import sys

import pandas as pd
import pytest


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "convert_cworld_IS_to_meltron_IS.py"


def write_bedgraph(root, sample, chrom, distance, rows):
    path = (
        root
        / "cworld_insulation"
        / sample
        / f"{sample}__{chrom}--is{distance + 1}--nt0--ids{distance // 2 + 1}--ss1--immean.insulation.bedGraph"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "".join(f"{chrom}\t{start}\t{end}\t{value}\n" for start, end, value in rows)
    path.write_text("track type=bedGraph name=IS\n" + body)
    return path


def write_manifest(root, sample, chrom, entries):
    path = root / "cworld_insulation" / sample / f"{sample}__{chrom}.is_manifest.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{distance}\t{filename}\n" for distance, filename in entries))
    return path


def run_converter(root, chroms, manifests, expected_bin_size, distances="32000"):
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--sample", "S1",
            "--chroms", *chroms,
            "--manifests", *(str(path) for path in manifests),
            "--outdir", str(root),
            "--is-distances", distances,
            "--expected-bin-size", str(expected_bin_size),
            "--log", str(root / "converter.log"),
        ],
        text=True,
        capture_output=True,
        check=False,
    )


@pytest.mark.parametrize("width", [4000, 10000, 7000])
def test_converter_exports_resolution_aligned_starts(tmp_path, width):
    rows = [(0, width, 0.25), (width, width * 2, 0.75)]
    bedgraph = write_bedgraph(tmp_path, "S1", "HC04_A01", 32000, rows)
    manifest = write_manifest(tmp_path, "S1", "HC04_A01", [(32000, bedgraph.name)])

    result = run_converter(tmp_path, ["HC04_A01"], [manifest], width)

    assert result.returncode == 0, result.stderr
    frame = pd.read_csv(tmp_path / "IS" / "S1" / "S1_IS.tsv", sep="\t")
    assert frame["viewpoint"].tolist() == [0, width]
    assert all(viewpoint % width == 0 for viewpoint in frame["viewpoint"])
    assert math.gcd(*frame["viewpoint"].tolist()) == width


def test_converter_rejects_width_that_disagrees_with_config(tmp_path):
    bedgraph = write_bedgraph(
        tmp_path,
        "S1",
        "HC04_A01",
        32000,
        [(0, 10000, 0.25), (10000, 20000, 0.75)],
    )
    manifest = write_manifest(tmp_path, "S1", "HC04_A01", [(32000, bedgraph.name)])

    result = run_converter(tmp_path, ["HC04_A01"], [manifest], 4000)

    assert result.returncode != 0
    assert "expected bin size 4000" in result.stderr
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()


def test_converter_rejects_off_grid_full_bin(tmp_path):
    bedgraph = write_bedgraph(
        tmp_path,
        "S1",
        "HC04_A01",
        32000,
        [(2000, 6000, 0.25), (6000, 10000, 0.75)],
    )
    manifest = write_manifest(tmp_path, "S1", "HC04_A01", [(32000, bedgraph.name)])

    result = run_converter(tmp_path, ["HC04_A01"], [manifest], 4000)

    assert result.returncode != 0
    assert "not aligned to bin size 4000" in result.stderr
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()


def test_converter_ignores_only_terminal_incomplete_bin(tmp_path):
    bedgraph = write_bedgraph(
        tmp_path,
        "S1",
        "HC04_A01",
        32000,
        [(0, 4000, 0.25), (4000, 8000, 0.75), (8000, 9000, "NA")],
    )
    manifest = write_manifest(tmp_path, "S1", "HC04_A01", [(32000, bedgraph.name)])

    result = run_converter(tmp_path, ["HC04_A01"], [manifest], 4000)

    assert result.returncode == 0, result.stderr
    frame = pd.read_csv(tmp_path / "IS" / "S1" / "S1_IS.tsv", sep="\t")
    assert frame["viewpoint"].tolist() == [0, 4000]


def test_converter_rejects_duplicate_output_key(tmp_path):
    rows = [(0, 4000, 0.25)]
    first_bedgraph = write_bedgraph(tmp_path, "S1", "HC04_A01", 32000, rows)
    second_bedgraph = write_bedgraph(tmp_path, "S1", "HC04_A02", 32000, rows)
    second_bedgraph.write_text(
        "track type=bedGraph name=IS\nHC04_A01\t0\t4000\t0.75\n"
    )
    first_manifest = write_manifest(
        tmp_path, "S1", "HC04_A01", [(32000, first_bedgraph.name)]
    )
    second_manifest = write_manifest(
        tmp_path, "S1", "HC04_A02", [(32000, second_bedgraph.name)]
    )

    result = run_converter(
        tmp_path,
        ["HC04_A01", "HC04_A02"],
        [first_manifest, second_manifest],
        4000,
    )

    assert result.returncode != 0
    assert "duplicate Meltron coordinates" in result.stderr
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()


def test_converter_preserves_manifest_distance_despite_cworld_suffix(tmp_path):
    bedgraph = write_bedgraph(
        tmp_path,
        "S1",
        "HC04_A01",
        64000,
        [(0, 4000, 0.25), (4000, 8000, 0.75)],
    )
    manifest = write_manifest(tmp_path, "S1", "HC04_A01", [(32000, bedgraph.name)])

    result = run_converter(tmp_path, ["HC04_A01"], [manifest], 4000)

    assert result.returncode == 0, result.stderr
    frame = pd.read_csv(tmp_path / "IS" / "S1" / "S1_IS.tsv", sep="\t")
    assert set(frame["IS_distance"]) == {32000}


def test_converter_rejects_duplicate_distance_in_manifest(tmp_path):
    bedgraph = write_bedgraph(
        tmp_path,
        "S1",
        "HC04_A01",
        32000,
        [(0, 4000, 0.25), (4000, 8000, 0.75)],
    )
    manifest = write_manifest(
        tmp_path,
        "S1",
        "HC04_A01",
        [(32000, bedgraph.name), (32000, bedgraph.name)],
    )

    result = run_converter(tmp_path, ["HC04_A01"], [manifest], 4000)

    assert result.returncode != 0
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()


def test_converter_rejects_missing_manifest_target(tmp_path):
    manifest = write_manifest(
        tmp_path,
        "S1",
        "HC04_A01",
        [(32000, "missing.bedGraph")],
    )

    result = run_converter(tmp_path, ["HC04_A01"], [manifest], 4000)

    assert result.returncode != 0
    assert "manifest target missing" in result.stderr
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()


def test_converter_rejects_unequal_chromosome_and_manifest_lists(tmp_path):
    bedgraph = write_bedgraph(
        tmp_path,
        "S1",
        "HC04_A01",
        32000,
        [(0, 4000, 0.25), (4000, 8000, 0.75)],
    )
    manifest = write_manifest(tmp_path, "S1", "HC04_A01", [(32000, bedgraph.name)])

    result = run_converter(tmp_path, ["HC04_A01", "HC04_A02"], [manifest], 4000)

    assert result.returncode != 0
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()


def test_converter_never_falls_back_to_legacy_bedgraph_glob(tmp_path):
    write_bedgraph(
        tmp_path,
        "S1",
        "HC04_A01",
        32000,
        [(0, 4000, 0.25), (4000, 8000, 0.75)],
    )

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
    bedgraph = write_bedgraph(
        tmp_path,
        "S1",
        "HC04_A01",
        32000,
        [(0, 4000, 0.25), (4000, 8000, 0.75)],
    )
    manifest = write_manifest(tmp_path, "S1", "HC04_A01", [(32000, bedgraph.name)])

    result = run_converter(
        tmp_path,
        ["HC04_A01", "HC04_A02"],
        [manifest, manifest],
        4000,
    )

    assert result.returncode != 0
    assert "duplicate manifest path" in result.stderr
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()


def test_converter_rejects_empty_manifest_frame(tmp_path):
    bedgraph = write_bedgraph(tmp_path, "S1", "HC04_A01", 32000, [])
    bedgraph.write_text("track type=bedGraph name=IS\n")
    manifest = write_manifest(tmp_path, "S1", "HC04_A01", [(32000, bedgraph.name)])

    result = run_converter(tmp_path, ["HC04_A01"], [manifest], 4000)

    assert result.returncode != 0
    assert "empty frame" in result.stderr
    assert not (tmp_path / "IS" / "S1" / "S1_IS.tsv").exists()
