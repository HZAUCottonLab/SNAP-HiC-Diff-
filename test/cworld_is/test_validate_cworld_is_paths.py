import gzip
import importlib.util
import subprocess
import sys
import types
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "validate_cworld_is.py"


def load_validator(monkeypatch):
    fake_cooler = types.ModuleType("cooler")
    fake_cooler.Cooler = None
    monkeypatch.setitem(sys.modules, "cooler", fake_cooler)

    spec = importlib.util.spec_from_file_location("validate_cworld_is_under_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def option_value(command, option):
    return command[command.index(option) + 1]


def test_validator_uses_converter_cworld_layout_for_single_sample(monkeypatch, tmp_path):
    validator = load_validator(monkeypatch)
    outdir = tmp_path / "validation"
    sample = "S1"
    chrom = "HC04_A01"
    distance = 32000
    cworld_dir = tmp_path / "cworld"
    perl_script = cworld_dir / "scripts" / "perl" / "matrix2insulation.pl"
    perl_script.parent.mkdir(parents=True)
    perl_script.touch()

    expected_sample_dir = outdir / "cworld_insulation" / sample
    expected_prefix = expected_sample_dir / f"{sample}__{chrom}"
    cworld_command = None
    converter_command = None
    expected_parent_created_before_cworld = False
    real_run = subprocess.run

    def fake_run(command, **kwargs):
        nonlocal cworld_command
        nonlocal converter_command
        nonlocal expected_parent_created_before_cworld

        command = [str(part) for part in command]
        if command[1] == "scripts/convert_cool_to_cworld.py":
            matrix_path = Path(option_value(command, "--out"))
            matrix_path.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(matrix_path, "wt") as matrix:
                matrix.write(f"\tHiC_1|aaa|{chrom}|0\n")
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        if command[0] == "perl":
            cworld_command = command
            expected_parent_created_before_cworld = expected_sample_dir.is_dir()
            prefix = Path(option_value(command, "-o"))
            bedgraph_path = prefix.with_name(
                f"{prefix.name}--is64001--nt0--ids32001--ss1--immean.insulation.bedGraph"
            )
            bedgraph_path.write_text(
                "track type=bedGraph name=IS\n"
                f"{chrom}\t0\t16000\t0.25\n"
            )
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        if command[1] == "scripts/convert_cworld_IS_to_meltron_IS.py":
            converter_command = command
            return real_run(command, **kwargs)

        raise AssertionError(f"Unexpected command: {command}")

    class FakeBins:
        def __getitem__(self, _):
            return validator.pd.DataFrame({"chrom": [chrom]})

    class FakeCooler:
        def __init__(self, _):
            pass

        def bins(self):
            return FakeBins()

    monkeypatch.chdir(REPO)
    monkeypatch.setattr(validator.subprocess, "run", fake_run)
    monkeypatch.setattr(validator.cooler, "Cooler", FakeCooler)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--cool", str(tmp_path / "input.cool"),
            "--sample", sample,
            "--chrom", chrom,
            "--is-distance", str(distance),
            "--cworld-dir", str(cworld_dir),
            "--outdir", str(outdir),
        ],
    )

    run_error = None
    try:
        validator.run()
    except subprocess.CalledProcessError as error:
        run_error = error

    assert option_value(cworld_command, "-o") == str(expected_prefix)
    assert expected_parent_created_before_cworld
    assert option_value(converter_command, "--manifests") == str(
        expected_sample_dir / f"{sample}__{chrom}.is_manifest.tsv"
    )
    assert option_value(converter_command, "--outdir") == str(outdir)
    assert run_error is None
    assert (outdir / "IS" / sample / f"{sample}_IS.tsv").exists()
