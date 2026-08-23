import importlib.util
from pathlib import Path
import sys
import textwrap

import pytest


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from tad_errors import TadError


def load_runner():
    path = REPO / "scripts" / "run_3dnetmod.py"
    spec = importlib.util.spec_from_file_location("run_3dnetmod", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def write_fake_stages(code_dir, fail_stage=None):
    for stage_name in (
        "3DNetMod_preprocess.py",
        "3DNetMod_GPS_MMCP.py",
        "3DNetMod_HSVM.py",
    ):
        stage = code_dir / stage_name
        stage.write_text(textwrap.dedent(f"""\
            import sys
            from pathlib import Path
            with (Path.cwd() / "stages.log").open("a") as handle:
                handle.write(Path(__file__).name + "\\n")
            if Path(__file__).name == {fail_stage!r}:
                raise SystemExit(7)
            if Path(__file__).name == "3DNetMod_HSVM.py":
                output = Path.cwd() / "FINAL_DOMAIN_CALLS"
                output.mkdir()
                (output / "calls_master_consistent_domains.txt").write_text(
                    "A02\\t4000\\t12000\\n"
                )
        """))


def test_stages_run_in_order_and_publish_canonical_calls(tmp_path):
    runner = load_runner()
    code_dir = tmp_path / "code"
    code_dir.mkdir()
    write_fake_stages(code_dir)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    settings = run_dir / "settings_S.raw.txt"
    settings.write_text("sample_1\tS\n")
    bed = tmp_path / "S.4kb.raw.3dnetmod.bed"
    bedpe = tmp_path / "S.4kb.raw.3dnetmod.bedpe"

    runner.run_stages(settings, code_dir, run_dir, sample="S", signal="raw", resolution=4000)
    runner.collect_final_calls(
        run_dir, bed, bedpe, {"A02": 20000},
        sample="S", signal="raw", resolution=4000,
    )

    assert (run_dir / "stages.log").read_text().splitlines() == [
        "3DNetMod_preprocess.py",
        "3DNetMod_GPS_MMCP.py",
        "3DNetMod_HSVM.py",
    ]
    assert bed.read_text() == "A02\t4000\t12000\n"
    assert bedpe.read_text() == "A02\t4000\t12000\tA02\t4000\t12000\n"


def test_stage_failure_names_real_context(tmp_path):
    runner = load_runner()
    code_dir = tmp_path / "code"
    code_dir.mkdir()
    write_fake_stages(code_dir, fail_stage="3DNetMod_GPS_MMCP.py")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    settings = run_dir / "settings_dynamic.balanced.txt"
    settings.write_text("sample_1\tdynamic\n")

    with pytest.raises(TadError) as error:
        runner.run_stages(
            settings, code_dir, run_dir,
            sample="dynamic", signal="balanced", resolution=8000,
        )
    message = str(error.value)
    assert message.startswith("ERROR[TAD_STAGE_FAILED]")
    assert "sample=dynamic" in message
    assert "signal=balanced" in message
    assert "resolution=8000" in message
    assert "stage=3DNetMod_GPS_MMCP.py" in message
    assert "returncode=7" in message


@pytest.mark.parametrize(
    "line,code",
    [
        ("A02\t-1\t4000\n", "TAD_INVALID_COORDINATE"),
        ("A02\t4000\t24000\n", "TAD_INVALID_COORDINATE"),
        ("unknown\t0\t4000\n", "TAD_INVALID_COORDINATE"),
    ],
)
def test_invalid_or_empty_calls_publish_nothing(tmp_path, line, code):
    runner = load_runner()
    final = tmp_path / "FINAL_DOMAIN_CALLS"
    final.mkdir()
    (final / "calls_master_consistent_domains.txt").write_text(line)
    bed = tmp_path / "calls.bed"
    bedpe = tmp_path / "calls.bedpe"
    with pytest.raises(TadError, match=rf"ERROR\[{code}\]"):
        runner.collect_final_calls(
            tmp_path, bed, bedpe, {"A02": 20000},
            sample="S", signal="raw", resolution=4000,
        )
    assert not bed.exists()
    assert not bedpe.exists()


def test_empty_calls_publish_empty_outputs(tmp_path):
    runner = load_runner()
    final = tmp_path / "FINAL_DOMAIN_CALLS"
    final.mkdir()
    (final / "calls_master_consistent_domains.txt").write_text("# no domains\n")
    bed = tmp_path / "calls.bed"
    bedpe = tmp_path / "calls.bedpe"

    runner.collect_final_calls(
        tmp_path, bed, bedpe, {"A02": 20000},
        sample="S", signal="raw", resolution=4000,
    )

    assert bed.read_text() == ""
    assert bedpe.read_text() == ""
