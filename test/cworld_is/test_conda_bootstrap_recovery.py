import os
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
CONDA_HELPER = REPO / "scripts" / "conda_envs.sh"


def make_executable(path: Path, content: str) -> None:
    path.write_text(content)
    path.chmod(0o755)


def test_bootstrap_creates_snakemake_compatible_prefixes_with_mamba(tmp_path: Path) -> None:
    root = tmp_path / "project"
    env_dir = root / "env"
    env_dir.mkdir(parents=True)
    for spec in ("SNAP-HiC-Diff.yaml", "cworld.yaml", "meltron.yaml", "tad.yaml"):
        (env_dir / spec).write_text(f"name: {spec}\n")

    prefix_dir = tmp_path / "conda-prefix"
    mamba_log = tmp_path / "mamba.log"
    fake_snakemake = tmp_path / "fake_snakemake.py"
    fake_mamba = tmp_path / "fake_mamba.py"

    make_executable(
        fake_snakemake,
        """#!/usr/bin/env python3
import os
import shutil
import sys
from pathlib import Path

root = Path(os.environ["SNAP_HIC_DIFF_ROOT"])
prefix_dir = Path(os.environ["CONDA_PREFIX_DIR"])
state = prefix_dir / "snakemake.calls"

def make_prefix(prefix, tools):
    (prefix / "conda-meta").mkdir(parents=True, exist_ok=True)
    (prefix / "conda-meta" / "history").write_text("complete\\n")
    bin_dir = prefix / "bin"
    bin_dir.mkdir(exist_ok=True)
    for tool in tools:
        executable = bin_dir / tool
        executable.write_text("#!/usr/bin/env bash\\nexit 0\\n")
        executable.chmod(0o755)

prefix_dir.mkdir(parents=True, exist_ok=True)
if not state.exists():
    shutil.copyfile(root / "env" / "meltron.yaml", prefix_dir / "meltron-cache.yaml")
    (prefix_dir / "meltron-cache").mkdir()
    state.touch()
    sys.exit(137)

for spec, tools in {
    "SNAP-HiC-Diff.yaml": ("python", "hictk"),
    "cworld.yaml": ("python", "perl", "bedtools"),
}.items():
    cache = prefix_dir / f"{spec}.cache.yaml"
    shutil.copyfile(root / "env" / spec, cache)
    make_prefix(Path(str(cache)[:-5]), tools)
sys.exit(0)
""",
    )
    make_executable(
        fake_mamba,
        """#!/usr/bin/env python3
import os
import sys
from pathlib import Path

args = sys.argv[1:]
prefix = Path(args[args.index("--prefix") + 1])
(prefix / "conda-meta").mkdir(parents=True, exist_ok=True)
(prefix / "conda-meta" / "history").write_text("complete\\n")
bin_dir = prefix / "bin"
bin_dir.mkdir(exist_ok=True)
rscript = bin_dir / "Rscript"
rscript.write_text("#!/usr/bin/env bash\\nexit 0\\n")
rscript.chmod(0o755)
spec = Path(args[args.index("--file") + 1]).name
for tool in {
    "SNAP-HiC-Diff.yaml": ("python", "hictk"),
    "cworld.yaml": ("python", "perl", "bedtools"),
    "meltron.yaml": ("Rscript", "python"),
    "tad.yaml": ("python",),
}[spec]:
    executable = bin_dir / tool
    executable.write_text("#!/usr/bin/env bash")
    executable.chmod(0o755)
Path(os.environ["MAMBA_LOG"]).write_text("\\n".join(args))
if spec == "meltron.yaml":
    for package in ("argparser", "data.table", "DescTools"):
        description = prefix / "lib" / "R" / "library" / package / "DESCRIPTION"
        description.parent.mkdir(parents=True, exist_ok=True)
        description.write_text(f"Package: {package}\\n")

""",
    )

    environment = {
        **os.environ,
        "SNAP_HIC_DIFF_ROOT": str(root),
        "CONDA_PREFIX_DIR": str(prefix_dir),
        "SNAKEMAKE_BIN": str(fake_snakemake),
        "MAMBA_BIN": str(fake_mamba),
        "MAMBA_LOG": str(mamba_log),
    }
    bootstrap = subprocess.run(
        ["bash", str(CONDA_HELPER), "bootstrap"],
        text=True,
        capture_output=True,
        env=environment,
        check=False,
    )

    assert bootstrap.returncode == 0, bootstrap.stderr
    assert not (prefix_dir / "snakemake.calls").exists()
    assert mamba_log.read_text().splitlines()[:3] == ["env", "create", "--yes"]
    marker = prefix_dir / "offline.ready.tsv"
    assert len(marker.read_text().splitlines()) == 4

    verify = subprocess.run(
        ["bash", str(CONDA_HELPER), "verify"],
        text=True,
        capture_output=True,
        env=environment,
        check=False,
    )
    assert verify.returncode == 0, verify.stderr
