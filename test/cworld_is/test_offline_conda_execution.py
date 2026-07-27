from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
CWORLD_ENV = REPO / "env" / "cworld.yaml"


def dependency_lines(path):
    return {
        line.strip()[2:]
        for line in path.read_text().splitlines()
        if line.strip().startswith("- ")
    }


def test_cworld_environment_uses_only_the_required_arm_runtime():
    dependencies = dependency_lines(CWORLD_ENV)

    assert {
        "python>=3.10",
        "cooler>=0.10",
        "numpy",
        "pandas",
        "perl=5.32",
        "perl-gd>=2.57",
        "zlib",
    } <= dependencies
    assert not {
        dependency
        for dependency in dependencies
        if dependency.startswith("perl-file-slurp")
        or dependency.startswith("perl-statistics-basic")
        or dependency.startswith("perl-math-statistics")
        or dependency.startswith("perl-set-intervaltree")
        or dependency.startswith("perl-list-moreutils")
        or dependency.startswith("perl-list-util")
        or dependency.startswith("perl-data-dumper")
        or dependency.startswith("perl-file-copy-recursive")
        or dependency.startswith("perl-text-csv")
    }
