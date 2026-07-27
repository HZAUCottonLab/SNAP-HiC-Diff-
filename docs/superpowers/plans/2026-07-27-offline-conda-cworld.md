# Offline Conda and CWorld Execution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Allow the complete SNAP-HiC-Diff workflow to run through dsub on offline ARM compute nodes by prebuilding every Conda environment on a login node and forbidding runtime package downloads.

**Architecture:** A single `scripts/conda_envs.sh` owns the three Snakemake environment specifications, creates them using Snakemake's hash-derived cache, validates runtime executables, and writes a checksum-bound ready marker. `run_snakemake.sh` validates that marker before it starts the controller; `config/dsub_submit.sh` injects the same Conda cache and offline setting into every rule-job wrapper.

**Tech Stack:** Bash, Snakemake 8 cluster-generic executor, Conda/libmamba, pytest, Conda `linux-aarch64` packages, dsub.

## Global Constraints

- Target platform is `linux-aarch64`; do not add x86-only packages.
- The CWorld runtime uses vendored `env/cworld-dekker`; do not reuse an environment path from `HiC_analysis_pipeline`.
- `perl-file-slurp` and the other unused non-core Perl dependencies must not appear in `env/cworld.yaml`.
- All three specs (`env/SNAP-HiC-Diff.yaml`, `env/cworld.yaml`, `env/meltron.yaml`) must have completed prefixes before dsub submission.
- Compute nodes must run with `CONDA_OFFLINE=true`; only a login-node bootstrap may contact Conda channels.
- Do not submit `dsub -s run_snakemake.sh` as part of implementation verification.

## File Map

| File | Responsibility |
| --- | --- |
| `env/cworld.yaml` | Minimal ARM-solvable Python/Perl runtime for conversion and CWorld insulation. |
| `scripts/conda_envs.sh` | Create, validate, mark, and verify the three hash-derived Snakemake Conda prefixes. |
| `run_snakemake.sh` | Reject an unprepared/stale cache before running the dsub controller in offline mode. |
| `config/dsub_submit.sh` | Inject cache and offline exports into every generated rule-job script. |
| `README.md` | Login-node bootstrap and dsub submission instructions. |
| `test/cworld_is/test_offline_conda_execution.py` | Regression coverage for ARM-safe CWorld dependencies and offline script contracts. |

---

### Task 1: Make the CWorld Conda specification solvable on ARM

**Files:**
- Modify: `env/cworld.yaml`
- Create: `test/cworld_is/test_offline_conda_execution.py`

**Interfaces:**
- Consumes: `env/cworld-dekker/lib/cworld/dekker.pm` and `scripts/convert_cool_to_cworld.py`.
- Produces: an environment specification that `conda env create --dry-run --file env/cworld.yaml` can solve on `linux-aarch64`.

- [ ] **Step 1: Write the failing dependency regression test**

```python
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
```

- [ ] **Step 2: Run the test and verify the current specification fails**

Run:

```bash
pytest -q test/cworld_is/test_offline_conda_execution.py::test_cworld_environment_uses_only_the_required_arm_runtime
```

Expected: failure because `perl-file-slurp>=9999.34` is still present and `perl=5.32` is absent.

- [ ] **Step 3: Replace the CWorld dependency list with the minimal runtime**

Replace `env/cworld.yaml` dependencies with:

```yaml
dependencies:
  - python>=3.10
  - cooler>=0.10
  - numpy
  - pandas
  - perl=5.32
  - perl-gd>=2.57
  - zlib
```

Keep the existing channels and environment name. Do not add the broader package set used by the reference project: this workflow neither builds CWorld from source nor runs its R/image utilities.

- [ ] **Step 4: Verify the test and the ARM solver**

Run:

```bash
pytest -q test/cworld_is/test_offline_conda_execution.py::test_cworld_environment_uses_only_the_required_arm_runtime
conda env create --dry-run --file env/cworld.yaml
```

Expected: the pytest target passes and Conda reports a solvable `linux-aarch64` transaction without `PackagesNotFoundError`.

- [ ] **Step 5: Commit the isolated dependency repair**

```bash
git add env/cworld.yaml test/cworld_is/test_offline_conda_execution.py
git commit -m "fix(cworld): use arm-solvable runtime dependencies"
```

### Task 2: Add a reusable login-node environment bootstrap and ready marker

**Files:**
- Create: `scripts/conda_envs.sh`
- Modify: `test/cworld_is/test_offline_conda_execution.py`
- Modify: `README.md`

**Interfaces:**
- Consumes: an action argument, `bootstrap` or `verify`; optional `SNAP_HIC_DIFF_ROOT`, `SNAKEMAKE_BIN`, and `CONDA_PREFIX_DIR` overrides for testing.
- Produces: `.snakemake/conda/offline.ready.tsv`, with one tab-separated `spec path checksum` record for each of the three source specifications.

- [ ] **Step 1: Extend the test file with the expected bootstrap contract**

```python
CONDA_HELPER = REPO / "scripts" / "conda_envs.sh"


def test_conda_helper_declares_all_specs_marker_and_runtime_checks():
    script = CONDA_HELPER.read_text()

    for spec in (
        "env/SNAP-HiC-Diff.yaml",
        "env/cworld.yaml",
        "env/meltron.yaml",
    ):
        assert spec in script
    assert 'ACTION="${1:-}"' in script
    assert 'offline.ready.tsv' in script
    assert '--conda-create-envs-only' in script
    assert 'conda-meta/history' in script
    assert 'perl -I "$ROOT/env/cworld-dekker/lib" -Mcworld::dekker' in script
    assert 'Rscript' in script
```

- [ ] **Step 2: Run the test and verify it fails because the helper is absent**

Run:

```bash
pytest -q test/cworld_is/test_offline_conda_execution.py::test_conda_helper_declares_all_specs_marker_and_runtime_checks
```

Expected: `FileNotFoundError` for `scripts/conda_envs.sh`.

- [ ] **Step 3: Implement `scripts/conda_envs.sh`**

Implement Bash functions with these exact responsibilities:

```bash
SPECS=(env/SNAP-HiC-Diff.yaml env/cworld.yaml env/meltron.yaml)
ROOT="${SNAP_HIC_DIFF_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
CONDA_PREFIX_DIR="${CONDA_PREFIX_DIR:-$ROOT/.snakemake/conda}"
READY_MARKER="$CONDA_PREFIX_DIR/offline.ready.tsv"
SNAKEMAKE_BIN="${SNAKEMAKE_BIN:-snakemake}"

find_completed_prefix() {
    local spec="$1" cached prefix
    for cached in "$CONDA_PREFIX_DIR"/*.yaml; do
        [[ -f "$cached" ]] || continue
        cmp -s "$ROOT/$spec" "$cached" || continue
        prefix="${cached%.yaml}"
        [[ -f "$prefix/conda-meta/history" ]] || continue
        printf '%s\n' "$prefix"
        return 0
    done
    return 1
}
```

`bootstrap` must remove any existing ready marker, invoke:

```bash
"$SNAKEMAKE_BIN" --use-conda --conda-create-envs-only \
    --conda-prefix "$CONDA_PREFIX_DIR" \
    -s "$ROOT/Snakefile" --configfile "$ROOT/config/config.yaml"
```

Then locate one complete prefix for every spec, run these validations, and atomically rename a temporary marker into place:

```bash
"$snap_prefix/bin/python" -c 'import cooler, cooltools, numpy, pandas'
"$snap_prefix/bin/hictk" --version
"$cworld_prefix/bin/python" -c 'import cooler, numpy, pandas'
"$cworld_prefix/bin/perl" -I "$ROOT/env/cworld-dekker/lib" -Mcworld::dekker -e 'print $cworld::dekker::VERSION, "\n"'
"$meltron_prefix/bin/Rscript" -e 'library(argparse); library(data.table); library(DescTools)'
```

Each marker record must contain the relative spec path, the SHA-256 of its source YAML, and its absolute prefix path. `verify` must parse records without `source`, require exactly the same three specs, recompute each SHA-256, verify `conda-meta/history`, and rerun the same runtime command for the matching environment. Any mismatch must exit nonzero and say:

```text
Run: bash scripts/conda_envs.sh bootstrap
```

- [ ] **Step 4: Add operator documentation in `README.md`**

Replace the stale local/HPC run snippets with:

```bash
# On an ARM login node with Conda channel access; run after any env/*.yaml change.
bash scripts/conda_envs.sh bootstrap

# Validate without contacting channels.
CONDA_OFFLINE=true bash scripts/conda_envs.sh verify

# Submit the already-prepared workflow controller.
dsub -s run_snakemake.sh
```

State explicitly that compute nodes are offline and must never be used to create Conda environments.

- [ ] **Step 5: Verify the helper contract and marker behavior**

Run:

```bash
pytest -q test/cworld_is/test_offline_conda_execution.py::test_conda_helper_declares_all_specs_marker_and_runtime_checks
bash scripts/conda_envs.sh bootstrap
CONDA_OFFLINE=true bash scripts/conda_envs.sh verify
```

Expected: the test passes; bootstrap creates three completed prefixes plus `.snakemake/conda/offline.ready.tsv`; offline verification completes without channel access.

- [ ] **Step 6: Commit the bootstrap boundary**

```bash
git add scripts/conda_envs.sh README.md test/cworld_is/test_offline_conda_execution.py
git commit -m "feat(conda): bootstrap offline workflow environments"
```

### Task 3: Enforce cache validation and offline Conda in dsub execution

**Files:**
- Modify: `run_snakemake.sh`
- Modify: `config/dsub_submit.sh`
- Modify: `test/cworld_is/test_offline_conda_execution.py`

**Interfaces:**
- Consumes: a valid `.snakemake/conda/offline.ready.tsv` made by `scripts/conda_envs.sh bootstrap`.
- Produces: controller and generated rule-job scripts with `CONDA_ENVS_PATH` set to the project cache and `CONDA_OFFLINE=true` before Snakemake/Conda can run.

- [ ] **Step 1: Add failing tests for controller and worker offline guards**

```python
RUNNER = REPO / "run_snakemake.sh"
SUBMITTER = REPO / "config" / "dsub_submit.sh"


def test_controller_validates_the_ready_marker_before_snakemake():
    script = RUNNER.read_text()

    assert 'bash "$WORKDIR/scripts/conda_envs.sh" verify' in script
    assert 'export CONDA_ENVS_PATH="$WORKDIR/.snakemake/conda"' in script
    assert 'export CONDA_OFFLINE=true' in script
    assert '--conda-prefix "$CONDA_ENVS_PATH"' in script
