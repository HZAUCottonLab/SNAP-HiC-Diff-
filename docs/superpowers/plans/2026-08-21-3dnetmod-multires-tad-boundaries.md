# Multi-resolution 3DNetMod TAD and Boundary Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add portable 4/8/20 kb 3DNetMod TAD calling for raw and balanced contacts, with per-resolution and multi-resolution merged boundary BED files for every manifest sample.

**Architecture:** Vendor the validated Python 3 3DNetMod closure from `HiC_analysis_pipeline`, then schedule independent `sample x resolution x signal` preparation and calling jobs. A focused boundary script expands TAD endpoints and applies strict BED-overlap union first within each resolution and then across resolutions; a dedicated Conda environment is built on the login node and consumed offline by workers.

**Tech Stack:** Snakemake 8, Python 3.11, Cooler/HDF5, NumPy below 2, SciPy, scikit-learn, Pandas, pytest, Conda/Mamba, Donau dsub.

## Global Constraints

- Implement the approved design in `docs/superpowers/specs/2026-08-21-3dnetmod-multires-tad-boundaries-design.md`.
- TAD resolutions are exactly 4000, 8000, and 20000 bp; signals are exactly `raw` and `balanced`.
- Raw values are finite positive pixel counts; balanced values are `count * weight(bin1) * weight(bin2)` and must be finite and positive.
- Raw and balanced runs use identical 3DNetMod D parameters, including the 4 kb override.
- Analyze every chromosome present in the mcool; default `excluded_chromosomes` is `[]`; default `badregionfilter` is `false`.
- A boundary at coordinate `p` and resolution `R` is `[max(0, p-R), min(chrom_length, p+R))`.
- Merge intervals only when `next_start < current_end`; endpoint adjacency is not overlap; merging is transitive and uses interval union.
- Merge raw and balanced boundaries separately and retain `support_resolutions`, `n_supporting_resolutions`, and `n_source_boundaries`.
- Stable error messages start with `ERROR[CODE]` and dynamically include real sample, signal, resolution, stage, return code, and path when applicable.
- Required codes are `TAD_CONFIG_INVALID`, `TAD_MISSING_RESOLUTION`, `TAD_MISSING_WEIGHT`, `TAD_EMPTY_INTERACTIONS`, `TAD_STAGE_FAILED`, `TAD_EMPTY_CALLS`, `TAD_INVALID_COORDINATE`, and `TAD_BOUNDARY_MERGE_FAILED`.
- Publish stable outputs atomically per file; retain run-local intermediates and logs after failure.
- Add every TAD, per-resolution boundary, and multi-resolution boundary result to `rule all`.
- Use 20 threads for each 3DNetMod call; boundary jobs use one thread.
- `env/tad.yaml` is created only on the login node; compute nodes stay offline and only consume the prepared prefix.
- The finished repository must not require `/share/home/xhhuang24/project/HiC_analysis_pipeline` at runtime.
- Preserve all unrelated worktree changes and stage only files named by the current task.

## File Map

**Create:**

- `env/tad.yaml` — isolated 3DNetMod runtime.
- `rules/tad.smk` — preparation, calling, boundary extraction, and boundary merge rules.
- `scripts/python3_3dnetmod/` — vendored algorithm closure and seed file.
- `scripts/tad_errors.py` — shared stable error representation.
- `scripts/tad_from_mcool.py` — signal-aware mcool exporter and settings renderer.
- `scripts/run_3dnetmod.py` — ordered stage runner and validated final collector.
- `scripts/tad_boundaries.py` — boundary expansion, overlap union, provenance, and CLI.
- `scripts/tad_config.py` — focused configuration validation.
- `test/tad/test_mcool_conversion.py` — 8 kb and `--balance` conversion contract.
- `test/tad/test_vendored_3dnetmod.py` — vendored closure contract.
- `test/tad/test_tad_from_mcool.py` — raw/balanced input and stable-error tests.
- `test/tad/test_tad_runner.py` — stage order, failure context, coordinate, and publication tests.
- `test/tad/test_tad_boundaries.py` — expansion and merge semantics.
- `test/tad/test_tad_config.py` — TAD configuration parser tests.
- `test/tad/test_tad_workflow.py` — rules, outputs, and DAG contract.
- `CONTEXT.md` — project vocabulary for TADs, signals, and boundary regions.
- `docs/adr/0001-3dnetmod-raw-and-balanced-signals.md` — scientific input decision.

**Modify:**

- `config/config.yaml:17-19,181-187` — add 8 kb conversion and the approved TAD section.
- `scripts/convert_cool_to_mcool.py:43-55` — balance newly generated mcools.
- `scripts/conda_envs.sh:12,74-117,161-183` — bootstrap and verify `env/tad.yaml`.
- `test/cworld_is/test_offline_conda_execution.py:4-93` — require the TAD environment and validation case.
- `test/cworld_is/test_conda_bootstrap_recovery.py:15-121` — include a fourth fake prefix and marker record.
- `Snakefile:4-27,64-72,123-140,146-203` — load TAD settings, define paths, include rules, and extend `rule all`.
- `README.md:15-55,85-135,207-227` — document login-node bootstrap, TAD outputs, and targeted execution.

---

### Task 1: Balance future mcool conversions and include 8 kb

**Files:**
- Create: `test/tad/test_mcool_conversion.py`
- Modify: `scripts/convert_cool_to_mcool.py:43-55`
- Modify: `config/config.yaml:17-19`

**Interfaces:**
- Consumes: existing `run_convert(cool_uri, mcool_out, resolutions, threads, log_file)`.
- Produces: new mcools with balanced 4/8/20/100 kb layers; later TAD rules rely on the 8 kb group and `bins/weight`.

- [ ] **Step 1: Write the failing conversion contract test**

```python
# test/tad/test_mcool_conversion.py
import importlib.util
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import yaml


REPO = Path(__file__).resolve().parents[2]


def load_converter():
    path = REPO / "scripts" / "convert_cool_to_mcool.py"
    spec = importlib.util.spec_from_file_location("convert_cool_to_mcool", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_conversion_builds_balanced_4_8_20_100kb_mcool():
    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    assert config["convert"]["resolutions"] == [4000, 8000, 20000, 100000]

    converter = load_converter()
    with patch.object(converter.subprocess, "run") as run:
        converter.run_convert(
            "input.cool",
            "output.mcool",
            config["convert"]["resolutions"],
            20,
            StringIO(),
        )

    run.assert_called_once_with(
        [
            "cooler", "zoomify", "input.cool",
            "-o", "output.mcool",
            "-r", "4000,8000,20000,100000",
            "--balance",
            "--nproc", "20",
        ],
        check=True,
    )
```

- [ ] **Step 2: Run the test and verify the missing behavior**

Run:

```bash
pytest test/tad/test_mcool_conversion.py -v
```

Expected: FAIL because `convert.resolutions` lacks 8000 and the command lacks `--balance`.

- [ ] **Step 3: Add 8 kb and the native Cooler balance flag**

Change the config block to:

```yaml
convert:
  resolutions: [4000, 8000, 20000, 100000]
  threads: 20
```

Change the command construction to:

```python
cmd = [
    "cooler", "zoomify",
    cool_uri,
    "-o", mcool_out,
    "-r", ",".join(map(str, resolutions)),
    "--balance",
    "--nproc", str(threads),
]
```

- [ ] **Step 4: Run the focused test**

Run:

```bash
pytest test/tad/test_mcool_conversion.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit the conversion contract**

```bash
git add config/config.yaml scripts/convert_cool_to_mcool.py test/tad/test_mcool_conversion.py
git commit -m "feat: balance 3DNetMod mcool resolutions"
```

### Task 2: Vendor the validated Python 3 3DNetMod closure

**Files:**
- Create: `scripts/python3_3dnetmod/`
- Create: `test/tad/test_vendored_3dnetmod.py`

**Interfaces:**
- Consumes: tracked source directory at HiC analysis commit `38c7c8d`.
- Produces: local entrypoints `3DNetMod_preprocess.py`, `3DNetMod_GPS_MMCP.py`, and `3DNetMod_HSVM.py`, plus every recursive import and seed resource used by Task 5.

- [ ] **Step 1: Write the failing vendored-closure test**

```python
# test/tad/test_vendored_3dnetmod.py
from pathlib import Path
import runpy


REPO = Path(__file__).resolve().parents[2]
CODE_DIR = REPO / "scripts" / "python3_3dnetmod"


def test_vendored_entrypoints_recursive_closure_and_seed_exist():
    required = {
        "3DNetMod_preprocess.py",
        "3DNetMod_GPS_MMCP.py",
        "3DNetMod_HSVM.py",
        "outdir_helper.py",
        "read_settings.py",
        "read_seeds.py",
        "gen_louvain_100_seeds.txt",
        "find_bad_regions_mean_var_simplified_per_chrom.py",
        "load_bad_regions_simplified_per_chrom.py",
        "find_max_gamma_simplified_per_chrom.py",
        "ThreeDNetMod_GPS_MMCP_per_chr_simplified_per_chrom.py",
    }
    assert required <= {path.name for path in CODE_DIR.iterdir()}


def test_vendored_code_has_no_runtime_reference_to_source_checkout():
    source = "\n".join(
        path.read_text(errors="replace")
        for path in CODE_DIR.glob("*.py")
    )
    assert "/share/home/xhhuang24/project/HiC_analysis_pipeline" not in source
    assert "PanTE/compare_mtx" not in source


def test_consensus_filename_is_deterministic_and_filesystem_safe():
    helper = runpy.run_path(str(CODE_DIR / "outdir_helper.py"))
    name = helper["consensus_filtered_filename"](
        "sample", "long_parameter_tag_" * 20, 0.4050000000000002
    )
    assert name == helper["consensus_filtered_filename"](
        "sample", "long_parameter_tag_" * 20, 0.4050000000000002
    )
    assert len(name.encode()) <= 255
```

- [ ] **Step 2: Run the test and verify the directory is absent**

Run:

```bash
pytest test/tad/test_vendored_3dnetmod.py -v
```

Expected: FAIL with `FileNotFoundError` for `scripts/python3_3dnetmod`.

- [ ] **Step 3: Verify the approved source commit and closure entrypoint exist**

Run:

```bash
git -C /share/home/xhhuang24/project/HiC_analysis_pipeline cat-file -e '38c7c8d^{commit}'
git -C /share/home/xhhuang24/project/HiC_analysis_pipeline ls-tree -r --name-only 38c7c8d scripts/python3_3dnetmod | grep -Fx scripts/python3_3dnetmod/3DNetMod_preprocess.py
```

Expected: both commands exit 0. Stop if the commit or entrypoint is absent.

- [ ] **Step 4: Export only the tracked closure from that commit**

Run:

```bash
test ! -e scripts/python3_3dnetmod
git -C /share/home/xhhuang24/project/HiC_analysis_pipeline archive 38c7c8d scripts/python3_3dnetmod | tar -x
```

Expected: `scripts/python3_3dnetmod/3DNetMod_preprocess.py` exists and the directory is about 988 KB. This excludes untracked files from the source checkout and does not edit algorithm files.

- [ ] **Step 5: Run the closure tests**

Run:

```bash
pytest test/tad/test_vendored_3dnetmod.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit the vendored algorithm as one mechanical change**

```bash
git add scripts/python3_3dnetmod test/tad/test_vendored_3dnetmod.py
git commit -m "feat: vendor Python 3 3DNetMod"
```

### Task 3: Add the login-node TAD Conda environment

**Files:**
- Create: `env/tad.yaml`
- Modify: `scripts/conda_envs.sh:12,74-117,161-183`
- Modify: `test/cworld_is/test_offline_conda_execution.py:4-93`
- Modify: `test/cworld_is/test_conda_bootstrap_recovery.py:15-121`

**Interfaces:**
- Consumes: existing `scripts/conda_envs.sh bootstrap|verify` cache protocol.
- Produces: one completed `env/tad.yaml` prefix and one fourth record in `offline.ready.tsv`; all TAD rules use `../env/tad.yaml`.

- [ ] **Step 1: Extend the environment contract tests first**

Add to `test/cworld_is/test_offline_conda_execution.py`:

```python
TAD_ENV = REPO / "env" / "tad.yaml"


def test_tad_environment_is_isolated_and_pinned():
    dependencies = dependency_lines(TAD_ENV)
    assert {
        "python=3.11",
        "cooler",
        "h5py",
        "numpy<2",
        "scipy",
        "matplotlib=3.8.*",
        "scikit-learn",
        "pandas",
        "pytest",
    } <= dependencies
```

Extend the spec loop in `test_conda_helper_declares_all_specs_marker_and_runtime_checks`:

```python
for spec in (
    "env/SNAP-HiC-Diff.yaml",
    "env/cworld.yaml",
    "env/meltron.yaml",
    "env/tad.yaml",
):
    assert spec in script
```

In `test_conda_bootstrap_recovery.py`, create all four fake specs:

```python
for spec in ("SNAP-HiC-Diff.yaml", "cworld.yaml", "meltron.yaml", "tad.yaml"):
    (env_dir / spec).write_text(f"name: {spec}\n")
```

Add the fake TAD prefix tool map entry:

```python
"tad.yaml": ("python",),
```

Change the marker assertion to:

```python
assert len(marker.read_text().splitlines()) == 4
```

- [ ] **Step 2: Run the environment tests and verify they fail**

Run:

```bash
pytest test/cworld_is/test_offline_conda_execution.py test/cworld_is/test_conda_bootstrap_recovery.py -v
```

Expected: FAIL because `env/tad.yaml` and its helper case do not exist.

- [ ] **Step 3: Add the dedicated environment specification**

```yaml
# env/tad.yaml
name: snap_hic_diff_tad
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - cooler
  - h5py
  - numpy<2
  - scipy
  - matplotlib=3.8.*
  - scikit-learn
  - pandas
  - pytest
```

- [ ] **Step 4: Register and validate the TAD prefix**

Change the specs array to:

```bash
SPECS=(env/SNAP-HiC-Diff.yaml env/cworld.yaml env/meltron.yaml env/tad.yaml)
```

Add this `validate_prefix` case before the catch-all case:

```bash
        env/tad.yaml)
            local tad_prefix="$prefix"
            [[ -x "$tad_prefix/bin/python" ]] || die "Missing Python in $tad_prefix"
            validate_python_package "$tad_prefix/bin/python" cooler
            validate_python_package "$tad_prefix/bin/python" h5py
            validate_python_package "$tad_prefix/bin/python" numpy
            validate_python_package "$tad_prefix/bin/python" scipy
            validate_python_package "$tad_prefix/bin/python" matplotlib
            validate_python_package "$tad_prefix/bin/python" scikit-learn
            validate_python_package "$tad_prefix/bin/python" pandas
            ;;
```

Extend the ready-marker whitelist exactly:

```bash
            env/SNAP-HiC-Diff.yaml|env/cworld.yaml|env/meltron.yaml|env/tad.yaml) ;;
```

- [ ] **Step 5: Run the static and fake-bootstrap tests**

Run:

```bash
bash -n scripts/conda_envs.sh
pytest test/cworld_is/test_offline_conda_execution.py test/cworld_is/test_conda_bootstrap_recovery.py -v
```

Expected: shell syntax succeeds and all tests PASS. Do not create the real environment in this task; Task 8 performs login-node bootstrap after all code is ready.

- [ ] **Step 6: Commit the isolated runtime**

```bash
git add env/tad.yaml scripts/conda_envs.sh test/cworld_is/test_offline_conda_execution.py test/cworld_is/test_conda_bootstrap_recovery.py
git commit -m "feat: add offline 3DNetMod environment"
```

### Task 4: Build signal-aware 3DNetMod inputs

**Files:**
- Create: `scripts/tad_errors.py`
- Create: `scripts/tad_from_mcool.py`
- Create: `test/tad/test_tad_from_mcool.py`

**Interfaces:**
- Consumes: `mcool_path`, integer resolution, sample, signal, D config JSON, processors, bad-region settings, and chromosome exclusions.
- Produces: `TadError(code, **context)`, global one-based BED bins, raw or balanced interaction triples, one signal-specific settings file, and one atomic ready manifest.

- [ ] **Step 1: Write focused input and error tests**

```python
# test/tad/test_tad_from_mcool.py
import importlib.util
from pathlib import Path
import sys

import cooler
import numpy as np
import pandas as pd
import pytest


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from tad_errors import TadError


def load_adapter():
    path = REPO / "scripts" / "tad_from_mcool.py"
    spec = importlib.util.spec_from_file_location("tad_from_mcool", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def tiny_mcool(tmp_path):
    path = tmp_path / "tiny.mcool"
    bins = pd.DataFrame({
        "chrom": ["A02", "A02", "A03"],
        "start": [0, 4000, 0],
        "end": [4000, 8000, 4000],
        "weight": [0.25, 0.5, np.nan],
    })
    pixels = pd.DataFrame({
        "bin1_id": [0, 0, 0, 1],
        "bin2_id": [0, 1, 2, 1],
        "count": [8.0, 4.0, 4.0, 0.0],
    })
    cooler.create_cooler(
        f"{path}::/resolutions/4000", bins, pixels,
        dtypes={"count": "float64"},
    )
    return path


def test_raw_values_and_balanced_values_with_nonfinite_weight_filter(tiny_mcool):
    adapter = load_adapter()
    clr = cooler.Cooler(f"{tiny_mcool}::/resolutions/4000")
    assert list(adapter.iter_pixels_as_triples(clr, "raw")) == [
        (1, 1, 8.0),
        (1, 2, 4.0),
        (1, 3, 4.0),
    ]
    assert list(adapter.iter_pixels_as_triples(clr, "balanced")) == [
        (1, 1, 0.5),
        (1, 2, 0.5),
    ]


def test_4kb_override_is_identical_for_both_signal_runs():
    adapter = load_adapter()
    d_config = {
        "size_threshold_bins": 4,
        "variance_thresholds": [45, 45, 60, 45, 45],
        "resolution_overrides": {
            "4000": {
                "size_threshold_bins": 10,
                "variance_thresholds": [20, 20, 30, 20, 20],
            }
        },
    }
    raw_effective = adapter.effective_d_config(d_config, 4000)
    balanced_effective = adapter.effective_d_config(d_config, 4000)
    assert raw_effective == balanced_effective
    assert raw_effective["size_threshold_bins"] == 10
    assert raw_effective["variance_thresholds"] == [20, 20, 30, 20, 20]


def test_balanced_mode_reports_dynamic_missing_weight_context(tmp_path):
    adapter = load_adapter()
    path = tmp_path / "unbalanced.mcool"
    bins = pd.DataFrame({
        "chrom": ["A02"], "start": [0], "end": [4000]
    })
    pixels = pd.DataFrame({"bin1_id": [0], "bin2_id": [0], "count": [1]})
    cooler.create_cooler(f"{path}::/resolutions/4000", bins, pixels)

    with pytest.raises(TadError) as error:
        adapter.validate_mcool_resolution(path, 4000, "sample-X", "balanced")
    assert str(error.value).startswith("ERROR[TAD_MISSING_WEIGHT]")
    assert "sample=sample-X" in str(error.value)
    assert "signal=balanced" in str(error.value)
    assert "resolution=4000" in str(error.value)


def test_missing_resolution_and_empty_contacts_have_stable_codes(tiny_mcool, tmp_path):
    adapter = load_adapter()
    with pytest.raises(TadError, match=r"ERROR\[TAD_MISSING_RESOLUTION\]"):
        adapter.validate_mcool_resolution(tiny_mcool, 8000, "S", "raw")

    clr = cooler.Cooler(f"{tiny_mcool}::/resolutions/4000")
    clr.pixels()[:]
    with pytest.raises(TadError, match=r"ERROR\[TAD_EMPTY_INTERACTIONS\]"):
        adapter.write_interactions(
            tmp_path / "empty.interactions",
            [],
            sample="S",
            signal="raw",
            resolution=4000,
        )


def test_global_bed_ids_are_one_based(tiny_mcool, tmp_path):
    adapter = load_adapter()
    clr = cooler.Cooler(f"{tiny_mcool}::/resolutions/4000")
    output = tmp_path / "S.bed"
    assert adapter.write_bed_from_bins(output, clr) == 3
    assert output.read_text() == (
        "A02\t0\t4000\t1\n"
        "A02\t4000\t8000\t2\n"
        "A03\t0\t4000\t3\n"
    )
```

- [ ] **Step 2: Run the tests and verify imports fail**

Run:

```bash
pytest test/tad/test_tad_from_mcool.py -v
```

Expected: collection FAIL because `scripts/tad_errors.py` and `scripts/tad_from_mcool.py` do not exist.

- [ ] **Step 3: Add the shared stable error type**

```python
# scripts/tad_errors.py
class TadError(RuntimeError):
    def __init__(self, code, **context):
        self.code = code
        self.context = context
        details = " ".join(
            f"{key}={value}" for key, value in context.items() if value is not None
        )
        super().__init__(f"ERROR[{code}] {details}".rstrip())
```

- [ ] **Step 4: Copy the validated adapter and make signal handling explicit**

Use `/share/home/xhhuang24/project/HiC_analysis_pipeline/scripts/tad_from_mcool.py` as the mechanical base, then add `sys` and `TadError` imports and replace the input functions with:

```python
import sys
from tad_errors import TadError

VALID_SIGNALS = {"raw", "balanced"}


def validate_mcool_resolution(mcool_path, resolution, sample, signal):
    bins_path = f"resolutions/{resolution}/bins"
    weight_path = f"{bins_path}/weight"
    with h5py.File(mcool_path, "r") as handle:
        if bins_path not in handle:
            raise TadError(
                "TAD_MISSING_RESOLUTION",
                sample=sample,
                signal=signal,
                resolution=resolution,
                path=mcool_path,
            )
        if signal == "balanced" and weight_path not in handle:
            raise TadError(
                "TAD_MISSING_WEIGHT",
                sample=sample,
                signal=signal,
                resolution=resolution,
                path=mcool_path,
            )


def iter_pixels_as_triples(clr, signal):
    if signal not in VALID_SIGNALS:
        raise TadError("TAD_CONFIG_INVALID", signal=signal)
    pixels = clr.pixels()[:]
    weights = None
    if signal == "balanced":
        weights = clr.bins()[:]["weight"].to_numpy(dtype=float)
    for bin1, bin2, count in zip(
        pixels["bin1_id"].to_numpy(),
        pixels["bin2_id"].to_numpy(),
        pixels["count"].to_numpy(dtype=float),
    ):
        if not np.isfinite(count) or count <= 0:
            continue
        value = count
        if weights is not None:
            value = count * weights[int(bin1)] * weights[int(bin2)]
        if np.isfinite(value) and value > 0:
            yield int(bin1) + 1, int(bin2) + 1, float(value)


def write_interactions(path, triples, *, sample, signal, resolution):
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w") as handle:
        for bin1, bin2, value in triples:
            handle.write(f"{bin1}\t{bin2}\t{value}\n")
            count += 1
    if count == 0:
        raise TadError(
            "TAD_EMPTY_INTERACTIONS",
            sample=sample,
            signal=signal,
            resolution=resolution,
            path=path,
        )
    return count
```

Add the required parser argument:

```python
parser.add_argument("--signal", required=True, choices=sorted(VALID_SIGNALS))
```

Use these signal-specific paths and calls in `main`:

```python
validate_mcool_resolution(mcool_path, args.resolution, args.sample, args.signal)
interactions_path = input_dir / f"{args.sample}.{args.signal}.interactions"
settings_path = run_dir / f"settings_{args.sample}.{args.signal}.txt"
n_interactions = write_interactions(
    interactions_path,
    iter_pixels_as_triples(clr, args.signal),
    sample=args.sample,
    signal=args.signal,
    resolution=args.resolution,
)
```

Keep the migrated `effective_d_config` and `render_settings` logic unchanged except that `counts_file` receives the signal-specific filename. Add `signal=<value>` to the ready manifest and publish it atomically:

```python
ready_text = "\n".join([
    f"sample={args.sample}",
    f"signal={args.signal}",
    f"mcool={mcool_path}",
    f"resolution={args.resolution}",
    f"run_dir={run_dir}",
    f"settings={settings_path}",
    f"n_bins={n_bins}",
    f"n_interactions={n_interactions}",
]) + "\n"
ready_tmp = ready_flag.with_name(ready_flag.name + ".tmp")
ready_tmp.write_text(ready_text)
ready_tmp.replace(ready_flag)
```

Move the migrated `main` body after argument parsing into `run_from_args(args)`; keep `build_parser()` as the only parser definition. Use this exact command-line boundary so expected and unexpected preparation failures remain searchable while imported functions still raise `TadError` in tests:

```python
def cli(argv=None):
    args = build_parser().parse_args(argv)
    try:
        run_from_args(args)
    except TadError as error:
        print(error, file=sys.stderr)
        return 1
    except Exception as error:
        wrapped = TadError(
            "TAD_STAGE_FAILED",
            sample=args.sample,
            signal=args.signal,
            resolution=args.resolution,
            stage="prepare",
            returncode="NA",
            path=Path(args.run_dir).resolve(),
            error=type(error).__name__,
        )
        print(wrapped, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
```

- [ ] **Step 5: Run the adapter tests**

Run:

```bash
pytest test/tad/test_tad_from_mcool.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit signal-aware input preparation**

```bash
git add scripts/tad_errors.py scripts/tad_from_mcool.py test/tad/test_tad_from_mcool.py
git commit -m "feat: prepare raw and balanced 3DNetMod inputs"
```

### Task 5: Run 3DNetMod stages and validate final TADs

**Files:**
- Create: `scripts/run_3dnetmod.py`
- Create: `test/tad/test_tad_runner.py`

**Interfaces:**
- Consumes: absolute settings path, vendored code directory, mcool, resolution, sample, signal, and output paths.
- Produces: ordered stage execution, `TAD_STAGE_FAILED` context, canonical sorted BED3/BEDPE6, and `TAD_EMPTY_CALLS` or `TAD_INVALID_COORDINATE` instead of silent files.

- [ ] **Step 1: Write runner tests with fake stages**

```python
# test/tad/test_tad_runner.py
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
        ("# no domains\n", "TAD_EMPTY_CALLS"),
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
```

- [ ] **Step 2: Run the tests and verify the runner is missing**

Run:

```bash
pytest test/tad/test_tad_runner.py -v
```

Expected: FAIL because `scripts/run_3dnetmod.py` does not exist.

- [ ] **Step 3: Copy the validated runner and add dynamic failure context**

Use `/share/home/xhhuang24/project/HiC_analysis_pipeline/scripts/run_3dnetmod.py` as the mechanical base. Import `cooler` and `TadError`, then replace stage execution with:

```python
def run_stages(settings, code_dir, run_dir, *, sample, signal, resolution):
    settings = Path(settings).resolve()
    code_dir = Path(code_dir).resolve()
    run_dir = Path(run_dir).resolve()
    env = build_env(code_dir, run_dir)
    for stage in STAGES:
        try:
            subprocess.run(
                [
                    str(Path(sys.executable).resolve()),
                    str((code_dir / stage).resolve()),
                    str(settings),
                ],
                cwd=run_dir,
                env=env,
                check=True,
            )
        except subprocess.CalledProcessError as error:
            raise TadError(
                "TAD_STAGE_FAILED",
                sample=sample,
                signal=signal,
                resolution=resolution,
                stage=stage,
                returncode=error.returncode,
                path=run_dir,
            ) from error
```

- [ ] **Step 4: Canonicalize, validate, and atomically publish calls**

Replace collection with these exact behaviors:

```python
def read_chrom_sizes(mcool, resolution):
    clr = cooler.Cooler(f"{Path(mcool).resolve()}::/resolutions/{resolution}")
    return {str(chrom): int(length) for chrom, length in clr.chromsizes.items()}


def atomic_write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text)
    temporary.replace(path)


def collect_final_calls(
    run_dir, bed_output, bedpe_output, chrom_sizes,
    *, sample, signal, resolution,
):
    run_dir = Path(run_dir).resolve()
    final_dir = run_dir / "FINAL_DOMAIN_CALLS"
    domain_files = sorted(final_dir.glob("*_master_consistent_domains.txt"))
    chrom_order = {chrom: index for index, chrom in enumerate(chrom_sizes)}
    calls = []
    for domain_file in domain_files:
        for raw_line in domain_file.read_text().splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            columns = line.split()
            try:
                chrom = columns[0]
                start = int(columns[1])
                end = int(columns[2])
            except (IndexError, ValueError) as error:
                raise TadError(
                    "TAD_INVALID_COORDINATE",
                    sample=sample, signal=signal, resolution=resolution,
                    path=domain_file, line=line,
                ) from error
            if (
                chrom not in chrom_sizes
                or start < 0
                or start >= end
                or end > chrom_sizes.get(chrom, -1)
            ):
                raise TadError(
                    "TAD_INVALID_COORDINATE",
                    sample=sample, signal=signal, resolution=resolution,
                    path=domain_file, line=line,
                )
            calls.append((chrom, start, end))
    if not calls:
        raise TadError(
            "TAD_EMPTY_CALLS",
            sample=sample, signal=signal, resolution=resolution,
            path=final_dir,
        )
    calls.sort(key=lambda row: (chrom_order[row[0]], row[1], row[2]))
    bed_text = "".join(f"{chrom}\t{start}\t{end}\n" for chrom, start, end in calls)
    bedpe_text = "".join(
        f"{chrom}\t{start}\t{end}\t{chrom}\t{start}\t{end}\n"
        for chrom, start, end in calls
    )
    atomic_write(bed_output, bed_text)
    atomic_write(bedpe_output, bedpe_text)
```

Replace the CLI entry point with the exact argument and error contract below:

```python
def build_parser():
    parser = ArgumentParser()
    parser.add_argument("--settings", required=True)
    parser.add_argument("--code-dir", required=True)
    parser.add_argument("--mcool", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--signal", required=True, choices=("raw", "balanced"))
    parser.add_argument("--resolution", required=True, type=int)
    parser.add_argument("--bed-output", required=True)
    parser.add_argument("--bedpe-output", required=True)
    return parser


def run_from_args(args):
    settings = Path(args.settings).resolve()
    run_dir = settings.parent
    run_stages(
        settings,
        args.code_dir,
        run_dir,
        sample=args.sample,
        signal=args.signal,
        resolution=args.resolution,
    )
    collect_final_calls(
        run_dir,
        args.bed_output,
        args.bedpe_output,
        read_chrom_sizes(args.mcool, args.resolution),
        sample=args.sample,
        signal=args.signal,
        resolution=args.resolution,
    )


def cli(argv=None):
    args = build_parser().parse_args(argv)
    try:
        run_from_args(args)
    except TadError as error:
        print(error, file=sys.stderr)
        return 1
    except Exception as error:
        wrapped = TadError(
            "TAD_STAGE_FAILED",
            sample=args.sample,
            signal=args.signal,
            resolution=args.resolution,
            stage="runner",
            returncode="NA",
            path=Path(args.settings).resolve().parent,
            error=type(error).__name__,
        )
        print(wrapped, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
```

- [ ] **Step 5: Run the runner tests**

Run:

```bash
pytest test/tad/test_tad_runner.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit the isolated stage runner**

```bash
git add scripts/run_3dnetmod.py test/tad/test_tad_runner.py
git commit -m "feat: run and validate 3DNetMod calls"
```

### Task 6: Generate and merge boundary regions

**Files:**
- Create: `scripts/tad_boundaries.py`
- Create: `test/tad/test_tad_boundaries.py`

**Interfaces:**
- Consumes: canonical TAD BED3, mcool chromosome order/length, integer resolution, sample, signal, and per-resolution boundary BED6 files.
- Produces: per-resolution BED6 and multi-resolution BED7 with deterministic IDs and provenance counts.

- [ ] **Step 1: Write the boundary semantics tests**

```python
# test/tad/test_tad_boundaries.py
from pathlib import Path
import sys

import pytest


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from tad_boundaries import (
    Interval,
    boundaries_from_tads,
    merge_intervals,
    write_multires_boundaries,
    write_resolution_boundaries,
)
from tad_errors import TadError


def interval(start, end, resolution, count=1):
    return Interval("A02", start, end, {resolution}, count)


def test_endpoint_expansion_clips_at_chromosome_edges():
    boundaries = boundaries_from_tads(
        [("A02", 4000, 12000)],
        {"A02": 16000},
        4000,
        sample="S",
        signal="raw",
    )
    assert boundaries == [
        interval(0, 8000, 4000),
        interval(8000, 16000, 4000),
    ]


def test_same_resolution_and_transitive_overlap_take_union():
    merged = merge_intervals(
        [
            interval(0, 8000, 4000),
            interval(4000, 12000, 4000),
            interval(10000, 20000, 8000),
        ],
        {"A02": 0},
    )
    assert merged == [Interval("A02", 0, 20000, {4000, 8000}, 3)]


def test_endpoint_adjacency_does_not_merge():
    assert merge_intervals(
        [interval(0, 8000, 4000), interval(8000, 16000, 8000)],
        {"A02": 0},
    ) == [interval(0, 8000, 4000), interval(8000, 16000, 8000)]


def test_invalid_tad_coordinate_has_dynamic_stable_error():
    with pytest.raises(TadError) as error:
        boundaries_from_tads(
            [("A02", -1, 4000)],
            {"A02": 16000},
            4000,
            sample="actual-sample",
            signal="balanced",
        )
    assert str(error.value).startswith("ERROR[TAD_INVALID_COORDINATE]")
    assert "sample=actual-sample" in str(error.value)
    assert "signal=balanced" in str(error.value)


def test_provenance_counts_survive_second_merge_pass():
    merged = merge_intervals(
        [
            Interval("A02", 0, 12000, {4000}, 2),
            Interval("A02", 10000, 20000, {8000, 20000}, 4),
        ],
        {"A02": 0},
    )
    assert merged == [Interval("A02", 0, 20000, {4000, 8000, 20000}, 6)]


def test_bed6_and_bed7_contracts_have_no_header(tmp_path):
    interval_row = Interval("A02", 8000, 20000, {4000, 8000}, 3)
    bed6 = tmp_path / "resolution.bed"
    bed7 = tmp_path / "multires.bed"
    write_resolution_boundaries(bed6, [interval_row], 4000)
    write_multires_boundaries(bed7, [interval_row])
    assert bed6.read_text() == "A02\t8000\t20000\tboundary_000001\t4000\t3\n"
    assert bed7.read_text() == (
        "A02\t8000\t20000\tboundary_000001\t4000,8000\t2\t3\n"
    )
```

- [ ] **Step 2: Run the tests and verify the module is missing**

Run:

```bash
pytest test/tad/test_tad_boundaries.py -v
```

Expected: collection FAIL because `scripts/tad_boundaries.py` does not exist.

- [ ] **Step 3: Implement the interval model, validation, expansion, and strict overlap union**

```python
# scripts/tad_boundaries.py
from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
import sys

import cooler

from tad_errors import TadError


@dataclass(eq=True)
class Interval:
    chrom: str
    start: int
    end: int
    resolutions: set
    n_source_boundaries: int


def load_chrom_sizes(mcool, resolution):
    clr = cooler.Cooler(f"{Path(mcool).resolve()}::/resolutions/{resolution}")
    return {str(chrom): int(length) for chrom, length in clr.chromsizes.items()}


def read_tads(path, *, sample, signal, resolution):
    rows = []
    for line_number, raw_line in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        columns = line.split()
        try:
            rows.append((columns[0], int(columns[1]), int(columns[2])))
        except (IndexError, ValueError) as error:
            raise TadError(
                "TAD_INVALID_COORDINATE",
                sample=sample,
                signal=signal,
                resolution=resolution,
                path=path,
                line_number=line_number,
                line=line,
            ) from error
    return rows


def boundaries_from_tads(tads, chrom_sizes, resolution, *, sample, signal):
    boundaries = []
    for chrom, start, end in tads:
        if (
            chrom not in chrom_sizes
            or start < 0
            or start >= end
            or end > chrom_sizes.get(chrom, -1)
        ):
            raise TadError(
                "TAD_INVALID_COORDINATE",
                sample=sample, signal=signal, resolution=resolution,
                chrom=chrom, start=start, end=end,
            )
        chrom_length = chrom_sizes[chrom]
        for position in (start, end):
            boundaries.append(Interval(
                chrom,
                max(0, position - resolution),
                min(chrom_length, position + resolution),
                {resolution},
                1,
            ))
    return boundaries


def merge_intervals(intervals, chrom_order):
    ordered = sorted(
        intervals,
        key=lambda item: (chrom_order[item.chrom], item.start, item.end),
    )
    merged = []
    for item in ordered:
        if (
            merged
            and merged[-1].chrom == item.chrom
            and item.start < merged[-1].end
        ):
            merged[-1].end = max(merged[-1].end, item.end)
            merged[-1].resolutions.update(item.resolutions)
            merged[-1].n_source_boundaries += item.n_source_boundaries
        else:
            merged.append(Interval(
                item.chrom, item.start, item.end,
                set(item.resolutions), item.n_source_boundaries,
            ))
    return merged
```

- [ ] **Step 4: Add deterministic BED parsing and atomic writers**

```python
def atomic_write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text)
    temporary.replace(path)


def write_resolution_boundaries(path, intervals, resolution):
    lines = []
    for index, item in enumerate(intervals, 1):
        lines.append(
            f"{item.chrom}\t{item.start}\t{item.end}\t"
            f"boundary_{index:06d}\t{resolution}\t"
            f"{item.n_source_boundaries}\n"
        )
    atomic_write(path, "".join(lines))


def read_resolution_boundaries(path, *, sample, signal):
    intervals = []
    for line_number, line in enumerate(Path(path).read_text().splitlines(), 1):
        columns = line.split("\t")
        try:
            chrom, start, end = columns[0], int(columns[1]), int(columns[2])
            resolution = int(columns[4])
            count = int(columns[5])
        except (IndexError, ValueError) as error:
            raise TadError(
                "TAD_BOUNDARY_MERGE_FAILED",
                sample=sample,
                signal=signal,
                path=path,
                line_number=line_number,
                line=line,
            ) from error
        intervals.append(Interval(chrom, start, end, {resolution}, count))
    return intervals


def write_multires_boundaries(path, intervals):
    lines = []
    for index, item in enumerate(intervals, 1):
        resolutions = sorted(item.resolutions)
        lines.append(
            f"{item.chrom}\t{item.start}\t{item.end}\t"
            f"boundary_{index:06d}\t"
            f"{','.join(map(str, resolutions))}\t{len(resolutions)}\t"
            f"{item.n_source_boundaries}\n"
        )
    atomic_write(path, "".join(lines))
```

Add the exact two-subcommand CLI below. Both paths load chromosome order from the mcool; `extract` performs the same-resolution merge before writing BED6, while `merge` validates every BED6 row against the chromosome sizes before the multi-resolution union.

```python
def build_parser():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    extract = subparsers.add_parser("extract")
    extract.add_argument("--tad", required=True)
    extract.add_argument("--mcool", required=True)
    extract.add_argument("--resolution", required=True, type=int)
    extract.add_argument("--sample", required=True)
    extract.add_argument("--signal", required=True, choices=("raw", "balanced"))
    extract.add_argument("--output", required=True)

    merge = subparsers.add_parser("merge")
    merge.add_argument("--inputs", required=True, nargs="+")
    merge.add_argument("--mcool", required=True)
    merge.add_argument("--chrom-resolution", required=True, type=int)
    merge.add_argument("--sample", required=True)
    merge.add_argument("--signal", required=True, choices=("raw", "balanced"))
    merge.add_argument("--output", required=True)
    return parser


def run_from_args(args):
    chrom_resolution = (
        args.resolution if args.command == "extract" else args.chrom_resolution
    )
    chrom_sizes = load_chrom_sizes(args.mcool, chrom_resolution)
    chrom_order = {chrom: index for index, chrom in enumerate(chrom_sizes)}
    if args.command == "extract":
        intervals = boundaries_from_tads(
            read_tads(
                args.tad,
                sample=args.sample,
                signal=args.signal,
                resolution=args.resolution,
            ),
            chrom_sizes,
            args.resolution,
            sample=args.sample,
            signal=args.signal,
        )
        write_resolution_boundaries(
            args.output,
            merge_intervals(intervals, chrom_order),
            args.resolution,
        )
        return

    intervals = []
    for input_path in args.inputs:
        intervals.extend(
            read_resolution_boundaries(
                input_path,
                sample=args.sample,
                signal=args.signal,
            )
        )
    for item in intervals:
        if (
            item.chrom not in chrom_sizes
            or item.start < 0
            or item.start >= item.end
            or item.end > chrom_sizes.get(item.chrom, -1)
        ):
            raise TadError(
                "TAD_BOUNDARY_MERGE_FAILED",
                sample=args.sample,
                signal=args.signal,
                path=args.output,
                chrom=item.chrom,
                start=item.start,
                end=item.end,
            )
    write_multires_boundaries(
        args.output,
        merge_intervals(intervals, chrom_order),
    )


def cli(argv=None):
    args = build_parser().parse_args(argv)
    try:
        run_from_args(args)
    except TadError as error:
        print(error, file=sys.stderr)
        return 1
    except Exception as error:
        wrapped = TadError(
            "TAD_BOUNDARY_MERGE_FAILED",
            sample=args.sample,
            signal=args.signal,
            path=args.output,
            error=type(error).__name__,
        )
        print(wrapped, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
```

- [ ] **Step 5: Run boundary tests**

Run:

```bash
pytest test/tad/test_tad_boundaries.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit boundary extraction and merging**

```bash
git add scripts/tad_boundaries.py test/tad/test_tad_boundaries.py
git commit -m "feat: merge multiresolution TAD boundaries"
```

### Task 7: Integrate configuration, rules, paths, and the default DAG

**Files:**
- Create: `scripts/tad_config.py`
- Create: `rules/tad.smk`
- Create: `test/tad/test_tad_config.py`
- Create: `test/tad/test_tad_workflow.py`
- Modify: `config/config.yaml:181-187`
- Modify: `Snakefile:4-27,64-72,123-140,146-203`

**Interfaces:**
- Consumes: all scripts and environment from Tasks 2-6 plus the active manifest and mcools.
- Produces: validated `TadSettings`, stable path helpers, four independently retryable rules, and all raw/balanced outputs in `rule all`.

- [ ] **Step 1: Write configuration and workflow tests first**

```python
# test/tad/test_tad_config.py
from copy import deepcopy
from pathlib import Path
import sys

import pytest
import yaml


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from tad_config import TadConfigError, parse_tad_settings


def repository_config():
    return yaml.safe_load((REPO / "config" / "config.yaml").read_text())


def test_repository_tad_contract_is_exact():
    settings = parse_tad_settings(repository_config())
    assert settings.resolutions == (4000, 8000, 20000)
    assert settings.signals == ("raw", "balanced")
    assert settings.threads == 20
    assert settings.excluded_chromosomes == ()
    assert settings.badregionfilter is False
    assert settings.parameter_set == "D"
    assert settings.d["resolution_overrides"][4000]["size_threshold_bins"] == 10


@pytest.mark.parametrize(
    "mutate,match",
    [
        (lambda config: config["tad"].update(resolutions=[4000, 20000]), "resolutions"),
        (lambda config: config["tad"].update(signals=["raw"]), "signals"),
        (lambda config: config["tad"].update(threads=0), "threads"),
        (lambda config: config["tad"].update(parameter_set="A"), "parameter_set"),
        (lambda config: config["tad"]["d"].pop("size_bins"), "size_bins"),
    ],
)
def test_invalid_tad_configuration_has_stable_parse_error(mutate, match):
    config = deepcopy(repository_config())
    mutate(config)
    with pytest.raises(TadConfigError, match=match) as error:
        parse_tad_settings(config)
    assert str(error.value).startswith("ERROR[TAD_CONFIG_INVALID]")
```

```python
# test/tad/test_tad_workflow.py
from pathlib import Path
import subprocess


REPO = Path(__file__).resolve().parents[2]


def test_rule_contract_contains_all_signal_resolution_stages():
    rules = (REPO / "rules" / "tad.smk").read_text()
    snakefile = (REPO / "Snakefile").read_text()
    for rule in (
        "rule prepare_tad:",
        "rule call_tad:",
        "rule extract_tad_boundaries:",
        "rule merge_tad_boundaries:",
    ):
        assert rule in rules
    assert rules.count('../env/tad.yaml') == 4
    assert 'include: "rules/tad.smk"' in snakefile
    for token in (
        "tad_bed_path(sample, resolution, signal)",
        "tad_bedpe_path(sample, resolution, signal)",
        "tad_boundary_path(sample, resolution, signal)",
        "tad_merged_boundary_path(sample, signal)",
    ):
        assert token in snakefile


def test_two_signal_three_resolution_dag_builds():
    completed = subprocess.run(
        [
            "snakemake", "-n", "-p",
            "--configfile", "config/config.yaml",
            "--",
            "results/tad_boundaries/H1394.multires.raw.boundaries.bed",
            "results/tad_boundaries/H1394.multires.balanced.boundaries.bed",
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
    )
    plan = completed.stdout + completed.stderr
    assert completed.returncode == 0, plan
    for rule in (
        "prepare_tad", "call_tad",
        "extract_tad_boundaries", "merge_tad_boundaries",
    ):
        assert rule in plan
    for signal in ("raw", "balanced"):
        for resolution in ("4kb", "8kb", "20kb"):
            assert f"H1394.{resolution}.{signal}.3dnetmod.bed" in plan
```

- [ ] **Step 2: Run the tests and verify configuration/rules are absent**

Run:

```bash
pytest test/tad/test_tad_config.py test/tad/test_tad_workflow.py -v
```

Expected: collection FAIL for `tad_config` and workflow FAIL for missing `rules/tad.smk`.

- [ ] **Step 3: Implement focused TAD configuration validation**

```python
# scripts/tad_config.py
from dataclasses import dataclass


class TadConfigError(ValueError):
    pass


REQUIRED_D_KEYS = {
    "region_size_bins", "overlap_bins", "logged", "scale", "plateau",
    "chaosfilter", "chaos_pct", "max_gamma", "num_part", "plots",
    "pctile_threshold", "pct_value", "final_consistent_domains",
    "size_threshold_bins", "boundary_buffer_bins", "size_bins",
    "variance_thresholds", "resolution_overrides",
}


@dataclass(frozen=True)
class TadSettings:
    resolutions: tuple
    signals: tuple
    threads: int
    parameter_set: str
    run_root: str
    code_dir: str
    badregionfilter: bool
    bad_region_bed: object
    excluded_chromosomes: tuple
    d: dict


def parse_tad_settings(config):
    try:
        section = config["tad"]
        resolutions = tuple(section["resolutions"])
        signals = tuple(section["signals"])
        threads = section["threads"]
        parameter_set = section["parameter_set"]
        d = section["d"]
    except (KeyError, TypeError) as error:
        raise TadConfigError(f"ERROR[TAD_CONFIG_INVALID] missing={error}") from error

    if set(resolutions) != {4000, 8000, 20000} or len(resolutions) != 3:
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] tad.resolutions must be 4000,8000,20000")
    if any(resolution not in config["convert"]["resolutions"] for resolution in resolutions):
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] tad.resolutions require convert.resolutions")
    if set(signals) != {"raw", "balanced"} or len(signals) != 2:
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] tad.signals must be raw,balanced")
    if isinstance(threads, bool) or not isinstance(threads, int) or threads <= 0:
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] tad.threads must be positive")
    if parameter_set != "D":
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] tad.parameter_set must be D")
    missing = REQUIRED_D_KEYS.difference(d)
    if missing:
        raise TadConfigError(
            "ERROR[TAD_CONFIG_INVALID] tad.d missing=" + ",".join(sorted(missing))
        )
    if len(d["size_bins"]) != 5 or len(d["variance_thresholds"]) != 5:
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] D arrays must contain five values")
    excluded = section.get("excluded_chromosomes", [])
    if not isinstance(excluded, list):
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] excluded_chromosomes must be a list")
    badregionfilter = section.get("badregionfilter", False)
    bad_region_bed = section.get("bad_region_bed")
    if badregionfilter and not bad_region_bed:
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] bad_region_bed is required")
    return TadSettings(
        resolutions=tuple(sorted(resolutions)),
        signals=tuple(signal for signal in ("raw", "balanced") if signal in signals),
        threads=threads,
        parameter_set=parameter_set,
        run_root=section["run_root"],
        code_dir=section["code_dir"],
        badregionfilter=badregionfilter,
        bad_region_bed=bad_region_bed,
        excluded_chromosomes=tuple(excluded),
        d=d,
    )
```

- [ ] **Step 4: Add the approved TAD configuration**

Append the complete approved block exactly as written:

```yaml
tad:
  resolutions: [4000, 8000, 20000]
  signals: [raw, balanced]
  threads: 20
  parameter_set: D
  run_root: results/tad_runs
  code_dir: scripts/python3_3dnetmod
  badregionfilter: false
  bad_region_bed: null
  excluded_chromosomes: []
  d:
    region_size_bins: 150
    overlap_bins: 100
    logged: true
    scale: genomewide
    plateau: 3
    chaosfilter: true
    chaos_pct: 0.9
    max_gamma: 3.2
    num_part: 20
    plots: false
    pctile_threshold: 0
    pct_value: 0
    final_consistent_domains: true
    size_threshold_bins: 4
    boundary_buffer_bins: 8
    size_bins: [100, 200, 400, 750, 3000]
    variance_thresholds: [45, 45, 60, 45, 45]
    resolution_overrides:
      4000:
        size_threshold_bins: 10
        variance_thresholds: [20, 20, 30, 20, 20]
```

- [ ] **Step 5: Add Snakefile settings and stable path helpers**

Import `json` and the parser, then define:

```python
from tad_config import parse_tad_settings

TAD = parse_tad_settings(config)
TAD_RESOLUTIONS = TAD.resolutions
TAD_SIGNALS = TAD.signals
TAD_LABELS = {resolution: f"{resolution // 1000}kb" for resolution in TAD_RESOLUTIONS}
TAD_RESOLUTION_BY_LABEL = {label: resolution for resolution, label in TAD_LABELS.items()}


def tad_run_dir(sample, resolution, signal):
    return f"{TAD.run_root}/{sample}/{signal}/{TAD_LABELS[resolution]}/D"


def tad_ready_path(sample, resolution, signal):
    return f"{OUTDIR}/tad_prepare/{sample}.{TAD_LABELS[resolution]}.{signal}.ready"


def tad_bed_path(sample, resolution, signal):
    return f"{OUTDIR}/tads/{sample}.{TAD_LABELS[resolution]}.{signal}.3dnetmod.bed"


def tad_bedpe_path(sample, resolution, signal):
    return f"{OUTDIR}/tads/{sample}.{TAD_LABELS[resolution]}.{signal}.3dnetmod.bedpe"


def tad_boundary_path(sample, resolution, signal):
    return f"{OUTDIR}/tad_boundaries/{sample}.{TAD_LABELS[resolution]}.{signal}.boundaries.bed"


def tad_merged_boundary_path(sample, signal):
    return f"{OUTDIR}/tad_boundaries/{sample}.multires.{signal}.boundaries.bed"
```

Define `TAD_D_CONFIG_JSON = json.dumps(TAD.d, sort_keys=True)` and include `rules/tad.smk` after `rules/is_change.smk`.

- [ ] **Step 6: Add four Snakemake rules**

Create `rules/tad.smk` with this rule interface and the repository's required field order:

```python
wildcard_constraints:
    resolution_label="|".join(TAD_RESOLUTION_BY_LABEL),
    signal="|".join(TAD_SIGNALS)


rule prepare_tad:
    input:
        manifest="data/manifest_mcool.tsv",
        mcool=lambda wc: get_sample_mcool(wc.sample),
    output:
        ready=OUTDIR + "/tad_prepare/{sample}.{resolution_label}.{signal}.ready",
    log:
        OUTDIR + "/logs/tad_prepare/{sample}.{resolution_label}.{signal}.log"
    benchmark:
        OUTDIR + "/benchmark/tad_prepare/{sample}.{resolution_label}.{signal}.tsv"
    threads:
        TAD.threads
    params:
        resolution=lambda wc: TAD_RESOLUTION_BY_LABEL[wc.resolution_label],
        run_dir=lambda wc: tad_run_dir(
            wc.sample, TAD_RESOLUTION_BY_LABEL[wc.resolution_label], wc.signal
        ),
        d_config_json=TAD_D_CONFIG_JSON,
        excluded_chromosomes=",".join(TAD.excluded_chromosomes),
        bad_region_bed=TAD.bad_region_bed or "None",
        badregionfilter=str(TAD.badregionfilter).lower(),
    conda:
        "../env/tad.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.ready})" "$(dirname {log})" "$(dirname {benchmark})"
        python scripts/tad_from_mcool.py \
            --mcool {input.mcool:q} \
            --sample {wildcards.sample:q} \
            --signal {wildcards.signal:q} \
            --resolution {params.resolution} \
            --run-dir {params.run_dir:q} \
            --ready-flag {output.ready:q} \
            --processors {threads} \
            --badregionfilter {params.badregionfilter:q} \
            --bad-region-bed {params.bad_region_bed:q} \
            --excluded-chromosomes {params.excluded_chromosomes:q} \
            --d-config-json {params.d_config_json:q} > {log:q} 2>&1
        """


rule call_tad:
    input:
        ready=lambda wc: tad_ready_path(
            wc.sample, TAD_RESOLUTION_BY_LABEL[wc.resolution_label], wc.signal
        ),
        mcool=lambda wc: get_sample_mcool(wc.sample),
    output:
        bed=OUTDIR + "/tads/{sample}.{resolution_label}.{signal}.3dnetmod.bed",
        bedpe=OUTDIR + "/tads/{sample}.{resolution_label}.{signal}.3dnetmod.bedpe",
    log:
        OUTDIR + "/logs/tad_call/{sample}.{resolution_label}.{signal}.log"
    benchmark:
        OUTDIR + "/benchmark/tad_call/{sample}.{resolution_label}.{signal}.tsv"
    threads:
        TAD.threads
    params:
        resolution=lambda wc: TAD_RESOLUTION_BY_LABEL[wc.resolution_label],
        settings=lambda wc: (
            f"{tad_run_dir(wc.sample, TAD_RESOLUTION_BY_LABEL[wc.resolution_label], wc.signal)}"
            f"/settings_{wc.sample}.{wc.signal}.txt"
        ),
        code_dir=TAD.code_dir,
    conda:
        "../env/tad.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bed})" "$(dirname {log})" "$(dirname {benchmark})"
        python scripts/run_3dnetmod.py \
            --settings {params.settings:q} \
            --code-dir {params.code_dir:q} \
            --mcool {input.mcool:q} \
            --sample {wildcards.sample:q} \
            --signal {wildcards.signal:q} \
            --resolution {params.resolution} \
            --bed-output {output.bed:q} \
            --bedpe-output {output.bedpe:q} > {log:q} 2>&1
        """


rule extract_tad_boundaries:
    input:
        tad=OUTDIR + "/tads/{sample}.{resolution_label}.{signal}.3dnetmod.bed",
        mcool=lambda wc: get_sample_mcool(wc.sample),
    output:
        bed=OUTDIR + "/tad_boundaries/{sample}.{resolution_label}.{signal}.boundaries.bed",
    log:
        OUTDIR + "/logs/tad_boundaries/{sample}.{resolution_label}.{signal}.log"
    benchmark:
        OUTDIR + "/benchmark/tad_boundaries/{sample}.{resolution_label}.{signal}.tsv"
    threads: 1
    params:
        resolution=lambda wc: TAD_RESOLUTION_BY_LABEL[wc.resolution_label],
    conda:
        "../env/tad.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bed})" "$(dirname {log})" "$(dirname {benchmark})"
        python scripts/tad_boundaries.py extract \
            --tad {input.tad:q} --mcool {input.mcool:q} \
            --resolution {params.resolution} \
            --sample {wildcards.sample:q} --signal {wildcards.signal:q} \
            --output {output.bed:q} > {log:q} 2>&1
        """


rule merge_tad_boundaries:
    input:
        boundaries=lambda wc: [
            tad_boundary_path(wc.sample, resolution, wc.signal)
            for resolution in TAD_RESOLUTIONS
        ],
        mcool=lambda wc: get_sample_mcool(wc.sample),
    output:
        bed=OUTDIR + "/tad_boundaries/{sample}.multires.{signal}.boundaries.bed",
    log:
        OUTDIR + "/logs/tad_boundaries/{sample}.multires.{signal}.log"
    benchmark:
        OUTDIR + "/benchmark/tad_boundaries/{sample}.multires.{signal}.tsv"
    threads: 1
    conda:
        "../env/tad.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bed})" "$(dirname {log})" "$(dirname {benchmark})"
        python scripts/tad_boundaries.py merge \
            --inputs {input.boundaries:q} --mcool {input.mcool:q} \
            --chrom-resolution 4000 \
            --sample {wildcards.sample:q} --signal {wildcards.signal:q} \
            --output {output.bed:q} > {log:q} 2>&1
        """
```

- [ ] **Step 7: Extend `rule all` with every stable output**

Define the flattened output list immediately before `rule all`:

```python
TAD_ALL_OUTPUTS = [
    output
    for sample in SAMPLES
    for resolution in TAD_RESOLUTIONS
    for signal in TAD_SIGNALS
    for output in (
        tad_bed_path(sample, resolution, signal),
        tad_bedpe_path(sample, resolution, signal),
        tad_boundary_path(sample, resolution, signal),
    )
] + [
    tad_merged_boundary_path(sample, signal)
    for sample in SAMPLES
    for signal in TAD_SIGNALS
]
```

Then append `TAD_ALL_OUTPUTS,` as the final `input:` entry after the existing IS-change figure manifest expansion. Do not nest another list inside `TAD_ALL_OUTPUTS`.

- [ ] **Step 8: Run config, workflow, and impacted scheduler tests**

Run:

```bash
pytest test/tad/test_tad_config.py test/tad/test_tad_workflow.py test/cworld_is/test_compute_is_scheduling.py -v
snakemake -n -p --configfile config/config.yaml -- results/tad_boundaries/H1394.multires.raw.boundaries.bed results/tad_boundaries/H1394.multires.balanced.boundaries.bed
```

Expected: all tests PASS; dry-run exits 0 and lists six each of `prepare_tad`, `call_tad`, and `extract_tad_boundaries`, plus two `merge_tad_boundaries` jobs unless existing outputs are already current.

- [ ] **Step 9: Commit the workflow integration**

```bash
git add config/config.yaml Snakefile rules/tad.smk scripts/tad_config.py test/tad/test_tad_config.py test/tad/test_tad_workflow.py
git commit -m "feat: schedule multiresolution 3DNetMod TADs"
```

### Task 8: Document, bootstrap, and verify the complete feature

**Files:**
- Create: `CONTEXT.md`
- Create: `docs/adr/0001-3dnetmod-raw-and-balanced-signals.md`
- Modify: `README.md:15-55,85-135,207-227`
- Test: all files in `test/tad/` plus impacted offline and scheduler tests.

**Interfaces:**
- Consumes: completed workflow from Tasks 1-7.
- Produces: user-facing operating instructions, recorded scientific definitions, a prepared login-node environment, a clean full verification run, and one real-sample smoke result.

- [ ] **Step 1: Add the project vocabulary**

```markdown
# SNAP-HiC-Diff Domain Context

## TAD analysis

Per-sample identification of topologically associating domain intervals from a
Hi-C contact map at 4 kb, 8 kb, or 20 kb.

## Raw TAD signal

Finite positive mcool pixel counts used directly as 3DNetMod interaction
values.

## Balanced TAD signal

Finite positive values calculated as pixel count multiplied by both endpoint
bin weights. Balanced calls are separate results, not replacements for raw
calls.

## TAD boundary region

A TAD endpoint expanded left and right by one active resolution and clipped to
the chromosome. Left and right endpoints are not distinguished.

## Multi-resolution boundary region

The transitive interval union of actually overlapping 4 kb, 8 kb, and 20 kb
boundary regions for one sample and one signal. Endpoint adjacency alone does
not merge regions.
```

- [ ] **Step 2: Record the scientific-method ADR**

```markdown
# ADR-0001: Publish raw and balanced 3DNetMod calls separately

## Status

Accepted 2026-08-21.

## Decision

Run 3DNetMod at 4 kb, 8 kb, and 20 kb twice per sample with identical D
parameters. Raw runs use finite positive pixel counts. Balanced runs use
`count * weight(bin1) * weight(bin2)` and skip non-finite or non-positive
results. Boundary merging occurs within each signal and never combines raw
with balanced calls.

## Consequences

The workflow produces two scientifically distinguishable result families and
216 independent calls for the current 36 samples. Existing mcools remain
read-only; missing weights fail rather than trigger implicit rebalancing.
```

- [ ] **Step 3: Add concise operating documentation to README**

Document all of the following with the exact commands and paths:

````markdown
## 3DNetMod TADs and boundary regions

The default workflow runs raw and balanced 3DNetMod D calls at 4 kb, 8 kb,
and 20 kb for every manifest sample. Raw calls use positive pixel counts;
balanced calls use `count * weight1 * weight2`. No chromosome is excluded by
default.

Build or refresh all Conda prefixes on the login node before submission:

```bash
bash scripts/conda_envs.sh bootstrap
CONDA_OFFLINE=true bash scripts/conda_envs.sh verify
```

Compute nodes do not create environments or access Conda channels.

Stable outputs are:

```text
results/tads/{sample}.{4kb|8kb|20kb}.{raw|balanced}.3dnetmod.bed
results/tads/{sample}.{4kb|8kb|20kb}.{raw|balanced}.3dnetmod.bedpe
results/tad_boundaries/{sample}.{4kb|8kb|20kb}.{raw|balanced}.boundaries.bed
results/tad_boundaries/{sample}.multires.{raw|balanced}.boundaries.bed
```

Run one sample by targeting its two merged boundary files; Snakemake will pull
the six required TAD calls:

```bash
snakemake --use-conda --conda-prefix "$PWD/.snakemake/conda" \
  --configfile config/config.yaml -- \
  results/tad_boundaries/H1394.multires.raw.boundaries.bed \
  results/tad_boundaries/H1394.multires.balanced.boundaries.bed
```
````

- [ ] **Step 4: Run all focused and impacted tests**

Run:

```bash
pytest test/tad test/cworld_is/test_offline_conda_execution.py test/cworld_is/test_conda_bootstrap_recovery.py test/cworld_is/test_compute_is_scheduling.py -v
bash -n scripts/conda_envs.sh config/dsub_submit.sh run_snakemake.sh
snakemake -n -p --configfile config/config.yaml -- results/tad_boundaries/H1394.multires.raw.boundaries.bed results/tad_boundaries/H1394.multires.balanced.boundaries.bed
```

Expected: pytest PASS, shell syntax PASS, and Snakemake dry-run exits 0 with both signals and all three resolutions.

- [ ] **Step 5: Bootstrap the real TAD environment on the login node**

Run only where Conda channels are reachable:

```bash
bash scripts/conda_envs.sh bootstrap
CONDA_OFFLINE=true bash scripts/conda_envs.sh verify
```

Expected: `offline.ready.tsv` has four validated records including `env/tad.yaml`; verification prints one `Verified env/tad.yaml -> ...` line.

- [ ] **Step 6: Submit the one-sample six-call smoke target**

Run from the project root after bootstrap:

```bash
snakemake \
  --executor cluster-generic \
  --use-conda \
  --conda-prefix "$PWD/.snakemake/conda" \
  -j 1000 \
  -s Snakefile \
  --configfile config/config.yaml \
  --rerun-incomplete \
  --keep-going \
  --latency-wait 180 \
  --cluster-generic-submit-cmd "bash config/dsub_submit.sh {threads} {rulename}" \
  --cluster-generic-cancel-cmd dkill \
  --cluster-generic-cancel-nargs 1 \
  results/tad_boundaries/H1394.multires.raw.boundaries.bed \
  results/tad_boundaries/H1394.multires.balanced.boundaries.bed
```

Expected: six independent `call_tad` jobs complete, followed by six resolution-specific boundary jobs and two merge jobs.

- [ ] **Step 7: Verify every smoke artifact is non-empty**

Run:

```bash
for signal in raw balanced; do
  for resolution in 4 8 20; do
    test -s "results/tads/H1394.${resolution}kb.${signal}.3dnetmod.bed"
    test -s "results/tads/H1394.${resolution}kb.${signal}.3dnetmod.bedpe"
    test -s "results/tad_boundaries/H1394.${resolution}kb.${signal}.boundaries.bed"
  done
  test -s "results/tad_boundaries/H1394.multires.${signal}.boundaries.bed"
done
```

Expected: exit 0. On failure, search `results/logs/tad_*` for `ERROR[TAD_` and preserve the corresponding `results/tad_runs/H1394/` directory.

- [ ] **Step 8: Commit documentation after verification**

```bash
git add CONTEXT.md docs/adr/0001-3dnetmod-raw-and-balanced-signals.md README.md
git commit -m "docs: document 3DNetMod TAD workflow"
```

- [ ] **Step 9: Run final repository checks without staging generated results**

Run:

```bash
git diff --check
git status --short
```

Expected: no formatting errors; generated `results/`, logs, and run directories remain untracked or ignored and are not staged.
