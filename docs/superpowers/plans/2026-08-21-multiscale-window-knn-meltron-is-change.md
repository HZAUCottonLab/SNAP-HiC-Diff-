# Multi-window KNN, Meltron, and IS-change Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run five resolution-derived genomic window scales through window comparison, KNN calibration, Meltron scoring, cross-scale union, and the existing IS-change caller in one portable Snakemake invocation.

**Architecture:** Parse window arithmetic and background presets once into immutable settings, then expand every window-dependent rule over a `window_bins` wildcard beneath `results/window_scales/r<resolution>/w<window_bins>/`. Shared balanced matrices, O/E inputs, and ten-distance CWorld IS remain single-copy. Per-scale KNN/Meltron candidates are joined with scale-qualified source identities, concatenated into an evidence ledger, merged by genomic overlap/adjacency, and passed to the unchanged IS-change caller beneath `results/multiscale/r<resolution>/`.

**Tech Stack:** Snakemake 8.30, Python 3.12, pandas, NumPy, scikit-learn BallTree/joblib, pytest, YAML, CWorld Perl tools, and Meltron R/Python wrappers.

**Design:** `docs/superpowers/specs/2026-08-21-multiscale-window-knn-meltron-is-change-design.md`

## Global Constraints

- Analyse one matrix resolution per run; default `resolution: 4000`.
- Default windows are `[25, 50, 75, 100, 125]` bins.
- Derive `step_bins = ceil(window_bins * 0.20)`; 3.4 bins is 4 bins.
- Derive all base-pair values from `resolution`; do not encode 4 kb, 100 kb, 200 kb, 300 kb, 400 kb, or 500 kb in workflow logic.
- Require exactly one configured KNN background preset for every active `(resolution_bp, window_bins, step_bins)` tuple.
- Never invoke or copy code from `background_pool`; a missing or invalid preset is a hard error.
- Preserve five raw KNN empirical P values, the any-value `< 0.05` rule, and `meltron_best_score >= 5.0`.
- Do not add multiple-window correction or require support from two genomic window scales.
- Preserve the existing requirement for at least two CWorld IS scale bands in primary IS-change calls.
- Compute the ten CWorld IS distances once per sample and share them across all window scales.
- Qualify every window-dependent result, log, benchmark, and figure path by `r<resolution>/w<window_bins>`.
- Do not delete, alias, or consume legacy unqualified single-window outputs.
- Treat a zero-candidate scale and an all-zero-candidate run as valid header-only outputs; never publish an aggregate after an incomplete scale failure.
- Implement with TDD and one commit per task. During execution, dispatch one fresh implementation subagent per task, followed by specification and code-quality review gates.
- Production dry-run is not accepted until all five configured background files exist. Synthetic 4 kb and 10 kb fixtures provide code-level DAG verification before those files are delivered.

## File Structure

**Create:**

- `scripts/window_analysis_config.py` — immutable scale/preset models and configuration preflight.
- `scripts/validate_knn_background.py` — one-time background schema, coordinate, row-count, and provenance validator.
- `test/multiscale_window/test_window_analysis_config.py` — arithmetic and preset-selection tests.
- `test/multiscale_window/test_validate_knn_background.py` — background contract tests.
- `test/multiscale_window/test_window_generation.py` — scale metadata and coordinate construction tests.
- `test/multiscale_window/test_knn_scale_provenance.py` — KNN and heatmap path/provenance tests.
- `test/multiscale_window/test_multiscale_candidates.py` — cross-scale selection and union tests.
- `test/multiscale_window/test_workflow_dag.py` — synthetic 4 kb/10 kb Snakemake dry-run tests.

**Modify:**

- `config/config.yaml` — replace fixed bp window/step and single background with five bin-derived presets.
- `Snakefile` — parse scales, define scale/multiscale roots, expand final targets over `window_bins`.
- `rules/build_windows.smk` — one window-coordinate table per scale.
- `rules/compute_gap_qc.smk` — one QC checkpoint per `sample x window_bins`.
- `rules/sliding_window.smk` — one HDF5 checkpoint per `sample x window_bins`.
- `rules/window_compare.smk` — scale-qualified covariate, Gaussian, TW, and merged metrics.
- `rules/knn_empirical_pvalue.smk` — matching background audit dependency and scale-qualified KNN output.
- `rules/knn_significant_heatmap.smk` — scale-qualified scan/plot paths.
- `rules/meltron.smk` — one Meltron job per `comparison x window_bins` with shared IS.
- `rules/is_change.smk` — five-scale evidence/union and multiscale IS-change output paths.
- `scripts/build_windows.py` — emit configured scale metadata.
- `scripts/compute_knn_empirical_pvalue.py` — accept explicit config/scale arguments and preserve scale identity.
- `scripts/plot_knn_significant_heatmap.py` — accept shared and scale-specific result roots instead of hard-coded `results/` paths.
- `scripts/prepare_meltron_intervals.py` — use canonical raw KNN candidates and carry source identity.
- `scripts/merge_meltron_scores.py` — restore all source identity columns from the interval map.
- `scripts/is_change_candidates.py` — select scale-qualified windows and merge their cross-scale union.
- `scripts/build_is_change_candidates.py` — accept five aligned KNN/Meltron inputs and write three atomic outputs.
- `test/cworld_is/test_prepare_meltron_intervals.py` — scale identity and empty-schema regression.
- `test/cworld_is/test_merge_meltron_scores.py` — scale identity restoration regression.
- `test/is_change/test_is_change_candidates.py` — scale-qualified joins and aggregate schema.
- `test/is_change/test_is_change_workflow.py` — new multiscale handoff paths.
- `test/cworld_is/test_compute_is_scheduling.py` — scale-qualified Meltron target while retaining shared IS counts.
- `README.md` — configuration, background responsibility, formulas, thresholds, outputs, and resolution examples.

---

### Task 1: Parse window scales and exact background presets

**Files:**

- Create: `scripts/window_analysis_config.py`
- Create: `test/multiscale_window/test_window_analysis_config.py`
- Modify: `config/config.yaml:21-27,78-103`

**Interfaces:**

- Produces: `WindowScale(resolution_bp, window_bins, window_size_bp, step_bins, step_size_bp, tag)`.
- Produces: `BackgroundPreset(resolution_bp, window_bins, step_bins, file)`.
- Produces: `WindowAnalysisSettings(resolution_bp, step_fraction, scales, backgrounds_by_window_bins)`.
- Produces: `parse_window_analysis_settings(config: dict, check_files: bool = True) -> WindowAnalysisSettings`.
- Later tasks import only these public models and parser.

- [ ] **Step 1: Write failing arithmetic and configuration tests**

```python
from pathlib import Path
import sys

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from window_analysis_config import WindowAnalysisConfigError, parse_window_analysis_settings


def make_config(tmp_path, resolution=4000, bins=(25, 50, 75, 100, 125), fraction=0.20):
    presets = []
    for window_bins in bins:
        path = tmp_path / f"r{resolution}_w{window_bins}.tsv"
        path.write_text("comparison\tchrom\twindow_start\twindow_end\n")
        step_bins = __import__("math").ceil(window_bins * fraction)
        presets.append({
            "resolution_bp": resolution,
            "window_bins": window_bins,
            "step_bins": step_bins,
            "file": str(path),
        })
    return {
        "resolution": resolution,
        "window_analysis": {"window_bins": list(bins), "step_fraction": fraction},
        "knn_empirical_pvalue": {"background_presets": presets},
    }


def test_default_4kb_window_math(tmp_path):
    settings = parse_window_analysis_settings(make_config(tmp_path))
    assert [s.window_size_bp for s in settings.scales] == [100000, 200000, 300000, 400000, 500000]
    assert [s.step_bins for s in settings.scales] == [5, 10, 15, 20, 25]
    assert [s.step_size_bp for s in settings.scales] == [20000, 40000, 60000, 80000, 100000]


def test_fractional_step_rounds_up(tmp_path):
    settings = parse_window_analysis_settings(make_config(tmp_path, bins=(17,)))
    assert settings.scales[0].step_bins == 4


def test_10kb_changes_only_physical_sizes(tmp_path):
    settings = parse_window_analysis_settings(make_config(tmp_path, resolution=10000))
    assert [s.window_size_bp for s in settings.scales] == [250000, 500000, 750000, 1000000, 1250000]
    assert [s.step_bins for s in settings.scales] == [5, 10, 15, 20, 25]


def test_missing_exact_preset_names_expected_tuple(tmp_path):
    config = make_config(tmp_path)
    config["knn_empirical_pvalue"]["background_presets"].pop()
    with pytest.raises(WindowAnalysisConfigError, match=r"r4000/w125.*step_bins=25"):
        parse_window_analysis_settings(config)


def test_legacy_only_window_config_is_rejected():
    with pytest.raises(WindowAnalysisConfigError, match="window_analysis"):
        parse_window_analysis_settings({"resolution": 4000, "window_size": 300000, "step_size": 32000})
```

- [ ] **Step 2: Run the tests and verify the module is missing**

Run: `python -m pytest -q test/multiscale_window/test_window_analysis_config.py`

Expected: collection fails with `ModuleNotFoundError: No module named 'window_analysis_config'`.

- [ ] **Step 3: Implement immutable scale parsing**

Implement `scripts/window_analysis_config.py` with these exact public definitions and checks:

```python
from dataclasses import dataclass
from math import ceil, isfinite
from pathlib import Path


class WindowAnalysisConfigError(ValueError):
    pass


@dataclass(frozen=True)
class WindowScale:
    resolution_bp: int
    window_bins: int
    window_size_bp: int
    step_bins: int
    step_size_bp: int

    @property
    def tag(self):
        return f"r{self.resolution_bp}/w{self.window_bins}"


@dataclass(frozen=True)
class BackgroundPreset:
    resolution_bp: int
    window_bins: int
    step_bins: int
    file: Path


@dataclass(frozen=True)
class WindowAnalysisSettings:
    resolution_bp: int
    step_fraction: float
    scales: tuple[WindowScale, ...]
    backgrounds_by_window_bins: dict[int, BackgroundPreset]

    def scale(self, window_bins):
        value = int(window_bins)
        return next(scale for scale in self.scales if scale.window_bins == value)


def parse_window_analysis_settings(config, check_files=True):
    if "window_analysis" not in config:
        raise WindowAnalysisConfigError(
            "window_analysis.window_bins and window_analysis.step_fraction are required; "
            "legacy window_size/step_size are not supported"
        )
    resolution = config.get("resolution")
    bins = config["window_analysis"].get("window_bins")
    fraction = config["window_analysis"].get("step_fraction")
    if isinstance(resolution, bool) or not isinstance(resolution, int) or resolution <= 0:
        raise WindowAnalysisConfigError("resolution must be a positive integer")
    if not isinstance(bins, list) or not bins or any(isinstance(v, bool) or not isinstance(v, int) or v <= 0 for v in bins):
        raise WindowAnalysisConfigError("window_analysis.window_bins must be a non-empty list of positive integers")
    if len(set(bins)) != len(bins):
        raise WindowAnalysisConfigError("window_analysis.window_bins contains duplicates")
    if isinstance(fraction, bool) or not isinstance(fraction, (int, float)) or not isfinite(float(fraction)) or fraction <= 0:
        raise WindowAnalysisConfigError("window_analysis.step_fraction must be finite and positive")
    scales = tuple(
        WindowScale(resolution, window_bins, resolution * window_bins,
                    ceil(window_bins * float(fraction)),
                    resolution * ceil(window_bins * float(fraction)))
        for window_bins in bins
    )
    records = config.get("knn_empirical_pvalue", {}).get("background_presets", [])
    active = {}
    for record in records:
        key = (record.get("resolution_bp"), record.get("window_bins"), record.get("step_bins"))
        if key[0] != resolution:
            continue
        if key in {(resolution, scale.window_bins, scale.step_bins) for scale in scales}:
            if key[1] in active:
                raise WindowAnalysisConfigError(f"duplicate background preset for r{resolution}/w{key[1]}")
            path = Path(record.get("file", ""))
            if check_files and not path.is_file():
                raise WindowAnalysisConfigError(f"background file missing for r{resolution}/w{key[1]}: {path}")
            active[key[1]] = BackgroundPreset(*key, path)
    for scale in scales:
        if scale.window_bins not in active:
            raise WindowAnalysisConfigError(
                f"missing background preset for r{resolution}/w{scale.window_bins} step_bins={scale.step_bins}"
            )
    return WindowAnalysisSettings(resolution, float(fraction), scales, active)
```

- [ ] **Step 4: Replace fixed window and background configuration**

Set `config/config.yaml` to:

```yaml
resolution: 4000
window_analysis:
  window_bins: [25, 50, 75, 100, 125]
  step_fraction: 0.20

knn_empirical_pvalue:
  background_presets:
    - {resolution_bp: 4000, window_bins: 25, step_bins: 5, file: "/share/home/xhhuang24/project/PanTE/background_pool/results/window_scales/r4000/w25/window_compare/HiC_window_background.tsv"}
    - {resolution_bp: 4000, window_bins: 50, step_bins: 10, file: "/share/home/xhhuang24/project/PanTE/background_pool/results/window_scales/r4000/w50/window_compare/HiC_window_background.tsv"}
    - {resolution_bp: 4000, window_bins: 75, step_bins: 15, file: "/share/home/xhhuang24/project/PanTE/background_pool/results/window_scales/r4000/w75/window_compare/HiC_window_background.tsv"}
    - {resolution_bp: 4000, window_bins: 100, step_bins: 20, file: "/share/home/xhhuang24/project/PanTE/background_pool/results/window_scales/r4000/w100/window_compare/HiC_window_background.tsv"}
    - {resolution_bp: 4000, window_bins: 125, step_bins: 25, file: "/share/home/xhhuang24/project/PanTE/background_pool/results/window_scales/r4000/w125/window_compare/HiC_window_background.tsv"}
  knn_k: 10000
```

Keep the existing `knn_metric`, scaler, score columns, derived scores, directions, and `n_jobs` unchanged below these entries. Remove top-level `window_size`, `step_size`, and singular `background_file`.

- [ ] **Step 5: Run focused tests**

Run: `python -m pytest -q test/multiscale_window/test_window_analysis_config.py`

Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add config/config.yaml scripts/window_analysis_config.py test/multiscale_window/test_window_analysis_config.py
git commit -m "feat: define resolution-derived window scales"
```

---

### Task 2: Validate each external KNN background table

**Files:**

- Create: `scripts/validate_knn_background.py`
- Create: `test/multiscale_window/test_validate_knn_background.py`

**Interfaces:**

- Consumes: `BackgroundPreset` and `WindowScale` values from Task 1.
- Produces: `validate_background_table(path, resolution_bp, window_bins, step_bins, knn_k) -> dict`.
- Produces CLI audit JSON with `path`, expected tuple, total/valid row counts, complete/terminal-partial counts, comparisons, and chromosomes.
- Task 5 makes the KNN rule depend on this audit.

- [ ] **Step 1: Write failing contract tests**

Create fixtures with columns:

```python
REQUIRED = [
    "comparison", "chrom", "window_start", "window_end",
    "EC", "IB", "SP", "MSE", "SSIM", "Spearman", "tw_pvalue", "NCC",
]


def background_row(comparison, start, end):
    return {
        "comparison": comparison, "chrom": "chr1",
        "window_start": start, "window_end": end,
        "EC": 1.0, "IB": 2.0, "SP": 3.0, "MSE": 4.0,
        "SSIM": 0.8, "Spearman": 0.7, "tw_pvalue": 0.1, "NCC": 0.6,
    }
```

Cover these assertions:

- the same coordinates under comparisons `rep1` and `rep2` are valid;
- a duplicate `(comparison, chrom, window_start, window_end)` fails;
- missing any required column names that column;
- fewer than `knn_k` finite feature/score rows fails;
- a complete 100 kb window with 20 kb steps passes for `r4000/w25`;
- a non-terminal 96 kb span fails;
- a final 52 kb partial window passes only when it is the last interval in its `comparison x chrom` group;
- an off-grid start or start difference not divisible by 20 kb fails;
- CLI failure leaves a previous audit file unchanged.

- [ ] **Step 2: Run tests and verify failure**

Run: `python -m pytest -q test/multiscale_window/test_validate_knn_background.py`

Expected: collection fails because `validate_knn_background` does not exist.

- [ ] **Step 3: Implement validation and atomic audit output**

Implement the public function with this algorithm:

```python
REQUIRED_COLUMNS = (
    "comparison", "chrom", "window_start", "window_end",
    "EC", "IB", "SP", "MSE", "SSIM", "Spearman", "tw_pvalue", "NCC",
)
NUMERIC_COLUMNS = REQUIRED_COLUMNS[2:]


def validate_background_table(path, resolution_bp, window_bins, step_bins, knn_k):
    frame = pd.read_csv(path, sep="\t", usecols=lambda c: c in REQUIRED_COLUMNS)
    missing = sorted(set(REQUIRED_COLUMNS).difference(frame.columns))
    if missing:
        raise ValueError(f"background missing columns: {', '.join(missing)}")
    for column in NUMERIC_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    logical_keys = ["comparison", "chrom", "window_start", "window_end"]
    if frame.duplicated(logical_keys).any():
        raise ValueError("background has duplicate comparison/chrom/window coordinates")
    expected_span = resolution_bp * window_bins
    expected_step = resolution_bp * step_bins
    if (frame["window_start"] < 0).any() or (frame["window_start"] % resolution_bp != 0).any():
        raise ValueError(f"background starts must align to {resolution_bp} bp")
    complete_count = 0
    partial_count = 0
    for keys, group in frame.groupby(["comparison", "chrom"], sort=False):
        group = group.sort_values(["window_start", "window_end"], kind="stable")
        differences = group["window_start"].diff().dropna()
        if (differences <= 0).any() or (differences % expected_step != 0).any():
            raise ValueError(f"background start grid mismatch for {keys}")
        spans = group["window_end"] - group["window_start"]
        complete = spans.eq(expected_span)
        partial = spans.lt(expected_span) & spans.ge(expected_span / 2)
        if (~complete & ~partial).any():
            raise ValueError(f"background window span mismatch for {keys}")
        if partial.any() and (partial.sum() != 1 or not partial.iloc[-1]):
            raise ValueError(f"partial background window must be terminal for {keys}")
        complete_count += int(complete.sum())
        partial_count += int(partial.sum())
    finite = np.isfinite(frame[list(NUMERIC_COLUMNS)].to_numpy(dtype=float)).all(axis=1)
    valid_rows = int(finite.sum())
    if valid_rows < knn_k:
        raise ValueError(f"background has {valid_rows} valid rows, fewer than knn_k={knn_k}")
    return {
        "path": str(Path(path).resolve()), "resolution_bp": resolution_bp,
        "window_bins": window_bins, "window_size_bp": expected_span,
        "step_bins": step_bins, "step_size_bp": expected_step,
        "total_rows": len(frame), "valid_rows": valid_rows,
        "complete_rows": complete_count, "terminal_partial_rows": partial_count,
        "comparisons": int(frame["comparison"].nunique()),
        "chromosomes": int(frame["chrom"].nunique()),
    }
```

The CLI accepts `--background`, `--resolution-bp`, `--window-bins`, `--step-bins`, `--knn-k`, and `--audit`; write JSON through a temporary file in the destination directory and `os.replace` it only after validation succeeds.

- [ ] **Step 4: Run focused tests**

Run: `python -m pytest -q test/multiscale_window/test_validate_knn_background.py`

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add scripts/validate_knn_background.py test/multiscale_window/test_validate_knn_background.py
git commit -m "feat: validate scale-specific KNN backgrounds"
```

---

### Task 3: Expand window construction, QC, and smoothing by scale

**Files:**

- Modify: `Snakefile:17-27,108-120`
- Modify: `rules/build_windows.smk`
- Modify: `rules/compute_gap_qc.smk`
- Modify: `rules/sliding_window.smk`
- Modify: `scripts/build_windows.py:13-21,59-95,98-142`
- Create: `test/multiscale_window/test_window_generation.py`

**Interfaces:**

- Consumes: `WINDOW_ANALYSIS = parse_window_analysis_settings(config)` and `WINDOW_BINS`.
- Produces: `SCALE_ROOT = OUTDIR + "/window_scales/r<resolution>/w{window_bins}"`.
- Produces per-scale window parquet, QC TSVs, HDF5s, logs, and benchmarks.
- Shared balanced and O/E files stay under existing `results/balanced` and `results/oe` paths.

- [ ] **Step 1: Write failing scale-expansion tests**

Test `build_windows()` with 4 kb/25 bins/5 step bins and assert every row contains:

```python
{
    "resolution_bp": 4000,
    "window_bins": 25,
    "window_size_bp": 100000,
    "step_bins": 5,
    "step_size_bp": 20000,
}
```

Also inspect the three rule files and assert all output/log/benchmark paths contain both `window_scales/r` and `{window_bins}`, while balanced and O/E inputs remain unqualified.

- [ ] **Step 2: Run tests and verify failure**

Run: `python -m pytest -q test/multiscale_window/test_window_generation.py`

Expected: failures show missing scale metadata and missing `window_bins` wildcards.

- [ ] **Step 3: Parse settings once in `Snakefile`**

Replace `WINDOW_SIZE` and `STEP_SIZE` with:

```python
from window_analysis_config import parse_window_analysis_settings

# Mapping completeness is a global preflight check. File existence is left to
# the scale-specific validation rule so unrelated shared-IS targets remain runnable.
WINDOW_ANALYSIS = parse_window_analysis_settings(config, check_files=False)
RESOLUTION = WINDOW_ANALYSIS.resolution_bp
WINDOW_BINS = [scale.window_bins for scale in WINDOW_ANALYSIS.scales]
SCALE_ROOT = "{}/window_scales/r{}/w{{window_bins}}".format(OUTDIR, RESOLUTION)
MULTISCALE_ROOT = "{}/multiscale/r{}".format(OUTDIR, RESOLUTION)


def get_window_scale(window_bins):
    return WINDOW_ANALYSIS.scale(int(window_bins))


def get_background_preset(window_bins):
    return WINDOW_ANALYSIS.backgrounds_by_window_bins[int(window_bins)]
```

- [ ] **Step 4: Add scale metadata to `build_windows.py`**

Add required CLI arguments `--window-bins` and `--step-bins`. Pass them into `build_windows` and append constant `resolution_bp`, `window_bins`, `window_size_bp`, `step_bins`, and `step_size_bp` columns to every row, including terminal partial windows. Keep existing partial-window eligibility unchanged.

- [ ] **Step 5: Rewrite the three rules over `{window_bins}`**

Use `get_window_scale(w.window_bins)` for physical params. Required outputs are:

```text
<SCALE_ROOT>/windows/windows.parquet
<SCALE_ROOT>/windows/{sample}/{sample}__{chrom}_gaps.tsv
<SCALE_ROOT>/windows/{sample}/{sample}__{chrom}_window_qc.tsv
<SCALE_ROOT>/windows/{sample}/{sample}__{chrom}_{stage}.h5
<SCALE_ROOT>/logs/gap_qc/{sample}.done
<SCALE_ROOT>/logs/windows/{sample}.done
```

Change checkpoint lookups to include both `sample=w.sample` and `window_bins=w.window_bins`. Pass derived `window_size_bp` and `step_size_bp` to existing QC/sliding CLIs; no physical window constant remains in the rules.

- [ ] **Step 6: Run focused and existing window tests**

Run: `python -m pytest -q test/multiscale_window/test_window_generation.py`

Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add Snakefile rules/build_windows.smk rules/compute_gap_qc.smk rules/sliding_window.smk scripts/build_windows.py test/multiscale_window/test_window_generation.py
git commit -m "feat: expand window preprocessing by scale"
```

---

### Task 4: Qualify all window-comparison metrics by scale

**Files:**

- Modify: `rules/window_compare.smk`
- Modify: `test/multiscale_window/test_window_generation.py`

**Interfaces:**

- Consumes: per-scale QC/HDF5s and shared balanced coolers.
- Produces one covariate, Gaussian, TW, and merged metrics table per `comparison x window_bins`.
- Leaves metric schemas unchanged until Task 5 adds scale provenance to KNN output.

- [ ] **Step 1: Add failing rule-path assertions**

Assert every output/log/benchmark in `rules/window_compare.smk` is rooted below `SCALE_ROOT`; assert comparison rules expose `{comparison}` and `{window_bins}`; assert Gaussian inputs point to scale-specific HDF5s and TW still receives global `OUTDIR` for balanced coolers.

- [ ] **Step 2: Run the test and verify old paths fail**

Run: `python -m pytest -q test/multiscale_window/test_window_generation.py -k window_compare`

Expected: failures name legacy `results/window_compare` and `results/windows` paths.

- [ ] **Step 3: Rewrite comparison paths and checkpoint lookups**

Use these output forms for every comparison rule:

```text
<SCALE_ROOT>/window_compare/{comparison}/{comparison}_covariates.tsv
<SCALE_ROOT>/window_compare/{comparison}/{comparison}_gaussian_metrics.tsv
<SCALE_ROOT>/window_compare/{comparison}/{comparison}_tw.tsv
<SCALE_ROOT>/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv
```

Update Gaussian checkpoint calls to `checkpoints.sliding_window.get(sample=..., window_bins=w.window_bins)`. Construct explicit HDF5 paths from the same scale root. Keep `params.outdir = OUTDIR` only for TW balanced-cooler lookup. Make `merge_all_window_pair_metrics` scale-qualified or remove it if no target depends on it; do not allow its filesystem scan to mix scales.

- [ ] **Step 4: Run tests**

Run: `python -m pytest -q test/multiscale_window/test_window_generation.py`

Expected: all assertions pass.

- [ ] **Step 5: Commit**

```bash
git add rules/window_compare.smk test/multiscale_window/test_window_generation.py
git commit -m "feat: isolate comparison metrics by window scale"
```

---

### Task 5: Wire matching KNN backgrounds, provenance, and heatmaps

**Files:**

- Modify: `rules/knn_empirical_pvalue.smk`
- Modify: `rules/knn_significant_heatmap.smk`
- Modify: `Snakefile:108-120`
- Modify: `scripts/compute_knn_empirical_pvalue.py:165-227`
- Modify: `scripts/plot_knn_significant_heatmap.py:160-190,302-310,480-560`
- Create: `test/multiscale_window/test_knn_scale_provenance.py`

**Interfaces:**

- Consumes: Task 2 validator and Task 4 merged metrics.
- Produces KNN rows with `resolution_bp`, `window_bins`, `window_size_bp`, `step_bins`, and `step_size_bp`.
- Produces one background audit per scale and scale-qualified heatmap outputs.

- [ ] **Step 1: Write failing KNN provenance tests**

Run `compute_knn_empirical_pvalue.py` on two query rows and a small background with `knn_k: 2`; assert the five scale columns are present and constant, and duplicate `(resolution_bp, window_bins, chrom, window_start, window_end)` query identities fail. Add rule-text assertions that KNN depends on `background_validation.json` and no singular `background_file` remains.

Add a plotting-path unit test that monkeypatches HDF5 loading and proves the CLI uses explicit `--shared-results-root` for balanced coolers and `--window-results-root` for scale-specific HDF5s.

- [ ] **Step 2: Run tests and verify failure**

Run: `python -m pytest -q test/multiscale_window/test_knn_scale_provenance.py`

Expected: missing CLI options, provenance columns, and scale-qualified rule paths.

- [ ] **Step 3: Add the background validation rule**

Create `rule validate_knn_background` in `rules/knn_empirical_pvalue.smk`:

```text
input.background = str(get_background_preset(w.window_bins).file)
output.audit = <SCALE_ROOT>/logs/knn_background/background_validation.json
params = scale resolution/window/step values and knn_k
shell = python scripts/validate_knn_background.py ... --audit {output.audit}
```

Make `rule knn_empirical_pvalue` consume both the external background and audit, then write under `<SCALE_ROOT>/knn_empirical_pvalue/{comparison}/`.

- [ ] **Step 4: Preserve scale identity in KNN output**

Add CLI options `--config`, `--resolution-bp`, `--window-bins`, `--window-size-bp`, `--step-bins`, and `--step-size-bp`. Load KNN options from `args.config`, verify query coordinate identity is complete and unique, append those five constant columns before writing, and order key columns as:

```python
[
    "resolution_bp", "window_bins", "window_size_bp", "step_bins", "step_size_bp",
    "chrom", "window_id", "window_start", "window_end", "comparison",
]
```

Keep BallTree scaling, scores, empirical P-value formula, and `n_jobs` unchanged.

- [ ] **Step 5: Qualify heatmap paths and remove script hard-coding**

Add `window_bins` to checkpoint/rule wildcards and helper calls. Put scans, plots, logs, and benchmarks under the scale root. Replace script-level `PROJECT_ROOT/results/...` lookups with CLI-supplied roots:

```text
--shared-results-root results
--window-results-root results/window_scales/r4000/w25
```

Balanced coolers resolve below the shared root; `ln_oe` and Gaussian HDF5s resolve below the scale root.

- [ ] **Step 6: Run focused tests**

Run: `python -m pytest -q test/multiscale_window/test_validate_knn_background.py test/multiscale_window/test_knn_scale_provenance.py`

Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add Snakefile rules/knn_empirical_pvalue.smk rules/knn_significant_heatmap.smk scripts/compute_knn_empirical_pvalue.py scripts/plot_knn_significant_heatmap.py test/multiscale_window/test_knn_scale_provenance.py
git commit -m "feat: calibrate KNN per window scale"
```

---

### Task 6: Run Meltron per scale and retain source-window identity

**Files:**

- Modify: `rules/meltron.smk`
- Modify: `scripts/prepare_meltron_intervals.py`
- Modify: `scripts/merge_meltron_scores.py`
- Modify: `test/cworld_is/test_prepare_meltron_intervals.py`
- Modify: `test/cworld_is/test_merge_meltron_scores.py`
- Modify: `test/cworld_is/test_compute_is_scheduling.py`

**Interfaces:**

- Consumes: scale-qualified KNN output and shared `results/IS/<sample>/<sample>_IS.tsv`.
- Produces one interval map and merged Meltron table per `comparison x window_bins`.
- The interval map is the authoritative adapter between source window coordinates and Meltron algorithm intervals.

- [ ] **Step 1: Write failing source-identity tests**

Extend preparer fixtures with all five scale columns. Assert IDs are scale-qualified, for example `r4000_w25_HC04_A01_17`, and interval maps contain:

```python
[
    "resolution_bp", "window_bins", "window_size_bp", "step_bins", "step_size_bp",
    "chrom", "interval_start", "interval_end", "window_start", "window_end",
    "window_id", "ID", "n_significant_metrics",
]
```

Assert empty maps keep the identical schema. Extend merge tests to require those scale/source columns in normal and both-directions-empty outputs.

- [ ] **Step 2: Run tests and verify failure**

Run: `python -m pytest -q test/cworld_is/test_prepare_meltron_intervals.py test/cworld_is/test_merge_meltron_scores.py`

Expected: scale identity columns and scale-qualified IDs are missing.

- [ ] **Step 3: Select canonical raw KNN candidates and carry identity**

Change the preparer to require exactly the five `*_empirical_pvalue` columns, calculate `n_significant_metrics = (p < 0.05).sum(axis=1)`, and retain rows meeting `meltron.min_significant` (default 1). This removes derived BH-FDR as the Meltron entry gate while preserving the configured numeric threshold and canonical raw-P-value semantics.

Keep Meltron's strict four-column interval file. Store both source `window_start/window_end` and algorithm-facing `interval_start/interval_end` in the map. Build `ID` from resolution, window bins, chromosome, and window ID.

- [ ] **Step 4: Restore provenance in merged Meltron output**

Join the interval map with `validate="one_to_one"` on `ID` and carry every scale/source field into output. In the empty result branch, emit the complete stable schema. Preserve directional score normalization, NA columns, best-score logic, and `meltron_is_significant` unchanged.

- [ ] **Step 5: Qualify the Meltron rule**

Add `{window_bins}` to KNN input and every Meltron/log/benchmark output under the scale root. Keep `is_a` and `is_b` at shared `results/IS` paths. Pass `--raw-pvalue-threshold 0.05` and keep `--min-significant {params.min_significant}`. Update scheduling tests to target one scale and assert CWorld still schedules only the two samples' per-chromosome IS jobs, not five copies.

In the scheduling test, create five temporary background preset files and a temporary config passed through `--configfile`; never depend on the production background paths being populated during a unit test.

- [ ] **Step 6: Run Meltron and CWorld regressions**

Run: `python -m pytest -q test/cworld_is/test_prepare_meltron_intervals.py test/cworld_is/test_merge_meltron_scores.py test/cworld_is/test_compute_is_scheduling.py`

Expected: all tests pass and shared IS scheduling counts remain unchanged.

- [ ] **Step 7: Commit**

```bash
git add rules/meltron.smk scripts/prepare_meltron_intervals.py scripts/merge_meltron_scores.py test/cworld_is/test_prepare_meltron_intervals.py test/cworld_is/test_merge_meltron_scores.py test/cworld_is/test_compute_is_scheduling.py
git commit -m "feat: run Meltron for each window scale"
```

---

### Task 7: Build the cross-scale candidate ledger and genomic union

**Files:**

- Modify: `scripts/is_change_candidates.py`
- Modify: `scripts/build_is_change_candidates.py`
- Modify: `rules/is_change.smk`
- Modify: `test/is_change/test_is_change_candidates.py`
- Modify: `test/is_change/test_is_change_workflow.py`
- Create: `test/multiscale_window/test_multiscale_candidates.py`

**Interfaces:**

- Consumes: aligned lists of five KNN and five Meltron tables.
- Produces: `all_scale_candidates.tsv`, `candidate_union.tsv`, and `candidate_union.bed`.
- Produces region columns required by the unchanged caller: `region_id`, `chrom`, `region_start`, `region_end`.

- [ ] **Step 1: Write failing multiscale selection tests**

Use one passing w25 window, one overlapping passing w50 window, one touching passing w75 window, one KNN-only w100 window, and one Meltron-only w125 window. Assert:

- only the first three enter the ledger;
- one-window-scale support remains eligible;
- overlap and BED adjacency merge;
- `support_window_bins == "25,50,75"` and `support_scale_count == 3`;
- `source_window_count == 3`;
- the maximum Meltron score and all five minimum empirical P values are correct;
- source identifiers are prefixed by scale;
- chromosome boundaries never merge;
- duplicate scale-qualified identities fail;
- all-empty inputs write stable empty schemas;
- CLI failure preserves all previous output files.

- [ ] **Step 2: Run tests and verify failure**

Run: `python -m pytest -q test/multiscale_window/test_multiscale_candidates.py test/is_change/test_is_change_candidates.py`

Expected: current two-column identity and single-input CLI fail the new assertions.

- [ ] **Step 3: Implement scale-qualified selection**

Set:

```python
SOURCE_IDENTITY_KEYS = [
    "resolution_bp", "window_bins", "chrom", "window_start", "window_end",
]
```

Require both KNN and Meltron to contain that identity; Meltron's separate `start/end` remain audit fields. Join one-to-one on `SOURCE_IDENTITY_KEYS`, require KNN any raw P `< settings.raw_knn_pvalue_threshold`, and require `meltron_best_score >= settings.meltron_score_threshold`.

- [ ] **Step 4: Implement evidence concatenation and region union**

Implement the union deterministically:

1. Build `p_columns` from the five columns ending in `_empirical_pvalue`; reject any other count.
2. Concatenate the per-scale selected frames, or return the complete empty output schema when every scale is empty.
3. Stable-sort by `chrom, window_start, window_end, window_bins`.
4. Sweep once through the sorted rows. Start a new region at a chromosome change or only when `next.window_start > current.region_end`; therefore overlapping and BED-adjacent half-open intervals merge.
5. For every active region, update `region_end=max(region_end, window_end)`, collect unique `window_bins`, append scale-qualified member IDs, increment `source_window_count`, maximize `meltron_best_score`, and minimize each empirical P value.
6. Finalize the active region before opening the next one and once more after iteration; never merge across chromosomes.
7. Assign deterministic one-based `region_id` values after the sweep and serialize sorted scale support as comma-separated `support_window_bins`.
8. Set `support_scale_count` from the number of unique supported scales; do not filter out support from a single scale.

Prefix `member_window_ids` with `r{resolution_bp}_w{window_bins}` and preserve Meltron `ID` values in `member_ids`.
Return `region_id, chrom, region_start, region_end, support_window_bins, support_scale_count, source_window_count, member_window_ids, member_ids, best_meltron_score`, followed by the five `min_*_empirical_pvalue` audit columns.

- [ ] **Step 5: Expand the CLI to three atomic outputs**

Accept `--knn-enriched` and `--meltron-merged` with `nargs="+"`, plus `--expected-window-bins`, `--out-all`, `--out-regions`, `--out-bed`, `--config`, and `--log`. Verify list lengths and observed scale set before selection. Write all outputs to temporary siblings, then replace final paths only after all three serialize successfully.

- [ ] **Step 6: Rewire IS-change rules beneath the multiscale root**

The candidate rule expands all configured scales and writes:

```text
results/multiscale/r4000/{comparison}/all_scale_candidates.tsv
results/multiscale/r4000/{comparison}/candidate_union.tsv
results/multiscale/r4000/{comparison}/candidate_union.bed
```

Move caller, plot, log, and figure outputs beneath `results/multiscale/r4000/{comparison}/is_change/`. Keep shared IS/cool inputs and `call_is_change_loci.py` scoring unchanged.

- [ ] **Step 7: Run candidate and IS-change regressions**

Run: `python -m pytest -q test/multiscale_window/test_multiscale_candidates.py test/is_change`

Expected: all tests pass.

- [ ] **Step 8: Commit**

```bash
git add rules/is_change.smk scripts/is_change_candidates.py scripts/build_is_change_candidates.py test/is_change/test_is_change_candidates.py test/is_change/test_is_change_workflow.py test/multiscale_window/test_multiscale_candidates.py
git commit -m "feat: merge candidates across window scales"
```

---

### Task 8: Expand final targets, document the contract, and verify portability

**Files:**

- Modify: `Snakefile:146-202`
- Create: `test/multiscale_window/test_workflow_dag.py`
- Modify: `README.md:40-70,84-115,150-215`

**Interfaces:**

- Consumes every earlier task.
- Produces the final five-scale DAG and user-facing portable configuration contract.

- [ ] **Step 1: Write failing 4 kb and 10 kb DAG tests**

In pytest, create five valid small background TSVs and a temporary config for each resolution. Invoke Snakemake dry-run against the final `candidate_union.tsv` target. Assert the 4 kb plan contains `w25`, `w50`, `w75`, `w100`, and `w125`; assert the 10 kb plan uses the same bin tags; assert `compute_IS_chrom` counts are not multiplied by five; assert missing one preset fails with the exact expected tuple.

Also assert no final target contains legacy roots `results/window_compare/`, `results/knn_empirical_pvalue/`, `results/meltron/`, or `results/is_change/`.

- [ ] **Step 2: Run DAG tests and verify failure**

Run: `python -m pytest -q test/multiscale_window/test_workflow_dag.py`

Expected: `rule all` still exposes legacy single-window targets.

- [ ] **Step 3: Expand `rule all` over scales**

Use `expand(..., window_bins=WINDOW_BINS, comparison=get_window_compare_pairs())` for window coordinate, gap/sliding markers, comparison stages, KNN, heatmaps, and Meltron. Keep balanced, expected, O/E, and IS targets single-copy. Make final aggregate targets only:

```text
<MULTISCALE_ROOT>/{comparison}/candidate_union.tsv
<MULTISCALE_ROOT>/{comparison}/is_change/{comparison}_is_change_loci.tsv
<MULTISCALE_ROOT>/{comparison}/is_change/figures/manifest.tsv
```

- [ ] **Step 4: Update README**

Document:

- `window_analysis.window_bins` and ceiling step formula;
- the exact 4 kb and 10 kb tables;
- five required background preset records and required TSV columns;
- `background_pool` as producer and SNAP-HiC-Diff as strict consumer;
- missing-preset and invalid-table failure behavior;
- per-scale KNN raw-P-value and Meltron score thresholds;
- no cross-window correction and no minimum window-scale support;
- the distinction between genomic window support and CWorld IS scale-band support;
- scale-qualified outputs, evidence ledger, candidate union, and IS-change outputs;
- shared CWorld IS computation and CPU-driven Donau memory behavior.

- [ ] **Step 5: Run focused tests**

Run: `python -m pytest -q test/multiscale_window`

Expected: all multiscale tests pass.

- [ ] **Step 6: Run all relevant regressions**

Run: `python -m pytest -q test/is_change test/cworld_is`

Expected baseline: at least the pre-implementation `187 passed`; no failures.

- [ ] **Step 7: Run repository-wide tests**

Run: `python -m pytest -q`

Expected: all collected tests pass; existing matplotlib deprecation warnings are acceptable.

- [ ] **Step 8: Verify the real production preflight**

Run: `CONDA_OFFLINE=true snakemake -n -p --use-conda --conda-prefix "$PWD/.snakemake/conda" --configfile config/config.yaml`

Expected when all five configured background files exist: exit 0 and five scale branches in the job summary.

Expected before the background producer delivers all five files: a deliberate error naming the first missing exact `r4000/w<bins>` preset path. Do not weaken validation to make this command pass.

- [ ] **Step 9: Inspect final scope and commit**

Run: `git diff --check`

Run: `git status --short`

Expected: only files listed in this plan are modified; no runtime `results/`, `tmp/`, `.snakemake/`, `.codegraph/`, or scheduler logs are staged.

```bash
git add Snakefile README.md test/multiscale_window/test_workflow_dag.py
git commit -m "docs: document multiscale window workflow"
```

---

## Final Acceptance Gate

- [ ] Five background preset files exist and pass their audit rules.
- [ ] Synthetic 4 kb and 10 kb DAG tests pass.
- [ ] Every comparison expands five KNN and five Meltron branches.
- [ ] CWorld IS remains one per sample/chromosome and one merge per sample.
- [ ] A one-scale-supported candidate survives the cross-scale union.
- [ ] Evidence ledger and region union preserve scale, score, and five-P-value provenance.
- [ ] Existing IS-change multi-band thresholds and classifications remain unchanged.
- [ ] Focused, regression, and repository-wide tests pass.
- [ ] Real `snakemake -n` passes once the five external backgrounds are present.
- [ ] Implementation branch contains only task commits and no unrelated runtime artifacts.
