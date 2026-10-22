# Cross-Layer Thinking Guide

> **Purpose**: Think through data flow across layers before implementing.

---

## The Problem

**Most bugs happen at layer boundaries**, not within layers.

Common cross-layer bugs:
- API returns format A, frontend expects format B
- Database stores X, service transforms to Y, but loses data
- Multiple layers implement the same logic differently

---

## Before Implementing Cross-Layer Features

### Step 1: Map the Data Flow

Draw out how data moves:

```
Source → Transform → Store → Retrieve → Transform → Display
```

For each arrow, ask:
- What format is the data in?
- What could go wrong?
- Who is responsible for validation?

### Step 2: Identify Boundaries

| Boundary | Common Issues |
|----------|---------------|
| API ↔ Service | Type mismatches, missing fields |
| Service ↔ Database | Format conversions, null handling |
| Backend ↔ Frontend | Serialization, date formats |
| Component ↔ Component | Props shape changes |
| Local shell ↔ Remote cluster | Working directory, relative paths, executor placeholder mismatches |

### Step 3: Define Contracts

For each boundary:
- What is the exact input format?
- What is the exact output format?
- What errors can occur?

---

## Common Cross-Layer Mistakes

### Mistake 1: Implicit Format Assumptions

**Bad**: Assuming date format without checking

**Good**: Explicit format conversion at boundaries

### Mistake 2: Scattered Validation

**Bad**: Validating the same thing in multiple layers

**Good**: Validate once at the entry point

### Mistake 3: Leaky Abstractions

**Bad**: Component knows about database schema

**Good**: Each layer only knows its neighbors

---

## Checklist for Cross-Layer Features

Before implementation:
- [ ] Mapped the complete data flow
- [ ] Identified all layer boundaries
- [ ] Defined format at each boundary
- [ ] Decided where validation happens

After implementation:
- [ ] Tested with edge cases (null, empty, invalid)
- [ ] Verified error handling at each boundary
- [ ] Checked data survives round-trip

---

## When to Create Flow Documentation

Create detailed flow docs when:
- Feature spans 3+ layers
- Multiple teams are involved
- Data format is complex
- Feature has caused bugs before

---

## Cluster Execution Cross-Layer Pitfalls (Snakemake + dsub)

When a Snakemake pipeline runs on a remote HPC cluster via dsub, there are **four execution layers** with implicit contracts between them:

```
Local shell (snakemake CLI)
    → Snakemake jobscript (Python, generated on remote node)
    → Bash wrapper (dsub_submit.sh)
    → dsub job (executed on compute node)
```

Each layer can introduce failures the others don't expect:

### Layer 1 → 2: Snakemake → jobscript

- **Contract**: Snakemake generates a Python script and passes the path to the submit command.
- **Contract**: Snakemake generates a Python script and passes its path to the submit command as a **quoted positional argument** (auto-appended by the plugin).
- **Failure mode**: Including `{jobscript}` in the `--cluster-generic-submit-cmd` template causes `NameError` — `{jobscript}` is NOT a valid Python `.format()` variable in Snakemake v8. The jobscript path arrives as `$3`, not as a `{placeholder}`.
- **Detection**: Run `snakemake -j1` on the cluster executor (not dry-run); watch for `NameError` or bash `set -u` failures.

### Layer 2 → 3: jobscript → bash wrapper

- **Contract**: The bash wrapper receives three positional args: `{threads}`, `{rulename}`, and the jobscript path.
- **Failure mode**: If `set -u` is enabled and an arg is missing, the script exits. The wrapper must use `${4:?missing jobscript}` syntax.
- **Detection**: Check that dsub logs show the wrapper received all three args.

### Layer 3 → 4: wrapper → dsub

- **Contract**: dsub runs the wrapper with `--cwd` set to the project root.
- **Failure mode**: Relative paths in `configfile:` or Snakefile `include:` resolve from the wrong directory. `set -e` causes immediate exit on any non-zero return.
- **Detection**: Always use `--directory "$PROJECT_ROOT"` in Snakemake invocation and absolute paths for `--configfile`.

### Layer 4 → file system

- **Contract**: Remote nodes have the same filesystem view as the submit node.
- **Failure mode**: Symlinks, NFS mounts, or `$HOME` expansion differ on remote nodes. Config files not copied to shared filesystem.
- **Detection**: Use absolute paths throughout; verify all referenced files exist on a shared filesystem accessible from compute nodes.

### Checklist for Cluster Pipeline Changes

Before submitting:
- [ ] Used `--directory "$PROJECT_ROOT"` (not just `--cwd` in dsub)
- [ ] Used absolute path for `--configfile` (not relative `config/config.yaml`)
- [ ] Verified submit command does NOT include `{jobscript}` in the template (plugin auto-appends it; including it causes `NameError`)
- [ ] Tested with `snakemake -j1 --executor cluster-generic` on a real job (not `-n`)
- [ ] Checked dsub logs from a compute node for the actual error

---

## Case Study: `{jobscript}` NameError in cluster-generic (Apr 2026)

### Symptom

Two independent `FileNotFoundError` failures on consecutive dsub job runs:

```
# Run 1 — wrong working directory
FileNotFoundError: '/share/home/xhhuang24/project/SNAP-HiC-Diff/config/config.yaml'
# (note: missing "SNAP-HiC-Diff" subdirectory in path)

# Run 2 — after fixing directory, but {jobscript} still present
NameError: The name 'jobscript' is unknown in this context.
```

### Root Cause

`run_snakemake.sh` included `{jobscript}` in the `--cluster-generic-submit-cmd` template:

```bash
# WRONG — causes NameError
--cluster-generic-submit-cmd 'bash dsub_submit.sh {threads} {rulename} {jobscript}'
```

Snakemake's cluster-generic plugin uses Python's `str.format()` with a **fixed variable set**: `{threads}`, `{rulename}`, `{rule}`, `{jobid}`, `{profile}`, `{resources.*}`, `{attempt}`. The `{jobscript}` variable is **not in that set**. The plugin additionally auto-appends the jobscript path as a quoted positional argument to the shell command — it arrives as `$3`, bypassing the template entirely.

### Why the Bug Persisted

1. **Spec was updated but implementation was not**: `snakemake-guidelines.md` already documented the correct pattern, and `cross-layer-thinking-guide.md` already had a checklist item for it. But `run_snakemake.sh` was never synchronized.
2. **Local dry-run cannot detect it**: `snakemake -n` never exercises the cluster-generic path.
3. **Two independent fixes needed**: First the `--directory` flag was needed to fix `FileNotFoundError`, then the `{jobscript}` removal was needed to fix `NameError`. They were separate bugs in the same command template.

### Fix Applied

```bash
# CORRECT — positional args only, plugin auto-appends jobscript as $3
--cluster-generic-submit-cmd 'bash dsub_submit.sh {threads} {rulename}'
```

The `dsub_submit.sh` wrapper captures the jobscript path with:

```bash
_jobscript="${3:?missing jobscript}"
```

### Prevention

- **Code review checklist**: When changing `--cluster-generic-submit-cmd`, verify it only uses the known placeholder set.
- **CI/audit check**: A grep for `{jobscript}` in all shell scripts should be part of the pipeline health check.
- **Spec-implementation sync**: Before deploying spec changes, run `grep -r "cluster-generic-submit-cmd" *.sh` to find all call sites that must be updated.

---

## Case Study: `checkpoint()` at Module Level + Missing conda Declarations (May 2026)

### Symptom

```
NameError: name 'checkpoint' is not defined
```

Snakemake DAG construction failed with `checkpoint` not defined, then subsequent fixes introduced more errors.

### Root Cause

Two independent bugs:

1. **`checkpoint()` at module level**: `checkpoint()` is a Snakemake DSL function only available inside rule execution contexts (input functions, params, shell blocks). It was called at module level (top of `.smk` file), where Snakemake's rule infrastructure is not initialized.
2. **Missing `conda:` on all rules**: None of the 12 rules had `conda:` declarations. `--conda-create-envs-only` built the DAG but created no environments.

### Key Distinction: `checkpoint` vs `checkpoints`

| Name | Type | Where to use |
|------|------|-------------|
| `checkpoint` | Singular function, used inside `output:` | Inside rule's `output:` to declare a checkpoint output |
| `checkpoints` | Plural module-level object | Inside input functions (lambdas or named functions) |

### Why the Bug Persisted

1. **Four failed fix attempts** before finding the correct pattern. Each attempt revealed a new misunderstanding:
   - Attempt 1: Plain path string — unresolved `{chunk_id}` wildcard remained
   - Attempt 2: Input function reading manifest — couldn't resolve output wildcards
   - Attempt 3: `ancient(checkpoint(...))` — `ancient()` also evaluates at module level
   - Attempt 4: `checkpoints` inside lambda input function — finally worked
2. **Second independent bug**: `--conda-create-envs-only` building DAG without creating environments was not immediately recognized as a missing `conda:` issue.

### Fix Applied

```python
# Correct: checkpoints inside input function
def get_metric_chunks(wildcards):
    manifest_path = checkpoints.chunk_windows.get(pair_id=wildcards.pair_id).output.manifest
    df = pd.read_csv(manifest_path, sep="\t")
    return [f"{OUTDIR}/metrics/{wildcards.pair_id}/chunk_{cid:03d}_metrics.tsv"
            for cid in df["chunk_id"].tolist()]

rule aggregate_metrics:
    input:
        get_metric_chunks  # No parentheses
```

And for every rule:

```python
rule my_rule:
    input: ...
    output: ...
    conda:
        "../env/SNAP-HiC-Diff.yaml"  # Required for --conda-create-envs-only
    script:
        "scripts/my_rule.py"
```

### Prevention

- **Lint rule**: `grep -n "checkpoint(" rules/*.smk` before committing; `checkpoint()` should only appear inside rule `output:` blocks, never at module level.
- **Conda checklist**: Every new rule must have `conda:` before first commit.
- **Template enforcement**: The rule template in Section 17 must list `conda` as a required field.

---

## Case Study: 1D/2D Mask Semantic Confusion in Hi-C Gap Interpolation (May 2026)

### Symptom

Two independent failures on real data, not caught by unit tests:

1. **Without gaps**: 13/25 positions in a 5×5 matrix were incorrectly zeroed.
2. **With gaps**: interpolation source came from NaN region instead of valid values.

### Root Cause

`gap_bin_mask` is a **1D** array (per-bin gap indicator), but the 2D Hi-C matrix interpretation is:

```
gap_bin_mask[k] = True
        ↓
row k AND column k are gap regions
        ↓
matrix position (i, j) is gap iff gap_bin_mask[i] == True OR gap_bin_mask[j] == True
```

**Bug #1** (first attempt): Used `gap_bin_mask` directly in a 2D loop as if it were 2D — `gap_bin_mask[i, j]` is invalid for a 1D array. The correct construction requires broadcasting and run-length analysis.

**Bug #2** (second attempt): After fixing the mask construction, used `scipy.interpolate.RegularGridInterpolator`. The gap row's NaN values were set to 0 as placeholder — RGI's `method='nearest'` then returned 0 as the nearest value (treating 0 as the nearest neighbor, not as an invalid marker).

### Why Unit Tests Passed But Integration Failed

Unit tests used uniform matrices (`np.ones()`) where original value == interpolation value == 1.0. The bug was invisible because the "wrong" value happened to equal the "right" value.

**Lesson**: Test with heterogeneous values, not uniform ones. At minimum, test with a gradient or known pattern that makes "right vs wrong" distinguishable.

### Correct Fix

```python
# Step 1: Convert 1D gap_bin_mask → 2D masks via run-length analysis
def _build_2d_gap_masks(gap_bin_mask, n):
    """1D → 2D: gap_bin_mask[k]=True → rows [k] AND cols [k] are gap."""
    fill_cand = np.zeros((n, n), dtype=bool)
    hard_mask = np.zeros((n, n), dtype=bool)
    j = 0
    while j < n:
        if gap_bin_mask[j]:
            start = j
            while j < n and gap_bin_mask[j]:
                j += 1
            run_len = j - start
            target = fill_cand if run_len <= 2 else hard_mask
            target[start:j, :] = True   # row side
            target[:, start:j] = True   # column side
        else:
            j += 1
    return fill_cand, hard_mask

# Step 2: Explicit nearest neighbor (don't rely on RGI treating 0 as valid)
def _interpolate_2d_nearest(fill_input, target_mask, fallback=0.0):
    """For each target, find closest non-zero source."""
    source_mask = (fill_input != 0.0)  # zeros = placeholder, not valid
    ...
```

### Prevention

- **Docstring contract**: Every function that takes a 1D mask for a 2D matrix must document the exact 2D interpretation of each bit.
- **Non-uniform test data**: Never test matrix functions with all-ones or all-zeros matrices.
- **Integration test on real data**: Unit tests are necessary but not sufficient; always add at least one real-data smoke test.

---

## Case Study: Pandas Merge Column Collision + Metadata Loss (May 2026)

### Symptom

Three sequential failures on the same pipeline:

1. `ValueError: The column label 'window_start' is not unique.`
2. `KeyError: "['chrom'] not in index"`
3. `KeyError: "['comparison', 'sample_a', 'sample_b'] not in index"`

### Root Cause

The merge pattern was copy-pasted across three files. Each fix addressed one symptom without propagating to the other two:

| Submission | Fix attempted | Symptom addressed |
|-----------|---------------|-----------------|
| 36521934 | Added `and c not in keys` filter | ValueError fixed in merge |
| 36522171 | `set_index(keys)` + `right_index=True` | KeyError in gauss/tw |
| 36522386 | Fixed bracket syntax | KeyError on missing metadata |

### Layer Model

```
Test fixture (pytest)
    ↓ calls
Scripts 1-4 (compute_window_pair_covariates, compare_gaussian_smoothed_windows,
              compare_iced_windows_tw, merge_window_pair_metrics)
    ↓ produce/consume
TSV files (covariates, gaussian_metrics, tw_metrics)
    ↓ join via
pandas.merge() ← BREAKS at three points
```

### Prevention

- **P0**: After fixing a pattern used in N files, grep for same pattern in all N files before declaring done.
- **P0**: Never use `os.system()` in tests — always `subprocess.run()` (exit code semantics differ).
- **P0**: When fixing merge logic, check all scripts that consume the same output format.
- **P1**: Extract shared merge logic to `scripts/utils/df_transforms.py`.
- **P1**: Add integration tests that call `main()` end-to-end for all pipeline scripts.
- **P2**: Add warning filter in tests for known acceptable RuntimeWarnings (e.g., `np.log(0)`).

---

## Case Study: Chromosome Naming Mismatch — Working Tree vs Git HEAD (Jun 2026)

### Symptom

36 个 `extract_and_balance_sample` 样本全部失败，同一错误：

```
[ERROR] Chromosome Ghir_A01 not found in /path/to/L282.mcool::/resolutions/4000
exit 1
```

### Root Cause

**Working tree 变更对 dsub 不可见。**

`config/chroms.yaml` committed 版本使用 `Ghir_*`，working tree 改为 `HC04_*`（mcool 实际染色体名），但未 commit。

```
Working tree:  chroms.yaml = HC04_A01   (← mcool 实际命名)
Git HEAD:      chroms.yaml = Ghir_A01   (← dsub 从此构建)
dsub job:      读取 Git HEAD → for CHR in Ghir_A01...
                  → awk '$1 == Ghir_A01' → no match → exit 1
```

`rules/extract_chrom.smk` shell block 用 `for CHR in {params.chroms}`，本身正确（参数化），但参数值来自 committed `chroms.yaml`。

### Four Execution Layers (Re-enforced)

```
Local shell (git commit)
    → Snakemake DAG (reads committed config)
    → dsub job (remote node, committed files only)
    → hictk / awk (reads actual mcool chromosome names)
```

**Layer 3→4 contract**: `hictk dump --table chroms` 输出的染色体名必须和 shell `for CHR in ...` 的值完全一致。

### Why the Bug Persisted

1. **Working tree 变更对 dsub 不可见**：dsub 从 git HEAD 构建环境，未 commit 变更完全忽略
2. **染色体命名约定未文档化**：pipeline 没有记录每批数据的染色体命名来源
3. **两次修改未原子 commit**：commit `c83f43b` 同时改了 `chroms.yaml` 和 `extract_chrom.smk`，但都基于旧命名体系
4. **无染色体名预验证**：Snakemake DAG 构建阶段不检查 `hictk dump` 能否匹配染色体名

### Prevention

| Priority | Action |
|----------|--------|
| P0 | 提交 dsub job 前：运行 `git status` 确认 working tree 无相关变更 |
| P0 | 同一 commit 中修改 config 和对应 rule，确保两者同步 |
| P0 | 提交后立即验证：取第一个样本，验证 `hictk dump --table chroms` 的输出列名是否在 `chroms.yaml` 中 |
| P1 | 在 manifest TSV 中增加染色体命名列（`chrom_naming: HC04`）供 pipeline 自动推断 |
| P1 | 建立数据源清单，记录每批 rawdata 的染色体命名约定 |
| P2 | Snakemake rule 在 shell block 开头加染色体名预检查 |

**Pre-submission checklist for config changes:**

```
Before `dsub -s run_snakemake.sh`:
- [ ] git status → no uncommitted changes in config/ or rules/
- [ ] Verify config values match actual data:
      hictk dump {input} --table chroms | head -1  → must match chroms.yaml[0]
- [ ] snakemake -n on cluster (not local dry-run)
```
