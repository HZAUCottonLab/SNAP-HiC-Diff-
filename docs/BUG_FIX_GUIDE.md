# Bug Fix Guide — SNAP-HiC-Diff Pipeline

> Practical reference for diagnosing and fixing common pipeline failures. Each entry covers symptom, root cause, fix, and verification.

---

## B-01: Chromosome Naming Mismatch

**Severity**: CRITICAL
**Symptoms**: All 36 `extract_and_balance_sample` jobs fail with `exit 1`, log shows:

```
[ERROR] Chromosome Ghir_A01 not found in /path/to/sample.mcool::/resolutions/4000
```

**Root Cause**: `config/chroms.yaml` uses a different chromosome naming convention than the actual mcool files.

This project has used two naming conventions:

| Convention | Example | Used by |
|-----------|---------|---------|
| `Ghir_*` | `Ghir_A01` | Old `Ghir_Pyb` mcool files |
| `HC04_*` | `HC04_A01` | Current `compare_mtx_v2` mcool files |

The rule shell block loops over `for CHR in {params.chroms}`, and `hictk dump --table chroms` outputs the actual names from the mcool. If they don't match, awk produces an empty file → `exit 1`.

**Why dsub amplifies this**: dsub runs from git HEAD, NOT the working tree. A working-tree-only change to `chroms.yaml` is invisible to dsub.

**Fix**:

```bash
# Step 1: Check what names the mcool files actually use
hictk dump /path/to/your/sample.mcool::/resolutions/4000 --table chroms | head -3

# Step 2: Verify chroms.yaml matches
cat config/chroms.yaml

# Step 3: If they differ, update chroms.yaml AND commit
# (dsub won't see working-tree-only changes)
git add config/chroms.yaml
git commit -m "fix: align chroms.yaml to actual mcool chromosome naming"

# Step 4: Clean stale outputs from failed run
rm -f results/logs/extract_balance/*.log results/logs/extract_balance/*.done
rm -f tmp/*__Ghir_*   # stale chrom.sizes files

# Step 5: Re-submit
dsub -s run_snakemake.sh
```

**Verification**:

```bash
# Check the first job's output paths in the .err log
grep "output:" job_XXXX.err | head -1
# Should contain HC04_A01, not Ghir_A01
```

**Prevention**: Before any dsub submission:

```bash
git status   # Must be clean for relevant files
# OR commit all changes first

# Verify chromosome names match
hictk dump "$(head -1 data/manifest_mcool.tsv | cut -f3)::/resolutions/4000" \
    --table chroms | head -1
# Output must match line 1 of config/chroms.yaml
```

---

## B-02: Snakemake checkpoint Missing Upstream Dependencies

**Severity**: HIGH
**Symptoms**: `FileNotFoundError` for `.npz` files during `sliding_window` execution:

```
FileNotFoundError: results/oe/H1398_rep2__HC04_A01_oe.npz
```

**Root Cause**: A `checkpoint` rule declares only some of its upstream inputs. When `checkpoint` and a producer rule run in parallel, the producer may not finish before the checkpoint starts.

**Fix**: Ensure every `checkpoint` explicitly declares all upstream dependencies:

```python
# WRONG: missing oe_files dependency
checkpoint sliding_window:
    input:
        gaps = expand(..., chrom=CHROMS),
        script = "scripts/sliding_window.py"

# CORRECT: both gaps and oe_files declared
checkpoint sliding_window:
    input:
        gaps = expand(..., chrom=CHROMS),
        oe_files = expand(OUTDIR + "/oe/{{sample}}__{chrom}_oe.npz", chrom=CHROMS),
        script = "scripts/sliding_window.py"
```

**Prevention**: When adding an input to a rule that a `checkpoint` depends on, always update the `checkpoint`'s `input:` block.

---

## B-03: Pandas Merge Column Collision

**Severity**: MEDIUM
**Symptoms**: Sequential failures on `compare_gaussian_smoothed_windows`, `compare_iced_windows_tw`:

```
ValueError: The column label 'window_start' is not unique
KeyError: "['chrom'] not in index"
KeyError: "['comparison', 'sample_a', 'sample_b'] not in index"
```

**Root Cause**: Three independent merge points all had the same pattern bug. After fixing one, the others still failed.

**Fix**: The merge pattern was copy-pasted across four files. After fixing one merge, grep for the same pattern in all related scripts:

```bash
# After fixing merge logic in one script, check all consumers
grep -n "pd.merge\|merge(" scripts/compare_*.py scripts/merge_*.py
```

**Prevention**: P0 rule — after fixing a pattern in N files, verify the same pattern across all N files before declaring done.

---

## B-04: checkpoint() Called at Module Level

**Severity**: HIGH
**Symptoms**:

```
NameError: name 'checkpoint' is not defined
```

**Root Cause**: `checkpoint()` is a Snakemake DSL function available only inside rule execution contexts. Calling it at the top of a `.smk` file (module level) fails because Snakemake's rule infrastructure is not initialized at parse time.

**Fix**: Always use `checkpoints` (plural, module-level object) inside input functions:

```python
# WRONG: checkpoint at module level
input:
    ancient(checkpoint(OUTDIR + "/chunks/{pair_id}/manifest.tsv")),

# CORRECT: checkpoints inside input function
def get_chunks(wildcards):
    return checkpoints.chunk_windows.get(
        pair_id=wildcards.pair_id
    ).output.manifest

rule aggregate:
    input:
        get_chunks  # no parentheses — pass the function object
```

---

## B-05: Cluster-generic {jobscript} Placeholder Trap

**Severity**: HIGH
**Symptoms**:

```
NameError: The name 'jobscript' is unknown in this context
FileNotFoundError: config/config.yaml (missing subdirectory in path)
```

**Root Cause**: Two independent issues in `--cluster-generic-submit-cmd`:

1. `{jobscript}` is NOT a valid Snakemake placeholder (the plugin auto-appends the jobscript path as a positional argument `$3`)
2. Relative paths in `--configfile` resolve from the wrong directory on remote nodes

**Fix**:

```bash
# CORRECT: positional args only, plugin auto-appends jobscript
--cluster-generic-submit-cmd 'bash dsub_submit.sh {threads} {rulename}'

# WRONG: {jobscript} causes NameError
--cluster-generic-submit-cmd 'bash dsub_submit.sh {threads} {rulename} {jobscript}'
```

Always use `--directory "$PROJECT_ROOT"` and absolute paths for `--configfile`:

```bash
PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
snakemake \
    --directory "$PROJECT_ROOT" \
    --configfile "$PROJECT_ROOT/config/config.yaml" \
    ...
```

**Valid placeholders**: `{threads}`, `{rule}`, `{rulename}`, `{jobid}`, `{profile}`, `{resources.*}`, `{attempt}`

**Prevention**: Before each deployment, `grep -rn '{jobscript}' *.sh` to detect this pattern.

---

## B-06: Shell Directive + config Object NameError

**Severity**: MEDIUM
**Symptoms**:

```
NameError: name 'config' is not defined
```

**Root Cause**: A Python script uses `config.get()` inside a `main()` function, but `config` is only defined inside an `if "snakemake" in dir()` block. When the rule uses `shell:` directive (not `script:`), `snakemake` is not injected → the `if` block is skipped → `config` is never defined.

**Fix**: Always define `config = {}` in `main()` before any `config.get()` calls:

```python
def main(args=None):
    if args is None:
        args = sys.argv[1:]
    config = {}   # Define before any config.get() calls
    ap = argparse.ArgumentParser()
    mcool_path = config.get("rawdata_dir", "")  # safe now
```

**Prevention**: If a script uses `snakemake` or `config`, use `script:` directive instead of `shell:`.

---

## Pre-Submission Checklist

Before running `dsub -s run_snakemake.sh`:

- [ ] `git status` — no uncommitted changes in `config/`, `rules/`, `scripts/`
- [ ] Chromosome naming verified: `hictk dump {sample}.mcool --table chroms | head -1` matches `config/chroms.yaml`
- [ ] Snakemake dry-run on cluster executor (not local `-n`):
  ```bash
  snakemake --executor cluster-generic -n -j1  # test one job on cluster
  ```
- [ ] All `checkpoint` rules have complete `input:` declarations
- [ ] No `{jobscript}` placeholder in `run_snakemake.sh`
- [ ] `--configfile` uses absolute path
- [ ] `--directory` set to project root
