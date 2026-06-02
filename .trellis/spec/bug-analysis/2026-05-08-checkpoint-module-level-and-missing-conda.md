# Bug Analysis: checkpoint() at Module Level + Missing conda Declarations (2026-05-08)

## 1. Root Cause Category

- **Category**: B — Cross-Layer Contract (Snakemake Python DSL ↔ Rule Evaluation Context)
- **Specific Cause**: `checkpoint()` is a Python function callable only **inside** rule execution contexts (input functions, params, shell blocks). At module level (top of `.smk` file), Snakemake's rule infrastructure is not initialized, so `checkpoint` is not in scope.

## 2. Why Fixes Failed

1. **First attempt**: Changed `checkpoint(...)` to a plain path string — partially worked but failed because `score_metrics` rule also had an unresolved `{chunk_id}` wildcard in its output with no mechanism to expand it from the checkpoint.
2. **Second attempt**: Added `get_score_input_chunks()` input function that reads the manifest — failed because input functions can't resolve output wildcards; the `{chunk_id}` in `score_metrics` output remained unresolvable.
3. **Third attempt**: Used `ancient(checkpoint(...))` — failed because `ancient()` is a Python function that evaluates its argument at module level, and `checkpoint()` is still not callable there.
4. **Final fix**: Used `lambda wildcards: checkpoints.chunk_windows.get(...).output.manifest` inside the rule's `input:` directive. The `checkpoints` module-level object is available for use inside input functions. The `{chunk_id}` in output was resolved because the `checkpoint` directive in `chunk_windows` tells Snakemake to expand `{chunk_id}` from the checkpoint's output.

## 3. Prevention Mechanisms

| Priority | Mechanism | Specific Action | Status |
|----------|-----------|----------------|--------|
| P0 | Documentation | Add checkpoint usage rules to `snakemake-guidelines.md` Section 10 (Input functions) | DONE |
| P0 | Documentation | Add explicit "module-level vs rule-context" distinction for `checkpoint()` | DONE |
| P1 | Code Review Checklist | Require `conda:` declaration on every rule; verify via `grep -l "rule " *.smk \| xargs grep -L "conda:"` | TODO |
| P1 | Spec Sync | After creating any new rule, ensure it has `conda:` before first commit | TODO |

## 4. Systematic Expansion

- **Similar Issues**: All Snakemake functions that are only valid inside rule context (`checkpoints`, `expand` with unresolved wildcards, `ancient`) will fail at module level. This pattern applies to any new pipeline that introduces checkpoints.
- **Design Improvement**: Consider a Snakemake linting rule that checks every `.smk` file for: (a) missing `conda:` directive on all non-comment rules, (b) `checkpoint()` calls that are not inside input functions.
- **Process Improvement**: When adding a new rule, the template must include `conda:` as a required field (already in `snakemake-guidelines.md` Section 17 but not enforced).

## 5. Knowledge Capture

### Bug Analysis Document Created

- [x] Created `.trellis/spec/bug-analysis/2026-05-08-checkpoint-module-level-and-missing-conda.md` (this file)

### Spec Updates

- [x] Update `snakemake-guidelines.md` — add Section 10.5 checkpoint usage rules
- [x] Update `cross-layer-thinking-guide.md` — add script: vs shell: injection distinction
- [ ] Create issue for Snakemake linting rule
- [ ] Create issue for spec sync enforcement

---

## Summary of 6 Bugs Fixed in One Session

### Bug 1: `checkpoint()` called at module level in `score_metrics.smk`

**Error**: `NameError: name 'checkpoint' is not defined`

**Root Cause**: `checkpoint()` is a Snakemake DSL function only available inside rule contexts (input functions, params, shell blocks). It was called at module level (top of `.smk` file), where Snakemake's rule infrastructure is not initialized.

**Fix**: Replace `ancient(checkpoint(...))` with `lambda wildcards: checkpoints.chunk_windows.get(pair_id=wildcards.pair_id).output.manifest` inside the rule's `input:` directive.

### Bug 2: Unresolved `{chunk_id}` wildcard in `score_metrics` output

**Error**: DAG construction failed with unresolved wildcard error.

**Root Cause**: `score_metrics` output pattern included `{chunk_id}` but there was no input file with that wildcard, and no `checkpoint` directive to tell Snakemake how to expand it.

**Fix**: The `checkpoint chunk_windows` rule with `output: manifest = OUTDIR + "/chunks/{pair_id}/chunk_manifest.tsv"` provides the mechanism. Snakemake uses the checkpoint's output to determine `{chunk_id}` values dynamically.

### Bug 3: `ancient()` wrapping `checkpoints.chunk_windows.get()` — invalid Python

**Error**: `NameError: name 'checkpoint' is not defined` (persisted after first fix)

**Root Cause**: `ancient()` evaluates its argument at module level. Even though the value was `checkpoints.chunk_windows.get(...)`, Python tried to evaluate `checkpoints.chunk_windows` at module load time.

**Fix**: Do not use `ancient()` with checkpoint output. The checkpoint itself acts as the dependency marker.

### Bug 4: Missing `conda:` declarations on all rules

**Symptom**: `--conda-create-envs-only` built the DAG successfully but `.snakemake/conda/` remained empty. No conda environments were created.

**Root Cause**: None of the 12 rules in `rules/*.smk` had a `conda:` directive. Snakemake only creates environments for rules that declare them.

**Fix**: Added `conda: "../env/SNAP-HiC-Diff.yaml"` to all rules.

### Bug 5: `extract_expected.smk` EOF / triple-quote nesting error

**Error**: Python `SyntaxError: EOF in multi-line string`

**Root Cause**: The `shell:` block contained a triple-quoted Python string (`'''...'''`) with single quotes inside. Python's tokenizer saw the single quote inside the string as a string delimiter, causing premature EOF.

**Fix**: Extracted the Python script to `scripts/extract_expected.py` and changed the rule to use `script:` directive instead of `shell:`.

### Bug 6: `score_metrics` shell block using `{input.balanced_a}` with lambda inputs

**Symptom**: After converting mcool paths to lambda functions for wildcard derivation, the shell block still referenced `{input.balanced_a}` which is invalid — lambda input functions produce dynamic paths, not named attributes.

**Fix**: Replaced `{input.balanced_a}` etc. with shell variable construction: `SAMPLE_A="${{PAIR_ID%%_vs_*}}"`, then construct full paths in shell. Alternatively, convert to `script:` directive which provides cleaner Snakemake object access.
