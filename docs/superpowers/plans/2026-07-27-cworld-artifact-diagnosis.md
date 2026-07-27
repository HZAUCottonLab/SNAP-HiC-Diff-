# CWorld Artifact Diagnosis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Retain one-chromosome CWorld smoke artifacts and expose Perl diagnostics needed to identify an absent bedGraph.

**Architecture:** The smoke harness owns the shared run directory. The validator records CWorld stdout, stderr, and prefix artifacts, while keeping its current fail-closed bedGraph check. Production Snakemake scheduling is unchanged.

**Tech Stack:** Bash, Python, pytest, CWorld Dekker Perl, dsub.

## Global Constraints

- Use L021 × HC04_D08 × 32000 bp only.
- Write retained artifacts below `$PROJ/tmp/cworld_smoke/`.
- Request four dsub CPUs and no explicit memory limit.

---

### Task 1: Retained diagnostic boundary

**Files:**
- Modify: `test/cworld_is/test_compute_is_oom_script.py`
- Modify: `test/compute_is_oom_20260609/run_test.sh`
- Modify: `scripts/validate_cworld_is.py`

- [ ] Write a failing test asserting the shared run directory, `-v`,
  `cworld_stdout.log`, `cworld_stderr.log`, and `cworld_artifacts.txt`.
- [ ] Run the focused pytest test and confirm it fails.
- [ ] Add the minimum harness and validator behavior.
- [ ] Run syntax plus focused smoke/validator tests and confirm they pass.
- [ ] Commit only the diagnostic implementation files.

### Task 2: Retained single-chromosome smoke

**Files:** no source changes.

- [ ] Submit the existing L021 smoke with code, data, and Conda roots explicit.
- [ ] Inspect retained matrix header, stdout, stderr, prefix artifact list,
  bedGraph, and final Meltron TSV.
- [ ] Report the exact first failing boundary before any production-rule fix.
