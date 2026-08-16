# IS-change Candidate Identity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let every window that passes the declared KNN-plus-Meltron gates enter IS-change region construction even when Meltron reports its right endpoint one 4 kb bin beyond the original KNN window.

**Architecture:** Candidate membership is established by the stable `chrom + window_id` identity key, not by an interval endpoint representation. KNN `window_start/window_end` remain the canonical BED-style core coordinates for merging and IS-change reporting; Meltron coordinates are audited for a shared start but never replace those KNN coordinates.

**Tech Stack:** Python 3, pandas, pytest, Snakemake, existing SNAP-HiC-Diff Conda environment.

## Global Constraints

- Retain the configured raw KNN P-value gate and `meltron_best_score >= 5.0` gate.
- Retain KNN `[window_start, window_end)` as the reportable candidate core interval.
- Retain `context_bp = 300000`; flanking context supports scoring only and cannot create reportable loci outside the merged core.
- Do not alter IS values, missing-bin masking, scale bands, event thresholds, or TAD/SV interpretation.
- Reject duplicate Meltron identity keys, duplicate KNN identity keys, Meltron keys absent from KNN, and Meltron rows whose start differs from the KNN start.
- Allow and log Meltron--KNN endpoint differences; do not subtract 4 kb from the Meltron source table.

---

## File structure

- Modify: `scripts/is_change_candidates.py` — validate stable keys, perform the identity join, preserve canonical KNN coordinates, and return a compact join audit.
- Modify: `scripts/build_is_change_candidates.py` — log the audit and the number of selected windows before writing merged regions.
- Modify: `test/is_change/test_is_change_candidates.py` — prove the endpoint-offset case, malformed-key failures, and CLI audit logging.
- Create: no production files.

### Task 1: Implement identity-based candidate selection

**Files:**
- Modify: `scripts/is_change_candidates.py:1-27`
- Test: `test/is_change/test_is_change_candidates.py:53-91`

**Interfaces:**
- Consumes: `knn` with `chrom`, `window_id`, `window_start`, `window_end`, and five `*_empirical_pvalue` columns; `meltron` with `chrom`, `window_id`, `start`, `end`, `ID`, and `meltron_best_score`.
- Produces: `select_final_windows(knn, meltron, settings) -> tuple[pd.DataFrame, dict]`, where the DataFrame retains KNN interval columns and the audit contains `matched_windows`, `selected_windows`, and `meltron_end_delta_bp`.

- [ ] **Step 1: Write the failing endpoint-offset test**

  Replace `test_selects_raw_knn_and_meltron_intersection` with a test whose first Meltron row uses the same ID and start as KNN but has a 4 kb larger end:

  ```python
  meltron = pd.DataFrame([
      meltron_row(1, 0, 304000, 5.0),
      meltron_row(2, 300000, 604000, 6.0),
      meltron_row(3, 600000, 904000, 6.0),
  ])

  selected, audit = select_final_windows(knn, meltron, settings)

  assert selected[["window_id", "window_start", "window_end"]].to_dict("records") == [
      {"window_id": 1, "window_start": 0, "window_end": 300000}
  ]
  assert audit == {
      "matched_windows": 3,
      "selected_windows": 1,
      "meltron_end_delta_bp": {4000: 3},
  }
  ```

- [ ] **Step 2: Run the new test and verify RED**

  Run:

  ```bash
  pytest -q test/is_change/test_is_change_candidates.py -k endpoint
  ```

  Expected: failure because the current coordinate join returns no selected row
  and `select_final_windows` returns only one value.

- [ ] **Step 3: Add required-key, uniqueness, and start validation**

  In `scripts/is_change_candidates.py`, replace `KNN_COORDS` with explicit
  identity and coordinate constants and add these helpers:

  ```python
  IDENTITY_KEYS = ["chrom", "window_id"]
  KNN_COORDS = ["chrom", "window_start", "window_end"]

  def _require_columns(frame, required, label):
      missing = set(required).difference(frame.columns)
      if missing:
          raise ValueError(f"{label} missing columns: {', '.join(sorted(missing))}")

  def _require_unique_identity(frame, label):
      if frame.duplicated(IDENTITY_KEYS).any():
          raise ValueError(f"{label} has duplicate chrom/window_id keys")
  ```

  Require `chrom/window_id/window_start/window_end` in KNN and
  `chrom/window_id/start/end/ID/meltron_best_score` in Meltron.  Check both
  identity key sets for duplicates.  Before the final merge, left-join Meltron
  keys to KNN keys with `validate="one_to_one"`; if any Meltron key is absent,
  raise `ValueError("Meltron windows missing KNN chrom/window_id keys")`.

- [ ] **Step 4: Replace the coordinate join with the identity join**

  Keep raw-P calculation unchanged.  Use KNN coordinates as the only output
  geometry and calculate an explicit audit:

  ```python
  merged = knn.merge(
      meltron,
      on=IDENTITY_KEYS,
      how="inner",
      validate="one_to_one",
      suffixes=("", "_meltron"),
  )
  if not merged["window_start"].eq(merged["start"]).all():
      raise ValueError("KNN and Meltron start coordinates disagree for matched windows")

  audit = {
      "matched_windows": len(merged),
      "meltron_end_delta_bp": {
          int(delta): int(count)
          for delta, count in merged["end"].sub(merged["window_end"])
          .value_counts().sort_index().items()
      },
  }
  selected = merged.loc[
      merged["knn_any_raw_pvalue_pass"]
      & (merged["meltron_best_score"] >= settings.meltron_score_threshold)
  ].sort_values(KNN_COORDS, kind="stable").reset_index(drop=True)
  audit["selected_windows"] = len(selected)
  return selected, audit
  ```

- [ ] **Step 5: Run the endpoint test and verify GREEN**

  Run:

  ```bash
  pytest -q test/is_change/test_is_change_candidates.py -k endpoint
  ```

  Expected: `1 passed`.

- [ ] **Step 6: Add failure tests for unsafe identity mappings**

  Add an exact-coordinate compatibility test and three failure tests:

  ```python
  def test_selects_exact_coordinate_identity_match(settings):
      knn = pd.DataFrame([knn_row(1, 0, 300000, MSE_empirical_pvalue=0.01)])
      meltron = pd.DataFrame([meltron_row(1, 0, 300000, 6.0)])
      selected, audit = select_final_windows(knn, meltron, settings)
      assert selected["window_id"].tolist() == [1]
      assert audit["meltron_end_delta_bp"] == {0: 1}

  def test_rejects_duplicate_knn_identity(settings):
      knn = pd.DataFrame([knn_row(1, 0, 300000), knn_row(1, 0, 300000)])
      meltron = pd.DataFrame([meltron_row(1, 0, 304000, 6.0)])
      with pytest.raises(ValueError, match="KNN has duplicate chrom/window_id keys"):
          select_final_windows(knn, meltron, settings)

  def test_rejects_meltron_identity_absent_from_knn(settings):
      knn = pd.DataFrame([knn_row(1, 0, 300000)])
      meltron = pd.DataFrame([meltron_row(2, 0, 304000, 6.0)])
      with pytest.raises(ValueError, match="Meltron windows missing KNN"):
          select_final_windows(knn, meltron, settings)

  def test_rejects_shifted_meltron_start(settings):
      knn = pd.DataFrame([knn_row(1, 0, 300000)])
      meltron = pd.DataFrame([meltron_row(1, 4000, 304000, 6.0)])
      with pytest.raises(ValueError, match="start coordinates disagree"):
          select_final_windows(knn, meltron, settings)
  ```

- [ ] **Step 7: Run the full candidate unit-test module**

  Run:

  ```bash
  pytest -q test/is_change/test_is_change_candidates.py
  ```

  Expected: all tests pass, including existing merge and atomic-output tests.

- [ ] **Step 8: Commit Task 1**

  ```bash
  git add scripts/is_change_candidates.py test/is_change/test_is_change_candidates.py
  git commit -m "fix: join IS-change candidates by window identity"
  ```

### Task 2: Log the candidate-interface audit

**Files:**
- Modify: `scripts/build_is_change_candidates.py:39-48`
- Test: `test/is_change/test_is_change_candidates.py:169-220`

**Interfaces:**
- Consumes: `(selected, audit)` from `select_final_windows`.
- Produces: candidate-build log lines that record matched rows, selected rows,
  and the Meltron-minus-KNN endpoint-delta distribution.

- [ ] **Step 1: Write the failing CLI audit assertion**

  In `test_cli_writes_sorted_region_schema`, change the fixture Meltron ends
  from 300 kb to 304 kb, then add:

  ```python
  log_text = log.read_text()
  assert "Matched 3 KNN/Meltron windows by chrom/window_id" in log_text
  assert "Selected 3 windows after KNN and Meltron thresholds" in log_text
  assert "Meltron minus KNN end delta(s) in bp: {4000: 3}" in log_text
  ```

- [ ] **Step 2: Run the CLI audit test and verify RED**

  Run:

  ```bash
  pytest -q test/is_change/test_is_change_candidates.py -k cli_writes_sorted
  ```

  Expected: failure because the CLI still expects a single return value and
  does not write the audit lines.

- [ ] **Step 3: Consume the audit and write deterministic log messages**

  In `scripts/build_is_change_candidates.py`, replace the one-value assignment
  and final log with:

  ```python
  selected, audit = select_final_windows(
      pd.read_csv(args.knn_enriched, sep="\t"),
      pd.read_csv(args.meltron_merged, sep="\t"),
      settings,
  )
  log(
      f"Matched {audit['matched_windows']} KNN/Meltron windows by chrom/window_id",
      args.log,
  )
  log(
      f"Selected {audit['selected_windows']} windows after KNN and Meltron thresholds",
      args.log,
  )
  log(f"Meltron minus KNN end delta(s) in bp: {audit['meltron_end_delta_bp']}", args.log)
  ```

  Keep the existing atomic output write and final `Wrote ... candidate regions`
  log unchanged.

- [ ] **Step 4: Run the CLI audit test and verify GREEN**

  Run:

  ```bash
  pytest -q test/is_change/test_is_change_candidates.py -k cli_writes_sorted
  ```

  Expected: `1 passed`.

- [ ] **Step 5: Run the complete IS-change test suite**

  Run:

  ```bash
  pytest -q test/is_change
  ```

  Expected: all tests pass; any pre-existing third-party Matplotlib warnings
  are reported separately from test failures.

- [ ] **Step 6: Commit Task 2**

  ```bash
  git add scripts/build_is_change_candidates.py test/is_change/test_is_change_candidates.py
  git commit -m "feat: audit IS-change candidate coordinate conventions"
  ```

### Task 3: Verify the repaired interface on `H1398_vs_L134`

**Files:**
- Modify: no source files.
- Verify: real KNN, Meltron, IS, and cooler inputs already present under the
  main worktree `results/` directory.

**Interfaces:**
- Consumes: branch versions of the candidate builder and caller; existing
  `H1398_vs_L134` data.
- Produces: a temporary candidate table, temporary IS-change outputs, and an
  evidence-backed count of selected windows, merged regions, and called loci.

- [ ] **Step 1: Build candidates with real inputs into a temporary directory**

  From the isolated worktree, run the branch script against absolute paths in
  the main worktree, writing only to a `mktemp -d` directory:

  ```bash
  MAIN_WORKTREE=/share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff
  CONDA_PREFIX_DIR=/share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff/.snakemake/conda/59271960b794774df1c934b1e7affef8_
  candidate_verify_dir=$(mktemp -d /tmp/is-change-candidates.XXXXXX)
  conda run -p "$CONDA_PREFIX_DIR" python scripts/build_is_change_candidates.py \
    --knn-enriched "$MAIN_WORKTREE/results/knn_empirical_pvalue/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv" \
    --meltron-merged "$MAIN_WORKTREE/results/meltron/H1398_vs_L134/H1398_vs_L134_meltron_merged.tsv" \
    --config config/config.yaml \
    --out "$candidate_verify_dir/candidate_regions.tsv" \
    --log "$candidate_verify_dir/candidates.log"
  ```

- [ ] **Step 2: Verify the candidate counts and provenance**

  Run:

  ```bash
  wc -l "$candidate_verify_dir/candidate_regions.tsv"
  grep -E 'Matched|Selected|end delta|Wrote' "$candidate_verify_dir/candidates.log"
  ```

  Expected: 2,622 selected windows, `{4000: 13370}` endpoint audit, and 1,061
  merged candidate rows plus the TSV header.

- [ ] **Step 3: Replay IS-change calling with real inputs to temporary outputs**

  Run the branch `scripts/call_is_change_loci.py` with the temporary candidate
  table, the two absolute IS-track paths, all absolute H1398 and L134 cooler
  shards, and temporary `primary.tsv`, `preliminary.tsv`, `masked.bed`, and
  `caller.log` paths.  Preserve the current configuration and do not write into
  the main `results/is_change/` directory.


  ```bash
  conda run -p "$CONDA_PREFIX_DIR" python scripts/call_is_change_loci.py \
    --candidate-regions "$candidate_verify_dir/candidate_regions.tsv" \
    --is-a "$MAIN_WORKTREE/results/IS/H1398/H1398_IS.tsv" \
    --is-b "$MAIN_WORKTREE/results/IS/L134/L134_IS.tsv" \
    --cool-a "$MAIN_WORKTREE"/results/balanced/H1398/*_iced.cool \
    --cool-b "$MAIN_WORKTREE"/results/balanced/L134/*_iced.cool \
    --config config/config.yaml \
    --out-primary "$candidate_verify_dir/primary.tsv" \
    --out-preliminary "$candidate_verify_dir/preliminary.tsv" \
    --out-masked "$candidate_verify_dir/masked.bed" \
    --log "$candidate_verify_dir/caller.log"
  ```

- [ ] **Step 4: Inspect data-derived output counts and failure logs**

  Run:

  ```bash
  wc -l "$candidate_verify_dir/primary.tsv" \
        "$candidate_verify_dir/preliminary.tsv" \
        "$candidate_verify_dir/masked.bed"
  grep -n -E 'Traceback|ValueError|Wrote' "$candidate_verify_dir/caller.log"
  ```

  Expected: no traceback or `ValueError`; report the observed primary,
  preliminary, and masked counts without precommitting to a biological result.

- [ ] **Step 5: Commit verification-only metadata only if source files changed**

  Do not commit temporary outputs.  If Tasks 1 and 2 already committed all
  source and test changes, leave this task without an additional commit.
