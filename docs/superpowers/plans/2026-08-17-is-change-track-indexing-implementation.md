# IS-change Track Indexing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove repeated whole-track IS scans so the repaired candidate set can be called efficiently without changing IS-change results.

**Architecture:** Keep `read_is_track` as the validation boundary, then partition each validated track once into a `(chrom, IS_distance)` lookup. The caller uses those small `start`/`IS` tables in its existing context-grid merges. Every pytest and real-data command runs in a Donau `dsub` job, never on the login node.

**Tech Stack:** Python 3.12, pandas, cooler, pytest, Donau `dsub`, existing Snakemake Conda environment.

## Global Constraints

- Retain KNN-plus-Meltron identity repair and KNN `[window_start, window_end)` candidate cores.
- Do not change IS values, masks, context size, scale bands, scoring, thresholds, clustering, ranking, schemas, or atomic writes.
- Index only after `read_is_track` succeeds, preserving its malformed-input failures.
- Absent chromosome/scale combinations must retain the current all-missing merge behaviour.
- Run all tests and real-data verification through `dsub` on compute nodes.
- Start compute-node validation with `cpu=2` on `arm` (10 GB) and record `MEM_MAX_SUM` and CPU use before changing production resources.
- Do not commit generated logs, temporary outputs, or `__pycache__`.

---

## File structure

- Modify: `scripts/call_is_change_loci.py` — create and consume IS lookup tables.
- Modify: `test/is_change/test_call_is_change_loci.py` — test the partitioning contract.
- Create: `test/is_change/run_is_change_tests.dsub` — full compute-node pytest runner.
- Create: `test/is_change/run_is_change_realdata.dsub` — real candidate/caller replay runner.
- Modify: `docs/superpowers/specs/2026-08-16-is-change-candidate-identity-design.md` — correct observed counts.
- Modify: `docs/superpowers/plans/2026-08-16-is-change-candidate-identity-implementation.md` — correct replay expectations.

### Task 1: Index validated IS tracks once

**Files:**
- Modify: `scripts/call_is_change_loci.py:18-180,246-256`
- Modify: `test/is_change/test_call_is_change_loci.py:13-205`
- Create: `test/is_change/run_is_change_tests.dsub`

**Interfaces:**
- Consumes: a validated `pd.DataFrame` with `chrom`, `start`, `IS_distance`, and `IS`.
- Produces: `index_is_track(track) -> dict[tuple[object, int], pd.DataFrame]`, where each value contains only `start` and `IS`; `track_values(index, chrom, scale) -> pd.DataFrame` returns an empty table with those two columns for an absent key.
- `main`, `output_rows`, and `_event_rows_for_region` consume the lookup dictionaries. `context_frame`, scoring, masking, and output schemas retain their current values.

- [ ] **Step 1: Write the failing lookup test and compute-node test runner**

Create `test/is_change/run_is_change_tests.dsub`:

```bash
#!/bin/bash
#DSUB -n is_change_tests
#DSUB -q arm
#DSUB -aa
#DSUB -R 'cpu=2'
#DSUB -o logs/cluster/is_change_tests_%J.out
set -euo pipefail
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="$ROOT/.snakemake/conda/59271960b794774df1c934b1e7affef8_/bin/python"
cd "$ROOT"
"$PYTHON" -m pytest -q test/is_change
```

Make it executable. In `test/is_change/test_call_is_change_loci.py`, import
`index_is_track` and `track_values`, then add:

```python
def test_index_is_track_partitions_chromosome_and_scale():
    track = pd.DataFrame([
        {"chrom": "chr1", "start": 0, "IS_distance": 32000, "IS": 1.0},
        {"chrom": "chr1", "start": 4000, "IS_distance": 32000, "IS": 2.0},
        {"chrom": "chr1", "start": 0, "IS_distance": 60000, "IS": 3.0},
        {"chrom": "chr2", "start": 0, "IS_distance": 32000, "IS": 4.0},
    ])
    index = index_is_track(track)
    assert set(index) == {("chr1", 32000), ("chr1", 60000), ("chr2", 32000)}
    assert index[("chr1", 32000)].to_dict("records") == [
        {"start": 0, "IS": 1.0}, {"start": 4000, "IS": 2.0},
    ]
    missing = track_values(index, "chr2", 60000)
    assert missing.columns.tolist() == ["start", "IS"]
    assert missing.empty
```

- [ ] **Step 2: Verify RED on a compute node**

Run `bash -n test/is_change/run_is_change_tests.dsub`, then submit:

```bash
dsub -s test/is_change/run_is_change_tests.dsub
```

Inspect `djob -ll <jobid>` and `logs/cluster/is_change_tests_<jobid>.out`.
Expected: failure because `index_is_track` is not importable. Do not run pytest locally.

- [ ] **Step 3: Implement the lookup and replace global filters**

Near the output constants in `scripts/call_is_change_loci.py`, add:

```python
TRACK_VALUE_COLUMNS = ["start", "IS"]

def index_is_track(track):
    return {
        (chrom, int(scale)): group[TRACK_VALUE_COLUMNS].copy().reset_index(drop=True)
        for (chrom, scale), group in track.groupby(
            ["chrom", "IS_distance"], sort=False, observed=True
        )
    }

def track_values(track_index, chrom, scale):
    return track_index.get(
        (chrom, scale), pd.DataFrame(columns=TRACK_VALUE_COLUMNS)
    )
```

In `_event_rows_for_region`, replace each full-track filter with:

```python
for track_index, column in ((track_a, "is_a"), (track_b, "is_b")):
    values = track_values(track_index, chrom, scale)
    frame = frame.merge(
        values.rename(columns={"IS": column}), on="start", how="left"
    )
```

In `main`, replace the two raw assignments with:

```python
track_a = index_is_track(read_is_track(args.is_a, settings))
track_b = index_is_track(read_is_track(args.is_b, settings))
```

Do not change any later `output_rows`, schema assertion, atomic-write, or log operation.

- [ ] **Step 4: Verify GREEN on a compute node and commit**

Submit the unchanged dsub runner again:

```bash
dsub -s test/is_change/run_is_change_tests.dsub
```

Expected: `57 passed` with only known third-party warnings. Then commit:

```bash
git add scripts/call_is_change_loci.py test/is_change/test_call_is_change_loci.py \
  test/is_change/run_is_change_tests.dsub
git commit -m "perf: index IS tracks by chromosome and scale"
```

### Task 2: Replay real data and correct the evidence records

**Files:**
- Create: `test/is_change/run_is_change_realdata.dsub`
- Modify: `docs/superpowers/specs/2026-08-16-is-change-candidate-identity-design.md:13-45,61-62`
- Modify: `docs/superpowers/plans/2026-08-16-is-change-candidate-identity-implementation.md:224-285`

**Interfaces:**
- Consumes: Task 1 caller plus main-worktree H1398/L134 KNN, Meltron, IS, and cooler data.
- Produces: only temporary candidate/caller tables and scheduler-log evidence; it never overwrites main-worktree results.
- Documents: 7,329 selected windows and 1,364 candidate regions as actual `H1398_vs_L134` evidence.

- [ ] **Step 1: Add the real-data dsub runner**

Create `test/is_change/run_is_change_realdata.dsub`:

```bash
#!/bin/bash
#DSUB -n is_change_realdata
#DSUB -q arm
#DSUB -aa
#DSUB -R 'cpu=2'
#DSUB -o logs/cluster/is_change_realdata_%J.out
set -euo pipefail
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MAIN_WORKTREE=/share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff
PYTHON="$MAIN_WORKTREE/.snakemake/conda/59271960b794774df1c934b1e7affef8_/bin/python"
OUTDIR="$(mktemp -d "${TMPDIR:-/tmp}/is-change-realdata.XXXXXX")"
cd "$ROOT"
"$PYTHON" scripts/build_is_change_candidates.py \
  --knn-enriched "$MAIN_WORKTREE/results/knn_empirical_pvalue/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv" \
  --meltron-merged "$MAIN_WORKTREE/results/meltron/H1398_vs_L134/H1398_vs_L134_meltron_merged.tsv" \
  --config config/config.yaml --out "$OUTDIR/candidate_regions.tsv" --log "$OUTDIR/candidates.log"
wc -l "$OUTDIR/candidate_regions.tsv"
grep -E 'Matched|Selected|end delta|Wrote' "$OUTDIR/candidates.log"
"$PYTHON" scripts/call_is_change_loci.py \
  --candidate-regions "$OUTDIR/candidate_regions.tsv" \
  --is-a "$MAIN_WORKTREE/results/IS/H1398/H1398_IS.tsv" \
  --is-b "$MAIN_WORKTREE/results/IS/L134/L134_IS.tsv" \
  --cool-a "$MAIN_WORKTREE"/results/balanced/H1398/*_iced.cool \
  --cool-b "$MAIN_WORKTREE"/results/balanced/L134/*_iced.cool \
  --config config/config.yaml \
  --out-primary "$OUTDIR/primary.tsv" --out-preliminary "$OUTDIR/preliminary.tsv" \
  --out-masked "$OUTDIR/masked.bed" --log "$OUTDIR/caller.log"
wc -l "$OUTDIR/primary.tsv" "$OUTDIR/preliminary.tsv" "$OUTDIR/masked.bed"
cat "$OUTDIR/caller.log"
```

Make it executable, check `bash -n`, then submit:

```bash
dsub -s test/is_change/run_is_change_realdata.dsub
```

- [ ] **Step 2: Verify the replay and scheduler allocation**

Inspect the terminal job and log:

```bash
djob -ll <jobid>
sed -n '1,240p' logs/cluster/is_change_realdata_<jobid>.out
```

Expected candidate evidence: 13,370 identity matches, 7,329 selected windows,
endpoint audit `{4000: 13370}`, and 1,365 TSV lines (header plus 1,364 regions).
Record primary, preliminary, masked counts, `MEM_MAX_SUM`, and `CPU_UTIL_AVG`.
Fail on a traceback, `ValueError`, or non-success job status; do not prescribe a biological locus count.

- [ ] **Step 3: Correct stale candidate-count evidence and commit**

In both existing candidate-identity documents, replace every H1398/L134 claim
of 2,622 selected windows or approximately 1,061 regions with 7,329 and
1,364 respectively. Keep identity, endpoint, coordinate, threshold, and
scientific semantics unchanged. Then commit:

```bash
git add test/is_change/run_is_change_realdata.dsub \
  docs/superpowers/specs/2026-08-16-is-change-candidate-identity-design.md \
  docs/superpowers/plans/2026-08-16-is-change-candidate-identity-implementation.md
git commit -m "test: add IS-change compute replay validation"
```

### Task 3: Whole-branch review and integration handoff

**Files:**
- Verify: all commits on `codex/is-change-candidate-identity`
- Verify: terminal dsub logs from Tasks 1 and 2

**Interfaces:**
- Consumes: reviewed task commits and compute-node evidence.
- Produces: a reviewed branch; it does not merge or delete any worktree.

- [ ] **Step 1: Create an independent review package**

```bash
BASE=$(git merge-base master HEAD)
bash /share/home/xhhuang24/.codex/skills/subagent-driven-development/scripts/review-package \
  "$BASE" HEAD .superpowers/sdd/review-track-indexing-final.diff
```

Use a fresh reviewer to check scientific semantics, lookup completeness,
missing-key behaviour, ordering, test coverage, and generated-file scope.

- [ ] **Step 2: Recheck branch and dsub evidence**

```bash
git diff --check "$BASE"..HEAD
git status --short
djob -ll <realdata-jobid>
```

Expected: no whitespace errors; only intended commits; successful compute-node
test and replay jobs. Treat third-party deprecation warnings separately from failures.

- [ ] **Step 3: Offer integration choice**

Report the commit range, dsub evidence, real-data counts, and reviewer result.
Offer merge to `master`, a pull request, or retaining the isolated branch. Do
not merge or delete the worktree without the user's choice.
