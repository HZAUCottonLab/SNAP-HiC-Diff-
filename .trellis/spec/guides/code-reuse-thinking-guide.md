# Code Reuse Thinking Guide

> **Purpose**: Stop and think before creating new code - does it already exist?

---

## The Problem

**Duplicated code is the #1 source of inconsistency bugs.**

When you copy-paste or rewrite existing logic:
- Bug fixes don't propagate
- Behavior diverges over time
- Codebase becomes harder to understand

---

## Before Writing New Code

### Step 1: Search First

```bash
# Search for similar function names
grep -r "functionName" .

# Search for similar logic
grep -r "keyword" .
```

### Step 2: Ask These Questions

| Question | If Yes... |
|----------|-----------|
| Does a similar function exist? | Use or extend it |
| Is this pattern used elsewhere? | Follow the existing pattern |
| Could this be a shared utility? | Create it in the right place |
| Am I copying code from another file? | **STOP** - extract to shared |

---

## Common Duplication Patterns

### Pattern 1: Copy-Paste Functions

**Bad**: Copying a validation function to another file

**Good**: Extract to shared utilities, import where needed

### Pattern 2: Similar Components

**Bad**: Creating a new component that's 80% similar to existing

**Good**: Extend existing component with props/variants

### Pattern 3: Repeated Constants

**Bad**: Defining the same constant in multiple files

**Good**: Single source of truth, import everywhere

---

## When to Abstract

**Abstract when**:
- Same code appears 3+ times
- Logic is complex enough to have bugs
- Multiple people might need this

**Don't abstract when**:
- Only used once
- Trivial one-liner
- Abstraction would be more complex than duplication

---

## After Batch Modifications

When you've made similar changes to multiple files:

1. **Review**: Did you catch all instances?
2. **Search**: Run grep to find any missed
3. **Consider**: Should this be abstracted?

---

## Gotcha: Asymmetric Mechanisms Producing Same Output

**Problem**: When two different mechanisms must produce the same file set (e.g., recursive directory copy for init vs. manual `files.set()` for update), structural changes (renaming, moving, adding subdirectories) only propagate through the automatic mechanism. The manual one silently drifts.

**Symptom**: Init works perfectly, but update creates files at wrong paths or misses files entirely.

**Prevention checklist**:
- [ ] When migrating directory structures, search for ALL code paths that reference the old structure
- [ ] If one path is auto-derived (glob/copy) and another is manually listed, the manual one needs updating
- [ ] Add a regression test that compares outputs from both mechanisms

---

## Gotcha: Copy-Pasted Merge Patterns Lose Sync

**Problem**: Three scripts (`merge`, `gauss`, `tw`) all needed the same data-join operation: covariates TSV + per-window metrics → merged TSV with standard column order. The merge logic was copy-pasted into each script, not shared. Fixing a bug in one didn't propagate to the others.

**Symptom**: Bug fixed in one file recurs identically in the other two.

**Root cause categories**:
- **C. Change Propagation Failure**: A fix to `merge_window_pair_metrics.py`'s merge was not applied to `compare_gaussian_smoothed_windows.py` and `compare_iced_windows_tw.py`.
- **D. Test Coverage Gap**: Each script had unit tests for individual functions, but no integration test calling `main()` end-to-end on real data. The bug only surfaced in integration.

**When to abstract**:
- Same `pd.merge` pattern appears 3+ times
- Same column ordering logic appears 3+ times
- Same metadata injection pattern appears 3+ times

**Correct pattern** (extract to shared utility):

```python
# scripts/utils/df_transforms.py
def merge_metrics_with_covariates(df_cov, df_metrics, comparison, sample_a, sample_b):
    """Standard join: covariates (left) + metrics (right) on window_id."""
    keys = ["chrom", "window_start", "window_end", "window_id"]
    metric_cols = ["window_id"] + [
        c for c in df_metrics.columns
        if c not in ["comparison", "sample_a", "sample_b", "chrom"]
        and c not in keys
    ]
    df_right = df_metrics.set_index(keys)[metric_cols]
    df_out = pd.merge(
        df_cov.set_index(keys),
        df_right,
        left_index=True,
        right_index=True,
        how="inner"
    ).reset_index()
    df_out["comparison"] = comparison
    df_out["sample_a"] = sample_a
    df_out["sample_b"] = sample_b
    return df_out
```

**After any multi-file pattern fix**:
- [ ] Grep for the same pattern in all other files before declaring done
- [ ] Ask: should this become a shared utility?
- [ ] If three or more files share the same logic, extract it

---

## Gotcha: Statistical Test Bounds Must Be Evidence-Based

**Problem**: Unit test assumed TW statistic for random matrix always falls within `[-5, 5]`. Empirical result: `-23.5`. The tight bound was based on intuition, not literature or empirical data.

**Prevention**:
- [ ] For statistical tests, cite expected distribution bounds from literature
- [ ] Use empirical ranges (e.g., `(-100, 100)`) rather than round numbers
- [ ] Test finiteness and range before testing specific values

---

## Checklist Before Commit

- [ ] Searched for existing similar code
- [ ] No copy-pasted logic that should be shared
- [ ] Constants defined in one place
- [ ] Similar patterns follow same structure
