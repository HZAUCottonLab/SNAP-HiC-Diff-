# Plan 03-02 Summary: Snakefile Update

**Phase:** 03-knn-significant-heatmap
**Plan:** 02
**Status:** COMPLETE
**Wave:** 1

## Tasks Completed

### Task 1: Update Snakefile

**Files modified:** `Snakefile`

**Changes:**
1. Added `include: "rules/knn_significant_heatmap.smk"` in the include block
2. Added `_all_knn_heatmap_pngs()` helper function (reads manifest TSV, expands all PNG paths)
3. Replaced `rule all` checkpoint manifest target with `_all_knn_heatmap_pngs()` call

**Verification:**
```bash
snakemake -n 2>&1 | grep "plot_knn" | tail -3
# plot_knn_significant_heatmap   1
# plot_knn_significant_window  395
# total  396
```

## Deviations from Plan

1. **`rule all` target** — plan specified checkpoint `.done` file; actual uses Python helper `_all_knn_heatmap_pngs()` to directly enumerate all PNG outputs
2. **`include` placement** — added after `knn_empirical_pvalue.smk` include
