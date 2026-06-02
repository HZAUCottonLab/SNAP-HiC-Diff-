# Plan: 01-knn-empirical-pvalue

**Wave 1 — Foundation (fully parallel)**

| # | Plan | Objective | Type | Depends |
|---|------|----------|------|--------|
| 01-knn-config | KNN Config | Add `knn_empirical_pvalue` section to `config/config.yaml` | execute | - |
| 01-knn-script | KNN Script | Create `scripts/compute_knn_empirical_pvalue.py` | execute | - |
| 01-knn-rule | KNN Rule | Add `rule knn_empirical_pvalue` to `rules/window_compare.smk` | execute | 01-knn-config, 01-knn-script |
| 01-knn-test | KNN Dry-run | Dry-run Snakemake to verify rule graph | execute | 01-knn-rule |

**Wave 2 — Integration (after Wave 1)**

| # | Plan | Objective | Type | Depends |
|---|------|----------|------|--------|
| 01-knn-all-rule | Update rule all | Add `{comparison}_knn_enriched.tsv` to `rule all` input | execute | 01-knn-rule |

---

## Plan: 01-knn-config

```yaml
wave: 1
depends_on: []
files_modified:
  - config/config.yaml
```

<read_first>
- config/config.yaml
</read_first>

<acceptance_criteria>
- `config/config.yaml` contains a top-level key `knn_empirical_pvalue`
- The section contains all keys: `background_file`, `knn_k`, `knn_metric`, `scaler`, `score_columns`, `derive_scores`, `pvalue_direction`, `compute_fdr`
- All 5 score columns are listed: `MSE`, `d_SSIM`, `d_Spearman`, `tw_score`, `d_NCC`
- `derive_scores` maps `d_SSIM`, `d_Spearman`, `tw_score`, `d_NCC` to their formulas
- `pvalue_direction` for all 5 scores is `"greater"`
- `background_file` points to `/share/home/xhhuang24/project/PanTE/background_pool/results/window_compare/HiC_window_background.tsv`
</acceptance_criteria>

<action>
Add the following YAML block to `config/config.yaml` under the existing `window_compare:` section (after `tw: n_jobs: 80`):

```yaml
# KNN Empirical P-Value Parameters
knn_empirical_pvalue:
  background_file: "/share/home/xhhuang24/project/PanTE/background_pool/results/window_compare/HiC_window_background.tsv"
  knn_k: 10000
  knn_metric: "euclidean"
  scaler: "robust"
  score_columns:
    - MSE
    - d_SSIM
    - d_Spearman
    - tw_score
    - d_NCC
  derive_scores:
    d_SSIM: "1 - SSIM"
    d_Spearman: "1 - Spearman"
    tw_score: "-log10(max(tw_pvalue, 1e-30))"
    d_NCC: "1 - NCC"
  pvalue_direction:
    MSE: "greater"
    d_SSIM: "greater"
    d_Spearman: "greater"
    tw_score: "greater"
    d_NCC: "greater"
  compute_fdr: false
```

Insert as a sibling of `window_compare:`, not inside it — both are top-level keys under the root config dict.
</action>

---

## Plan: 01-knn-script

```yaml
wave: 1
depends_on: []
files_modified:
  - scripts/compute_knn_empirical_pvalue.py
```

<read_first>
- scripts/compute_window_pair_covariates.py
- scripts/merge_window_pair_metrics.py
- .cursor/rules/010-snakemake-rule-format.mdc
</read_first>

<acceptance_criteria>
- `scripts/compute_knn_empirical_pvalue.py` exists and is executable
- `python scripts/compute_knn_empirical_pvalue.py --help` outputs usage with all required arguments
- Script accepts: `--comparison`, `--input`, `--background`, `--out`, `--log`
- Config is loaded internally via `yaml.safe_load("config/config.yaml")["knn_empirical_pvalue"]` — no `--config` CLI flag
- Background load: reads `EC`, `IB`, `SP`, `MSE`, `SSIM`, `Spearman`, `tw_pvalue`, `NCC` from TSV, drops rows with NA/inf/-inf
- Scaler: `robust` (median + IQR) fitted only on background rows for `EC`, `IB`, `SP` and each score column
- BallTree: built from background rows, standardized using fitted scaler, queried with k=10000 per comparison window
- Per-window query: batched with `joblib.Parallel` (n_jobs from threads param) using `delayed(tree.query)` to parallelize across all windows in a comparison — NOT a sequential for-loop
- Empirical P-value: right-tail `(count(bg_value >= query_value) + 1) / (knn_k + 1)` for each score column
- Output: `chrom`, `window_id`, `window_start`, `window_end`, `comparison`, then for each of 5 scores: `{score}_empirical_pvalue`, `bg_{score}_median`, then `knn_valid_k`, `knn_dist_median`, `knn_dist_max`, `match_quality_flag`
- Log: prints background row count, EC/IB/SP median and IQR, BallTree parameters, K value, P-value directions, output path
</acceptance_criteria>

<action>
Create `scripts/compute_knn_empirical_pvalue.py` with the following structure:

```
#!/usr/bin/env python3
"""
KNN-based empirical P-value computation for window pair metrics.

Usage:
    python scripts/compute_knn_empirical_pvalue.py \
        --comparison H1398_vs_L134 \
        --input results/window_compare/H1398_vs_L134/H1398_vs_L134_window_pair_metrics.tsv \
        --background /share/home/xhhuang24/project/PanTE/background_pool/results/window_compare/HiC_window_background.tsv \
        --out results/window_compare/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv \
        --log logs/window_compare/H1398_vs_L134/knn_empirical_pvalue.log
"""

import argparse
import yaml
import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree

def derive_scores(df, config):
    """Derive score columns from raw metrics."""
    df = df.copy()
    derive = config.get("derive_scores", {})
    for col, formula in derive.items():
        if "max(tw_pvalue" in formula:
            df[col] = -np.log10(np.maximum(df["tw_pvalue"], 1e-30))
        elif "1 - " in formula:
            src = formula.split("1 - ")[1]
            df[col] = 1.0 - df[src]
    return df

def robust_scaler_fit(data):
    """Fit robust scaler: median + IQR."""
    med = np.nanmedian(data)
    q75, q25 = np.nanpercentile(data, [75, 25])
    iqr = q75 - q25
    return med, iqr

def robust_scaler_transform(data, med, iqr):
    """Transform using fitted scaler."""
    if iqr == 0:
        return np.zeros_like(data, dtype=float)
    return (data - med) / iqr

def compute_knn_pvalues(query_df, bg_df, knn_k, n_jobs, log_file):
    """Compute KNN-based empirical P-values."""
    # Fit scaler on background for EC, IB, SP
    ec_med, ec_iqr = robust_scaler_fit(bg_df["EC"].values)
    ib_med, ib_iqr = robust_scaler_fit(bg_df["IB"].values)
    sp_med, sp_iqr = robust_scaler_fit(bg_df["SP"].values)

    # Build BallTree from background
    bg_knn_features = np.column_stack([
        robust_scaler_transform(bg_df["EC"].values, ec_med, ec_iqr),
        robust_scaler_transform(bg_df["IB"].values, ib_med, ib_iqr),
        robust_scaler_transform(bg_df["SP"].values, sp_med, sp_iqr),
    ])
    tree = BallTree(bg_knn_features, metric="euclidean")

    # Fit scalers for each score column on background
    score_cols = ["MSE", "d_SSIM", "d_Spearman", "tw_score", "d_NCC"]
    bg_medians = {}
    bg_iqrs = {}
    for col in score_cols:
        m, i = robust_scaler_fit(bg_df[col].values)
        bg_medians[col] = m
        bg_iqrs[col] = i

    # Pre-compute scaled query features for all windows at once
    query_ec = robust_scaler_transform(query_df["EC"].values, ec_med, ec_iqr)
    query_ib = robust_scaler_transform(query_df["IB"].values, ib_med, ib_iqr)
    query_sp = robust_scaler_transform(query_df["SP"].values, sp_med, sp_iqr)
    query_features = np.column_stack([query_ec, query_ib, query_sp])

    # Batch KNN query with joblib
    n_windows = len(query_df)
    log(f"[INFO] Running {n_windows} KNN queries with {n_jobs} workers...", log_file)
    all_dists, all_indices = tree.query(query_features, k=knn_k)

    # Compute P-values per window (vectorized where possible, fallback to loop for edge cases)
    results = []
    for i in range(n_windows):
        dists_i = all_dists[i]
        idx_i = all_indices[i]
        valid_k = len(idx_i)

        row_result = {
            "knn_valid_k": valid_k,
            "knn_dist_median": float(np.median(dists_i)) if valid_k > 0 else float("nan"),
            "knn_dist_max": float(dists_i[-1]) if valid_k > 0 else float("nan"),
        }
        for col in score_cols:
            bg_vals = bg_df[col].values[idx_i]
            q_val = query_df[col].values[i]
            count = int(np.sum(bg_vals >= q_val))
            pval = (count + 1) / (valid_k + 1)
            row_result[f"{col}_empirical_pvalue"] = pval
            row_result[f"bg_{col}_median"] = float(np.nanmedian(bg_vals))

        results.append(row_result)

    return pd.DataFrame(results)

def main():
    # Load config internally (script always reads config/config.yaml relative to cwd)
    with open("config/config.yaml") as f:
        knn_cfg = yaml.safe_load(f)["knn_empirical_pvalue"]
    # ... argparse + logic ...
    # Load input, derive scores, load background
    # Drop NA/inf/-inf rows
    # Compute KNN P-values
    # Select only necessary columns for output (not full merge)
    # Write output
```

See acceptance_criteria for exact requirements on each component.
</action>

---

## Plan: 01-knn-rule

```yaml
wave: 1
depends_on:
  - 01-knn-config
  - 01-knn-script
files_modified:
  - rules/window_compare.smk
```

<read_first>
- rules/window_compare.smk
- scripts/compute_knn_empirical_pvalue.py
- config/config.yaml
- .cursor/rules/011-snakemake-wildcards.mdc
</read_first>

<acceptance_criteria>
- `rules/window_compare.smk` contains `rule knn_empirical_pvalue`
- Rule input: `OUTDIR + "/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv"` (no lambda/expand needed — single file per comparison)
- Rule output columns: `chrom`, `window_id`, `window_start`, `window_end`, `comparison`, `MSE_empirical_pvalue`, `bg_MSE_median`, `d_SSIM_empirical_pvalue`, `bg_d_SSIM_median`, `d_Spearman_empirical_pvalue`, `bg_d_Spearman_median`, `tw_score_empirical_pvalue`, `bg_tw_score_median`, `d_NCC_empirical_pvalue`, `bg_d_NCC_median`, `knn_valid_k`, `knn_dist_median`, `knn_dist_max`, `match_quality_flag`
- Rule log: `OUTDIR + "/logs/window_compare/{comparison}/knn_empirical_pvalue.log"`
- Rule benchmark: `OUTDIR + "/benchmark/window_compare/{comparison}/knn_empirical_pvalue.tsv"`
- Rule threads: `config["parallel"]["threads_per_job"]`
- Rule params: reads `knn_k`, `knn_metric`, `scaler`, `score_columns`, `derive_scores`, `pvalue_direction`, `compute_fdr`, `background_file` from `config["knn_empirical_pvalue"]`
- Shell block uses `{wildcards.comparison}` extraction
- Rule fields in correct order: input → output → log → benchmark → threads → params → conda → shell
- Rule placed after `rule merge_window_pair_metrics` in the file
</acceptance_criteria>

<action>
Append the following rule to `rules/window_compare.smk` (after the closing of `rule merge_all_window_pair_metrics`):

```python
# ==============================================================================
# Knn_empirical_pvalue
# ==============================================================================

rule knn_empirical_pvalue:
    """Compute KNN-based empirical P-values for window pair metrics."""
    input:
        script = SCRIPT_DIR + "/scripts/compute_knn_empirical_pvalue.py",
        merged = OUTDIR + "/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv",
    output:
        enriched = OUTDIR + "/window_compare/{comparison}/{comparison}_knn_enriched.tsv",
    log:
        OUTDIR + "/logs/window_compare/{comparison}/knn_empirical_pvalue.log",
    benchmark:
        OUTDIR + "/benchmark/window_compare/{comparison}/knn_empirical_pvalue.tsv",
    threads: config["parallel"]["threads_per_job"]
    params:
        background_file = config["knn_empirical_pvalue"]["background_file"],
        knn_k = config["knn_empirical_pvalue"]["knn_k"],
        knn_metric = config["knn_empirical_pvalue"]["knn_metric"],
        scaler = config["knn_empirical_pvalue"]["scaler"],
        score_columns = config["knn_empirical_pvalue"]["score_columns"],
        derive_scores = config["knn_empirical_pvalue"]["derive_scores"],
        pvalue_direction = config["knn_empirical_pvalue"]["pvalue_direction"],
        compute_fdr = config["knn_empirical_pvalue"]["compute_fdr"],
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        COMPARISON="{wildcards.comparison}"

        mkdir -p "$(dirname "{output.enriched}")" "$(dirname "{log}")"

        python {input.script} \
            --comparison "$COMPARISON" \
            --input "{input.merged}" \
            --background "{params.background_file}" \
            --out "{output.enriched}" \
            --log "{log}" \
            > "{log}" 2>&1
        """
```

Note: the Python script receives all knn_empirical_pvalue params by loading config.yaml directly
(the config file path is implicit: always `config/config.yaml` relative to Snakefile cwd).
No `--config` CLI flag needed — the script uses yaml.safe_load() on `config/config.yaml` internally.
The `--background` flag is passed the hardcoded absolute path from params.
</action>

---

## Plan: 01-knn-test

```yaml
wave: 1
depends_on:
  - 01-knn-rule
files_modified: []
```

<read_first>
- Snakefile
- rules/window_compare.smk
</read_first>

<acceptance_criteria>
- `snakemake --dryrun` exits 0 and shows `knn_empirical_pvalue` jobs in the DAG
- At least one comparison (e.g. `H1398_vs_L134`) shows as a pending job for `knn_empirical_pvalue`
- No missing input errors for the `knn_empirical_pvalue` rule
</acceptance_criteria>

<action>
Run `cd /share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff && snakemake --dryrun 2>&1 | grep -E "knn_empirical_pvalue|missing"` to verify the rule is recognized and no inputs are missing.
</action>

---

## Plan: 01-knn-all-rule

```yaml
wave: 2
depends_on:
  - 01-knn-rule
files_modified:
  - Snakefile
```

<read_first>
- Snakefile
- rules/window_compare.smk
</read_first>

<acceptance_criteria>
- `Snakefile` `rule all` input includes `expand(OUTDIR + "/window_compare/{comparison}/{comparison}_knn_enriched.tsv", comparison=get_window_compare_pairs())`
- The new input is added after the existing `window_compare` entry (line ~123)
- `snakemake --dryrun` still exits 0 after the change
</acceptance_criteria>

<action>
In `Snakefile`, after the existing `window_compare` entry in `rule all`:

```python
        # 7. Window compare (per comparison pair)
        expand(
            OUTDIR + "/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv",
            comparison=get_window_compare_pairs()
        ),
        # 8. KNN empirical P-value (per comparison pair)
        expand(
            OUTDIR + "/window_compare/{comparison}/{comparison}_knn_enriched.tsv",
            comparison=get_window_compare_pairs()
        ),
```

Update the comment count from `# 7.` to `# 7.` and add `# 8.`.
</action>

---

## Verification

<must_haves>

- [ ] `config/config.yaml` contains `knn_empirical_pvalue` section with all 5 score columns
- [ ] `scripts/compute_knn_empirical_pvalue.py` exists and `python --help` works
- [ ] `rules/window_compare.smk` contains `rule knn_empirical_pvalue` with correct field order
- [ ] `Snakefile` `rule all` includes `{comparison}_knn_enriched.tsv`
- [ ] `snakemake --dryrun` exits 0 with no missing input errors for `knn_empirical_pvalue`

</must_haves>
