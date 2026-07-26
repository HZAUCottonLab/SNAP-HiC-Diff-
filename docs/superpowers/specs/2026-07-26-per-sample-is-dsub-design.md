# Per-Sample CWorld IS dsub Scheduling Design

## Goal

For `N` samples with balanced chromosome-level `.cool` inputs, submit exactly `N` CWorld insulation-score worker jobs through Snakemake's dsub executor. Each worker processes every configured chromosome for one sample and writes the final Meltron-format IS TSV.

The outer `dsub -s run_snakemake.sh` controller remains a separate lightweight dsub task. Thus, once upstream balanced files are available, a 10-sample IS run has one controller plus at most 10 `compute_IS` child jobs.

## Context and root cause

The committed workflow already models `compute_IS` as one checkpoint per `{sample}`. The uncommitted CWorld adaptation changed the checkpoint to `{sample} × {chrom}`, added `compute_IS_bedgraphs` as an output-only rule, and made `merge_IS` a separate rule.

`compute_IS_bedgraphs` has no `shell`, `run`, or `script` directive. Snakemake therefore executes it as a local dependency-tracking rule. This is not where IS computation runs. The observed `38487865.err` records its `MissingOutputException`: the CWorld worker wrote its `.done` marker even after warning that expected bedGraph files were missing.

## Chosen architecture

### One CWorld worker per sample

Keep `checkpoint compute_IS`, but restore its sole wildcard to `{sample}`. Its input waits for `extract_and_balance_sample` for the same sample; its output consists of:

- `results/logs/compute_IS/{sample}.done`
- `results/IS/{sample}/{sample}_IS.tsv`

The job loops over `CHROMS` internally. For each chromosome it reads `results/balanced/{sample}/{sample}__{chrom}_iced.cool`, converts that matrix to CWorld format, and runs `matrix2insulation.pl` at every configured insulation distance. When all chromosome-distance bedGraphs are present, the same dsub job invokes `convert_cworld_IS_to_meltron_IS.py` to write the final TSV, then writes the done marker.

### Remove the local output-only layer

Remove `compute_IS_bedgraphs` and `merge_IS`. The final IS TSV is a declared output of `compute_IS`, so Meltron retains its existing input path without a synthetic local rule or an additional dsub job.

### Fail closed on missing CWorld output

For every chromosome and insulation distance, the worker must verify the expected CWorld bedGraph exists after `matrix2insulation.pl` completes. If it is absent, the shell exits non-zero and does not create either final output. The conversion step must likewise complete before Snakemake can observe the final TSV and `.done` marker.

### Keep cluster responsibility unchanged

`run_snakemake.sh` remains the dsub-submitted Snakemake controller and keeps its cluster-generic executor. `config/dsub_submit.sh` remains responsible for submitting each executable rule as an individual dsub child job. No rule is added to `localrules`.

## Files and responsibilities

| File | Change |
| --- | --- |
| `rules/compute_IS.smk` | Convert the CWorld implementation to one executable checkpoint per sample; inline final TSV conversion; validate every bedGraph; remove local declaration and merge rules. |
| `Snakefile` | Keep the two final per-sample IS targets but ensure they point at `compute_IS` outputs rather than a removed merge stage. |
| `config/config.yaml` | Retain CWorld insulation-distance and per-worker thread configuration; document that `compute_is.threads` is the dsub CPU request for one sample worker. |
| `scripts/convert_cworld_IS_to_meltron_IS.py` | Retain as the conversion boundary; validate its documented CLI contract against the new one-job flow. |
| `scripts/validate_cworld_is.py` | Use as the output-integrity validator for the smoke fixture. |
| `docs/SNAKEMAKE_PIPELINE_MAP.md`, `docs/CODEBASE_OVERVIEW.md` | Update the IS section to describe one dsub job per sample and remove the obsolete local placeholder / merge stage. |

## Acceptance criteria

1. A dry-run targeting 10 explicit IS TSV outputs schedules exactly 10 `compute_IS` jobs (one distinct `{sample}` per target) and no `compute_IS_bedgraphs` or `merge_IS` jobs.
2. Each `compute_IS` job requests `config["compute_is"]["threads"]` CPUs through the existing dsub wrapper; the controller is not counted as an IS worker.
3. A missing CWorld bedGraph exits the per-sample job non-zero and leaves both `results/IS/{sample}/{sample}_IS.tsv` and its `.done` marker absent.
4. A valid single-sample smoke fixture produces all chromosome-distance bedGraphs, one final Meltron-format IS TSV, and one `.done` marker.
5. `snakemake -n --until compute_IS` and the focused smoke validation complete successfully; the existing full-pipeline dry-run remains DAG-valid.
6. Documentation names the execution granularity correctly: `N` samples produce `N` IS worker jobs after upstream balancing.

## Non-goals

- Do not alter dsub queue selection, the outer controller submission mechanism, or the unrelated Hi-C comparison stages.
- Do not repartition the upstream balance checkpoint.
- Do not submit, cancel, or modify the currently running dsub job during this work.
