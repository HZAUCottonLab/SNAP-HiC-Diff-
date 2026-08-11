# Meltron 4 kb per-chromosome CWorld insulation design

> **Historical design (superseded for resolution portability).** This document
> records the original 4 kb implementation design and its H1398/L134 4 kb
> execution context. The current CWorld-to-Meltron resolution contract is the
> [resolution-portable design](2026-08-10-meltron-resolution-portability-design.md):
> any integer source resolution of at least 1,000 bp is supported when the
> configured scale and IDS constraints are satisfied.

## Purpose

This historical design generated the insulation-score (IS) input required by
Meltron using CWorld at a single source resolution of **4,000 bp**, while
making the expensive work parallel at the `sample × chromosome` level. One
scheduler task computed all ten requested IS scales for one chromosome of one
sample; a lightweight task then merged every chromosome's results into the
existing per-sample Meltron input table.

That historical implementation preserved Meltron's expected input interface
(`chrom`, `viewpoint`, `IS_distance`, `IS`) and intentionally prevented CWorld
runs at the other then-available mcool resolutions (20 kb and 100 kb). It does
not limit the current portable contract linked above.

## Scope and non-goals

In scope:

- CWorld insulation-score generation for Meltron.
- The Snakemake topology, scheduler resources, intermediate manifests, merge
  validation, cleanup helper, tests, and README documentation.

Out of scope:

- Changing the downstream Meltron comparison interface.
- Calculating CWorld IS from 20 kb or 100 kb matrices.
- Changing the chosen historical scale series unless the configuration is
  deliberately edited and passes validation.

## Scientific and technical constraints

Meltron consumes long-format IS tables and combines values across IS scales.
Its upstream helper defaults to ten scales and requires scales to be compatible
with the matrix resolution.  The project will retain its historical ten-scale
series, adapted exactly to the 4 kb resolution rather than using Meltron's
default 100 kb--1 Mb series.  See [Meltron's IS conversion
helper](https://github.com/pombo-lab/Meltron/blob/main/command_line_apps/long_matrix_to_IS.R)
and [Meltron's command-line application](https://github.com/pombo-lab/Meltron/blob/main/command_line_apps/Meltron.R).

The required scales are:

| Index | IS distance (bp) | IS distance (4 kb bins) | CWorld `--ids` (bp) |
| --- | ---: | ---: | ---: |
| 1 | 32,000 | 8 | 16,000 |
| 2 | 60,000 | 15 | 30,000 |
| 3 | 92,000 | 23 | 46,000 |
| 4 | 120,000 | 30 | 60,000 |
| 5 | 152,000 | 38 | 76,000 |
| 6 | 180,000 | 45 | 90,000 |
| 7 | 212,000 | 53 | 106,000 |
| 8 | 240,000 | 60 | 120,000 |
| 9 | 272,000 | 68 | 136,000 |
| 10 | 300,000 | 75 | 150,000 |

For every requested scale, the following invariants are mandatory:

1. the source `.cool` has bin size exactly 4,000 bp;
2. the IS distance is positive and divisible by 4,000 bp;
3. `--ids` equals `--is / 2`;
4. there are exactly ten distinct scales;
5. the requested scale recorded in the manifest and final TSV is the source of
   truth.  CWorld's output filename may display an apparent `+1` value (for
   example `is32001`) and must not redefine the scale.

## Configuration contract

All resolution, scale, and scheduler choices belong in Snakemake configuration,
not hard-coded in a shell command or rule body.  The intended configuration
shape is:

```yaml
meltron:
  cworld_is:
    source_resolution_bp: 4000
    is_distances_bp: [32000, 60000, 92000, 120000, 152000,
                      180000, 212000, 240000, 272000, 300000]
    ids_divisor: 2
    chrom_job_threads: 4
    merge_job_threads: 4
```

`chrom_job_threads: 4` is the default resource request for every `compute_IS`
chromosome job.  On this dsub/Donau installation CPU count determines the
memory allocation, so the workflow must request CPUs and must not add a
separate memory setting.  The value remains configurable for later tuning.

## Workflow topology

For a sample with 26 target chromosomes, the workflow must form this graph:

```text
sample × chromosome (26 jobs per sample)
  compute_IS_chrom: validate 4 kb .cool → calculate 10 sequential scales
                                ↓
                    per-chromosome manifest + 10 bedGraphs
                                ↓
all chromosome manifests for one sample
  merge_IS: validate complete matrix → <sample>_IS.tsv
                                ↓
Meltron comparison rules (unchanged input interface)
```

`compute_IS_chrom` is the dsub-submitted, resource-intensive job.  It runs one
chromosome's ten scales sequentially in the same job, so it is not one job per
scale.  `merge_IS` is a separate, lightweight per-sample task.  This gives
`26 × sample_count` CWorld tasks, plus one merge task per sample; for ten
samples, 260 CWorld tasks and 10 merge tasks are expected.

Failure isolation is required: one failed chromosome prevents only that
sample's merge and its downstream Meltron comparisons.  Independent
chromosomes and samples remain schedulable.

## Inputs and outputs

Each chromosome task consumes only its 4 kb chromosome-specific balanced cool
input, conceptually:

```text
results/balanced/<sample>/<sample>__<chromosome>_iced.cool
```

The task must inspect the actual cooler metadata before CWorld conversion and
fail before calculation if its bin size differs from
`source_resolution_bp`.

Intermediate CWorld products stay under:

```text
results/cworld_insulation/<sample>/
```

For each chromosome, a manifest maps every requested IS distance to precisely
one newly generated CWorld bedGraph.  The rule must reject stale, missing, or
ambiguous files instead of accepting a filename glob blindly.

The per-sample merge creates exactly the existing downstream contract:

```text
results/IS/<sample>/<sample>_IS.tsv
```

That TSV has four tab-separated columns, in this order:

```text
chrom  viewpoint  IS_distance  IS
```

There is no requirement to create a separate final IS TSV per chromosome; the
per-chromosome manifests and bedGraphs are the parallel-work products, and the
merged TSV is the Meltron input.

## Safe cleanup and restart

The restart helper accepts explicit sample names only.  It must default to a
dry run and require an explicit `--execute` opt-in before removing anything.
For every selected sample, the cleanup scope is limited to these four product
classes:

1. `results/cworld_insulation/<sample>/`
2. `results/IS/<sample>/`
3. `results/logs/compute_IS/<sample>*`
4. `results/benchmark/compute_IS/<sample>*`

It must never delete broad repository paths, unrelated samples, source inputs,
or Conda environments.  After an executed cleanup, the user can resubmit the
full workflow with `dsub -s run_snakemake.sh`; no existing controller job is
cancelled by the helper.

## Historical validation and acceptance criteria

These criteria governed the original fixed-4 kb implementation and are
preserved as an execution record. The current requirements are in the
[resolution-portable design](2026-08-10-meltron-resolution-portability-design.md).

That historical implementation was acceptable only when all of the following
were demonstrated:

1. Configuration validation rejected a non-4 kb source resolution, a scale not
   divisible by 4 kb, a wrong scale count, duplicate scales, or an inconsistent
   `ids_divisor`.
2. Automated tests assert the exact ten values and their expected `--ids`
   values listed above.
3. A Snakemake dry run for one sample exposes 26 `compute_IS_chrom` jobs and
   one `merge_IS` job (or the equivalent count for the configured chromosome
   set), not one sample-wide CWorld job and not jobs for other mcool
   resolutions.
4. A chromosome task validates its input cooler is 4 kb, generates exactly ten
   fresh bedGraphs, and writes a complete manifest mapping the requested
   distances to those files.
5. The merge task requires every configured chromosome × every one of the ten
   scales before creating output.  It emits the required four-column table and
   retains requested distances such as `32000`, not CWorld filename artifacts
   such as `32001`.
6. A deliberately failed chromosome blocks only its own sample merge and
   downstream Meltron calculation.
7. The cleanup helper is tested in both dry-run and `--execute` modes against
   explicit sample fixtures.
8. Existing Meltron rules continue to consume
   `results/IS/<sample>/<sample>_IS.tsv` without an interface change.

## Historical README requirements

The original README update described the fixed-4 kb execution configuration,
not the current portable contract. It covered:

- the scientific purpose: CWorld IS is produced solely as Meltron input;
- the fixed 4 kb source-resolution configuration and the exact ten historical
  scales;
- the `--ids = --is / 2` rule and the fact that CWorld output names may show a
  `+1` artifact;
- scheduler topology: one 4-CPU dsub task per sample × chromosome, ten scales
  sequentially within that task, followed by one merge task per sample;
- CPU-driven memory allocation on this dsub system and how to adjust the
  configurable CPU value;
- intermediate and final output locations, including the required final TSV
  schema;
- expected task counts, failure isolation, and the safe cleanup/restart
  command with dry-run first.

## Implementation sequence

1. Move the CWorld IS settings into the documented configuration contract and
   validate them early.
2. Refactor the sample-wide checkpoint into chromosome-level CWorld jobs plus
   a per-sample merge job, retaining the final Meltron input path.
3. Add the direct cooler-resolution guard, manifest validation, and safe
   explicit-sample cleanup helper.
4. Update tests, dry-run assertions, and README documentation.
5. Run static checks, targeted tests, and a representative Snakemake dry run
   before any full dsub resubmission.
