# Multi-resolution 3DNetMod TAD and Boundary Design

## Status

Approved on 2026-08-21.

## Goal

Add a self-contained 3DNetMod TAD capability to SNAP-HiC-Diff. For every
sample, the default workflow will call TADs at 4 kb, 8 kb, and 20 kb from both
raw and balanced contact values. It will publish per-resolution TAD and
boundary files plus one multi-resolution merged-boundary file for each signal.

The implementation will migrate the already validated Python 3 3DNetMod
workflow from `/share/home/xhhuang24/project/HiC_analysis_pipeline`; the
finished workflow must not depend on that checkout at runtime.

## Scope

- Vendor the complete Python 3 3DNetMod dependency closure and seed resource.
- Run the canonical 3DNetMod D parameter set at 4 kb, 8 kb, and 20 kb.
- Run raw and balanced analyses independently with identical D parameters.
- Generate and merge boundary intervals within and across resolutions.
- Include every stable final output in `rule all`.
- Build the dedicated TAD Conda environment on the login node and consume it
  offline on compute nodes.
- Preserve signal-specific run directories, intermediate results, and logs for
  audit and failure diagnosis.

Out of scope:

- Comparing or merging raw and balanced calls into a consensus result.
- Rebalancing an existing input mcool in place.
- Running alternative 3DNetMod parameter grids.
- Filtering chromosomes by name.
- Enabling a bad-region mask without a project-specific BED.
- Launching the complete 216-call production run as part of implementation
  verification.

## Existing Data Contract

The active manifest contains 36 mcool inputs. All currently contain 4 kb,
8 kb, and 20 kb resolution groups and a `bins/weight` dataset at each of those
resolutions. The workflow must still validate these properties for every
sample rather than relying on this observation.

The current conversion configuration omits 8 kb, and `cooler zoomify` is
invoked without balancing. Future `.cool` inputs therefore require both
changes:

- add 8000 to `convert.resolutions`;
- add `--balance` when building a new mcool.

Existing mcool inputs remain read-only. If an existing mcool lacks a required
resolution or weight dataset, the corresponding TAD job fails with a stable
error code instead of modifying the file.

## Scientific Input Definitions

For a pixel with raw contact count `c` and endpoint-bin weights `w1` and `w2`:

```text
raw_contact      = c
balanced_contact = c * w1 * w2
```

Raw mode exports finite positive counts. Balanced mode requires the weight
dataset and exports finite positive balanced values; pixels with non-finite
weights or results are skipped. A signal with no valid exported interactions
is an error.

Both modes use the same D parameters, including the existing 4 kb override.
Only the interaction values, input filename, signal label, and run directory
differ. All chromosomes present in the selected mcool resolution are analyzed;
`excluded_chromosomes` defaults to an empty list. `badregionfilter` defaults to
false and `bad_region_bed` defaults to null.

## Configuration Contract

The existing configuration gains 8 kb conversion and a TAD section following
the repository's component-local thread convention:

```yaml
convert:
  resolutions: [4000, 8000, 10000, 20000, 100000]

tad:
  resolutions: [4000, 8000, 20000]
  signals: [raw, balanced]
  threads: 20
  parameter_set: D
  run_root: results/tad_runs
  code_dir: scripts/python3_3dnetmod
  badregionfilter: false
  bad_region_bed: null
  excluded_chromosomes: []
  d:
    region_size_bins: 150
    overlap_bins: 100
    logged: true
    scale: genomewide
    plateau: 3
    chaosfilter: true
    chaos_pct: 0.9
    max_gamma: 3.2
    num_part: 20
    plots: false
    pctile_threshold: 0
    pct_value: 0
    final_consistent_domains: true
    size_threshold_bins: 4
    boundary_buffer_bins: 8
    size_bins: [100, 200, 400, 750, 3000]
    variance_thresholds: [45, 45, 60, 45, 45]
    resolution_overrides:
      4000:
        size_threshold_bins: 10
        variance_thresholds: [20, 20, 30, 20, 20]
```

Resolution and signal lists must be non-empty, duplicate-free, and restricted
to the supported values. TAD resolutions must exist in `convert.resolutions`.
The parameter set must be D, threads must be positive, and all required D
fields and five-element arrays must validate before the DAG schedules a
compute-heavy job.

## Components and Data Flow

### Snakemake orchestration

`rules/tad.smk` defines independently retryable work over
`sample x resolution x signal`:

1. prepare signal-specific 3DNetMod input and settings;
2. call 3DNetMod and publish TAD BED/BEDPE;
3. extract and merge boundaries for one resolution;
4. merge the three resolution-specific boundary sets for one sample and
   signal.

All TAD, per-resolution boundary, and multi-resolution boundary outputs are
included in `rule all`. With the current 36 samples, two signals, and three
resolutions, the DAG contains 216 independent 3DNetMod calls. Users may still
request a single output target to run a subset.

### Input adapter

`scripts/tad_from_mcool.py` is migrated and extended with a required signal
argument. It validates the mcool layer, exports global one-based bin IDs,
writes raw or balanced interaction triples, renders one D settings file, and
records a ready manifest. Raw and balanced artifacts never share a mutable run
directory.

### Vendored algorithm and runner

`scripts/python3_3dnetmod/` contains the complete migrated Python 3 code,
recursive imports, seeds, and static resources. `scripts/run_3dnetmod.py`
invokes these stages in order:

1. `3DNetMod_preprocess.py`;
2. `3DNetMod_GPS_MMCP.py`;
3. `3DNetMod_HSVM.py`.

The stable run directory is:

```text
results/tad_runs/{sample}/{signal}/{resolution}kb/D/
```

Failures preserve this directory and its stage output. Successful collection
publishes canonical BED3 and BEDPE6 files only after validating the final
domain coordinates.

### Boundary processor

`scripts/tad_boundaries.py` performs both resolution-specific extraction and
multi-resolution merging through small, directly tested functions. Keeping
the two operations in one script avoids duplicate interval parsing and merge
implementations while Snakemake retains separate rules and outputs.

The data flow is:

```text
manifest mcool
  -> signal-specific interactions and D settings
  -> 3DNetMod TAD BED/BEDPE
  -> per-resolution merged boundary BED
  -> per-signal multi-resolution merged boundary BED
```

## Output Contract

For sample `S`, resolution label `Rkb`, and signal `raw` or `balanced`:

```text
results/tads/S.Rkb.signal.3dnetmod.bed
results/tads/S.Rkb.signal.3dnetmod.bedpe
results/tad_boundaries/S.Rkb.signal.boundaries.bed
results/tad_boundaries/S.multires.signal.boundaries.bed
```

TAD BED is canonical BED3. TAD BEDPE repeats the same interval in both BEDPE
halves, matching the migrated runner's existing contract.

Resolution-specific boundary BED has six tab-separated columns and no header:

```text
chrom  start  end  boundary_id  resolution  n_source_boundaries
```

Multi-resolution boundary BED has seven tab-separated columns and no header:

```text
chrom  start  end  boundary_id  support_resolutions
n_supporting_resolutions  n_source_boundaries
```

`support_resolutions` is a comma-separated ascending list of distinct integer
base-pair resolutions. `n_source_boundaries` counts original TAD endpoints,
not already merged intervals. Boundary IDs are deterministic after sorting by
the chromosome order stored in the mcool, then start and end.

## Boundary Semantics

For a TAD `[tad_start, tad_end)` called at resolution `R`, the two endpoint
coordinates are expanded symmetrically by one resolution:

```text
[max(0, tad_start - R), min(chrom_length, tad_start + R))
[max(0, tad_end   - R), min(chrom_length, tad_end   + R))
```

Left and right endpoints are intentionally not distinguished. Chromosome
length and order come from the active mcool.

Intervals are merged independently by chromosome. Resolution-specific files
first merge all overlapping intervals at that resolution. The final file then
merges the 4 kb, 8 kb, and 20 kb interval sets for the same signal. Merging is
transitive and takes the interval union: the output start is the minimum start
and the output end is the maximum end of the connected overlap component.

BED half-open semantics are strict:

```text
[0, 8000) and [4000, 12000) -> [0, 12000)
[0, 8000) and [8000, 16000) -> remain separate
```

Therefore the merge condition is `next_start < current_end`, never less than
or equal. Transitive overlap may intentionally produce a broad boundary region;
the design does not impose an unrequested maximum width or split heuristic.
Raw and balanced boundaries are never merged with each other.

## Validation, Errors, and Atomic Publication

Each error message begins with a stable code and dynamically reports the
actual sample, signal, resolution, stage, and path where applicable. Required
codes are:

```text
TAD_CONFIG_INVALID
TAD_MISSING_RESOLUTION
TAD_MISSING_WEIGHT
TAD_EMPTY_INTERACTIONS
TAD_STAGE_FAILED
TAD_EMPTY_CALLS
TAD_INVALID_COORDINATE
TAD_BOUNDARY_MERGE_FAILED
```

No sample name is hard-coded in diagnostic output. A stage failure includes
the stage name and process return code. Unknown chromosomes, non-integer
coordinates, non-positive intervals, chromosome-boundary violations, missing
outputs, and empty final calls fail explicitly.

Final TAD and boundary files are written to sibling temporary files, validated,
and atomically replaced only after success. A failed job leaves its run-local
intermediates and log but never publishes a partial stable output.

Signal- and resolution-specific logs are:

```text
results/logs/tad_prepare/{sample}.{resolution}kb.{signal}.log
results/logs/tad_call/{sample}.{resolution}kb.{signal}.log
results/logs/tad_boundaries/{sample}.{resolution}kb.{signal}.log
results/logs/tad_boundaries/{sample}.multires.{signal}.log
```

## Environment and Scheduling

`env/tad.yaml` isolates the migrated compatibility set from the existing main
analysis environment. It pins Python 3.11, NumPy below 2, Matplotlib 3.8, and
declares Cooler, HDF5, SciPy, scikit-learn, Pandas, and pytest.

The environment specification is added to `scripts/conda_envs.sh`. Bootstrap
creates and validates it only on the login node:

```text
bash scripts/conda_envs.sh bootstrap
```

Compute nodes remain offline and execute only the existing verification path
before Snakemake activates the completed prefix. No worker downloads packages
or creates environments.

Each 3DNetMod call uses 20 threads. The existing dsub wrapper submits those
threads to the arm queue, whose CPU-linked allocation provides approximately
100 GB for the job. Raw and balanced calls are separate jobs, so either can be
retried without rerunning the other.

## Verification

Unit and contract tests cover:

- raw positive-count export;
- balanced contact calculation and non-finite weight filtering;
- stable errors for missing resolution, missing weight, and empty contacts;
- identical effective D settings for raw and balanced signals;
- the 4 kb D override;
- global one-based bin IDs across chromosomes;
- vendored entrypoints, recursive imports, and seed resources;
- ordered stage execution, stage-failure context, and atomic final collection;
- boundary expansion and chromosome-end clipping;
- duplicate, same-resolution, cross-resolution, and transitive overlap merge;
- non-merging endpoint adjacency;
- deterministic IDs and provenance counts;
- exact output and log path construction;
- inclusion of all sample/resolution/signal outputs in `rule all`;
- dedicated TAD environment bootstrap and offline verification.

A Snakemake fixture dry-run must build the complete raw/balanced three-
resolution DAG without submitting work. Final acceptance runs one real sample
through all six TAD combinations and verifies non-empty TADs, six
resolution-specific boundary files, and two multi-resolution boundary files.
The complete 36-sample production run is left to the normal workflow launch.

## Acceptance Criteria

The feature is complete when:

1. the repository runs without access to the source HiC analysis checkout;
2. login-node bootstrap creates and verifies the dedicated TAD environment;
3. Snakemake dry-run schedules raw and balanced 4/8/20 kb calls for every
   manifest sample and includes their stable outputs in `rule all`;
4. the test suite covers signal calculation, boundary semantics, stable errors,
   vendored closure, and workflow contracts;
5. a real-sample smoke test publishes valid, non-empty TAD and boundary
   artifacts for both signals;
6. failures remain searchable by stable dynamic error code and never publish
   partial stable outputs.
