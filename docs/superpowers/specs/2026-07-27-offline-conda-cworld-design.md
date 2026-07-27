# Offline Conda and CWorld execution design

**Date:** 2026-07-27

## Goal

Run the complete SNAP-HiC-Diff workflow through dsub on ARM compute nodes that
have no network access.  All Conda environments must be created on a login
node before submission; neither the dsub controller job nor any submitted rule
job may download packages.

## Evidence and root cause

- Job `38519427` failed before submitting workflow rules because Snakemake
  attempted to create `env/cworld.yaml` on an offline compute node.  The node
  inherited a proxy pointing to `127.0.0.1:47890`, where no proxy was running.
- The same environment cannot be created on the login node as written:
  `perl-file-slurp>=9999.34` has no candidate for `linux-aarch64`; the
  available Bioconda package version is at most `9999.32`.
- The vendored CWorld Perl code uses only Perl core modules and the vendored
  `cworld::dekker` library for this insulation path.  It does not import
  `File::Slurp` or the other non-core Perl modules currently listed.
- There are three Conda specifications (`SNAP-HiC-Diff.yaml`, `cworld.yaml`,
  and `meltron.yaml`) and no completed prefixes under `.snakemake/conda`.
  Preparing only CWorld would defer the same offline failure to later rules.

## Chosen design

### 1. Make the CWorld specification ARM-solvable and minimal

Keep the Python runtime packages needed by `convert_cool_to_cworld.py`
(`python`, `cooler`, `numpy`, and `pandas`) and the Perl runtime used by the
vendored `matrix2insulation.pl` (`perl=5.32`, `perl-gd`, and `zlib`).  Remove
the unavailable and unused `perl-file-slurp` constraint, plus unreferenced
non-core Perl package constraints.  Pinning Perl to 5.32 matches the
linux-aarch64 package family and avoids an unconstrained solver choice.

### 2. Prebuild every workflow environment on the login node

Add a project bootstrap command that runs Snakemake's own
`--conda-create-envs-only` path for the default DAG, so each environment lands
in the exact hash-derived prefix Snakemake later activates.  The bootstrap
will:

1. create or refresh all three environments under `.snakemake/conda`;
2. verify each prefix has `conda-meta/history` and its cached YAML matches the
   source specification checksum;
3. perform runtime checks appropriate to the environment, including Python
   imports for CWorld and a Perl import of `cworld::dekker` using the vendored
   library; and
4. write a single ready marker containing the source checksums and validated
   prefix paths.

An incomplete prefix is never treated as ready.  A changed environment YAML
invalidates the marker and requires rerunning bootstrap on a login node.

### 3. Enforce offline execution at dsub time

`run_snakemake.sh` will validate the ready marker before starting Snakemake,
set `CONDA_ENVS_PATH` to the project cache, and export `CONDA_OFFLINE=true`.
`config/dsub_submit.sh` will set the same two variables in every generated
rule-job wrapper.  A missing, corrupt, or stale prebuilt environment will fail
immediately with instructions to run bootstrap on a login node instead of
attempting a network operation on a compute node.

The implementation will not reuse the reference project's environment paths:
that would create an undeclared cross-project dependency and risk incompatible
package versions.  Its checked bootstrap/marker/offline pattern is reused, not
its full package set.

## Data flow

```text
login node: bootstrap all 3 Conda envs -> validate -> ready marker
                                             |
                                             v
dsub controller: verify marker + CONDA_OFFLINE=true
                                             |
                                             v
per-rule dsub wrapper: same cache + offline -> activate existing prefix only
```

## Error handling

- A solver failure on the login node leaves no ready marker and reports the
  offending environment specification.
- A missing/stale marker or invalid prefix makes the submitted job exit before
  DAG execution, with the bootstrap command shown.
- `CONDA_OFFLINE=true` ensures a missed prebuild cannot silently trigger a
  package download attempt on a compute node.

## Acceptance criteria

1. `env/cworld.yaml` solves on `linux-aarch64` without
   `perl-file-slurp`.
2. The bootstrap creates and validates all three Conda prefixes and writes a
   marker whose checksums match current YAML files.
3. With `CONDA_OFFLINE=true`, the default Snakemake dry-run and CWorld runtime
   checks succeed without channel access.
4. A deliberately absent or stale marker causes `run_snakemake.sh` to fail
