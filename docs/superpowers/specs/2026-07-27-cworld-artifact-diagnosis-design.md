# CWorld Artifact Diagnosis Design

## Goal

Make the single-chromosome CWorld smoke reproducible and inspectable so a
missing bedGraph can be diagnosed without relying on a compute node's local
`/tmp` directory.

## Accepted requirements

- Retain smoke artifacts below the feature checkout's `tmp/` directory.
- Test one sample, one chromosome, and one insulation window.
- Verify the CWorld matrix header before Perl is invoked.
- Preserve CWorld stdout, stderr, and every generated file.
- Request four dsub CPUs; do not add an explicit memory limit.

## Design

The smoke harness uses
`$PROJ/tmp/cworld_smoke/L021_HC04_D08_is32000`. The validator receives that
directory as `--outdir`, invokes Perl in verbose mode, writes both captured
streams to files, and records all files sharing the CWorld output prefix before
it fails closed for an absent bedGraph.

## Acceptance criteria

1. The smoke CWorld matrix and validation tree are retained below `tmp/`.
2. A regression test asserts the shared path and diagnostic logs.
3. The smoke log contains the matrix header check, Perl stderr, and artifact
   listing even when bedGraph is absent.
