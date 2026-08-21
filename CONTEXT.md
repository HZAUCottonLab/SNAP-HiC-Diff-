# SNAP-HiC-Diff Domain Context

## TAD analysis

Per-sample identification of topologically associating domain intervals from a
Hi-C contact map at 4 kb, 8 kb, or 20 kb.

## Raw TAD signal

Finite positive mcool pixel counts used directly as 3DNetMod interaction
values.

## Balanced TAD signal

Finite positive values calculated as pixel count multiplied by both endpoint
bin weights. Balanced calls are separate results, not replacements for raw
calls.

## TAD boundary region

A TAD endpoint expanded left and right by one active resolution and clipped to
the chromosome. Left and right endpoints are not distinguished.

## Multi-resolution boundary region

The transitive interval union of actually overlapping 4 kb, 8 kb, and 20 kb
boundary regions for one sample and one signal. Endpoint adjacency alone does
not merge regions.
