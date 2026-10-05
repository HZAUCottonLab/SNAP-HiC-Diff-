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

## Unified boundary region

A downstream-only, provenance-preserving union of actually overlapping raw and
balanced multi-resolution boundary regions for one sample. The unified view is
used only for IS-change annotation and never overwrites or aliases its source
outputs.

## Boundary position state

The downstream comparison classification of a unified boundary between H1398
(A) and L134 (B) on the same reference genome and the same 4 kb grid. The
approved states are `conserved`, `gained_in_B`, `lost_in_B`, `shifted`,
`complex`, and `ambiguous`; each retains source normalization provenance.

## Shifted boundary

A shifted boundary is an ordered A→B match only when
`J < abs(center_B - center_A) <= Dmax`. A matched boundary within J is
`conserved`. The shifted classification preserves both positions and does not
replace either raw or balanced upstream boundary output.

## IS boundary strength

A downstream IS-change annotation measure for a unified boundary, evaluated
against local reference boundaries across the configured IS scales. It retains
the boundary's signal and resolution support.
