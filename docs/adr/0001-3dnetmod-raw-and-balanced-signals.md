# ADR-0001: Publish raw and balanced 3DNetMod calls separately

## Status

Accepted 2026-08-21.

## Decision

Run 3DNetMod at 4 kb, 8 kb, and 20 kb twice per sample with identical D
parameters. Raw runs use finite positive pixel counts. Balanced runs use
`count * weight(bin1) * weight(bin2)` and skip non-finite or non-positive
results. Boundary merging occurs within each signal and never combines raw
with balanced calls.

## Consequences

The workflow produces two scientifically distinguishable result families and
216 independent calls for the current 36 samples. Existing mcools remain
read-only; missing weights fail rather than trigger implicit rebalancing.
