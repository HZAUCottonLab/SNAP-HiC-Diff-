# ADR-0002: Build a downstream unified TAD-boundary view for IS annotation

## Status

Accepted 2026-09-03.

## Decision

Keep the raw and balanced 3DNetMod TAD and multi-resolution boundary outputs
separate as required by ADR-0001. Add a downstream, provenance-preserving union
view used only for IS-change annotation. Overlapping raw/balanced boundaries
may become `both`; unmatched boundaries remain `raw_only` or `balanced_only`.
The union view never overwrites or aliases either source output.

## Consequences

Scientific users can audit normalization sensitivity while the comparison
stage receives one nonduplicated boundary collection per sample. Downstream
evidence tiers retain signal and resolution support, and complex overlap
components remain explicit rather than being forced into one boundary.
