# IS-change plot: aligned Hi-C heatmap panels

## Goal

Extend each primary IS-change diagnostic plot with two sample-specific Hi-C heatmaps above the existing multi-scale IS curves. The heatmaps and all IS panels must use exactly the same genomic x-coordinate interval so an event can be located consistently in every panel.

## Confirmed decisions

- Layout is vertical: sample A heatmap, sample B heatmap, then the existing ten IS-scale panels.
- Both heatmaps show the complete symmetric contact matrix for the existing plot interval: the candidate region plus the configured 300 kb context on both sides, rounded to 4 kb bins.
- All panels share that interval as their x-axis. Only the bottom IS panel shows x tick labels.
- Heatmaps use the per-chromosome `*_iced.cool` inputs with balanced contacts, transformed as `log1p(balanced contact)`.
- The two heatmaps within one figure use one shared color range: `vmin=0` and `vmax` is the 99.5th percentile of their pooled, finite `log1p(balanced contact)` values. A shared colorbar is labelled `log1p(balanced contact)`.
- Full matrices are rendered in fixed-height panels. Their genomic coordinates remain correct, while the y direction is visually compressed to keep the 12-panel diagnostic figure usable.
- Non-finite contact cells and cells involving unavailable bins are shown in light gray.
- Existing red primary-event markers are retained. Point events receive vertical and horizontal summit lines on both matrices; broad-shift events receive translucent red vertical and horizontal bands spanning the called interval. Existing IS-panel markers are unchanged.

## Data flow

1. The current plotter reads the primary-event TSV, candidate-region TSV, two IS tracks, and the per-chromosome iced cooler inputs.
2. It creates a chromosome-to-cooler mapping for each sample and validates that the primary region chromosome is available once per sample.
3. For a plotted region, it computes the existing bin-aligned plot interval once and uses it for every IS scale and both matrix fetches.
4. It fetches only the two local cis matrices for that interval, applies balancing and `log1p`, derives the specified joint shared color limit, and renders the two heatmaps before the IS axes.
5. It writes the same PNG and manifest schema as before; no caller, candidate, or event TSV schema changes.

## Robustness requirements

- Fail clearly if an input cooler does not provide the required chromosome/interval or if a matrix cannot be rendered with any finite non-negative balanced contacts.
- Do not write a partial PNG when rendering fails; retain the existing atomic output behavior.
- Keep existing IS-track and balanced-bin validation, unavailable-bin shading, primary-only selection, ranking, and manifest semantics.
- Do not add a difference heatmap, additional event-calling logic, new output tables, or separate per-sample color scales.

## Validation

Add synthetic-cooler coverage that verifies:

1. both heatmaps fetch the same genomic interval as the IS panels;
2. their matrix extents and color limits are identical within a figure;
3. sample labels and summit/broad-shift overlays render without error;
4. unavailable matrix cells use the configured missing-data rendering;
5. existing primary-only rejection and manifest behavior remain intact.

Run the full IS-change test suite through dsub on an ARM compute node. A real H1398-vs-L134 replay should visually confirm that the A05 event coordinate aligns between both heatmaps and IS curves.

## Acceptance criteria

For every diagnostic PNG, a reader can draw a vertical line from any heatmap x-coordinate to every IS panel and obtain the same genomic coordinate. The two heatmaps are quantitatively comparable by color, event markers fall at the same x-coordinate as the IS event markers, and failures do not replace prior valid figures.
