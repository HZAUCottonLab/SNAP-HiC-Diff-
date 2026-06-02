# Journal - xhhuang24

---

## Session 1: Trellis init — port knowledge from background_pool

**Date**: 2026-05-27
**Task**: trellis-init
**Branch**: —

### Summary

Ported reusable knowledge from `background_pool/docs/` and `background_pool/.trellis/` into SNAP-HiC-Diff project.

### Main Changes

- Created `.trellis/` directory with spec, scripts, workspace, and runtime
- Ported `docs/dsub_hpc_standard.md` (HPC cluster execution guide)
- Ported `docs/snakemake/wildcards.md` (adapted for project patterns)
- Ported `docs/templates/submit-backfill-ncc.sh`
- Ported `.trellis/spec/pipeline/` (snakemake-guidelines, dsub-guidelines, cluster-submission-ssot, hictk, hicexplorer-knowledge, pipeline index)
- Ported `.trellis/spec/bug-analysis/` (7 bug analysis docs)
- Ported `.trellis/spec/guides/` (3d-genomics, code-reuse, cross-layer thinking guides)
- Ported `.trellis/scripts/` (all Python scripts)
- Ported `.trellis/.runtime/` and `.trellis/config.yaml`
- Created `.trellis/.developer` (identity: xhhuang24)
- Created `.trellis/.version` (1.0.0)

### Non-Ported (selection rationale)

- `docs/project_design.md` — overlaps with existing `CODEBASE_OVERVIEW.md`
- `docs/rule_index.md` — full rule index, already covered by `RULE_SCRIPT_CONFIG_MAP.md`
- `docs/snakemake_v8_30_0_standard.md` — superseded by `.trellis/spec/pipeline/snakemake-guidelines.md`
- `.trellis/spec/backend/` — no backend layer in this project
- `.trellis/spec/frontend/` — no frontend layer in this project
- `.trellis/tasks/` — task records are developer-specific, not reusable
- `.trellis/workspace/` (except template) — already created for xhhuang24

### Status

[OK] **Completed**

### Next Steps

- Verify trellis scripts work: `python3 ./.trellis/scripts/get_developer.py`
- Run `codegraph init -i` to build CodeGraph index
