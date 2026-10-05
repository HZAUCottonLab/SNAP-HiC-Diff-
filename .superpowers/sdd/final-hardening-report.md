# Final hardening evidence

- 2026-08-23: archived every pre-existing 3DNetMod run workspace before successful preparation, retaining it under `.previous_attempts/prepare-<uuid>` while creating a fresh active workspace.
  - Regression: changed D configuration archives old `FINAL_DOMAIN_CALLS` and GPS diagnostics before the next prepare.
  - Regression: vendored source contract prohibits `taskset` masks and fixed pool sizes; each active pool is constrained by `settings["processors"]` and its stage workload.
  - Validation: fixed D scalars, five-value arrays, and 4-kb overrides are validated before compute.
  - Evidence: `pytest -q test/tad` — 63 passed; changed-module `py_compile`, `git diff --check`, and no-affinity source scan passed.
