# Final Review Fix Report

## Scope

- Fixed the CWorld layout contract in `scripts/validate_cworld_is.py`.
- Added `test/cworld_is/test_validate_cworld_is_paths.py`.
- Strengthened `test/cworld_is/test_compute_is_scheduling.py` to count every requested CWorld invocation.
- Did not modify production scheduling code, cluster configuration, execution tracking, or caches.

## RED

Command:

```text
pytest -q test/cworld_is/test_validate_cworld_is_paths.py::test_validator_uses_converter_cworld_layout_for_single_sample
```

Output:

```text
FAILED ...
assert '/tmp/.../validation/S1__HC04_A01' == '/tmp/.../validation/cworld_insulation/S1/S1__HC04_A01'
1 failed in 8.75s
```

The focused test calls `validate_cworld_is.run()` with only the Cooler and external conversion/CWorld boundaries mocked. It records the CWorld `-o` argument, runs the real converter subprocess, verifies the converter receives the root output directory, and requires its expected input layout to be discoverable.

## GREEN

Commands and outputs:

```text
pytest -q test/cworld_is/test_validate_cworld_is_paths.py::test_validator_uses_converter_cworld_layout_for_single_sample
1 passed in 2.59s

pytest -q test/cworld_is/test_convert_cworld_is.py test/cworld_is/test_compute_is_scheduling.py
3 passed in 6.97s

pytest -q test/cworld_is/test_validate_cworld_is_paths.py test/cworld_is/test_convert_cworld_is.py test/cworld_is/test_compute_is_scheduling.py
4 passed in 11.56s
```

## Change Summary

- The validator now creates `<outdir>/cworld_insulation/<sample>/` before invoking CWorld.
- Its CWorld matrix and prefix now use `<outdir>/cworld_insulation/<sample>/<sample>__<chrom>`.
- The validator still passes root `<outdir>` to `convert_cworld_IS_to_meltron_IS.py`.
- The scheduling dry-run test now asserts exactly ten `matrix2insulation.pl` appearances for the ten requested samples and `compute_IS=10`.

## Self-Review

- The matrix/prefix and converter discovery layout now agree.
- The regression test observes the validator's actual `run()` command construction and uses the actual converter for discovery and output creation.
- The parent directory is asserted to exist before the CWorld invocation.
- The change is limited to the requested validator and tests; no production scheduling behavior changed.
- `git diff --check` passed before final verification.

## Commit

`fix(is): align cworld validation layout`

## dsub

No `dsub` command was run. The scheduling test uses Snakemake dry-run mode only.
