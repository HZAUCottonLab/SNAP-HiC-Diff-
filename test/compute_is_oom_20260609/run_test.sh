#!/bin/bash
#DSUB -q arm
#DSUB -aa
#DSUB -n test-compute-is
#DSUB -e %J.err
#DSUB -o %J.out
#DSUB -R "cpu=4"

set -euo pipefail

PROJ=/share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff
export CONDA_ENVS_PATH="${PROJ}/.snakemake/conda"

CONDA_PREFIX="${CONDA_ENVS_PATH}/59271960b794774df1c934b1e7affef8_"
export PATH="${CONDA_PREFIX}/bin:$PATH"

echo "=== Test Start: $(date) ==="
echo "Python: $(which python3)"

cd "$PROJ"

# Test on L021 (17935 bins - small but real)
python3 "$PROJ/scripts/convert_cool_to_cworld.py" \
    --cool "$PROJ/results/balanced/L021/L021__HC04_D08_iced.cool" \
    --chrom HC04_D08 \
    --out /tmp/L021__HC04_D08_test.matrix.gz \
    --log /tmp/test_convert.log

if [ ! -f /tmp/L021__HC04_D08_test.matrix.gz ]; then
    echo "ERROR: output file not created"
    exit 1
fi

SIZE=$(stat -c%s /tmp/L021__HC04_D08_test.matrix.gz)
if [ "$SIZE" -eq 0 ]; then
    echo "ERROR: output file is empty"
    exit 1
fi

echo "Output size: $SIZE bytes"

# Verify gzip integrity
if ! gzip -t /tmp/L021__HC04_D08_test.matrix.gz 2>&1; then
    echo "ERROR: gzip integrity check failed"
    exit 1
fi

# Count lines
LINES=$(zcat /tmp/L021__HC04_D08_test.matrix.gz | wc -l)
echo "Output lines: $LINES (expected 17936 = 1 header + 17935 data rows)"
if [ "$LINES" -ne 17936 ]; then
    echo "ERROR: unexpected line count"
    exit 1
fi

echo "=== Test End: $(date) ==="
echo "TEST PASSED"
