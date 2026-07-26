#!/bin/bash
#DSUB -q arm
#DSUB -aa
#DSUB -n test-compute-is
#DSUB -e %J.err
#DSUB -o %J.out
#DSUB -R "cpu=4"

set -euo pipefail

PROJ="${PROJ:-/share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff}"
DATA_PROJ="${DATA_PROJ:-$PROJ}"
ENV_PROJ="${ENV_PROJ:-$PROJ}"
export CONDA_ENVS_PATH="${ENV_PROJ}/.snakemake/conda"

CONDA_PREFIX="${CONDA_ENVS_PATH}/59271960b794774df1c934b1e7affef8_"
export PATH="${CONDA_PREFIX}/bin:$PATH"

echo "=== Test Start: $(date) ==="
echo "Python: $(which python3)"

cd "$PROJ"

# Test on L021 (17935 bins - small but real)
python3 "$PROJ/scripts/convert_cool_to_cworld.py" \
    --cool "$DATA_PROJ/results/balanced/L021/L021__HC04_D08_iced.cool" \
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

python3 "$PROJ/scripts/validate_cworld_is.py" \
    --cool "$DATA_PROJ/results/balanced/L021/L021__HC04_D08_iced.cool" \
    --sample L021 \
    --chrom HC04_D08 \
    --is-distance 32000 \
    --cworld-dir "$PROJ/env/cworld-dekker" \
    --outdir /tmp/cworld_validate_L021_D08

MELTRON_TSV=/tmp/cworld_validate_L021_D08/IS/L021/L021_IS.tsv
if [ ! -s "$MELTRON_TSV" ]; then
    echo "ERROR: Meltron TSV is missing or empty: $MELTRON_TSV"
    exit 1
fi

EXPECTED_HEADER=$'chrom\tviewpoint\tIS_distance\tIS'
ACTUAL_HEADER="$(head -n 1 "$MELTRON_TSV")"
if [ "$ACTUAL_HEADER" != "$EXPECTED_HEADER" ]; then
    echo "ERROR: unexpected Meltron TSV columns: $ACTUAL_HEADER"
    exit 1
fi

echo "Meltron TSV columns: $ACTUAL_HEADER"

echo "=== Test End: $(date) ==="
echo "TEST PASSED"
