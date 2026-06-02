#!/bin/bash
#DSUB -q arm
#DSUB -aa
#DSUB -n <job-name>
#DSUB -e %J.err
#DSUB -o %J.out
#DSUB -R "cpu=10"

set -euo pipefail

# =============================================================================
# Disable nested parallelism — CRITICAL for HPC efficiency
# =============================================================================
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# =============================================================================
# Configure paths (project-specific)
# =============================================================================
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# =============================================================================
# Activate conda environment
# =============================================================================
source conda activate snakemake-snaphic

# =============================================================================
# Your commands here
# =============================================================================
python3 -m pytest tests/ -v
