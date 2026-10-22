#!/bin/bash
#DSUB -R 'cpu=10'
#DSUB -e %J.err
#DSUB -o %J.out
#DSUB -n SNAPHiCDiff_background_build
#DSUB -aa
#DSUB -q arm

# SNAP-HiC-Diff Background Build Pipeline
# Submission script for Snakemake with dsub cluster executor

set -euo pipefail

WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "${WORKDIR}"

MELTRON_CLI="${WORKDIR}/env/Meltron/command_line_apps/Meltron.R"
if [[ ! -f "${MELTRON_CLI}" ]]; then
    echo "ERROR: Meltron source is missing: ${MELTRON_CLI}" >&2
    echo "Initialize it with: git submodule update --init --recursive" >&2
    exit 1
fi

SNAKEFILE="${SNAKEFILE:-Snakefile}"
CONFIG="${CONFIG:-config/config.yaml}"
DSUB_SUBMIT_CMD="${DSUB_SUBMIT_CMD:-config/dsub_submit.sh}"

# Conda environments are built on the login node; compute nodes must stay offline.
export CONDA_PREFIX_DIR="${WORKDIR}/.snakemake/conda"
bash "${WORKDIR}/scripts/conda_envs.sh" verify
export CONDA_OFFLINE=true

# Create log directories
mkdir -p logs/cluster

# NOTE: Snakemake --use-conda handles environment activation automatically.
# No need to manually activate conda environments.

# Run snakemake with cluster-generic executor
snakemake \
    --executor cluster-generic \
    --use-conda \
    --conda-prefix "${CONDA_PREFIX_DIR}" \
    -j 1000 \
    -s "${SNAKEFILE}" \
    --configfile "${CONFIG}" \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 180 \
    --cluster-generic-submit-cmd "bash ${DSUB_SUBMIT_CMD} {threads} {rulename}" \
    --cluster-generic-status-cmd "bash config/dsub_status.sh" \
    --cluster-generic-cancel-cmd dkill \
    --cluster-generic-cancel-nargs 1 \
    all
