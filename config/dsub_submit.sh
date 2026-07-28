#!/usr/bin/env bash
# SNAP-HiC-Diff Background Build Pipeline
# dsub submit wrapper for Snakemake --executor cluster-generic

set -euo pipefail

_threads="${1:?missing threads}"
_rulename="${2:?missing rulename}"
_jobscript="${3:?missing jobscript}"

# Disable nested parallelism so each ProcessPoolExecutor worker uses exactly 1 core
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=256
export OMP_MAX_ACTIVE_LEVELS=1

# Resolve jobscript to absolute path
if [[ "$_jobscript" = /* ]]; then
    _jobscript_abs="$_jobscript"
else
    _jobscript_abs="$(pwd)/$_jobscript"
fi
_wrapper_dir="$(dirname "$_jobscript_abs")"
_wrapper="${_wrapper_dir}/.dsub_wrapper_$$.sh"

# Build wrapper: exports + original script + self-deletion
{
    printf '#!/usr/bin/env bash\n'
    printf 'set -euo pipefail\n'
    printf 'export CONDA_OFFLINE=true\n'
    cat "$_jobscript"
    printf 'rm -f "%s"\n' "$_wrapper"
} > "$_wrapper"
chmod +x "$_wrapper"

# Get project root (parent of workflow/profile)
PROJECT_ROOT="$(cd "${_wrapper_dir}/../.." && pwd)"

# Set conda environments path for offline cluster (environment created on login node)
export CONDA_ENVS_PATH="${PROJECT_ROOT}/.snakemake/conda"

# Submit to dsub
msg=$(dsub \
    --cwd "${PROJECT_ROOT}" \
    -n "${_rulename}" \
    -aa \
    -q arm \
    -R "cpu=${_threads}" \
    -o "logs/cluster/${_rulename}_%J.out" \
    -e "logs/cluster/${_rulename}_%J.err" \
    bash "$_wrapper"
)

# Parse job ID from dsub output
if [[ "$msg" =~ \<([0-9]+)\> ]]; then
    echo "${BASH_REMATCH[1]}"
else
    echo "Cannot parse dsub job id from: $msg" >&2
    exit 1
fi
