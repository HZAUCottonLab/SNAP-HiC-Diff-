# dsub HPC Standard

## Overview

This document defines the standards for running the SNAP-HiC-Diff pipeline on HPC clusters using dsub. All job submissions must follow these conventions.

---

## 1. dsub Basic Usage

### Command Syntax

```bash
dsub [options] --command "your command"
dsub [options] --script script.sh
```

### Essential Flags

| Flag | Description | Example |
|------|-------------|---------|
| `--script` | Shell script to execute | `--script job.sh` |
| `--name` | Job name | `--name snakemake-validate` |
| `-q`, `--queue` | Queue name | `-q x86` |
| `-aa`, `--any-arch` | Accept any architecture | `-aa` |
| `-e`, `--env` | Environment variable | `-e MY_VAR=value` |
| `-i`, `--input` | Input file | `-i /data/input.txt` |
| `-o`, `--output` | Output file | `-o /results/output.txt` |
| `-u`, `--user` | Username | `-u $USER` |
| `--dry-run` | Print without executing | `--dry-run` |

### Example Basic Usage

```bash
#!/bin/bash
# Basic dsub job submission

dsub \
    --script my_script.sh \
    --name my-job \
    --queue x86 \
    -aa \
    --env MY_DATA=/data/sample1 \
    --env OUTPUT_DIR=/results \
    --input /data/reference.fasta \
    --output /results/output.txt \
    --log /logs/my-job.log
```

---

## 2. DSUB Header Comments

Shell scripts executed by dsub should include `#DSUB` comments for metadata:

```bash
#!/bin/bash
#DSUB -n my-job                  # Job name
#DSUB -q x86                     # Queue
#DSUB -aa                        # Any architecture
#DSUB -e SNEMAKE_THREADS         # Environment variable
#DSUB -e SNEMAKE_MEM_MB
#DSUB -e OUTDIR
#DSUB -o %o.log                  # stdout
#DSUB -e %e.log                  # stderr
#DSUB -R "cpu=4"                 # Resource requirements
#DSUB -R "mem=8192"              # Memory requirements

set -euo pipefail

# Your commands here
echo "Job started at $(date)"
```

### Header Comment Reference

| Comment | Meaning |
|---------|---------|
| `#DSUB -n name` | Job name |
| `#DSUB -q x86` | Queue selection |
| `#DSUB -aa` | Any architecture allowed |
| `#DSUB -e VAR` | Export environment variable |
| `#DSUB -o path` | stdout file |
| `#DSUB -e path` | stderr file |
| `#DSUB -R "key=value"` | Resource specification |

---

## 3. Wrapper Script Pattern

Wrapper scripts must disable nested parallelism and clean up after themselves.

### Template Wrapper Script

```bash
#!/bin/bash
#DSUB -n snakemake-wrapper
#DSUB -q x86
#DSUB -aa
#DSUB -o %o.log
#DSUB -e %e.log

set -euo pipefail

# =============================================================================
# Disable nested parallelism — CRITICAL for HPC efficiency
# =============================================================================
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# Capture arguments
SCRIPT="$1"
shift

# Execute the actual script
echo "Executing: $SCRIPT $@"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "Start time: $(date)"

bash "$SCRIPT" "$@"

EXIT_CODE=$?

echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"

# Self-delete wrapper after successful execution
if [ $EXIT_CODE -eq 0 ]; then
    rm -f "$SCRIPT"
fi

exit $EXIT_CODE
```

### Why Disable Nested Parallelism?

```
Snakemake (1 thread per rule)
    └── Python script (N threads via joblib)
            └── Libraries (OpenBLAS, MKL, OpenMP)
                    └── Each using M threads
```

Without disabling, you get `N × M` threads per rule, causing:
- CPU oversubscription
- Memory exhaustion
- Poor performance
- Job failures

---

## 4. Snakemake Integration

### Profile Configuration

The dsub profile is located at `profiles/dsub/`:

```
profiles/dsub/
├── config.yaml          # Snakemake profile configuration
└── dsub_submit.sh       # Submission command template
```

### config.yaml

```yaml
executor: cluster-generic
cluster-generic-submit-cmd: ./profiles/dsub/dsub_submit.sh
jobs: 100
latency-wait: 60
keep-going: true
rerun-incomplete: true
```

### dsub_submit.sh

```bash
#!/bin/bash
# =============================================================================
# dsub_submit.sh — dsub submission command for Snakemake
#
# Usage: dsub_submit.sh <script.sh> <rule> <jobid> <threads> <mem_mb> <log_file>
# =============================================================================

set -euo pipefail

SCRIPT="$1"
RULE="$2"
JOBID="$3"
THREADS="${4:-1}"
MEM_MB="${5:-4096}"
LOG_FILE="${6:-/tmp/snakemake.log}"

# Create unique job name
JOBNAME="smk-${RULE}-${JOBID}"

# Disable nested parallelism
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Parse input/output parameters from additional args
INPUT_FILES=""
OUTPUT_FILES=""
ENV_VARS=""

shift 6
while [ $# -gt 0 ]; do
    case "$1" in
        --input)
            INPUT_FILES="$INPUT_FILES --input $2"
            shift 2
            ;;
        --output)
            OUTPUT_FILES="$OUTPUT_FILES --output $2"
            shift 2
            ;;
        --env)
            ENV_VARS="$ENV_VARS --env $2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

# Calculate CPU and memory requirements
CPUS=$(( (THREADS + 1) / 2 ))  # Round up

# Submit to dsub
dsub \
    --script "$SCRIPT" \
    --name "$JOBNAME" \
    --queue x86 \
    --any-arch \
    --env SNEMAKE_THREADS="$THREADS" \
    --env SNEMAKE_MEM_MB="$MEM_MB" \
    --env SNEMAKE_LOG="$LOG_FILE" \
    --env OMP_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env OPENBLAS_NUM_THREADS=1 \
    $ENV_VARS \
    $INPUT_FILES \
    $OUTPUT_FILES \
    --log "$LOG_FILE" \
    --machine-type n2-standard-32 \
    --boot-disk-size 100 \
    --disk-size 500

echo "Submitted job: $JOBNAME"
```

### Running Snakemake with dsub

```bash
# Activate environment
conda activate snakemake-snaphic

# Run with dsub profile
snakemake \
    --profile profiles/dsub \
    --configfile config.yaml \
    --jobs 100 \
    --latency-wait 60 \
    --keep-going \
    --rerun-incomplete

# Dry run first
snakemake -n -p --profile profiles/dsub --configfile config.yaml

# Lint check
snakemake --lint --profile profiles/dsub
```

---

## 5. dsub Commands Reference

### dsub — Submit Jobs

```bash
# Basic job
dsub --script script.sh --name my-job --queue x86 -aa

# With input/output
dsub \
    --script process.sh \
    --name data-process \
    --queue x86 \
    -aa \
    --input /data/input1.txt \
    --input /data/input2.txt \
    --output /results/output.txt

# With environment variables
dsub \
    --script train.py \
    --name model-train \
    --queue x86 \
    -aa \
    --env MODEL_TYPE=transformer \
    --env EPOCHS=100 \
    --env BATCH_SIZE=32

# With resource specifications
dsub \
    --script compute.sh \
    --name compute-job \
    --queue x86 \
    -aa \
    --env THREADS=8 \
    --env MEMORY_GB=32 \
    --machine-type n2-standard-32

# Dry run
dsub --dry-run --script script.sh --name test-job --queue x86 -aa
```

### djob — Job Management

```bash
# List all jobs for user
djob -u $USER

# List jobs with status
djob -u $USER -s

# Show job details
djob -j JOB_ID

# Show job history
djob -u $USER -H

# Filter by status
djob -u $USER -s | grep -E "RUNNING|PENDING|FAILURE"
```

### dkill — Job Cancellation

```bash
# Kill single job
dkill -j JOB_ID

# Kill multiple jobs
dkill -j JOB_ID_1 -j JOB_ID_2

# Kill all jobs matching pattern
dkill -u $USER -t pattern

# Kill all running jobs
dkill -u $USER -s RUNNING
```

### dpeek — Log Retrieval

```bash
# Retrieve stdout
dpeek -j JOB_ID

# Retrieve stderr
dpeek -j JOB_ID -e

# Tail last 100 lines
dpeek -j JOB_ID | tail -100

# Save to file
dpeek -j JOB_ID > job_output.txt
```

---

## 6. Job Monitoring

### Status Check Workflow

```bash
# 1. List your jobs
djob -u $USER -s

# Sample output:
# JOB_ID              NAME                    STATUS
# 12345               smk-validate-abc        RUNNING
# 12346               smk-extract-def         PENDING
# 12347               smk-windows-ghi         RUNNING

# 2. Get details for specific job
djob -j 12345

# 3. Check logs
dpeek -j 12345

# 4. If failed, get error
dpeek -j 12345 -e
```

### Common Status Values

| Status | Meaning |
|--------|---------|
| `PENDING` | Waiting for resources |
| `RUNNING` | Currently executing |
| `SUCCESS` | Completed successfully |
| `FAILURE` | Completed with errors |
| `CANCELED` | User or admin canceled |

### Batch Monitoring Script

```bash
#!/bin/bash
# monitor_jobs.sh — Monitor running jobs

USER="${1:-$USER}"
INTERVAL="${2:-30}"

echo "Monitoring jobs for $USER (Ctrl+C to exit)"
echo "=========================================="
echo ""

while true; do
    clear
    echo "Updated: $(date)"
    echo ""
    djob -u $USER -s
    echo ""
    echo "Running: $(djob -u $USER -s | grep RUNNING | wc -l)"
    echo "Pending: $(djob -u $USER -s | grep PENDING | wc -l)"
    echo "Failed:  $(djob -u $USER -s | grep FAILURE | wc -l)"
    sleep "$INTERVAL"
done
```

---

## 7. Log Retrieval

### Manual Log Access

```bash
# Get stdout log
dpeek -j JOB_ID

# Get stderr log
dpeek -j JOB_ID -e

# Get both combined
(dpeek -j JOB_ID 2>&1) | head -100

# Search in logs
dpeek -j JOB_ID | grep -i error

# Count errors
dpeek -j JOB_ID | grep -ci error
```

### Log File Locations

By default, logs are written to:
- stdout: Specified by `-o` flag or `%o.log`
- stderr: Specified by `-e` flag or `%e.log`

```bash
# Custom log locations
dsub --script job.sh --name my-job \
    -o /logs/my-job.stdout \
    -e /logs/my-job.stderr
```

---

## 8. Failure Handling

### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | General error |
| 2 | Misuse of command |
| 126 | Command not executable |
| 127 | Command not found |
| 130 | Terminated by signal (Ctrl+C) |

### Failure Investigation

```bash
# 1. Check job status
djob -j JOB_ID

# 2. Get error log
dpeek -j JOB_ID -e

# 3. Check common issues
dpeek -j JOB_ID | grep -iE "error|exception|traceback|failed"

# 4. Check resource usage
dpeek -j JOB_ID | grep -iE "memory|time|cpu"
```

### Resume After Failure

```bash
# Clean up failed jobs
dkill -u $USER -s FAILURE
dkill -u $USER -s CANCELED

# Resume Snakemake
snakemake \
    --profile profiles/dsub \
    --configfile config.yaml \
    --jobs 100 \
    --rerun-incomplete \
    --keep-going

# Or clean and restart
snakemake \
    --profile profiles/dsub \
    --configfile config.yaml \
    --unlock
snakemake \
    --profile profiles/dsub \
    --configfile config.yaml \
    --jobs 100 \
    --rerun-incomplete
```

### Automatic Retry Script

```bash
#!/bin/bash
# retry_failed.sh — Retry failed Snakemake jobs

MAX_RETRIES=3
RETRY_COUNT=0

while [ $RETRY_COUNT -lt $MAX_RETRIES ]; do
    echo "Attempt $((RETRY_COUNT + 1)) of $MAX_RETRIES"
    
    snakemake \
        --profile profiles/dsub \
        --configfile config.yaml \
        --jobs 100 \
        --rerun-incomplete \
        --keep-going
    
    EXIT_CODE=$?
    
    if [ $EXIT_CODE -eq 0 ]; then
        echo "Success!"
        exit 0
    fi
    
    RETRY_COUNT=$((RETRY_COUNT + 1))
    
    if [ $RETRY_COUNT -lt $MAX_RETRIES ]; then
        echo "Failed. Waiting 60 seconds before retry..."
        sleep 60
    fi
done

echo "All retries exhausted"
exit 1
```

---

## 9. Resource Specification

### CPU-Based Model

dsub uses a CPU-based resource model:

```bash
# Request 4 CPUs
dsub --script job.sh --name my-job -R "cpu=4"

# Request 8 CPUs
dsub --script job.sh --name my-job -R "cpu=8"

# Auto-calculate from Snakemake threads
CPUS=$(( (THREADS + 1) / 2 ))
```

### Memory Proportional to CPU

| CPUs | Recommended Memory |
|------|-------------------|
| 1 | 4 GB |
| 2 | 8 GB |
| 4 | 16 GB |
| 8 | 32 GB |
| 16 | 64 GB |

```bash
# For 8 threads with 32GB memory
dsub --script job.sh --name my-job \
    -R "cpu=8" \
    -R "mem=32000"
```

### Machine Types

```bash
# Standard machine types
--machine-type n2-standard-32    # 32 CPUs, 128 GB RAM
--machine-type n2-standard-16    # 16 CPUs, 64 GB RAM
--machine-type n2-standard-8     # 8 CPUs, 32 GB RAM
--machine-type n2-standard-4     # 4 CPUs, 16 GB RAM
```

---

## 10. Queue Selection

### x86 Queue

The primary queue for x86_64 architecture:

```bash
# Basic x86 queue
dsub --script job.sh --name my-job -q x86

# With any-arch flag (faster scheduling)
dsub --script job.sh --name my-job -q x86 -aa
```

### Priority Levels

```bash
# Standard priority (default)
dsub --script job.sh --name my-job -q x86

# High priority (if authorized)
dsub --script job.sh --name my-job -q x86 --priority high
```

---

## 11. Common Issues

### Path Resolution

```bash
# PROBLEM: Relative paths don't work in dsub
--input data/sample.txt    # WRONG - relative path

# SOLUTION: Use absolute paths
--input /full/path/to/data/sample.txt

# SOLUTION: Pass via environment variable
--env DATA_DIR=/full/path/to/data \
--input %(DATA_DIR)s/sample.txt
```

### Environment Variables

```bash
# PROBLEM: Variables not passed correctly
--env PATH=/usr/bin

# SOLUTION: Use explicit key=value
--env MY_VAR=value

# PROBLEM: Quoting issues
--env VAR="value with spaces"

# SOLUTION: Escape inner quotes
--env 'VAR=value with spaces'
```

### Nested Parallelism

```bash
# PROBLEM: Jobs use too many threads
OMP_NUM_THREADS not set  # Uses all available cores

# SOLUTION: Always disable in wrapper
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# SOLUTION: Verify in script
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
```

### File Permissions

```bash
# PROBLEM: Cannot write to output
Permission denied: /results/output.txt

# SOLUTION: Ensure parent directory exists
mkdir -p /results

# SOLUTION: Use user-specific paths
--env OUTDIR=/home/$USER/results

# SOLUTION: Set proper permissions
chmod 755 /results
```

### Job Name Conflicts

```bash
# PROBLEM: Duplicate job names cause issues
--name snakemake-job

# SOLUTION: Include unique identifier
--name "smk-$(date +%s)-$$"
--name "smk-$RULE-$JOBID"
```

---

## Summary Checklist

- [ ] All wrapper scripts disable nested parallelism
- [ ] All scripts have proper `#DSUB` headers
- [ ] All paths are absolute
- [ ] All environment variables are explicitly passed
- [ ] Memory is proportional to CPU count
- [ ] Job names include unique identifiers
- [ ] Logs are captured to accessible locations
- [ ] Failed jobs are investigated before retry
- [ ] `snakemake -n -p` passes with dsub profile
