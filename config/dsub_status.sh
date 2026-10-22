#!/usr/bin/env bash
# Translate a Donau job state for Snakemake's cluster-generic executor.

set -euo pipefail

job_id="${1:?missing job id}"
state="$(djob "${job_id}" | awk 'NR == 2 {print $3}')"

case "${state}" in
    SUCCEEDED) printf 'success\n' ;;
    FAILED|CANCELED|CANCELLED|ENDED) printf 'failed\n' ;;
    *) printf 'running\n' ;;
esac
