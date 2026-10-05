# Nutures
# time: 2026-01-07
# Purpose: Central output directory for Python3 3DNetMod; overridable via PY3_OUTDIR.

import os
import hashlib


def consensus_filtered_filename(region, tag, gamma):
    digest = hashlib.sha256(
        "{}|{}|{}".format(region, tag, gamma).encode()
    ).hexdigest()[:16]
    return "Consensus_filtered_{}_{}_gamma_{}.csv".format(
        region, digest, gamma
    )


# When PY3_OUTDIR is set (e.g. by run_python3_pipeline.sh), all output goes there.
# Default "output" keeps legacy behavior (repo-root/output when run from repo root).
OUTDIR = os.environ.get("PY3_OUTDIR", "output")
