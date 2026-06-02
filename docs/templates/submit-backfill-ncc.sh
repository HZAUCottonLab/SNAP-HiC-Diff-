#!/bin/bash
#DSUB -q arm
#DSUB -aa
#DSUB -n backfill-ncc
#DSUB -e %J.err
#DSUB -o %J.out
#DSUB -R "cpu=10"

set -euo pipefail

cd /share/home/xhhuang24/project/PanTE/background_pool
chmod 777 scripts/backfill_ncc_and_rebuild_background.sh
bash scripts/backfill_ncc_and_rebuild_background.sh
