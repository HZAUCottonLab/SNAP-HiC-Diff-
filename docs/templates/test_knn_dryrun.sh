#!/bin/bash
#DSUB -q arm
#DSUB -aa
#DSUB -n test-knn-dryrun
#DSUB -e %J.err
#DSUB -o %J.out
#DSUB -R "cpu=10"

cd /share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff
source ~/.bashrc
conda activate snap-hic-diff
snakemake --dryrun 2>&1 | head -50
