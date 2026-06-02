#!/bin/bash
#DSUB -q arm
#DSUB -aa
#DSUB -n knn-viz
#DSUB -e %J.err
#DSUB -o %J.out
#DSUB -R "cpu=10"

cd /share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff
python3 jupyter/knn_viz.py --comparison H1398_vs_L134
