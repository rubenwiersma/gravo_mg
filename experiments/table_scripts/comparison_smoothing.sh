#!/bin/bash
export TOL=1e-4
export ALPHA=1e-3 # Sets the smoothing coefficient
python experiments/python/comparisons.py --label noef_smoothing_all --in_dir ../../data/final --sig06 --amg --direct --tolerance $TOL --tau $ALPHA
python experiments/python/comparisons.py --label noef_smoothing_nonmanifold --in_dir ../../data/nonmanifold --nonmanifold --robust --sig06 --amg --direct --tolerance $TOL --tau $ALPHA
python experiments/python/comparisons.py --label noef_smoothing_pointcloud --in_dir ../../data/pointcloud --pointcloud --nested --sig06 --amg --direct --tolerance $TOL --tau $ALPHA