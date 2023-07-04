#!/bin/bash
export TOL=1e-4
export ETA=1e-6 # Sets the regularizer in the Poisson problem
python experiments/python/comparisons.py --label noef_poisson_all --poisson --in_dir ../../data/final --sig06 --amg --direct --tolerance $TOL --tau $ETA
python experiments/python/comparisons.py --label noef_poisson_nonmanifold --poisson --in_dir ../../data/nonmanifold --nonmanifold --robust --sig06 --amg --direct --tolerance $TOL --tau $ETA
python experiments/python/comparisons.py --label noef_poisson_pointcloud --poisson --in_dir ../../data/pointcloud --pointcloud --nested --sig06 --amg --direct --tolerance $TOL --tau $ETA