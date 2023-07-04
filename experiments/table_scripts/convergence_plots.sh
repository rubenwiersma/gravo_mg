#!/bin/bash
export TOL=1e-12
export ALPHA=1e-3 # Sets the smoothing coefficient
python experiments/python/comparisons.py --label convergence_smoothing --in_dir ../../data/main --sig06 --tolerance $TOL --tau $ALPHA