#!/bin/bash
export TOL=1e-4
export TAU=1e-3
python experiments/python/comparisons.py --label ablation_weighting_baseline --in_dir ../../data/ablations --nosig21 --tau $TAU --tolerance $TOL
python experiments/python/comparisons.py --label ablation_weighting_uniform --in_dir ../../data/ablations --nosig21 --tau $TAU --tolerance $TOL --no_names --weighting uniform
python experiments/python/comparisons.py --label ablation_weighting_invdist --in_dir ../../data/ablations --nosig21 --tau $TAU --tolerance $TOL --no_names --weighting invdist
python experiments/python/comparisons.py --label ablation_weighting_nested --in_dir ../../data/ablations --nosig21 --tau $TAU --tolerance $TOL --no_names --nested