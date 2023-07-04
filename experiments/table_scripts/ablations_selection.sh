#!/bin/bash
export TOL=1e-4
export TAU=1e-3
python experiments/python/comparisons.py --label ablation_selection_baseline --in_dir ../../data/ablations --nosig21 --tau $TAU --tolerance $TOL
python experiments/python/comparisons.py --label ablation_selection_2closest --in_dir ../../data/ablations --nosig21 --tau $TAU --tolerance $TOL --ablation --ablation_n 2
python experiments/python/comparisons.py --label ablation_selection_3closest --in_dir ../../data/ablations --nosig21 --tau $TAU --tolerance $TOL --ablation --ablation_n 3
python experiments/python/comparisons.py --label ablation_selection_3random --in_dir ../../data/ablations --nosig21 --tau $TAU --tolerance $TOL --ablation --ablation_n 3 --ablation_random
python experiments/python/comparisons.py --label ablation_selection_4closest --in_dir ../../data/ablations --nosig21 --tau $TAU --tolerance $TOL --ablation --ablation_n 4
python experiments/python/comparisons.py --label ablation_selection_alltri --in_dir ../../data/ablations --nested --nosig21 --tau $TAU --tolerance $TOL --all_triangles 