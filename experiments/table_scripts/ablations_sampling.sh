#!/bin/bash
export TOL=1e-4
python experiments/python/comparisons.py --label ablation_sampling_baseline --in_dir ../../data/ablations --nosig21 --tolerance $TOL
python experiments/python/comparisons.py --label ablation_sampling_random --in_dir ../../data/ablations --nosig21 --tolerance $TOL --no_names --sampling random
python experiments/python/comparisons.py --label ablation_sampling_pds --in_dir ../../data/ablations --nosig21 --tolerance $TOL --no_names --sampling poissondisk
python experiments/python/comparisons.py --label ablation_sampling_mis --in_dir ../../data/ablations --nosig21 --tolerance $TOL --no_names --sampling mis
python experiments/python/comparisons.py --label ablation_sampling_fps --in_dir ../../data/ablations --nosig21 --tolerance $TOL --no_names --sampling fps