# Gravo MG: Graph Voronoi Multigrid
[[Paper]](https://graphics.tudelft.nl/~klaus/papers/Gravo_MG.pdf) [[Project page]](https://rubenwiersma.nl/gravomg) [[C++ library]](https://github.com/rubenwiersma/gravo_mg_cpp) [[Python Bindings]](https://github.com/rubenwiersma/gravo_mg_python)

![](https://rubenwiersma.nl/assets/img/publications/gravomg/teaser_gravomg.png)

Repository for **"A Fast Geometric Multigrid Method for Curved Surfaces"**, published at **SIGGRAPH 2023**
<br />
by Ruben Wiersma, [Ahmad Nasikun](https://github.com/a-nasikun) (equal contribution); Elmar Eisemann; Klaus Hildebrandt.

If you need to solve linear systems on meshes or point clouds with >50.000 vertices, Gravo MG is what you need. Gravo MG is a fast geometric multigrid method that quickly computes a hierarchy used in an iterative multigrid solver. The use of graph Voronoi cells enables fast convergence, while allowing for fast construction.

## Standalone C++ and Python libraries
This repository contains many components that are not strictly relevant for using Gravo MG (comparisons, related work, etc.). It is recommended to use the following packages, instead of this repository:
- [Gravo MG C++ library](https://github.com/rubenwiersma/gravo_mg_cpp)
- [Gravo MG Python bindings](https://github.com/rubenwiersma/gravo_mg_python)

For Python, that means you can use
```
pip install gravomg
```
and start solving. More info can be found in the [Python binding repository](https://github.com/rubenwiersma/gravo_mg_python).

## Replicating our results
The tables in our paper are created using the scripts in `experiments/table_scripts`. These scripts require you to build and install a Python package using pip.

### Setting up the environment
1. Clone this repository, including the submodules (required to pull in Pybind11), and change into the directory:
```bash
git clone --recursive https://github.com/rubenwiersma/gravo_mg.git
cd gravo_mg
```

If you forgot to use the `--recursive` flag, you can pull the required submodules with
```bash
git submodule update --init --recursive
```

2. Create a Conda environment with the necessary requirements and activate environment:
```bash
conda env create -f environment.yml
conda activate gravomg
```
3. Install the `gravomg` Python package:
```bash
pip install -e ./gravomg_bindings
```

This builds the Gravo MG C++ library and wraps it in a Python binding. The `-e` flag installs the package in the current folder, rather than your default Conda folder.

### Running the experiments
First, [download the data](https://surfdrive.surf.nl/files/index.php/s/gOAGyWdSVJVPrBb), unzip, and place it in the `gravo_mg` folder. The result should be a folder `data` in the root of this repository.

You can run each experiment from the experiments folder, e.g.:
```bash
sh experiments/table_scripts/comparison_poisson.sh
```

The output is written in `out/timing` and a formatted LaTeX table in `out/latex`.

### Pardiso
To compute the timings for Pardiso, make sure that you have [installed Intel's OneMKL on your machine](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html). Our library automatically looks for OneMKL in `/usr/include/mkl` and enables it for use when compiling. You can change the search path in `gravomg/CMakeLists.txt`.

If you've already installed the `gravomg` package before installing OneMKL, make sure that you reinstall `gravomg`.

## Acknowledgements

### Code
Our code uses Derek Liu and colleagues' code for [Surface Multigrid via Intrinsic Prolongation](https://github.com/HTDerekLiu/surface_multigrid_code) to replicate their results and as a basis for the main solver routine.

### Data
You can download the data used in the paper here: [download the data](https://surfdrive.surf.nl/files/index.php/s/gOAGyWdSVJVPrBb).

[Please find the source and author of each mesh used in our experiments in the linked sheet.](https://docs.google.com/spreadsheets/d/1s5ogLIqmCHthTtyOcgc1SADOlBfXtgP1vG7Qh-oaVdk/edit?usp=sharing)

## Citations
Please cite our paper if this code contributes to an academic publication:

```bib
@Article{WiersmaNasikun2023GravoMG,
author = {Ruben Wiersma, Ahmad Nasikun, Elmar Eisemann, Klaus Hildebrandt},
journal = {SIGGRAPH 2023},
title = {A Fast Geometric Multigrid Method for Curved Surfaces},
year = {2023},
month = jul,
number = {4},
volume = {41},
doi = {10.1145/3588432.3591502},
publisher = {ACM}
}
```
