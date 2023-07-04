# Gravo MG: Graph Voronoi Multigrid
[[Paper]](https://graphics.tudelft.nl/~klaus/papers/Gravo_MG.pdf) [[Project page]](https://rubenwiersma.nl/gravomg)
![](https://rubenwiersma.nl/assets/img/publications/gravomg/teaser_gravomg.png)

Repository for **SIGGRAPH 2023** paper **"A Fast Geometric Multigrid Method for Curved Surfaces"** by Ruben Wiersma, [Ahmad Nasikun](https://github.com/a-nasikun) (equal contribution); Elmar Eisemann; Klaus Hildebrandt.

This repository contains code for Gravo MG: an approach to solve linear systems on curved surfaces quickly. We make use of a geometric multigrid method that solves the system by iterating on several levels in a hierarchy. Gravo MG is a fast way to compute this hierarchy that keeps solving times low.

## Coming soon: Standalone C++ and Python libraries
This repository contains the code to replicate the SIGGRAPH 2023 paper, including ablations and comparisons. That means the repository contains many components that are not strictly relevant for using Gravo MG. We plan to release a smaller, more lightweight library in the coming weeks.

# Replicating our results
The tables in our paper are created using the scripts in `experiments/table_scripts`. These scripts require you to build and install a Python package, using pip.

## Setting up the environment
0. Clone this repository, including the submodules (required to pull in Pybind11), and change into the directory:
```
$ git clone --recurse-submodules https://github.com/rubenwiersma/gravo-mg.git
$ cd gravo-mg
```
1. Create a Conda environment with the necessary requirements and activate environment:
```
$ conda env create -f environment.yml
$ conda activate gravomg
```
2. Install the `gravomg` Python package:
```
$ pip install -e ./gravomg_bindings
```

This builds the Gravo MG C++ library and wraps it in a Python binding.

## Running the experiments
First, [download the data](https://surfdrive.surf.nl/files/index.php/s/gOAGyWdSVJVPrBb), unzip, and place it in the `gravo-mg` folder. The result should be a folder `data` in the root of this repository.

You can run each experiment from the experiments folder, e.g.:
```
$ sh experiments/table_scripts/comparison_poisson.sh
```

The output is written in `out/timing` and a formatted LaTeX table in `out/latex`.

## Pardiso
To compute the timings for Pardiso, make sure that you have [installed Pardiso on your machine](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html). Our library automatically looks for Pardiso in `/usr/include/mkl` and enables it for use when compiling. You can change the search path in `gravomg/CMakeLists.txt`.

If you've already installed the `gravomg` package before installing OneMKL, make sure that you reinstall `gravomg`.

# Acknowledgements

## Code
Our code uses Derek Liu and colleagues' code for [Surface Multigrid via Intrinsic Prolongation](https://github.com/HTDerekLiu/surface_multigrid_code) to replicate their results and as a basis for the main solver routine.

## Data
[Please find the source and author of each mesh used in our experiments in the linked sheet.](https://docs.google.com/spreadsheets/d/1s5ogLIqmCHthTtyOcgc1SADOlBfXtgP1vG7Qh-oaVdk/edit?usp=sharing)