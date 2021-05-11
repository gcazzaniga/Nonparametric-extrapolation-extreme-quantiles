# Nonparametric-extrapolation-extreme-quantiles

This repository contains the code to compute non parametric extrapolation of extreme quantiles, from synthetic dataset. In the code are implemented three nonparametric techniques for extrapolation: Hutson (2002), Kernel Density Estimation, and Scholz (2005). For any question or comment concerning the code please contact Fabiola Banfi (fabiola.banfi@polimi.it) or Greta Cazzaniga (greta.cazzaniga@polimi.it). 

The code for generation of synthetic data and extrapolation of extreme quantiles is written in R. Packages required to run the code include `climod`, `kdensity`, `rstudioapi`, `rtop`, and `reshape2`.
The code, which plot figures is instead written in Python 3 and dependencies required are `pyreadr`, `numpy`, `pandas`, and `matplotlib`. Moreover an installation of jupyter is needed to run the *.ipynb notebook.

## Repository structure:
- `0_RunCode.R` script contains 
- `1_main.R` script
- `2_functions.R` script
- `CreateFigures.ipynb`
