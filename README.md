[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15883770.svg)](https://doi.org/10.5281/zenodo.15883770)


# Code repository for _Model-data synthesis of benthic isotopes suggests a warmer Miocene Climatic Optimum_

This repository gathers notebooks for the manuscript _Model-data synthesis of benthic isotopes suggests a warmer Miocene Climatic Optimum_.

The code is tested with Python 3.13, and the package [`x4c`](https://ncar.github.io/x4c/) is required to perform essential analysis and the corresponding visualization.

## Repository Structure
- `notebooks`: the directory that includes the Jupyter notebooks
    - `BayesianEst.py`: the probabilistic proxy system model (PSM) proposed in this study, which is used to generate the `df_CO2pdf.csv` under the `data` directory
    - `SIFig.03.ipynb`: the notebook that performs probabilistic proxy system modeling using `BayesianEst.py` and generates SI Fig. 3 in the Supplementary Information
    - `Fig.01.ipynb`: the notebook that performs analysis and generates Fig. 1 in the main text
    - `Fig.02.ipynb`: the notebook that performs analysis and generates Fig. 2 in the main text
    - `Fig.03.ipynb`: the notebook that performs analysis and generates Fig. 3 in the main text
    - `Fig.04.ipynb`: the notebook that performs analysis and generates Fig. 4 in the main text
    - `Fig.05.ipynb`: the notebook that performs analysis and generates Fig. 5 in the main text
- `data`: the directory that includes the auxiliary data for the analysis and visualization (note that the raw iCESM1.3 output are stored on NSF NCAR's Glade system)
- `figs`: the folder that includes the figures in the main text


## How to cite this repository
This repository can be cited with DOI: [10.5281/zenodo.15883770](https://doi.org/10.5281/zenodo.15883770)

