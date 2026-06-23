[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15883769.svg)](https://doi.org/10.5281/zenodo.15883769)


# Code repository for _Model-data synthesis of benthic isotopes suggests a warmer Miocene Climatic Optimum_

This repository gathers notebooks for the manuscript _Model-data synthesis of benthic isotopes suggests a warmer Miocene Climatic Optimum_.

The code is tested with Python 3.13, and the package [`x4c`](https://ncar.github.io/x4c/) is required to perform essential analysis and the corresponding visualization.

## Repository Structure
- `notebooks`: the directory that includes the Jupyter notebooks and the utility scripts
    - [`BayesianEst.py`](https://github.com/fzhu2e/paper-MCO_iCESM/blob/main/notebooks/BayesianEst.py): the probabilistic proxy system model (PSM) proposed in this study, which is used to generate the `df_CO2pdf.csv` under the `data` directory
    - [`export_MLE.ipynb`](https://github.com/fzhu2e/paper-MCO_iCESM/blob/main/notebooks/export_MLE.ipynb): the notebook that performs interpolation of the raw iCESM1.3 outputs and exports variables for the MLE case
    - [`Fig.01.ipynb`](https://github.com/fzhu2e/paper-MCO_iCESM/blob/main/notebooks/Fig.01.ipynb): the notebook that performs analysis and generates Fig. 1 in the main text
    - [`Fig.02.ipynb`](https://github.com/fzhu2e/paper-MCO_iCESM/blob/main/notebooks/Fig.02.ipynb): the notebook that performs analysis and generates Fig. 2 in the main text
    - [`Fig.03.ipynb`](https://github.com/fzhu2e/paper-MCO_iCESM/blob/main/notebooks/Fig.03.ipynb): the notebook that performs analysis and generates Fig. 3 in the main text
    - [`Fig.04.ipynb`](https://github.com/fzhu2e/paper-MCO_iCESM/blob/main/notebooks/Fig.04.ipynb): the notebook that performs analysis and generates Fig. 4 in the main text
    - [`Fig.05.ipynb`](https://github.com/fzhu2e/paper-MCO_iCESM/blob/main/notebooks/Fig.05.ipynb): the notebook that performs analysis and generates Fig. 5 in the main text
    - [`SIFig.03.ipynb`](https://github.com/fzhu2e/paper-MCO_iCESM/blob/main/notebooks/SIFig.03.ipynb): the notebook that performs probabilistic proxy system modeling using `BayesianEst.py` and generates SI Fig. 3 in the Supplementary Information
    - [`SIFig.05.ipynb`](https://github.com/fzhu2e/paper-MCO_iCESM/blob/main/notebooks/SIFig.05.ipynb): the notebook that performs the sensitivity test on site selection and generates SI Fig. 5 in the Supplementary Information
- `data`: the directory that includes the auxiliary data for the analysis and visualization (note that the raw iCESM1.3 output are stored on NSF NCAR's Glade system)
    - [`iCESM1.3_1.5x_TEMP_clim_eq.nc`](https://github.com/fzhu2e/paper-MCO_iCESM/raw/refs/heads/main/data/iCESM1.3_1.5x_TEMP_clim_eq.nc): the equilibrium state of TEMP (ocean temperature) of the 1.5xCO2 case
    - [`iCESM1.3_1.5x_SALT_clim_eq.nc`](https://github.com/fzhu2e/paper-MCO_iCESM/raw/refs/heads/main/data/iCESM1.3_1.5x_SALT_clim_eq.nc): the equilibrium state of SALT (ocean salinity) of the 1.5xCO2 case
    - [`iCESM1.3_1.5x_d18Osw_clim_eq.nc`](https://github.com/fzhu2e/paper-MCO_iCESM/raw/refs/heads/main/data/iCESM1.3_1.5x_d18Osw_clim_eq.nc): the equilibrium state of d18Osw (sea-water d18O) of the 1.5xCO2 case; the global volume mean is set to -0.55 permil
    - [`iCESM1.3_3x_TEMP_clim_eq.nc`](https://github.com/fzhu2e/paper-MCO_iCESM/raw/refs/heads/main/data/iCESM1.3_3x_TEMP_clim_eq.nc): the equilibrium state of TEMP (ocean temperature) of the 3xCO2 case
    - [`iCESM1.3_3x_SALT_clim_eq.nc`](https://github.com/fzhu2e/paper-MCO_iCESM/raw/refs/heads/main/data/iCESM1.3_3x_SALT_clim_eq.nc): the equilibrium state of SALT (ocean salinity) of the 3xCO2 case
    - [`iCESM1.3_3x_d18Osw_clim_eq.nc`](https://github.com/fzhu2e/paper-MCO_iCESM/raw/refs/heads/main/data/iCESM1.3_3x_d18Osw_clim_eq.nc): the equilibrium state of d18Osw (sea-water d18O) of the 3xCO2 case; the global volume mean is set to -0.55 permil
    - [`iCESM1.3_MLE_TEMP_clim_eq.nc`](https://github.com/fzhu2e/paper-MCO_iCESM/raw/refs/heads/main/data/iCESM1.3_MLE_TEMP_clim_eq.nc): the equilibrium state of TEMP (ocean temperature) of the MLE case
    - [`iCESM1.3_MLE_SALT_clim_eq.nc`](https://github.com/fzhu2e/paper-MCO_iCESM/raw/refs/heads/main/data/iCESM1.3_MLE_SALT_clim_eq.nc): the equilibrium state of SALT (ocean salinity) of the MLE case
    - [`iCESM1.3_MLE_d18Osw_clim_eq.nc`](https://github.com/fzhu2e/paper-MCO_iCESM/raw/refs/heads/main/data/iCESM1.3_MLE_d18Osw_clim_eq.nc): the equilibrium state of d18Osw (sea-water d18O) of the MLE case
    - [`iCESM1.3_MLE_SST_climo.nc`](https://github.com/fzhu2e/paper-MCO_iCESM/raw/refs/heads/main/data/iCESM1.3_MLE_SST_climo.nc): the monthly climatology of SST (sea-water temperature) of the MLE case
    - [`iCESM1.3_MLE_aice_climo.nc`](https://github.com/fzhu2e/paper-MCO_iCESM/raw/refs/heads/main/data/iCESM1.3_MLE_aice_climo.nc): the monthly climatology of aice (sea-ice area) of the MLE case
    - [`iCESM1.3_MLE_R18O_climo.nc`](https://github.com/fzhu2e/paper-MCO_iCESM/raw/refs/heads/main/data/iCESM1.3_MLE_R18O_climo.nc): the monthly climatology of d18Osw (sea-water d18O) of the MLE case; the global volume mean is set to -0.55 permil
    - [`iCESM1.3_MLE_RHDO_climo.nc`](https://github.com/fzhu2e/paper-MCO_iCESM/raw/refs/heads/main/data/iCESM1.3_MLE_RHDO_climo.nc): the monthly climatology of dDsw (sea-water dD) of the MLE case; the global volume mean is reset accordingly
- `CESM_configs`: the directory that includes the CESM configuration related files, including the namelists and the necessary MCO topography and bathymetry files, etc.
- `pygplates_data`: the directory that includes the static files for rotating present-day locations to paleo-locations using [PyGPlates](https://github.com/GPlates/GPlates)


## How to cite this repository
This repository can be cited with DOI: [10.5281/zenodo.15883769](https://doi.org/10.5281/zenodo.15883769)

