# Exposure and sensitivity to biological invasions
<!-- badges: start -->
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://choosealicense.com/licenses/mit/)
<!-- badges: end -->

Scripts and data for reproducing the results obtained by Marino et al. in the paper **Exposure and sensitivity of terrestrial vertebrates to biological invasions worldwide**.


## Description of the Data
All information related to data are detailed in [`Data/README.md`](https://github.com/claramarino/BIVA_framework/main/Data/README.md).

## How to use the scripts

### 1. Data preparation

As all raw data cannot be uploaded on github, the scripts for computing the metrics are based on already derived and clean data. 
The scripts [`BIVA_10_Exposure_metrics.R`](https://github.com/claramarino/BIVA_framework/main/BIVA_10_Exposure_metrics.R) and [`BIVA_13_Sensitivity_metrics.R`](https://github.com/claramarino/BIVA_framework/main/BIVA_13_Sensitivity_metrics.R) are used for computing the exposure and sensitivity metrics, respectively, using data from the [`Data/derived-data/`](https://github.com/claramarino/BIVA_framework/main/Data/derived-data/) folder. The script [`BIVA_14_Completeness_Sensitivity.R`](https://github.com/claramarino/BIVA_framework/main/BIVA_14_Completeness_Sensitivity.R) is used for deriving the completeness for sensitivity. All the outputs of those 

### 2. Analyses and visualization

The scripts from BIVA_20 to BIVA_23 reproduce all the analyses conducted in the paper. All the data used in these analyses are detailed in [`Data/README.md`](https://github.com/claramarino/BIVA_framework/main/Data/README.md). 

## Citation

Please use the following citation:

> Marino C (2024) Data and codes for Exposure and sensitivity of terrestrial vertebrates to biological invasions worldwide. URL: https://github.com/claramarino/BIVA_framework.
