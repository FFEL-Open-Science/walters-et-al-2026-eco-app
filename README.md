# Description

This repository contains data and codes associated with the article [**Disentangling the contributions of density dependence and independence to population growth rates**]()

**IMPORTANT:** You will need to set up Git LFS when cloning this repository.
This is required so you can download output files (`out_fit_dat.rds`), which can be up to about 2GB in size.
Once you have setup Git LFS (see instructions [here](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage)), issue the command `git lfs pull` in an RStudio terminal, so the output file is downloaded, rather than just the pointer file.

[![DOI](https://zenodo.org/badge/1126025607.svg)](https://doi.org/10.5281/zenodo.18140770)

## Data and file structure

The following data files were used in the analyses.
Missing data are identified as NA.

### Raw data

Data files show in this section are located in the folder `data/raw`.

`age_length.csv`: Data file containing length-at age data for Stellako River rainbow trout collected in 1988 and 2017.
The following variables are included: - *year*: Year of sample.
- *date*: Date the sample was taken.
- *sample_id*: Unique ID assigned to sample.
- *fl*: Fork length of the fish at capture in mm. - *sex*: Sex of the fish.
- *maturity*: Maturity of the fish at capture (mature, immature, spent) - *age_capture*: Age of the fish at capture.
- *1 to 9*: Column number denotes the age and the cell values are the estimated length.

`conditions.csv`: Data file containing visibility conditions for each snorkel survey from 1988 to 2022.
Some values are missing from years they were not collected.
The following variables are included: - *date*: Date the sample was taken.
- *visibility*: Measured visibility in m.
- *discharge*: Stellako River discharge in m3/s.
- *temperature*: Stellako River water temperature in °C.\
- *stage*: Water depth in m.
- *comment*: Relevant comments of river conditions.

`lane_reach.csv`: Date file describing the lanes and reach surveyed for each snorkel survey day from 1988 to 2022.
The following variables are included: - *year*: Year of survey.
- *date*: Date the survey was completed.
- *survey_day*: The relative day the reach was surveyed in a year.
- *nlanes*: The number of lanes surveyed.
- *reach*: The reach surveyed.

`rb_counts.csv`: Data file with the snorkel survey counts.
The following variables are included: - *year*: Year of survey.
- *date*: Date the survey was completed.
- *survey_day*: The relative day the reach was surveyed in a year.
- *size_class*: Rainbow trout size class.
- *count*: The number of rainbow trout in a size class counted during the snorkel survey day.

`sk_nuseds.csv`: Data file with the annual abundance of spawning sockeye salmon returns to the Stellako River.
Created from data downloaded from [DFO NuSEDS](https://open.canada.ca/data/en/dataset/c48669a3-045b-400d-b730-48aafe8c5ee6). The following variables are included: - *year*: Year of sample.
- *abundance*: The total abundance of spawning sockeye salmon returns.

`wsc_water_temperature_2011-2024.csv`: The Water Survey of Canada (WSC) data on water temperature in the Stellako River, station no. 08JB002, obtained via email request.
- *dttm*: The date and time of the sample.
- *temperature*: Stellako River water temperature in °C.

`era5_air_temperature_1980-1994.nc`: Data file containing 1980-1994 daily mean air temperature data in the region containing the Stellako River, which was defined by longitudes 125.0°W – 124.9°W and latitudes 54.0°N – 54.1°N.
The data were downloaded from the Copernicus Climate Data Store in the NetCDF format (.nc) and processed into a data frame using R code `codes/process/02-process_environmental_data.R`.

`era5_air_temperature_1995-2008.nc`: Data file containing 1995-2008 daily mean air temperature data in the region containing the Stellako River, which was defined by longitudes 125.0°W – 124.9°W and latitudes 54.0°N – 54.1°N.
The data were downloaded from the Copernicus Climate Data Store in the NetCDF format (.nc) and processed into a data frame using R code `codes/process/02-process_environmental_data.R`.

`era5_air_temperature_2009-2022.nc`: Data file containing 2009-2022 daily mean air temperature data in the region containing the Stellako River, which was defined by longitudes 125.0°W – 124.9°W and latitudes 54.0°N – 54.1°N.
The data were downloaded from the Copernicus Climate Data Store in the NetCDF format (.nc) and processed into a data frame using R code `codes/process/02-process_environmental_data.R`.

### Processed data

Data files shown in this section were processed from raw data using codes included in the folder `codes/process` (see their description below).
All are in the `.rds` format and located in the folder `data/processed`.

`count_data.rds`: This data file contains an array of snorkel count records where rows denote survey day, columns denote year of the survey, and the slice denotes the size class.
This data file was produced with code `codes/process/01-process_count_data.R`

`daily_air_temperature_era5.rds`: This data file contains daily mean air temperature from June to September of the study years in the region containing the Stellako River, which was defined by longitudes 125.0°W – 124.9°W and latitudes 54.0°N – 54.1°N.
The data were obtained from ERA5 using R package *ecmwfr* (also available at the [Climate Data Store](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview)).
This data file was produced with code `codes/process/01-process_environmental_data.R`.
The following variables are included: - *date*: Date of the record.
- *air_temp*: Mean air temperature in °C.

`daily_discharge.rds`: This data file contains daily mean discharge of the Stellako River between 1988-2022.
The data were obtained from Water Survey of Canada (WSC) using R package *tidyhydat* (also available at the [WSC website](https://wateroffice.ec.gc.ca/search/historical_e.html)).
This data file was produced with code `codes/process/01-process_environmental_data.R`.
The following variables are included: - *date*: Date of the record.
- *discharge*: Mean discharge in $\text{m}^3 / \text{s}$

`daily_temperature.rds`: This data file contains daily mean water temperature in the Stellako River between 2011-2022.
Hourly water temperature data were obtained upon request from Water Survey of Canada (WSC).
As WSC does not QA/QC water temperature data, we applied QA/QC approaches described by Boyer et al (2017) before computing mean daily water temperature.
These data are not used in the analyses, but are included here because we used it to determine the correlation between air and water temperature.
This data file was produced with code `codes/process/01-process_environmental_data.R`.
The following variables are included: - *date*: Date of the record.
- *discharge*: Mean water temperature in °C

Boyer, C., A. St-Hilaire, and D.
Boivin.
2017.
Quality assured water temperature data from the Water Survey of Canada.
Adapted from report prepared to Water Survey of Canada.

`reach_length.rds`: This data file contains the length of the snorkel survey reaches and was produced with code `codes/process/03-process_reach_length.R`.
The following variables are included: - *site*: Unique ID of the site defining the start of the reach.
- *rm*: River meters at the site.
- *lon*: Longitude at the site.
- *lat*: Latitude of the site.
- *reach*: Number of the reach starting at the site.
- *length*: Length of the reach in meters.

### GIS data

-   `stellako_reaches.kml`: This data file contain coordinates of the points defining the snorkel survey reaches in the Stellako River.
    This file is used by R code `codes/process/03-process_reach_length.R`

-   `stellako_network_riverdist.rds`: This file contains a `rivernetwork` object produced by functions in the package `riverdist`.
    This file is produced and used by R code `codes/process/03-process_reach_length.R`.

### Analysis outputs

Data files listed in this section are produced by codes fitting the N-mixture and von Bertalanffy models to the data and applying the transient life table response experiments (tLTRE) to the N-mixture model output.
All are in the `.rds` format and located in the folder `outputs`.

-   `output_ltre_var_vr.rds`: This file contains the outputs of the tLTRE analysis decomposing variability in population growth rate among vital rates and population structure components.
    It is produced by R function `codes/ltre/00-function_ltre_var_vr.R`.

-   `output_ltre_var_ep.rds`: This file contains the outputs of the tLTRE analysis decomposing variability in population growth rate among the linear predictors of vital rates and population structure components.
    It is produced by R function `codes/ltre/00-function_ltre_var_ep.R`.

-   `output_ltre_diff.rds`: This file contains the outputs of the tLTRE analysis decomposing sequential differences in population growth rate among linear predictors of vital rates and population structure.
    This output is later used to aggregate decompositions at the level of linear predictors to the level of vital rates components.
    It is produced by R function `codes/ltre/00-function_ltre_diff.R`.

-   `output_ltre_lamg.rds`: This file contains the outputs of the tLTRE analysis decomposing differences in the geometric mean of population growth rate between two periods among vital rates and population structure components.
    It is produced by R function `codes/ltre/00-function_ltre_lamg.R`.

-   `spl_dat_beta-binomial_dd_all.rds`: This file is an R list containing the outputs of the integrated, dynamic N-mixture model presented in the manuscript (i.e., using the beta-binomial distribution for the observation model and density dependence in both productivity, $\gamma$ , and survival probability, $\phi$) (in list element `spl`), as well as the data (in list element `dat`) and constants (in list element `cnt`) sent to NIMBLE.
    It is produced by R code `codes/nmix/01-fit_nmix_model.R`.

-   `spl_dat_vbgf.rds`: This file is an R list containing the outputs of the von Bertalanffy growth model (in list element `spl`), as well as the data (in list element `dat`) and constants (in list element `cnt`) sent to NIMBLE.
    It is produced by R code `codes/vbgf/01-fit_vbgf_model.R`.

## Code/Software

The analyses were conducted using the R package nimble 1.3.0 and nimbleEcology 0.4.1 in R 4.4.3.
Data processing as well as summaries and visualization of results were conducted in R 4.4.3 using packages bayesplot 1.12.0, bcdata 0.5.1, cowplot 1.1.3, ecmwfr 1.5.0, ggfan 0.1.4, grid 4.4.3, gridExtra 2.3, HDInterval 0.2.4, MCMCvis 0.16.3, ncdf4 1.24, parallel 4.4.3, RColorBrewer 1.1.3, riverdist 0.17.1, sf 1.0.21, tidyhydat 0.7.1, tidyverse 2.0.0, truncnorm 1.0.9, and zoo 1.8.13

NOTE: Function `dBetaBinom_One`, used in the model code `00-nmix_model.R`, was removed from recent versions of package nimbleEcology ($\ge$ 0.5.0). Please make sure you install nimbleEcology 0.4.1 from the CRAN archive.

### Data processing

The following files contain code to process the data and are located in folder `codes/process`:

-   `01-process_count_data.R`: R code used to process raw snorkel survey count data into an array.
    It produces file: `data/processed/count_data.rds`.

-   `02-process_environmental_data.R`: R code used to process Stellako River discharge and water temperature data from Water Survey of Canada, and air temperature data from ERA5.
    It produces files: `data/processed/daily_discharge.rds`, `data/processed/daily_temperature.rds`, and `data/processed/daily_air_temperature_era5.rds`

-   `03-process_reach_length.R`: R code used to compute the length of snorkel survey reaches using stream network data available at the BC Data Catalogue (via package `bcdata`).
    It produces file: `data/processed/reach_length.rds`

### N-mixture model

The following files contain code to fit the N-mixture model to the data and check the model outputs and are located in folder `codes/nmix`:

-   `00-nmix_model.R`: R code containing the integrated dynamic N-mixture model in the BUGS dialect used by NIMBLE.

-   `00-function_prepare_data.R`: R function used to read in and further process the datasets and constants into a list to be used by NIMBLE.
    This function is called from `codes/nmix/01-fit_nmix_model.R`.

-   `00-function_fit_nmix_workflow.R`: R function used to generate initial values and setup the model to be run in NIMBLE.
    This function is called from `codes/nmix/01-fit_nmix_model.R`.

-   `01-fit_nmix_model.R`: R code used as a wrapper to execute functions `00-function_prepare_data.R` and `00-function_fit_nmix_workflow.R` and setup MCMC runs.
    The output file is saved in `outputs/spl_dat_X_dd_Y.rds`, where X is the distribution used in the observation model (beta-binomial, binomial, or poisson), and Y is the vital rate with density dependence (gamma, phi, or all).

-   `02-check_nmix_fit.qmd`: Quarto file used to generate summaries, statistics, and plots to assess MCMC convergence and model fit.

### tLTRE analysis

The following files contain code to conduct the transient life table response experiments (tLTRE) and are located in folder `codes/ltre`:

-   `00-function_get_samples.R`: R function used to get all MCMC samples from parameters used in the tLTRE analysis.
    This function is called from functions `codes/ltre/00-function_ltre_var_vr.R`, `codes/ltre/00-function_ltre_var_ep.R`, `codes/ltre/00-function_ltre_diff.R`, and `codes/ltre/00-function_ltre_lamg.R`.

-   `00-function_extract_sample.R`: R function used to extract individual MCMC samples and compute linear predictors.
    This function is called from functions `codes/ltre/00-function_ltre_var_vr.R`, `codes/ltre/00-function_ltre_var_ep.R`, `codes/ltre/00-function_ltre_diff.R`, and `codes/ltre/00-function_ltre_lamg.R`.

-   `00-function_compute_sensitivities.R`: R function used to compute sensitivities of population growth rates to changes in vital rates and population structure.
    This function is called from functions `codes/ltre/00-function_ltre_var_vr.R`, `codes/ltre/00-function_ltre_var_ep.R`, and `codes/ltre/00-function_ltre_diff.R`.

-   `00-function_compute_log_diff.R`: R function used to compute log-difference of vital rates and population structure between two periods of equal length.
    This function is called from function `codes/ltre/00-function_ltre_lamg.R`.

-   `00-function_compute_mat_derivatives.R`: R function used to compute derivatives of the transition matrix with respect to vital rates.
    This function is called from function `codes/ltre/00-function_ltre_lamg.R`.

-   `00-function_compute_rte.R`: R function used to compute real time elasticities.
    This function is called from function `codes/ltre/00-function_ltre_lamg.R`.

-   `00-function_ltre_var_vr.R`: R function used to conduct tLTRE decomposing variability in population growth rate among vital rates and population structure components.
    This function is called from R code `codes/ltre/01-ltre_analysis.R`.

-   `00-function_ltre_var_ep.R`: R function used to conduct tLTRE decomposing variability in population growth rate among linear predictors of vital rates and population structure components.
    This function is called from R code `codes/ltre/01-ltre_analysis.R`.

-   `00-function_ltre_diff.R`: R function used to conduct tLTRE decomposing sequential differences in population growth rate among linear predictors of vital rates and population structure.
    The output is later used to aggregate decompositions at the level of linear predictors to the level of vital rates components.
    This function is called from R code `codes/ltre/01-ltre_analysis.R`.

-   `00-function_ltre_lamg.R`: R function used to conduct tLTRE decomposing differences in the geometric mean of population growth rate between two periods among vital rates and population structure components.
    This function is called from R code `codes/ltre/01-ltre_analysis.R`.

-   `01-ltre_analysis.R`: R code used as a wrapper to execute functions `00-function_ltre_var_vr.R`, `00-function_ltre_var_ep.R`, `00-function_ltre_diff.R`, and `00-function_ltre_lamg.R`.
    The output files are saved in `outputs/output_ltre_var_vr.rds`, `outputs/output_ltre_var_vep.rds`, `outputs/output_ltre_diff.rds`, and `outputs/output_ltre_lamg.rds`

### von Bertalanffy growth model

The following files contain code to fit the von Bertalanffy growth model to the length-at-age data and check the model outputs and are located in folder `codes/vbgf`:

-   `00-vbgf_model.R`: R code containing the von Bertalanffy growth model in the BUGS dialect used by NIMBLE.

-   `01-fit_vbgf_model.R`: R code used to setup MCMC runs and fit the von Bertalanffy growth model to the data.
    The output file is saved in `outputs/spl_dat_vbgf.rds`.

-   `02-check_vbgf_fit.qmd`: Quarto file used to generate summaries, statistics, and plots to assess MCMC convergence and fit of the von Bertalanffy growth model.

### Figures

The following files contain code to generate figures presented in the manuscript and in the appendices and are located in folder `codes/figures`

-   `main-figures.qmd`: Quarto file used to generate figures 3-9 presented in the manuscript.

-   `appendix-s1-figure.qmd`: Quarto file used to generate figures presented in Appendix S1.

-   `appendix-s2-figure.qmd`: Quarto file used to generate figures presented in Appendix S2.

## Important instructions

If the user is not cloning the GitHub repository with these files, make sure to organize the files in a project folder with the structure below, otherwise the user will have to manually edit file paths in the .R and .qmd files.

```         
├── codes
│   ├── figures
│   │   ├── appendix-s1-figures.qmd
│   │   ├── appendix-s2-figures.qmd
│   │   └── main-figures.qmd
│   ├── ltre
│   │   ├── 00-function_compute_log_diff.R
│   │   ├── 00-function_compute_mat_derivatives.R
│   │   ├── 00-function_compute_rte.R
│   │   ├── 00-function_compute_sensitivities.R
│   │   ├── 00-function_extract_sample.R
│   │   ├── 00-function_get_samples.R
│   │   ├── 00-function_ltre_diff.R
│   │   ├── 00-function_ltre_lamg.R
│   │   ├── 00-function_ltre_var_ep.R
│   │   ├── 00-function_ltre_var_vr.R
│   │   └── 01-ltre_analysis.R
│   ├── nmix
│   │   ├── 00-function_fit_nmix_workflow.R
│   │   ├── 00-function_prepare_data.R
│   │   ├── 00-nmix_model.R
│   │   ├── 01-fit_nmix_model.R
│   │   └── 02-check_nmix_fit.qmd
│   ├── process
│   │   ├── 01-process_count_data.R
│   │   ├── 02-process_environmental_data.R
│   │   └── 03-process_reach_length.R
│   └── vbgf
│       ├── 00-vbgf_model.R
│       ├── 01-fit_vbgf_model.R
│       └── 02-check_vbgf_fit.qmd
├── data
│   ├── gis
│   │   ├── stellako_network_riverdist.rds
│   │   └── stellako_reaches.kml
│   ├── processed
│   │   ├── count_data.rds
│   │   ├── daily_air_temperature_era5.rds
│   │   ├── daily_discharge.rds
│   │   ├── daily_temperature.rds
│   │   ├── daily_water_temperature.rds
│   │   └── reach_lengths.rds
│   └── raw
│       ├── age_length.csv
│       ├── conditions.csv
│       ├── era5_air_temperature_1980-1994.nc
│       ├── era5_air_temperature_1995-2008.nc
│       ├── era5_air_temperature_2009-2022.nc
│       ├── lane_reach.csv
│       ├── rb_counts.csv
│       ├── sk_nuseds.csv
│       └── wsc_water_temperature_2011-2024.csv
├── manuscript-repo.Rproj
└── outputs
    ├── output_ltre_diff.rds
    ├── output_ltre_lamg.rds
    ├── output_ltre_var_ep.rds
    ├── output_ltre_var_vr.rds
    ├── spl_dat_nmix_beta-binomial_dd_all.rds
    └── spl_dat_vbgf.rds
```

## License

This project is licensed under the following terms:

-   **Code**: The code in this repository is licensed under the [GNU General Public License v3.0 (GPL-3.0)](https://www.gnu.org/licenses/gpl-3.0.html). This means you are free to use, modify, and distribute the code, but any derivative work must also be released under the same license.
-   **Data**: The data in this repository is licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/). You are free to share and adapt the data, provided that you give appropriate credit.

By using this repository, you agree to comply with the terms of these licenses.
