# Readme
Code and data originally to reproduce analyses presented in Seaman et al, A global investment case for hepatitis B elimination – a modelling study, 2023 (Under Review). The code has now been expanded to perform additional regional/country analyses, including work for the Vaccine Impact Modelling Consortium (VIMC).

The hepatitis B model is run using Python (≥v3.7). Make sure you have the packages numpy, pandas, matplotlib and Atomica ( https://github.com/atomicateam/atomica) installed before you run the model.

Full details on implementation of compartmental models in Atomica can be found at: atomica.tools

## applications
THe applications folder contains a list of all country and regional applications currently available, normally characteristed by a prefix of the project it's under (currently either region or vimc), and the specific region or country of the application. 

Within each application folder, there is an associated databook and calibration file for each application. Regional analyses also contain additional databooks which have datapoints to represent different scenario targets (2030 - 2050 for S1 - S3 respectively). 

The VIMC analysis folder also contains specific functions and scripts designed for running and plotting VIMC analyses (vimc_functions and vimc_plotting).

## cost and agg data
Data files relevant to running the `model_results` function.

## documents
Miscellaneous files collected from functions for meetings. 

## frameworks
A folder containing all the different versions of the framework at the moment. The current list of frameworks include:
- hbv_v14_gamma_max.xlsx (regional analysis)
- hbv_v14_gamma_vimc.xlsx (vimc analysis)

## hbv_functions
A python module containing all relevant functions for all hbv analyses. Functions include: collecting model results, editing scenarios, generating calibrations, and creating calibration error tables. Specific details of each function is listed in the `README.md` in `hbv_functions`. 

## output_plots
A list of plots created by the `epi_plots` function.

## scripts
A list of example scripts for hbv analyses. These include how to use functions in `hbv_functions` to run model results, calibrating models, plotting results, and generating error tables. 

