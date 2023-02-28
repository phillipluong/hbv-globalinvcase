# Readme
Code and data to reproduce analyses presented in Seaman et al, A global investment case for hepatitis B elimination – a modelling study, 2023 (Under Review).

The hepatitis B model is run using Python (≥v3.7). Make sure you have the packages numpy, pandas, matplotlib and Atomica ( https://github.com/atomicateam/atomica) installed before you run the model.

Full details on implementation of compartmental models in Atomica can be found at: atomica.tools

## hbv_v14_gamma_mav.xlsx
Framework file which contains descriptions of model parameters, model compartment details, and the transition matrix which determines progression within the model. 

## Databooks
The databooks folder contains a series of 24 databooks, each an individual scenario (baseline, 2030 target [S1], 2040 target [S2], 2050 target [S3]) in a given WHO region (AFR, AMR, EMR, EUR, SEAR, WPR). Databooks contain parameter values.

Calibration and cost and agg data folders include spreadsheets used for model calibration (completed within the hbvincase_fullmodelrun.py script), plotting outputs, and economic inputs. Any economic inputs (except YLD weights) can be ignored, as these represent unused parameters within the model (a reminant of model development).
