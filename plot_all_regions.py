import os
wd= 'C:/Users/iamph/Documents/GitHub/hbv-globalinvcase'#set own working directory
os.chdir(wd)

import numpy as np
import atomica as at
import pandas as pd

import hbv_functions as hbv

F=at.ProjectFramework("hbv_v14_gamma_mav.xlsx") #updated to include maternal antivirals
runs=2 #number of model simulations
regions = ['AFR', 'AMR', 'EMR', 'EUR', 'SEAR', 'WPR']

for ct in regions:
    bl_runs, bl_cent = hbv.model_results(F, 'regional', f"{ct}_db_mav.xlsx", f"{ct}_calib.xlsx", "Status Quo", runs)
    s1_runs, s1_cent = hbv.model_results(F, 'regional', f"{ct}_db_s1_mav.xlsx", f"{ct}_calib.xlsx", "S1: 2040 Target", runs)
    s2_runs, s2_cent = hbv.model_results(F, 'regional', f"{ct}_db_s2_mav.xlsx", f"{ct}_calib.xlsx", "S2: 2040 Target", runs)
    s3_runs, s3_cent = hbv.model_results(F, 'regional', f"{ct}_db_s3_mav.xlsx", f"{ct}_calib.xlsx", "S3: 2050 Target", runs)

    hbv.epi_plots(bl_cent, s1_cent, s2_cent,s3_cent, bl_runs, s1_runs, s2_runs, s3_runs, wd, "/cost and agg data/pop_calib.xlsx", ct)
