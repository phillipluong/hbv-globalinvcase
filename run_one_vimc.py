import os
wd= 'C:/Users/Phil/Documents/GitHub/hbv-globalinvcase'#set own working directory
os.chdir(wd)

import numpy as np
import atomica as at
import pandas as pd

import hbv_functions as hbv

## Run this separately
F=at.ProjectFramework("hbv_v14_gamma_vimc.xlsx") #updated to include maternal antivirals
runs=20 #number of model simulations

afr_bl_runs, afr_bl_cent = hbv.model_results(F, 'vimc', "vimc_AFR_db_v1_3.xlsx", "vimc_AFR_calib_v1_2.xlsx", "vimc", runs)

## Model results: it runs, but there are couple of errors from the outputs
## There will be some value taking a look into these errors and disussing these with Chris!
'''
## Is there a need for econ analysis for vimc? not at the moment
afr_econ= hbv.econ_analysis (afr_bl_cent, afr_s1_cent, afr_s2_cent, afr_s3_cent,afr_bl_runs, afr_s1_runs, afr_s2_runs, afr_s3_runs, "AFR", wd, "cost and agg data/costs.xlsx", runs, 0.03, 0.03)

hbv.epi_plots(afr_bl_cent, afr_s1_cent, afr_s2_cent,afr_s3_cent, afr_bl_runs, afr_s1_runs, afr_s2_runs, afr_s3_runs, wd, "cost and agg data/pop_calib.xlsx", "AFR")

hbv.econ_plots (afr_econ, "AFR")
'''

### ------ Still to add onto the code

'''
#Generate weights for prevalence aggregation (central)
afr_bl_cent["prev_w"]=afr_bl_cent["prev"]*afr_bl_cent["pop"]
afr_s1_cent["prev_w"]=afr_s1_cent["prev"]*afr_s1_cent["pop"]
afr_s2_cent["prev_w"]=afr_s2_cent["prev"]*afr_s2_cent["pop"]
afr_s3_cent["prev_w"]=afr_s3_cent["prev"]*afr_s3_cent["pop"]

afr_bl_cent["prev_u5_w"]=afr_bl_cent["prev_u5"]*afr_bl_cent["pop_u5"]
afr_s1_cent["prev_u5_w"]=afr_s1_cent["prev_u5"]*afr_s1_cent["pop_u5"]
afr_s2_cent["prev_u5_w"]=afr_s2_cent["prev_u5"]*afr_s2_cent["pop_u5"]
afr_s3_cent["prev_u5_w"]=afr_s3_cent["prev_u5"]*afr_s3_cent["pop_u5"]

#Generate weights for prevalence aggregation (sampled)
afr_bl_runs["prev_w"]=afr_bl_runs["prev"]*afr_bl_runs["pop"]
afr_s1_runs["prev_w"]=afr_s1_runs["prev"]*afr_s1_runs["pop"]
afr_s2_runs["prev_w"]=afr_s2_runs["prev"]*afr_s2_runs["pop"]
afr_s3_runs["prev_w"]=afr_s3_runs["prev"]*afr_s3_runs["pop"]

afr_bl_runs["prev_u5_w"]=afr_bl_runs["prev_u5"]*afr_bl_runs["pop_u5"]
afr_s1_runs["prev_u5_w"]=afr_s1_runs["prev_u5"]*afr_s1_runs["pop_u5"]
afr_s2_runs["prev_u5_w"]=afr_s2_runs["prev_u5"]*afr_s2_runs["pop_u5"]
afr_s3_runs["prev_u5_w"]=afr_s3_runs["prev_u5"]*afr_s3_runs["pop_u5"]

## MTCT worst case
afr_mtctw_bl, afr_mtctw_s1= par_sens(F, "AFR_db_mav.xlsx", "AFR_db_s1_mav.xlsx", "AFR_calib.xlsx", mtct_wc)

## MTCT best case
afr_mtctb_bl, afr_mtctb_s1=  par_sens(F, "AFR_db_mav.xlsx", "AFR_db_s1_mav.xlsx", "AFR_calib.xlsx", mtct_bc)

## Mortality worst case
afr_mortw_bl, afr_mortw_s1= par_sens(F, "AFR_db_mav.xlsx", "AFR_db_s1_mav.xlsx", "AFR_calib.xlsx", mort_wc)

## Mortality best case
afr_mortb_bl, afr_mortb_s1=  par_sens(F, "AFR_db_mav.xlsx", "AFR_db_s1_mav.xlsx", "AFR_calib.xlsx", mort_bc)

## Treatment worst case (dis prog)
afr_trtw_bl, afr_trtw_s1= par_sens(F, "AFR_db_mav.xlsx", "AFR_db_s1_mav.xlsx", "AFR_calib.xlsx", tex_wc)

## Treatment best case (dis prog)
afr_trtb_bl, afr_trtb_s1=  par_sens(F, "AFR_db_mav.xlsx", "AFR_db_s1_mav.xlsx", "AFR_calib.xlsx", tex_bc)

## MTCT worst case
afr_mtctw_s1=parsens_econ(afr_mtctw_bl, afr_mtctw_s1, "AFR", wd, "cost and agg data/costs.xlsx", 0.03)

## MTCT best case
afr_mtctb_s1=parsens_econ(afr_mtctb_bl, afr_mtctb_s1, "AFR", wd, "cost and agg data/costs.xlsx", 0.03)

## Mortality worst case
afr_mortw_s1=parsens_econ(afr_mortw_bl, afr_mortw_s1, "AFR", wd, "cost and agg data/costs.xlsx", 0.03)

## Mortality best case
afr_mortb_s1=parsens_econ(afr_mortb_bl, afr_mortb_s1, "AFR", wd, "cost and agg data/costs.xlsx", 0.03)

## Treatment worst case
afr_trtw_s1=parsens_econ(afr_trtw_bl, afr_trtw_s1, "AFR", wd, "cost and agg data/costs.xlsx", 0.03)

## Treatment best case
afr_trtb_s1=parsens_econ(afr_trtb_bl, afr_trtb_s1, "AFR", wd, "cost and agg data/costs.xlsx", 0.03)

afr_esens_outs=sens_nebicer(afr_econ_sens, afr_econ)

afr_sens_tab=sens_tables(afr_s1_cent, afr_econ, afr_mtctw_s1, afr_mtctb_s1, afr_mortw_s1, afr_mortb_s1, afr_trtw_s1, afr_trtb_s1, afr_econ_sens, afr_esens_outs)
'''