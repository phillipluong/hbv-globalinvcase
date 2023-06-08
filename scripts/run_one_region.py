import os
# wd= 'C:/Users/iamph/Documents/GitHub/hbv-globalinvcase'#set own working directory
# os.chdir(wd)

import atomica as at

import hbv_functions as hbv
hbv.get_gitlab_folder()

## Run this separately
F=at.ProjectFramework("frameworks/hbv_v14_gamma_mav.xlsx") #updated to include maternal antivirals
runs=2 #number of model simulations
ct = 'EMR'

# bl_runs, bl_cent = hbv.model_results(F, 'regional', f"{ct}_db_mav.xlsx", f"{ct}_calib.xlsx", "Status Quo", runs)
s1_runs, s1_cent = hbv.model_results(F, f'region_{ct.lower()}', f"{ct}_db_s1_mav.xlsx", f"{ct}_calib.xlsx", "S1: 2040 Target", runs)
# s2_runs, s2_cent = hbv.model_results(F, 'regional', f"{ct}_db_s2_mav.xlsx", f"{ct}_calib.xlsx", "S2: 2040 Target", runs)
# s3_runs, s3_cent = hbv.model_results(F, 'regional', f"{ct}_db_s3_mav.xlsx", f"{ct}_calib.xlsx", "S3: 2050 Target", runs)

# econ= hbv.econ_analysis (bl_cent, s1_cent, s2_cent, s3_cent,bl_runs, s1_runs, s2_runs, s3_runs, "AFR", wd, "cost and agg data/costs.xlsx", runs, 0.03, 0.03)

# hbv.epi_plots(bl_cent, s1_cent, s2_cent,s3_cent, bl_runs, s1_runs, s2_runs, s3_runs, wd, "cost and agg data/pop_calib.xlsx", ct

# hbv.econ_plots (econ, ct)

### ------ Still to add onto the code

'''
#Generate weights for prevalence aggregation (central)
bl_cent["prev_w"]=bl_cent["prev"]*bl_cent["pop"]
s1_cent["prev_w"]=s1_cent["prev"]*s1_cent["pop"]
s2_cent["prev_w"]=s2_cent["prev"]*s2_cent["pop"]
s3_cent["prev_w"]=s3_cent["prev"]*s3_cent["pop"]

bl_cent["prev_u5_w"]=bl_cent["prev_u5"]*bl_cent["pop_u5"]
s1_cent["prev_u5_w"]=s1_cent["prev_u5"]*s1_cent["pop_u5"]
s2_cent["prev_u5_w"]=s2_cent["prev_u5"]*s2_cent["pop_u5"]
s3_cent["prev_u5_w"]=s3_cent["prev_u5"]*s3_cent["pop_u5"]

#Generate weights for prevalence aggregation (sampled)
bl_runs["prev_w"]=bl_runs["prev"]*bl_runs["pop"]
s1_runs["prev_w"]=s1_runs["prev"]*s1_runs["pop"]
s2_runs["prev_w"]=s2_runs["prev"]*s2_runs["pop"]
s3_runs["prev_w"]=s3_runs["prev"]*s3_runs["pop"]

bl_runs["prev_u5_w"]=bl_runs["prev_u5"]*bl_runs["pop_u5"]
s1_runs["prev_u5_w"]=s1_runs["prev_u5"]*s1_runs["pop_u5"]
s2_runs["prev_u5_w"]=s2_runs["prev_u5"]*s2_runs["pop_u5"]
s3_runs["prev_u5_w"]=s3_runs["prev_u5"]*s3_runs["pop_u5"]

## MTCT worst case
mtctw_bl, mtctw_s1= par_sens(F, "db_mav.xlsx", "db_s1_mav.xlsx", "calib.xlsx", mtct_wc)

## MTCT best case
mtctb_bl, mtctb_s1=  par_sens(F, "db_mav.xlsx", "db_s1_mav.xlsx", "calib.xlsx", mtct_bc)

## Mortality worst case
mortw_bl, mortw_s1= par_sens(F, "db_mav.xlsx", "db_s1_mav.xlsx", "calib.xlsx", mort_wc)

## Mortality best case
mortb_bl, mortb_s1=  par_sens(F, "db_mav.xlsx", "db_s1_mav.xlsx", "calib.xlsx", mort_bc)

## Treatment worst case (dis prog)
trtw_bl, trtw_s1= par_sens(F, "db_mav.xlsx", "db_s1_mav.xlsx", "calib.xlsx", tex_wc)

## Treatment best case (dis prog)
trtb_bl, trtb_s1=  par_sens(F, "db_mav.xlsx", "db_s1_mav.xlsx", "calib.xlsx", tex_bc)

## MTCT worst case
mtctw_s1=parsens_econ(mtctw_bl, mtctw_s1, "AFR", wd, "cost and agg data/costs.xlsx", 0.03)

## MTCT best case
mtctb_s1=parsens_econ(mtctb_bl, mtctb_s1, "AFR", wd, "cost and agg data/costs.xlsx", 0.03)

## Mortality worst case
mortw_s1=parsens_econ(mortw_bl, mortw_s1, "AFR", wd, "cost and agg data/costs.xlsx", 0.03)

## Mortality best case
mortb_s1=parsens_econ(mortb_bl, mortb_s1, "AFR", wd, "cost and agg data/costs.xlsx", 0.03)

## Treatment worst case
trtw_s1=parsens_econ(trtw_bl, trtw_s1, "AFR", wd, "cost and agg data/costs.xlsx", 0.03)

## Treatment best case
trtb_s1=parsens_econ(trtb_bl, trtb_s1, "AFR", wd, "cost and agg data/costs.xlsx", 0.03)

esens_outs=sens_nebicer(econ_sens, econ)

sens_tab=sens_tables(s1_cent, econ, mtctw_s1, mtctb_s1, mortw_s1, mortb_s1, trtw_s1, trtb_s1, econ_sens, esens_outs)
'''