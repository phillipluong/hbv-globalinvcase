## Generate Required Data for Global Epi Outcomes (central)
global_bl_cent={}
global_bl_cent["prev"]=(afr_bl_cent["prev_w"]+amr_bl_cent["prev_w"]+emr_bl_cent["prev_w"]+eur_bl_cent["prev_w"]+sear_bl_cent["prev_w"]+wpr_bl_cent["prev_w"])/(afr_bl_cent["pop"]+amr_bl_cent["pop"]+emr_bl_cent["pop"]+eur_bl_cent["pop"]+sear_bl_cent["pop"]+wpr_bl_cent["pop"])
global_bl_cent["prev_u5"]=(afr_bl_cent["prev_u5_w"]+amr_bl_cent["prev_u5_w"]+emr_bl_cent["prev_u5_w"]+eur_bl_cent["prev_u5_w"]+sear_bl_cent["prev_u5_w"]+wpr_bl_cent["prev_u5_w"])/(afr_bl_cent["pop_u5"]+amr_bl_cent["pop_u5"]+emr_bl_cent["pop_u5"]+eur_bl_cent["pop_u5"]+sear_bl_cent["pop_u5"]+wpr_bl_cent["pop_u5"])
global_bl_cent["hcc_inc"]=afr_bl_cent["hcc_inc"]+amr_bl_cent["hcc_inc"]+emr_bl_cent["hcc_inc"]+eur_bl_cent["hcc_inc"]+sear_bl_cent["hcc_inc"]+wpr_bl_cent["hcc_inc"]
global_bl_cent["mort"]=afr_bl_cent["mort"]+amr_bl_cent["mort"]+emr_bl_cent["mort"]+eur_bl_cent["mort"]+sear_bl_cent["mort"]+wpr_bl_cent["mort"]
global_bl_cent["chb_inc"]=afr_bl_cent["chb_inc"]+amr_bl_cent["chb_inc"]+emr_bl_cent["chb_inc"]+eur_bl_cent["chb_inc"]+sear_bl_cent["chb_inc"]+wpr_bl_cent["chb_inc"]

global_s1_cent={}
global_s1_cent["prev"]=(afr_s1_cent["prev_w"]+amr_s1_cent["prev_w"]+emr_s1_cent["prev_w"]+eur_s1_cent["prev_w"]+sear_s1_cent["prev_w"]+wpr_s1_cent["prev_w"])/(afr_s1_cent["pop"]+amr_s1_cent["pop"]+emr_s1_cent["pop"]+eur_s1_cent["pop"]+sear_s1_cent["pop"]+wpr_s1_cent["pop"])
global_s1_cent["prev_u5"]=(afr_s1_cent["prev_u5_w"]+amr_s1_cent["prev_u5_w"]+emr_s1_cent["prev_u5_w"]+eur_s1_cent["prev_u5_w"]+sear_s1_cent["prev_u5_w"]+wpr_s1_cent["prev_u5_w"])/(afr_s1_cent["pop_u5"]+amr_s1_cent["pop_u5"]+emr_s1_cent["pop_u5"]+eur_s1_cent["pop_u5"]+sear_s1_cent["pop_u5"]+wpr_s1_cent["pop_u5"])
global_s1_cent["hcc_inc"]=afr_s1_cent["hcc_inc"]+amr_s1_cent["hcc_inc"]+emr_s1_cent["hcc_inc"]+eur_s1_cent["hcc_inc"]+sear_s1_cent["hcc_inc"]+wpr_s1_cent["hcc_inc"]
global_s1_cent["mort"]=afr_s1_cent["mort"]+amr_s1_cent["mort"]+emr_s1_cent["mort"]+eur_s1_cent["mort"]+sear_s1_cent["mort"]+wpr_s1_cent["mort"]
global_s1_cent["chb_inc"]=afr_s1_cent["chb_inc"]+amr_s1_cent["chb_inc"]+emr_s1_cent["chb_inc"]+eur_s1_cent["chb_inc"]+sear_s1_cent["chb_inc"]+wpr_s1_cent["chb_inc"]

global_s2_cent={}
global_s2_cent["prev"]=(afr_s2_cent["prev_w"]+amr_s2_cent["prev_w"]+emr_s2_cent["prev_w"]+eur_s2_cent["prev_w"]+sear_s2_cent["prev_w"]+wpr_s2_cent["prev_w"])/(afr_s2_cent["pop"]+amr_s2_cent["pop"]+emr_s2_cent["pop"]+eur_s2_cent["pop"]+sear_s2_cent["pop"]+wpr_s2_cent["pop"])
global_s2_cent["prev_u5"]=(afr_s2_cent["prev_u5_w"]+amr_s2_cent["prev_u5_w"]+emr_s2_cent["prev_u5_w"]+eur_s2_cent["prev_u5_w"]+sear_s2_cent["prev_u5_w"]+wpr_s2_cent["prev_u5_w"])/(afr_s2_cent["pop_u5"]+amr_s2_cent["pop_u5"]+emr_s2_cent["pop_u5"]+eur_s2_cent["pop_u5"]+sear_s2_cent["pop_u5"]+wpr_s2_cent["pop_u5"])
global_s2_cent["hcc_inc"]=afr_s2_cent["hcc_inc"]+amr_s2_cent["hcc_inc"]+emr_s2_cent["hcc_inc"]+eur_s2_cent["hcc_inc"]+sear_s2_cent["hcc_inc"]+wpr_s2_cent["hcc_inc"]
global_s2_cent["mort"]=afr_s2_cent["mort"]+amr_s2_cent["mort"]+emr_s2_cent["mort"]+eur_s2_cent["mort"]+sear_s2_cent["mort"]+wpr_s2_cent["mort"]
global_s2_cent["chb_inc"]=afr_s2_cent["chb_inc"]+amr_s2_cent["chb_inc"]+emr_s2_cent["chb_inc"]+eur_s2_cent["chb_inc"]+sear_s2_cent["chb_inc"]+wpr_s2_cent["chb_inc"]

global_s3_cent={}
global_s3_cent["prev"]=(afr_s3_cent["prev_w"]+amr_s3_cent["prev_w"]+emr_s3_cent["prev_w"]+eur_s3_cent["prev_w"]+sear_s3_cent["prev_w"]+wpr_s3_cent["prev_w"])/(afr_s3_cent["pop"]+amr_s3_cent["pop"]+emr_s3_cent["pop"]+eur_s3_cent["pop"]+sear_s3_cent["pop"]+wpr_s3_cent["pop"])
global_s3_cent["prev_u5"]=(afr_s3_cent["prev_u5_w"]+amr_s3_cent["prev_u5_w"]+emr_s3_cent["prev_u5_w"]+eur_s3_cent["prev_u5_w"]+sear_s3_cent["prev_u5_w"]+wpr_s3_cent["prev_u5_w"])/(afr_s3_cent["pop_u5"]+amr_s3_cent["pop_u5"]+emr_s3_cent["pop_u5"]+eur_s3_cent["pop_u5"]+sear_s3_cent["pop_u5"]+wpr_s3_cent["pop_u5"])
global_s3_cent["hcc_inc"]=afr_s3_cent["hcc_inc"]+amr_s3_cent["hcc_inc"]+emr_s3_cent["hcc_inc"]+eur_s3_cent["hcc_inc"]+sear_s3_cent["hcc_inc"]+wpr_s3_cent["hcc_inc"]
global_s3_cent["mort"]=afr_s3_cent["mort"]+amr_s3_cent["mort"]+emr_s3_cent["mort"]+eur_s3_cent["mort"]+sear_s3_cent["mort"]+wpr_s3_cent["mort"]
global_s3_cent["chb_inc"]=afr_s3_cent["chb_inc"]+amr_s3_cent["chb_inc"]+emr_s3_cent["chb_inc"]+eur_s3_cent["chb_inc"]+sear_s3_cent["chb_inc"]+wpr_s3_cent["chb_inc"]

## Generate Required Data for Global Epi Outcomes (sampled)
global_bl_runs={}
global_bl_runs["prev"]=(afr_bl_runs["prev_w"]+amr_bl_runs["prev_w"]+emr_bl_runs["prev_w"]+eur_bl_runs["prev_w"]+sear_bl_runs["prev_w"]+wpr_bl_runs["prev_w"])/(afr_bl_runs["pop"]+amr_bl_runs["pop"]+emr_bl_runs["pop"]+eur_bl_runs["pop"]+sear_bl_runs["pop"]+wpr_bl_runs["pop"])
global_bl_runs["prev_u5"]=(afr_bl_runs["prev_u5_w"]+amr_bl_runs["prev_u5_w"]+emr_bl_runs["prev_u5_w"]+eur_bl_runs["prev_u5_w"]+sear_bl_runs["prev_u5_w"]+wpr_bl_runs["prev_u5_w"])/(afr_bl_runs["pop_u5"]+amr_bl_runs["pop_u5"]+emr_bl_runs["pop_u5"]+eur_bl_runs["pop_u5"]+sear_bl_runs["pop_u5"]+wpr_bl_runs["pop_u5"])
global_bl_runs["hcc_inc"]=afr_bl_runs["hcc_inc"]+amr_bl_runs["hcc_inc"]+emr_bl_runs["hcc_inc"]+eur_bl_runs["hcc_inc"]+sear_bl_runs["hcc_inc"]+wpr_bl_runs["hcc_inc"]
global_bl_runs["mort"]=afr_bl_runs["mort"]+amr_bl_runs["mort"]+emr_bl_runs["mort"]+eur_bl_runs["mort"]+sear_bl_runs["mort"]+wpr_bl_runs["mort"]
global_bl_runs["chb_inc"]=afr_bl_runs["chb_inc"]+amr_bl_runs["chb_inc"]+emr_bl_runs["chb_inc"]+eur_bl_runs["chb_inc"]+sear_bl_runs["chb_inc"]+wpr_bl_runs["chb_inc"]

global_s1_runs={}
global_s1_runs["prev"]=(afr_s1_runs["prev_w"]+amr_s1_runs["prev_w"]+emr_s1_runs["prev_w"]+eur_s1_runs["prev_w"]+sear_s1_runs["prev_w"]+wpr_s1_runs["prev_w"])/(afr_s1_runs["pop"]+amr_s1_runs["pop"]+emr_s1_runs["pop"]+eur_s1_runs["pop"]+sear_s1_runs["pop"]+wpr_s1_runs["pop"])
global_s1_runs["prev_u5"]=(afr_s1_runs["prev_u5_w"]+amr_s1_runs["prev_u5_w"]+emr_s1_runs["prev_u5_w"]+eur_s1_runs["prev_u5_w"]+sear_s1_runs["prev_u5_w"]+wpr_s1_runs["prev_u5_w"])/(afr_s1_runs["pop_u5"]+amr_s1_runs["pop_u5"]+emr_s1_runs["pop_u5"]+eur_s1_runs["pop_u5"]+sear_s1_runs["pop_u5"]+wpr_s1_runs["pop_u5"])
global_s1_runs["hcc_inc"]=afr_s1_runs["hcc_inc"]+amr_s1_runs["hcc_inc"]+emr_s1_runs["hcc_inc"]+eur_s1_runs["hcc_inc"]+sear_s1_runs["hcc_inc"]+wpr_s1_runs["hcc_inc"]
global_s1_runs["mort"]=afr_s1_runs["mort"]+amr_s1_runs["mort"]+emr_s1_runs["mort"]+eur_s1_runs["mort"]+sear_s1_runs["mort"]+wpr_s1_runs["mort"]
global_s1_runs["chb_inc"]=afr_s1_runs["chb_inc"]+amr_s1_runs["chb_inc"]+emr_s1_runs["chb_inc"]+eur_s1_runs["chb_inc"]+sear_s1_runs["chb_inc"]+wpr_s1_runs["chb_inc"]

global_s2_runs={}
global_s2_runs["prev"]=(afr_s2_runs["prev_w"]+amr_s2_runs["prev_w"]+emr_s2_runs["prev_w"]+eur_s2_runs["prev_w"]+sear_s2_runs["prev_w"]+wpr_s2_runs["prev_w"])/(afr_s2_runs["pop"]+amr_s2_runs["pop"]+emr_s2_runs["pop"]+eur_s2_runs["pop"]+sear_s2_runs["pop"]+wpr_s2_runs["pop"])
global_s2_runs["prev_u5"]=(afr_s2_runs["prev_u5_w"]+amr_s2_runs["prev_u5_w"]+emr_s2_runs["prev_u5_w"]+eur_s2_runs["prev_u5_w"]+sear_s2_runs["prev_u5_w"]+wpr_s2_runs["prev_u5_w"])/(afr_s2_runs["pop_u5"]+amr_s2_runs["pop_u5"]+emr_s2_runs["pop_u5"]+eur_s2_runs["pop_u5"]+sear_s2_runs["pop_u5"]+wpr_s2_runs["pop_u5"])
global_s2_runs["hcc_inc"]=afr_s2_runs["hcc_inc"]+amr_s2_runs["hcc_inc"]+emr_s2_runs["hcc_inc"]+eur_s2_runs["hcc_inc"]+sear_s2_runs["hcc_inc"]+wpr_s2_runs["hcc_inc"]
global_s2_runs["mort"]=afr_s2_runs["mort"]+amr_s2_runs["mort"]+emr_s2_runs["mort"]+eur_s2_runs["mort"]+sear_s2_runs["mort"]+wpr_s2_runs["mort"]
global_s2_runs["chb_inc"]=afr_s2_runs["chb_inc"]+amr_s2_runs["chb_inc"]+emr_s2_runs["chb_inc"]+eur_s2_runs["chb_inc"]+sear_s2_runs["chb_inc"]+wpr_s2_runs["chb_inc"]

global_s3_runs={}
global_s3_runs["prev"]=(afr_s3_runs["prev_w"]+amr_s3_runs["prev_w"]+emr_s3_runs["prev_w"]+eur_s3_runs["prev_w"]+sear_s3_runs["prev_w"]+wpr_s3_runs["prev_w"])/(afr_s3_runs["pop"]+amr_s3_runs["pop"]+emr_s3_runs["pop"]+eur_s3_runs["pop"]+sear_s3_runs["pop"]+wpr_s3_runs["pop"])
global_s3_runs["prev_u5"]=(afr_s3_runs["prev_u5_w"]+amr_s3_runs["prev_u5_w"]+emr_s3_runs["prev_u5_w"]+eur_s3_runs["prev_u5_w"]+sear_s3_runs["prev_u5_w"]+wpr_s3_runs["prev_u5_w"])/(afr_s3_runs["pop_u5"]+amr_s3_runs["pop_u5"]+emr_s3_runs["pop_u5"]+eur_s3_runs["pop_u5"]+sear_s3_runs["pop_u5"]+wpr_s3_runs["pop_u5"])
global_s3_runs["hcc_inc"]=afr_s3_runs["hcc_inc"]+amr_s3_runs["hcc_inc"]+emr_s3_runs["hcc_inc"]+eur_s3_runs["hcc_inc"]+sear_s3_runs["hcc_inc"]+wpr_s3_runs["hcc_inc"]
global_s3_runs["mort"]=afr_s3_runs["mort"]+amr_s3_runs["mort"]+emr_s3_runs["mort"]+eur_s3_runs["mort"]+sear_s3_runs["mort"]+wpr_s3_runs["mort"]
global_s3_runs["chb_inc"]=afr_s3_runs["chb_inc"]+amr_s3_runs["chb_inc"]+emr_s3_runs["chb_inc"]+eur_s3_runs["chb_inc"]+sear_s3_runs["chb_inc"]+wpr_s3_runs["chb_inc"]


epi_plots(global_bl_cent, global_s1_cent, global_s2_cent, global_s3_cent,global_bl_runs, global_s1_runs, global_s2_runs, global_s3_runs, wd, "cost and agg data/pop_calib.xlsx", "global")

## Generate Required Data for Global Economic Outcomes
global_econ={}

## Add Intervention Costs
global_econ["cbl_intc"]=afr_econ["cbl_intc"]+amr_econ["cbl_intc"]+emr_econ["cbl_intc"]+eur_econ["cbl_intc"]+sear_econ["cbl_intc"]+wpr_econ["cbl_intc"]
global_econ["cs1_intc"]=afr_econ["cs1_intc"]+amr_econ["cs1_intc"]+emr_econ["cs1_intc"]+eur_econ["cs1_intc"]+sear_econ["cs1_intc"]+wpr_econ["cs1_intc"]
global_econ["cs2_intc"]=afr_econ["cs2_intc"]+amr_econ["cs2_intc"]+emr_econ["cs2_intc"]+eur_econ["cs2_intc"]+sear_econ["cs2_intc"]+wpr_econ["cs2_intc"]
global_econ["cs3_intc"]=afr_econ["cs3_intc"]+amr_econ["cs3_intc"]+emr_econ["cs3_intc"]+eur_econ["cs3_intc"]+sear_econ["cs3_intc"]+wpr_econ["cs3_intc"]

global_econ["bl_intc"]=afr_econ["bl_intc"]+amr_econ["bl_intc"]+emr_econ["bl_intc"]+eur_econ["bl_intc"]+sear_econ["bl_intc"]+wpr_econ["bl_intc"]
global_econ["s1_intc"]=afr_econ["s1_intc"]+amr_econ["s1_intc"]+emr_econ["s1_intc"]+eur_econ["s1_intc"]+sear_econ["s1_intc"]+wpr_econ["s1_intc"]
global_econ["s2_intc"]=afr_econ["s2_intc"]+amr_econ["s2_intc"]+emr_econ["s2_intc"]+eur_econ["s2_intc"]+sear_econ["s2_intc"]+wpr_econ["s2_intc"]
global_econ["s3_intc"]=afr_econ["s3_intc"]+amr_econ["s3_intc"]+emr_econ["s3_intc"]+eur_econ["s3_intc"]+sear_econ["s3_intc"]+wpr_econ["s3_intc"]
## Add medical costs
global_econ["cbl_medc"]=afr_econ["cbl_medc"]+amr_econ["cbl_medc"]+emr_econ["cbl_medc"]+eur_econ["cbl_medc"]+sear_econ["cbl_medc"]+wpr_econ["cbl_medc"]
global_econ["cs1_medc"]=afr_econ["cs1_medc"]+amr_econ["cs1_medc"]+emr_econ["cs1_medc"]+eur_econ["cs1_medc"]+sear_econ["cs1_medc"]+wpr_econ["cs1_medc"]
global_econ["cs2_medc"]=afr_econ["cs2_medc"]+amr_econ["cs2_medc"]+emr_econ["cs2_medc"]+eur_econ["cs2_medc"]+sear_econ["cs2_medc"]+wpr_econ["cs2_medc"]
global_econ["cs3_medc"]=afr_econ["cs3_medc"]+amr_econ["cs3_medc"]+emr_econ["cs3_medc"]+eur_econ["cs3_medc"]+sear_econ["cs3_medc"]+wpr_econ["cs3_medc"]

global_econ["bl_medc"]=afr_econ["bl_medc"]+amr_econ["bl_medc"]+emr_econ["bl_medc"]+eur_econ["bl_medc"]+sear_econ["bl_medc"]+wpr_econ["bl_medc"]
global_econ["s1_medc"]=afr_econ["s1_medc"]+amr_econ["s1_medc"]+emr_econ["s1_medc"]+eur_econ["s1_medc"]+sear_econ["s1_medc"]+wpr_econ["s1_medc"]
global_econ["s2_medc"]=afr_econ["s2_medc"]+amr_econ["s2_medc"]+emr_econ["s2_medc"]+eur_econ["s2_medc"]+sear_econ["s2_medc"]+wpr_econ["s2_medc"]
global_econ["s3_medc"]=afr_econ["s3_medc"]+amr_econ["s3_medc"]+emr_econ["s3_medc"]+eur_econ["s3_medc"]+sear_econ["s3_medc"]+wpr_econ["s3_medc"]
## Add productivity Loss
global_econ["cbl_prod"]=afr_econ["cbl_prod"]+amr_econ["cbl_prod"]+emr_econ["cbl_prod"]+eur_econ["cbl_prod"]+sear_econ["cbl_prod"]+wpr_econ["cbl_prod"]
global_econ["cs1_prod"]=afr_econ["cs1_prod"]+amr_econ["cs1_prod"]+emr_econ["cs1_prod"]+eur_econ["cs1_prod"]+sear_econ["cs1_prod"]+wpr_econ["cs1_prod"]
global_econ["cs2_prod"]=afr_econ["cs2_prod"]+amr_econ["cs2_prod"]+emr_econ["cs2_prod"]+eur_econ["cs2_prod"]+sear_econ["cs2_prod"]+wpr_econ["cs2_prod"]
global_econ["cs3_prod"]=afr_econ["cs3_prod"]+amr_econ["cs3_prod"]+emr_econ["cs3_prod"]+eur_econ["cs3_prod"]+sear_econ["cs3_prod"]+wpr_econ["cs3_prod"]

global_econ["bl_prod"]=afr_econ["bl_prod"]+amr_econ["bl_prod"]+emr_econ["bl_prod"]+eur_econ["bl_prod"]+sear_econ["bl_prod"]+wpr_econ["bl_prod"]
global_econ["s1_prod"]=afr_econ["s1_prod"]+amr_econ["s1_prod"]+emr_econ["s1_prod"]+eur_econ["s1_prod"]+sear_econ["s1_prod"]+wpr_econ["s1_prod"]
global_econ["s2_prod"]=afr_econ["s2_prod"]+amr_econ["s2_prod"]+emr_econ["s2_prod"]+eur_econ["s2_prod"]+sear_econ["s2_prod"]+wpr_econ["s2_prod"]
global_econ["s3_prod"]=afr_econ["s3_prod"]+amr_econ["s3_prod"]+emr_econ["s3_prod"]+eur_econ["s3_prod"]+sear_econ["s3_prod"]+wpr_econ["s3_prod"]

global_econ["cbl_daly"]=afr_econ["cbl_daly"]+amr_econ["cbl_daly"]+emr_econ["cbl_daly"]+eur_econ["cbl_daly"]+sear_econ["cbl_daly"]+wpr_econ["cbl_daly"]
global_econ["cs1_daly"]=afr_econ["cs1_daly"]+amr_econ["cs1_daly"]+emr_econ["cs1_daly"]+eur_econ["cs1_daly"]+sear_econ["cs1_daly"]+wpr_econ["cs1_daly"]
global_econ["cs2_daly"]=afr_econ["cs2_daly"]+amr_econ["cs2_daly"]+emr_econ["cs2_daly"]+eur_econ["cs2_daly"]+sear_econ["cs2_daly"]+wpr_econ["cs2_daly"]
global_econ["cs3_daly"]=afr_econ["cs3_daly"]+amr_econ["cs3_daly"]+emr_econ["cs3_daly"]+eur_econ["cs3_daly"]+sear_econ["cs3_daly"]+wpr_econ["cs3_daly"]

global_econ["bl_daly"]=afr_econ["bl_daly"]+amr_econ["bl_daly"]+emr_econ["bl_daly"]+eur_econ["bl_daly"]+sear_econ["bl_daly"]+wpr_econ["bl_daly"]
global_econ["s1_daly"]=afr_econ["s1_daly"]+amr_econ["s1_daly"]+emr_econ["s1_daly"]+eur_econ["s1_daly"]+sear_econ["s1_daly"]+wpr_econ["s1_daly"]
global_econ["s2_daly"]=afr_econ["s2_daly"]+amr_econ["s2_daly"]+emr_econ["s2_daly"]+eur_econ["s2_daly"]+sear_econ["s2_daly"]+wpr_econ["s2_daly"]
global_econ["s3_daly"]=afr_econ["s3_daly"]+amr_econ["s3_daly"]+emr_econ["s3_daly"]+eur_econ["s3_daly"]+sear_econ["s3_daly"]+wpr_econ["s3_daly"]

global_econ["cbl_dirc"]=afr_econ["cbl_dirc"]+amr_econ["cbl_dirc"]+emr_econ["cbl_dirc"]+eur_econ["cbl_dirc"]+sear_econ["cbl_dirc"]+wpr_econ["cbl_dirc"]
global_econ["cs1_dirc"]=afr_econ["cs1_dirc"]+amr_econ["cs1_dirc"]+emr_econ["cs1_dirc"]+eur_econ["cs1_dirc"]+sear_econ["cs1_dirc"]+wpr_econ["cs1_dirc"]
global_econ["cs2_dirc"]=afr_econ["cs2_dirc"]+amr_econ["cs2_dirc"]+emr_econ["cs2_dirc"]+eur_econ["cs2_dirc"]+sear_econ["cs2_dirc"]+wpr_econ["cs2_dirc"]
global_econ["cs3_dirc"]=afr_econ["cs3_dirc"]+amr_econ["cs3_dirc"]+emr_econ["cs3_dirc"]+eur_econ["cs3_dirc"]+sear_econ["cs3_dirc"]+wpr_econ["cs3_dirc"]

global_econ["bl_dirc"]=afr_econ["bl_dirc"]+amr_econ["bl_dirc"]+emr_econ["bl_dirc"]+eur_econ["bl_dirc"]+sear_econ["bl_dirc"]+wpr_econ["bl_dirc"]
global_econ["s1_dirc"]=afr_econ["s1_dirc"]+amr_econ["s1_dirc"]+emr_econ["s1_dirc"]+eur_econ["s1_dirc"]+sear_econ["s1_dirc"]+wpr_econ["s1_dirc"]
global_econ["s2_dirc"]=afr_econ["s2_dirc"]+amr_econ["s2_dirc"]+emr_econ["s2_dirc"]+eur_econ["s2_dirc"]+sear_econ["s2_dirc"]+wpr_econ["s2_dirc"]
global_econ["s3_dirc"]=afr_econ["s3_dirc"]+amr_econ["s3_dirc"]+emr_econ["s3_dirc"]+eur_econ["s3_dirc"]+sear_econ["s3_dirc"]+wpr_econ["s3_dirc"]

global_econ["cs1_icer"]=-((afr_econ["cbl_cumdc"]+amr_econ["cbl_cumdc"]+emr_econ["cbl_cumdc"]+eur_econ["cbl_cumdc"]+sear_econ["cbl_cumdc"]+wpr_econ["cbl_cumdc"])-(afr_econ["cs1_cumdc"]+amr_econ["cs1_cumdc"]+emr_econ["cs1_cumdc"]+eur_econ["cs1_cumdc"]+sear_econ["cs1_cumdc"]+wpr_econ["cs1_cumdc"]))/((afr_econ["cbl_cdal"]+amr_econ["cbl_cdal"]+emr_econ["cbl_cdal"]+eur_econ["cbl_cdal"]+sear_econ["cbl_cdal"]+wpr_econ["cbl_cdal"])-(afr_econ["cs1_cdal"]+amr_econ["cs1_cdal"]+emr_econ["cs1_cdal"]+eur_econ["cs1_cdal"]+sear_econ["cs1_cdal"]+wpr_econ["cs1_cdal"]))
global_econ["cs2_icer"]=-((afr_econ["cbl_cumdc"]+amr_econ["cbl_cumdc"]+emr_econ["cbl_cumdc"]+eur_econ["cbl_cumdc"]+sear_econ["cbl_cumdc"]+wpr_econ["cbl_cumdc"])-(afr_econ["cs2_cumdc"]+amr_econ["cs2_cumdc"]+emr_econ["cs2_cumdc"]+eur_econ["cs2_cumdc"]+sear_econ["cs2_cumdc"]+wpr_econ["cs2_cumdc"]))/((afr_econ["cbl_cdal"]+amr_econ["cbl_cdal"]+emr_econ["cbl_cdal"]+eur_econ["cbl_cdal"]+sear_econ["cbl_cdal"]+wpr_econ["cbl_cdal"])-(afr_econ["cs2_cdal"]+amr_econ["cs2_cdal"]+emr_econ["cs2_cdal"]+eur_econ["cs2_cdal"]+sear_econ["cs2_cdal"]+wpr_econ["cs2_cdal"]))
global_econ["cs3_icer"]=-((afr_econ["cbl_cumdc"]+amr_econ["cbl_cumdc"]+emr_econ["cbl_cumdc"]+eur_econ["cbl_cumdc"]+sear_econ["cbl_cumdc"]+wpr_econ["cbl_cumdc"])-(afr_econ["cs3_cumdc"]+amr_econ["cs3_cumdc"]+emr_econ["cs3_cumdc"]+eur_econ["cs3_cumdc"]+sear_econ["cs3_cumdc"]+wpr_econ["cs3_cumdc"]))/((afr_econ["cbl_cdal"]+amr_econ["cbl_cdal"]+emr_econ["cbl_cdal"]+eur_econ["cbl_cdal"]+sear_econ["cbl_cdal"]+wpr_econ["cbl_cdal"])-(afr_econ["cs3_cdal"]+amr_econ["cs3_cdal"]+emr_econ["cs3_cdal"]+eur_econ["cs3_cdal"]+sear_econ["cs3_cdal"]+wpr_econ["cs3_cdal"]))

global_econ["s1_icer"]=-((afr_econ["bl_cumdc"]+amr_econ["bl_cumdc"]+emr_econ["bl_cumdc"]+eur_econ["bl_cumdc"]+sear_econ["bl_cumdc"]+wpr_econ["bl_cumdc"])-(afr_econ["s1_cumdc"]+amr_econ["s1_cumdc"]+emr_econ["s1_cumdc"]+eur_econ["s1_cumdc"]+sear_econ["s1_cumdc"]+wpr_econ["s1_cumdc"]))/((afr_econ["bl_cdal"]+amr_econ["bl_cdal"]+emr_econ["bl_cdal"]+eur_econ["bl_cdal"]+sear_econ["bl_cdal"]+wpr_econ["bl_cdal"])-(afr_econ["s1_cdal"]+amr_econ["s1_cdal"]+emr_econ["s1_cdal"]+eur_econ["s1_cdal"]+sear_econ["s1_cdal"]+wpr_econ["s1_cdal"]))
global_econ["s2_icer"]=-((afr_econ["bl_cumdc"]+amr_econ["bl_cumdc"]+emr_econ["bl_cumdc"]+eur_econ["bl_cumdc"]+sear_econ["bl_cumdc"]+wpr_econ["bl_cumdc"])-(afr_econ["s2_cumdc"]+amr_econ["s2_cumdc"]+emr_econ["s2_cumdc"]+eur_econ["s2_cumdc"]+sear_econ["s2_cumdc"]+wpr_econ["s2_cumdc"]))/((afr_econ["bl_cdal"]+amr_econ["bl_cdal"]+emr_econ["bl_cdal"]+eur_econ["bl_cdal"]+sear_econ["bl_cdal"]+wpr_econ["bl_cdal"])-(afr_econ["s2_cdal"]+amr_econ["s2_cdal"]+emr_econ["s2_cdal"]+eur_econ["s2_cdal"]+sear_econ["s2_cdal"]+wpr_econ["s2_cdal"]))
global_econ["s3_icer"]=-((afr_econ["bl_cumdc"]+amr_econ["bl_cumdc"]+emr_econ["bl_cumdc"]+eur_econ["bl_cumdc"]+sear_econ["bl_cumdc"]+wpr_econ["bl_cumdc"])-(afr_econ["s3_cumdc"]+amr_econ["s3_cumdc"]+emr_econ["s3_cumdc"]+eur_econ["s3_cumdc"]+sear_econ["s3_cumdc"]+wpr_econ["s3_cumdc"]))/((afr_econ["bl_cdal"]+amr_econ["bl_cdal"]+emr_econ["bl_cdal"]+eur_econ["bl_cdal"]+sear_econ["bl_cdal"]+wpr_econ["bl_cdal"])-(afr_econ["s3_cdal"]+amr_econ["s3_cdal"]+emr_econ["s3_cdal"]+eur_econ["s3_cdal"]+sear_econ["s3_cdal"]+wpr_econ["s3_cdal"]))

global_econ["cs1_neb"]=afr_econ["cs1_neb"]+amr_econ["cs1_neb"]+emr_econ["cs1_neb"]+eur_econ["cs1_neb"]+sear_econ["cs1_neb"]+wpr_econ["cs1_neb"]
global_econ["cs2_neb"]=afr_econ["cs2_neb"]+amr_econ["cs2_neb"]+emr_econ["cs2_neb"]+eur_econ["cs2_neb"]+sear_econ["cs2_neb"]+wpr_econ["cs2_neb"]
global_econ["cs3_neb"]=afr_econ["cs3_neb"]+amr_econ["cs3_neb"]+emr_econ["cs3_neb"]+eur_econ["cs3_neb"]+sear_econ["cs3_neb"]+wpr_econ["cs3_neb"]

global_econ["s1_neb"]=afr_econ["s1_neb"]+amr_econ["s1_neb"]+emr_econ["s1_neb"]+eur_econ["s1_neb"]+sear_econ["s1_neb"]+wpr_econ["s1_neb"]
global_econ["s2_neb"]=afr_econ["s2_neb"]+amr_econ["s2_neb"]+emr_econ["s2_neb"]+eur_econ["s2_neb"]+sear_econ["s2_neb"]+wpr_econ["s2_neb"]
global_econ["s3_neb"]=afr_econ["s3_neb"]+amr_econ["s3_neb"]+emr_econ["s3_neb"]+eur_econ["s3_neb"]+sear_econ["s3_neb"]+wpr_econ["s3_neb"]

econ_plots (global_econ, "global")

#%% Table 2 Creation

def summary_stats_median(bl, s1, s2, s3, econ,lb,ub):
    """Returns summary stats for each region, can be
    collated post-hoc for table 2"""

    import numpy as np

    error=[50,lb,ub] #median and bounds
    bl_tab_vars=["chb_inc", "hcc_inc", "mort", "bl_daly", "bl_intc", "bl_medc", "bl_prod"]
    s1_tab_vars=["chb_inc", "hcc_inc", "mort", "s1_daly", "s1_intc", "s1_medc", "s1_prod"]
    s2_tab_vars=["chb_inc", "hcc_inc", "mort", "s2_daly", "s2_intc", "s2_medc", "s2_prod"]
    s3_tab_vars=["chb_inc", "hcc_inc", "mort", "s3_daly", "s3_intc", "s3_medc", "s3_prod"]

    bl_sum=np.zeros((8,3))
    s1_sum=np.zeros((8,3))
    s2_sum=np.zeros((8,3))
    s3_sum=np.zeros((8,3))


    for idx,val in enumerate(bl_tab_vars):
        for i, err in enumerate(error):
            if idx==0:      #change this when re-running the model with CHB incidence added!
                bl_sum[idx,i]=np.percentile(np.sum(bl[val][32:61,:],axis=0),err)
            elif idx >0 and idx <3:
                bl_sum[idx,i]=np.percentile(np.sum(bl[val][32:61,:],axis=0),err)
            elif idx>=3:
                bl_sum[idx,i]=np.percentile(np.sum(econ[val][32:61,:],axis=0),err)

    for idx,val in enumerate(s1_tab_vars):
        for i, err in enumerate(error):
            if idx==0:      #change this when re-running the model with CHB incidence added!
                s1_sum[idx,i]=np.percentile(np.sum(s1[val][32:61,:],axis=0),err)
            elif idx >0 and idx <3:
                s1_sum[idx,i]=np.percentile(np.sum(s1[val][32:61,:],axis=0),err)
            elif idx>=3:
                s1_sum[idx,i]=np.percentile(np.sum(econ[val][32:61,:],axis=0),err)

    for idx,val in enumerate(s2_tab_vars):
        for i, err in enumerate(error):
            if idx==0:      #change this when re-running the model with CHB incidence added!
                s2_sum[idx,i]=np.percentile(np.sum(s2[val][32:61,:],axis=0),err)
            elif idx >0 and idx <3:
                s2_sum[idx,i]=np.percentile(np.sum(s2[val][32:61,:],axis=0),err)
            elif idx>=3:
                s2_sum[idx,i]=np.percentile(np.sum(econ[val][32:61,:],axis=0),err)

    for idx,val in enumerate(s3_tab_vars):
        for i, err in enumerate(error):
            if idx==0:      #change this when re-running the model with CHB incidence added!
                s3_sum[idx,i]=np.percentile(np.sum(s3[val][32:61,:],axis=0),err)
            elif idx >0 and idx <3:
                s3_sum[idx,i]=np.percentile(np.sum(s3[val][32:61,:],axis=0),err)
            elif idx>=3:
                s3_sum[idx,i]=np.percentile(np.sum(econ[val][32:61,:],axis=0),err)

    #Add years (with bounds) net economic benefit is reached
    s1_neb_med=np.zeros((109,2))
    s1_neb_med[:,0]=np.arange(1990, 2099, 1)
    s1_neb_med[:,1]=np.percentile(econ["s1_neb"],50, axis=1)
    year=[]
    for idx,val in enumerate(s1_neb_med[:,1]):
        if val >0:
            year.append(s1_neb_med[idx,0])
    s1_sum[7,0]=min(year)

    s1_neb_lb=np.zeros((109,2))
    s1_neb_lb[:,0]=np.arange(1990, 2099, 1)
    s1_neb_lb[:,1]=np.percentile(econ["s1_neb"],97.5, axis=1)
    year=[]
    for idx,val in enumerate(s1_neb_lb[:,1]):
        if val >0:
            year.append(s1_neb_lb[idx,0])
    s1_sum[7,1]=min(year)

    s1_neb_ub=np.zeros((109,2))
    s1_neb_ub[:,0]=np.arange(1990, 2099, 1)
    s1_neb_ub[:,1]=np.percentile(econ["s1_neb"],2.5, axis=1)
    year=[]
    for idx,val in enumerate(s1_neb_ub[:,1]):
        if val >0:
            year.append(s1_neb_ub[idx,0])
    s1_sum[7,2]=min(year)


    s2_neb_med=np.zeros((109,2))
    s2_neb_med[:,0]=np.arange(1990, 2099, 1)
    s2_neb_med[:,1]=np.percentile(econ["s2_neb"],50, axis=1)
    year=[]
    for idx,val in enumerate(s2_neb_med[:,1]):
        if val >0:
            year.append(s2_neb_med[idx,0])
    s2_sum[7,0]=min(year)

    s2_neb_lb=np.zeros((109,2))
    s2_neb_lb[:,0]=np.arange(1990, 2099, 1)
    s2_neb_lb[:,1]=np.percentile(econ["s2_neb"],97.5, axis=1)
    year=[]
    for idx,val in enumerate(s2_neb_lb[:,1]):
        if val >0:
            year.append(s2_neb_lb[idx,0])
    s2_sum[7,1]=min(year)

    s2_neb_ub=np.zeros((109,2))
    s2_neb_ub[:,0]=np.arange(1990, 2099, 1)
    s2_neb_ub[:,1]=np.percentile(econ["s2_neb"],2.5, axis=1)
    year=[]
    for idx,val in enumerate(s2_neb_ub[:,1]):
        if val >0:
            year.append(s2_neb_ub[idx,0])
    s2_sum[7,2]=min(year)


    s3_neb_med=np.zeros((109,2))
    s3_neb_med[:,0]=np.arange(1990, 2099, 1)
    s3_neb_med[:,1]=np.percentile(econ["s3_neb"],50, axis=1)
    year=[]
    for idx,val in enumerate(s3_neb_med[:,1]):
        if val >0:
            year.append(s3_neb_med[idx,0])
    s3_sum[7,0]=min(year)


    s3_neb_lb=np.zeros((109,2))
    s3_neb_lb[:,0]=np.arange(1990, 2099, 1)
    s3_neb_lb[:,1]=np.percentile(econ["s3_neb"],97.5, axis=1)
    year=[]
    for idx,val in enumerate(s3_neb_lb[:,1]):
        if val >0:
            year.append(s3_neb_lb[idx,0])
    s3_sum[7,1]=min(year)

    s3_neb_ub=np.zeros((109,2))
    s3_neb_ub[:,0]=np.arange(1990, 2099, 1)
    s3_neb_ub[:,1]=np.percentile(econ["s3_neb"],2.5, axis=1)
    year=[]
    for idx,val in enumerate(s3_neb_ub[:,1]):
        if val >0:
            year.append(s3_neb_ub[idx,0])
    s3_sum[7,2]=min(year)

    t2_sum=np.concatenate([bl_sum,s1_sum, s2_sum, s3_sum], axis=1)

    return t2_sum

tab_2=np.zeros((7,8,12))
tab_2[0,:,:]=summary_stats_median(afr_bl_runs, afr_s1_runs, afr_s2_runs, afr_s3_runs, afr_econ,2.5,97.5)
tab_2[1,:,:]=summary_stats_median(amr_bl_runs, amr_s1_runs, amr_s2_runs, amr_s3_runs, amr_econ,2.5,97.5)
tab_2[2,:,:]=summary_stats_median(emr_bl_runs, emr_s1_runs, emr_s2_runs, emr_s3_runs, emr_econ,2.5,97.5)
tab_2[3,:,:]=summary_stats_median(eur_bl_runs, eur_s1_runs, eur_s2_runs, eur_s3_runs, eur_econ,2.5,97.5)
tab_2[4,:,:]=summary_stats_median(sear_bl_runs, sear_s1_runs, sear_s2_runs, sear_s3_runs, sear_econ,2.5,97.5)
tab_2[5,:,:]=summary_stats_median(wpr_bl_runs, wpr_s1_runs, wpr_s2_runs, wpr_s3_runs, wpr_econ,2.5,97.5)
tab_2[6,:,:]=summary_stats_median(global_bl_runs, global_s1_runs, global_s2_runs, global_s3_runs, global_econ,2.5,97.5)

##Replace medians with central estimates

## AFR
year=np.arange(1990, 2099,1)

tab_2[0,0,0]=np.sum(afr_bl_cent["chb_inc"][32:61])
tab_2[0,1,0]=np.sum(afr_bl_cent["hcc_inc"][32:61])
tab_2[0,2,0]=np.sum(afr_bl_cent["mort"][32:61])
tab_2[0,3,0]=np.sum(afr_econ["cbl_daly"][32:61])
tab_2[0,4,0]=np.sum(afr_econ["cbl_intc"][32:61])
tab_2[0,5,0]=np.sum(afr_econ["cbl_medc"][32:61])
tab_2[0,6,0]=np.sum(afr_econ["cbl_prod"][32:61])

tab_2[0,0,3]=np.sum(afr_s1_cent["chb_inc"][32:61])
tab_2[0,1,3]=np.sum(afr_s1_cent["hcc_inc"][32:61])
tab_2[0,2,3]=np.sum(afr_s1_cent["mort"][32:61])
tab_2[0,3,3]=np.sum(afr_econ["cs1_daly"][32:61])
tab_2[0,4,3]=np.sum(afr_econ["cs1_intc"][32:61])
tab_2[0,5,3]=np.sum(afr_econ["cs1_medc"][32:61])
tab_2[0,6,3]=np.sum(afr_econ["cs1_prod"][32:61])
year_idx=[]
for idx,val in enumerate(afr_econ["cs1_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[0,7,3]=year[yneb_i]

tab_2[0,0,6]=np.sum(afr_s2_cent["chb_inc"][32:61])
tab_2[0,1,6]=np.sum(afr_s2_cent["hcc_inc"][32:61])
tab_2[0,2,6]=np.sum(afr_s2_cent["mort"][32:61])
tab_2[0,3,6]=np.sum(afr_econ["cs2_daly"][32:61])
tab_2[0,4,6]=np.sum(afr_econ["cs2_intc"][32:61])
tab_2[0,5,6]=np.sum(afr_econ["cs2_medc"][32:61])
tab_2[0,6,6]=np.sum(afr_econ["cs2_prod"][32:61])
year_idx=[]
for idx,val in enumerate(afr_econ["cs2_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[0,7,6]=year[yneb_i]

tab_2[0,0,9]=np.sum(afr_s3_cent["chb_inc"][32:61])
tab_2[0,1,9]=np.sum(afr_s3_cent["hcc_inc"][32:61])
tab_2[0,2,9]=np.sum(afr_s3_cent["mort"][32:61])
tab_2[0,3,9]=np.sum(afr_econ["cs3_daly"][32:61])
tab_2[0,4,9]=np.sum(afr_econ["cs3_intc"][32:61])
tab_2[0,5,9]=np.sum(afr_econ["cs3_medc"][32:61])
tab_2[0,6,9]=np.sum(afr_econ["cs3_prod"][32:61])
year_idx=[]
for idx,val in enumerate(afr_econ["cs3_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[0,7,9]=year[yneb_i]


## AMR
tab_2[1,0,0]=np.sum(amr_bl_cent["chb_inc"][32:61])
tab_2[1,1,0]=np.sum(amr_bl_cent["hcc_inc"][32:61])
tab_2[1,2,0]=np.sum(amr_bl_cent["mort"][32:61])
tab_2[1,3,0]=np.sum(amr_econ["cbl_daly"][32:61])
tab_2[1,4,0]=np.sum(amr_econ["cbl_intc"][32:61])
tab_2[1,5,0]=np.sum(amr_econ["cbl_medc"][32:61])
tab_2[1,6,0]=np.sum(amr_econ["cbl_prod"][32:61])

tab_2[1,0,3]=np.sum(amr_s1_cent["chb_inc"][32:61])
tab_2[1,1,3]=np.sum(amr_s1_cent["hcc_inc"][32:61])
tab_2[1,2,3]=np.sum(amr_s1_cent["mort"][32:61])
tab_2[1,3,3]=np.sum(amr_econ["cs1_daly"][32:61])
tab_2[1,4,3]=np.sum(amr_econ["cs1_intc"][32:61])
tab_2[1,5,3]=np.sum(amr_econ["cs1_medc"][32:61])
tab_2[1,6,3]=np.sum(amr_econ["cs1_prod"][32:61])
year_idx=[]
for idx,val in enumerate(amr_econ["cs1_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[1,7,3]=year[yneb_i]

tab_2[1,0,6]=np.sum(amr_s2_cent["chb_inc"][32:61])
tab_2[1,1,6]=np.sum(amr_s2_cent["hcc_inc"][32:61])
tab_2[1,2,6]=np.sum(amr_s2_cent["mort"][32:61])
tab_2[1,3,6]=np.sum(amr_econ["cs2_daly"][32:61])
tab_2[1,4,6]=np.sum(amr_econ["cs2_intc"][32:61])
tab_2[1,5,6]=np.sum(amr_econ["cs2_medc"][32:61])
tab_2[1,6,6]=np.sum(amr_econ["cs2_prod"][32:61])
year_idx=[]
for idx,val in enumerate(amr_econ["cs2_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[1,7,6]=year[yneb_i]

tab_2[1,0,9]=np.sum(amr_s3_cent["chb_inc"][32:61])
tab_2[1,1,9]=np.sum(amr_s3_cent["hcc_inc"][32:61])
tab_2[1,2,9]=np.sum(amr_s3_cent["mort"][32:61])
tab_2[1,3,9]=np.sum(amr_econ["cs3_daly"][32:61])
tab_2[1,4,9]=np.sum(amr_econ["cs3_intc"][32:61])
tab_2[1,5,9]=np.sum(amr_econ["cs3_medc"][32:61])
tab_2[1,6,9]=np.sum(amr_econ["cs3_prod"][32:61])
year_idx=[]
for idx,val in enumerate(amr_econ["cs3_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[1,7,9]=year[yneb_i]

## EMR
tab_2[2,0,0]=np.sum(emr_bl_cent["chb_inc"][32:61])
tab_2[2,1,0]=np.sum(emr_bl_cent["hcc_inc"][32:61])
tab_2[2,2,0]=np.sum(emr_bl_cent["mort"][32:61])
tab_2[2,3,0]=np.sum(emr_econ["cbl_daly"][32:61])
tab_2[2,4,0]=np.sum(emr_econ["cbl_intc"][32:61])
tab_2[2,5,0]=np.sum(emr_econ["cbl_medc"][32:61])
tab_2[2,6,0]=np.sum(emr_econ["cbl_prod"][32:61])

tab_2[2,0,3]=np.sum(emr_s1_cent["chb_inc"][32:61])
tab_2[2,1,3]=np.sum(emr_s1_cent["hcc_inc"][32:61])
tab_2[2,2,3]=np.sum(emr_s1_cent["mort"][32:61])
tab_2[2,3,3]=np.sum(emr_econ["cs1_daly"][32:61])
tab_2[2,4,3]=np.sum(emr_econ["cs1_intc"][32:61])
tab_2[2,5,3]=np.sum(emr_econ["cs1_medc"][32:61])
tab_2[2,6,3]=np.sum(emr_econ["cs1_prod"][32:61])
year_idx=[]
for idx,val in enumerate(emr_econ["cs1_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[2,7,3]=year[yneb_i]

tab_2[2,0,6]=np.sum(emr_s2_cent["chb_inc"][32:61])
tab_2[2,1,6]=np.sum(emr_s2_cent["hcc_inc"][32:61])
tab_2[2,2,6]=np.sum(emr_s2_cent["mort"][32:61])
tab_2[2,3,6]=np.sum(emr_econ["cs2_daly"][32:61])
tab_2[2,4,6]=np.sum(emr_econ["cs2_intc"][32:61])
tab_2[2,5,6]=np.sum(emr_econ["cs2_medc"][32:61])
tab_2[2,6,6]=np.sum(emr_econ["cs2_prod"][32:61])
year_idx=[]
for idx,val in enumerate(emr_econ["cs2_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[2,7,6]=year[yneb_i]

tab_2[2,0,9]=np.sum(emr_s3_cent["chb_inc"][32:61])
tab_2[2,1,9]=np.sum(emr_s3_cent["hcc_inc"][32:61])
tab_2[2,2,9]=np.sum(emr_s3_cent["mort"][32:61])
tab_2[2,3,9]=np.sum(emr_econ["cs3_daly"][32:61])
tab_2[2,4,9]=np.sum(emr_econ["cs3_intc"][32:61])
tab_2[2,5,9]=np.sum(emr_econ["cs3_medc"][32:61])
tab_2[2,6,9]=np.sum(emr_econ["cs3_prod"][32:61])
year_idx=[]
for idx,val in enumerate(emr_econ["cs3_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[2,7,9]=year[yneb_i]

## EUR
tab_2[3,1,0]=np.sum(eur_bl_cent["hcc_inc"][32:61])
tab_2[3,2,0]=np.sum(eur_bl_cent["mort"][32:61])
tab_2[3,3,0]=np.sum(eur_econ["cbl_daly"][32:61])
tab_2[3,4,0]=np.sum(eur_econ["cbl_intc"][32:61])
tab_2[3,5,0]=np.sum(eur_econ["cbl_medc"][32:61])
tab_2[3,6,0]=np.sum(eur_econ["cbl_prod"][32:61])

tab_2[3,0,3]=np.sum(eur_s1_cent["chb_inc"][32:61])
tab_2[3,1,3]=np.sum(eur_s1_cent["hcc_inc"][32:61])
tab_2[3,2,3]=np.sum(eur_s1_cent["mort"][32:61])
tab_2[3,3,3]=np.sum(eur_econ["cs1_daly"][32:61])
tab_2[3,4,3]=np.sum(eur_econ["cs1_intc"][32:61])
tab_2[3,5,3]=np.sum(eur_econ["cs1_medc"][32:61])
tab_2[3,6,3]=np.sum(eur_econ["cs1_prod"][32:61])
year_idx=[]
for idx,val in enumerate(eur_econ["cs1_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[3,7,3]=year[yneb_i]

tab_2[3,0,6]=np.sum(eur_s2_cent["chb_inc"][32:61])
tab_2[3,1,6]=np.sum(eur_s2_cent["hcc_inc"][32:61])
tab_2[3,2,6]=np.sum(eur_s2_cent["mort"][32:61])
tab_2[3,3,6]=np.sum(eur_econ["cs2_daly"][32:61])
tab_2[3,4,6]=np.sum(eur_econ["cs2_intc"][32:61])
tab_2[3,5,6]=np.sum(eur_econ["cs2_medc"][32:61])
tab_2[3,6,6]=np.sum(eur_econ["cs2_prod"][32:61])
year_idx=[]
for idx,val in enumerate(eur_econ["cs2_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[3,7,6]=year[yneb_i]

tab_2[3,0,9]=np.sum(eur_s3_cent["chb_inc"][32:61])
tab_2[3,1,9]=np.sum(eur_s3_cent["hcc_inc"][32:61])
tab_2[3,2,9]=np.sum(eur_s3_cent["mort"][32:61])
tab_2[3,3,9]=np.sum(eur_econ["cs3_daly"][32:61])
tab_2[3,4,9]=np.sum(eur_econ["cs3_intc"][32:61])
tab_2[3,5,9]=np.sum(eur_econ["cs3_medc"][32:61])
tab_2[3,6,9]=np.sum(eur_econ["cs3_prod"][32:61])
year_idx=[]
for idx,val in enumerate(eur_econ["cs3_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[3,7,9]=year[yneb_i]

## SEAR
tab_2[4,1,0]=np.sum(sear_bl_cent["hcc_inc"][32:61])
tab_2[4,2,0]=np.sum(sear_bl_cent["mort"][32:61])
tab_2[4,3,0]=np.sum(sear_econ["cbl_daly"][32:61])
tab_2[4,4,0]=np.sum(sear_econ["cbl_intc"][32:61])
tab_2[4,5,0]=np.sum(sear_econ["cbl_medc"][32:61])
tab_2[4,6,0]=np.sum(sear_econ["cbl_prod"][32:61])

tab_2[4,0,3]=np.sum(sear_s1_cent["chb_inc"][32:61])
tab_2[4,1,3]=np.sum(sear_s1_cent["hcc_inc"][32:61])
tab_2[4,2,3]=np.sum(sear_s1_cent["mort"][32:61])
tab_2[4,3,3]=np.sum(sear_econ["cs1_daly"][32:61])
tab_2[4,4,3]=np.sum(sear_econ["cs1_intc"][32:61])
tab_2[4,5,3]=np.sum(sear_econ["cs1_medc"][32:61])
tab_2[4,6,3]=np.sum(sear_econ["cs1_prod"][32:61])
year_idx=[]
for idx,val in enumerate(sear_econ["cs1_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[4,7,3]=year[yneb_i]

tab_2[4,0,6]=np.sum(sear_s2_cent["chb_inc"][32:61])
tab_2[4,1,6]=np.sum(sear_s2_cent["hcc_inc"][32:61])
tab_2[4,2,6]=np.sum(sear_s2_cent["mort"][32:61])
tab_2[4,3,6]=np.sum(sear_econ["cs2_daly"][32:61])
tab_2[4,4,6]=np.sum(sear_econ["cs2_intc"][32:61])
tab_2[4,5,6]=np.sum(sear_econ["cs2_medc"][32:61])
tab_2[4,6,6]=np.sum(sear_econ["cs2_prod"][32:61])
year_idx=[]
for idx,val in enumerate(sear_econ["cs2_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[4,7,6]=year[yneb_i]

tab_2[4,0,9]=np.sum(sear_s3_cent["chb_inc"][32:61])
tab_2[4,1,9]=np.sum(sear_s3_cent["hcc_inc"][32:61])
tab_2[4,2,9]=np.sum(sear_s3_cent["mort"][32:61])
tab_2[4,3,9]=np.sum(sear_econ["cs3_daly"][32:61])
tab_2[4,4,9]=np.sum(sear_econ["cs3_intc"][32:61])
tab_2[4,5,9]=np.sum(sear_econ["cs3_medc"][32:61])
tab_2[4,6,9]=np.sum(sear_econ["cs3_prod"][32:61])
year_idx=[]
for idx,val in enumerate(sear_econ["cs3_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[4,7,9]=year[yneb_i]

## WPR
tab_2[5,1,0]=np.sum(wpr_bl_cent["hcc_inc"][32:61])
tab_2[5,2,0]=np.sum(wpr_bl_cent["mort"][32:61])
tab_2[5,3,0]=np.sum(wpr_econ["cbl_daly"][32:61])
tab_2[5,4,0]=np.sum(wpr_econ["cbl_intc"][32:61])
tab_2[5,5,0]=np.sum(wpr_econ["cbl_medc"][32:61])
tab_2[5,6,0]=np.sum(wpr_econ["cbl_prod"][32:61])

tab_2[5,0,3]=np.sum(wpr_s1_cent["chb_inc"][32:61])
tab_2[5,1,3]=np.sum(wpr_s1_cent["hcc_inc"][32:61])
tab_2[5,2,3]=np.sum(wpr_s1_cent["mort"][32:61])
tab_2[5,3,3]=np.sum(wpr_econ["cs1_daly"][32:61])
tab_2[5,4,3]=np.sum(wpr_econ["cs1_intc"][32:61])
tab_2[5,5,3]=np.sum(wpr_econ["cs1_medc"][32:61])
tab_2[5,6,3]=np.sum(wpr_econ["cs1_prod"][32:61])
year_idx=[]
for idx,val in enumerate(wpr_econ["cs1_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[5,7,3]=year[yneb_i]

tab_2[5,0,6]=np.sum(wpr_s2_cent["chb_inc"][32:61])
tab_2[5,1,6]=np.sum(wpr_s2_cent["hcc_inc"][32:61])
tab_2[5,2,6]=np.sum(wpr_s2_cent["mort"][32:61])
tab_2[5,3,6]=np.sum(wpr_econ["cs2_daly"][32:61])
tab_2[5,4,6]=np.sum(wpr_econ["cs2_intc"][32:61])
tab_2[5,5,6]=np.sum(wpr_econ["cs2_medc"][32:61])
tab_2[5,6,6]=np.sum(wpr_econ["cs2_prod"][32:61])
year_idx=[]
for idx,val in enumerate(wpr_econ["cs2_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[5,7,6]=year[yneb_i]

tab_2[5,0,9]=np.sum(wpr_s3_cent["chb_inc"][32:61])
tab_2[5,1,9]=np.sum(wpr_s3_cent["hcc_inc"][32:61])
tab_2[5,2,9]=np.sum(wpr_s3_cent["mort"][32:61])
tab_2[5,3,9]=np.sum(wpr_econ["cs3_daly"][32:61])
tab_2[5,4,9]=np.sum(wpr_econ["cs3_intc"][32:61])
tab_2[5,5,9]=np.sum(wpr_econ["cs3_medc"][32:61])
tab_2[5,6,9]=np.sum(wpr_econ["cs3_prod"][32:61])
year_idx=[]
for idx,val in enumerate(wpr_econ["cs3_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[5,7,9]=year[yneb_i]

## Global
tab_2[6,1,0]=np.sum(global_bl_cent["hcc_inc"][32:61])
tab_2[6,2,0]=np.sum(global_bl_cent["mort"][32:61])
tab_2[6,3,0]=np.sum(global_econ["cbl_daly"][32:61])
tab_2[6,4,0]=np.sum(global_econ["cbl_intc"][32:61])
tab_2[6,5,0]=np.sum(global_econ["cbl_medc"][32:61])
tab_2[6,6,0]=np.sum(global_econ["cbl_prod"][32:61])

tab_2[6,0,3]=np.sum(global_s1_cent["chb_inc"][32:61])
tab_2[6,1,3]=np.sum(global_s1_cent["hcc_inc"][32:61])
tab_2[6,2,3]=np.sum(global_s1_cent["mort"][32:61])
tab_2[6,3,3]=np.sum(global_econ["cs1_daly"][32:61])
tab_2[6,4,3]=np.sum(global_econ["cs1_intc"][32:61])
tab_2[6,5,3]=np.sum(global_econ["cs1_medc"][32:61])
tab_2[6,6,3]=np.sum(global_econ["cs1_prod"][32:61])
year_idx=[]
for idx,val in enumerate(global_econ["cs1_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[6,7,3]=year[yneb_i]

tab_2[6,0,6]=np.sum(global_s2_cent["chb_inc"][32:61])
tab_2[6,1,6]=np.sum(global_s2_cent["hcc_inc"][32:61])
tab_2[6,2,6]=np.sum(global_s2_cent["mort"][32:61])
tab_2[6,3,6]=np.sum(global_econ["cs2_daly"][32:61])
tab_2[6,4,6]=np.sum(global_econ["cs2_intc"][32:61])
tab_2[6,5,6]=np.sum(global_econ["cs2_medc"][32:61])
tab_2[6,6,6]=np.sum(global_econ["cs2_prod"][32:61])
year_idx=[]
for idx,val in enumerate(global_econ["cs2_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[6,7,6]=year[yneb_i]

tab_2[6,0,9]=np.sum(global_s3_cent["chb_inc"][32:61])
tab_2[6,1,9]=np.sum(global_s3_cent["hcc_inc"][32:61])
tab_2[6,2,9]=np.sum(global_s3_cent["mort"][32:61])
tab_2[6,3,9]=np.sum(global_econ["cs3_daly"][32:61])
tab_2[6,4,9]=np.sum(global_econ["cs3_intc"][32:61])
tab_2[6,5,9]=np.sum(global_econ["cs3_medc"][32:61])
tab_2[6,6,9]=np.sum(global_econ["cs3_prod"][32:61])
year_idx=[]
for idx,val in enumerate(global_econ["cs3_neb"]):
    if val>0:
        year_idx.append(idx)
yneb_i=min(year_idx)
tab_2[6,7,9]=year[yneb_i]


## Save Table 2 as a spreadsheet
writer=pd.ExcelWriter("Table 2_central.xlsx", engine="xlsxwriter")

region=["AFR", "AMR", "EMR", "EUR", "SEAR", "WPR", "Global"]

for i in range(len(reg)):
    df=pd.DataFrame(tab_2[i,:,:])
    df.to_excel(writer, sheet_name=region[i])

writer.save()

## Global Outcomes for parameter scenarios
global_mtctw_s1={}
global_mtctb_s1={}
global_mortw_s1={}
global_mortb_s1={}
global_trtw_s1={}
global_trtb_s1={}

for key, val in afr_mtctw_s1.items():
    global_mtctw_s1[key]=afr_mtctw_s1[key]+amr_mtctw_s1[key]+emr_mtctw_s1[key]+eur_mtctw_s1[key]+sear_mtctw_s1[key]+wpr_mtctw_s1[key]
    global_mtctb_s1[key]=afr_mtctb_s1[key]+amr_mtctb_s1[key]+emr_mtctb_s1[key]+eur_mtctb_s1[key]+sear_mtctb_s1[key]+wpr_mtctb_s1[key]
    global_mortw_s1[key]=afr_mortw_s1[key]+amr_mortw_s1[key]+emr_mortw_s1[key]+eur_mortw_s1[key]+sear_mortw_s1[key]+wpr_mortw_s1[key]
    global_mortb_s1[key]=afr_mortb_s1[key]+amr_mortb_s1[key]+emr_mortb_s1[key]+eur_mortb_s1[key]+sear_mortb_s1[key]+wpr_mortb_s1[key]
    global_trtw_s1[key]=afr_trtw_s1[key]+amr_trtw_s1[key]+emr_trtw_s1[key]+eur_trtw_s1[key]+sear_trtw_s1[key]+wpr_trtw_s1[key]
    global_trtb_s1[key]=afr_trtb_s1[key]+amr_trtb_s1[key]+emr_trtb_s1[key]+eur_trtb_s1[key]+sear_trtb_s1[key]+wpr_trtb_s1[key]

global_mtctw_s1["icer"]=-(global_mtctw_s1["cdirc_bl"]-global_mtctw_s1["cdirc_s1"])/(global_mtctw_s1["cdaly_bl"]-global_mtctw_s1["cdaly_s1"])
global_mtctb_s1["icer"]=-(global_mtctb_s1["cdirc_bl"]-global_mtctb_s1["cdirc_s1"])/(global_mtctb_s1["cdaly_bl"]-global_mtctb_s1["cdaly_s1"])
global_mortw_s1["icer"]=-(global_mortw_s1["cdirc_bl"]-global_mortw_s1["cdirc_s1"])/(global_mortw_s1["cdaly_bl"]-global_mortw_s1["cdaly_s1"])
global_mortb_s1["icer"]=-(global_mortb_s1["cdirc_bl"]-global_mortb_s1["cdirc_s1"])/(global_mortb_s1["cdaly_bl"]-global_mortb_s1["cdaly_s1"])
global_trtw_s1["icer"]=-(global_trtw_s1["cdirc_bl"]-global_trtw_s1["cdirc_s1"])/(global_trtw_s1["cdaly_bl"]-global_trtw_s1["cdaly_s1"])
global_trtb_s1["icer"]=-(global_trtb_s1["cdirc_bl"]-global_trtb_s1["cdirc_s1"])/(global_trtb_s1["cdaly_bl"]-global_trtb_s1["cdaly_s1"])

afr_econ_sens= econsens_econ(afr_bl_cent, afr_s1_cent, afr_econ, "AFR", "cost and agg data/costs.xlsx", 0.03)
amr_econ_sens= econsens_econ(amr_bl_cent, amr_s1_cent, amr_econ, "AMR", "cost and agg data/costs.xlsx", 0.03)
emr_econ_sens= econsens_econ(emr_bl_cent, emr_s1_cent, emr_econ, "EMR", "cost and agg data/costs.xlsx", 0.03)
eur_econ_sens= econsens_econ(eur_bl_cent, eur_s1_cent, eur_econ, "EUR", "cost and agg data/costs.xlsx", 0.03)
sear_econ_sens= econsens_econ(sear_bl_cent, sear_s1_cent, sear_econ, "SEAR", "cost and agg data/costs.xlsx", 0.03)
wpr_econ_sens= econsens_econ(wpr_bl_cent, wpr_s1_cent, wpr_econ, "WPR", "cost and agg data/costs.xlsx", 0.03)

## Global Outcomes for economic scenarios
global_econ_sens={}

for key, val in afr_econ_sens.items():
    global_econ_sens[key]=afr_econ_sens[key]+amr_econ_sens[key]+emr_econ_sens[key]+eur_econ_sens[key]+sear_econ_sens[key]+wpr_econ_sens[key]