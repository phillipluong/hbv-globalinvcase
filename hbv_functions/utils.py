import atomica as at
import functools
import hbv_functions as hbv
import numpy as np
import pandas as pd
import seaborn as sns
import os
from collections import defaultdict
from sklearn.metrics import mean_absolute_error, mean_squared_error, mean_absolute_percentage_error

import matplotlib as mpl
mpl.use('TkAgg')

def model_results(F, region, data, calib, res_name, runs):

    """ Runs the model and simulations, returns central and sampled data sets for analysis

    :param F:        Project Framework
    :param region:   Type of region the model we're looking at (str)
    :param data:     Project databook name
    :param calib:    Project Calibration
    :param res_name: Name of result (str)
    :param runs:     Number of runs in the model
    :returns:        store_runs, central_est
    """

    import atomica as at
    import numpy as np

    ## Run the model
    P=at.Project(framework=F, databook=f"applications/{region}/"+data, sim_start=1990, sim_end=2099, sim_dt=0.25, do_run=False)
    cal=P.make_parset()
    cal.load_calibration(f"applications/{region}/"+calib)

    # Central Estimate (can be used if needed, but median of ~100 runs converges very well [i.e., no meaningful difference])
    res=P.run_sim(parset=cal, result_name = res_name) # can also expand code here to check calibrations if needed

    #Probabilistic sensitivity analysis
    np.random.seed(25012023)
    psa=P.parsets[0]
    # psa_sample=psa.sample()
    psa_res=P.run_sampled_sims(cal, n_samples=runs)

    #Initiate dictionary for required results
    central_est={}
    store_runs={}
    out_res=[{"pop":["alive", "total", "sum"]}, #total population size
             {"pop_u5":["alive", {"Under 5":["0-4M", "0-4F"]}, "sum"]}, #total under 5 population size
             {"prev":["prev", "total", "weighted"]}, #total population prevalence
             {"chb_pop":["chb_pop", "total", "weighted"]}, #total population living with CHB
             {"prev_u5":["prev", {"Under 5":["0-4M", "0-4F"]}, "weighted"]}, #under 5y prevalence
             {"mort":[":dd_hbv", "total", "sum"]}, #Total HBV mortality
             {"hcc_inc": ["flw_hcc", "total", "sum"]}, #HCC incidence
             {"chb_inc": ["tot_inc", "total", "sum"]}, #CHB incidence
             {"hbe_preg": ["eag_ott", "15-49F", "sum"]}, #HBeAg prevalence in pregnant women
             {"yld": ["yld", "total", "sum"]}, #years lost to disability
             {"births": ["b_rate",{"Under 5":["0-4M", "0-4F"]}, "sum" ]}, #total births
             {"bd_cov": ["bd", {"Under 5":["0-4M", "0-4F"]}, "weighted"]}, #birth dose coverage
             {"hb3_cov": ["hb3", {"Under 5":["0-4M", "0-4F"]}, "weighted"]}, #hb3 coverage
             {"dx_rate":["tot_dx", "total", "sum"]}, #annual (incident) diagnoses
             {"tx_cov": ["treat", "total", "sum"]}, #total treatment coverage
             {"dm_dx": [[{"dm_dx":"it_dx+icl_dx+ict_dx+ie_dx+cc_dx+dc_dx+hcc_dx"}], "total", "sum"]}, #disease management, no treatment (for costs)
             {"dm_tx":[[{"dm_tx":"icl_tx+ict_tx+ie_tx+cc_tx+dc_tx+hcc_tx"}], "total", "sum"]}, #disease management, treatment (for costs)
             {"pop_hcc": [[{"pop_hcc":"cc_dx+cc_tx+dc_dx+dc_tx"}], "total", "sum"]}, #population HCC surveillance (for costs)
             {"tgt_hcc": [[{"tgt_hcc":"it_dx+icl_dx+icl_tx+ict_dx+ict_tx+ie_dx+ie_tx"}], {"50+":["50-69M", "50-69F", "70+M", "70+F"]}, "sum"]}, #targeted HCC surveillance (for costs)
             {"tgt_hcc_b": [[{"tgt_hcc":"it_dx+icl_dx+icl_tx+ict_dx+ict_tx+ie_dx+ie_tx"}], {"40+":["15-49M", "15-49F"]}, "sum"]}, #targeted HCC surveillance (for costs (40-49))
             {"hsp_tx": [[{"hsp_tx":"dc_tx+hcc_tx"}], "total", "sum"]}, #hospitalizations among treated (for costs)
             {"hsp_utx":[[{"hsp_utx":"dc+dc_dx+hcc+hcc_dx"}], "total", "sum"]},#hospitalization among untreated (for costs)
             {"dx_prop": ["diag_cov", "total", "weighted"]}, #diagnosis coverage
             {"tx_prop": ["treat_cov", "total", "weighted"]}, #treatment coverage (among eligible)
             {"mav_n": ["mav_births",{"Under 5":["0-4M", "0-4F"]}, "sum" ]}, #number of births recieving mAVs and HBIG,
             {"prg_scr": ["preg_scr_num", "15-49F", "sum"]}, #number of pregnant women screened for HBsAg,
             {"prg_hrs":["preg_scr_num", "15-49F", "sum"]}] #number of pregnant women screened for HBeAg

    for i in out_res:
        for key, val in i.items():
            df=at.PlotData(res, outputs=val[0], pops=val[1], pop_aggregation=val[2], t_bins=1).series[0].vals
            central_est[key]=df


    # Loop to get all outputs via ensemble method
    for i in out_res:
        for key,val in i.items():
            mapping_function = lambda x: at.PlotData(x,outputs=val[0], pops=val[1], pop_aggregation=val[2], t_bins=1)
            ensemble = at.Ensemble(mapping_function=mapping_function)
            ensemble.update(psa_res)
            df=pd.DataFrame([d.series[0].vals for d in ensemble.samples])
            store_runs[key]=np.array(df).T

    return store_runs, central_est

#%% Economic Outcome Calculations -> convert this to analyse a single set of results
def econ_analysis (cent, cent_s1, cent_s2, cent_s3, res, res_s1, res_s2, res_s3, reg, wd, cost_data, runs, cdisc_rate, hdisc_rate):

    #Discounting Array
    if cdisc_rate >1:
        cdisc_rate=cdisc_rate/100
    else:
        cdisc_rate=cdisc_rate

    if hdisc_rate >1:
        hdisc_rate=hdisc_rate/100
    else:
        hdisc_rate=hdisc_rate

    discount=np.zeros((len(np.arange(1990,2099,1)),4))
    discount[:,0]=np.arange(1990,2099,1)

    for idx,val in enumerate(discount[:,0]):
        if val <= 2021:
            discount[idx,1]=1 #consumption discount
            discount[idx,2]=0
            discount[idx,3]=1 #health discount
        else:
            discount[idx,1:3]=(1-cdisc_rate)**(val-2022)
            discount[idx,3]=(1-hdisc_rate)**(val-2022)



    # Vaccination Costing (seems to be working fine, needs more than 30 runs, and extraction of total births needs checking)
    vax_costs=pd.read_excel("cost and agg data/costs.xlsx", sheet_name="vax")

    np.random.seed(25012023)

    bd_vax=vax_costs[reg].iloc[0]
    hb3_vax=vax_costs[reg].iloc[1]

    bd_vax_samp=np.random.triangular(vax_costs[reg].iloc[6],vax_costs[reg].iloc[0],vax_costs[reg].iloc[7], runs)
    hb3_vax_samp=np.random.triangular(vax_costs[reg].iloc[8],vax_costs[reg].iloc[1],vax_costs[reg].iloc[9 ], runs)

    cbd_cost_bl, chb3_cost_bl, cbd_cost_s1, chb3_cost_s1,cbd_cost_s2, chb3_cost_s2,cbd_cost_s3, chb3_cost_s3=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))
    bd_cost_bl,hb3_cost_bl,bd_cost_s1,hb3_cost_s1,bd_cost_s2,hb3_cost_s2,bd_cost_s3,hb3_cost_s3=np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109, runs)),np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109,runs))

    cbd_cost_bl[:,0]=cent["bd_cov"][:]*cent["births"][:]*bd_vax*discount[:,1]
    chb3_cost_bl[:,0]=cent["hb3_cov"][:]*cent["births"][:]*hb3_vax*discount[:,1]
    cbd_cost_s1[:,0]=cent_s1["bd_cov"][:]*cent_s1["births"][:]*bd_vax*discount[:,1]
    chb3_cost_s1[:,0]=cent_s1["hb3_cov"][:]*cent_s1["births"][:]*hb3_vax*discount[:,1]
    cbd_cost_s2[:,0]=cent_s2["bd_cov"][:]*cent_s2["births"][:]*bd_vax*discount[:,1]
    chb3_cost_s2[:,0]=cent_s2["hb3_cov"][:]*cent_s2["births"][:]*hb3_vax*discount[:,1]
    cbd_cost_s3[:,0]=cent_s3["bd_cov"][:]*cent_s3["births"][:]*bd_vax*discount[:,1]
    chb3_cost_s3[:,0]=cent_s3["hb3_cov"][:]*cent_s3["births"][:]*hb3_vax*discount[:,1]


    for run in range(runs):
        bd_cost_bl[:,run]=res["bd_cov"][:,run]*res["births"][:,run]*bd_vax_samp[run]*discount[:,1]
        hb3_cost_bl[:,run]=res["hb3_cov"][:,run]*res["births"][:,run]*hb3_vax_samp[run]*discount[:,1]
        bd_cost_s1[:,run]=res_s1["bd_cov"][:,run]*res_s1["births"][:,run]*bd_vax_samp[run]*discount[:,1]
        hb3_cost_s1[:,run]=res_s1["hb3_cov"][:,run]*res_s1["births"][:,run]*hb3_vax_samp[run]*discount[:,1]
        bd_cost_s2[:,run]=res_s2["bd_cov"][:,run]*res_s2["births"][:,run]*bd_vax_samp[run]*discount[:,1]
        hb3_cost_s2[:,run]=res_s2["hb3_cov"][:,run]*res_s2["births"][:,run]*hb3_vax_samp[run]*discount[:,1]
        bd_cost_s3[:,run]=res_s3["bd_cov"][:,run]*res_s3["births"][:,run]*bd_vax_samp[run]*discount[:,1]
        hb3_cost_s3[:,run]=res_s3["hb3_cov"][:,run]*res_s3["births"][:,run]*hb3_vax_samp[run]*discount[:,1]

    care_costs=pd.read_excel(cost_data, sheet_name="care")
    dx_cost=care_costs[reg].iloc[0]
    hrs_cost=vax_costs[reg].iloc[10]
    mav_cost=vax_costs[reg].iloc[11]

    #msc- mAV and screening costs (plus HBIG)
    cmsc_cost_bl, cmsc_cost_s1,cmsc_cost_s2,cmsc_cost_s3=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))
    msc_cost_bl, msc_cost_s1,msc_cost_s2,msc_cost_s3=np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109, runs)),np.zeros((109,runs))

    cmsc_cost_bl[:,0]=(cent["prg_scr"][:]*dx_cost*discount[:,1])#(cent["mav_n"][:]*mav_cost*discount[:,1])+(cent["prg_hrs"][:]*hrs_cost*discount[:,1])
    cmsc_cost_s1[:,0]=(cent_s1["mav_n"][:]*mav_cost*discount[:,1])+(cent_s1["prg_scr"][:]*dx_cost*discount[:,1])+(cent_s1["prg_hrs"][:]*hrs_cost*discount[:,2])
    cmsc_cost_s2[:,0]=(cent_s2["mav_n"][:]*mav_cost*discount[:,1])+(cent_s2["prg_scr"][:]*dx_cost*discount[:,1])+(cent_s2["prg_hrs"][:]*hrs_cost*discount[:,2])
    cmsc_cost_s3[:,0]=(cent_s3["mav_n"][:]*mav_cost*discount[:,1])+(cent_s3["prg_scr"][:]*dx_cost*discount[:,1])+(cent_s3["prg_hrs"][:]*hrs_cost*discount[:,2])

    for run in range(runs):
        msc_cost_bl[:,run]=(res["prg_scr"][:,run]*dx_cost*discount[:,1])#(res["mav_n"][:,run]*mav_cost*discount[:,1])+(res["prg_hrs"][:,run]*hrs_cost*discount[:,1])
        msc_cost_s1[:,run]=(res_s1["mav_n"][:,run]*mav_cost*discount[:,1])+(res_s1["prg_scr"][:,run]*dx_cost*discount[:,1])+(res_s1["prg_hrs"][:,run]*hrs_cost*discount[:,2])
        msc_cost_s2[:,run]=(res_s2["mav_n"][:,run]*mav_cost*discount[:,1])+(res_s2["prg_scr"][:,run]*dx_cost*discount[:,1])+(res_s2["prg_hrs"][:,run]*hrs_cost*discount[:,2])
        msc_cost_s3[:,run]=(res_s3["mav_n"][:,run]*mav_cost*discount[:,1])+(res_s3["prg_scr"][:,run]*dx_cost*discount[:,1])+(res_s3["prg_hrs"][:,run]*hrs_cost*discount[:,2])

    # Diagnosis Costs (note three approaches coded, final one [dx_costb] used for analysis [proportional negative testing costs])
    cdx_cost_bl, cdx_cost_s1, cdx_cost_s2, cdx_cost_s3 = np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))
    dx_cost_bl, dx_cost_s1, dx_cost_s2, dx_cost_s3 = np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109,runs))

    cdx_cost_bl[:,0]=cent["dx_rate"][:]*dx_cost*discount[:,1]
    cdx_cost_s1[:,0]=cent_s1["dx_rate"][:]*dx_cost*discount[:,1]
    cdx_cost_s2[:,0]=cent_s2["dx_rate"][:]*dx_cost*discount[:,1]
    cdx_cost_s3[:,0]=cent_s3["dx_rate"][:]*dx_cost*discount[:,1]

    for run in range(runs):
        dx_cost_bl[:,run]=res["dx_rate"][:,run]*dx_cost*discount[:,1]
        dx_cost_s1[:,run]=res_s1["dx_rate"][:, run]*dx_cost*discount[:,1]
        dx_cost_s2[:,run]=res_s2["dx_rate"][:, run]*dx_cost*discount[:,1]
        dx_cost_s3[:,run]=res_s3["dx_rate"][:, run]*dx_cost*discount[:,1]


    ## Dx Cost_alt1 --> 1:5 ratio pos to neg
    cdx_costa_bl, cdx_costa_s1, cdx_costa_s2, cdx_costa_s3 = np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))
    dx_costa_bl, dx_costa_s1, dx_costa_s2, dx_costa_s3 = np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109,runs))

    cdx_costa_bl[:,0]=cent["dx_rate"][:]*(dx_cost*5)*discount[:,1]
    cdx_costa_s1[:,0]=cent_s1["dx_rate"][:]*(dx_cost*5)*discount[:,1]
    cdx_costa_s2[:,0]=cent_s2["dx_rate"][:]*(dx_cost*5)*discount[:,1]
    cdx_costa_s3[:,0]=cent_s3["dx_rate"][:]*(dx_cost*5)*discount[:,1]

    for run in range(runs):
        dx_costa_bl[:,run]=res["dx_rate"][:,run]*(dx_cost*5)*discount[:,1]
        dx_costa_s1[:,run]=res_s1["dx_rate"][:, run]*(dx_cost*5)*discount[:,1]
        dx_costa_s2[:,run]=res_s2["dx_rate"][:, run]*(dx_cost*5)*discount[:,1]
        dx_costa_s3[:,run]=res_s3["dx_rate"][:, run]*(dx_cost*5)*discount[:,1]

    ## Dx Cost_alt2 --> Proportional negative tests
    cbl_dx_inc, cs1_dx_inc,cs2_dx_inc,cs3_dx_inc=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))
    bl_dx_inc, s1_dx_inc,s2_dx_inc,s3_dx_inc=np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109,runs))

    for i in range(len(res["dx_prop"])):
        if i < 1:
            cbl_dx_inc[i,0]=cent["dx_prop"][i]
            cs1_dx_inc[i,0]=cent_s1["dx_prop"][i]
            cs2_dx_inc[i,0]=cent_s2["dx_prop"][i]
            cs3_dx_inc[i,0]=cent_s3["dx_prop"][i]
        else:
            cbl_dx_inc[i,0]=max(cent["dx_prop"][i]-cent["dx_prop"][i-1],0)
            cs1_dx_inc[i,0]=max(cent_s1["dx_prop"][i]-cent_s1["dx_prop"][i-1],0)
            cs2_dx_inc[i,0]=max(cent_s2["dx_prop"][i]-cent_s2["dx_prop"][i-1],0)
            cs3_dx_inc[i,0]=max(cent_s3["dx_prop"][i]-cent_s3["dx_prop"][i-1],0)

    for i in range(len(res["dx_prop"])):
        for run in range(runs):
            if i < 1:
                bl_dx_inc[i,run]=res["dx_prop"][i,run]
                s1_dx_inc[i,run]=res_s1["dx_prop"][i,run]
                s2_dx_inc[i,run]=res_s2["dx_prop"][i,run]
                s3_dx_inc[i,run]=res_s3["dx_prop"][i,run]
            else:
                bl_dx_inc[i,run]=max(res["dx_prop"][i,run]-res["dx_prop"][i-1,run],0)
                s1_dx_inc[i,run]=max(res_s1["dx_prop"][i,run]-res_s1["dx_prop"][i-1,run],0)
                s2_dx_inc[i,run]=max(res_s2["dx_prop"][i,run]-res_s2["dx_prop"][i-1,run],0)
                s3_dx_inc[i,run]=max(res_s3["dx_prop"][i,run]-res_s3["dx_prop"][i-1,run],0)

    cdx_costb_bl, cdx_costb_s1, cdx_costb_s2, cdx_costb_s3 = np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))
    dx_costb_bl, dx_costb_s1, dx_costb_s2, dx_costb_s3 = np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109,runs))

    for yr in range(len(cdx_costb_bl)):
        cdx_costb_bl[yr,0]=(cent["dx_rate"][yr]*dx_cost*discount[yr,1])+(dx_cost*cbl_dx_inc[yr,0]*(cent["pop"][yr]-cent["dx_rate"][yr])*discount[yr,1])
        cdx_costb_s1[yr,0]=(cent_s1["dx_rate"][yr]*dx_cost*discount[yr,1])+(dx_cost*cs1_dx_inc[yr,0]*(cent_s1["pop"][yr]-cent_s1["dx_rate"][yr])*discount[yr,1])
        cdx_costb_s2[yr,0]=(cent_s2["dx_rate"][yr]*dx_cost*discount[yr,1])+(dx_cost*cs2_dx_inc[yr,0]*(cent_s2["pop"][yr]-cent_s2["dx_rate"][yr])*discount[yr,1])
        cdx_costb_s3[yr,0]=(cent_s3["dx_rate"][yr]*dx_cost*discount[yr,1])+(dx_cost*cs3_dx_inc[yr,0]*(cent_s3["pop"][yr]-cent_s3["dx_rate"][yr])*discount[yr,1])

    for run in range(runs):
        for yr in range(len(dx_costb_bl)):
            dx_costb_bl[yr,run]=(res["dx_rate"][yr,run]*dx_cost*discount[yr,1])+(dx_cost*bl_dx_inc[yr,run]*(res["pop"][yr,run]-res["dx_rate"][yr,run])*discount[yr,1])
            dx_costb_s1[yr,run]=(res_s1["dx_rate"][yr, run]*dx_cost*discount[yr,1])+(dx_cost*s1_dx_inc[yr,run]*(res_s1["pop"][yr,run]-res_s1["dx_rate"][yr,run])*discount[yr,1])
            dx_costb_s2[yr,run]=(res_s2["dx_rate"][yr, run]*dx_cost*discount[yr,1])+(dx_cost*s2_dx_inc[yr,run]*(res_s2["pop"][yr,run]-res_s2["dx_rate"][yr,run])*discount[yr,1])
            dx_costb_s3[yr,run]=(res_s3["dx_rate"][yr, run]*dx_cost*discount[yr,1])+(dx_cost*s3_dx_inc[yr,run]*(res_s3["pop"][yr,run]-res_s3["dx_rate"][yr,run])*discount[yr,1])

    # Treatment Costs
    tx_cost=care_costs[reg].iloc[3]

    ctx_cost_bl, ctx_cost_s1, ctx_cost_s2, ctx_cost_s3 = np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))
    tx_cost_bl, tx_cost_s1, tx_cost_s2, tx_cost_s3 = np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109,runs)),np.zeros((109,runs))

    ctx_cost_bl[:,0]=cent["tx_cov"][:]*tx_cost*discount[:,1]
    ctx_cost_s1[:,0]=cent_s1["tx_cov"][:]*tx_cost*discount[:,1]
    ctx_cost_s2[:,0]=cent_s2["tx_cov"][:]*tx_cost*discount[:,1]
    ctx_cost_s3[:,0]=cent_s3["tx_cov"][:]*tx_cost*discount[:,1]

    for run in range(runs):
        tx_cost_bl[:,run]=res["tx_cov"][:,run]*tx_cost*discount[:,1]
        tx_cost_s1[:, run]=res_s1["tx_cov"][:, run]*tx_cost*discount[:,1]
        tx_cost_s2[:, run]=res_s2["tx_cov"][:, run]*tx_cost*discount[:,1]
        tx_cost_s3[:, run]=res_s3["tx_cov"][:, run]*tx_cost*discount[:,1]

    ## Direct Medical Costs (already discounted; central)
    dmc_cent_bl=cbd_cost_bl+chb3_cost_bl+cdx_costb_bl+ctx_cost_bl+cmsc_cost_bl
    dmc_cent_s1=cbd_cost_s1+chb3_cost_s1+cdx_costb_s1+ctx_cost_s1+cmsc_cost_s1
    dmc_cent_s2=cbd_cost_s2+chb3_cost_s2+cdx_costb_s2+ctx_cost_s2+cmsc_cost_s2
    dmc_cent_s3=cbd_cost_s3+chb3_cost_s3+cdx_costb_s3+ctx_cost_s3+cmsc_cost_s3

    ## Direct Medical Costs (already discounted; sampled)
    dmc_psa_bl=bd_cost_bl+hb3_cost_bl+dx_costb_bl+tx_cost_bl+msc_cost_bl
    dmc_psa_s1=bd_cost_s1+hb3_cost_s1+dx_costb_s1+tx_cost_s1+msc_cost_s1
    dmc_psa_s2=bd_cost_s2+hb3_cost_s2+dx_costb_s2+tx_cost_s2+msc_cost_s2
    dmc_psa_s3=bd_cost_s3+hb3_cost_s3+dx_costb_s3+tx_cost_s3+msc_cost_s3

    ## Disease Progression Management costs
    util=0.25  #care utilization assumption
    tx_hosp=0.5  #treatment preventing hospitalisation assumption

    ## Indirect Disease Management Costs (assume 50% HR covered under UHC; 25% hospital utilization)
    dx_dmc=care_costs[reg].iloc[4] #cost of disease management (diagnosed)
    tx_dmc=care_costs[reg].iloc[4] #cost of disease management (on treatment)
    hosp_cost=care_costs[reg].iloc[7] #cost of hospitalisation
    hcc_cost=care_costs[reg].iloc[6] #HCC surveillance costs
    hcc_prp=care_costs[reg].iloc[8] #proportion of CHB in 15-49y who are 40-49y (HCC screening)

    cmc_cost_bl, cmc_cost_s1, cmc_cost_s2, cmc_cost_s3= np.zeros((109,1)),np.zeros((109, 1)),np.zeros((109, 1)),np.zeros((109, 1))
    chosp_cost_bl, chosp_cost_s1, chosp_cost_s2, chosp_cost_s3=np.zeros((109, 1)),np.zeros((109,1)),np.zeros((109, 1)),np.zeros((109, 1))
    chcc_cost_bl, chcc_cost_s1, chcc_cost_s2, chcc_cost_s3=np.zeros((109, 1)),np.zeros((109, 1)),np.zeros((109, 1)),np.zeros((109, 1))

    ## Disease Management Costs
    cmc_cost_bl[:, 0]=((cent["dm_dx"][:]*dx_dmc*util)+(cent["dm_tx"][:]*tx_dmc*util))*discount[:,1]
    cmc_cost_s1[:, 0]=((cent_s1["dm_dx"][:]*dx_dmc*util)+(cent_s1["dm_tx"][:]*tx_dmc*util))*discount[:,1]
    cmc_cost_s2[:, 0]=((cent_s2["dm_dx"][:]*dx_dmc*util)+(cent_s2["dm_tx"][:]*tx_dmc*util))*discount[:,1]
    cmc_cost_s3[:, 0]=((cent_s3["dm_dx"][:]*dx_dmc*util)+(cent_s3["dm_tx"][:]*tx_dmc*util))*discount[:,1]
    ## Hospitalisation Costs
    chosp_cost_bl[:, 0]=((cent["hsp_utx"][:]*hosp_cost*util)+(cent["hsp_tx"][:]*hosp_cost*util*tx_hosp))*discount[:,1]
    chosp_cost_s1[:, 0]=((cent_s1["hsp_utx"][:]*hosp_cost*util)+(cent_s1["hsp_tx"][:]*hosp_cost*util*tx_hosp))*discount[:,1]
    chosp_cost_s2[:, 0]=((cent_s2["hsp_utx"][:]*hosp_cost*util)+(cent_s2["hsp_tx"][:]*hosp_cost*util*tx_hosp))*discount[:,1]
    chosp_cost_s3[:, 0]=((cent_s3["hsp_utx"][:]*hosp_cost*util)+(cent_s3["hsp_tx"][:]*hosp_cost*util*tx_hosp))*discount[:,1]
    ## HCC Surveillance Costs (assume 50% HR covered under UHC; 25% resource utilization)
    chcc_cost_bl[:, 0]=((cent["pop_hcc"][:]*hcc_cost*util)+(cent["tgt_hcc"][:]*hcc_cost*util)+(cent["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util))*discount[:,1]
    chcc_cost_s1[:, 0]=((cent_s1["pop_hcc"][:]*hcc_cost*util)+(cent_s1["tgt_hcc"][:]*hcc_cost*util)+(cent_s1["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util))*discount[:,1]
    chcc_cost_s2[:, 0]=((cent_s2["pop_hcc"][:]*hcc_cost*util)+(cent_s2["tgt_hcc"][:]*hcc_cost*util)+(cent_s2["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util))*discount[:,1]
    chcc_cost_s3[:, 0]=((cent_s3["pop_hcc"][:]*hcc_cost*util)+(cent_s3["tgt_hcc"][:]*hcc_cost*util)+(cent_s3["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util))*discount[:,1]


    mc_cost_bl, mc_cost_s1, mc_cost_s2, mc_cost_s3= np.zeros((109,runs)),np.zeros((109, runs)),np.zeros((109, runs)),np.zeros((109, runs))
    hosp_cost_bl, hosp_cost_s1, hosp_cost_s2, hosp_cost_s3=np.zeros((109, runs)),np.zeros((109, runs)),np.zeros((109, runs)),np.zeros((109, runs))
    hcc_cost_bl, hcc_cost_s1, hcc_cost_s2, hcc_cost_s3=np.zeros((109, runs)),np.zeros((109, runs)),np.zeros((109, runs)),np.zeros((109, runs))

    for run in range(runs):
        ## Disease Management Costs
        mc_cost_bl[:, run]=((res["dm_dx"][:, run]*dx_dmc*util)+(res["dm_tx"][:, run]*tx_dmc*util))*discount[:,1]
        mc_cost_s1[:, run]=((res_s1["dm_dx"][:, run]*dx_dmc*util)+(res_s1["dm_tx"][:, run]*tx_dmc*util))*discount[:,1]
        mc_cost_s2[:, run]=((res_s2["dm_dx"][:, run]*dx_dmc*util)+(res_s2["dm_tx"][:, run]*tx_dmc*util))*discount[:,1]
        mc_cost_s3[:, run]=((res_s3["dm_dx"][:, run]*dx_dmc*util)+(res_s3["dm_tx"][:, run]*tx_dmc*util))*discount[:,1]
        ## Hospitalisation Costs
        hosp_cost_bl[:, run]=((res["hsp_utx"][:,run]*hosp_cost*util)+(res["hsp_tx"][:, run]*hosp_cost*util*tx_hosp))*discount[:,1]
        hosp_cost_s1[:, run]=((res_s1["hsp_utx"][:,run]*hosp_cost*util)+(res_s1["hsp_tx"][:, run]*hosp_cost*util*tx_hosp))*discount[:,1]
        hosp_cost_s2[:, run]=((res_s2["hsp_utx"][:,run]*hosp_cost*util)+(res_s2["hsp_tx"][:, run]*hosp_cost*util*tx_hosp))*discount[:,1]
        hosp_cost_s3[:, run]=((res_s3["hsp_utx"][:,run]*hosp_cost*util)+(res_s3["hsp_tx"][:, run]*hosp_cost*util*tx_hosp))*discount[:,1]
        ## HCC Surveillance Costs (assume 50% HR covered under UHC; 25% resource utilization)
        hcc_cost_bl[:, run]=((res["pop_hcc"][:, run]*hcc_cost*util)+(res["tgt_hcc"][:, run]*hcc_cost*util)+(res["tgt_hcc_b"][:,run]*hcc_prp*hcc_cost*util))*discount[:,1]
        hcc_cost_s1[:, run]=((res_s1["pop_hcc"][:, run]*hcc_cost*util)+(res_s1["tgt_hcc"][:, run]*hcc_cost*util)+(res_s1["tgt_hcc_b"][:,run]*hcc_prp*hcc_cost*util))*discount[:,1]
        hcc_cost_s2[:, run]=((res_s2["pop_hcc"][:, run]*hcc_cost*util)+(res_s2["tgt_hcc"][:, run]*hcc_cost*util)+(res_s2["tgt_hcc_b"][:,run]*hcc_prp*hcc_cost*util))*discount[:,1]
        hcc_cost_s3[:, run]=((res_s3["pop_hcc"][:, run]*hcc_cost*util)+(res_s3["tgt_hcc"][:, run]*hcc_cost*util)+(res_s3["tgt_hcc_b"][:,run]*hcc_prp*hcc_cost*util))*discount[:,1]

    ## Indirect Medical Costs (discounted; central)
    imc_cent_bl=cmc_cost_bl+chosp_cost_bl+chcc_cost_bl
    imc_cent_s1=cmc_cost_s1+chosp_cost_s1+chcc_cost_s1
    imc_cent_s2=cmc_cost_s2+chosp_cost_s2+chcc_cost_s2
    imc_cent_s3=cmc_cost_s3+chosp_cost_s3+chcc_cost_s3

    ## Indirect Medical Costs (discounted; sampled)
    imc_psa_bl=mc_cost_bl+hosp_cost_bl+hcc_cost_bl
    imc_psa_s1=mc_cost_s1+hosp_cost_s1+hcc_cost_s1
    imc_psa_s2=mc_cost_s2+hosp_cost_s2+hcc_cost_s2
    imc_psa_s3=mc_cost_s3+hosp_cost_s3+hcc_cost_s3


    ## Productivity Loss (calculate for baseline, check outcomes, then repeat)
    prod_costs=pd.read_excel(cost_data, sheet_name="emp_gdp_lex")
    etp_ratio=prod_costs[reg].iloc[0]
    gdp=prod_costs[reg].iloc[1]
    life_exp=prod_costs[reg].iloc[2]

    ## GDP growth
    gdp_grw=np.zeros((len(np.arange(1990,2099,1)),4))
    gdp_grw[:,0]=np.arange(1990,2099,1)
    gdp_grw[:,1:4]=gdp

    for i,val in enumerate(gdp_grw[:,0]):
        if val>2022:
            gdp_grw[i,1]=gdp_grw[i-1,1]*1.00
            gdp_grw[i,2]=gdp_grw[i-1,2]*1.015
            gdp_grw[i,3]=gdp_grw[i-1,3]*1.03

    age_of_deaths=np.array([0.01, 0.031, 0.253, 0.341, 0.365])
    prop_leaving_age_categories=np.array([1/15, 1/15, 1/20, 1/15])
    all_cause_mort=np.array([0.003, 0.0013, 0.0022, 0.0103, (1/life_exp)])

    ## Baseline (central)
    cbl_deaths=np.zeros((len(cent["mort"]),2))
    cbl_deaths[:,0]=np.arange(1990,2099,1)
    cbl_deaths[:,1]=cent["mort"]

    for idx,val in enumerate(cbl_deaths[:,0]):
        if val < 2022:
            cbl_deaths[idx,1]=0

    ghosts_cbl=np.zeros((len(cbl_deaths), len(age_of_deaths)))
    ghosts_cbl[0,:]=cbl_deaths[0,1]*age_of_deaths

    for t in range(1,len(cbl_deaths)):
        ppl_who_age=ghosts_cbl[t,0:len(prop_leaving_age_categories)]*prop_leaving_age_categories
        ghosts_cbl[t,0]=max(0, ghosts_cbl[t-1,0]-ppl_who_age[0]-all_cause_mort[0]*ghosts_cbl[t-1,0])
        ghosts_cbl[t,1]=max(0, ghosts_cbl[t-1,1]-ppl_who_age[1]+ppl_who_age[0]-all_cause_mort[1]*ghosts_cbl[t-1,1])
        ghosts_cbl[t,2]=max(0, ghosts_cbl[t-1,2]-ppl_who_age[2]+ppl_who_age[1]-all_cause_mort[2]*ghosts_cbl[t-1,2])
        ghosts_cbl[t,3]=max(0, ghosts_cbl[t-1,3]-ppl_who_age[3]+ppl_who_age[2]-all_cause_mort[3]*ghosts_cbl[t-1,3])
        ghosts_cbl[t,4]=max(0, ghosts_cbl[t-1,4]+ppl_who_age[3]-all_cause_mort[4]*ghosts_cbl[t-1,4])

        ghosts_cbl[t,:]= ghosts_cbl[t,:]+cbl_deaths[t,1]*age_of_deaths

    ## Baseline (sampled)
    bl_deaths=np.zeros((len(res["mort"]), runs+1))
    bl_deaths[:,0]=np.arange(1990,2099,1)
    bl_deaths[:,1:]=res["mort"]

    for idx,val in enumerate(bl_deaths[:,0]):
        if val < 2022:
            bl_deaths[idx,1:]=0

    ghosts_bl=np.zeros((runs,len(bl_deaths),len(age_of_deaths)))

    for run in range(runs):
        ghosts_bl[run,0,:]=bl_deaths[0,run+1]*age_of_deaths

    for run in range(runs):
        for t in range(1,len(bl_deaths)):
            ppl_who_age=ghosts_bl[run,t,0:len(prop_leaving_age_categories)]*prop_leaving_age_categories
            ghosts_bl[run,t,0]=max(0, ghosts_bl[run,t-1,0]-ppl_who_age[0]-all_cause_mort[0]*ghosts_bl[run,t-1,0])
            ghosts_bl[run,t,1]=max(0, ghosts_bl[run,t-1,1]-ppl_who_age[1]+ppl_who_age[0]-all_cause_mort[1]*ghosts_bl[run,t-1,1])
            ghosts_bl[run,t,2]=max(0, ghosts_bl[run,t-1,2]-ppl_who_age[2]+ppl_who_age[1]-all_cause_mort[2]*ghosts_bl[run,t-1,2])
            ghosts_bl[run,t,3]=max(0, ghosts_bl[run,t-1,3]-ppl_who_age[3]+ppl_who_age[2]-all_cause_mort[3]*ghosts_bl[run,t-1,3])
            ghosts_bl[run,t,4]=max(0, ghosts_bl[run,t-1,4]+ppl_who_age[3]-all_cause_mort[4]*ghosts_bl[run,t-1,4])

            ghosts_bl[run,t,:]= ghosts_bl[run,t,:]+bl_deaths[t,run+1]*age_of_deaths

    ## Scenario 1 (central)
    cs1_deaths=np.zeros((len(cent_s1["mort"]),2))
    cs1_deaths[:,0]=np.arange(1990,2099,1)
    cs1_deaths[:,1]=cent_s1["mort"]

    for idx,val in enumerate(cs1_deaths[:,0]):
        if val < 2022:
            cs1_deaths[idx,1]=0

    ghosts_cs1=np.zeros((len(cs1_deaths), len(age_of_deaths)))
    ghosts_cs1[0,:]=cs1_deaths[0,1]*age_of_deaths

    for t in range(1,len(cs1_deaths)):
        ppl_who_age=ghosts_cs1[t,0:len(prop_leaving_age_categories)]*prop_leaving_age_categories
        ghosts_cs1[t,0]=max(0, ghosts_cs1[t-1,0]-ppl_who_age[0]-all_cause_mort[0]*ghosts_cs1[t-1,0])
        ghosts_cs1[t,1]=max(0, ghosts_cs1[t-1,1]-ppl_who_age[1]+ppl_who_age[0]-all_cause_mort[1]*ghosts_cs1[t-1,1])
        ghosts_cs1[t,2]=max(0, ghosts_cs1[t-1,2]-ppl_who_age[2]+ppl_who_age[1]-all_cause_mort[2]*ghosts_cs1[t-1,2])
        ghosts_cs1[t,3]=max(0, ghosts_cs1[t-1,3]-ppl_who_age[3]+ppl_who_age[2]-all_cause_mort[3]*ghosts_cs1[t-1,3])
        ghosts_cs1[t,4]=max(0, ghosts_cs1[t-1,4]+ppl_who_age[3]-all_cause_mort[4]*ghosts_cs1[t-1,4])

        ghosts_cs1[t,:]= ghosts_cs1[t,:]+cs1_deaths[t,1]*age_of_deaths


    ## Scenario 1 (sampled)
    s1_deaths=np.zeros((len(res_s1["mort"]), runs+1))
    s1_deaths[:,0]=np.arange(1990,2099,1)
    s1_deaths[:,1:]=res_s1["mort"]

    for idx,val in enumerate(s1_deaths[:,0]):
        if val < 2022:
            s1_deaths[idx,1:]=0

    ghosts_s1=np.zeros((runs,len(s1_deaths),len(age_of_deaths)))

    for run in range(runs):
        ghosts_s1[run,0,:]=s1_deaths[0,run+1]*age_of_deaths

    for run in range(runs):
        for t in range(1,len(s1_deaths)):
            ppl_who_age=ghosts_s1[run,t,0:len(prop_leaving_age_categories)]*prop_leaving_age_categories
            ghosts_s1[run,t,0]=max(0, ghosts_s1[run,t-1,0]-ppl_who_age[0]-all_cause_mort[0]*ghosts_s1[run,t-1,0])
            ghosts_s1[run,t,1]=max(0, ghosts_s1[run,t-1,1]-ppl_who_age[1]+ppl_who_age[0]-all_cause_mort[1]*ghosts_s1[run,t-1,1])
            ghosts_s1[run,t,2]=max(0, ghosts_s1[run,t-1,2]-ppl_who_age[2]+ppl_who_age[1]-all_cause_mort[2]*ghosts_s1[run,t-1,2])
            ghosts_s1[run,t,3]=max(0, ghosts_s1[run,t-1,3]-ppl_who_age[3]+ppl_who_age[2]-all_cause_mort[3]*ghosts_s1[run,t-1,3])
            ghosts_s1[run,t,4]=max(0, ghosts_s1[run,t-1,4]+ppl_who_age[3]-all_cause_mort[4]*ghosts_s1[run,t-1,4])

            ghosts_s1[run,t,:]= ghosts_s1[run,t,:]+s1_deaths[t,run+1]*age_of_deaths

    ## Scenario 2 (central)
    cs2_deaths=np.zeros((len(cent_s2["mort"]),2))
    cs2_deaths[:,0]=np.arange(1990,2099,1)
    cs2_deaths[:,1]=cent_s2["mort"]

    for idx,val in enumerate(cs2_deaths[:,0]):
        if val < 2022:
            cs2_deaths[idx,1]=0

    ghosts_cs2=np.zeros((len(cs2_deaths), len(age_of_deaths)))
    ghosts_cs2[0,:]=cs2_deaths[0,1]*age_of_deaths

    for t in range(1,len(cs2_deaths)):
        ppl_who_age=ghosts_cs2[t,0:len(prop_leaving_age_categories)]*prop_leaving_age_categories
        ghosts_cs2[t,0]=max(0, ghosts_cs2[t-1,0]-ppl_who_age[0]-all_cause_mort[0]*ghosts_cs2[t-1,0])
        ghosts_cs2[t,1]=max(0, ghosts_cs2[t-1,1]-ppl_who_age[1]+ppl_who_age[0]-all_cause_mort[1]*ghosts_cs2[t-1,1])
        ghosts_cs2[t,2]=max(0, ghosts_cs2[t-1,2]-ppl_who_age[2]+ppl_who_age[1]-all_cause_mort[2]*ghosts_cs2[t-1,2])
        ghosts_cs2[t,3]=max(0, ghosts_cs2[t-1,3]-ppl_who_age[3]+ppl_who_age[2]-all_cause_mort[3]*ghosts_cs2[t-1,3])
        ghosts_cs2[t,4]=max(0, ghosts_cs2[t-1,4]+ppl_who_age[3]-all_cause_mort[4]*ghosts_cs2[t-1,4])

        ghosts_cs2[t,:]= ghosts_cs2[t,:]+cs2_deaths[t,1]*age_of_deaths

    ## Scenario 2 (sampled)
    s2_deaths=np.zeros((len(res_s2["mort"]), runs+1))
    s2_deaths[:,0]=np.arange(1990,2099,1)
    s2_deaths[:,1:]=res_s2["mort"]

    for idx,val in enumerate(s2_deaths[:,0]):
        if val < 2022:
            s2_deaths[idx,1:]=0

    ghosts_s2=np.zeros((runs,len(s2_deaths),len(age_of_deaths)))

    for run in range(runs):
        ghosts_s2[run,0,:]=s2_deaths[0,run+1]*age_of_deaths

    for run in range(runs):
        for t in range(1,len(s2_deaths)):
            ppl_who_age=ghosts_s2[run,t,0:len(prop_leaving_age_categories)]*prop_leaving_age_categories
            ghosts_s2[run,t,0]=max(0, ghosts_s2[run,t-1,0]-ppl_who_age[0]-all_cause_mort[0]*ghosts_s2[run,t-1,0])
            ghosts_s2[run,t,1]=max(0, ghosts_s2[run,t-1,1]-ppl_who_age[1]+ppl_who_age[0]-all_cause_mort[1]*ghosts_s2[run,t-1,1])
            ghosts_s2[run,t,2]=max(0, ghosts_s2[run,t-1,2]-ppl_who_age[2]+ppl_who_age[1]-all_cause_mort[2]*ghosts_s2[run,t-1,2])
            ghosts_s2[run,t,3]=max(0, ghosts_s2[run,t-1,3]-ppl_who_age[3]+ppl_who_age[2]-all_cause_mort[3]*ghosts_s2[run,t-1,3])
            ghosts_s2[run,t,4]=max(0, ghosts_s2[run,t-1,4]+ppl_who_age[3]-all_cause_mort[4]*ghosts_s2[run,t-1,4])

            ghosts_s2[run,t,:]= ghosts_s2[run,t,:]+s2_deaths[t,run+1]*age_of_deaths

    ## Scenario 3 (central)
    cs3_deaths=np.zeros((len(cent_s3["mort"]),2))
    cs3_deaths[:,0]=np.arange(1990,2099,1)
    cs3_deaths[:,1]=cent_s3["mort"]

    for idx,val in enumerate(cs3_deaths[:,0]):
        if val < 2022:
            cs3_deaths[idx,1]=0

    ghosts_cs3=np.zeros((len(cs3_deaths), len(age_of_deaths)))
    ghosts_cs3[0,:]=cs3_deaths[0,1]*age_of_deaths

    for t in range(1,len(cs3_deaths)):
        ppl_who_age=ghosts_cs3[t,0:len(prop_leaving_age_categories)]*prop_leaving_age_categories
        ghosts_cs3[t,0]=max(0, ghosts_cs3[t-1,0]-ppl_who_age[0]-all_cause_mort[0]*ghosts_cs3[t-1,0])
        ghosts_cs3[t,1]=max(0, ghosts_cs3[t-1,1]-ppl_who_age[1]+ppl_who_age[0]-all_cause_mort[1]*ghosts_cs3[t-1,1])
        ghosts_cs3[t,2]=max(0, ghosts_cs3[t-1,2]-ppl_who_age[2]+ppl_who_age[1]-all_cause_mort[2]*ghosts_cs3[t-1,2])
        ghosts_cs3[t,3]=max(0, ghosts_cs3[t-1,3]-ppl_who_age[3]+ppl_who_age[2]-all_cause_mort[3]*ghosts_cs3[t-1,3])
        ghosts_cs3[t,4]=max(0, ghosts_cs3[t-1,4]+ppl_who_age[3]-all_cause_mort[4]*ghosts_cs3[t-1,4])

        ghosts_cs3[t,:]= ghosts_cs3[t,:]+cs3_deaths[t,1]*age_of_deaths

    ## Scenario 3 (sampled)
    s3_deaths=np.zeros((len(res_s3["mort"]), runs+1))
    s3_deaths[:,0]=np.arange(1990,2099,1)
    s3_deaths[:,1:]=res_s3["mort"]

    for idx,val in enumerate(s3_deaths[:,0]):
        if val < 2022:
            s3_deaths[idx,1:]=0

    ghosts_s3=np.zeros((runs,len(s3_deaths),len(age_of_deaths)))

    for run in range(runs):
        ghosts_s3[run,0,:]=s3_deaths[0,run+1]*age_of_deaths

    for run in range(runs):
        for t in range(1,len(s3_deaths)):
            ppl_who_age=ghosts_s3[run,t,0:len(prop_leaving_age_categories)]*prop_leaving_age_categories
            ghosts_s3[run,t,0]=max(0, ghosts_s3[run,t-1,0]-ppl_who_age[0]-all_cause_mort[0]*ghosts_s3[run,t-1,0])
            ghosts_s3[run,t,1]=max(0, ghosts_s3[run,t-1,1]-ppl_who_age[1]+ppl_who_age[0]-all_cause_mort[1]*ghosts_s3[run,t-1,1])
            ghosts_s3[run,t,2]=max(0, ghosts_s3[run,t-1,2]-ppl_who_age[2]+ppl_who_age[1]-all_cause_mort[2]*ghosts_s3[run,t-1,2])
            ghosts_s3[run,t,3]=max(0, ghosts_s3[run,t-1,3]-ppl_who_age[3]+ppl_who_age[2]-all_cause_mort[3]*ghosts_s3[run,t-1,3])
            ghosts_s3[run,t,4]=max(0, ghosts_s3[run,t-1,4]+ppl_who_age[3]-all_cause_mort[4]*ghosts_s3[run,t-1,4])

            ghosts_s3[run,t,:]= ghosts_s3[run,t,:]+s3_deaths[t,run+1]*age_of_deaths


    ## DALYs
    cbl_yll=np.sum(ghosts_cbl[:,0:5], axis=1)*discount[:,3]
    cs1_yll=np.sum(ghosts_cs1[:,0:5], axis=1)*discount[:,3]
    cs2_yll=np.sum(ghosts_cs2[:,0:5], axis=1)*discount[:,3]
    cs3_yll=np.sum(ghosts_cs3[:,0:5], axis=1)*discount[:,3]

    bl_yll=(np.sum(ghosts_bl[:,:,0:5], axis=2).T)
    s1_yll=(np.sum(ghosts_s1[:,:,0:5], axis=2).T)
    s2_yll=(np.sum(ghosts_s2[:,:,0:5], axis=2).T)
    s3_yll=(np.sum(ghosts_s3[:,:,0:5], axis=2).T)

    for idx in range(len(bl_yll)):
        bl_yll[idx,:]=bl_yll[idx,:]*discount[idx,3]
        s1_yll[idx,:]=s1_yll[idx,:]*discount[idx,3]
        s2_yll[idx,:]=s2_yll[idx,:]*discount[idx,3]
        s3_yll[idx,:]=s3_yll[idx,:]*discount[idx,3]

    cbl_dalys=cbl_yll+cent["yld"]
    cs1_dalys=cs1_yll+cent_s1["yld"]
    cs2_dalys=cs2_yll+cent_s2["yld"]
    cs3_dalys=cs3_yll+cent_s3["yld"]

    bl_dalys=bl_yll+res["yld"]
    s1_dalys=s1_yll+res_s1["yld"]
    s2_dalys=s2_yll+res_s2["yld"]
    s3_dalys=s3_yll+res_s3["yld"]

    ## Productivity
    cbl_prod=np.sum(ghosts_cbl[:,0:4], axis=1)*gdp_grw[:,2]*etp_ratio*discount[:,3]
    cs1_prod=np.sum(ghosts_cs1[:,0:4], axis=1)*gdp_grw[:,2]*etp_ratio*discount[:,3]
    cs2_prod=np.sum(ghosts_cs2[:,0:4], axis=1)*gdp_grw[:,2]*etp_ratio*discount[:,3]
    cs3_prod=np.sum(ghosts_cs3[:,0:4], axis=1)*gdp_grw[:,2]*etp_ratio*discount[:,3]

    bl_prod=(np.sum(ghosts_bl[:,:,0:4], axis=2).T)
    s1_prod=(np.sum(ghosts_s1[:,:,0:4], axis=2).T)
    s2_prod=(np.sum(ghosts_s2[:,:,0:4], axis=2).T)
    s3_prod=(np.sum(ghosts_s3[:,:,0:4], axis=2).T)

    for idx in range(len(bl_yll)):
        bl_prod[idx,:]=bl_prod[idx,:]*gdp_grw[idx,2]*etp_ratio*discount[idx,3]
        s1_prod[idx,:]=s1_prod[idx,:]*gdp_grw[idx,2]*etp_ratio*discount[idx,3]
        s2_prod[idx,:]=s2_prod[idx,:]*gdp_grw[idx,2]*etp_ratio*discount[idx,3]
        s3_prod[idx,:]=s3_prod[idx,:]*gdp_grw[idx,2]*etp_ratio*discount[idx,3]


    ##(Note: Both Measures Require Cumulative Costs, only instantaneous calculated so far)
    ## Net Economic Benefit (up to here for central estimates)
    cbl_ins_tc=cbl_prod[:]+dmc_cent_bl[:,0]+imc_cent_bl[:,0]
    cs1_ins_tc=cs1_prod[:]+dmc_cent_s1[:,0]+imc_cent_s1[:,0]
    cs2_ins_tc=cs2_prod[:]+dmc_cent_s2[:,0]+imc_cent_s2[:,0]
    cs3_ins_tc=cs3_prod[:]+dmc_cent_s3[:,0]+imc_cent_s3[:,0]

    cbl_cum_tc=np.zeros(np.shape(cbl_ins_tc))
    cs1_cum_tc=np.zeros(np.shape(cs1_ins_tc))
    cs2_cum_tc=np.zeros(np.shape(cs2_ins_tc))
    cs3_cum_tc=np.zeros(np.shape(cs3_ins_tc))

    for i in range(len(cbl_ins_tc)):
        if i < 1:
            cbl_cum_tc[i]= cbl_ins_tc[i]
            cs1_cum_tc[i]= cs1_ins_tc[i]
            cs2_cum_tc[i]= cs2_ins_tc[i]
            cs3_cum_tc[i]= cs3_ins_tc[i]
        else:
            cbl_cum_tc[i]= cbl_cum_tc[i-1]+cbl_ins_tc[i]
            cs1_cum_tc[i]= cs1_cum_tc[i-1]+cs1_ins_tc[i]
            cs2_cum_tc[i]= cs2_cum_tc[i-1]+cs2_ins_tc[i]
            cs3_cum_tc[i]= cs3_cum_tc[i-1]+cs3_ins_tc[i]

    bl_ins_tc=bl_prod+dmc_psa_bl+imc_psa_bl
    s1_ins_tc=s1_prod+dmc_psa_s1+imc_psa_s1
    s2_ins_tc=s2_prod+dmc_psa_s2+imc_psa_s2
    s3_ins_tc=s3_prod+dmc_psa_s3+imc_psa_s3

    bl_cum_tc=np.zeros(np.shape(bl_ins_tc))
    s1_cum_tc=np.zeros(np.shape(s1_ins_tc))
    s2_cum_tc=np.zeros(np.shape(s2_ins_tc))
    s3_cum_tc=np.zeros(np.shape(s3_ins_tc))

    for i in range(len(bl_ins_tc)):
        if i < 1:
            bl_cum_tc[i,:]=bl_ins_tc[i,:]
            s1_cum_tc[i,:]=s1_ins_tc[i,:]
            s2_cum_tc[i,:]=s2_ins_tc[i,:]
            s3_cum_tc[i,:]=s3_ins_tc[i,:]
        else:
            bl_cum_tc[i,:]=bl_cum_tc[i-1,:]+bl_ins_tc[i,:]
            s1_cum_tc[i,:]=s1_cum_tc[i-1,:]+s1_ins_tc[i,:]
            s2_cum_tc[i,:]=s2_cum_tc[i-1,:]+s2_ins_tc[i,:]
            s3_cum_tc[i,:]=s3_cum_tc[i-1,:]+s3_ins_tc[i,:]

    cs1_neb=cbl_cum_tc-cs1_cum_tc
    cs2_neb=cbl_cum_tc-cs2_cum_tc
    cs3_neb=cbl_cum_tc-cs3_cum_tc

    s1_neb=bl_cum_tc-s1_cum_tc
    s2_neb=bl_cum_tc-s2_cum_tc
    s3_neb=bl_cum_tc-s3_cum_tc

    ## ICERs
    cbl_ins_dc=dmc_cent_bl[:,0]+imc_cent_bl[:,0]
    cs1_ins_dc=dmc_cent_s1[:,0]+imc_cent_s1[:,0]
    cs2_ins_dc=dmc_cent_s2[:,0]+imc_cent_s2[:,0]
    cs3_ins_dc=dmc_cent_s3[:,0]+imc_cent_s3[:,0]

    cbl_cum_dc=np.zeros(np.shape(cbl_ins_dc))
    cs1_cum_dc=np.zeros(np.shape(cs1_ins_dc))
    cs2_cum_dc=np.zeros(np.shape(cs2_ins_dc))
    cs3_cum_dc=np.zeros(np.shape(cs3_ins_dc))

    cbl_cum_daly=np.zeros(np.shape(cbl_dalys))
    cs1_cum_daly=np.zeros(np.shape(cs1_dalys))
    cs2_cum_daly=np.zeros(np.shape(cs2_dalys))
    cs3_cum_daly=np.zeros(np.shape(cs3_dalys))

    for i in range(len(cbl_ins_dc)):
        if i < 1:
            cbl_cum_dc[i]=cbl_ins_dc[i]
            cs1_cum_dc[i]=cs1_ins_dc[i]
            cs2_cum_dc[i]=cs2_ins_dc[i]
            cs3_cum_dc[i]=cs3_ins_dc[i]

            cbl_cum_daly[i]=cbl_dalys[i]
            cs1_cum_daly[i]=cs1_dalys[i]
            cs2_cum_daly[i]=cs2_dalys[i]
            cs3_cum_daly[i]=cs3_dalys[i]
        else:
            cbl_cum_dc[i]=cbl_cum_dc[i-1]+cbl_ins_dc[i]
            cs1_cum_dc[i]=cs1_cum_dc[i-1]+cs1_ins_dc[i]
            cs2_cum_dc[i]=cs2_cum_dc[i-1]+cs2_ins_dc[i]
            cs3_cum_dc[i]=cs3_cum_dc[i-1]+cs3_ins_dc[i]

            cbl_cum_daly[i]=cbl_cum_daly[i-1]+cbl_dalys[i]
            cs1_cum_daly[i]=cs1_cum_daly[i-1]+cs1_dalys[i]
            cs2_cum_daly[i]=cs2_cum_daly[i-1]+cs2_dalys[i]
            cs3_cum_daly[i]=cs3_cum_daly[i-1]+cs3_dalys[i]

    cs1_icer=-(cbl_cum_dc-cs1_cum_dc)/(cbl_cum_daly-cs1_cum_daly)
    cs2_icer=-(cbl_cum_dc-cs2_cum_dc)/(cbl_cum_daly-cs2_cum_daly)
    cs3_icer=-(cbl_cum_dc-cs3_cum_dc)/(cbl_cum_daly-cs3_cum_daly)

    cs1_icer=np.nan_to_num(cs1_icer, 0)
    cs2_icer=np.nan_to_num(cs2_icer, 0)
    cs3_icer=np.nan_to_num(cs3_icer, 0)

    bl_ins_dc=dmc_psa_bl+imc_psa_bl
    s1_ins_dc=dmc_psa_s1+imc_psa_s1
    s2_ins_dc=dmc_psa_s2+imc_psa_s2
    s3_ins_dc=dmc_psa_s3+imc_psa_s3

    bl_cum_dc=np.zeros(np.shape(bl_ins_dc))
    s1_cum_dc=np.zeros(np.shape(s1_ins_dc))
    s2_cum_dc=np.zeros(np.shape(s2_ins_dc))
    s3_cum_dc=np.zeros(np.shape(s3_ins_dc))

    bl_cum_daly=np.zeros(np.shape(bl_dalys))
    s1_cum_daly=np.zeros(np.shape(s1_dalys))
    s2_cum_daly=np.zeros(np.shape(s2_dalys))
    s3_cum_daly=np.zeros(np.shape(s3_dalys))

    for i in range(len(bl_ins_dc)):
        if i < 1:
            bl_cum_dc[i,:]=bl_ins_dc[i,:]
            s1_cum_dc[i,:]=s1_ins_dc[i,:]
            s2_cum_dc[i,:]=s2_ins_dc[i,:]
            s3_cum_dc[i,:]=s3_ins_dc[i,:]

            bl_cum_daly[i,:]=bl_dalys[i,:]
            s1_cum_daly[i,:]=s1_dalys[i,:]
            s2_cum_daly[i,:]=s2_dalys[i,:]
            s3_cum_daly[i,:]=s3_dalys[i,:]
        else:
            bl_cum_dc[i,:]=bl_cum_dc[i-1,:]+bl_ins_dc[i,:]
            s1_cum_dc[i,:]=s1_cum_dc[i-1,:]+s1_ins_dc[i,:]
            s2_cum_dc[i,:]=s2_cum_dc[i-1,:]+s2_ins_dc[i,:]
            s3_cum_dc[i,:]=s3_cum_dc[i-1,:]+s3_ins_dc[i,:]

            bl_cum_daly[i,:]=bl_cum_daly[i-1,:]+bl_dalys[i,:]
            s1_cum_daly[i,:]=s1_cum_daly[i-1,:]+s1_dalys[i,:]
            s2_cum_daly[i,:]=s2_cum_daly[i-1,:]+s2_dalys[i,:]
            s3_cum_daly[i,:]=s3_cum_daly[i-1,:]+s3_dalys[i,:]

    s1_icer=-(bl_cum_dc-s1_cum_dc)/(bl_cum_daly-s1_cum_daly)
    s2_icer=-(bl_cum_dc-s2_cum_dc)/(bl_cum_daly-s2_cum_daly)
    s3_icer=-(bl_cum_dc-s3_cum_dc)/(bl_cum_daly-s3_cum_daly)

    s1_icer=np.nan_to_num(s1_icer, 0)
    s2_icer=np.nan_to_num(s2_icer, 0)
    s3_icer=np.nan_to_num(s3_icer, 0)

    econ_dict={}
    #Year by year DALYs
    econ_dict["bl_daly"]=bl_dalys
    econ_dict["s1_daly"]=s1_dalys
    econ_dict["s2_daly"]=s2_dalys
    econ_dict["s3_daly"]=s3_dalys

    econ_dict["cbl_daly"]=cbl_dalys
    econ_dict["cs1_daly"]=cs1_dalys
    econ_dict["cs2_daly"]=cs2_dalys
    econ_dict["cs3_daly"]=cs3_dalys

    #Year by year intervention costs
    econ_dict["bl_intc"]=dmc_psa_bl
    econ_dict["s1_intc"]=dmc_psa_s1
    econ_dict["s2_intc"]=dmc_psa_s2
    econ_dict["s3_intc"]=dmc_psa_s3

    econ_dict["cbl_intc"]=dmc_cent_bl
    econ_dict["cs1_intc"]=dmc_cent_s1
    econ_dict["cs2_intc"]=dmc_cent_s2
    econ_dict["cs3_intc"]=dmc_cent_s3

    #Year by year medical costs
    econ_dict["bl_medc"]=imc_psa_bl
    econ_dict["s1_medc"]=imc_psa_s1
    econ_dict["s2_medc"]=imc_psa_s2
    econ_dict["s3_medc"]=imc_psa_s3

    econ_dict["cbl_medc"]=imc_cent_bl
    econ_dict["cs1_medc"]=imc_cent_s1
    econ_dict["cs2_medc"]=imc_cent_s2
    econ_dict["cs3_medc"]=imc_cent_s3

    #Year by year productivity costs
    econ_dict["bl_prod"]=bl_prod
    econ_dict["s1_prod"]=s1_prod
    econ_dict["s2_prod"]=s2_prod
    econ_dict["s3_prod"]=s3_prod

    econ_dict["cbl_prod"]=cbl_prod
    econ_dict["cs1_prod"]=cs1_prod
    econ_dict["cs2_prod"]=cs2_prod
    econ_dict["cs3_prod"]=cs3_prod

   #Year by year direct costs (for plot)
    econ_dict["bl_dirc"]=bl_ins_dc
    econ_dict["s1_dirc"]=s1_ins_dc
    econ_dict["s2_dirc"]=s2_ins_dc
    econ_dict["s3_dirc"]=s3_ins_dc

    econ_dict["cbl_dirc"]=cbl_ins_dc
    econ_dict["cs1_dirc"]=cs1_ins_dc
    econ_dict["cs2_dirc"]=cs2_ins_dc
    econ_dict["cs3_dirc"]=cs3_ins_dc

    #Year by year total costs (just in case)
    econ_dict["bl_tcos"]=bl_ins_tc
    econ_dict["s1_tcos"]=s1_ins_tc
    econ_dict["s2_tcos"]=s2_ins_tc
    econ_dict["s3_tcos"]=s3_ins_tc

    econ_dict["cbl_tcos"]=cbl_ins_tc
    econ_dict["cs1_tcos"]=cs1_ins_tc
    econ_dict["cs2_tcos"]=cs2_ins_tc
    econ_dict["cs3_tcos"]=cs3_ins_tc

    #Year by year Net Economic Benefit
    econ_dict["s1_neb"]=s1_neb
    econ_dict["s2_neb"]=s2_neb
    econ_dict["s3_neb"]=s3_neb

    econ_dict["cs1_neb"]=cs1_neb
    econ_dict["cs2_neb"]=cs2_neb
    econ_dict["cs3_neb"]=cs3_neb

    # Year by year ICERs
    econ_dict["s1_icer"]=s1_icer
    econ_dict["s2_icer"]=s2_icer
    econ_dict["s3_icer"]=s3_icer

    econ_dict["cs1_icer"]=cs1_icer
    econ_dict["cs2_icer"]=cs2_icer
    econ_dict["cs3_icer"]=cs3_icer

    #Cumulative Direct Costs (for global ICERs)
    econ_dict["bl_cumdc"]=bl_cum_dc
    econ_dict["s1_cumdc"]=s1_cum_dc
    econ_dict["s2_cumdc"]=s2_cum_dc
    econ_dict["s3_cumdc"]=s3_cum_dc

    econ_dict["cbl_cumdc"]=cbl_cum_dc
    econ_dict["cs1_cumdc"]=cs1_cum_dc
    econ_dict["cs2_cumdc"]=cs2_cum_dc
    econ_dict["cs3_cumdc"]=cs3_cum_dc


    #Cumulative DALYs (for ICERS)
    econ_dict["bl_cdal"]=bl_cum_daly
    econ_dict["s1_cdal"]=s1_cum_daly
    econ_dict["s2_cdal"]=s2_cum_daly
    econ_dict["s3_cdal"]=s3_cum_daly

    econ_dict["cbl_cdal"]=cbl_cum_daly
    econ_dict["cs1_cdal"]=cs1_cum_daly
    econ_dict["cs2_cdal"]=cs2_cum_daly
    econ_dict["cs3_cdal"]=cs3_cum_daly

    return econ_dict

#%% Intervention Coverage and Epidemiological Outcomes Plots -> convert this to plot a single set of results -> have another function to plot together
def epi_plots(cbl, cs1, cs2, cs3, bl, s1, s2, s3, wd, calib, reg):

    from matplotlib import pyplot as plt
    import matplotlib.ticker as mtick
    import seaborn as sns
    import numpy as np

    plot_bd=pd.read_excel(wd+calib, sheet_name="bd_vax")
    plot_hb3=pd.read_excel(wd+calib, sheet_name="hb3_vax")
    plot_prev=pd.read_excel(wd+calib, sheet_name="pop_prv")
    plot_prev_u5=pd.read_excel(wd+calib, sheet_name="u5_prev")
    plot_death=pd.read_excel(wd+calib, sheet_name="t_deaths")
    plot_hcc=pd.read_excel(wd+calib, sheet_name="hcc_inc")

    colours=["#17202A","#66c2a5", "#fc8d62", "#8da0cb"] #black for baseline, then follow qualitative Set2

    sns.set_theme(style="whitegrid")

    if reg!="global":
    ## Intervention Coverage Plots
        cov_plot=plt.figure(figsize=(15,15)) # remove
    #Birth Dose Vaccine Coverage (caveat that coverage of mAVs and HBIG in the period 2022-2050 matches coverage of HepB-BD)
        bd=cov_plot.add_subplot(2,2,1)
        bd.plot(np.arange(1990,2099,1).T, cbl["bd_cov"]*1e2, color=colours[0], linestyle="dashed", alpha=0.8)
        #bd.plot(np.arange(1990,2099,1).T, np.percentile(bl["bd_cov"],50, axis=1)*1e2, color=colours[0], linestyle="dashed", alpha=0.8)
        bd.fill_between(np.arange(1990,2099,1).T, np.percentile(bl["bd_cov"],2.5, axis=1)*1e2, np.percentile(bl["bd_cov"],97.5, axis=1)*1e2,color=colours[0], alpha=0.2)
        bd.plot(np.arange(1990,2099,1).T, cs1["bd_cov"]*1e2, color=colours[1], alpha=0.8)
        #bd.plot(np.arange(1990,2099,1).T, np.percentile(s1["bd_cov"],50, axis=1)*1e2, color=colours[1], alpha=0.8)
        bd.fill_between(np.arange(1990,2099,1).T, np.percentile(s1["bd_cov"],2.5, axis=1)*1e2, np.percentile(s1["bd_cov"],97.5, axis=1)*1e2,color=colours[1], alpha=0.2)
        bd.plot(np.arange(1990,2099,1).T, cs2["bd_cov"]*1e2, color=colours[2], alpha=0.8)
        #bd.plot(np.arange(1990,2099,1).T, np.percentile(s2["bd_cov"],50, axis=1)*1e2, color=colours[2], alpha=0.8)
        bd.fill_between(np.arange(1990,2099,1).T, np.percentile(s2["bd_cov"],2.5, axis=1)*1e2, np.percentile(s2["bd_cov"],97.5, axis=1)*1e2,color=colours[2], alpha=0.2)
        bd.plot(np.arange(1990,2099,1).T, cs3["bd_cov"]*1e2, color=colours[3], alpha=0.8)
        #bd.plot(np.arange(1990,2099,1).T, np.percentile(s3["bd_cov"],50, axis=1)*1e2, color=colours[3], alpha=0.8)
        bd.fill_between(np.arange(1990,2099,1).T, np.percentile(s3["bd_cov"],2.5, axis=1)*1e2, np.percentile(s3["bd_cov"],97.5, axis=1)*1e2,color=colours[3], alpha=0.2)
        bd.scatter(plot_bd["year"], plot_bd[reg]*1e2)
        bd.set_xlim(2000,2050)
        bd.set_ylim([0,100])
        bd.yaxis.set_major_formatter(mtick.PercentFormatter())
        bd.tick_params(axis="both", which="major", labelsize=14)
        bd.hlines(y=90, xmin=2000, xmax=2050, color="red", alpha=0.6, linestyle=':')
        bd.set_title("Hepatitis B Birth Dose Vaccine", fontdict={"fontsize":20, "fontweight":"roman"})
        bd.set_ylabel("Coverage", fontsize=16)

    #HepB3 coverage
        hb3=cov_plot.add_subplot(2,2,2)
        hb3.plot(np.arange(1990,2099,1).T, cbl["hb3_cov"]*1e2, color=colours[0], linestyle="dashed", alpha=0.8)
        hb3.plot(np.arange(1990,2099,1).T, cbl["hb3_cov"]*1e2, color=colours[0], linestyle="dashed", alpha=0.8)
        #hb3.plot(np.arange(1990,2099,1).T, np.percentile(bl["hb3_cov"],50, axis=1)*1e2, color=colours[0], linestyle="dashed", alpha=0.8)
        hb3.fill_between(np.arange(1990,2099,1).T, np.percentile(bl["hb3_cov"],2.5, axis=1)*1e2, np.percentile(bl["hb3_cov"],97.5, axis=1)*1e2,color=colours[0], alpha=0.2)
        hb3.plot(np.arange(1990,2099,1).T, cs1["hb3_cov"]*1e2, color=colours[1], alpha=0.8)
        #hb3.plot(np.arange(1990,2099,1).T, np.percentile(s1["hb3_cov"],50, axis=1)*1e2, color=colours[1], alpha=0.8)
        hb3.fill_between(np.arange(1990,2099,1).T, np.percentile(s1["hb3_cov"],2.5, axis=1)*1e2, np.percentile(s1["hb3_cov"],97.5, axis=1)*1e2,color=colours[1], alpha=0.2)
        hb3.plot(np.arange(1990,2099,1).T, cs2["hb3_cov"]*1e2, color=colours[2], alpha=0.8)
        #hb3.plot(np.arange(1990,2099,1).T, np.percentile(s2["hb3_cov"],50, axis=1)*1e2, color=colours[2], alpha=0.8)
        hb3.fill_between(np.arange(1990,2099,1).T, np.percentile(s2["hb3_cov"],2.5, axis=1)*1e2, np.percentile(s2["hb3_cov"],97.5, axis=1)*1e2,color=colours[2], alpha=0.2)
        hb3.plot(np.arange(1990,2099,1).T, cs3["hb3_cov"]*1e2, color=colours[3], alpha=0.8)
        #hb3.plot(np.arange(1990,2099,1).T, np.percentile(s3["hb3_cov"],50, axis=1)*1e2, color=colours[3], alpha=0.8)
        hb3.fill_between(np.arange(1990,2099,1).T, np.percentile(s3["hb3_cov"],2.5, axis=1)*1e2, np.percentile(s3["hb3_cov"],97.5, axis=1)*1e2,color=colours[3], alpha=0.2)
        hb3.scatter(plot_bd["year"], plot_hb3[reg]*1e2)
        hb3.set_xlim(2000,2050)
        hb3.set_ylim([0,100])
        hb3.yaxis.set_major_formatter(mtick.PercentFormatter())
        hb3.hlines(y=90, xmin=2000, xmax=2050, color="red", alpha=0.6, linestyle=':')
        hb3.tick_params(axis="both", which="major", labelsize=14)
        hb3.set_title("Infant (three-dose HepB) Vaccine", fontdict={"fontsize":20, "fontweight":"roman"})
        hb3.set_ylabel("Coverage", fontsize=16)

    #Diagnosis Coverage
        dx=cov_plot.add_subplot(2,2,3)
        dx.plot(np.arange(1990,2099,1).T, cbl["dx_prop"]*1e2, color=colours[0], linestyle="dashed", alpha=0.8)
        #dx.plot(np.arange(1990,2099,1).T, np.percentile(bl["dx_prop"],50, axis=1)*1e2, color=colours[0], linestyle="dashed", alpha=0.8)
        dx.fill_between(np.arange(1990,2099,1).T, np.percentile(bl["dx_prop"],2.5, axis=1)*1e2, np.percentile(bl["dx_prop"],97.5, axis=1)*1e2,color=colours[0], alpha=0.2)
        dx.plot(np.arange(1990,2099,1).T, cs1["dx_prop"]*1e2, color=colours[1], alpha=0.8)
        #dx.plot(np.arange(1990,2099,1).T, np.percentile(s1["dx_prop"],50, axis=1)*1e2, color=colours[1], alpha=0.8)
        dx.fill_between(np.arange(1990,2099,1).T, np.percentile(s1["dx_prop"],2.5, axis=1)*1e2, np.percentile(s1["dx_prop"],97.5, axis=1)*1e2,color=colours[1], alpha=0.2)
        dx.plot(np.arange(1990,2099,1).T, cs2["dx_prop"]*1e2, color=colours[2], alpha=0.8)
        #dx.plot(np.arange(1990,2099,1).T, np.percentile(s2["dx_prop"],50, axis=1)*1e2, color=colours[2], alpha=0.8)
        dx.fill_between(np.arange(1990,2099,1).T, np.percentile(s2["dx_prop"],2.5, axis=1)*1e2, np.percentile(s2["dx_prop"],97.5, axis=1)*1e2,color=colours[2], alpha=0.2)
        dx.plot(np.arange(1990,2099,1).T, cs3["dx_prop"]*1e2, color=colours[3], alpha=0.8)
        #dx.plot(np.arange(1990,2099,1).T, np.percentile(s3["dx_prop"],50, axis=1)*1e2, color=colours[3], alpha=0.8)
        dx.fill_between(np.arange(1990,2099,1).T, np.percentile(s3["dx_prop"],2.5, axis=1)*1e2, np.percentile(s3["dx_prop"],97.5, axis=1)*1e2,color=colours[3], alpha=0.2)
        dx.hlines(y=90, xmin=2000, xmax=2050, color="red", alpha=0.6, linestyle=':')
        #dx.scatter(2021, plot_care[reg].iloc[3]*1e2, color="blue",marker="X", s=50)
        dx.tick_params(axis="both", which="major", labelsize=14)
        dx.set_xlim([2000,2050])
        dx.set_ylim([0,100])
        dx.yaxis.set_major_formatter(mtick.PercentFormatter())
        dx.set_title("Chronic Hepatitis B (Diagnosed)", fontdict={"fontsize":20, "fontweight":"roman"})
        dx.set_ylabel("Coverage", fontsize=16)

    #Treatment Coverage
        tx=cov_plot.add_subplot(2,2,4)
        tx.plot(np.arange(1990,2099,1).T, cbl["tx_prop"]*1e2, color=colours[0], linestyle="dashed", alpha=0.8)
        tx.plot(np.arange(1990,2099,1).T, cs1["tx_prop"]*1e2, color=colours[1], alpha=0.8)
        tx.plot(np.arange(1990,2099,1).T, cs2["tx_prop"]*1e2, color=colours[2], alpha=0.8)
        tx.plot(np.arange(1990,2099,1).T, cs3["tx_prop"]*1e2, color=colours[3], alpha=0.8)
        #tx.plot(np.arange(1990,2099,1).T, np.percentile(bl["tx_prop"],50, axis=1)*1e2, color=colours[0], linestyle="dashed", alpha=0.8)
        #tx.plot(np.arange(1990,2099,1).T, np.percentile(s1["tx_prop"],50, axis=1)*1e2, color=colours[1], alpha=0.8)
        #tx.plot(np.arange(1990,2099,1).T, np.percentile(s2["tx_prop"],50, axis=1)*1e2, color=colours[2], alpha=0.8)
        #tx.plot(np.arange(1990,2099,1).T, np.percentile(s3["tx_prop"],50, axis=1)*1e2, color=colours[3], alpha=0.8)
        tx.fill_between(np.arange(1990,2099,1).T, np.percentile(bl["tx_prop"],2.5, axis=1)*1e2, np.percentile(bl["tx_prop"],97.5, axis=1)*1e2,color=colours[0], alpha=0.2)
        tx.fill_between(np.arange(1990,2099,1).T, np.percentile(s1["tx_prop"],2.5, axis=1)*1e2, np.percentile(s1["tx_prop"],97.5, axis=1)*1e2,color=colours[1], alpha=0.2)
        tx.fill_between(np.arange(1990,2099,1).T, np.percentile(s2["tx_prop"],2.5, axis=1)*1e2, np.percentile(s2["tx_prop"],97.5, axis=1)*1e2,color=colours[2], alpha=0.2)
        tx.fill_between(np.arange(1990,2099,1).T, np.percentile(s3["tx_prop"],2.5, axis=1)*1e2, np.percentile(s3["tx_prop"],97.5, axis=1)*1e2,color=colours[3], alpha=0.2)
        tx.hlines(y=80, xmin=2000, xmax=2050, color="red", alpha=0.6, linestyle=':')
        #tx.scatter(2021, plot_care[reg].iloc[4]*1e2, color="blue", marker="X", s=50) #add marker for calibration tx coverage
        #tx.scatter(2021, np.percentile((bl["tx_prop"][32,:],50)*1e2, color="red", marker="X",s=50)  #add modelled equivalent of treatment(n)/diag(n) in 2021
        tx.set_xlim([2000,2050])
        tx.set_ylim([0,100])
        tx.yaxis.set_major_formatter(mtick.PercentFormatter())
        tx.tick_params(axis="both", which="major", labelsize=14)
        tx.set_title("Chronic Hepatitis B (Eligible on Treatment)", fontdict={"fontsize":20, "fontweight":"roman"})
        tx.set_ylabel("Coverage", fontsize=16)
        tx.legend(("Status Quo", "Coverage Targets Achieved in 2030", "Coverage Targets Achieved in 2040", "Coverage Targets Achieved in 2050"), fontsize="large", frameon=False, loc="upper left")

        #plt.subplots_adjust(top=0.96, bottom=0.08, left=0.055, right=0.959, hspace=0.35, wspace=0.35)
        plt.tight_layout()
        plt.savefig("output_plots/cov_plots/"+reg+"_coverage plot_PSA.png", dpi=300)

    ## Epidemiological/Calibration Plots
    psa_plot=plt.figure(figsize=(15,15))

    ## Population Prevalence
    prev=psa_plot.add_subplot(2,2,1)
    prev.plot(np.arange(1990,2099,1).T, cbl["prev"]*1e2,color=colours[0], linestyle="dashed", alpha=0.8)
    #prev.plot(np.arange(1990,2099,1).T, np.percentile(bl["prev"],50, axis=1)*1e2,color=colours[0], linestyle="dashed", alpha=0.8)
    prev.fill_between(np.arange(1990,2099,1).T,np.percentile(bl["prev"],2.5, axis=1)*1e2,np.percentile(bl["prev"],97.5, axis=1)*1e2, color=colours[0], alpha=0.2)
    prev.plot(np.arange(1990,2099,1).T, cs1["prev"]*1e2,color=colours[1], alpha=0.8)
    #prev.plot(np.arange(1990,2099,1).T, np.percentile(s1["prev"],50, axis=1)*1e2,color=colours[1], alpha=0.8)
    prev.fill_between(np.arange(1990,2099,1).T,np.percentile(s1["prev"],2.5, axis=1)*1e2,np.percentile(s1["prev"],97.5, axis=1)*1e2, color=colours[1], alpha=0.2)
    prev.plot(np.arange(1990,2099,1).T, cs2["prev"]*1e2,color=colours[2], alpha=0.8)
    #prev.plot(np.arange(1990,2099,1).T, np.percentile(s2["prev"],50, axis=1)*1e2,color=colours[2], alpha=0.8)
    prev.fill_between(np.arange(1990,2099,1).T,np.percentile(s2["prev"],2.5, axis=1)*1e2,np.percentile(s2["prev"],97.5, axis=1)*1e2, color=colours[2], alpha=0.2)
    prev.plot(np.arange(1990,2099,1).T, cs3["prev"]*1e2,color=colours[3], alpha=0.8)
    #prev.plot(np.arange(1990,2099,1).T, np.percentile(s3["prev"],50, axis=1)*1e2,color=colours[3], alpha=0.8)
    prev.fill_between(np.arange(1990,2099,1).T,np.percentile(s3["prev"],2.5, axis=1)*1e2,np.percentile(s3["prev"],97.5, axis=1)*1e2, color=colours[3], alpha=0.2)
    prev.scatter(plot_prev["year"], plot_prev[reg]*1e2)
    prev.set_ylim([0,(max(np.percentile(bl["prev"],97.5, axis=1)*1e2))+0.5])
    prev.set_xlim(2000,2050)
    prev.yaxis.set_major_formatter(mtick.PercentFormatter())
    prev.tick_params(axis="both", which="major", labelsize=14)
    prev.set_title("Population Prevalence", fontdict={"fontsize":20, "fontweight":"roman"})
    prev.set_ylabel("HBsAg seroprevalence", fontsize=16)
    #plt.grid()

    ## Under5y Prevalence
    prev_u5=psa_plot.add_subplot(2,2,2)
    prev_u5.plot(np.arange(1990,2099,1).T, cbl["prev_u5"]*1e2,color=colours[0], linestyle="dashed", alpha=0.8)
    #prev_u5.plot(np.arange(1990,2099,1).T, np.percentile(bl["prev_u5"],50, axis=1)*1e2,color=colours[0], linestyle="dashed", alpha=0.8)
    prev_u5.fill_between(np.arange(1990,2099,1).T,np.percentile(bl["prev_u5"],2.5, axis=1)*1e2,np.percentile(bl["prev_u5"],97.5, axis=1)*1e2, color=colours[0], alpha=0.2)
    prev_u5.plot(np.arange(1990,2099,1).T, cs1["prev_u5"]*1e2,color=colours[1], alpha=0.8)
    #prev_u5.plot(np.arange(1990,2099,1).T, np.percentile(s1["prev_u5"],50, axis=1)*1e2,color=colours[1], alpha=0.8)
    prev_u5.fill_between(np.arange(1990,2099,1).T,np.percentile(s1["prev_u5"],2.5, axis=1)*1e2,np.percentile(s1["prev_u5"],97.5, axis=1)*1e2, color=colours[1], alpha=0.2)
    prev_u5.plot(np.arange(1990,2099,1).T, cs2["prev_u5"]*1e2,color=colours[2], alpha=0.8)
    #prev_u5.plot(np.arange(1990,2099,1).T, np.percentile(s2["prev_u5"],50, axis=1)*1e2,color=colours[2], alpha=0.8)
    prev_u5.fill_between(np.arange(1990,2099,1).T,np.percentile(s2["prev_u5"],2.5, axis=1)*1e2,np.percentile(s2["prev_u5"],97.5, axis=1)*1e2, color=colours[2], alpha=0.2)
    prev_u5.plot(np.arange(1990,2099,1).T, cs3["prev_u5"]*1e2,color=colours[3], alpha=0.8)
    #prev_u5.plot(np.arange(1990,2099,1).T, np.percentile(s3["prev_u5"],50, axis=1)*1e2,color=colours[3], alpha=0.8)
    prev_u5.fill_between(np.arange(1990,2099,1).T,np.percentile(s3["prev_u5"],2.5, axis=1)*1e2,np.percentile(s3["prev_u5"],97.5, axis=1)*1e2, color=colours[3], alpha=0.2)
    prev_u5.scatter(plot_prev_u5["year"], plot_prev_u5[reg]*1e2)
    prev_u5.set_ylim([0,(max(np.percentile(bl["prev_u5"],97.5, axis=1)*1e2))+0.5])
    prev_u5.set_xlim(2000,2050)
    prev_u5.yaxis.set_major_formatter(mtick.PercentFormatter())
    prev_u5.tick_params(axis="both", which="major", labelsize=14)
    prev_u5.set_title("Under 5y Prevalence", fontdict={"fontsize":20, "fontweight":"roman"})
    prev_u5.set_ylabel("HBsAg seroprevalence", fontsize=16)
    #plt.grid()

    ## HBV Mortality
    mort=psa_plot.add_subplot(2,2,3)
    mort.plot(np.arange(1990,2099,1).T,cbl["mort"]/1e3,color=colours[0], linestyle="dashed", alpha=0.8)
    #mort.plot(np.arange(1990,2099,1).T, np.percentile(bl["mort"],50, axis=1)/1e3,color=colours[0], linestyle="dashed", alpha=0.8)
    mort.fill_between(np.arange(1990,2099,1).T,np.percentile(bl["mort"],2.5, axis=1)/1e3,np.percentile(bl["mort"],97.5, axis=1)/1e3, color=colours[0], alpha=0.2)
    mort.plot(np.arange(1990,2099,1).T, cs1["mort"]/1e3,color=colours[1], alpha=0.8)
    #mort.plot(np.arange(1990,2099,1).T, np.percentile(s1["mort"],50, axis=1)/1e3,color=colours[1], alpha=0.8)
    mort.fill_between(np.arange(1990,2099,1).T,np.percentile(s1["mort"],2.5, axis=1)/1e3,np.percentile(s1["mort"],97.5, axis=1)/1e3, color=colours[1], alpha=0.2)
    mort.plot(np.arange(1990,2099,1).T, cs2["mort"]/1e3,color=colours[2], alpha=0.8)
    #mort.plot(np.arange(1990,2099,1).T, np.percentile(s2["mort"],50, axis=1)/1e3,color=colours[2], alpha=0.8)
    mort.fill_between(np.arange(1990,2099,1).T,np.percentile(s2["mort"],2.5, axis=1)/1e3,np.percentile(s2["mort"],97.5, axis=1)/1e3, color=colours[2], alpha=0.2)
    mort.plot(np.arange(1990,2099,1).T, cs3["mort"]/1e3,color=colours[3], alpha=0.8)
    #mort.plot(np.arange(1990,2099,1).T, np.percentile(s3["mort"],50, axis=1)/1e3,color=colours[3], alpha=0.8)
    mort.fill_between(np.arange(1990,2099,1).T,np.percentile(s3["mort"],2.5, axis=1)/1e3,np.percentile(s3["mort"],97.5, axis=1)/1e3, color=colours[3], alpha=0.2)
    mort.scatter(plot_death["year"], plot_death[reg]/1e3)
    mort.set_xlim(2000,2050)
    #mort.set_ylim([0, np.percentile(bl["mort"][32:61,:],97.5)+500])
    mort.tick_params(axis="both", which="major", labelsize=14)
    mort.set_title("Hepatitis B Attributable Mortality", fontdict={"fontsize":20, "fontweight":"roman"})
    mort.set_ylabel("Annual Deaths (thousand)", fontsize=16)
    #plt.gca().yaxis.set_major_formatter(plt.matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
    #plt.grid()

    ## HCC incidence
    hcc=psa_plot.add_subplot(2,2,4)
    hcc.plot(np.arange(1990,2099,1).T, cbl["hcc_inc"]/1e3,color=colours[0], linestyle="dashed", alpha=0.8)
    hcc.plot(np.arange(1990,2099,1).T, cs1["hcc_inc"]/1e3,color=colours[1], alpha=0.8)
    hcc.plot(np.arange(1990,2099,1).T, cs2["hcc_inc"]/1e3,color=colours[2], alpha=0.8)
    hcc.plot(np.arange(1990,2099,1).T, cs3["hcc_inc"]/1e3,color=colours[3], alpha=0.8)
    #hcc.plot(np.arange(1990,2099,1).T, np.percentile(bl["hcc_inc"],50, axis=1)/1e3,color=colours[0], linestyle="dashed", alpha=0.8)
    #hcc.plot(np.arange(1990,2099,1).T, np.percentile(s1["hcc_inc"],50, axis=1)/1e3,color=colours[1], alpha=0.8)
    #hcc.plot(np.arange(1990,2099,1).T, np.percentile(s2["hcc_inc"],50, axis=1)/1e3,color=colours[2], alpha=0.8)
    #hcc.plot(np.arange(1990,2099,1).T, np.percentile(s3["hcc_inc"],50, axis=1)/1e3,color=colours[3], alpha=0.8)
    hcc.fill_between(np.arange(1990,2099,1).T,np.percentile(bl["hcc_inc"],2.5, axis=1)/1e3,np.percentile(bl["hcc_inc"],97.5, axis=1)/1e3, color=colours[0], alpha=0.2)
    hcc.fill_between(np.arange(1990,2099,1).T,np.percentile(s1["hcc_inc"],2.5, axis=1)/1e3,np.percentile(s1["hcc_inc"],97.5, axis=1)/1e3, color=colours[1], alpha=0.2)
    hcc.fill_between(np.arange(1990,2099,1).T,np.percentile(s2["hcc_inc"],2.5, axis=1)/1e3,np.percentile(s2["hcc_inc"],97.5, axis=1)/1e3, color=colours[2], alpha=0.2)
    hcc.fill_between(np.arange(1990,2099,1).T,np.percentile(s3["hcc_inc"],2.5, axis=1)/1e3,np.percentile(s3["hcc_inc"],97.5, axis=1)/1e3, color=colours[3], alpha=0.2)
    hcc.scatter(plot_hcc["year"], plot_hcc[reg]/1e3)
    hcc.set_xlim(2000,2050)
    #hcc.set_ylim([0, np.percentile(bl["hcc_inc"][32:61,:],97.5)+500])
    hcc.tick_params(axis="both", which="major", labelsize=14)
    hcc.set_title("Hepatocellular Carcinoma Incidence", fontdict={"fontsize":20, "fontweight":"roman"})
    hcc.set_ylabel("Annual Cases (thousand)", fontsize=16)
    plt.gca().yaxis.set_major_formatter(plt.matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
    hcc.legend(("Status Quo", "Coverage Targets Achieved in 2030", "Coverage Targets Achieved in 2040", "Coverage Targets Achieved in 2050"), fontsize="large", frameon=False)
    #plt.grid()

    #plt.subplots_adjust(top=0.96, bottom=0.04, left=0.065, right=0.98, hspace=0.35, wspace=0.35)
    plt.tight_layout()

    plt.savefig("output_plots/epi_plots/"+reg+"_epi plot_PSA.png", dpi=500)

    if reg=="global":
        plt.savefig("output_plots/epi_plots/"+reg+"_epi plot_PSA.pdf", dpi=500)

    return print("Figures Generated!")

#%% Plot Economic Outcomes
def econ_plots (econ_dict, reg):

    from matplotlib import pyplot as plt
    import matplotlib.ticker as mtick
    import seaborn as sns
    import numpy as np

    colours=["#17202A","#66c2a5", "#fc8d62", "#8da0cb"] #black for baseline, then follow qualitative Set2

    sns.set_theme(style="whitegrid")

    econ_plot=plt.figure(figsize=(15,15))
    ## DALYs
    daly=econ_plot.add_subplot(2,2,1)
    daly.plot(np.arange(1990,2099,1).T, econ_dict["cbl_daly"]/1e6,color=colours[0], linestyle="dashed", alpha=0.8)
    daly.plot(np.arange(1990,2099,1).T, econ_dict["cs1_daly"]/1e6,color=colours[1], alpha=0.8)
    daly.plot(np.arange(1990,2099,1).T, econ_dict["cs2_daly"]/1e6,color=colours[2], alpha=0.8)
    daly.plot(np.arange(1990,2099,1).T, econ_dict["cs3_daly"]/1e6,color=colours[3], alpha=0.8)
    # daly.plot(np.arange(1990,2099,1).T, np.percentile(econ_dict["bl_daly"],50, axis=1)/1e6,color=colours[0], linestyle="dashed", alpha=0.8)
    # daly.plot(np.arange(1990,2099,1).T, np.percentile(econ_dict["s1_daly"],50, axis=1)/1e6,color=colours[1], alpha=0.8)
    # daly.plot(np.arange(1990,2099,1).T, np.percentile(econ_dict["s2_daly"],50, axis=1)/1e6,color=colours[2], alpha=0.8)
    # daly.plot(np.arange(1990,2099,1).T, np.percentile(econ_dict["s3_daly"],50, axis=1)/1e6,color=colours[3], alpha=0.8)
    daly.fill_between(np.arange(1990,2099,1).T,np.percentile(econ_dict["bl_daly"],2.5, axis=1)/1e6,np.percentile(econ_dict["bl_daly"],97.5, axis=1)/1e6, color=colours[0], alpha=0.2)
    daly.fill_between(np.arange(1990,2099,1).T,np.percentile(econ_dict["s1_daly"],2.5, axis=1)/1e6,np.percentile(econ_dict["s1_daly"],97.5, axis=1)/1e6, color=colours[1], alpha=0.2)
    daly.fill_between(np.arange(1990,2099,1).T,np.percentile(econ_dict["s2_daly"],2.5, axis=1)/1e6,np.percentile(econ_dict["s2_daly"],97.5, axis=1)/1e6, color=colours[2], alpha=0.2)
    daly.fill_between(np.arange(1990,2099,1).T,np.percentile(econ_dict["s3_daly"],2.5, axis=1)/1e6,np.percentile(econ_dict["s3_daly"],97.5, axis=1)/1e6, color=colours[3], alpha=0.2)
    daly.set_xlim(2022,2050)
    daly.set_title("Disability Adjusted Life Years", fontdict={"fontsize":20, "fontweight":"roman"})
    daly.set_ylabel("DALYs (million)", fontsize=16)
    daly.tick_params(axis="both", which="major", labelsize=14)
    daly.legend(("Status Quo", "Coverage Targets Achieved in 2030", "Coverage Targets Achieved in 2040", "Coverage Targets Achieved in 2050"), fontsize="large", frameon=False)


    ## Direct Costs
    dc=econ_plot.add_subplot(2,2,2)
    dc.plot(np.arange(1990,2099,1).T, econ_dict["cbl_dirc"]/1e9,color=colours[0], linestyle="dashed", alpha=0.8)
    #dc.plot(np.arange(1990,2099,1).T, np.percentile(econ_dict["bl_dirc"],50, axis=1)/1e9,color=colours[0], linestyle="dashed", alpha=0.8)
    dc.fill_between(np.arange(1990,2099,1).T,np.percentile(econ_dict["bl_dirc"],2.5, axis=1)/1e9,np.percentile(econ_dict["bl_dirc"],97.5, axis=1)/1e9, color=colours[0], alpha=0.2)
    dc.plot(np.arange(1990,2099,1).T, econ_dict["cs1_dirc"]/1e9,color=colours[1], alpha=0.8)
    #dc.plot(np.arange(1990,2099,1).T, np.percentile(econ_dict["s1_dirc"],50, axis=1)/1e9,color=colours[1], alpha=0.8)
    dc.fill_between(np.arange(1990,2099,1).T,np.percentile(econ_dict["s1_dirc"],2.5, axis=1)/1e9,np.percentile(econ_dict["s1_dirc"],97.5, axis=1)/1e9, color=colours[1], alpha=0.2)
    dc.plot(np.arange(1990,2099,1).T, econ_dict["cs2_dirc"]/1e9,color=colours[2], alpha=0.8)
    #dc.plot(np.arange(1990,2099,1).T, np.percentile(econ_dict["s2_dirc"],50, axis=1)/1e9,color=colours[2], alpha=0.8)
    dc.fill_between(np.arange(1990,2099,1).T,np.percentile(econ_dict["s2_dirc"],2.5, axis=1)/1e9,np.percentile(econ_dict["s2_dirc"],97.5, axis=1)/1e9, color=colours[2], alpha=0.2)
    dc.plot(np.arange(1990,2099,1).T, econ_dict["cs3_dirc"]/1e9,color=colours[3], alpha=0.8)
    #dc.plot(np.arange(1990,2099,1).T, np.percentile(econ_dict["s3_dirc"],50, axis=1)/1e9,color=colours[3], alpha=0.8)
    dc.fill_between(np.arange(1990,2099,1).T,np.percentile(econ_dict["s3_dirc"],2.5, axis=1)/1e9,np.percentile(econ_dict["s3_dirc"],97.5, axis=1)/1e9, color=colours[3], alpha=0.2)
    dc.set_xlim(2022,2050)
    dc.tick_params(axis="both", which="major", labelsize=14)
    dc.set_title("Healthcare Costs", fontdict={"fontsize":20, "fontweight":"roman"})
    dc.set_ylabel("US$ (billion)", fontsize=16)


    ## ICERs
    icer=econ_plot.add_subplot(2,2,3)
    icer.plot(np.arange(2025,2051,1).T, econ_dict["cs1_icer"][35:61]/1e3,color=colours[1], alpha=0.8)
    #icer.plot(np.arange(2025,2051,1).T, np.percentile(econ_dict["s1_icer"][35:61,:],50, axis=1)/1e3,color=colours[1], alpha=0.8)
    icer.fill_between(np.arange(2025,2051,1).T,np.percentile(econ_dict["s1_icer"][35:61,:],2.5, axis=1)/1e3,np.percentile(econ_dict["s1_icer"][35:61,:],97.5, axis=1)/1e3, color=colours[1], alpha=0.2)
    icer.plot(np.arange(2025,2051,1).T, econ_dict["cs2_icer"][35:61]/1e3,color=colours[2], alpha=0.8)
    #icer.plot(np.arange(2025,2051,1).T, np.percentile(econ_dict["s2_icer"][35:61,:],50, axis=1)/1e3,color=colours[2], alpha=0.8)
    icer.fill_between(np.arange(2025,2051,1).T,np.percentile(econ_dict["s2_icer"][35:61,:],2.5, axis=1)/1e3,np.percentile(econ_dict["s2_icer"][35:61,:],97.5, axis=1)/1e3, color=colours[2], alpha=0.2)
    icer.plot(np.arange(2025,2051,1).T, econ_dict["cs3_icer"][35:61]/1e3,color=colours[3], alpha=0.8)
    #icer.plot(np.arange(2025,2051,1).T, np.percentile(econ_dict["s3_icer"][35:61,:],50, axis=1)/1e3,color=colours[3], alpha=0.8)
    icer.fill_between(np.arange(2025,2051,1).T,np.percentile(econ_dict["s3_icer"][35:61,:],2.5, axis=1)/1e3,np.percentile(econ_dict["s3_icer"][35:61,:],97.5, axis=1)/1e3, color=colours[3], alpha=0.2)
    icer.set_xlim(2025,2050)
    icer.tick_params(axis="both", which="major", labelsize=14)
    #icer.set_ylim(0,8e4)
    icer.set_title("ICER", fontdict={"fontsize":20, "fontweight":"roman"})
    icer.set_ylabel("US$ (thousand) per DALY averted", fontsize=16)

    ## Net Economic Benefit
    neb=econ_plot.add_subplot(2,2,4)
    neb.plot(np.arange(2022,2051,1).T, econ_dict["cs1_neb"][32:61]/1e9,color=colours[1], alpha=0.8)
    #neb.plot(np.arange(2022,2051,1).T, np.percentile(econ_dict["s1_neb"][32:61,:],50, axis=1)/1e9,color=colours[1], alpha=0.8)
    neb.fill_between(np.arange(2022,2051,1).T,np.percentile(econ_dict["s1_neb"][32:61,:],2.5, axis=1)/1e9,np.percentile(econ_dict["s1_neb"][32:61,:],97.5, axis=1)/1e9, color=colours[1], alpha=0.2)
    neb.plot(np.arange(2022,2051,1).T,econ_dict["cs2_neb"][32:61]/1e9,color=colours[2], alpha=0.8)
    #neb.plot(np.arange(2022,2051,1).T, np.percentile(econ_dict["s2_neb"][32:61,:],50, axis=1)/1e9,color=colours[2], alpha=0.8)
    neb.fill_between(np.arange(2022,2051,1).T,np.percentile(econ_dict["s2_neb"][32:61,:],2.5, axis=1)/1e9,np.percentile(econ_dict["s2_neb"][32:61,:],97.5, axis=1)/1e9, color=colours[2], alpha=0.2)
    neb.plot(np.arange(2022,2051,1).T, econ_dict["cs3_neb"][32:61]/1e9,color=colours[3], alpha=0.8)
    #neb.plot(np.arange(2022,2051,1).T, np.percentile(econ_dict["s3_neb"][32:61,:],50, axis=1)/1e9,color=colours[3], alpha=0.8)
    neb.fill_between(np.arange(2022,2051,1).T,np.percentile(econ_dict["s3_neb"][32:61,:],2.5, axis=1)/1e9,np.percentile(econ_dict["s3_neb"][32:61,:],97.5, axis=1)/1e9, color=colours[3], alpha=0.2)
    neb.set_xlim(2022,2050)
    neb.tick_params(axis="both", which="major", labelsize=14)
    neb.set_title("Net Economic Benefit", fontdict={"fontsize":20, "fontweight":"roman"})
    neb.set_ylabel("US$ (billion)", fontsize=16)

    #plt.subplots_adjust(top=0.96, bottom=0.04, left=0.065, right=0.98, hspace=0.35, wspace=0.35)
    plt.tight_layout()

    if reg=="global":
        plt.savefig("econ_plots/"+reg+"_econ plot PSA.pdf", dpi=500)

    plt.savefig("econ_plots/"+reg+"_econ plot PSA.png", dpi=300)



    return print("Figure Generated!")

def scenario_edit(D, par, cov_dict = defaultdict(lambda: []), t_dict = defaultdict(lambda: []), assump_dict = defaultdict(lambda: None)):
    '''
    A single edit of one variable in a databook, driven by dictionaries
    :param D:             Project Databook to be edited (atomica Databook)
    :param par:           Specific parameter that will be edited (str)
    :param cov_dict:      Dictionary of new time series values (dict; keys -> population, values -> list of new input variables)
    :param t_dict:        Dictionary of new time series year values (dict; keys -> population, values -> list of relative years of input variables)
    :param assump_dict:   Dictionary of new assumption values (dict; keys -> population, values -> assumption)
    :return D:			  New Edited databook (databook)
    '''
    assert cov_dict.keys() == t_dict.keys(), 'cov and t dict must have the same pops!'
    tmp = D.tdve[par].ts['0-4M'].units
    
    for pop in cov_dict.keys():
        D.tdve[par].ts[pop] = at.TimeSeries(t_dict[pop], cov_dict[pop], units = tmp)
        
    for pop in assump_dict.keys():
        D.tdve[par].ts[pop] = at.TimeSeries([], [], units = tmp, assumption = assump_dict[pop])
        
    return D

def no_vax_scenario(db_loc, fw_loc = "hbv_v14_gamma_mav.xlsx", keep_db = False):
    '''
	Introducing the same scenario (from a databook) to simulate a simple no-vaccination scenario
	:param db_loc: location of the databook to be edited (str)
	:param fw_loc: location of the atomica framework (str)
	:param keep_db: whether or not you want to keep the new databook (boolean)
	:return D: New Edited databook (databook)
    '''
    # adjust the two vaccination components to zero (bd, hb3)
    assump_dict = {}

    for sx in ['M', 'F']:
        assump_dict[f'0-4{sx}'] = 0
    
    F = at.ProjectFramework("hbv_v14_gamma_mav.xlsx")
    D = at.ProjectData.from_spreadsheet(db_loc, framework = F)
    D = scenario_edit(D, 'bd', assump_dict = assump_dict)
    D = scenario_edit(D, 'hb3', assump_dict = assump_dict)
    D.save('tmp_databook.xlsx')
    D = at.ProjectData.from_spreadsheet('tmp_databook.xlsx', framework = F) # We need to reload the databook
    
    
    
    # Remove the temporary databook (you can choose to keep it too by commenting this out)
    if keep_db == False:
        os.remove('tmp_databook.xlsx')
    # save new db and reload it to return it!
    return D
    
def no_treat_scenario(db_loc, fw_loc = "hbv_v14_gamma_mav.xlsx", keep_db = False):
    """
	Introducing the same scenario (from a databook) to simulate a simple no-treatment scenario
	:param db_loc: location of the databook to be edited (str)
	:param fw_loc: location of the atomica framework (str)
	:param keep_db: whether or not you want to keep the new databook (boolean)
	:return D: New Edited databook (databook)
    """
    F = at.ProjectFramework("hbv_v14_gamma_mav.xlsx")
    D = at.ProjectData.from_spreadsheet(db_loc, framework = F)
    assump_dict = {}
    
    for pop in D.pops:
        assump_dict[pop] = 0
    
    D = scenario_edit(D, 't_cov_net', assump_dict = assump_dict)
    D = scenario_edit(D, 't_cov_hb', assump_dict = assump_dict)
    D.save('tmp_databook.xlsx')
    D = at.ProjectData.from_spreadsheet('tmp_databook.xlsx', framework = F)
    
    # Remove the temporary databook (you can choose to keep it too by commenting this out)
    if keep_db == False:
        os.remove('tmp_databook.xlsx')

    return D # save new db and reload it to return it!

def pop_calib(P,cal, maxtime = 300):
    '''
    Assuming acm and alive are still in the framework
    :param P: atomica Project to be calibrated
    :param cal: existing atomica Project calibration
    :param maxtime: the maximum time allowed for the calibration
    :return cal: updated population calibration
    '''
    return P.calibrate(max_time=maxtime, parset=cal, adjustables=["acm"], measurables=["alive"], default_min_scale=0.5, default_max_scale=2)

def epi_calib(P, cal,maxtime = 100, ranges = 4):
    '''
    Calibrating general epidemiology of the model (acute, chronic, cirrhosis, and hcc)

    '''
    for i in range(ranges): #lets iterate a few more times
        print(f'RUN NUMBER {i+1}, acute and chronic')
        cal = P.calibrate(max_time=maxtime, parset=cal, adjustables=["sag_ix", "ad_pop_sus", "ch_pop_sus", "eag_ix"], measurables=["prev"], default_min_scale=0.1, default_max_scale=5)
        cal = P.calibrate(max_time=maxtime, parset=cal, adjustables=["m_acu"], measurables=["cl_acu"], default_min_scale=0.1, default_max_scale=5)
        cal = P.calibrate(max_time=maxtime, parset=cal, adjustables=["it_icl", "icl_ict"], measurables=["eag_ott"], default_min_scale=0.1, default_max_scale=3)

        print(f'RUN NUMBER {i+1}, cirr and hcc')
        # Testing this new parameter (lower and upper bounds can be adjusted)
        cal = P.calibrate(max_time=maxtime, parset=cal, adjustables=["cc_dc", "dc_hcc"], measurables=["cl_cir", 'cl_hcc'], default_min_scale=0.1, default_max_scale=3)
        cal = P.calibrate(max_time=maxtime, parset=cal, adjustables=["m_dc"], measurables=["cl_cir"], default_min_scale=0.1, default_max_scale=5)
        cal = P.calibrate(max_time=maxtime, parset=cal, adjustables=["m_hcc"], measurables=["cl_hcc"], default_min_scale=0.1, default_max_scale=5)
    return cal

def highlight_results(s, lb, ub, props=""):
    return np.where((lb <= s) & (s < ub), props, "")

def determine_cal_quality(df):
    """
    Created by: Phillip Luong
    Last Updated: 26/04/2023
    Taken from cat-typhoid
    Needs highlight_results, helps determine how good the calibrations are after generate_output_df
    """

    df = df.style.apply(highlight_results, lb=1, ub=1000, props="color:white;background-color:black",  axis=0).apply(highlight_results, lb=0.6, ub=1, props="color:white;background-color:red",axis=0).apply(highlight_results, lb=0.4, ub=0.6, props="color:black;background-color:gold",  axis=0).apply(highlight_results, lb=0, ub=0.4, props="color:black;background-color:mediumseagreen",  axis=0)
    return df

def metric(data, est, mt):
    """
    LAST EDIT: 10/08/21
    by: Phillip Luong

    :param data: An array or list of data points to calibrate to
    :param est: An array or list of the estimated points
    :param mt: 'mae', 'mse', 'rmse', 'mape'
    :return: Numercal value for the error
    """
    if mt == "mae":
        return mean_absolute_error(data, est, sample_weight=data)
    if mt == "mse":
        return mean_squared_error(data, est, sample_weight=data)
    if mt == "rmse":
        return mean_squared_error(data, est, sample_weight=data, squared=False)
    if mt == "mape":
        return mean_absolute_percentage_error(data, est, sample_weight=data)

def error_table(regions, fw_loc = 'frameworks/hbv_v14_gamma_mav.xlsx'):
	'''
	THIS ONLY WORKS WITH WHO REGIONS RIGHT NOW
	TODO: comment on how parts of the funtion work
	'''
	edf = pd.DataFrame()
	years = [2010, 2015, 2019]
	pars = ['chb_pop', 'cl_acu', 'cl_cir', 'cl_hcc']
	F = at.ProjectFramework(fw_loc)
	for ct in regions:
	    print(f"Current region: {ct}")

	    D = at.ProjectData.from_spreadsheet(f"applications/region_{ct.lower()}/{ct}_db_mav.xlsx", framework=F)
	    P= at.Project(framework=F, databook=D, sim_dt=0.25, sim_start=1990, sim_end=2099, do_run=False)
	    cal = P.make_parset()
	    cal.load_calibration(f"applications/region_{ct.lower()}/{ct}_calib.xlsx")
	    res = P.run_sim(parset=cal)

	    for par in pars:
	        data = []
	        est = []

	        if par == 'chb_pop':
	            tmp_yr = years[1:]
	        else:
	            tmp_yr = years

	        for yr in tmp_yr:
	            tmp = [P.data.tdve[par].ts[x].interpolate(2000 + yr) for x in range(10)]  # data
	            data = [i[0] for i in tmp]
	            est = [x.vals[(yr-2000) * 4] for x in res.get_variable(par)]  # estimate

	            edf.loc[f"{ct}", par] = metric(data, est, "mape")
	    
	edf = edf.apply(pd.to_numeric).round(4)
	return edf