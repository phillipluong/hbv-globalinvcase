## Code to reproduce analyses presented input
## A global investment case for hepatitis B elimination – a modelling study
## Seaman et al, 2023

#%% Import pre-requisites
import os
wd= 'C:/Users/Phil/Documents/GitHub/hbv-globalinvcase'#set own working directory
os.chdir(wd)

import atomica as at
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def model_results(F, data, calib, res_name, runs):

    """ Runs the model and simulations, returns central and sampled data sets for analysis"""

    import atomica as at
    import numpy as np

    ## Run the model
    D=at.ProjectData.from_spreadsheet("databooks/"+data, framework=F)
    P=at.Project(framework=F, databook="databooks/"+data, sim_start=1990, sim_end=2099, sim_dt=0.25, do_run=False)
    cal=P.make_parset()
    cal.load_calibration("calibrations/"+calib)

    # Central Estimate (can be used if needed, but median of ~100 runs converges very well [i.e., no meaningful difference])
    res=P.run_sim(parset=cal, result_name = res_name) # can also expand code here to check calibrations if needed

    #Probabilistic sensitivity analysis
    np.random.seed(25012023)
    psa=P.parsets[0]
    psa_sample=psa.sample()
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

#%% Economic Outcome Calculations
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

#%% Intervention Coverage and Epidemiological Outcomes Plots
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
        cov_plot=plt.figure(figsize=(15,15))
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

    #Diagnsosis Coverage
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
        plt.savefig("cov_plots/"+reg+"_coverage plot_PSA.png", dpi=300)

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

    plt.savefig("epi_plots/"+reg+"_epi plot PSA.png", dpi=500)

    if reg=="global":
        plt.savefig("epi_plots/"+reg+"_epi plot PSA.pdf", dpi=500)

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

def par_sens(F, db_bl, db_s1, calib, par_dict):

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

    D_bl=at.ProjectData.from_spreadsheet("databooks/"+db_bl, framework=F)
    P_bl=at.Project(framework=F, databook="databooks/"+db_bl, sim_start=1990, sim_end=2099, sim_dt=0.25, do_run=False)
    cal_bl=P_bl.make_parset()
    cal_bl.load_calibration("calibrations/"+calib)

    bl=P_bl.run_sim(parset=cal_bl, result_name = "Status Quo")

    par_scen_bl=at.ParameterScenario(name="scen_bl")

    for i in par_dict:
        for key,val in i.items():
            par_scen_bl.add(val[0], val[1], [val[2], val[3], val[4]], [val[5], val[6], val[7]])

    par_scen_bl_res=par_scen_bl.run(P_bl, parset=cal_bl)

    central_est_bl={}
    for i in out_res:
        for key, val in i.items():
            df=at.PlotData(par_scen_bl_res, outputs=val[0], pops=val[1], pop_aggregation=val[2], t_bins=1).series[0].vals
            central_est_bl[key]=df


    D_s1=at.ProjectData.from_spreadsheet("databooks/"+db_s1, framework=F)
    P_s1=at.Project(framework=F, databook="databooks/"+db_s1, sim_start=1990, sim_end=2099, sim_dt=0.25, do_run=False)
    cal_s1=P_s1.make_parset()
    cal_s1.load_calibration("calibrations/"+calib)

    s1=P_s1.run_sim(parset=cal_s1, result_name = "Scenario 1: 2030 Target")

    par_scen_s1=at.ParameterScenario(name="scen_s1")

    for i in par_dict:
        for key,val in i.items():
            par_scen_s1.add(val[0], val[1], [val[2], val[3], val[4]], [val[5], val[6], val[7]])

    par_scen_s1_res=par_scen_s1.run(P_s1, parset=cal_s1)

    central_est_s1={}
    for i in out_res:
        for key, val in i.items():
            df=at.PlotData(par_scen_s1_res, outputs=val[0], pops=val[1], pop_aggregation=val[2], t_bins=1).series[0].vals
            central_est_s1[key]=df

    return central_est_bl, central_est_s1

def parsens_econ(psens_bl, psens_s1, reg, wd, cost_data, disc_rate):

    # psens_bl=amr_mtctw_bl
    # psens_s1=amr_mtctw_s1
    # reg="AMR"
    # disc_rate=0.03
    # cost_data="cost and agg data/costs.xlsx"

    #Discounting Array
    if disc_rate >1:
        disc_rate=disc_rate/100
    else:
        disc_rate=disc_rate

    discount=np.zeros((len(np.arange(1990,2099,1)),3))
    discount[:,0]=np.arange(1990,2099,1)

    for idx,val in enumerate(discount[:,0]):
        if val <= 2021:
            discount[idx,1]=1
            discount[idx,2]=0 #use this for high risk screening as a check
        else:
            discount[idx,1:]=(1-disc_rate)**(val-2022)


    vax_costs=pd.read_excel(cost_data, sheet_name="vax")
    bd_vax=vax_costs[reg].iloc[0]
    hb3_vax=vax_costs[reg].iloc[1]
    hrs_cost=vax_costs[reg].iloc[10]
    mav_cost=vax_costs[reg].iloc[11]

    care_costs=pd.read_excel(cost_data, sheet_name="care")
    dx_cost=care_costs[reg].iloc[0]
    tx_cost=care_costs[reg].iloc[3]
    dx_dmc=care_costs[reg].iloc[4] #cost of disease management (diagnosed)
    tx_dmc=care_costs[reg].iloc[4] #cost of disease management (on treatment)
    hosp_cost=care_costs[reg].iloc[7] #cost of hospitalisation
    hcc_cost=care_costs[reg].iloc[6] #HCC surveillance costs
    hcc_prp=care_costs[reg].iloc[8] #proportion of CHB in 15-49y who are 40-49y (HCC screening)

    ## Vaccines and mAV

    bd_cost_bl, hb3_cost_bl, bd_cost_s1, hb3_cost_s1=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))
    bd_cost_bl[:,0]=psens_bl["bd_cov"][:]*psens_bl["births"][:]*bd_vax*discount[:,1]
    hb3_cost_bl[:,0]=psens_bl["hb3_cov"][:]*psens_bl["births"][:]*hb3_vax*discount[:,1]
    bd_cost_s1[:,0]=psens_s1["bd_cov"][:]*psens_s1["births"][:]*bd_vax*discount[:,1]
    hb3_cost_s1[:,0]=psens_s1["hb3_cov"][:]*psens_s1["births"][:]*hb3_vax*discount[:,1]

    bl_cost_mav, s1_cost_mav=np.zeros((109,1)),np.zeros((109,1))
    bl_cost_mav[:,0]=(psens_bl["prg_scr"][:]*dx_cost*discount[:,1])#(cent["mav_n"][:]*mav_cost*discount[:,1])+(cent["prg_hrs"][:]*hrs_cost*discount[:,1])
    s1_cost_mav[:,0]=(psens_s1["mav_n"][:]*mav_cost*discount[:,1])+(psens_s1["prg_scr"][:]*dx_cost*discount[:,1])+(psens_s1["prg_hrs"][:]*hrs_cost*discount[:,2])

    ## Diagnoses
    bl_dx_inc, s1_dx_inc=np.zeros((109,1)),np.zeros((109,1))

    for i in range(len(psens_bl["dx_prop"])):
        if i < 1:
            bl_dx_inc[i,0]=psens_bl["dx_prop"][i]
            s1_dx_inc[i,0]=psens_s1["dx_prop"][i]
        else:
            bl_dx_inc[i,0]=max(psens_bl["dx_prop"][i]-psens_bl["dx_prop"][i-1],0)
            s1_dx_inc[i,0]=max(psens_s1["dx_prop"][i]-psens_s1["dx_prop"][i-1],0)

    dx_costb_bl, dx_costb_s1=np.zeros((109,1)),np.zeros((109,1))

    for yr in range(len(bl_dx_inc)):
        dx_costb_bl[yr,0]=(psens_bl["dx_rate"][yr]*dx_cost*discount[yr,1])+(dx_cost*bl_dx_inc[yr,0]*(psens_bl["pop"][yr]-psens_bl["dx_rate"][yr])*discount[yr,1])
        dx_costb_s1[yr,0]=(psens_s1["dx_rate"][yr]*dx_cost*discount[yr,1])+(dx_cost*s1_dx_inc[yr,0]*(psens_s1["pop"][yr]-psens_s1["dx_rate"][yr])*discount[yr,1])

    ## Treatment
    bl_cost_tx, s1_cost_tx=np.zeros((109,1)),np.zeros((109,1))

    bl_cost_tx[:,0]=psens_bl["tx_cov"][:]*tx_cost*discount[:,1]
    s1_cost_tx[:,0]=psens_s1["tx_cov"][:]*tx_cost*discount[:,1]

    dmc_bl=bd_cost_bl+hb3_cost_bl+bl_cost_mav+dx_costb_bl+bl_cost_tx
    dmc_s1=bd_cost_s1+hb3_cost_s1+s1_cost_mav+dx_costb_s1+s1_cost_tx

    ## Disease Management
    util=0.25
    tx_hosp=0.5

    manc_bl, manc_s1=np.zeros((109,1)),np.zeros((109,1))
    hospc_bl, hospc_s1=np.zeros((109,1)),np.zeros((109,1))
    hccs_bl, hccs_s1=np.zeros((109,1)),np.zeros((109,1))

    manc_bl[:,0]=((psens_bl["dm_dx"][:]*dx_dmc*util)+(psens_bl["dm_tx"][:]*tx_dmc*util))*discount[:,1]
    manc_s1[:,0]=((psens_s1["dm_dx"][:]*dx_dmc*util)+(psens_s1["dm_tx"][:]*tx_dmc*util))*discount[:,1]

    hospc_bl[:, 0]=((psens_bl["hsp_utx"][:]*hosp_cost*util)+(psens_bl["hsp_tx"][:]*hosp_cost*util*tx_hosp))*discount[:,1]
    hospc_s1[:, 0]=((psens_s1["hsp_utx"][:]*hosp_cost*util)+(psens_s1["hsp_tx"][:]*hosp_cost*util*tx_hosp))*discount[:,1]

    hccs_bl[:, 0]=((psens_bl["pop_hcc"][:]*hcc_cost*util)+(psens_bl["tgt_hcc"][:]*hcc_cost*util)+(psens_bl["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util))*discount[:,1]
    hccs_s1[:, 0]=((psens_s1["pop_hcc"][:]*hcc_cost*util)+(psens_s1["tgt_hcc"][:]*hcc_cost*util)+(psens_s1["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util))*discount[:,1]

    imc_bl=manc_bl+hospc_bl+hccs_bl
    imc_s1=manc_s1+hospc_s1+hccs_s1

    ## Productivity, YLL, YPLL, DALYs
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

    #Baseline
    cbl_deaths=np.zeros((len(psens_bl["mort"]),2))
    cbl_deaths[:,0]=np.arange(1990,2099,1)
    cbl_deaths[:,1]=psens_bl["mort"]

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

    ## Scenario 1 (central)
    cs1_deaths=np.zeros((len(psens_s1["mort"]),2))
    cs1_deaths[:,0]=np.arange(1990,2099,1)
    cs1_deaths[:,1]=psens_s1["mort"]

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

    #DALYs
    bl_yll=np.sum(ghosts_cbl[:,0:5], axis=1)*discount[:,1]
    s1_yll=np.sum(ghosts_cs1[:,0:5], axis=1)*discount[:,1]

    bl_dalys=bl_yll+psens_bl["yld"]
    s1_dalys=s1_yll+psens_s1["yld"]

    #Productivity
    bl_prod=np.sum(ghosts_cbl[:,0:4], axis=1)*discount[:,1]*etp_ratio*gdp_grw[:,2]
    s1_prod=np.sum(ghosts_cs1[:,0:4], axis=1)*discount[:,1]*etp_ratio*gdp_grw[:,2]

    ## NEB
    bl_tc_ins=bl_prod[:]+dmc_bl[:,0]+imc_bl[:,0]
    s1_tc_ins=s1_prod[:]+dmc_s1[:,0]+imc_s1[:,0]

    bl_cum_tc=np.zeros(np.shape(bl_tc_ins))
    s1_cum_tc=np.zeros(np.shape(s1_tc_ins))

    for i in range(len(bl_tc_ins)):
        if i < 1:
            bl_cum_tc[i]=bl_tc_ins[i]
            s1_cum_tc[i]=s1_tc_ins[i]
        else:
            bl_cum_tc[i]=bl_cum_tc[i-1]+bl_tc_ins[i]
            s1_cum_tc[i]=s1_cum_tc[i-1]+s1_tc_ins[i]

    s1_neb=bl_cum_tc-s1_cum_tc

    ## ICER
    bl_dc_ins=dmc_bl[:,0]+imc_bl[:,0]
    s1_dc_ins=dmc_s1[:,0]+imc_s1[:,0]

    bl_cum_dc=np.zeros(np.shape(bl_dc_ins))
    s1_cum_dc=np.zeros(np.shape(s1_dc_ins))
    bl_cum_daly=np.zeros(np.shape(bl_dc_ins))
    s1_cum_daly=np.zeros(np.shape(s1_dc_ins))

    for i in range(len(bl_dc_ins)):
        if i <1:
            bl_cum_dc[i]=bl_dc_ins[i]
            bl_cum_daly[i]=bl_dalys[i]

            s1_cum_dc[i]=s1_dc_ins[i]
            s1_cum_daly[i]=s1_dalys[i]
        else:
            bl_cum_dc[i]=bl_dc_ins[i]+bl_cum_dc[i-1]
            bl_cum_daly[i]=bl_dalys[i]+bl_cum_daly[i-1]

            s1_cum_dc[i]=s1_dc_ins[i]+s1_cum_dc[i-1]
            s1_cum_daly[i]=s1_dalys[i]+s1_cum_daly[i-1]

    s1_icer=-(bl_cum_dc-s1_cum_dc)/(bl_cum_daly-s1_cum_daly)
    s1_icer=np.nan_to_num(s1_icer,0)

    ## Return results needed for tabulation (within the input dictionaries for simplicity)


    psens_s1["cdaly_bl"]=bl_cum_daly
    psens_s1["cdaly_s1"]=s1_cum_daly

    psens_s1["cdirc_bl"]=bl_cum_dc
    psens_s1["cdirc_s1"]=s1_cum_dc

    psens_s1["prod"]=s1_prod
    psens_s1["icer"]=s1_icer
    psens_s1["neb"]=s1_neb

    return psens_s1

def econsens_econ(bl_cent, cent_s1, econ, reg, cost_data, disc_rate):
    import numpy as np
    import pandas as pd

    #Discounting Array
    if disc_rate >1:
        disc_rate=disc_rate/100
    else:
        disc_rate=disc_rate

    discount=np.zeros((len(np.arange(1990,2099,1)),3))
    discount[:,0]=np.arange(1990,2099,1)

    for idx,val in enumerate(discount[:,0]):
        if val <= 2021:
            discount[idx,1]=1
            discount[idx,2]=0 #use this for high risk screening as a check
        else:
            discount[idx,1:]=(1-disc_rate)**(val-2022)

    vax_costs=pd.read_excel(cost_data, sheet_name="vax")
    care_costs=pd.read_excel(cost_data, sheet_name="care")


    """Vaccine Costs (just add/substract costs from relevant parts)"""
    #Current
    bd_pe=vax_costs[reg].iloc[0]
    hb3_pe=vax_costs[reg].iloc[1]

    bl_bd_pe, s1_bd_pe, bl_hb3_pe, s1_hb3_pe=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))

    bl_bd_pe[:,0]=bl_cent["bd_cov"][:]*bl_cent["births"][:]*bd_pe*discount[:,1]
    bl_hb3_pe[:,0]=bl_cent["hb3_cov"][:]*bl_cent["births"][:]*hb3_pe*discount[:,1]
    s1_bd_pe[:,0]=cent_s1["bd_cov"][:]*cent_s1["births"][:]*bd_pe*discount[:,1]
    s1_hb3_pe[:,0]=cent_s1["hb3_cov"][:]*cent_s1["births"][:]*hb3_pe*discount[:,1]

    #Lower Bound
    bd_lb=vax_costs[reg].iloc[6]
    hb3_lb=vax_costs[reg].iloc[8]

    bl_bd_lb, s1_bd_lb, bl_hb3_lb, s1_hb3_lb=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))

    bl_bd_lb[:,0]=bl_cent["bd_cov"][:]*bl_cent["births"][:]*bd_lb*discount[:,1]
    bl_hb3_lb[:,0]=bl_cent["hb3_cov"][:]*bl_cent["births"][:]*hb3_lb*discount[:,1]
    s1_bd_lb[:,0]=cent_s1["bd_cov"][:]*cent_s1["births"][:]*bd_lb*discount[:,1]
    s1_hb3_lb[:,0]=cent_s1["hb3_cov"][:]*cent_s1["births"][:]*hb3_lb*discount[:,1]

    #Upper Bound
    bd_ub=vax_costs[reg].iloc[7]
    hb3_ub=vax_costs[reg].iloc[9]

    bl_bd_ub, s1_bd_ub, bl_hb3_ub, s1_hb3_ub=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))

    bl_bd_ub[:,0]=bl_cent["bd_cov"][:]*bl_cent["births"][:]*bd_ub*discount[:,1]
    bl_hb3_ub[:,0]=bl_cent["hb3_cov"][:]*bl_cent["births"][:]*hb3_ub*discount[:,1]
    s1_bd_ub[:,0]=cent_s1["bd_cov"][:]*cent_s1["births"][:]*bd_ub*discount[:,1]
    s1_hb3_ub[:,0]=cent_s1["hb3_cov"][:]*cent_s1["births"][:]*hb3_ub*discount[:,1]

    # Difference (to be added to direct costs)

    bl_bd_lbd=bl_bd_lb-bl_bd_pe
    bl_hb3_lbd=bl_hb3_lb-bl_hb3_pe
    s1_bd_lbd=s1_bd_lb-s1_bd_pe
    s1_hb3_lbd=s1_hb3_lb-s1_hb3_pe

    bl_bd_ubd=bl_bd_ub-bl_bd_pe
    bl_hb3_ubd=bl_hb3_ub-bl_hb3_pe
    s1_bd_ubd=s1_bd_ub-s1_bd_pe
    s1_hb3_ubd=s1_hb3_ub-s1_hb3_pe

    #total direct costs for tabulation, NEB, and ICER
    vcost_lb_bl=econ["cbl_dirc"]+(bl_bd_lbd[:,0]+bl_hb3_lbd[:,0])
    vcost_ub_bl=econ["cbl_dirc"]+(bl_bd_ubd[:,0]+bl_hb3_ubd[:,0])

    vcost_lb_s1=econ["cs1_dirc"]+(s1_bd_lbd[:,0]+s1_hb3_lbd[:,0])
    vcost_ub_s1=econ["cs1_dirc"]+(s1_bd_ubd[:,0]+s1_hb3_ubd[:,0])

    """Treatment Costs (add mAVs)"""
    #current
    tx_pe=care_costs[reg].iloc[3]
    mav_cost=vax_costs[reg].iloc[11]

    tx_bl_pe=bl_cent["tx_cov"]*tx_pe*discount[:,1]
    mav_bl_s1=cent_s1["mav_n"][:]*mav_cost*discount[:,1]
    tx_s1_pe=cent_s1["tx_cov"]*tx_pe*discount[:,1]

    #halved
    tx_bl_lb=tx_bl_pe*0.5
    mav_s1_lb=mav_bl_s1*0.5
    tx_s1_lb=tx_s1_pe*0.5

    #doubled
    tx_bl_ub=tx_bl_pe*2
    mav_s1_ub=mav_bl_s1*2
    tx_s1_ub=tx_s1_pe*2

    # Difference (to be added to direct costs)
    tx_bl_lbd=tx_bl_lb-tx_bl_pe
    tx_bl_ubd=tx_bl_ub-tx_bl_pe

    tx_s1_lbd=(tx_s1_lb+mav_s1_lb)-(tx_s1_pe+mav_bl_s1)
    tx_s1_ubd=(tx_s1_ub+mav_s1_ub)-(tx_s1_pe+mav_bl_s1)

    #total direct costs for tabulation, NEB, and ICER
    tcost_lb_bl=econ["cbl_dirc"]+tx_bl_lbd
    tcost_ub_bl=econ["cbl_dirc"]+tx_bl_ubd

    tcost_lb_s1=econ["cs1_dirc"]+tx_s1_lbd
    tcost_ub_s1=econ["cs1_dirc"]+tx_s1_ubd

    """mAV with HBIG"""
    hrs_cost=vax_costs[reg].iloc[10]
    dx_cost=care_costs[reg].iloc[0]
    mav_hbig=vax_costs[reg].iloc[12]

    mv_bl_pe, mv_s1_pe, mv_bl_hbig, mv_s1_hbig=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))
    #current (no HBIG)
    mv_bl_pe[:,0]=(bl_cent["prg_scr"][:]*dx_cost*discount[:,1])#(cent["mav_n"][:]*mav_cost*discount[:,1])+(cent["prg_hrs"][:]*hrs_cost*discount[:,1])
    mv_s1_pe[:,0]=(cent_s1["mav_n"][:]*mav_cost*discount[:,1])+(cent_s1["prg_scr"][:]*dx_cost*discount[:,1])+(cent_s1["prg_hrs"][:]*hrs_cost*discount[:,2])
    # upper bound (with HBIG)
    mv_bl_hbig[:,0]=(bl_cent["prg_scr"][:]*dx_cost*discount[:,1])#(cent["mav_n"][:]*mav_cost*discount[:,1])+(cent["prg_hrs"][:]*hrs_cost*discount[:,1])
    mv_s1_hbig[:,0]=(cent_s1["mav_n"][:]*mav_hbig*discount[:,1])+(cent_s1["prg_scr"][:]*dx_cost*discount[:,1])+(cent_s1["prg_hrs"][:]*hrs_cost*discount[:,2])
    # Difference (to be added to direct costs)
    bl_mvd=mv_bl_hbig-mv_bl_pe
    s1_mvd=mv_s1_hbig-mv_s1_pe

    #total direct costs for tabulation, NEB, and ICER
    hbig_bl=econ["cbl_dirc"]+bl_mvd[:,0]
    hbig_s1=econ["cs1_dirc"]+s1_mvd[:,0]

    """HCC surveillance costs"""
    #current
    hcc_cost=care_costs[reg].iloc[6] #HCC surveillance costs
    hcc_prp=care_costs[reg].iloc[8] #proportion of CHB in 15-49y who are 40-49y (HCC screening)
    util=0.25

    hccs_bl_pe, hccs_s1_pe, hccs_bl_lb, hccs_s1_lb, hccs_bl_ub, hccs_s1_ub=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))

    hccs_bl_pe[:, 0]=((bl_cent["pop_hcc"][:]*hcc_cost*util)+(bl_cent["tgt_hcc"][:]*hcc_cost*util)+(bl_cent["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util))*discount[:,1]
    hccs_s1_pe[:, 0]=((cent_s1["pop_hcc"][:]*hcc_cost*util)+(cent_s1["tgt_hcc"][:]*hcc_cost*util)+(cent_s1["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util))*discount[:,1]

    #halved
    hccs_bl_lb=hccs_bl_pe*0.5
    hccs_s1_lb=hccs_s1_pe*0.5

    #doubled
    hccs_bl_ub=hccs_bl_pe*2
    hccs_s1_ub=hccs_s1_pe*2

    # Difference (to be added to indirect costs)
    bl_hccs_lbd=hccs_bl_lb-hccs_bl_pe
    bl_hccs_ubd=hccs_bl_ub-hccs_bl_pe

    s1_hccs_lbd=hccs_s1_lb-hccs_s1_pe
    s1_hccs_ubd=hccs_s1_ub-hccs_s1_pe

    #total direct costs for tabulation, NEB, and ICER
    hccs_lb_bl=econ["cbl_dirc"]+bl_hccs_lbd[:,0]
    hccs_ub_bl=econ["cbl_dirc"]+bl_hccs_ubd[:,0]

    hccs_lb_s1=econ["cs1_dirc"]+s1_hccs_lbd[:,0]
    hccs_ub_s1=econ["cs1_dirc"]+s1_hccs_ubd[:,0]

    """Hospitalization costs"""
    #current
    hosp_cost=care_costs[reg].iloc[7]
    tx_hosp=0.5

    hospc_bl_pe, hospc_s1_pe, hospc_bl_lb, hospc_s1_lb, hospc_bl_ub, hospc_s1_ub=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))

    hospc_bl_pe[:, 0]=((bl_cent["hsp_utx"][:]*hosp_cost*util)+(bl_cent["hsp_tx"][:]*hosp_cost*util*tx_hosp))*discount[:,1]
    hospc_s1_pe[:, 0]=((cent_s1["hsp_utx"][:]*hosp_cost*util)+(cent_s1["hsp_tx"][:]*hosp_cost*util*tx_hosp))*discount[:,1]

    #halved
    hospc_bl_lb=hospc_bl_pe*0.5
    hospc_s1_lb=hospc_s1_pe*0.5

    #doubled
    hospc_bl_ub=hospc_bl_pe*2
    hospc_s1_ub=hospc_s1_pe*2

    # Difference (to be added to indirect costs)
    bl_hospc_lbd=hospc_bl_lb-hospc_bl_pe
    bl_hospc_ubd=hospc_bl_ub-hospc_bl_pe

    s1_hospc_lbd=hospc_s1_lb-hospc_s1_pe
    s1_hospc_ubd=hospc_s1_ub-hospc_s1_pe

    #total direct costs for tabulation, NEB, and ICER
    hospc_lb_bl=econ["cbl_dirc"]+bl_hospc_lbd[:,0]
    hospc_ub_bl=econ["cbl_dirc"]+bl_hospc_ubd[:,0]

    hospc_lb_s1=econ["cs1_dirc"]+s1_hospc_lbd[:,0]
    hospc_ub_s1=econ["cs1_dirc"]+s1_hospc_ubd[:,0]

    """UHC assumptions (vaccination, diagnosis and high-risk screening)"""
    bd_lb_hr=vax_costs[reg].iloc[2]
    hb3_lb_hr=vax_costs[reg].iloc[3]
    dx_lb_hr=care_costs[reg].iloc[1]
    hrs_lb_hr=vax_costs[reg].iloc[13]

    bd_ub_hr=vax_costs[reg].iloc[4]
    hb3_ub_hr=vax_costs[reg].iloc[5]
    dx_ub_hr=care_costs[reg].iloc[2]
    hrs_ub_hr=vax_costs[reg].iloc[14]

    #current
    bl_bd_pe, s1_bd_pe, bl_hb3_pe, s1_hb3_pe=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))

    bl_bd_pe[:,0]=bl_cent["bd_cov"][:]*bl_cent["births"][:]*bd_pe*discount[:,1]
    bl_hb3_pe[:,0]=bl_cent["hb3_cov"][:]*bl_cent["births"][:]*hb3_pe*discount[:,1]
    s1_bd_pe[:,0]=cent_s1["bd_cov"][:]*cent_s1["births"][:]*bd_pe*discount[:,1]
    s1_hb3_pe[:,0]=cent_s1["hb3_cov"][:]*cent_s1["births"][:]*hb3_pe*discount[:,1]

    mv_bl_pe, mv_s1_pe=np.zeros((109,1)),np.zeros((109,1))

    mv_bl_pe[:,0]=(bl_cent["prg_scr"][:]*dx_cost*discount[:,1])#(cent["mav_n"][:]*mav_cost*discount[:,1])+(cent["prg_hrs"][:]*hrs_cost*discount[:,1])
    mv_s1_pe[:,0]=(cent_s1["mav_n"][:]*mav_cost*discount[:,1])+(cent_s1["prg_scr"][:]*dx_cost*discount[:,1])+(cent_s1["prg_hrs"][:]*hrs_cost*discount[:,2])

    bl_dx_inc, s1_dx_inc=np.zeros((109,1)),np.zeros((109,1))

    for i in range(len(bl_cent["dx_prop"])):
        if i < 1:
            bl_dx_inc[i,0]=bl_cent["dx_prop"][i]
            s1_dx_inc[i,0]=cent_s1["dx_prop"][i]
        else:
            bl_dx_inc[i,0]=max(bl_cent["dx_prop"][i]-bl_cent["dx_prop"][i-1],0)
            s1_dx_inc[i,0]=max(cent_s1["dx_prop"][i]-cent_s1["dx_prop"][i-1],0)

    bl_dxc_pe, s1_dxc_pe=np.zeros((109,1)),np.zeros((109,1))

    for yr in range(len(bl_dx_inc)):
        bl_dxc_pe[yr,0]=(bl_cent["dx_rate"][yr]*dx_cost*discount[yr,1])+(dx_cost*bl_dx_inc[yr,0]*(bl_cent["pop"][yr]-bl_cent["dx_rate"][yr])*discount[yr,1])
        s1_dxc_pe[yr,0]=(cent_s1["dx_rate"][yr]*dx_cost*discount[yr,1])+(dx_cost*s1_dx_inc[yr,0]*(cent_s1["pop"][yr]-cent_s1["dx_rate"][yr])*discount[yr,1])

    #0%
    bl_bd_ulb, s1_bd_ulb, bl_hb3_ulb, s1_hb3_ulb=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))

    bl_bd_ulb[:,0]=bl_cent["bd_cov"][:]*bl_cent["births"][:]*bd_lb_hr*discount[:,1]
    bl_hb3_ulb[:,0]=bl_cent["hb3_cov"][:]*bl_cent["births"][:]*hb3_lb_hr*discount[:,1]
    s1_bd_ulb[:,0]=cent_s1["bd_cov"][:]*cent_s1["births"][:]*bd_lb_hr*discount[:,1]
    s1_hb3_ulb[:,0]=cent_s1["hb3_cov"][:]*cent_s1["births"][:]*hb3_lb_hr*discount[:,1]

    mv_bl_ulb, mv_s1_ulb=np.zeros((109,1)),np.zeros((109,1))

    mv_bl_ulb[:,0]=(bl_cent["prg_scr"][:]*dx_lb_hr*discount[:,1])#(cent["mav_n"][:]*mav_cost*discount[:,1])+(cent["prg_hrs"][:]*hrs_cost*discount[:,1])
    mv_s1_ulb[:,0]=(cent_s1["mav_n"][:]*mav_cost*discount[:,1])+(cent_s1["prg_scr"][:]*dx_lb_hr*discount[:,1])+(cent_s1["prg_hrs"][:]*hrs_lb_hr*discount[:,2])

    bl_dxc_ulb, s1_dxc_ulb=np.zeros((109,1)),np.zeros((109,1))

    for yr in range(len(bl_dx_inc)):
        bl_dxc_ulb[yr,0]=(bl_cent["dx_rate"][yr]*dx_lb_hr*discount[yr,1])+(dx_lb_hr*bl_dx_inc[yr,0]*(bl_cent["pop"][yr]-bl_cent["dx_rate"][yr])*discount[yr,1])
        s1_dxc_ulb[yr,0]=(cent_s1["dx_rate"][yr]*dx_lb_hr*discount[yr,1])+(dx_lb_hr*s1_dx_inc[yr,0]*(cent_s1["pop"][yr]-cent_s1["dx_rate"][yr])*discount[yr,1])

    #100%
    bl_bd_uub, s1_bd_uub, bl_hb3_uub, s1_hb3_uub=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))

    bl_bd_uub[:,0]=bl_cent["bd_cov"][:]*bl_cent["births"][:]*bd_ub_hr*discount[:,1]
    bl_hb3_uub[:,0]=bl_cent["hb3_cov"][:]*bl_cent["births"][:]*hb3_ub_hr*discount[:,1]
    s1_bd_uub[:,0]=cent_s1["bd_cov"][:]*cent_s1["births"][:]*bd_ub_hr*discount[:,1]
    s1_hb3_uub[:,0]=cent_s1["hb3_cov"][:]*cent_s1["births"][:]*hb3_ub_hr*discount[:,1]

    mv_bl_uub, mv_s1_uub=np.zeros((109,1)),np.zeros((109,1))

    mv_bl_uub[:,0]=(bl_cent["prg_scr"][:]*dx_ub_hr*discount[:,1])#(cent["mav_n"][:]*mav_cost*discount[:,1])+(cent["prg_hrs"][:]*hrs_cost*discount[:,1])
    mv_s1_uub[:,0]=(cent_s1["mav_n"][:]*mav_cost*discount[:,1])+(cent_s1["prg_scr"][:]*dx_ub_hr*discount[:,1])+(cent_s1["prg_hrs"][:]*hrs_ub_hr*discount[:,2])

    bl_dxc_uub, s1_dxc_uub=np.zeros((109,1)),np.zeros((109,1))

    for yr in range(len(bl_dx_inc)):
        bl_dxc_uub[yr,0]=(bl_cent["dx_rate"][yr]*dx_ub_hr*discount[yr,1])+(dx_ub_hr*bl_dx_inc[yr,0]*(bl_cent["pop"][yr]-bl_cent["dx_rate"][yr])*discount[yr,1])
        s1_dxc_uub[yr,0]=(cent_s1["dx_rate"][yr]*dx_ub_hr*discount[yr,1])+(dx_ub_hr*s1_dx_inc[yr,0]*(cent_s1["pop"][yr]-cent_s1["dx_rate"][yr])*discount[yr,1])

    # Difference (to be added to direct (dx, mav, hepB-bd) costs)
    bl_lb_uhcd=(bl_bd_ulb+bl_hb3_ulb+mv_bl_ulb+bl_dxc_ulb)-(bl_bd_pe+bl_hb3_pe+mv_bl_pe+bl_dxc_pe)
    bl_ub_uhcd=(bl_bd_uub+bl_hb3_uub+mv_bl_uub+bl_dxc_uub)-(bl_bd_pe+bl_hb3_pe+mv_bl_pe+bl_dxc_pe)

    s1_lb_uhcd=(s1_bd_ulb+s1_hb3_ulb+mv_s1_ulb+s1_dxc_ulb)-(s1_bd_pe+s1_hb3_pe+mv_s1_pe+s1_dxc_pe)
    s1_ub_uhcd=(s1_bd_uub+s1_hb3_uub+mv_s1_uub+s1_dxc_uub)-(s1_bd_pe+s1_hb3_pe+mv_s1_pe+s1_dxc_pe)

    #total direct costs for tabulation, NEB, and ICER
    uhc_lb_bl=econ["cbl_dirc"]+bl_lb_uhcd[:,0]
    uhc_ub_bl=econ["cbl_dirc"]+bl_ub_uhcd[:,0]

    uhc_lb_s1=econ["cs1_dirc"]+s1_lb_uhcd[:,0]
    uhc_ub_s1=econ["cs1_dirc"]+s1_ub_uhcd[:,0]

    """Utilization assumptions (HCC surveillance and hospitalisation)"""
    util_lb=0
    util_ub=0.5

    dx_dmc=care_costs[reg].iloc[4] #cost of disease management (diagnosed)
    tx_dmc=care_costs[reg].iloc[4] #cost of disease management (on treatment)

    #current

    manu_bl_pe, manu_s1_pe=np.zeros((109,1)), np.zeros((109,1))

    manu_bl_pe[:,0]=((bl_cent["dm_dx"][:]*dx_dmc*util)+(bl_cent["dm_tx"][:]*tx_dmc*util))*discount[:,1]
    manu_s1_pe[:,0]=((cent_s1["dm_dx"][:]*dx_dmc*util)+(cent_s1["dm_tx"][:]*tx_dmc*util))*discount[:,1]


    hospu_bl_pe, hospu_s1_pe=np.zeros((109,1)),np.zeros((109,1))

    hospu_bl_pe[:, 0]=((bl_cent["hsp_utx"][:]*hosp_cost*util)+(bl_cent["hsp_tx"][:]*hosp_cost*util*tx_hosp))*discount[:,1]
    hospu_s1_pe[:, 0]=((cent_s1["hsp_utx"][:]*hosp_cost*util)+(cent_s1["hsp_tx"][:]*hosp_cost*util*tx_hosp))*discount[:,1]

    hccu_bl_pe, hccu_s1_pe=np.zeros((109,1)),np.zeros((109,1))

    hccu_bl_pe[:, 0]=((bl_cent["pop_hcc"][:]*hcc_cost*util)+(bl_cent["tgt_hcc"][:]*hcc_cost*util)+(bl_cent["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util))*discount[:,1]
    hccu_s1_pe[:, 0]=((cent_s1["pop_hcc"][:]*hcc_cost*util)+(cent_s1["tgt_hcc"][:]*hcc_cost*util)+(cent_s1["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util))*discount[:,1]

    #0%
    manu_bl_lb, manu_s1_lb=np.zeros((109,1)), np.zeros((109,1))

    manu_bl_lb[:,0]=((bl_cent["dm_dx"][:]*dx_dmc*util_lb)+(bl_cent["dm_tx"][:]*tx_dmc*util_lb))*discount[:,1]
    manu_s1_lb[:,0]=((cent_s1["dm_dx"][:]*dx_dmc*util_lb)+(cent_s1["dm_tx"][:]*tx_dmc*util_lb))*discount[:,1]

    hospu_bl_lb, hospu_s1_lb=np.zeros((109,1)),np.zeros((109,1))

    hospu_bl_lb[:, 0]=((bl_cent["hsp_utx"][:]*hosp_cost*util_lb)+(bl_cent["hsp_tx"][:]*hosp_cost*util_lb*tx_hosp))*discount[:,1]
    hospu_s1_lb[:, 0]=((cent_s1["hsp_utx"][:]*hosp_cost*util_lb)+(cent_s1["hsp_tx"][:]*hosp_cost*util_lb*tx_hosp))*discount[:,1]

    hccu_bl_lb, hccu_s1_lb=np.zeros((109,1)),np.zeros((109,1))

    hccu_bl_lb[:, 0]=((bl_cent["pop_hcc"][:]*hcc_cost*util_lb)+(bl_cent["tgt_hcc"][:]*hcc_cost*util_lb)+(bl_cent["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util_lb))*discount[:,1]
    hccu_s1_lb[:, 0]=((cent_s1["pop_hcc"][:]*hcc_cost*util_lb)+(cent_s1["tgt_hcc"][:]*hcc_cost*util_lb)+(cent_s1["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util_lb))*discount[:,1]

    #50%
    manu_bl_ub, manu_s1_ub=np.zeros((109,1)), np.zeros((109,1))

    manu_bl_ub[:,0]=((bl_cent["dm_dx"][:]*dx_dmc*util_ub)+(bl_cent["dm_tx"][:]*tx_dmc*util_ub))*discount[:,1]
    manu_s1_ub[:,0]=((cent_s1["dm_dx"][:]*dx_dmc*util_ub)+(cent_s1["dm_tx"][:]*tx_dmc*util_ub))*discount[:,1]


    hospu_bl_ub, hospu_s1_ub=np.zeros((109,1)),np.zeros((109,1))

    hospu_bl_ub[:, 0]=((bl_cent["hsp_utx"][:]*hosp_cost*util_ub)+(bl_cent["hsp_tx"][:]*hosp_cost*util_ub*tx_hosp))*discount[:,1]
    hospu_s1_ub[:, 0]=((cent_s1["hsp_utx"][:]*hosp_cost*util_ub)+(cent_s1["hsp_tx"][:]*hosp_cost*util_ub*tx_hosp))*discount[:,1]

    hccu_bl_ub, hccu_s1_ub=np.zeros((109,1)),np.zeros((109,1))

    hccu_bl_ub[:, 0]=((bl_cent["pop_hcc"][:]*hcc_cost*util_ub)+(bl_cent["tgt_hcc"][:]*hcc_cost*util_ub)+(bl_cent["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util_ub))*discount[:,1]
    hccu_s1_ub[:, 0]=((cent_s1["pop_hcc"][:]*hcc_cost*util_ub)+(cent_s1["tgt_hcc"][:]*hcc_cost*util_ub)+(cent_s1["tgt_hcc_b"][:]*hcc_prp*hcc_cost*util_ub))*discount[:,1]
    # Difference (to be added to indirect (hcc surveillance, dis man, hosp) costs)
    bl_util_lbd=(hospu_bl_lb+hccu_bl_lb+manu_bl_lb)-(hospu_bl_pe+hccu_bl_pe+manu_bl_pe)
    bl_util_ubd=(hospu_bl_ub+hccu_bl_ub+manu_bl_ub)-(hospu_bl_pe+hccu_bl_pe+manu_bl_pe)

    s1_util_lbd=(hospu_s1_lb+hccu_s1_lb+manu_s1_lb)-(hospu_s1_pe+hccu_s1_pe+manu_bl_pe)
    s1_util_ubd=(hospu_s1_ub+hccu_s1_ub+manu_s1_ub)-(hospu_s1_pe+hccu_s1_pe+manu_bl_pe)

    #total direct costs for tabulation, NEB, and ICER
    util_lb_bl=econ["cbl_dirc"]+bl_util_lbd[:,0]
    util_ub_bl=econ["cbl_dirc"]+bl_util_ubd[:,0]

    util_lb_s1=econ["cs1_dirc"]+s1_util_lbd[:,0]
    util_ub_s1=econ["cs1_dirc"]+s1_util_ubd[:,0]

    """Hospitalisation (tx) assumption"""
    tx_hosp_lb=1
    tx_hosp_ub=0.25

    hospt_bl_pe, hospt_s1_pe, hospt_bl_lb, hospt_s1_lb, hospt_bl_ub, hospt_s1_ub=np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1)),np.zeros((109,1))

    #current
    hospt_bl_pe[:, 0]=((bl_cent["hsp_utx"][:]*hosp_cost*util)+(bl_cent["hsp_tx"][:]*hosp_cost*util*tx_hosp))*discount[:,1]
    hospt_s1_pe[:, 0]=((cent_s1["hsp_utx"][:]*hosp_cost*util)+(cent_s1["hsp_tx"][:]*hosp_cost*util*tx_hosp))*discount[:,1]
    #0%
    hospt_bl_lb[:, 0]=((bl_cent["hsp_utx"][:]*hosp_cost*util)+(bl_cent["hsp_tx"][:]*hosp_cost*util*tx_hosp_lb))*discount[:,1]
    hospt_s1_lb[:, 0]=((cent_s1["hsp_utx"][:]*hosp_cost*util)+(cent_s1["hsp_tx"][:]*hosp_cost*util*tx_hosp_lb))*discount[:,1]
    #75%
    hospt_bl_ub[:, 0]=((bl_cent["hsp_utx"][:]*hosp_cost*util)+(bl_cent["hsp_tx"][:]*hosp_cost*util*tx_hosp_ub))*discount[:,1]
    hospt_s1_ub[:, 0]=((cent_s1["hsp_utx"][:]*hosp_cost*util)+(cent_s1["hsp_tx"][:]*hosp_cost*util*tx_hosp_ub))*discount[:,1]

    # Difference (to be added to indirect costs)
    bl_hospt_lbd=hospt_bl_lb-hospt_bl_pe
    bl_hospt_ubd=hospt_bl_ub-hospt_bl_pe

    s1_hospt_lbd=hospt_s1_lb-hospt_s1_pe
    s1_hospt_ubd=hospt_s1_ub-hospt_s1_pe

    #total direct costs for tabulation, NEB, and ICER
    hospt_lb_bl=econ["cbl_dirc"]+bl_hospt_lbd[:,0]
    hospt_ub_bl=econ["cbl_dirc"]+bl_hospt_ubd[:,0]

    hospt_lb_s1=econ["cs1_dirc"]+s1_hospt_lbd[:,0]
    hospt_ub_s1=econ["cs1_dirc"]+s1_hospt_ubd[:,0]

    res_dict={}
    res_dict["vax_lb_bl"]=vcost_lb_bl
    res_dict["vax_ub_bl"]=vcost_ub_bl
    res_dict["vax_lb_s1"]=vcost_lb_s1
    res_dict["vax_ub_s1"]=vcost_ub_s1
    res_dict["trt_lb_bl"]=tcost_lb_bl
    res_dict["trt_ub_bl"]=tcost_ub_bl
    res_dict["trt_lb_s1"]=tcost_lb_s1
    res_dict["trt_ub_s1"]=tcost_ub_s1
    res_dict["hbig_bl"]=hbig_bl
    res_dict["hbig_s1"]=hbig_s1
    res_dict["surv_lb_bl"]=hccs_lb_bl
    res_dict["surv_ub_bl"]=hccs_ub_bl
    res_dict["surv_lb_s1"]=hccs_lb_s1
    res_dict["surv_ub_s1"]=hccs_ub_s1
    res_dict["hosp_lb_bl"]=hospc_lb_bl
    res_dict["hosp_ub_bl"]=hospc_ub_bl
    res_dict["hosp_lb_s1"]=hospc_lb_s1
    res_dict["hosp_ub_s1"]=hospc_ub_s1
    res_dict["uhc_lb_bl"]=uhc_lb_bl
    res_dict["uhc_ub_bl"]=uhc_ub_bl
    res_dict["uhc_lb_s1"]=uhc_lb_s1
    res_dict["uhc_ub_s1"]=uhc_ub_s1
    res_dict["util_lb_bl"]=util_lb_bl
    res_dict["util_ub_bl"]=util_ub_bl
    res_dict["util_lb_s1"]=util_lb_s1
    res_dict["util_ub_s1"]=util_ub_s1
    res_dict["htx_lb_bl"]=hospt_lb_bl
    res_dict["htx_ub_bl"]=hospt_ub_bl
    res_dict["htx_lb_s1"]=hospt_lb_s1
    res_dict["htx_ub_s1"]=hospt_ub_s1

    return res_dict

def sens_nebicer(econ_sens, econ_main):

    import numpy as np
    econ_outs={}

    ## Vaccine Costs
    econ_outs["vax_lb_neb"]=(np.cumsum(econ_sens["vax_lb_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["vax_lb_s1"])+np.cumsum(econ_main["cs1_prod"]))
    econ_outs["vax_ub_neb"]=(np.cumsum(econ_sens["vax_ub_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["vax_ub_s1"])+np.cumsum(econ_main["cs1_prod"]))

    econ_outs["vax_lb_icer"]=-(np.cumsum(econ_sens["vax_lb_bl"])-np.cumsum(econ_sens["vax_lb_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))
    econ_outs["vax_ub_icer"]=-(np.cumsum(econ_sens["vax_ub_bl"])-np.cumsum(econ_sens["vax_ub_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))

    ## Treatment Costs
    econ_outs["tc_lb_neb"]=(np.cumsum(econ_sens["trt_lb_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["trt_lb_s1"])+np.cumsum(econ_main["cs1_prod"]))
    econ_outs["tc_ub_neb"]=(np.cumsum(econ_sens["trt_ub_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["trt_ub_s1"])+np.cumsum(econ_main["cs1_prod"]))

    econ_outs["tc_lb_icer"]=-(np.cumsum(econ_sens["trt_lb_bl"])-np.cumsum(econ_sens["trt_lb_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))
    econ_outs["tc_ub_icer"]=-(np.cumsum(econ_sens["trt_ub_bl"])-np.cumsum(econ_sens["trt_ub_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))

    ##HBIG
    econ_outs["hbig_neb"]=(np.cumsum(econ_sens["hbig_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["hbig_s1"])+np.cumsum(econ_main["cs1_prod"]))

    econ_outs["hbig_icer"]=-(np.cumsum(econ_sens["hbig_bl"])-np.cumsum(econ_sens["hbig_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))

    #HCC Surveillance
    econ_outs["surv_lb_neb"]=(np.cumsum(econ_sens["surv_lb_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["surv_lb_s1"])+np.cumsum(econ_main["cs1_prod"]))
    econ_outs["surv_ub_neb"]=(np.cumsum(econ_sens["surv_ub_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["surv_ub_s1"])+np.cumsum(econ_main["cs1_prod"]))

    econ_outs["surv_lb_icer"]=-(np.cumsum(econ_sens["surv_lb_bl"])-np.cumsum(econ_sens["surv_lb_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))
    econ_outs["surv_ub_icer"]=-(np.cumsum(econ_sens["surv_ub_bl"])-np.cumsum(econ_sens["surv_ub_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))

    #Hospital Costs
    econ_outs["hosp_lb_neb"]=(np.cumsum(econ_sens["hosp_lb_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["hosp_lb_s1"])+np.cumsum(econ_main["cs1_prod"]))
    econ_outs["hosp_ub_neb"]=(np.cumsum(econ_sens["hosp_ub_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["hosp_ub_s1"])+np.cumsum(econ_main["cs1_prod"]))

    econ_outs["hosp_lb_icer"]=-(np.cumsum(econ_sens["hosp_lb_bl"])-np.cumsum(econ_sens["hosp_lb_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))
    econ_outs["hosp_ub_icer"]=-(np.cumsum(econ_sens["hosp_ub_bl"])-np.cumsum(econ_sens["hosp_ub_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))

    #UHC assumpyion
    econ_outs["uhc_lb_neb"]=(np.cumsum(econ_sens["uhc_lb_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["uhc_lb_s1"])+np.cumsum(econ_main["cs1_prod"]))
    econ_outs["uhc_ub_neb"]=(np.cumsum(econ_sens["uhc_ub_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["uhc_ub_s1"])+np.cumsum(econ_main["cs1_prod"]))

    econ_outs["uhc_lb_icer"]=-(np.cumsum(econ_sens["uhc_lb_bl"])-np.cumsum(econ_sens["uhc_lb_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))
    econ_outs["uhc_ub_icer"]=-(np.cumsum(econ_sens["uhc_ub_bl"])-np.cumsum(econ_sens["uhc_ub_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))

    #Utilization assumption
    econ_outs["util_lb_neb"]=(np.cumsum(econ_sens["util_lb_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["util_lb_s1"])+np.cumsum(econ_main["cs1_prod"]))
    econ_outs["util_ub_neb"]=(np.cumsum(econ_sens["util_ub_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["util_ub_s1"])+np.cumsum(econ_main["cs1_prod"]))

    econ_outs["util_lb_icer"]=-(np.cumsum(econ_sens["util_lb_bl"])-np.cumsum(econ_sens["util_lb_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))
    econ_outs["util_ub_icer"]=-(np.cumsum(econ_sens["util_ub_bl"])-np.cumsum(econ_sens["util_ub_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))

    #Treatment effectiveness (hospitalizations)
    econ_outs["htx_lb_neb"]=(np.cumsum(econ_sens["htx_lb_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["htx_lb_s1"])+np.cumsum(econ_main["cs1_prod"]))
    econ_outs["htx_ub_neb"]=(np.cumsum(econ_sens["htx_ub_bl"])+np.cumsum(econ_main["cbl_prod"]))-(np.cumsum(econ_sens["htx_ub_s1"])+np.cumsum(econ_main["cs1_prod"]))

    econ_outs["htx_lb_icer"]=-(np.cumsum(econ_sens["htx_lb_bl"])-np.cumsum(econ_sens["htx_lb_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))
    econ_outs["htx_ub_icer"]=-(np.cumsum(econ_sens["htx_ub_bl"])-np.cumsum(econ_sens["htx_ub_s1"]))/(np.cumsum(econ_main["cbl_daly"])-np.cumsum(econ_main["cs1_daly"]))

    return econ_outs

def sens_tables(s1_cent, econ, mtctw, mtctb, mortw, mortb, txtw, txtb, esens, esens_outs):
    import pandas as pd
    import numpy as np
    year=np.arange(1990, 2099,1)


    table=pd.DataFrame(columns=["Scenario", "CHB Incidence (millions)", "HCC Incidence (millions)", "HBV attributable deaths (millions)", "DALYs (millions)",
                          "Healthcare Costs (US$ billions)", "Productivity Losses (US$, billions)", "ICER", "Net Economic Benefit (year achieved)"])

    table["Scenario"]=["Base Case", "MTCT worst case", "MTCT best case", "Mortality worst case", "Mortality best case", "Treatment Effectiveness worst case",
                       "Treatment Effectiveness best case", "Lower bound vaccine commodity cost", "Upper bound vaccine commodity cost", "Lower bound TDF cost",
                       "Upper bound TDF cost", "HBIG + maternal antivirals", "Lower bound HCC surveillance cost", "Upper bound HCC surveillance cost", "Lower bound hospitalisation cost",
                       "Upper bound hospitalisation cost", "Lower bound UHC absorption", "Upper bound UHC absorption", "Lower bound care utilization", "Upper bound care utilization",
                       "Lower bound treatment hospital prevention", "Upper bound treatment hospital prevention"]

    ## Base Case
    table.iloc[0,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[0,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[0,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[0,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[0,5]=round(np.sum(econ["cs1_dirc"][32:61]/1e9),2)
    table.iloc[0,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[0,7]=round(econ["cs1_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(econ["cs1_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[0,8]=year[yneb_i]

    ## MTCT w
    table.iloc[1,1]=round(np.sum(mtctw["chb_inc"][32:61])/1e6,2)
    table.iloc[1,2]=round(np.sum(mtctw["hcc_inc"][32:61])/1e6,2)
    table.iloc[1,3]=round(np.sum(mtctw["mort"][32:61])/1e6,2)
    table.iloc[1,4]=round((mtctw["cdaly_s1"][61]-mtctw["cdaly_s1"][32])/1e6,2)
    table.iloc[1,5]=round((mtctw["cdirc_s1"][61]-mtctw["cdirc_s1"][32])/1e9,2)
    table.iloc[1,6]=round(np.sum(mtctw["prod"][32:61]/1e9),2)
    table.iloc[1,7]=round(mtctw["icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(mtctw["neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[1,8]=year[yneb_i]
    ## MTCT b
    table.iloc[2,1]=round(np.sum(mtctb["chb_inc"][32:61])/1e6,2)
    table.iloc[2,2]=round(np.sum(mtctb["hcc_inc"][32:61])/1e6,2)
    table.iloc[2,3]=round(np.sum(mtctb["mort"][32:61])/1e6,2)
    table.iloc[2,4]=round((mtctb["cdaly_s1"][61]-mtctb["cdaly_s1"][32])/1e6,2)
    table.iloc[2,5]=round((mtctb["cdirc_s1"][61]-mtctb["cdirc_s1"][32])/1e9,2)
    table.iloc[2,6]=round(np.sum(mtctb["prod"][32:61]/1e9),2)
    table.iloc[2,7]=round(mtctb["icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(mtctb["neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[2,8]=year[yneb_i]
    ## mort w
    table.iloc[3,1]=round(np.sum(mortw["chb_inc"][32:61])/1e6,2)
    table.iloc[3,2]=round(np.sum(mortw["hcc_inc"][32:61])/1e6,2)
    table.iloc[3,3]=round(np.sum(mortw["mort"][32:61])/1e6,2)
    table.iloc[3,4]=round((mortw["cdaly_s1"][61]-mortw["cdaly_s1"][32])/1e6,2)
    table.iloc[3,5]=round((mortw["cdirc_s1"][61]-mortw["cdirc_s1"][32])/1e9,2)
    table.iloc[3,6]=round(np.sum(mortw["prod"][32:61]/1e9),2)
    table.iloc[3,7]=round(mortw["icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(mortw["neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[3,8]=year[yneb_i]
    ## mort b
    table.iloc[4,1]=round(np.sum(mortb["chb_inc"][32:61])/1e6,2)
    table.iloc[4,2]=round(np.sum(mortb["hcc_inc"][32:61])/1e6,2)
    table.iloc[4,3]=round(np.sum(mortb["mort"][32:61])/1e6,2)
    table.iloc[4,4]=round((mortb["cdaly_s1"][61]-mortb["cdaly_s1"][32])/1e6,2)
    table.iloc[4,5]=round((mortb["cdirc_s1"][61]-mortb["cdirc_s1"][32])/1e9,2)
    table.iloc[4,6]=round(np.sum(mortb["prod"][32:61]/1e9),2)
    table.iloc[4,7]=round(mortb["icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(mortb["neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[4,8]=year[yneb_i]
    ## tx w
    table.iloc[5,1]=round(np.sum(txtw["chb_inc"][32:61])/1e6,2)
    table.iloc[5,2]=round(np.sum(txtw["hcc_inc"][32:61])/1e6,2)
    table.iloc[5,3]=round(np.sum(txtw["mort"][32:61])/1e6,2)
    table.iloc[5,4]=round((txtw["cdaly_s1"][61]-txtw["cdaly_s1"][32])/1e6,2)
    table.iloc[5,5]=round((txtw["cdirc_s1"][61]-txtw["cdirc_s1"][32])/1e9,2)
    table.iloc[5,6]=round(np.sum(txtw["prod"][32:61]/1e9),2)
    table.iloc[5,7]=round(txtw["icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(txtw["neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[5,8]=year[yneb_i]
    ##tx b
    table.iloc[6,1]=round(np.sum(txtb["chb_inc"][32:61])/1e6,2)
    table.iloc[6,2]=round(np.sum(txtb["hcc_inc"][32:61])/1e6,2)
    table.iloc[6,3]=round(np.sum(txtb["mort"][32:61])/1e6,2)
    table.iloc[6,4]=round((txtb["cdaly_s1"][61]-txtb["cdaly_s1"][32])/1e6,2)
    table.iloc[6,5]=round((txtb["cdirc_s1"][61]-txtb["cdirc_s1"][32])/1e9,2)
    table.iloc[6,6]=round(np.sum(txtb["prod"][32:61]/1e9),2)
    table.iloc[6,7]=round(txtb["icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(txtb["neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[6,8]=year[yneb_i]
  ##tx b
    ## vax com cost lb
    table.iloc[7,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[7,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[7,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[7,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[7,5]=round(np.sum(esens["vax_lb_s1"][32:61]/1e9),2)
    table.iloc[7,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[7,7]=round(esens_outs["vax_lb_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["vax_lb_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[7,8]=year[yneb_i]

    ## vax com cost ub
    table.iloc[8,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[8,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[8,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[8,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[8,5]=round(np.sum(esens["vax_ub_s1"][32:61]/1e9),2)
    table.iloc[8,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[8,7]=round(esens_outs["vax_ub_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["vax_ub_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[8,8]=year[yneb_i]
    # tdf cost lb
    table.iloc[9,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[9,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[9,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[9,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[9,5]=round(np.sum(esens["trt_lb_s1"][32:61]/1e9),2)
    table.iloc[9,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[9,7]=round(esens_outs["tc_lb_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["tc_lb_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[9,8]=year[yneb_i]

    # tdf cost ub
    table.iloc[10,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[10,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[10,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[10,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[10,5]=round(np.sum(esens["trt_ub_s1"][32:61]/1e9),2)
    table.iloc[10,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[10,7]=round(esens_outs["tc_ub_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["tc_ub_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[10,8]=year[yneb_i]


    # hbig
    table.iloc[11,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[11,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[11,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[11,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[11,5]=round(np.sum(esens["hbig_s1"][32:61]/1e9),2)
    table.iloc[11,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[11,7]=round(esens_outs["hbig_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["hbig_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[11,8]=year[yneb_i]
    # hcc surveillance lb
    table.iloc[12,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[12,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[12,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[12,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[12,5]=round(np.sum(esens["surv_lb_s1"][32:61]/1e9),2)
    table.iloc[12,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[12,7]=round(esens_outs["surv_lb_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["surv_lb_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[12,8]=year[yneb_i]
    # hcc surveillance ub
    table.iloc[13,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[13,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[13,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[13,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[13,5]=round(np.sum(esens["surv_ub_s1"][32:61]/1e9),2)
    table.iloc[13,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[13,7]=round(esens_outs["surv_ub_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["surv_ub_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[13,8]=year[yneb_i]
    #hosp cost lb
    table.iloc[14,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[14,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[14,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[14,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[14,5]=round(np.sum(esens["hosp_lb_s1"][32:61]/1e9),2)
    table.iloc[14,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[14,7]=round(esens_outs["hosp_lb_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["hosp_lb_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[14,8]=year[yneb_i]
    #hosp cost ub
    table.iloc[15,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[15,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[15,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[15,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[15,5]=round(np.sum(esens["hosp_ub_s1"][32:61]/1e9),2)
    table.iloc[15,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[15,7]=round(esens_outs["hosp_ub_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["hosp_ub_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[15,8]=year[yneb_i]
    #uhc lb
    table.iloc[16,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[16,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[16,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[16,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[16,5]=round(np.sum(esens["uhc_lb_s1"][32:61]/1e9),2)
    table.iloc[16,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[16,7]=round(esens_outs["uhc_lb_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["uhc_lb_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[16,8]=year[yneb_i]
    #uhc ub
    table.iloc[17,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[17,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[17,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[17,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[17,5]=round(np.sum(esens["uhc_ub_s1"][32:61]/1e9),2)
    table.iloc[17,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[17,7]=round(esens_outs["uhc_ub_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["uhc_ub_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[17,8]=year[yneb_i]
    #util lb
    table.iloc[18,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[18,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[18,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[18,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[18,5]=round(np.sum(esens["util_lb_s1"][32:61]/1e9),2)
    table.iloc[18,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[18,7]=round(esens_outs["util_lb_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["util_lb_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[18,8]=year[yneb_i]
    #util ub
    table.iloc[19,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[19,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[19,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[19,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[19,5]=round(np.sum(esens["util_ub_s1"][32:61]/1e9),2)
    table.iloc[19,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[19,7]=round(esens_outs["util_ub_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["util_ub_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[19,8]=year[yneb_i]
    #tx hosp lb
    table.iloc[20,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[20,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[20,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[20,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[20,5]=round(np.sum(esens["htx_lb_s1"][32:61]/1e9),2)
    table.iloc[20,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[20,7]=round(esens_outs["htx_lb_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["htx_lb_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[20,8]=year[yneb_i]
    #tx hosp ub
    table.iloc[21,1]=round(np.sum(s1_cent["chb_inc"][32:61])/1e6,2)
    table.iloc[21,2]=round(np.sum(s1_cent["hcc_inc"][32:61])/1e6,2)
    table.iloc[21,3]=round(np.sum(s1_cent["mort"][32:61])/1e6,2)
    table.iloc[21,4]=round(np.sum(econ["cs1_daly"][32:61]/1e6),2)
    table.iloc[21,5]=round(np.sum(esens["htx_ub_s1"][32:61]/1e9),2)
    table.iloc[21,6]=round(np.sum(econ["cs1_prod"][32:61]/1e9),2)
    table.iloc[21,7]=round(esens_outs["htx_ub_icer"][60],2)
    year_idx=[]
    for idx,val in enumerate(esens_outs["htx_ub_neb"]):
        if val>0:
            year_idx.append(idx)
    yneb_i=min(year_idx)
    table.iloc[21,8]=year[yneb_i]

    return table