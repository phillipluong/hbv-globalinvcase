import os
wd= 'C:/Users/Phil/Documents/GitHub/hbv-globalinvcase'#set own working directory
os.chdir(wd)

import atomica as at
import matplotlib.pyplot as plt
import hbv_functions as hbv

# for region in ['AFR', 'AMR', 'EMR', 'EUR', 'SEAR', 'WPR']:
for region in ['AMR', 'WPR']:
    F = at.ProjectFramework('frameworks/hbv_v14_gamma_mav.xlsx')
    D = at.ProjectData.from_spreadsheet(f"applications/region_{region.lower()}/{region}_db_mav.xlsx", framework=F)
    P= at.Project(framework=F, databook=D, sim_dt=0.25, sim_start=1990, sim_end=2099, do_run=False)
    cal = P.make_parset()
    cal.load_calibration(f"applications/region_{region.lower()}/{region}_calib.xlsx")
    cal = hbv.pop_calib(P,cal)
    cal = hbv.epi_calib(P,cal)
    cal.save_calibration(f'{region}_calib.xlsx')