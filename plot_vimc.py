import os
wd= 'C:/Users/Phil/Documents/GitHub/hbv-globalinvcase'#set own working directory
os.chdir(wd)

import numpy as np
import atomica as at
import pandas as pd
import matplotlib.pyplot as plt

import hbv_functions as hbv

loc = f"valuations/vimc/"
fw = 'hbv_v14_gamma_vimc.xlsx'
db = loc+'vimc_AFR_db_v1_3.xlsx'
cl = loc+'vimc_AFR_calib_v1_2.xlsx'
P = at.Project(framework=fw, databook=db, sim_start=1990, sim_end=2101, sim_dt=0.25, do_run=False)
cal = P.make_parset()
cal.load_calibration(cl)
res = P.run_sim(parset=cal, result_name = 'Test')

d=at.PlotData([res], outputs=["alive"])
at.plot_series(d, data=P.data, axis="pops")

d=at.PlotData([res], outputs=["tot_inc"])
at.plot_series(d, data=P.data, axis="pops")

d=at.PlotData([res], outputs=["dalys"])
at.plot_series(d, data=P.data, axis="pops")

d=at.PlotData([res], outputs=["deaths"]) # Deaths is cl_acu, cl_cir, cl_hcc
at.plot_series(d, data=P.data, axis="pops")

plt.show()

'''
Other plots:
- Total number of cases
- Total Deaths
- Total DALYs
- Compare baseline vs novax
'''

