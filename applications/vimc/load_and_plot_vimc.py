import os
wd= 'C:/Users/Phil/Documents/GitHub/hbv-globalinvcase'#set own working directory
wd_v = wd + '/applications/vimc'
os.chdir(wd)

import numpy as np
import pandas as pd
import atomica as at

import hbv_functions as hbv
import vimc_functions as vimc
import vimc_plotting as vplt

fw = 'hbv_v14_gamma_vimc.xlsx'
db =  wd_v + '/vimc_AFR_db_v1_3.xlsx'
cl =  wd_v +'/vimc_AFR_calib_v1_2.xlsx'

in_df, cen_df, sto_df = vimc.get_all_results(fw = fw, db = db, cl = cl)

plt.figure()

plt.show()
