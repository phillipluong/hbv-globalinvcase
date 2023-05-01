import os
wd= 'C:/Users/Phil/Documents/GitHub/hbv-globalinvcase'#set own working directory
wd_v = wd + '/applications/vimc'
wd_data = 'C:/Users/Phil/Google Drive (phillip.luong@burnet.edu.au)/Projects/2022-12 VIMC HepB Modelling/2023-04_Wk1_result_archives/'
os.chdir(wd)

import pandas as pd
import matplotlib.pyplot as plt
import vimc_plotting as vplt

in_df = pd.read_csv(wd_data + 'input_results_novax_230130_2.csv')
cen_df = pd.read_csv(wd_data + 'central_results_novax_230130.csv')
sto_df = pd.read_csv(wd_data + 'stochastic_results_novax_230130.csv')

vplt.plot_all(cen_df, sto_df, 10)
vplt.plot_bounds(cen_df, sto_df, 10, 'cases')

plt.show()

