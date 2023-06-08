import os
# wd= 'C:/Users/Phil/Documents/GitHub/hbv-globalinvcase'#set own working directory
# os.chdir(wd)

import atomica as at
import matplotlib.pyplot as plt
import hbv_functions as hbv
hbv.get_gitlab_folder()

region = 'EUR'
F = at.ProjectFramework('frameworks/hbv_v14_gamma_mav.xlsx')
D = at.ProjectData.from_spreadsheet(f"applications/region_{region.lower()}/{region}_db_mav.xlsx", framework=F)
P= at.Project(framework=F, databook=D, sim_dt=0.25, sim_start=1990, sim_end=2099, do_run=False)
cal = P.make_parset()

cal = hbv.pop_calib(P,cal)
cal = hbv.epi_calib(P,cal)
res = P.run_sim(parset=cal, result_name=f'Sample {region} calibration')

plt.figure()
d=at.PlotData([res], outputs=['mtct_inf', 'j_acu:it' ,"j_acu:icl"], pops = 'total')
at.plot_series(d, data=P.data, plot_type = 'stacked')

d=at.PlotData([res], outputs=["prev"])
at.plot_series(d, data=P.data, axis="results")

d=at.PlotData([res], outputs=["cl_acu"])
at.plot_series(d, data=P.data, axis="results")

d=at.PlotData([res], outputs=["cl_cir"])
at.plot_series(d, data=P.data, axis="results")

d=at.PlotData([res], outputs=["cl_hcc"])
at.plot_series(d, data=P.data, axis="results")
plt.show()