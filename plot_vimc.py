import os
wd= 'C:/Users/iamph/Documents/GitHub/hbv-globalinvcase'#set own working directory
os.chdir(wd)

import numpy as np
import atomica as at
import pandas as pd

import hbv_functions as hbv

F = at.ProjectFramework()
D = at.ProjectData()
P = at.Project()
cal = P.make_parset()

d=at.PlotData([res], outputs=["alive"])
at.plot_series(d, data=P.data, axis="results")

d=at.PlotData([res], outputs=["cl_cir"])
at.plot_series(d, data=P.data, axis="results")

d=at.PlotData([res], outputs=["decomp"])
at.plot_series(d, data=P.data, axis="results")

d=at.PlotData([res], outputs=["prev"])
at.plot_series(d, data=P.data, axis="results")

d=at.PlotData([res], outputs=["cl_acu"])
at.plot_series(d, data=P.data, axis="results")

d=at.PlotData([res], outputs=["eag_ott"])
at.plot_series(d, data=P.data, axis="results")

'''
Other plots:
- Total number of cases
- Total Deaths
- Total DALYs
- Compare baseline vs novax
'''

