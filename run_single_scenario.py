import atomica as at
at.model.model_settings['tolerance'] = 2e-6
import pandas as pd
import sciris as sc
import os
import numpy as np
import matplotlib.pyplot as plt
import gcat as gc
from pathlib import Path
import itertools
import cat_syphilis as syph

result_dir = syph.root/'results_adjusted_2'
result_dir.mkdir(exist_ok=True)

eff = 'high'
dur = 25
st = 'pn'
sx = 'mf'

# scr = 'more'

# scenario = syph.scenario(eff, dur, st, sx, scr)
# scenario.run((2000, 2051))

# scen_label = '_'.join([str(x) for x in args])
# scen_label = f'{eff}_{dur}_{st}_{sx}_{scr}'


scenario = syph.scenario(eff, dur, st, sx)
scenario.run((2000, 2051))

# scen_label = '_'.join([str(x) for x in args])
scen_label = f'{eff}_{dur}_{st}_{sx}'

cost_df = gc.get_cost_df(scenario)
cost_df.to_excel(result_dir / f"cost_df_{scen_label}.xlsx", merge_cells=False)

epi_df = gc.get_epi_df(scenario)
epi_df.to_excel(result_dir / f"epi_df_{scen_label}.xlsx", merge_cells=False)

bubble_df = gc.get_bubble_rollup_df(scenario)
bubble_df.to_excel(result_dir / f"bubble_df_{scen_label}.xlsx")

df = gc.get_coverage_df(scenario, quantity='coverage_number')
df.to_excel(result_dir / f"coverage_number_{scen_label}.xlsx", merge_cells=False)

df = gc.get_coverage_df(scenario, quantity='coverage_fraction')
df.to_excel(result_dir / f"coverage_fraction_{scen_label}.xlsx", merge_cells=False)