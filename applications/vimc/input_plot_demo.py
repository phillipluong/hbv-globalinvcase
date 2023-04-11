import pandas as pd
import matplotlib.pyplot as plt

wd_data = 'C:/Users/Phil/Google Drive (phillip.luong@burnet.edu.au)/Projects/2022-12 VIMC HepB Modelling/2023-04_Wk1_result_archives/'


in_n = pd.read_csv(wd_data + 'input_results_vax_230130_2.csv')

# possible parameters
nathis = ['ci_p', 'm_acu', 'm_dc', 'm_hcc']
trtinv = ['te_dc_cc', 'te_icl_ict', 'te_icl_cc', 'te_cc_dc', 'te_cc_hcc', 'te_ie_cc', 'te_m_dc', 'te_m_hcc', 'te_ict_hcc', 'te_ie_hcc', 'te_icl_hcc', 'te_dc_hcc']
vacinv = ['eag_ve', 'sag_ve', 'hb3_ve', 'mav_ve']

# possible ages
ages = ['0-4', '5-14', '15-49', '50-69', '70+']

# Possible sexes
sexes = ['M', 'F']

# State an input parameter, age group, and sex that you would like to observe
input_par = 'te_m_dc'
age = '5-14'
sex = 'M'

in_n.boxplot(column = [f'{input_par}_{age}{sex}'])
plt.show()