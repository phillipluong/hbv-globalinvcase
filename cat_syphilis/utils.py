import atomica as at
import functools
import cat_syphilis as syph
import gcat as gc
import numpy as np

framework_version = 'v1_8'
databook_version = 'v1_2'
calibration_version = 'v1_2'

# Scenario dimensions
hiv_pops = ['p', 'pn'] # p: positive, pn: positive and negative
sex_pops = ['f', 'mf'] # f: female, fm: female and male
durations = [2, 5, 10, 25]
vacc_eff = ['direct','low', 'med', 'high']
scrn = ['norm','more']

# Name mappings
hiv_dict = {'p': 'HIV-positive only',
            'pn': 'Both HIV-positive and -negative'}
sex_dict = {'f': 'Females only',
            'mf': 'Females and Males'}
scrn_dict = {'norm': 'Normal Screening',
             'more': 'Increased Screening'}

# Vaccine efficacy levels

eff_acq = {'direct': 0, 'low': 0, 'med': 0, 'high': 0.2} # Efficacy against syphilis acquisition
cs_eff = {'direct': 0.5, 'low': 0.5, 'med': 0.75, 'high': 0.9} # Efficacy against stillbirth and neonatal death
gen_eff = {'direct': 0, 'low': 0.5, 'med': 0.75, 'high': 0.9} # General efficacy - against syphilis transmission (sexual, mother-to-child) and stillbirths/neonatal deaths
hiv_eff = {'direct': 1, 'low': 0.5, 'med': 0.75, 'high': 1}


@functools.lru_cache()
def load_products() -> list:
    products = {}
    for fname in (syph.root/"syp_products").iterdir():
        if not fname.stem.startswith("~"):
            products[fname.stem] = gc.Product.from_spreadsheet(fname)
    return products

def framework():
    return at.ProjectFramework(f'seilt_framework_{framework_version}.xlsx')

def load_country_data(country:str) -> at.ProjectData:
    f = framework()
    d = at.ProjectData.from_spreadsheet(syph.root / f"{country}_databook_{databook_version}.xlsx",framework=f)
    d.validate(f)
    return d


def project(country:str, calibrated = True, eff='high', scr='norm') -> at.Project:

    D = load_country_data(country)

    ## Reduced HIV susceptibility
    for pop in D.pops:
        if pop.endswith('+'):
            D.tdve['hiv_vac_eff'].ts[pop].assumption = hiv_eff[eff]

    if scr == 'more':
        # Scale up treatments
        for ts in D.tdve['trt_likelihood'].ts.values():
            ts.insert(2032, 1)
            ts.insert(2033,1.25)

    P = at.Project(framework=framework(), databook=D, do_run=False)
    P.settings = syph.settings # Use the MNCH standard project years/timestep
    P.country = country  # Add country to project (normally linked in by FE)
    P.data_specification = syph.data_specification

    if calibrated:
        P.parsets[0].load_calibration(f'{country}_cal_{calibration_version}.xlsx')

    return P


def scenario(eff='low', dur=None, st=None, sx=None, scr=None):
    """

    :param eff:
    :param dur:
    :param st:
    :param sx:
    :return:
    """

    # scenario = gc.Scenario(f"{dur} year duration, {eff} efficacy, {sex_dict[sx]}, {hiv_dict[st]}, {scrn_dict[scr]}", "")

    scenario = gc.Scenario(f"{dur} year duration, {eff} efficacy, {sex_dict[sx]}, {hiv_dict[st]}", "")

    scenario.products = [
        gc.Product.from_spreadsheet(f'syp_products/syphilis_proxy_vaccine.xlsx'),
        gc.Product.from_spreadsheet(f'syp_products/syphilis_vaccine_{dur}_{sx}_{st}_{eff}.xlsx'),
        gc.Product.from_spreadsheet(f'syp_products/syphilis_proxy_treatment.xlsx'),
        gc.Product.from_spreadsheet(f'syp_products/congenital_syphilis_treatment.xlsx')
    ]

    for country in syph.countries:
        scenario.projects.append(project(country, calibrated=True, eff=eff, scr=scr))

    scenario.tree = gc.Tree(products=scenario.products, countries=syph.countries, tvec=np.arange(2000, 2051))

    scenario.tree.load_spreadsheet(f'demandbooks/syphilis_demandbook_{dur}_{sx}_{st}_{eff}.xlsx')

    scenario.products = {x.name: x for x in scenario.products}

    return scenario