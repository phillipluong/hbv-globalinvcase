import atomica as at
import gcat as gc
import numpy as np

analysis_years = [2000, 2050]
dt=1 / 12
analysis_years[1] = analysis_years[0] + np.floor(((analysis_years[1] + 1 - dt * 1.1) - analysis_years[0]) / dt) * dt

# analysis_years[1] = analysis_years[0] + np.floor(((analysis_years[1] + 1 - dt * 1.1) - analysis_years[0]) / dt) * dt

settings = at.ProjectSettings(sim_start=analysis_years[0], sim_end=analysis_years[1], sim_dt=dt)  #: Atomica project settings

countries = {"UGA": 'Uganda'}

# Load projects
pops = [
    gc.dataspecification.PopulationSpecification('<1y-', "<1 year, M+F, HIV negative", age_min=0, age_max=0, male=True, female=True),
    gc.dataspecification.PopulationSpecification("<1y+", "<1y year M+F, HIV positive", age_min=0, age_max=0, male=True, female=True),
    gc.dataspecification.PopulationSpecification("1-9y-", "1-9 years, M+F, HIV negative", age_min=1, age_max=9, male=True, female=True),
    gc.dataspecification.PopulationSpecification("1-9y+", "1-9 years, M+F, HIV positive", age_min=1, age_max=9, male=True, female=True),
    gc.dataspecification.PopulationSpecification("10-17M-", "10-17 years, Male, HIV negative", age_min=10, age_max=17, male=True, female=False),
    gc.dataspecification.PopulationSpecification("10-17M+", "10-17 years, Male, HIV positive", age_min=10, age_max=17, male=True, female=False),
    gc.dataspecification.PopulationSpecification("10-17F-", "10-17 years, Female, HIV negative", age_min=10, age_max=17, male=False, female=True),
    gc.dataspecification.PopulationSpecification("10-17F+", "10-17 years, Female, HIV positive", age_min=10, age_max=17, male=False, female=True),
    gc.dataspecification.PopulationSpecification("18-29M-", "18-29 years, Male, HIV negative", age_min=18, age_max=29, male=True, female=False),
    gc.dataspecification.PopulationSpecification("18-29M+", "18-29 years, Male, HIV positive", age_min=18, age_max=29, male=True, female=False),
    gc.dataspecification.PopulationSpecification("18-29F-", "18-29 years, Female, HIV negative", age_min=18, age_max=29, male=False, female=True),
    gc.dataspecification.PopulationSpecification("18-29F+", "18-29 years, Female, HIV positive", age_min=18, age_max=29, male=False, female=True),
    gc.dataspecification.PopulationSpecification("30+M-", "30+ years, Male, HIV negative", age_min=30, age_max=np.inf, male=True, female=False),
    gc.dataspecification.PopulationSpecification("30+M+", "30+ years, Male, HIV positive", age_min=30, age_max=np.inf, male=True, female=False),
    gc.dataspecification.PopulationSpecification("30+F-", "30+ years, Female, HIV negative", age_min=30, age_max=np.inf, male=False, female=True),
    gc.dataspecification.PopulationSpecification("30+F+", "30+ years, Female, HIV positive", age_min=30, age_max=np.inf, male=False, female=True)
]

data_specification = gc.dataspecification.DataSpecification(pops=pops)


