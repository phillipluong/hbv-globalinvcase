import atomica as at
at.model.model_settings['tolerance'] = 2e-6

import sciris as sc
import numpy as np
import gcat as gc
from pathlib import Path
import cat_syphilis as syph

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



scenario = syph.scenario("high",5,"p","f", "norm")
scenario.run((2000, 2051))


# Need to plot

