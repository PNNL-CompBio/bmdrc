import numpy as np
import pandas as pd
import scipy.stats as stats
from astropy import stats as astrostats

import BMD_Analysis_Functions as baf 

__author__ = "David Degnan"

def fit_the_models(self, models, fit_threshold, BMD_Measurements):
    '''
    Fit the EPA recommended models to your dataset. 

    plate_groups: (pandas DataFrame) A plate_groups object with columns to indicate filters.
    See "filter modules" for more details. 
    
    
    '''