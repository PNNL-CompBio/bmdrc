import operator
import pandas as pd
import numpy as np
from abc import abstractmethod

from .model_fitting import fit_the_models, gen_response_curve, fits_table
from .output_modules import benchmark_dose, dose_table, report_binary
__author__ = "David Degnan"


class HalfClass(object):
    '''
    Class for fitting any data with a response from 0 to 1. 
    Does not contain the pre-processing & filtering options of BinaryClass.
    '''

    ###########################
    ## MODEL FITTING MODULES ##
    ###########################

    @abstractmethod
    def fit_models(self, gof_threshold = 0.1, aic_threshold = 2, model_selection = "lowest BMDL", diagnostic_mode = False):
        fit_the_models(self, gof_threshold, aic_threshold, model_selection, diagnostic_mode)

    @abstractmethod
    def response_curve(self, chemical_name, endpoint_name, model, steps = 10):
        gen_response_curve(self, chemical_name, endpoint_name, model, steps)

    ####################
    ## OUTPUT MODULES ##
    ####################

    @abstractmethod
    def output_benchmark_dose(self, path = None):
        benchmark_dose(self, path)

    @abstractmethod
    def output_dose_table(self, path = None):
        dose_table(self, path)

    @abstractmethod
    def output_fits_table(self, path = None):
        fits_table(self, path)

    @abstractmethod
    def report(self, out_folder, report_name = "Benchmark Dose Curves"):
        report_binary(self, out_folder, report_name)

class ResponseOnlyClass(HalfClass):
    '''
    Class for fitting any data with a response from 0 to 1. 
    Does not contain the pre-processing & filtering options of BinaryClass.

    df: (pandas DataFrame) A dataframe containing columns title chemical, 
    plate, well, concentration, endpoint (long format only), value (long format only).
    If the data is in wide format, all additional columns are assumed to be endpoints.

    chemical: (string) name of the column containing the chemical IDs, which
    should be strings

    concentration: (string) name of the column containing the concentrations, which
    should be numerics

    endpoint: (string) name of the column containing endpoints, which should be
    a string. 

    response: (string) name of the column containing the response values, which should 
    range from 0 to 1.
    '''

    # Define the input checking functions 
    def __init__(self, df, chemical, concentration, endpoint, response):
        self.df = df
        self.chemical = chemical
        self.concentration = concentration
        self.endpoint = endpoint
        self.response = response
        self.unacceptable = ["bmdrc.Well.ID", "bmdrc.num.tot", "bmdrc.num.nonna", "bmdrc.num.affected", \
                            "bmdrc.Plate.ID", "bmdrc.Endpoint.ID", "bmdrc.filter", "bmdrc.filter.reason", \
                            "bmdrc.frac.affected"]
    
    # Set property returning functions 
    df = property(operator.attrgetter('_df'))
    chemical = property(operator.attrgetter('_chemical'))
    concentration = property(operator.attrgetter('_concentration'))
    endpoint = property(operator.attrgetter('_endpoint'))
    response = property(operator.attrgetter('_response'))
    unacceptable = ["bmdrc.Well.ID", "bmdrc.num.tot", "bmdrc.num.nonna", "bmdrc.num.affected", \
                    "bmdrc.Plate.ID", "bmdrc.Endpoint.ID", "bmdrc.filter", "bmdrc.filter.reason", \
                    "bmdrc.frac.affected"]
    
    ################
    ## SET INPUTS ##
    ################

    @df.setter
    def df(self, theDF):
        if not isinstance(theDF, pd.DataFrame):
            raise Exception("df must be a pandas DataFrame.")
        if theDF.empty:
            raise Exception("df cannot be empty. Please provide a pandas DataFrame.")
        self._df = theDF

    @chemical.setter
    def chemical(self, chemicalname):
        if not isinstance(chemicalname, str):
            raise Exception("chemical must be a name of a column in df.")
        if not chemicalname in self._df.columns:
            raise Exception(chemicalname + " is not in the column names of df.")
        if chemicalname in self.unacceptable:
            raise Exception(chemicalname + " is not a permitted name. Please rename this column.")
        
        # Make this column a string
        self._df[chemicalname] = self._df[chemicalname].astype(str)

        if any([" " in x for x in self._df[chemicalname].unique()]):
            raise Exception("spaces are not permitted in chemical names. Check this column in your dataframe.")
        self._chemical = chemicalname

    @concentration.setter
    def concentration(self, concentrationname):
        if not isinstance(concentrationname, str):
            raise Exception("concentration must be a name of a column in df.")
        if not concentrationname in self._df.columns:
            raise Exception(concentrationname + " is not in the column names of df.")
        if concentrationname in self.unacceptable:
            raise Exception(concentrationname + " is not a permitted name. Please rename this column.")
        self._df[concentrationname] = pd.to_numeric(self._df[concentrationname])
        self._concentration = concentrationname

    @response.setter
    def response(self, responsename):
        if not isinstance(responsename, str):
            raise Exception("response must be a name of a column in df.")
        if not responsename in self._df.columns:
            raise Exception(responsename + " is not in the column name of df.")
        if responsename in self.unacceptable:
            raise Exception(responsename + " is not a permitted name. Please rename this column.")

