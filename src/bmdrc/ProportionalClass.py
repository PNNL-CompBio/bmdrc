import operator
import pandas as pd
from abc import abstractmethod

from .preprocessing import remove_endpoints
from .filtering import min_concentration, correlation_score
from .model_fitting import fit_the_models, gen_response_curve, fits_table
from .output_modules import benchmark_dose, dose_table, report_binary

__author__ = "David Degnan"


class ProportionalClass():
    '''
    Generates a bmdrc object from proportions (ranging from 0 to 1). Does not contain the pre-processing & filtering options of BinaryClass.
    Data must be in long format.

    Parameters
    ----------
    df
        A pandas dataframe containing columns with chemical, concentration, endpoint, and response information.

    chemical
        A string indicating the name of the column containing the chemical IDs, which should be strings

    concentration
        A string indicating the name of the column containing the concentrations, which should be numerics

    endpoint
        A string indicating the name of the column containing endpoints, which should be strings. 

    response
        A string indicating the name of the column containing the response values, which should range from 0 to 1.

    '''

    ############################
    ## PRE-PROCESSING MODULES ##
    ############################

    @abstractmethod
    def remove_endpoints(self, endpoint_name):
        remove_endpoints(self, endpoint_name)

    #######################
    ## FILTERING MODULES ##
    #######################

    @abstractmethod
    def filter_min_concentration(self, count = 3, apply = False, diagnostic_plot = False):
        min_concentration(self, count, apply, diagnostic_plot)

    @abstractmethod
    def filter_correlation_score(self, score = 0.2, apply = False, diagnostic_plot = False, direction = "below"):
        correlation_score(self, score, apply, diagnostic_plot, direction)

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

    #####################
    ## INIT DEFINITION ##
    #####################
    
    # Define the input checking functions 
    def __init__(self, df, chemical, concentration, endpoint, response):
        self.df = df
        self.chemical = chemical
        self.concentration = concentration
        self.endpoint = endpoint
        self.response = response
        self.plate = "plate"  # Placeholder to enable specific functions to run
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
        theDF["plate"] = "NoPlate" # There is no plate information to track
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

    @endpoint.setter
    def endpoint(self, endpointname):
        if not isinstance(endpointname, str):
            raise Exception("endpoint must be a name of a column in df.")
        if not endpointname in self._df.columns:
            raise Exception(endpointname + " is not in the column names of df.")
        if endpointname in self.unacceptable:
            raise Exception(endpointname + " is not a permitted name. Please rename this column.")
        
        # Make this column a string
        self._df[endpointname] = self._df[endpointname].astype(str)

        if any([" " in x for x in self._df[endpointname].unique()]):
            raise Exception("spaces are not permitted in endpoint names. Check this column in your dataframe.")
        self._endpoint = endpointname

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
        self._df[responsename] = pd.to_numeric(self._df[responsename])
        if min(self._df[responsename]) < 0 or max(self._df[responsename]) > 1:
            raise Exception("The response column must range in values from 0 to 1. Filter out NAs before using the ProportionalClass.")
        self._response = responsename
