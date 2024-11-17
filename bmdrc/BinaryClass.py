import operator
import pandas as pd
import numpy as np
from abc import abstractmethod

from .preprocessing import endpoint_combine, well_to_na, remove_endpoints
from .filtering import make_plate_groups, negative_control, min_concentration, correlation_score
from .model_fitting import fit_the_models, gen_response_curve, fits_table
from .output_modules import benchmark_dose, dose_table, report_binary

__author__ = "David Degnan"

## Guidance for all selections can be found here: https://www.epa.gov/sites/default/files/2015-01/documents/benchmark_dose_guidance.pdf
## Package built following the guidance from here: https://packaging.python.org/en/latest/tutorials/packaging-projects/

class DataClass(object):
    '''
    An abstract class for all bmdrc accepted datatypes
    '''
    
    ############################
    ## PRE-PROCESSING MODULES ##
    ############################

    @abstractmethod
    def set_well_to_na(self, endpoint_name, endpoint_value, except_endpoint = None):
        well_to_na(self, endpoint_name, endpoint_value, except_endpoint)

    @abstractmethod
    def combine_and_create_new_endpoints(self, endpoint_dict):
        endpoint_combine(self, endpoint_dict)

    @abstractmethod
    def remove_endpoints(self, endpoint_name):
        remove_endpoints(self, endpoint_name)

    #######################
    ## FILTERING MODULES ##
    #######################

    @abstractmethod
    def make_plate_groups(self):
        make_plate_groups(self)

    @abstractmethod
    def filter_negative_control(self, percentage = 50, apply = False, diagnostic_plot = False):
        negative_control(self, percentage, apply, diagnostic_plot)

    @abstractmethod
    def filter_min_concentration(self, count = 3, apply = False, diagnostic_plot = False):
        min_concentration(self, count, apply, diagnostic_plot)

    @abstractmethod
    def filter_correlation_score(self, score = 0.2, apply = False, diagnostic_plot = False):
        correlation_score(self, score, apply, diagnostic_plot)

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

class BinaryClass(DataClass):
    '''
    Generates a binary class object where input values are either a 0, 1, or NA.

    df: (pandas DataFrame) A dataframe containing columns title chemical, 
    plate, well, concentration, endpoint (long format only), value (long format only).
    If the data is in wide format, all additional columns are assumed to be endpoints.

    chemical: (string) name of the column containing the chemical IDs, which
    should be strings

    plate: (string) name of the column indicating the plate IDs, which should be
    strings

    well: (string) name of the column with the well IDs, which should be strings

    concentration: (string) name of the column containing the concentrations, which
    should be numerics

    endpoint: (string) name of the column containing endpoints, which should be
    a string. Not used if the data is in wide format. 

    value: (string) name of the column containing the binary values, which should 
    be 0 for absent, and 1 for present. Not used if the data is in wide format.

    format: (string) indicate whether the data is in 'long' or 'wide' format. Wide
    format requires only the chemical, plate, well, and concentration columns.
    The rest of the columns are assumed to be endpoints. Wide formats are then converted
    to the long format. 
    '''

    # Define the input checking functions 
    def __init__(self, df, chemical, plate, well, concentration, endpoint = None, value = None, format = "long"):
        self.df = df
        self.chemical = chemical
        self.plate = plate
        self.well = well
        self.concentration = concentration
        self.format = format
        self.endpoint = endpoint
        self.value = value
        self.unacceptable = ["bmdrc.Well.ID"]

    # Set property returning functions 
    df = property(operator.attrgetter('_df'))
    chemical = property(operator.attrgetter('_chemical'))
    plate = property(operator.attrgetter('_plate'))
    well = property(operator.attrgetter('_well'))
    concentration = property(operator.attrgetter('_concentration'))
    format = property(operator.attrgetter('_format'))
    endpoint = property(operator.attrgetter('_endpoint'))
    value = property(operator.attrgetter('_value'))
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
        self._chemical = chemicalname

    @plate.setter
    def plate(self, platename):
        if not isinstance(platename, str):
            raise Exception("plate must be a name of a column in df.")
        if not platename in self._df.columns:
            raise Exception(platename + " is not in the column names of df.")
        if platename in self.unacceptable:
            raise Exception(platename + " is not a permitted name. Please rename this column.")
        self._plate = platename
        
    @well.setter
    def well(self, wellname):
        if not isinstance(wellname, str):
            raise Exception("well must be a name of a column in df.")
        if not wellname in self._df.columns:
            raise Exception(wellname + " is not in the column names of df.")
        if wellname in self.unacceptable:
            raise Exception(wellname + " is not a permitted name. Please rename this column.")
        self._well = wellname
        
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

    # The format variable by default is long. If the data is wide, it needs to be
    # pivoted.
    @format.setter
    def format(self, long_or_wide):
        if not (long_or_wide == "wide" or long_or_wide == "long"):
            raise Exception("format must be 'long' or 'wide'.")
        if long_or_wide == "wide":
            self._df = self._df.melt(id_vars = [self._chemical, self._concentration, self._plate, self._well], var_name = "endpoint")
        self._format = long_or_wide
        
    @endpoint.setter
    def endpoint(self, endpointname):
        if self._format == "long":
            if not isinstance(endpointname, str):
                raise Exception("endpoint must be a name of a column in df.")
            if not endpointname in self._df.columns:
                raise Exception(endpointname + " is not in the column names of df.")
            if endpointname in self.unacceptable:
                raise Exception(endpointname + " is not a permitted name. Please rename this column.")
            self._endpoint = endpointname
        else:
            self._endpoint = "endpoint"
        
    @value.setter
    def value(self, valuename):
        if self._format == "long":
            if not isinstance(valuename, str):
                raise Exception("value must be a name of a column in df.")
            if not valuename in self._df.columns:
                raise Exception(valuename + " is not in the column names of df.")
            if valuename in self.unacceptable:
                raise Exception(valuename + " is not a permitted name. Please rename this column.")
            
            # Save the unique values in this list
            theValues = self._df[valuename].unique().tolist()

            # Remove NAs (these are acceptable)
            theValues = [int(val) for val in theValues if np.isnan(val) == False]

            if not set(theValues) == {0, 1}:
                raise Exception("The value column must be comprised of only zeroes, ones, and NA values.")
            self._value = valuename
        else:
            self._value = "value"
