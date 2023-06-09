import operator
import pandas as pd
import numpy as np

__author__ = "David Degnan"

class BinaryClass(object):
    '''
    Generates a binary class object where input values are either a 0 or a 1

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

    # Set property returning functions 
    df = property(operator.attrgetter('_df'))
    chemical = property(operator.attrgetter('_chemical'))
    plate = property(operator.attrgetter('_plate'))
    well = property(operator.attrgetter('_well'))
    concentration = property(operator.attrgetter('_concentration'))
    format = property(operator.attrgetter('_format'))
    endpoint = property(operator.attrgetter('_endpoint'))
    value = property(operator.attrgetter('_value'))

    # Now, ensure all other input is correct 

    @df.setter
    def df(self, theDF):
        if theDF.empty:
            raise Exception("df cannot be empty. Please provide a pandas DataFrame.")
        if not isinstance(theDF, pd.DataFrame):
            raise Exception("df must be a pandas DataFrame")
        self._df = theDF

    @chemical.setter
    def chemical(self, chemicalname):
        if not chemicalname: 
           raise Exception("chemical cannot be empty. Please enter the column name for \
                            the chemicals.")
        if not isinstance(chemicalname, str):
            raise Exception("chemical must be a name of a column in df.")
        if not chemicalname in self._df.columns:
            raise Exception(chemicalname + "is not in the column names of df")
        self._chemical = chemicalname

    @plate.setter
    def plate(self, platename):
        if not platename: 
           raise Exception("plate cannot be empty. Please enter the column name for \
                            the plate ids.")
        if not isinstance(platename, str):
            raise Exception("plate must be a name of a column in df.")
        if not platename in self._df.columns:
            raise Exception(platename + " is not in the column names of df")
        self._plate = platename
        
    @well.setter
    def well(self, wellname):
        if not wellname: 
           raise Exception("well cannot be empty. Please enter the column name for \
                            the well ids.")
        if not isinstance(wellname, str):
            raise Exception("well must be a name of a column in df.")
        if not wellname in self._df.columns:
            raise Exception(wellname + " is not in the column names of df")
        self._well = wellname
        
    @concentration.setter
    def concentration(self, concentrationname):
        if not concentrationname: 
           raise Exception("concentration cannot be empty. Please enter the column name for \
                            the concentration.")
        if not isinstance(concentrationname, str):
            raise Exception("concentration must be a name of a column in df.")
        if not concentrationname in self._df.columns:
            raise Exception(concentrationname + " is not in the column names of df")
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
            if not endpointname: 
                raise Exception("endpoint cannot be empty. Please enter the column name for \
                                the endpoint.")
            if not isinstance(endpointname, str):
                raise Exception("endpoint must be a name of a column in df.")
            if not endpointname in self._df.columns:
                raise Exception(endpointname + " is not in the column names of df")
            self._endpoint = endpointname
        else:
            self._endpoint = "endpoint"

        
    @value.setter
    def value(self, valuename):
        if self._format == "value":
            if not valuename: 
                raise Exception("value cannot be empty. Please enter the column name for \
                                    the value.")
            if not isinstance(valuename, str):
                raise Exception("value must be a name of a column in df.")
            if not valuename in self._df.columns:
                raise Exception(valuename + " is not in the column names of df")
            if not np.isin(self._df["value"].unique(), [0,1]).all():
                raise Exception("The value column must be comprised of only zeroes and ones.")
            self._value = valuename
        else:
            self._value = "value"