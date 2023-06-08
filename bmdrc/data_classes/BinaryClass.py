import operator
import pandas as pd

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
    def __init__(self, df, chemical, plate, well, concentration, endpoint, value): #format):
        self.df = df
        self.chemical = chemical
        self.plate = plate
        self.well = well
        self.concentration = concentration
        self.endpoint = endpoint
        self.value = value
        #self.format = format

    # Set property returning functions 
    df = property(operator.attrgetter('_df'))
    chemical = property(operator.attrgetter('_chemical'))
    plate = property(operator.attrgetter('_plate'))
    well = property(operator.attrgetter('_well'))
    concentration = property(operator.attrgetter('_concentration'))
    endpoint = property(operator.attrgetter('_endpoint'))
    value = property(operator.attrgetter('_value'))
    #format = property(operator.attrgetter('_format'))

    # Now, ensure each input is correct 

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
        
    @endpoint.setter
    def endpoint(self, endpoints):
    #    if not endpoints: 
    #        raise Exception("endpoint cannot be empty. Please provide a list \
    #                        of strings representing each endpoint.")
    #    if not isinstance(endpoints, str):
    #        raise Exception("endpoint must be a list of strings.")
        self._endpoint = endpoints
        
        
    @value.setter
    def value(self, values):
    #    if not values: 
     #       raise Exception("value cannot be empty. Please provide a list \
     #                       of zeroes and ones representing whether a condition \
     #                       is absent or present, respectively.")
     #   if all(val == 0 or val == 1 for val in values):
     #       raise Exception("value must be a list of zeroes or ones indicating \
     #                       whether a condition is absent or present, respectively.")
        self._value = values