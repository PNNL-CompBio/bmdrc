import operator
import pandas as pd
import numpy as np
import re
from abc import abstractmethod

from .BinaryClass import DataClass

class LPRClass(DataClass):
    '''
    Generate a binary class object from light photomotor response data,
    which must be in long format. 

    df: (pandas DataFrame) A dataframe containing columns with the chemical, 
    concentration, plate, well, time, and value. 

    chemical: (string) name of the column containing the chemical IDs, which
    should be strings

    plate: (string) name of the column indicating the plate IDs, which should be
    strings

    well: (string) name of the column with the well IDs, which should be strings

    concentration: (string) name of the column containing the concentrations, which
    should be numerics

    time: (string) name of the column containing time, which should be
    a string or integer. Strings should contain a number. 

    cycle_time (numeric): length of a cycle. Default is 30. 

    value: (string) name of the column containing the binary values, which should 
    be 0 for absent, and 1 for present. Not used if the light photomotor response 
    '''

    # Define the input checking functions 
    def __init__(self, df, chemical, plate, well, concentration, time, value, cycle_length = 30):
        self.df = df
        self.chemical = chemical
        self.plate = plate
        self.well = well
        self.concentration = concentration
        self.time = time
        self.value = value
        self.cycle_length = cycle_length

    # Set property returning functions 
    df = property(operator.attrgetter('_df'))
    chemical = property(operator.attrgetter('_chemical'))
    plate = property(operator.attrgetter('_plate'))
    well = property(operator.attrgetter('_well'))
    concentration = property(operator.attrgetter('_concentration'))
    time = property(operator.attrgetter('_time'))
    value = property(operator.attrgetter('_value'))
    cycle_length = property(operator.attrgetter("_cycle_length"))
    unacceptable = ["bmdrc.Well.ID", "bmdrc.num.tot", "bmdrc.num.nonna", "bmdrc.num.affected", \
                    "bmdrc.Plate.ID", "bmdrc.Endpoint.ID", "bmdrc.filter", "bmdrc.filter.reason", \
                    "bmdrc.frac.affected"]
    
    ################
    ## SET INPUTS ##
    ################

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
            raise Exception(chemicalname + " is not in the column names of df")
        if chemicalname in self.unacceptable:
            raise Exception(chemicalname + " is not a permitted name. Please rename this column.")
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
        if platename in self.unacceptable:
            raise Exception(platename + " is not a permitted name. Please rename this column.")
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
        if wellname in self.unacceptable:
            raise Exception(wellname + " is not a permitted name. Please rename this column.")
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
        if concentrationname in self.unacceptable:
            raise Exception(concentrationname + " is not a permitted name. Please rename this column.")
        self._df[concentrationname] = pd.to_numeric(self._df[concentrationname])
        self._concentration = concentrationname

    @time.setter
    def time(self, timename):
        if not timename: 
           raise Exception("time cannot be empty. Please enter the column name for \
                            the time.")
        if not isinstance(timename, str):
            raise Exception("time must be a name of a column in df.")
        if not timename in self._df.columns:
            raise Exception(timename + " is not in the column names of df")
        if timename in self.unacceptable:
            raise Exception(timename + " is not a permitted name. Please rename this column.")
        self._df[timename] = pd.to_numeric(re.sub("[^0-9]", "", self._df[timename]))
        self._time = timename

    @value.setter
    def value(self, valuename):
        if not valuename: 
            raise Exception("value cannot be empty. Please enter the column name for \
                                the value.")
        if not isinstance(valuename, str):
            raise Exception("value must be a name of a column in df.")
        if not valuename in self._df.columns:
            raise Exception(valuename + " is not in the column names of df")
        if valuename in self.unacceptable:
            raise Exception(valuename + " is not a permitted name. Please rename this column.")
        self._df[valuename] = pd.to_numeric(self._df[valuename])
        self._value = valuename

    @cycle_length.setter
    def cycle_length(self, cycle_length):
        if not cycle_length:
            raise Exception("a cycle_length value must be set.")
        if not isinstance(cycle_length, float):
            raise Exception("cycle_length should be a float")
        self._cycle_length = cycle_length

    




    

