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

    value: (string) name of the column containing the binary values, which should 
    be 0 for absent, and 1 for present. Not used if the light photomotor response 

    cycle_time: (numeric) length of a light or dark cycle. Default is 20. 
    The unit is a 6-second measure, so 20 six second measures is 2 minutes.

    cycle_cooldown: (numeric) length of time between cycles. Default is 10.
    The unit is a 6-second measure, so 10 six second measures is 1 minute. 

    starting_cycle: (string) either "light" or "dark" depending on whether
    the first measurement was a light or dark cycle. Default is "light". 

    samples_to_remove: (string) list of strings with plate and well ids to remove, written
    as "plate.id & well". For example if plate 2 and well H01 was to be removed, it should
    be written as "2 & H01"  # TO Update
    '''

    # Define the input checking functions. Include raw and transformed data.frames 
    def __init__(self, df, chemical, plate, well, concentration, time, value, cycle_length = 20.0, 
                 cycle_cooldown = 10.0, starting_cycle = "light", samples_to_remove = None):
        self.df = df
        self.chemical = chemical
        self.plate = plate
        self.well = well
        self.concentration = concentration
        self.time = time
        self.value = value
        self.cycle_length = cycle_length
        self.cycle_cooldown = cycle_cooldown
        self.starting_cycle = starting_cycle
        self.samples_to_remove = samples_to_remove
        self.add_cycles()

    # Set property returning functions 
    df = property(operator.attrgetter('_df'))
    chemical = property(operator.attrgetter('_chemical'))
    plate = property(operator.attrgetter('_plate'))
    well = property(operator.attrgetter('_well'))
    concentration = property(operator.attrgetter('_concentration'))
    time = property(operator.attrgetter('_time'))
    value = property(operator.attrgetter('_value'))
    cycle_length = property(operator.attrgetter('_cycle_length'))
    cycle_cooldown = property(operator.attrgetter('_cycle_cooldown'))
    starting_cycle = property(operator.attrgetter('_starting_cycle'))
    samples_to_remove = property(operator.attrgetter('_samples_to_remove'))
    cycles = property(operator.attrgetter('_cycle'))
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
        self._ori_df = theDF
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
        self._df[timename] = self._df[timename].str.extract('(\d+)', expand=False).astype(float)
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
        self._df[valuename] = self._df[valuename].astype(float)
        self._value = valuename

    @cycle_length.setter
    def cycle_length(self, cycle_length):
        if not cycle_length:
            raise Exception("a cycle_length value must be set.")
        if not isinstance(cycle_length, float):
            raise Exception("cycle_length should be a float")
        self._cycle_length = cycle_length

    @cycle_cooldown.setter
    def cycle_cooldown(self, cycle_cooldown):
        if not cycle_cooldown:
            raise Exception("a cycle_cooldown value must be set.")
        if not isinstance(cycle_cooldown, float):
            raise Exception("cycle_length should be a float")
        self._cycle_cooldown = cycle_cooldown

    @starting_cycle.setter
    def starting_cycle(self, starting_cycle):
        if not starting_cycle:
            raise Exception("starting_cycle must be either 'light' or 'dark'.")
        if not starting_cycle in ['light', 'dark']:
            raise Exception("starting_cycle must be either 'light' or 'dark'.")
        self._starting_cycle = starting_cycle

    @samples_to_remove.setter
    def samples_to_remove(self, samples_to_remove):
        if samples_to_remove:
            self._df = self._df[~(self._df["well"].isin(samples_to_remove))]
            self._samples_to_remove = samples_to_remove

    # Write LPR specific formatting function
    def add_cycles(self):
        '''Specific LPR function that adds cycle information following users setting cycle_time,
        cycle_cooldown, and samples_to_remove'''

        print("...defining cycles")

        # Unique and arrange times 
        cycle_info = pd.DataFrame(self._df[self._time].unique()).rename({0:self._time}, axis = 1)
        cycle_info[self._time] = cycle_info[self._time].astype(float)
        cycle_info = cycle_info.sort_values(by = [self._time])

        # Build cycle names and order. First, define all the needed variables to make this happen
        cycle_order = []
        first_count = 0
        gap_a_count = 0
        gap_b_count = 0
        second_count = 0
        cycle_count = 1
        if self._starting_cycle == "light":
            other_cycle = "dark"
        else:
            other_cycle = "light"

        # Cycle through the light, gap, dark, and then reset
        for pos in range(len(cycle_info)):
            if (first_count < self._cycle_length):
                cycle_order.append(self._starting_cycle + str(cycle_count))
                first_count += 1
            elif (gap_a_count < self._cycle_cooldown):
                cycle_order.append("gap_" + self._starting_cycle + str(cycle_count))
                gap_a_count += 1
            elif (second_count < self._cycle_length):
                cycle_order.append(other_cycle + str(cycle_count))
                second_count += 1
            elif (gap_b_count < self._cycle_cooldown):
                cycle_order.append("gap_" + other_cycle + str(cycle_count))
                gap_b_count += 1
            else:
                cycle_count += 1
                cycle_order.append(self._starting_cycle + str(cycle_count))
                first_count = 1
                gap_a_count = 0
                second_count = 0
                gap_b_count = 0
            
        # Add essential order information to cycle_info file
        cycle_info["cycle"] = cycle_order

        # Merge with data.frame
        self._cycles = cycle_info
        self._df = self._df.merge(cycle_info)






    




    

