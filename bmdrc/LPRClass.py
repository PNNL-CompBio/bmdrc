import operator
import pandas as pd
import numpy as np
import re
from abc import abstractmethod

from .BinaryClass import DataClass, BinaryClass

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

    '''

    # Define the input checking functions. Include raw and transformed data.frames 
    def __init__(self, df, chemical, plate, well, concentration, time, value, cycle_length = 20.0, 
                 cycle_cooldown = 10.0, starting_cycle = "light"):
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
        self.convert_LPR()


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
    cycles = property(operator.attrgetter('_cycle'))
    unacceptable = ["bmdrc.Well.ID", "bmdrc.num.tot", "bmdrc.num.nonna", "bmdrc.num.affected", \
                    "bmdrc.Plate.ID", "bmdrc.Endpoint.ID", "bmdrc.filter", "bmdrc.filter.reason", \
                    "bmdrc.frac.affected", "cycle"]
    
    ################
    ## SET INPUTS ##
    ################

    @df.setter
    def df(self, theDF):
        if not isinstance(theDF, pd.DataFrame):
            raise Exception("df must be a pandas DataFrame.")
        if theDF.empty:
            raise Exception("df cannot be empty. Please provide a pandas DataFrame.")
        self._ori_df = theDF
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

    @time.setter
    def time(self, timename):
        if not isinstance(timename, str):
            raise Exception("time must be a name of a column in df.")
        if not timename in self._df.columns:
            raise Exception(timename + " is not in the column names of df.")
        if timename in self.unacceptable:
            raise Exception(timename + " is not a permitted name. Please rename this column.")
        self._df[timename] = self._df[timename].str.extract('(\d+)', expand=False).astype(float)
        self._time = timename

    @value.setter
    def value(self, valuename):
        if not isinstance(valuename, str):
            raise Exception("value must be a name of a column in df.")
        if not valuename in self._df.columns:
            raise Exception(valuename + " is not in the column names of df.")
        if valuename in self.unacceptable:
            raise Exception(valuename + " is not a permitted name. Please rename this column.")
        self._df[valuename] = self._df[valuename].astype(float)
        self._value = valuename

    @cycle_length.setter
    def cycle_length(self, cycle_length):
        if not isinstance(cycle_length, float):
            raise Exception("cycle_length should be a float.")
        self._cycle_length = cycle_length

    @cycle_cooldown.setter
    def cycle_cooldown(self, cycle_cooldown):
        if not isinstance(cycle_cooldown, float):
            raise Exception("cycle_cooldown should be a float.")
        self._cycle_cooldown = cycle_cooldown

    @starting_cycle.setter
    def starting_cycle(self, starting_cycle):
        if not starting_cycle in ['light', 'dark']:
            raise Exception("starting_cycle must be either 'light' or 'dark'.")
        self._starting_cycle = starting_cycle

    # LPR-specific function: determines light and dark cycles
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
        self._max_cycle = cycle_count
        return(self._df.merge(cycle_info))

    # LPR-specific function: Converts continuous to dichotomous following: https://www.sciencedirect.com/science/article/pii/S2468111318300732
    def to_dichotomous(self, the_df, the_value):
        '''Specific LPR function that converts continuous AUC or MOV to dichotomous'''

        LPR_plateGroups = the_df[[self._chemical, self._plate, self._concentration, the_value]]
        LPR_zero = LPR_plateGroups[LPR_plateGroups[self._concentration] == 0].groupby([self._chemical, self._plate])

        # Pull quartile calculations
        rangeValues = LPR_zero.apply(lambda df: df[the_value].quantile(0.25)).reset_index().rename(columns = {0:"Q1"})
        rangeValues["Q3"] = LPR_zero.apply(lambda df: df[the_value].quantile(0.75)).reset_index()[0]

        # Add IQR and lower and upper bonds. 
        rangeValues["IQR"] = rangeValues["Q3"] - rangeValues["Q1"]
        rangeValues["Low"] = rangeValues["Q1"] - (1.5 * rangeValues["IQR"])
        rangeValues["High"] = rangeValues["Q3"] + (1.5 * rangeValues["IQR"])
        rangeValues = rangeValues[[self._chemical, self._plate, "Low", "High"]]

        LPR_plateGroups = pd.merge(LPR_plateGroups, rangeValues)
        LPR_plateGroups["result"] = (LPR_plateGroups[the_value] < 0) | (LPR_plateGroups[the_value] < LPR_plateGroups["Low"]) | (LPR_plateGroups[the_value] > LPR_plateGroups["High"])

        return(LPR_plateGroups["result"].astype(float))
    
    # LPR-specific function: Calculate AUC values
    def calculate_aucs(self, cycles):
        '''Specific LPR function for calculating AUC values'''

        print("...calculating AUC values")

        # Remove gaps
        aucs = cycles[~cycles["cycle"].str.contains("gap")].drop(labels = self._time, axis = 1)
        aucs = aucs.groupby(by = [self._chemical, self._concentration, self._plate, self._well, "cycle"]).sum().reset_index()

        # Initiate list to store all values
        store_aucs = []

        # Iterate through all cycles, subtracting dark from light
        for cycle_num in range(self._max_cycle):

            # Pull cycle information
            cycle_num = cycle_num + 1
            light_name = "light" + str(cycle_num)
            dark_name = "dark" + str(cycle_num)

            # Merge light and dark information
            to_calc_auc = pd.merge(
                aucs[aucs["cycle"] == light_name].rename(columns = {"value":"light"}).drop("cycle", axis = 1),
                aucs[aucs["cycle"] == dark_name].rename(columns = {"value":"dark"}).drop("cycle", axis = 1),
                how = "left"
            )

            # Sore calculate aucs 
            store_aucs.append([to_calc_auc["dark"][x] - to_calc_auc["light"][x] for x in range(len(to_calc_auc))])

        # Convert to a data.frame
        auc_values = pd.DataFrame(store_aucs).transpose()

        # Rename columns 
        for name in auc_values.columns:
            new_name = "AUC" + str(int(name) + 1)
            auc_values.rename(columns={name:new_name}, inplace = True)

        # Add additional columns that are needed
        auc_process = pd.concat([
            cycles[[self._chemical, self._concentration, self._plate, self._well]].drop_duplicates(),
            auc_values
        ], axis = 1)

        # Convert to dichotomous
        for x in range(self._max_cycle):
            value = "AUC" + str(x + 1)
            auc_process[value] = self.to_dichotomous(auc_process, value)

        # Return AUC
        return(auc_process)

    # LPR-specific function: Calculate MOV values
    def calculate_movs(self, cycles):
        '''Specific LPR function for calculating MOV values'''

        print("...calculating MOV values")

        # Determine gap beginning and end times 
        if self._starting_cycle == "light":
            self._light_gaps = [(self._cycle_length * (x + 1)) + (x * self._cycle_length * 2) for x in range(self._max_cycle)]
            self._dark_gaps = [(self._cycle_length * (x + 1) + self._cycle_cooldown) + (x * self._cycle_length * 2) for x in range(self._max_cycle)]
        else:
            self._light_gaps = [(self._cycle_length * (x + 1) + self._cycle_cooldown) + (x * self._cycle_length * 2) for x in range(self._max_cycle)]
            self._dark_gaps = [(self._cycle_length * (x + 1)) + (x * self._cycle_length * 2) for x in range(self._max_cycle)]

        # Select and sum values 
        movs = cycles[(cycles[self._time].isin(self._light_gaps)) | (cycles[self._time].isin(self._dark_gaps))]

        # Initiate list to store all calculate values 
        store_movs = []

        # Iterate through all cycles, subtracting dark from light
        for x in range(self._max_cycle):

            # Merge light and dark information
            to_calc_mov = pd.merge(
                movs[movs[self._time] == self._light_gaps[x]].rename(columns = {"value":"light"}).drop([self._time, "cycle"], axis = 1),
                movs[movs[self._time] == self._dark_gaps[x]].rename(columns = {"value":"dark"}).drop([self._time, "cycle"], axis = 1)
            )

            # Store calculated MOVs
            store_movs.append([to_calc_mov["dark"][x] - to_calc_mov["light"][x] for x in range(len(to_calc_mov))]) 

        # Make data.frame
        mov_values = pd.DataFrame(store_movs).transpose()

        # Rename columns
        for name in mov_values.columns:
            new_name = "MOV" + str(int(name) + 1)
            mov_values.rename(columns={name:new_name}, inplace = True)

        # Add additional columns that are needed
        mov_process = pd.concat([
            cycles[[self._chemical, self._concentration, self._plate, self._well]].drop_duplicates(),
            mov_values
        ], axis = 1)

        # Convert to dichotomous
        for x in range(self._max_cycle):
            value = "MOV" + str(x + 1)
            mov_process[value] = self.to_dichotomous(mov_process, value)

        # Return AUC
        return(mov_process)

    # LPR-specific: Convert LPR continuous to Dichotomous 
    def convert_LPR(self):
        '''Wrapper function for all LPR-specific functions for calculating cycles,
        AUC values, MOV values, and converting them to dichotomous values'''

        # Step 1: Make Cycle Information
        CycleInfo = self.add_cycles()

        # Step 2: Calculate AUC values
        AUCs = self.calculate_aucs(CycleInfo)

        # Step 3: Calculate MOV values
        MOVs = self.calculate_movs(CycleInfo)

        # Step 4: Merge the results
        NewValues = pd.merge(AUCs, MOVs)

        # Step 5: Make a binary class
        self._df = NewValues.melt(id_vars = [self._chemical, self._concentration, self._plate, self._well], var_name = "endpoint")
        self.endpoint = "endpoint"
        self._endpoint = "endpoint"
        self.value = "value"
        self._value = "value"
