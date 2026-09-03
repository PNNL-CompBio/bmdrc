import operator
from abc import abstractmethod

import numpy as np
import pandas as pd

from .BinaryClass import BinaryClass, DataClass

class LPRClass(DataClass):
    """
    Generates a bmdrc object from larval photomotor response data, which must be in long format.

    Parameters
    ----------
    df
        A pandas dataframe containing columns with the chemical, concentration, plate, well, time, and value information.
    chemical
        A string indicating the name of the column containing the chemical IDs, which should be strings
    plate
        A string indicating the name of the column indicating the plate IDs, which should be strings
    well
        A string indicating the name of the column with the well IDs, which should be strings
    concentration
        A string indicating the name of the column containing the concentrations, which should be numerics
    time
        A string indicating the name of the column containing time, which should be a string or integer. Strings should contain a number.
    value
        A string indicating the name of the column containing the binary values, which should be 0 for absent, and 1 for present. Not used if the light photomotor response
    cycle_length
        A numeric for the length of a light or dark cycle. Default is 20. The unit is a 6-second measure, so 20 six second measures is 2 minutes.
    cycle_cooldown
        A numeric for the length of time between cycles. Default is 10. The unit is a 6-second measure, so 10 six second measures is 1 minute.
    starting_cycle
        A string of either the "light" or "dark" cycle depending on whether the first measurement was a light or dark cycle. Default is "light".
    """

    # Define the input checking functions. Include raw and transformed data.frames
    def __init__(
        self,
        df,
        chemical,
        plate,
        well,
        concentration,
        time,
        value,
        cycle_length=20.0,
        cycle_cooldown=10.0,
        starting_cycle="light",
    ):
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
    df = property(operator.attrgetter("_df"))
    chemical = property(operator.attrgetter("_chemical"))
    plate = property(operator.attrgetter("_plate"))
    well = property(operator.attrgetter("_well"))
    concentration = property(operator.attrgetter("_concentration"))
    time = property(operator.attrgetter("_time"))
    value = property(operator.attrgetter("_value"))
    cycle_length = property(operator.attrgetter("_cycle_length"))
    cycle_cooldown = property(operator.attrgetter("_cycle_cooldown"))
    starting_cycle = property(operator.attrgetter("_starting_cycle"))
    cycles = property(operator.attrgetter("_cycle"))

    unacceptable = [
        "bmdrc.Well.ID",
        "bmdrc.num.tot",
        "bmdrc.num.nonna",
        "bmdrc.num.affected",
        "bmdrc.Plate.ID",
        "bmdrc.Endpoint.ID",
        "bmdrc.filter",
        "bmdrc.filter.reason",
        "bmdrc.frac.affected",
        "cycle",
    ]

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
            raise Exception(
                chemicalname + " is not a permitted name. Please rename this column."
            )
        self._df[chemicalname] = self._df[chemicalname].astype(str)
        self._chemical = chemicalname

    @plate.setter
    def plate(self, platename):
        if not isinstance(platename, str):
            raise Exception("plate must be a name of a column in df.")
        if not platename in self._df.columns:
            raise Exception(platename + " is not in the column names of df.")
        if platename in self.unacceptable:
            raise Exception(
                platename + " is not a permitted name. Please rename this column."
            )
        self._plate = platename

    @well.setter
    def well(self, wellname):
        if not isinstance(wellname, str):
            raise Exception("well must be a name of a column in df.")
        if not wellname in self._df.columns:
            raise Exception(wellname + " is not in the column names of df.")
        if wellname in self.unacceptable:
            raise Exception(
                wellname + " is not a permitted name. Please rename this column."
            )
        self._well = wellname

    @concentration.setter
    def concentration(self, concentrationname):
        if not isinstance(concentrationname, str):
            raise Exception("concentration must be a name of a column in df.")
        if not concentrationname in self._df.columns:
            raise Exception(concentrationname + " is not in the column names of df.")
        if concentrationname in self.unacceptable:
            raise Exception(
                concentrationname
                + " is not a permitted name. Please rename this column."
            )
        self._df[concentrationname] = pd.to_numeric(self._df[concentrationname])
        self._concentration = concentrationname

    @time.setter
    def time(self, timename):
        if not isinstance(timename, str):
            raise Exception("time must be a name of a column in df.")
        if not timename in self._df.columns:
            raise Exception(timename + " is not in the column names of df.")
        if timename in self.unacceptable:
            raise Exception(
                timename + " is not a permitted name. Please rename this column."
            )
        self._df[timename] = (
            self._df[timename].str.extract("(\d+)", expand=False).astype(float)
        )
        self._time = timename

    @value.setter
    def value(self, valuename):
        if not isinstance(valuename, str):
            raise Exception("value must be a name of a column in df.")
        if not valuename in self._df.columns:
            raise Exception(valuename + " is not in the column names of df.")
        if valuename in self.unacceptable:
            raise Exception(
                valuename + " is not a permitted name. Please rename this column."
            )
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
        if not starting_cycle in ["light", "dark"]:
            raise Exception("starting_cycle must be either 'light' or 'dark'.")
        self._starting_cycle = starting_cycle

    # LPR-specific function: determines light and dark cycles
    def add_cycles(self):
        """Specific LPR function that adds cycle information following users setting cycle_time,
        cycle_cooldown, and samples_to_remove"""

        print("...defining cycles")

        # Unique and arrange times
        cycle_info = pd.DataFrame(self._df[self._time].unique()).rename(
            {0: self._time}, axis=1
        )
        cycle_info[self._time] = cycle_info[self._time].astype(float)
        cycle_info = cycle_info.sort_values(by=[self._time])

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
            if first_count < self._cycle_length:
                cycle_order.append(self._starting_cycle + str(cycle_count))
                first_count += 1
            elif gap_a_count < self._cycle_cooldown:
                cycle_order.append("gap_" + self._starting_cycle + str(cycle_count))
                gap_a_count += 1
            elif second_count < self._cycle_length:
                cycle_order.append(other_cycle + str(cycle_count))
                second_count += 1
            elif gap_b_count < self._cycle_cooldown:
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

        return self._df.merge(cycle_info)

    def to_dichotomous(self, the_df, the_value):
        """
        Converts continuous AUC or MOV values to dichotomous (0/1) classification.

        Per the paper methodology:
        - Thresholds are computed from ONLY the normally responding (positive-valued) control fish.
        - Control fish are classified as abnormal ONLY if their value is negative (hypoactive).
        - Chemical-exposed fish are classified as abnormal if:
            (a) their value is negative (hypoactive), OR
            (b) their positive value is an outlier relative to the normal control distribution
                (above Q3 + 1.5*IQR or below Q1 - 1.5*IQR using Tukey's method).

        Returns a numpy array of 0/1 values aligned positionally with the input DataFrame.
        """

        working = the_df[[self._chemical, self._plate, self._concentration, the_value]].copy()

        # --------------------------------------------------------------------------
        # Step 1: Compute outlier thresholds from ONLY normally responding controls
        #         (positive endpoint values in the control group)
        # --------------------------------------------------------------------------
        positive_controls = working[
            (working[self._concentration] == 0) & (working[the_value] > 0)
        ]

        # If there are no positive controls for a chemical-plate, thresholds will be NaN
        # (handled gracefully by the left merge below)
        thresholds = (
            positive_controls.groupby([self._chemical, self._plate])[the_value]
            .agg(Q1=lambda x: x.quantile(0.25), Q3=lambda x: x.quantile(0.75))
            .reset_index()
        )
        thresholds["IQR"] = thresholds["Q3"] - thresholds["Q1"]
        thresholds["Low"] = thresholds["Q1"] - (1.5 * thresholds["IQR"])
        thresholds["High"] = thresholds["Q3"] + (1.5 * thresholds["IQR"])
        thresholds = thresholds[[self._chemical, self._plate, "Low", "High"]]

        # --------------------------------------------------------------------------
        # Step 2: Left merge to attach thresholds to all rows.
        #         Plates with no positive controls get NaN thresholds (only hypoactivity
        #         can be detected on those plates). Left merge preserves row order.
        # --------------------------------------------------------------------------
        merged = working.merge(
            thresholds, on=[self._chemical, self._plate], how="left"
        )

        # --------------------------------------------------------------------------
        # Step 3: Apply DIFFERENT classification rules for controls vs. exposed fish
        # --------------------------------------------------------------------------
        is_control = merged[self._concentration] == 0
        values = merged[the_value]
        low = merged["Low"]
        high = merged["High"]

        # Controls: abnormal ONLY if hypoactive (negative value)
        control_abnormal = values < 0

        # Chemical-exposed: abnormal if hypoactive OR positive outlier
        exposed_abnormal = (
            (values < 0)
            | (values > high)
            | ((values > 0) & (values < low))
        )

        # Combine: use control rule for controls, exposed rule for chemical groups
        result = np.where(is_control, control_abnormal, exposed_abnormal)

        # Return as integer array (0 = normal, 1 = abnormal) aligned positionally
        return result.astype(int)

    # LPR-specific function: Calculate AUC values
    def calculate_aucs(self, cycles):
        """Specific LPR function for calculating AUC values"""

        print("...calculating AUC values")

        # Remove gaps from the AUC calculation (only sum light and dark periods)
        aucs = cycles[~cycles["cycle"].str.contains("gap")].drop(
            labels=self._time, axis=1
        )

        # Sum movement values within each cycle phase for each fish
        aucs = (
            aucs.groupby(
                by=[
                    self._chemical,
                    self._concentration,
                    self._plate,
                    self._well,
                    "cycle",
                ]
            )
            .sum()
            .reset_index()
        )

        # Initiate list to store all values
        store_aucs = []

        # Iterate through all cycles, subtracting light AUC from dark AUC
        for cycle_num in range(self._max_cycle):

            cycle_num = cycle_num + 1
            light_name = "light" + str(cycle_num)
            dark_name = "dark" + str(cycle_num)

            # Merge light and dark sums for each fish
            to_calc_auc = pd.merge(
                aucs[aucs["cycle"] == light_name]
                .rename(columns={"value": "light"})
                .drop("cycle", axis=1),
                aucs[aucs["cycle"] == dark_name]
                .rename(columns={"value": "dark"})
                .drop("cycle", axis=1),
                how="left",
            )

            # AUC = dark_sum - light_sum (positive = more active in dark = normal)
            to_calc_auc["Cycle"] = "AUC" + str(cycle_num)
            to_calc_auc["AUC"] = to_calc_auc["dark"] - to_calc_auc["light"]
            store_aucs.append(to_calc_auc)

        # Concatenate all cycles
        auc_values = pd.concat(store_aucs).drop_duplicates().dropna()

        # Pivot so each cycle is its own column (AUC1, AUC2, etc.)
        auc_process = (
            auc_values[
                [self._chemical, self._plate, self._concentration, self._well, "Cycle", "AUC"]
            ]
            .pivot(
                index=[self._chemical, self._plate, self._concentration, self._well],
                columns="Cycle",
                values="AUC",
            )
            .reset_index()
        )

        # Convert each AUC column from continuous to dichotomous (0/1)
        for x in range(self._max_cycle):
            value = "AUC" + str(x + 1)
            auc_process[value] = self.to_dichotomous(auc_process, value)

        return auc_process

    # LPR-specific function: Calculate MOV values
    def calculate_movs(self, cycles):
        """Specific LPR function for calculating MOV values"""

        print("...calculating MOV values")

        # Full cycle length (light + gap + dark + gap)
        full_cycle = (self._cycle_length * 2) + (self._cycle_cooldown * 2)

        # --------------------------------------------------------------------------
        # Determine the exact time indices for the MOV calculation.
        # MOV = movement at first dark time point - movement at last light time point
        #
        # For starting_cycle == "light":
        #   Structure: light(0 to CL-1), gap(CL to CL+CC-1), dark(CL+CC to 2*CL+CC-1), gap(...)
        #   Last light time = (CL - 1) + x * full_cycle
        #   First dark time = (CL + CC) + x * full_cycle
        #
        # For starting_cycle == "dark":
        #   Structure: dark(0 to CL-1), gap(CL to CL+CC-1), light(CL+CC to 2*CL+CC-1), gap(...)
        #   Light-to-dark transition spans across cycles:
        #   Last light time = (2*CL + CC - 1) + x * full_cycle
        #   First dark time = full_cycle + x * full_cycle = (x+1) * full_cycle
        # --------------------------------------------------------------------------
        if self._starting_cycle == "light":
            self._last_light_times = [
                (self._cycle_length - 1) + (x * full_cycle)
                for x in range(self._max_cycle)
            ]
            self._first_dark_times = [
                (self._cycle_length + self._cycle_cooldown) + (x * full_cycle)
                for x in range(self._max_cycle)
            ]
        else:
            # Dark starts first: dark, gap, light, gap, dark, gap, light, gap, ...
            # Light-to-dark transitions cross cycle boundaries
            self._last_light_times = [
                (2 * self._cycle_length + self._cycle_cooldown - 1) + (x * full_cycle)
                for x in range(self._max_cycle)
            ]
            self._first_dark_times = [
                full_cycle + (x * full_cycle)
                for x in range(self._max_cycle)
            ]

        # Select only the rows at our target time points
        movs = cycles[
            (cycles[self._time].isin(self._last_light_times))
            | (cycles[self._time].isin(self._first_dark_times))
        ]

        # Initiate list to store all calculated values
        store_movs = []

        # Iterate through all cycles, computing dark - light at transition
        for x in range(self._max_cycle):

            light_rows = movs[movs[self._time] == self._last_light_times[x]].rename(
                columns={"value": "light"}
            ).drop([self._time, "cycle"], axis=1)

            dark_rows = movs[movs[self._time] == self._first_dark_times[x]].rename(
                columns={"value": "dark"}
            ).drop([self._time, "cycle"], axis=1)

            # Inner merge: only fish with both measurements
            to_calc_mov = pd.merge(light_rows, dark_rows)

            # MOV = dark_value - light_value (positive = increased activity at transition = normal)
            to_calc_mov["Cycle"] = "MOV" + str(x + 1)
            to_calc_mov["MOV"] = to_calc_mov["dark"] - to_calc_mov["light"]
            store_movs.append(to_calc_mov)

        # Concatenate all cycles
        mov_values = pd.concat(store_movs).drop_duplicates().dropna()

        # Pivot so each cycle is its own column (MOV1, MOV2, etc.)
        mov_process = (
            mov_values[
                [self._chemical, self._plate, self._concentration, self._well, "Cycle", "MOV"]
            ]
            .drop_duplicates()
            .reset_index(drop = True)
            .pivot(
                index=[self._chemical, self._plate, self._concentration, self._well],
                columns="Cycle",
                values="MOV",
            )
            .reset_index()
        )

        # Convert each MOV column from continuous to dichotomous (0/1)
        for x in range(self._max_cycle):
            value = "MOV" + str(x + 1)
            try:
                mov_process[value] = self.to_dichotomous(mov_process, value)
            except KeyError:
                pass

        return mov_process

    # LPR-specific: Convert LPR continuous to Dichotomous
    def convert_LPR(self):
        """Wrapper function for all LPR-specific functions for calculating cycles,
        AUC values, MOV values, and converting them to dichotomous values"""

        id_vars = [self._chemical, self._concentration, self._plate, self._well]

        # Step 1: Make Cycle Information
        CycleInfo = self.add_cycles()

        # Step 2: Calculate AUC values
        AUCs = self.calculate_aucs(CycleInfo)
        AUCs = AUCs.dropna(subset=id_vars)

        # Step 3: Calculate MOV values
        MOVs = self.calculate_movs(CycleInfo)
        MOVs = MOVs.dropna(subset=id_vars)

        # Step 4: Merge the results
        NewValues = pd.merge(AUCs, MOVs)

        # Step 5: Make a binary class
        self._df = NewValues.melt(
            id_vars=id_vars,
            var_name="endpoint",
        )
        self.endpoint = "endpoint"
        self._endpoint = "endpoint"
        self.value = "value"
        self._value = "value"