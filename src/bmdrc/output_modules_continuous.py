from .filtering import make_plate_groups
import numpy as np
import pandas as pd
from scipy import stats
from bmdrc import filtering 
import os
import json

def benchmark_dose(self, path: str):
    '''
    Calculate high level of statistics of benchmark dose fits. The Data_QC flag has is determined as follows:

    | Flag | Number of Non-Control Concentrations | Spearman Correlation | Goodness of Fit | BMD50 |
    | -- | -- | -- | -- | -- |
    | Not Fit | < 3 | < 0.2 | < 0.1 | Not within concentration range |
    | Moderate | >= 3 | 0.2 - 0.7 | >= 0.1 | Not within concentration range |
    | Good | >= 5 | > 0.7 | >= 0.1 | Within concentration range |

    Parameters
    ----------
    path
        The path to write the benchmark dose file to
    
    '''

    # Add plate groups (needed for DataQC_Flag)
    try:
        self.plate_groups
    except AttributeError:
        make_plate_groups(self)

        
    # Pull BMDS. 
    BMDS = self.bmds

    # If modeled, the data passed all filters.
    BMDS["Modeled_Flag"] = "Pass"
    BMDS.loc[BMDS["Model"].isna() | (BMDS["Model"] == "No model"), "Modeled_Flag"] = "Fail - no model fit"
    if hasattr(self, "failed_bmdl_gt_bmd10"):
        BMDS.loc[BMDS["bmdrc.Endpoint.ID"].isin(self.failed_bmdl_gt_bmd10), "Modeled_Flag"] = "Fail - poor fit (BMDL > BMD10)"

    # Add filtered data as needed
    if self.bmds_filtered is not None:

        # Pull filtered information
        BMDS_Filtered = self.bmds_filtered
        BMDS_Filtered["Modeled_Flag"] = "Fail - other filter"

        # Add where minimum concentration filter was the issue
        for row in range(len(BMDS_Filtered)):

            the_endpoint = BMDS_Filtered["bmdrc.Endpoint.ID"][row]
            the_reasons = self.plate_groups[self.plate_groups["bmdrc.Endpoint.ID"] == the_endpoint]["bmdrc.filter.reason"].unique().tolist()

            if " correlation_score_filter" in the_reasons:
                BMDS_Filtered["Modeled_Flag"][row] = "Fail - correlation score filter"

        # Remove endpoints whose models were already fit
        the_ids = BMDS["bmdrc.Endpoint.ID"].unique().tolist()
        BMDS_Filtered = BMDS_Filtered[BMDS_Filtered["bmdrc.Endpoint.ID"].isin(the_ids) == False]

        # Start final BMDS data frame 
        BMDS_Final = pd.concat([BMDS, BMDS_Filtered])

    else:
        BMDS_Final = BMDS

    # Add BMD10 and BMD50 flags
    BMDS_Final["BMD10_Flag"] = "Fail"
    BMDS_Final["BMD50_Flag"] = "Fail"
    BMDS_Final.loc[(BMDS_Final["BMD10"] >= BMDS_Final["Min_Dose"]) & (BMDS_Final["BMD10"] <= BMDS_Final["Max_Dose"]), "BMD10_Flag"] = "Pass"
    BMDS_Final.loc[(BMDS_Final["BMD50"] >= BMDS_Final["Min_Dose"]) & (BMDS_Final["BMD50"] <= BMDS_Final["Max_Dose"]), "BMD50_Flag"] = "Pass"

        ##################
    ## DATA QC FLAG ##
    ###################

    ## Add Number of Concentrations ## 
    PlateGroupsNonZero = self.plate_groups[self.plate_groups[self.concentration] != 0]

    # Get a count per concentration group
    ConcCount = PlateGroupsNonZero.loc[PlateGroupsNonZero["bmdrc.filter"] == "Keep", ["bmdrc.Endpoint.ID", self.concentration]].groupby("bmdrc.Endpoint.ID").nunique().reset_index().rename(columns = {self.concentration:"NumConc"})

    # Add to the benchmark dose table
    BMDS_Final = BMDS_Final.merge(ConcCount, left_on = "bmdrc.Endpoint.ID", right_on = "bmdrc.Endpoint.ID", how = "left")

    ## Add Spearman Correlation ## 

    CorScore = self.plate_groups

    # If the data is BinaryClass where plate and well information is available, do the following
    if hasattr(self, "value"):

        # First, only keep the values that aren't being filtered
        CorScore = CorScore.loc[CorScore["bmdrc.filter"] == "Keep", [self.concentration, "bmdrc.Endpoint.ID", "bmdrc.num.nonna", "bmdrc.num.affected"]]

        # Sum up counts
        CorScore = CorScore.groupby([self.concentration, "bmdrc.Endpoint.ID"]).sum().reset_index()

        # Calculate response
        CorScore["Response"] = CorScore["bmdrc.num.affected"] / CorScore["bmdrc.num.nonna"]

    else:

        # Calculate the response
        CorScore = CorScore.loc[CorScore["bmdrc.filter"] == "Keep", [self.concentration, "bmdrc.Endpoint.ID", self.response]].rename(columns = {self.response:"Response"})

    # Sort data.frame appropriately
    CorScore.sort_values(by = ["bmdrc.Endpoint.ID", self.concentration])

    # Calculate spearman correlations
    CorScore = CorScore[[self.concentration, "bmdrc.Endpoint.ID", "Response"]].groupby(["bmdrc.Endpoint.ID"]).corr(method = "spearman").unstack().iloc[:,1].reset_index()
    CorScore.columns = ["bmdrc.Endpoint.ID", "Spearman_Correlation"]

    # Add to the benchmark dose table
    BMDS_Final = BMDS_Final.merge(CorScore, left_on = "bmdrc.Endpoint.ID", right_on = "bmdrc.Endpoint.ID", how = "left")

    # Add Final Data QC Flag
    BMDS_Final["DataQC_Flag"] = "Not Fit"
    BMDS_Final.loc[(BMDS_Final["NumConc"] >= 3) &
                   (BMDS_Final["Spearman_Correlation"] >= 0.2), "DataQC_Flag"] = "Moderate"
    BMDS_Final.loc[(BMDS_Final["NumConc"] >= 5) & 
                   (BMDS_Final["Spearman_Correlation"] >= 0.7) &
                   (BMDS_Final["BMD50_Flag"] == "Pass"), "DataQC_Flag"] = "Good"
    
    # Fix cases where a models is not fit
    BMDS_Final.loc[BMDS_Final["Modeled_Flag"] != "Pass", "DataQC_Flag"] = "Not Fit"

    # Add columns for printing
    BMDS_Final["Chemical_ID"] = [x.split(" ")[0] for x in BMDS_Final["bmdrc.Endpoint.ID"].to_list()]
    BMDS_Final["End_Point"] = [x.split(" ")[1] for x in BMDS_Final["bmdrc.Endpoint.ID"].to_list()]

    BMDS_Final = BMDS_Final[["Chemical_ID", "End_Point", "Model", "BMD10", "BMDL", "BMD50", "AUC", "Min_Dose", "Max_Dose", "AUC_Norm", 
                             "Modeled_Flag", "DataQC_Flag", "BMD10_Flag", "BMD50_Flag", "NumConc", "Spearman_Correlation", "bmdrc.Endpoint.ID"]]
    
    # Save output table
    self.output_res_benchmark_dose = BMDS_Final

    # Write file if path is not none
    if path is not None:
        BMDS_Final.to_csv(path, header = True, index = False)

def dose_table(self, path: str):
    '''
    Calculate confidence intervals for each measured dose

    Parameters
    ----------
    path
        The path to write the dose table file to
    
    '''
        
    # Extract the specific dosages that were measured with their additional information
    dose_table = self.df[[self.chemical, self.endpoint, self.concentration, self.response]].groupby([self.chemical, self.endpoint, self.concentration]).agg(["mean", "sem", "size"]).reset_index()
    dose_table.columns = [self.chemical, self.endpoint, self.concentration, "mean", "sem", "size"]
    dose_table["dof"] = dose_table["size"] - 1
    
    # Add 95% confidence intervals
    dose_table["Low"] = np.nan
    dose_table["High"] = np.nan
    
    # Add confidence intervals
    for row in range(len(dose_table)):
        CI = stats.t.interval(0.95, dose_table["dof"][row], loc = dose_table["mean"][row], scale = dose_table["sem"][row])
        dose_table["Low"][row] = np.round(CI[0], 8) 
        dose_table["High"][row] = np.round(CI[1], 8) 
    
    # Select columns
    dose_table = dose_table[[self.chemical, self.endpoint, self.concentration, "mean", "Low", "High"]]
    
    # Rename columns
    dose_table = dose_table.rename({self.chemical: "Chemical_ID", 
                                    self.endpoint: "End_Point", 
                                    self.concentration: "Dose",
                                    "mean": "Response",
                                    "Low": "CI_Lo",
                                    "High": "CI_Hi"}, axis = 1)

    # Save output table
    self.output_res_dose_table = dose_table

    # Write file if path is not none
    if path is not None:
        dose_table.to_csv(path, header = True, index = False)
