import numpy as np
import pandas as pd
from scipy import stats
from bmdrc import filtering 
import os
import json

def benchmark_dose(self, path: str):
    '''
    Calculate high level of statistics of benchmark dose fits

    Parameters
    ----------
    path
        The path to write the benchmark dose file to
    
    '''
        
    # Pull BMDS. 
    BMDS = self.bmds

    # If modeled, the data passed all filters.
    BMDS["DataQC_Flag"] = "Pass"

    # Add filtered data as needed
    if self.bmds_filtered is not None:

        # Pull filtered information
        BMDS_Filtered = self.bmds_filtered
        BMDS_Filtered["DataQC_Flag"] = "Fail - other filter"

        # Add where minimum concentration filter was the issue
        for row in range(len(BMDS_Filtered)):

            the_endpoint = BMDS_Filtered["bmdrc.Endpoint.ID"][row]
            the_reasons = self.plate_groups[self.plate_groups["bmdrc.Endpoint.ID"] == the_endpoint]["bmdrc.filter.reason"].unique().tolist()

            if " correlation_score_filter" in the_reasons:
                BMDS_Filtered["DataQC_Flag"][row] = "Fail - correlation score filter"

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

    # Add BMD Analysis Flag
    BMDS_Final["BMD_Analysis_Flag"] = BMDS_Final.apply(
        lambda x: "Pass" if x["BMD10_Flag"] == "Pass" and x["BMD50_Flag"] == "Pass" else "Fail", axis=1
    )
    
    # Add columns for printing
    BMDS_Final["Chemical_ID"] = [x.split(" ")[0] for x in BMDS_Final["bmdrc.Endpoint.ID"].to_list()]
    BMDS_Final["End_Point"] = [x.split(" ")[1] for x in BMDS_Final["bmdrc.Endpoint.ID"].to_list()]

    BMDS_Final = BMDS_Final[["Chemical_ID", "End_Point", "Model", "BMD10", "BMDL", "BMD50", "AUC", "Min_Dose", "Max_Dose", "AUC_Norm", 
                             "DataQC_Flag", "BMD_Analysis_Flag", "BMD10_Flag", "BMD50_Flag", "bmdrc.Endpoint.ID"]]
    
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
