import numpy as np
import pandas as pd
from astropy import stats as astrostats
import os

__author__ = "David Degnan"

def benchmark_dose(self, path):
    
    def BMD_Range_Flag(id, BMD):
        
        # Get concentrations
        concs = self.plate_groups[self.plate_groups["bmdrc.Endpoint.ID"] == id][self.concentration]

        if (np.isnan(BMD)):
            return(np.nan)
        elif (BMD <= max(concs) and BMD >= min(concs)):
            return(0)
        else:
            return(1)
        
    # Pull BMDS. Flag of 0 is bad, and flag of 1 is good. 
    BMDS = self.bmds
    BMDS["DataQC_Flag"] = 1
    BMDS_Filtered = self.bmds_filtered
    BMDS_Filtered["DataQC_Flag"] = 0

    # Start final BMDS data frame 
    BMDS_Final = pd.concat([BMDS, BMDS_Filtered])

    # Add BMD10 and BMD50 flags
    BMDS_Final["BMD10_Flag"] = 0
    BMDS_Final["BMD50_Flag"] = 0
    BMDS_Final.loc[(BMDS_Final["BMD10"] >= BMDS_Final["Min_Dose"]) & (BMDS_Final["BMD10"] <= BMDS_Final["Max_Dose"]), "BMD10_Flag"] = 1
    BMDS_Final.loc[(BMDS_Final["BMD50"] >= BMDS_Final["Min_Dose"]) & (BMDS_Final["BMD50"] <= BMDS_Final["Max_Dose"]), "BMD50_Flag"] = 1

    # Add BMD Analysis Flag
    BMDS_Final["BMD_Analysis_Flag"] = BMDS_Final["BMD10_Flag"] + BMDS_Final["BMD50_Flag"]
    
    # Add columns for printing
    BMDS_Final["Chemical_ID"] = [x.split(" ")[0] for x in BMDS_Final["bmdrc.Endpoint.ID"].to_list()]
    BMDS_Final["End_Point"] = [x.split(" ")[1] for x in BMDS_Final["bmdrc.Endpoint.ID"].to_list()]

    BMDS_Final = BMDS_Final[["Chemical_ID", "End_Point", "Model", "BMD10", "BMDL", "BMD50", "AUC", "Min_Dose", "Max_Dose", "AUC_Norm", 
                "DataQC_Flag", "BMD_Analysis_Flag", "BMD10_Flag", "BMD50_Flag", "bmdrc.Endpoint.ID"]]
    
    # Arrange by analysis flag
    BMDS_Final = BMDS_Final.sort_values("BMD_Analysis_Flag", ascending = False)
    
    self.output_res_benchmark_dose = BMDS_Final

def report(self, outpath, title = "Benchmark Dose Response Curves", curve_plots = False):
    command = "quarto render /Users/degn400/Git_Repos/bmdrc/bmdrc/output_modules/bmdrc.qmd -o outpath -P thedata:self -P titleName:title -P curve_plots:curvePlots"
    os.system(command)