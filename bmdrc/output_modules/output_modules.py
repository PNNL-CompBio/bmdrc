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

def report(self, out_path):

    #####################################
    ## PULL INPUT DATA CHARACTERISTICS ##
    #####################################

    out_string = "# Benchmark Dose Curves\n\n" + \
    "## Input Data\n\n" + \
    "A **binary class** object was created using data in **" + str(self.format) + "** format." + \
    " The following column names were set:\n\n" + \
    "|Parameter|Column Name|\n" + \
    "|---------|-----------|\n" + \
    "|Chemical|" + str(self.chemical) + "|\n" + \
    "|Plate|" + str(self.plate) + "|\n" + \
    "|Well|" + str(self.well) + "|\n" + \
    "|Concentration|"  + str(self.concentration) + "|\n" + \
    "|Endpoint|"  + str(self.endpoint) + "|\n" + \
    "|Value|"  + str(self.value) + "|\n\n" + \
    "## Pre-Processing\n\n#### **Combine & Make New Endpoints**\n"

    ############################
    ## PRE-PROCESSING RESULTS ##
    ############################

    ## Combine & Make New Endpoints----------------------------------------------------------------------------- 

    try:

        # Trigger except if object does not exist
        test_combine_endpoints = self.report_combination

        out_string = out_string + "New endpoints were made using existing endpoints. See a summary below:\n\n"
        
        the_combined = "|New Endpoint Name|Combined Existing Endpoints|\n|---|---|\n"
        for key in self.report_combination:
            value_collapse = "|"
            for val in self.report_combination[key]:
                value_collapse = value_collapse + val + ", "
            the_combined = the_combined + "|" + key + value_collapse[0:(len(value_collapse)-2)] + "|\n"
        the_combined = the_combined + "\n"
        out_string = out_string + the_combined + "#### **Set Invalid Wells to NA**\n\n"

    except:
        out_string = out_string + "This step was not conducted.\n\n#### **Set Invalid Wells to NA**\n\n"

    # Set Invalid Wells to NA-----------------------------------------------------------------------------------
        
    try:
        
        # Trigger except if the object does not exist
        test_report_well_na = self.report_well_na

        out_string = out_string + "In some cases, like when a sample fish dies, many affected endpoints" + \
                     " need to be set to NA. Here, the 'Endpoint Name' column denotes the specific" + \
                     " endpoint that sets this rule. In this example, it could be MORT for mortality." + \
                     " Then, the endpoint value needs to be set, which in this case would be a 1 to" + \
                     " indicate sample fish that did die. All endpoints would then be set to NA except" + \
                     " for cases where the endpoint should not be affected, which are referred to as" + \
                     " 'Endpoint Exceptions.'\n\n"
        
        endpoints = "|Endpoint Name|Endpoint Value|Endpoint Exceptions|\n|---|---|---|\n"

        for el in range(len(self.report_well_na)):

            embedded_list = self.report_well_na[el]
            endpoints = endpoints + "|"

            for the_endpoint in embedded_list[0]:
                endpoints = endpoints + the_endpoint + ", "
            endpoints = endpoints[0:(len(endpoints)-2)] + "|"

            for the_value in embedded_list[1]:
                endpoints = endpoints + str(the_value) + ", "
            endpoints = endpoints[0:(len(endpoints)-2)] + "|"

            if embedded_list[2] is not None:
                for the_exception in embedded_list[2]:
                    endpoints = endpoints + str(the_exception) + ", "
                endpoints = endpoints[0:(len(endpoints)-2)] + "|\n"
            else:
                endpoints = endpoints + "None|\n"
            
        endpoints = endpoints + "\n#### **Remove Invalid Endpoints**\n\n"

        out_string = out_string + endpoints

    except:
        out_string = out_string + "This step was not conducted.\n\n#### **Remove Invalid Endpoints**\n\n"

    # Remove Invalid Endpoints---------------------------------------------------------------------------------
        
    try:

        # Trigger except if the object does not exist
        test_report_endpoint_removal = self.report_endpoint_removal

        the_removed = ""
        for removed in self.report_endpoint_removal:
            the_removed = the_removed + removed + ", "
        the_removed = the_removed[0:(len(the_removed)-2)]

        out_string = out_string + "The following endpoints were removed: " + the_removed + "\n\n"

    except:

        out_string = out_string + "This step was not conducted.\n\n"

    out_string = out_string + "## Filtering\n\n#### **Negative Control Filter**\n\n"

    #######################
    ## FILTERING RESULTS ##
    #######################

    # Negative Control Filter------------------------------------------------------------------------------------

    try:

        # Trigger except if object does not exist
        test_fnc = self.filter_negative_control_df

        out_string = out_string + "Unusually high responses in negative control samples were filtered." +\
                     " The response threshold was set to **" + str(self.filter_negative_control_percentage)  + "**. See a summary below:\n\n"
        
        fnc_table = "|Response|Count|Filter|\n|---|---|---|\n"

        for el in range(len(self.filter_negative_control_df)):
            row = self.filter_negative_control_df.loc[el]
            fnc_table = fnc_table + "|" + str(np.round(row["Response"], 4)) + "|" + \
                        str(row["Count"]) + "|" + row["Filter"] + "|\n"

        out_string = out_string + fnc_table + "\n#### **Minimum Concentration Filter**\n\n"

    except:
        out_string = out_string + "This step was not conducted.\n\n#### **Minimum Concentration Filter**\n\n"

    file = open(out_path, "w")
    file.write(out_string)
    file.close()