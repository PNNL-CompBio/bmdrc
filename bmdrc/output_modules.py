import numpy as np
import pandas as pd
from astropy import stats as astrostats
from bmdrc import filtering 
import os

__author__ = "David Degnan"

def benchmark_dose(self, path):
        
    # Pull BMDS. Flag of 0 is bad, and flag of 1 is good. 
    BMDS = self.bmds
    BMDS["DataQC_Flag"] = 1

    if self.bmds_filtered is not None:
        BMDS_Filtered = self.bmds_filtered
        BMDS_Filtered["DataQC_Flag"] = 0

        # Start final BMDS data frame 
        BMDS_Final = pd.concat([BMDS, BMDS_Filtered])

    else:
        BMDS_Final = BMDS

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
    
    # Save output table
    self.output_res_benchmark_dose = BMDS_Final

    # Write file if path is not none
    if path is not None:
        BMDS_Final.to_csv(path, header = True, index = False)

def dose_table(self, path):

    # Extract the specific dosages that were measured with their additional information
    dose_table = self.plate_groups[[self.chemical, self.endpoint, self.concentration, "bmdrc.num.affected", \
                                    "bmdrc.num.nonna", "bmdrc.Endpoint.ID"]]
    
    # Add up values 
    dose_table = dose_table.groupby([self.chemical, self.endpoint, self.concentration, "bmdrc.Endpoint.ID"]).sum().reset_index()

    # Add 95% confidence intervals
    dose_table["Low"] = np.nan
    dose_table["High"] = np.nan

    # Add confidence intervals 
    for row in range(len(dose_table)):
        NumAffected = dose_table["bmdrc.num.affected"][row]
        NumNonNa = dose_table["bmdrc.num.nonna"][row]
        if NumNonNa != 0:
            CI = astrostats.binom_conf_interval(NumAffected, NumNonNa, confidence_level = 0.95)
            dose_table["Low"][row] = np.round(CI[0], 8) 
            dose_table["High"][row] = np.round(CI[1], 8) 

    # Rename columns
    dose_table = dose_table.rename({self.chemical:"Chemical_ID", self.endpoint:"End_Point", self.concentration:"Dose", \
                                    "bmdrc.num.affected":"num.affected", "bmdrc.num.nonna":"num.nonna", \
                                    "bmdrc.Endpoint.ID":"ids", "Low":"CI_Lo", "High":"CI_Hi"}, axis = 1)
    
    # Save output table
    self.output_res_dose_table = dose_table

    # Write file if path is not none
    if path is not None:
        dose_table.to_csv(path, header = True, index = False)

def report_binary(self, out_folder, report_name):

    if os.path.isdir(out_folder) == False:
        os.mkdir(out_folder)

    #####################################
    ## PULL INPUT DATA CHARACTERISTICS ##
    #####################################

    if str(type(self)) != "<class 'bmdrc.LPRClass.LPRClass'>":

        out_string = "# " + str(report_name) + "\n\n" + \
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

    else: 

        out_string = "# " + str(report_name) + "\n\n" + \
        "## Input Data\n\n" + \
        "A **lpr class** object was created." + \
        " The following column names were set:\n\n" + \
        "|Parameter|Column Name|\n" + \
        "|---------|-----------|\n" + \
        "|Chemical|" + str(self.chemical) + "|\n" + \
        "|Plate|" + str(self.plate) + "|\n" + \
        "|Well|" + str(self.well) + "|\n" + \
        "|Concentration|"  + str(self.concentration) + "|\n" + \
        "|Time|"  + str(self.time) + "|\n" + \
        "|Value|"  + str(self.value) + "|\n" + \
        "|Cycle Length|" + str(self.cycle_length) + "|\n" + \
        "|Cycle Cooldown|" + str(self.cycle_cooldown) + "|\n" + \
        "|Starting Cycle|" + str(self.starting_cycle) + "|\n\n" + \
        "## Pre-Processing\n\n#### **Combine & Make New Endpoints**\n"

    ############################
    ## PRE-PROCESSING RESULTS ##
    ############################

    ## Combine & Make New Endpoints----------------------------------------------------------------------------- 

    try:

        out_string = out_string + "New endpoints were made using existing endpoints using 'or', which means that" + \
        " if there is any endpoints with a '1', this new endpoint will also have a '1', regardless of" + \
        " how many zeroes there are in the other endpoints. See a summary table of added endpoints below:\n\n"
        
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

        out_string = out_string + "Plates with unusually high responses in negative control samples were filtered." +\
                     " The response threshold was set to **" + str(self.filter_negative_control_thresh)  + "**. See a summary below:\n\n"
        
        # Make table
        fnc_table = "|Response|Number of Plates|Filter|\n|---|---|---|\n"

        for el in range(len(self.filter_negative_control_df)):
            row = self.filter_negative_control_df.loc[el]
            fnc_table = fnc_table + "|" + str(np.round(row["Response"], 4)) + "|" + \
                        str(row["Count"]) + "|" + row["Filter"] + "|\n"

        # Save plot
        self.filter_negative_control_plot.savefig(out_folder + "/" + "filter_negative_control.png")

        out_string = out_string + fnc_table + "\nAnd here is the plot:\n![Filter Negative Control](./filter_negative_control.png)\n"
        out_string = out_string +  "\n#### **Minimum Concentration Filter**\n\n"

    except:
        out_string = out_string + "This step was not conducted.\n\n#### **Minimum Concentration Filter**\n\n"

    # Minimum Concentration Filter--------------------------------------------------------------------------------
        
    try:

        out_string = out_string + "Endpoints with too few concentration measurements (non-NA) to model are removed." +\
                     " The minimum was set to **" + str(self.filter_min_concentration_thresh)  + "**. See a summary below:\n\n"
        
        # Make table
        mc_table = "|Number of Concentrations|Number of Endpoints|Filter|\n|---|---|---|\n"

        for el in range(len(self.filter_min_concentration_df)):
            row = self.filter_min_concentration_df.loc[el]
            mc_table = mc_table + "|" + str(row["NumConc"]) + "|" + \
                        str(row["Count"]) + "|" + row["Filter"] + "|\n"

        # Save plot
        self.filter_min_concentration_plot.savefig(out_folder + "/" + "filter_minimum_concentration.png")

        out_string = out_string + mc_table + "\nAnd here is the plot:\n![Filter Minimum Concentration](./filter_minimum_concentration.png)\n"
        out_string = out_string +  "\n#### **Correlation Score Filter**\n\n"

    except:

        out_string = out_string + "This step was not conducted.\n\n#### **Correlation Score Filter**\n\n"

    # Correlation Score Filter----------------------------------------------------------------------------------
        
    try:

        out_string = out_string + "Endpoints with little to no positive correlation with dose are unexpected" +\
                     " and should be removed. The correlation threshold was set to **" + str(self.filter_correlation_score_thresh)  + "**. See a summary below:\n\n"
        
        # Correlation Score Summary Table
        the_bins = [x/100 for x in range(-100, 105, 20)]
        correlation_score = self.filter_correlation_score_df
        counts, bins = np.histogram(correlation_score["Spearman"], bins = the_bins)
        cor_score_summary = pd.DataFrame([bins, counts]).transpose().rename({0:"CorrelationScoreBin", 1:"Count"}, axis = 1).loc[0:9]

        # Make table
        cs_table = "|Correlation Score Bin|Number of Endpoints|\n|---|---|\n"

        for el in range(len(cor_score_summary)):
            row = cor_score_summary.loc[el]
            cs_table = cs_table + "|" + str(row["CorrelationScoreBin"]) + "|" + \
                       str(np.round(row["Count"], 0)) + "|\n"

        # Save plot
        self.filter_correlation_score_plot.savefig(out_folder + "/" + "filter_correlation_score.png")

        out_string = out_string + cs_table + "\nAnd here is the plot:\n![Filter Correlation Score](./filter_correlation_score.png)\n\n## Model Fitting & Output Modules\n\n#### **Filter Summary**\n\n"

    except:

        out_string = out_string + "This step was not conducted.\n\n## Model Fitting & Output Modules\n\n#### **Filter Summary**\n\n"
        
    ###################
    ## MODEL FITTING ##
    ###################
        
    try:
    
        # Filter Summary-------------------------------------------------------------------------------------
            
        # Get removal and kept counts 
        removed = len(pd.unique(self.plate_groups[self.plate_groups["bmdrc.filter"] == "Remove"]["bmdrc.Endpoint.ID"]))
        kept = len(pd.unique(self.plate_groups[self.plate_groups["bmdrc.filter"] == "Keep"]["bmdrc.Endpoint.ID"]))       
        total = removed + kept

        out_string = out_string + "Overall, " + str(total) + " endpoint and chemical combinations were considered. " + str(kept) + \
                    " were deemed eligible for modeling, and " + str(removed) + " were not based on filtering selections explained" + \
                    " in the previous section.\n\n#### **Model Fitting Selections**\n\n"
        
        # Model Fitting Selections------------------------------------------------------------------------------

        out_string = out_string + "The following model fitting parameters were selected.\n\n|Parameter|Value|Parameter Description|\n" + \
                    "|---|---|---|\n|Goodness of Fit Threshold|" + str(self.model_fitting_gof_threshold) + "|Minimum p-value for fitting a model. Default is 0.1|\n" + \
                    "|Akaike Information Criterion (AIC) Threshold|" + str(self.model_fitting_aic_threshold) + "|Any models with an AIC within this value are considered" + \
                    " an equitable fit. Default is 2.\n" + \
                    "|Model Selection|" + self.model_fitting_model_selection + "|Either return one model with the lowest BMDL, or combine equivalent fits|\n\n#### **Model Quality Summary**\n\n"
        
        # Model Quality Summary---------------------------------------------------------------------------------

        out_string = out_string + "Below is a summary table of the number of endpoints with a high quality fit (a flag of 1, meaning that" + \
                    " the BMD10 value is within the range of measured doses) and those that are not high quality (a flag of 0).\n\n"

        # Make dataframe of flag counts, adding any missing flags
        dataqc_table = self.bmds[["Model", "DataQC_Flag"]].groupby("DataQC_Flag").count().reset_index().rename({"Model":"Count"}, axis = 1)
        if ((0 in dataqc_table["DataQC_Flag"].values.tolist()) == False):
            dataqc_table = pd.concat([dataqc_table, pd.DataFrame({"DataQC_Flag":[0], "Count":[0]})])
        if ((1 in dataqc_table["DataQC_Flag"].values.tolist()) == False):
            dataqc_table = pd.concat([dataqc_table, pd.DataFrame({"DataQC_Flag":[1], "Count":[0]})])

        # Add flag counts 
        out_string = out_string + "|Flag|Count|\n|---|---|\n|0|" + str(dataqc_table[dataqc_table["DataQC_Flag"] == 0]["Count"].values[0]) + "|\n" + \
                    "|1|" + str(dataqc_table[dataqc_table["DataQC_Flag"] == 1]["Count"].values[0]) + "|\n\n#### **Output Modules**\n\nBelow, see a table of" + \
                    " useful methods for extracting outputs from bmdrc.\n\n"
        
        # Add useful parameters 
        out_string = out_string + "|Method|Description|\n|---|---|\n|.bmds|Table of fitted benchmark dose values|\n" + \
                    "|.bmds_filtered|Table of filtered models not eligible for benchmark dose calculations|\n" + \
                    "|.output_res_benchmark_dose|Table of benchmark doses for all models, regardless of whether they were filtered or not|\n" + \
                    "|.p_value_df|Table of goodness of fit p-values for every eligible endpoint|\n" + \
                    "|.aic_df|Table of Akaike Information Criterion values for every eligible endpoint|\n" + \
                    "|.response_curve|Plot a benchmark dose curve for an endpoint|\n\n"
        
    except:

        out_string = out_string + "Model fits were not conducted."
        
    file = open(out_folder + "/" + report_name + ".md", "w")
    file.write(out_string)
    file.close()