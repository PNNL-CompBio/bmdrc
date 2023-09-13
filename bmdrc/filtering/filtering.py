import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

__author__ = "David Degnan"

def make_plate_groups(self):
    '''
    Support function for the filter modules. 
    Assign groups based on the chemical, concentration, plate, and endpoint.
    This step should be completed after all pre-processing steps are finished. 
    '''

    # If there is a bmdrc well id, drop it
    if "bmdrc.Well.ID" in self.df.columns.tolist():
        self.plate_groups = self.df.drop([self.well, "bmdrc.Well.ID"], 1).groupby(by = [self.chemical, self.concentration, self.plate, self.endpoint], as_index = False)
    else:
        self.plate_groups = self.df.drop([self.well], 1).groupby(by = [self.chemical, self.concentration, self.plate, self.endpoint], as_index = False)
        
    # Get the number of samples per group
    num_tot_samples = self.plate_groups.size().rename(columns = {"size": "bmdrc.num.tot"})

    # Get the number of non-na samples per groups
    num_nonna = self.plate_groups.count().rename(columns = {"value": "bmdrc.num.nonna"})

    # Get the number affected
    num_affected = self.plate_groups.sum().rename(columns = {"value": "bmdrc.num.affected"})

    # Merge to create missingness dataframe
    self.plate_groups = pd.merge(pd.merge(num_tot_samples, num_nonna), num_affected)

    # Create IDs of chemical.id, plate.id, and endpoint in plate_groups 
    self.plate_groups["bmdrc.Plate.ID"] = self.plate_groups[self.chemical].astype(str) + " " + \
                                            self.plate_groups[self.plate].astype(str) + " " + \
                                            self.plate_groups[self.endpoint].astype(str)
    
def check_apply_diagnostic(apply, diagnostic):
    '''
    Support function for the filter modules. 
    Check that apply and diagnostic are acceptable inputs.
    '''

    ##################
    ## CHECK INPUTS ##
    ##################

    # Check that apply is either True or False
    if type(apply) != bool:
        raise Exception("apply must be a True or False.")
    
    # Check that diagnostic is either plot, df, or both
    if isinstance(diagnostic, list):
        if all(diagnostic in ["plot", "df"]) == False:
            raise Exception("diagnostic must be either 'plot', 'df', or a list with 'plot' and 'df'")
    else:
        if diagnostic != "plot" and diagnostic != "df":
            raise Exception("diagnostic must be either 'plot', 'df', or a list with 'plot' and 'df'")
        
def negative_control_plot(neg_control_df):
    '''
    Support function for the filter modules. 
    Return the negative control diagnostic plot. 
    '''

    fig = plt.figure(figsize = (10, 5))

    colors = {'Keep':'steelblue', 'Filter':'firebrick'}
    color_choices = neg_control_df["Filter"].apply(lambda x: colors[x])
    labels = list(colors.keys())
    handles = [plt.Rectangle((0,0),1,1, color=colors[label]) for label in labels]

    plt.bar(x = [x for x in range(len(neg_control_df))], height = neg_control_df["Count"], 
            edgecolor = "black", tick_label = np.round(neg_control_df["Response"], 4),
            color = color_choices, label = colors)
    plt.title("Counts of proportional responses per plate, endpoint, and chemical group")
    plt.xlabel("Proportional response in negative controls")
    plt.ylabel("Count")
    plt.legend(handles, labels)

    return(fig)

def negative_control(self, percentage, apply, diagnostic):
    '''
    Filter to remove plates with unusually high expression in the controls. 

    percentage: (float) A value between 0 and 100 indicating the percentage of phenotypic
    expression in the controls that is permissable. Default is 50. 

    apply: (logical) Apply the filter. Default is False. 

    diagnostic: (list - string) If apply if False, see a diagnostic plot with "plot" or 
    a diagnostics data.frame with "df"
    '''

    ##################
    ## CHECK INPUTS ##
    ##################

    # Assert that percentage is numeric
    try:
        percentage = float(percentage)
    except ValueError:
        raise Exception("percentage must be a float.")
    
    # Check that percentage is in the right range 
    if percentage < 0 or percentage > 100:
        raise Exception("percentage must be between 0 and 100.")
    
    # Let the user know they picked a small value, but that it is acceptable. 
    if (percentage < 1):
        print("percentage should range from 0-100, as in 0-100%. This value will be a very small percentage.")
    
    # Check apply and diagnostic
    check_apply_diagnostic(apply, diagnostic)

    ##############################
    ## MAKE GROUPS IF NECESSARY ##
    ##############################

    try:
        self.plate_groups
    except AttributeError:
        make_plate_groups(self)
    
    ###############################
    ## CREATE DIAGNOSTIC SUMMARY ##
    ###############################

    # Extract negative controls 
    NegControls = self.plate_groups[self.plate_groups[self.concentration] == 0]

    # Calculate responses in negative controls
    NegControlRes = pd.DataFrame((NegControls["bmdrc.num.affected"] / NegControls["bmdrc.num.nonna"])).value_counts().rename_axis("Response").reset_index().rename(columns = {0:"Count"}).sort_values(by = ["Response"]).reset_index(drop=True)
    
    # Set all values to be kept 
    NegControlRes["Filter"] = "Keep"

    # Determine what values will be removed
    NegControlRes.loc[NegControlRes["Response"] > percentage/100, "Filter"] = "Filter"

    #######################
    ## RETURN DIAGNOSTIC ##
    #######################

    if isinstance(diagnostic, list) == False: 
        if diagnostic == "df":
            self.filter_negative_control_df = NegControlRes
        if diagnostic == "plot":
            self.filter_negative_control_plot = negative_control_plot(NegControlRes)



  


    