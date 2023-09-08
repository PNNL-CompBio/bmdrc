import pandas as pd
import numpy as np

__author__ = "David Degnan"

def well_to_na(self, endpoint_name, endpoint_value, except_endpoint = None):
    '''
    Remove any wells where a specific endpoint has a specific value. 
    Wells are set to NA.

    endpoint_name: (list - string) list of endpoints to remove 

    endpoint_value: (numeric) specific value or list of values that the endpoint needs to have to remove the well

    except_endpoint: (list - string) list of endpoints that should not have their wells affected
    '''

    ############################
    ## CHECK INPUT PARAMETERS ##
    ############################

    # Iterate through each endpoint to confirm it is a valid choice 
    for endpoint in endpoint_name:
        if endpoint in self.df[self.endpoint].unique() is False:
            raise Exception(endpoint + " is not in an endpoint in the DataClass object")
        
    # Convert endpoint value to a list
    if isinstance(endpoint_value, list) == False:
        endpoint_value = list(endpoint_value)

    # Confirm the value exists for at least one endpoint
    if any(self.df[self.df[self.endpoint].isin(endpoint_name)][self.value].isin(endpoint_value)) == False:
        raise Exception("None of the provided endpoint_names contain the provided endpoint_value")
    
    # Iterate through each except endpoint to confirm they are valid choices
    if except_endpoint is not None:
        for endpoint in except_endpoint:
            if endpoint in self.df[self.endpoint].unique() is False:
                raise Exception(endpoint + " is not in an endpoint in the DataClass object")
            
    ####################################
    ## SET WELLS TO NA TO REMOVE THEM ##
    ####################################

    # If it hasn't been made yet, make the bmdrc.Well.ID column
    if "bmdrc.Well.ID" not in self.df.columns.tolist():
        self.df["bmdrc.Well.ID"] = self.df[self.chemical].astype(str) + " " + self.df[self.concentration].astype(str) + " " + self.df[self.plate].astype(str) + " " + self.df[self.well].astype(str)
    
    # List wells to remove
    wells_rm = self.df[(self.df[self.endpoint].isin(endpoint_name)) & (self.df[self.value].isin(endpoint_value))]["bmdrc.Well.ID"]

    # Pull all endpoints
    all_endpoints = self.df[self.endpoint].unique().tolist()

    # Remove specific rows if applicable 
    if except_endpoint is not None:
        all_endpoints = [end for end in all_endpoints if end not in except_endpoint]

    # Set values to NA
    self.df[self.df["bmdrc.Well.ID"].isin(wells_rm)][self.value] = np.nan

    ################################
    ## ADD ATTRIBUTES FOR REPORTS ##
    ################################

    # Only add new inputs to the dictionary. 
    if hasattr(self, "report_well_removal"):
        self.report_well_removal.append([endpoint_name, endpoint_value, except_endpoint])
    else:
        self.report_combination = [endpoint_name, endpoint_value, except_endpoint]



def endpoint_combine(self, endpoint_dict):
    '''
    Combine endpoints and create new endpoints.
    For example, multiple 24 hour endpoints can be combined to create an "Any 24" endpoint.
    New endpoints are created with a binary or statement, meaning that if there is a 1 
    in any of the other endpoints, the resulting endpoint is a 1. Otherwise, it is 
    0 unless the other endpoints are all NA. Then the final value is NA.

    endpoint_dict: (dictionary) A dictionary where names are the new endpoint, and values are a list
    containing the endpoints to calculate these values from. 

    '''

    ############################
    ## CHECK INPUT PARAMETERS ##
    ############################
    
    # Assert that EndpointDictionary is a Dictionary
    if not isinstance(endpoint_dict, dict):
        raise Exception("EndpointDictionary is not a dict object.")
    
    # Assert that all values in the dictionary are endpoints in the DataClass object
    dict_unique_endpoints = list(set(sum(list(endpoint_dict.values()), [])))
    bc_unique_endpoints = list(set(self.df["endpoint"].tolist()))
    if not all([item in bc_unique_endpoints for item in dict_unique_endpoints]):
        raise Exception("All endpoints in EndpointDictionary must be in the DataClass endpoints.")
    
    #########################
    ## CREATE NEW ENDPOINT ##
    #########################

    # Define a small function to create new endpoints
    def new_endpoint(new_name, endpoints):
        sub_df = self.df[self.df[self.endpoint].isin(endpoints)].copy()
        sub_df[self.endpoint] = new_name
        sub_df = sub_df.groupby(by = [self.chemical, self.concentration, self.plate, self.well, self.endpoint], as_index = False).sum()
        sub_df[self.value].values[sub_df[self.value] > 1] = 1 
        return(sub_df)

    # Iterate through each dictionary entry 
    for NewEndpoint in endpoint_dict:

        if NewEndpoint in self.df[self.endpoint].unique():
            print(NewEndpoint + " is already an existing endpoint")
        else:
            self.df = pd.concat([self.df, new_endpoint(NewEndpoint, endpoint_dict[NewEndpoint])])

    ################################
    ## ADD ATTRIBUTES FOR REPORTS ##
    ################################

    # Only add new inputs to the dictionary. 
    if hasattr(self, "report_combination"):
        endpoint_dict.update(self.report_combination)
    
    # Set value 
    self.report_combination = endpoint_dict




    




    
        

