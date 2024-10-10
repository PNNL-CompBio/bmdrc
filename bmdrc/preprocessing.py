import pandas as pd
import numpy as np

__author__ = "David Degnan"

def well_to_na(self, endpoint_name, endpoint_value, except_endpoint):
    '''
    Remove any wells where a specific endpoint has a specific value. 
    Wells are set to NA.

    endpoint_name: (list - string) list of endpoints to remove 

    endpoint_value: (list - numeric) specific value or list of values that the endpoint needs to have to remove the well

    except_endpoint: (list - string) list of endpoints that should not have their wells affected
    '''

    ############################
    ## CHECK INPUT PARAMETERS ##
    ############################
    
    # Convert endpoint name to a list
    if isinstance(endpoint_name, list) == False:
        endpoint_name = [endpoint_name]

    # Iterate through each endpoint to confirm it is a valid choice 
    for endpoint in endpoint_name:
        if (endpoint in self.df[self.endpoint].unique().tolist()) == False:
            raise Exception(endpoint + " is not an endpoint in the DataClass object.")

    # Convert endpoint value to a list
    if isinstance(endpoint_value, list) == False:
        endpoint_value = [endpoint_value]
    
    # Iterate through each except endpoint to confirm they are valid choices
    if except_endpoint is not None:

        # Convert except endpoints to a list
        if isinstance(except_endpoint, list) == False:
            except_endpoint = [except_endpoint]

        # Check each endpoint
        for endpoint in except_endpoint:
            if (endpoint in self.df[self.endpoint].unique().tolist()) == False:
                raise Exception(endpoint + " is not an endpoint in the DataClass object.")
            
    ####################################
    ## SET WELLS TO NA TO REMOVE THEM ##
    ####################################

    # If it hasn't been made yet, make the bmdrc.Well.ID column
    if "bmdrc.Well.ID" not in self.df.columns.tolist():
        self.df["bmdrc.Well.ID"] = self.df[self.chemical].astype(str) + " " + self.df[self.concentration].astype(str) + " " + self.df[self.plate].astype(str) + " " + self.df[self.well].astype(str)
    
    # List wells to remove
    wells_rm = self.df[(self.df[self.endpoint].isin(endpoint_name)) & (self.df[self.value].isin(endpoint_value))]["bmdrc.Well.ID"]

    # Pull all endpoints as acceptalbe
    acc_endpoints = self.df[self.endpoint].unique().tolist()

    # Remove specific rows if applicable 
    if except_endpoint is not None:
        acc_endpoints = [end for end in acc_endpoints if end not in except_endpoint]

    # Set values to NA
    self.df.loc[self.df["bmdrc.Well.ID"].isin(wells_rm) & self.df[self.endpoint].isin(acc_endpoints), self.value] = np.nan

    ################################
    ## ADD ATTRIBUTES FOR REPORTS ##
    ################################

    # New attributes 
    new_attributes = [endpoint_name, endpoint_value, except_endpoint]

    # Only add new inputs to the dictionary. 
    if hasattr(self, "report_well_na"):
        self.report_well_na.append(new_attributes)
    else:
        self.report_well_na = [new_attributes]

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
    
    #########################
    ## CREATE NEW ENDPOINT ##
    #########################

    # Define a small function to create new endpoints
    def new_endpoint(new_name, endpoints):
        
        # Convert endpoints to list 
        if isinstance(endpoints, list) == False:
            endpoints = [endpoints]

        # If endpoints are not in the data.frame, trigger error 
        for endpoint in endpoints: 
            if (endpoint in self.df[self.endpoint].unique().tolist()) == False:
                raise Exception(endpoint + " is not an endpoint in the DataClass object.")

        # Combine endpoints 
        sub_df = self.df[self.df[self.endpoint].isin(endpoints)].copy()
        sub_df[self.endpoint] = new_name
        sub_df = sub_df.groupby(by = [self.chemical, self.concentration, self.plate, self.well, self.endpoint], as_index = False).sum()
        sub_df[self.value].values[sub_df[self.value] > 1] = 1 
        return(sub_df)

    # Iterate through each dictionary entry 
    for NewEndpoint in endpoint_dict:

        if NewEndpoint in self.df[self.endpoint].unique().tolist():
            raise Exception(NewEndpoint + " is already an existing endpoint")
        else:
            self.df = pd.concat([self.df, new_endpoint(NewEndpoint, endpoint_dict[NewEndpoint])])

    ################################
    ## ADD ATTRIBUTES FOR REPORTS ##
    ################################

    # Only add new inputs to the dictionary. 
    if hasattr(self, "report_combination"):
        self.report_combination = self.report_combination | endpoint_dict
    else:
        self.report_combination = endpoint_dict

def remove_endpoints(self, endpoint_name):
    '''
    Completely remove an endpoint or set of endpoints from the dataset

    endpoint_name: (list - string) A list of endpoints to remove
    '''

    ############################
    ## CHECK INPUT PARAMETERS ##
    ############################

    # Convert endpoint name to a list
    if isinstance(endpoint_name, list) == False:
        endpoint_name = [endpoint_name]

    # Iterate through each endpoint to confirm it is a valid choice 
    for endpoint in endpoint_name:
        if (endpoint in self.df[self.endpoint].unique().tolist()) == False:
            raise Exception(endpoint + " is not an endpoint in the DataClass object.")
        
    ######################
    ## REMOVE ENDPOINTS ##
    ######################

    self.df = self.df[self.df[self.endpoint].isin(endpoint_name) == False]

    ################################
    ## ADD ATTRIBUTES FOR REPORTS ##
    ################################

    # Only add new inputs to the dictionary. 
    if hasattr(self, "report_endpoint_removal"):
        self.report_endpoint_removal.extend(endpoint_name)
    else:
        self.report_endpoint_removal = endpoint_name
