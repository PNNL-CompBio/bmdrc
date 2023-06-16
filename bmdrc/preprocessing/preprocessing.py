
import pandas as pd
from bmdrc.input_data_classes.BinaryClass import BinaryClass as BC

__author__ = "David Degnan"

def combine_and_create_new_endpoints(BinaryClass, EndpointDictionary):
    '''
    Combine endpoints and create new endpoints for only the BinaryClass datatype.
    For example, multiple 24 hour endpoints can be combined to create an "Any 24" endpoint.
    New endpoints are created with a binary or statement, meaning that if there is a 1 
    in any of the other endpoints, the resulting endpoint is a 1. Otherwise, it is 
    0 unless the other endpoints are all NA. Then the final value is NA.
    
    BinaryClass: (bmdrc BinaryClass object) Build with the BinaryClass function

    EndpointDictionary: (dictionary) A dictionary where names are the new endpoint, and values are a list
    containing the endpoints to calculate these values from. 

    '''

    ############################
    ## CHECK INPUT PARAMETERS ##
    ############################

    # Assert that BinaryClass is the correct object type
    if not isinstance(BinaryClass, BC):
        raise Exception("BinaryClass is not a BinaryClass object.")
    
    # Assert that EndpointDictionary is a Dictionary
    if not isinstance(EndpointDictionary, dict):
        raise Exception("EndpointDictionary is not a dict object.")
    
    # Assert that all values in the dictionary are endpoints in the BinaryClass object
    dict_unique_endpoints = list(set(sum(list(EndpointDictionary.values()), [])))
    bc_unique_endpoints = list(set(BinaryClass.df["endpoint"].tolist()))
    if not all([item in bc_unique_endpoints for item in dict_unique_endpoints]):
        raise Exception("All endpoints in EndpointDictionary must be in the BinaryClass endpoints.")
    
    #########################
    ## CREATE NEW ENDPOINT ##
    #########################

    # Define a small function to create new endpoints
    def new_endpoint(endpoints, new_name):
        sub_df = df_morph[df_morph["endpoint"].isin(endpoints)]
        sub_df["endpoint"] = new_name
        sub_df = sub_df.groupby(by = ["chemical.id", "conc", "plate.id", "well", "endpoint"], as_index = False).sum()
        sub_df['value'].values[sub_df['value'] > 1] = 1 
        return(sub_df)

    # Iterate through each dictionary entry 
    for NewEndpoint, ExistingEndpoints in EndpointDictionary:
        pd.concat([BinaryClass.df, new_endpoint(NewEndpoint, ExistingEndpoints)])


    
        

