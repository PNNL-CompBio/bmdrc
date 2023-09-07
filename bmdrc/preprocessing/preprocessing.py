import pandas as pd
import ipdb

__author__ = "David Degnan"

def endpoint_combine(self, endpoint_dict):
    '''
    Combine endpoints and create new endpoints for only the BinaryClass datatype.
    For example, multiple 24 hour endpoints can be combined to create an "Any 24" endpoint.
    New endpoints are created with a binary or statement, meaning that if there is a 1 
    in any of the other endpoints, the resulting endpoint is a 1. Otherwise, it is 
    0 unless the other endpoints are all NA. Then the final value is NA.
    
    dataclass: (bmdrc dataclass object) Build with the DataClass function. Acceptable options are
    BinaryClass, LPRClass, and MovementClass

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


    
        

