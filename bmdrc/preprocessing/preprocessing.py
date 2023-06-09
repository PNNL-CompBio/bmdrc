
__author__ = "David Degnan"

def combine_endpoints(BinaryClass, EndpointDictionary):
    '''
    Combine endpoints and create new endpoints for only the BinaryClass datatype.
    For example, multiple 24 hour endpoints can be combined to create an "Any 24" endpoint.
    New endpoints are created with a binary or statement, meaning that if there is a 1 
    in any of the other endpoints, the resulting endpoint is a 1. Otherwise, it is 
    0 unless the other endpoints are all NA. Then the final value is NA.
    
    BinaryClass: (bmdrc BinaryClass object) 

    EndpointDictionary: (dictionary) A dictionary where names are the new endpoint, and values are a list
    containing the endpoints to calculate these values from.

    '''