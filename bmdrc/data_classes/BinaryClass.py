import operator
import pandas as pd
import numpy as np
import sys 

__author__ = "David Degnan"

class Binary(object):
    '''
    Generates a binary class object where input values are either a 0 or a 1, 
    which requires a designation for the following columns:

    1. chemicalID
    2. concentration
    3. plateID
    4. well
    5. endpoint
    6. value --> here, should be a 0 or a 1

    If the data is in long format, the "endpoint" and "value" designations will be
    created before the class call. 
    '''

    # Define the input checking functions 
    def __init__(self, chemicalID, concentration, plateID, well, endpoint, value):
        self.chemicalID = chemicalID
        self.concentration = concentration
        self.plateID = plateID
        self.well = well
        self.endpoint = endpoint
        self.value = value

    # Set property returning functions 
    chemicalID = property(operator.attrgetter('_chemicalID'))
    concentration = property(operator.attrgetter('_concentration'))
    plateID = property(operator.attrgetter('_plateID'))
    well = property(operator.attrgetter('_well'))
    endpoint = property(operator.attrgetter('_endpoint'))
    value = property(operator.attrgetter('_value'))

    # Now, ensure each input is correct 

    @chemicalID.setter
    def chemicalID(self, chemicalIDs):
        if not chemicalIDs: 
            raise Exception("chemicalID cannot be empty. Please provide a list \
                            of strings representing each chemical's ID")
        if not isinstance(chemicalIDs, str):
            raise Exception("chemicalID must be a list of strings.")