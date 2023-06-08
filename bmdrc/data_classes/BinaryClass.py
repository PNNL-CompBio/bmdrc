import operator

__author__ = "David Degnan"

class BinaryClass(object):
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
        self._chemicalID = chemicalIDs
        
    @concentration.setter
    def concentration(self, concentrations):
        if not concentrations: 
            raise Exception("concentration cannot be empty. Please provide a list \
                            of numerics representing each chemical's concentration")
        if not isinstance(concentrations, float):
            raise Exception("concentration must be a list of numerics.")
        self._concentration = concentrations
        
    @plateID.setter
    def plateID(self, plateIDs):
        if not plateIDs: 
            raise Exception("plateID cannot be empty. Please provide a list \
                            of strings representing each plateID")
        if not isinstance(plateIDs, str):
            raise Exception("plateID must be a list of strings.")
        self._plateID = plateIDs
        
    @well.setter
    def well(self, wells):
        if not wells: 
            raise Exception("well cannot be empty. Please provide a list \
                            of strings representing each well")
        if not isinstance(wells, str):
            raise Exception("well must be a list of strings.")
        self._well = wells
        
    @endpoint.setter
    def endpoint(self, endpoints):
        if not endpoints: 
            raise Exception("endpoint cannot be empty. Please provide a list \
                            of strings representing each endpoint.")
        if not isinstance(endpoints, str):
            raise Exception("endpoint must be a list of strings.")
        self._endpoint = endpoints
        
        
    @value.setter
    def value(self, values):
        if not values: 
            raise Exception("value cannot be empty. Please provide a list \
                            of zeroes and ones representing whether a condition \
                            is absent or present, respectively.")
        if all(val == 0 or val == 1 for val in values):
            raise Exception("value must be a list of zeroes or ones indicating \
                            whether a condition is absent or present, respectively.")
        self._value = values