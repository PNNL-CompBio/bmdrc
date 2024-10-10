import pandas as pd
import numpy as np

import pytest
import bmdrc
from bmdrc import BinaryClass

## How to calculate coverage (from within main package directory): 
# coverage run --source=bmdrc -m pytest -x tests/*
# coverage report
# coverage html

## Filtering tests ## 

# Save example data 
Long_Test = BinaryClass.BinaryClass(
    df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1), # Input is a pandas DataFrame
    chemical = "chemical.id", # The name of the chemical column 
    plate = "plate.id", # The name of the plate ID column
    well = "well", # The name of the column with well names
    concentration = "concentration", # The name of the concentration column
    endpoint = "endpoint", # The name of the column with endpoints
    value = "value", # The name of the column with values
    format = "long" # The format of the input data, either 'long' or 'wide' is accepted
)

# Run esential pre-processing steps
Long_Test.combine_and_create_new_endpoints({"ANY24":["NC24", "DP24", "SM24"], "ANY":["NC24", "DP24", "SM24", "JAW"]})
Long_Test.set_well_to_na(endpoint_name = "DNC", endpoint_value = 1, except_endpoint = ["ANY24"])
Long_Test.remove_endpoints("DNC")

