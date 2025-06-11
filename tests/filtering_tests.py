import pandas as pd
import numpy as np

import pytest
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

# Run essential pre-processing steps
Long_Test.combine_and_create_new_endpoints({"ANY24":["NC24", "DP24", "SM24"], "ANY":["NC24", "DP24", "SM24", "JAW"]})
Long_Test.set_well_to_na(endpoint_name = "DNC", endpoint_value = 1, except_endpoint = ["ANY24"])
Long_Test.remove_endpoints("DNC")

# Test the negative control filter
def test_negative_control_filter():

    # Small percentages are permitted
    Long_Test.filter_negative_control(percentage = 0.01)

    # Recommended default is 50% or higher
    Long_Test.filter_negative_control(apply = True)

    # Ensure the stored parameter is correct 
    assert Long_Test.filter_negative_control_thresh == 50

# Test the minimum concentration filter
def test_minimum_concentration_filter():

    New_Test = BinaryClass.BinaryClass(
        df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1), # Input is a pandas DataFrame
        chemical = "chemical.id", # The name of the chemical column 
        plate = "plate.id", # The name of the plate ID column
        well = "well", # The name of the column with well names
        concentration = "concentration", # The name of the concentration column
        endpoint = "endpoint", # The name of the column with endpoints
        value = "value", # The name of the column with values
        format = "long" # The format of the input data, either 'long' or 'wide' is accepted
    )

    # Run minimum concentration count
    Long_Test.filter_min_concentration(apply = True)

    # Run minimum concentration on a clean dataset
    New_Test.filter_min_concentration(apply = False, diagnostic_plot = False)

    # Tracked value should be the default of 3
    assert Long_Test.filter_min_concentration_thresh == 3

# Test the correlation score filter 
def test_correlation_score_filter():

    New_Test = BinaryClass.BinaryClass(
        df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1), # Input is a pandas DataFrame
        chemical = "chemical.id", # The name of the chemical column 
        plate = "plate.id", # The name of the plate ID column
        well = "well", # The name of the column with well names
        concentration = "concentration", # The name of the concentration column
        endpoint = "endpoint", # The name of the column with endpoints
        value = "value", # The name of the column with values
        format = "long" # The format of the input data, either 'long' or 'wide' is accepted
    )

    # Run correlation score test 
    Long_Test.filter_correlation_score(apply = True)

    # Run correlation score test on a clean dataset
    New_Test.filter_correlation_score(apply = False, diagnostic_plot = False)

    # Tracked value should be the default of 0.2
    assert Long_Test.filter_correlation_score_thresh == 0.2

    
