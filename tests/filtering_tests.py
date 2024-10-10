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

    # Percentage must be a number 
    with pytest.raises(Exception, match = "percentage must be a float."):
        Long_Test.filter_negative_control(percentage = "WaffleCone")

    # Percentage must be between 0 and 100
    with pytest.raises(Exception, match = "percentage must be between 0 and 100."):
        Long_Test.filter_negative_control(percentage = -1)

    # Small percentages are permitted
    Long_Test.filter_negative_control(percentage = 0.01)

    # Apply must be a true or false
    with pytest.raises(Exception, match = "apply must be a True or False."):
        Long_Test.filter_negative_control(apply = "Ahi Tuna")

    # Diagnostic plot must be true or false
    with pytest.raises(Exception, match = "diagnostic_plot must be a True or False."):
        Long_Test.filter_negative_control(diagnostic_plot = "Sashimi")

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

    # Count needs to be an integer
    with pytest.raises(Exception, match = "count must be an integer."):
        Long_Test.filter_min_concentration(count = "French Fries")

    # Count can't be a negative number
    with pytest.raises(Exception, match = "count must be at least 1."):
        Long_Test.filter_min_concentration(count = 0)

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

    # Score needs to be an integer
    with pytest.raises(Exception, match = "score must be an integer."):
        Long_Test.filter_correlation_score(score = "Carrots")

    # Score should be between -1 and 1
    with pytest.raises(Exception, match = "score must be larger than -1 or less than 1 to filter any values."):
        Long_Test.filter_correlation_score(score = 1.2)

    # Run correlation score test 
    Long_Test.filter_correlation_score(apply = True)

    # Run correlation score test on a clean dataset
    New_Test.filter_correlation_score(apply = False, diagnostic_plot = False)

    # Tracked value should be the default of 0.2
    assert Long_Test.filter_correlation_score_thresh == 0.2

    
