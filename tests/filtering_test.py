import pandas as pd
import numpy as np
import pytest
from bmdrc import BinaryClass
from bmdrc import ProportionalClass
from bmdrc import ContinuousClass
import matplotlib
matplotlib.use("Agg")

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=SyntaxWarning)

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
Long_Test.combine_and_create_new_endpoints({"ANY24":["DP24", "SM24", "JAW"]})

# Test the negative control filter
def test_negative_control_filter():

    # Small percentages are permitted
    Long_Test.filter_negative_control(percentage = 0.01)

    # Percentages greater than a 100 and less than 0 are not permitted, reset to default of 50%
    Long_Test.filter_negative_control(percentage = 150)

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

    # Count below 1 is not permitted, reset to default of 1
    New_Test.filter_min_concentration(count = 0, apply = False, diagnostic_plot = False)
    assert New_Test.filter_min_concentration_thresh == 1

    # A high count than the data creates a removal
    New_Test.filter_min_concentration(count = 10, apply = True)
    assert (New_Test.plate_groups["bmdrc.filter"] == "Remove").any()

    # Diagnostic plot should be created when requested
    New_Test.filter_min_concentration(apply = False, diagnostic_plot = True)

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

    # DP24 has a constant response across all concentration giving a NaN Spearman correlation
    assert not Long_Test.filter_correlation_score_df["Spearman"].isna().any()

    #  Test correlation score when = -5 and that it clipped to nearest bound. 
    New_Test.filter_correlation_score(score = -5, apply = False, diagnostic_plot = False)
    assert New_Test.filter_correlation_score_thresh == -1

    #  Test correlation score when = 5 and that it clipped to nearest bound. 
    New_Test.filter_correlation_score(score = 5, apply = False, diagnostic_plot = False)
    assert New_Test.filter_correlation_score_thresh == 1

    # Test correlation when direction = "above" the threshold
    New_Test.filter_correlation_score(score = 0.2, direction = "above", apply = True)
    
    # Between was being applied right after the above test, so we need to reset the data to a clean state
    Between_Test = BinaryClass.BinaryClass(  
        df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1), # Input is a pandas DataFrame
        chemical = "chemical.id", # The name of the chemical column 
        plate = "plate.id", # The name of the plate ID column
        well = "well", # The name of the column with well names
        concentration = "concentration", # The name of the concentration column
        endpoint = "endpoint", # The name of the column with endpoints
        value = "value", # The name of the column with values
        format = "long" # The format of the input data, either 'long' or 'wide' is accepted
    )
    
    # Test correlation when direction = "between" the threshold
    Between_Test.filter_correlation_score(score = 0.2, direction = "between", apply = True)

    # Test unrecognized direction raises a ValueError
    with pytest.raises(ValueError):
        New_Test.filter_correlation_score(score = 0.2, direction = "invalid_direction", apply = True)
    
    # Diagnostic plot should be created when requested
    New_Test.filter_correlation_score(apply = False, diagnostic_plot = True)

# Save example proportional data
Proportional_Test = ProportionalClass.ProportionalClass(
    df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1), # Input is a pandas DataFrame
    chemical = "chemical.id", # The name of the chemical column
    concentration = "concentration", # The name of the concentration column
    endpoint = "endpoint", # The name of the column with endpoints
    response = "value", # The name of the response column
)

# Save example continuous data
Continuous_Test = ContinuousClass.ContinuousClass(
    df = pd.read_csv("data/Continuous.txt", sep = "\t").drop("Notes", axis = 1), # Input is a pandas DataFrame
    chemical = "Chemical ID", # The name of the chemical column
    concentration = "Concentration_uM", # The name of the concentration column
    endpoint = "Endpoint", # The name of the column with endpoints
    response = "Measurement", # The name of the response column
)

# Test "no plate" branches of make_plate_groups for proportional and continuous data
def test_proportional_correlation_score_filter():

    # Proportional data
    Proportional_Test.filter_correlation_score(apply = True)
    assert hasattr(Proportional_Test, "plate_groups")

# Test continuous negativeontrol filter
def test_continuous_negative_control_filter():
    Continuous_Test.filter_negative_control(apply = False, diagnostic_plot = False)
    Continuous_Test.filter_negative_control(apply = True)
    Continuous_Test.filter_min_concentration(apply = False, diagnostic_plot = True)
    assert "Filter" in Continuous_Test.filter_min_concentration_df.columns