import pandas as pd
import pytest
from bmdrc import BinaryClass

## How to run test functions: 
# coverage run -m bmdrc.BinaryClass tests/binary_class_tests.py
# coverage report
# coverage html

## Binary Class Tests ## 

# Test to ensure the function runs without error
def make_BinaryClass():
    BinaryClass(
        df = pd.read_csv("../data/Binary_Simplified_Long.csv").drop("Notes", axis = 1), # Input is a pandas DataFrame
        chemical = "chemical.id", # The name of the chemical column 
        plate = "plate.id", # The name of the plate ID column
        well = "well", # The name of the column with well names
        concentration = "concentration", # The name of the concentration column
        endpoint = "endpoint", # The name of the column with endpoints
        value = "value", # The name of the column with values
        format = "long" # The format of the input data, either 'long' or 'wide' is accepted
    )
    
# Test errors 
def test_df():

    BinaryClass.BinaryClass(
        df = pd.DataFrame(),
        chemical = "chemical.id", # The name of the chemical column 
        plate = "plate.id", # The name of the plate ID column
        well = "well", # The name of the column with well names
        concentration = "concentration", # The name of the concentration column
        endpoint = "endpoint", # The name of the column with endpoints
        value = "value", # The name of the column with values
        format = "long" # The format of the input data, either 'long' or 'wide' is accepted
    )

    


