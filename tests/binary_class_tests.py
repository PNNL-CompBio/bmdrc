import pandas as pd
import pytest
import bmdrc
from bmdrc import BinaryClass

## How to calculate coverage (from within main package directory): 
# coverage run --source=bmdrc -m pytest -x tests/*
# coverage report
# coverage html

## Binary Class Tests ## 

# Test to ensure the function runs without error
def make_BinaryClass():

    LongTest = BinaryClass.BinaryClass(
        df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1), # Input is a pandas DataFrame
        chemical = "chemical.id", # The name of the chemical column 
        plate = "plate.id", # The name of the plate ID column
        well = "well", # The name of the column with well names
        concentration = "concentration", # The name of the concentration column
        endpoint = "endpoint", # The name of the column with endpoints
        value = "value", # The name of the column with values
        format = "long" # The format of the input data, either 'long' or 'wide' is accepted
    )

    WideTest = BinaryClass.BinaryClass(
        df = pd.read_csv("data/Binary_Morphology_Wide.csv"),
        chemical = "chemical.id",
        plate = "plate.id",
        well = "well",
        concentration = "conc",
        endpoint = "endpoint",
        value = "value",
        format = "wide"
    )
    
# Test wrong inputs 
def test_wrong_inputs():

    with pytest.raises(Exception, match = "df cannot be empty. Please provide a pandas DataFrame.") as test_info:
        BinaryClass.BinaryClass(
            df = pd.DataFrame(),
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = "concentration", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )
    assert test_info.type == Exception

def run_tests():
    make_BinaryClass()
    test_wrong_inputs()

run_tests()

