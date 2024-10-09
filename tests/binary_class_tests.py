import pandas as pd
import pytest
import bmdrc
from bmdrc import BinaryClass

## How to calculate coverage (from within main package directory): 
# coverage run --source=bmdrc -m pytest -x tests/*
# coverage report
# coverage html

## Binary Class Tests ## 

# Test to ensure long data runs without error 
def test_long_BinaryClass():

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
    assert isinstance(LongTest, bmdrc.BinaryClass.BinaryClass)

# Test to ensure wide data runs without error
def test_wide_BinaryClass():

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
    assert isinstance(WideTest, bmdrc.BinaryClass.BinaryClass)
    
# Test wrong inputs for data.frame 
def test_df():

    # The df must be a pandas DataFrame, no exceptions
    with pytest.raises(Exception, match = "df must be a pandas DataFrame"):
        BinaryClass.BinaryClass(
            df = "celery",
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = "concentration", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )
    
    # The df must be a pandas DataFrame with data in it 
    with pytest.raises(Exception, match = "df cannot be empty. Please provide a pandas DataFrame."):
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

# Test wrong inputs for chemicals
def test_chemical():

    # The chemical must be a string
    with pytest.raises(Exception, match = "chemical must be a name of a column in df."):
        BinaryClass.BinaryClass(
            df = pd.read_csv("data/Binary_Morphology_Wide.csv"),
            chemical = 3, 
            plate = "plate.id", 
            well = "well", 
            concentration = "concentration", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )

    # The chemical must be a name in the dataframe 
    with pytest.raises(Exception, match = "cantelope is not in the column names of df."):
        BinaryClass.BinaryClass(
            df = pd.read_csv("data/Binary_Morphology_Wide.csv"),
            chemical = "cantelope", 
            plate = "plate.id", 
            well = "well", 
            concentration = "concentration", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )
    
    # The chemical name must not be an unaceeptable name 
    with pytest.raises(Exception, match = "bmdrc.Well.ID is not a permitted name. Please rename this column."):

        df2 = pd.read_csv("data/Binary_Morphology_Wide.csv")
        df2 = df2.rename({"chemical.id":"bmdrc.Well.ID"}, axis = 1)

        BinaryClass.BinaryClass(
            df = df2,
            chemical = "bmdrc.Well.ID", 
            plate = "plate.id", 
            well = "well", 
            concentration = "concentration", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )