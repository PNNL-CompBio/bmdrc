import pandas as pd
import numpy as np

import pytest
from bmdrc import BinaryClass

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=SyntaxWarning)

## How to calculate coverage (from within main package directory): 
# coverage run --source=bmdrc -m pytest -x tests/*
# coverage report
# coverage html

## Preprocessing tests ## 

# Create data Tests
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

# Ensure well to NA works as expected 
def test_well_to_na():

    # Set do not count cases to NA
    Long_Test.set_well_to_na(endpoint_name = "DNC", endpoint_value = 1, except_endpoint = ["JAW"])

    # Set NC24 cases to NA as well
    Long_Test.set_well_to_na(endpoint_name = "NC24", endpoint_value = 1, except_endpoint = ["JAW"])

    # Well to NA dictionary should have both these values
    assert Long_Test.report_well_na == [[['DNC'], [1], ['JAW']], [['NC24'], [1], ['JAW']]]

    # All things with DNC and NC24 with ones should now be NA 
    IDs_Removed = set(Long_Test.df[(Long_Test.df[Long_Test.endpoint].isin(["DNC", "NC24"])) & (np.isnan(Long_Test.df[Long_Test.value]))]["bmdrc.Well.ID"])
    assert IDs_Removed == {'1 0.0 B 6', '1 0.0 C 6', '1 1.5 A 2', '1 1.5 B 6', '1 1.5 C 6', '1 3.0 B 1', '1 6.0 A 1'}

# Ensure combine endpoint_name works as expected 
def test_endpoint_combine(): 

    # Add new endpoint
    endpoint_dict = {"ANY24":["NC24", "DP24", "SM24"], 
                     "ANY":["NC24", "DP24", "SM24", "JAW"]}
    Long_Test.combine_and_create_new_endpoints(endpoint_dict)

    # Add another endpoint
    Long_Test.combine_and_create_new_endpoints({"THE24":["NC24", "SM24"], "DEL":["JAW"]})

    # Test the stored dictionary dictionary
    assert Long_Test.report_combination == {"ANY24":["NC24", "DP24", "SM24"], "ANY":["NC24", "DP24", "SM24", "JAW"], "THE24":["NC24", "SM24"], "DEL":["JAW"]}

# Ensure remove endpoints works as expected
def test_remove_endpoints():

    # Remove endpoint
    Long_Test.remove_endpoints(["DNC", "THE24"])

    # Remove more endpoints
    Long_Test.remove_endpoints(["DEL"])

    # Test the stored removed endpoints 
    assert Long_Test.report_endpoint_removal == ["DNC", "THE24", "DEL"]
#### removed ####
# Ensure value error if endpoint name is not in the dataset
def test_well_to_na_invalid_except_endpoint_raises():
    try:
        Long_Test.set_well_to_na(
            endpoint_name="DNC",
            endpoint_value=1,
            except_endpoint=["NOT A REAL ENDPOINT"]
        )
        assert False, "EXPECTED VALUEERROR NOT RAISED"
    except ValueError:
        pass

# Ensure scalar endpoint is converted to a list and works
def test_combine_endpoints_scalar_endpoints():
    Long_Test.combine_and_create_new_endpoints({"ANY SCALAR": "NC24"})

    assert "ANY SCALAR" in Long_Test.report_combination

# Ensure if any component endpoint does not exist in the dataframe, raise Exception
def test_combine_endpoints_invalid_endpoint_raises():
    try:
        Long_Test.combine_and_create_new_endpoints({
            "ANY_INVALID": ["NOT A REAL ENDPOINT"]
        })
        assert False, "EXCEPTION NOT RAISED FOR INVALID ENDPOINT"
    except Exception:
        pass

#Ensure if NewEnpoint exists, raise Exception
def test_combine_existing_endpoint_raises():
    try:
        Long_Test.combine_and_create_new_endpoints({
            "NC24": ["DP24"]  # endpoint name already exists
        })
        assert False, "EXPECTED EXCEPTION NOT RAISED FOR EXISTING ENDPOINT"
    except Exception:
        pass

# Ensure the removal of a non-existent endpoint raises ValueError
def test_remove_endpoints_invalid_endpoint_raises():
    try:
        Long_Test.remove_endpoints(["NOT A REAL ENDPOINT"])
        assert False, "EXPECTED VALUEERROR NOT RAISED"
    except ValueError:
        pass

# Ensure new BinaryClass for each test (prevents cross-test mutation)
def make_long_test():
    return BinaryClass.BinaryClass(
        df=pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis=1),
        chemical="chemical.id",
        plate="plate.id",
        well="well",
        concentration="concentration",
        endpoint="endpoint",
        value="value",
        format="long"
    )
# Ensure invalid endpoint raises an error when trying to set to NA
def test_well_to_na_endpoint_name_as_list():
    Long_Test = make_long_test()

    # pick a real endpoint from the data instead of hardcoding "DNC"
    endpoint = Long_Test.df[Long_Test.endpoint].unique().tolist()[0]

    Long_Test.set_well_to_na(endpoint_name=[endpoint], endpoint_value=1, except_endpoint=["JAW"])

    assert Long_Test.report_well_na[-1][0][0] == endpoint  # stored endpoint name

# Ensure a scalar invalid except endpoint reaches the validation loop and raises ValueError.
def test_well_to_na_invalid_except_endpoint_scalar_raises():
    Long_Test = make_long_test()
    endpoint = Long_Test.df[Long_Test.endpoint].unique().tolist()[0]

    with pytest.raises(ValueError, match = "is not an endpoint in the DataClass object"):
        Long_Test.set_well_to_na(endpoint_name = endpoint, endpoint_value = 1, except_endpoint = "NOT A REAL ENDPOINT")