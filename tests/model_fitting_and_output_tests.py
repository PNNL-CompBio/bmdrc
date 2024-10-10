import pandas as pd
import numpy as np

import pytest
from bmdrc import BinaryClass, LPRClass

## How to calculate coverage (from within main package directory): 
# coverage run --source=bmdrc -m pytest -x tests/*
# coverage report
# coverage html

## Model Fitting and Output Tests ## 

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

# Run essential pre-processing steps
Long_Test.combine_and_create_new_endpoints({"ANY24":["NC24", "DP24", "SM24"], "ANY":["NC24", "DP24", "SM24", "JAW"]})
Long_Test.set_well_to_na(endpoint_name = "DNC", endpoint_value = 1, except_endpoint = ["ANY24"])
Long_Test.remove_endpoints("DNC")

# Run essential filtering steps 
Long_Test.filter_negative_control(apply = True)
Long_Test.filter_min_concentration(apply = True)
Long_Test.filter_correlation_score(apply = True)

# Ensure model fits works appropriately 
def test_model_fits():

    # GOF threshold must be a float
    with pytest.raises(Exception, match = "gof_threshold must be a float."):
        Long_Test.fit_models(gof_threshold = "watermelon")

    # GOF threshold must be between 0 and 1
    with pytest.raises(Exception, match = "gof_threshold must be larger than 0 or less than 1."):
        Long_Test.fit_models(gof_threshold = -0.1)

    # AIC threshold must be a float 
    with pytest.raises(Exception, match = "aic_threshold must be a float."):
        Long_Test.fit_models(aic_threshold = "melon")

    # Model selection can only be lowest_BMDL for now
    with pytest.raises(Exception, match = "Currently only 'lowest BMDL' is supported for model_selection."):
        Long_Test.fit_models(model_selection = "do nothing")

    # Fit models 
    Long_Test.fit_models(gof_threshold = 0.1, aic_threshold = 2, model_selection = "lowest BMDL")

    # Fit an unusual case when no pre-processing or filter steps occur
    Unusual = BinaryClass.BinaryClass(
        df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1), # Input is a pandas DataFrame
        chemical = "chemical.id", # The name of the chemical column 
        plate = "plate.id", # The name of the plate ID column
        well = "well", # The name of the column with well names
        concentration = "concentration", # The name of the concentration column
        endpoint = "endpoint", # The name of the column with endpoints
        value = "value", # The name of the column with values
        format = "long" # The format of the input data, either 'long' or 'wide' is accepted
    )
    Unusual.fit_models()

# Test curve generation functions
def test_gen_response_curve():

    # Chemical name must be real
    with pytest.raises(Exception, match = "rabbit is not a recognized chemical_name."):
        Long_Test.response_curve(chemical_name = "rabbit", endpoint_name = "JAW", model = "logistic")

    # Endpoint name must be real
    with pytest.raises(Exception, match = "freedom is not a recognized endpoint_name."):
        Long_Test.response_curve(chemical_name = 2, endpoint_name = "freedom", model = "logistic")

    # Model must be of the acceptable options
    with pytest.raises(Exception, match = "linear is not an acceptable model option. Acceptable options are: logistic, gamma, weibull, log logistic, probit, log probit, multistage2, quantal linear."):
        Long_Test.response_curve(chemical_name = 2, endpoint_name = "JAW", model = "linear")

    # Generate a response curve of every type
    Long_Test.response_curve(chemical_name = 2, endpoint_name = "JAW", model = "logistic")
    Long_Test.response_curve(chemical_name = 2, endpoint_name = "JAW", model = "gamma")
    Long_Test.response_curve(chemical_name = 2, endpoint_name = "JAW", model = "weibull")
    Long_Test.response_curve(chemical_name = 2, endpoint_name = "JAW", model = "log logistic")
    Long_Test.response_curve(chemical_name = 2, endpoint_name = "JAW", model = "probit")
    Long_Test.response_curve(chemical_name = 2, endpoint_name = "JAW", model = "log probit")
    Long_Test.response_curve(chemical_name = 2, endpoint_name = "JAW", model = "multistage2")
    Long_Test.response_curve(chemical_name = 2, endpoint_name = "JAW", model = "quantal linear")

# Test output module
def test_outputs():
    
    # Output benchmark doses
    Long_Test.output_benchmark_dose()

    # Output reports
    Long_Test.report(out_folder = "./long_report_delete")

    # Build a report for data that underwent no processing
    Quick_Test = LPRClass.LPRClass(
        df = pd.read_csv("data/LPR_Long.csv"),
        chemical = "chemical.id",
        plate = "plate.id",
        well = "well",
        concentration = "conc",
        time = "variable",
        value = "value",
        cycle_length = 20.0,
        cycle_cooldown = 10.0, 
        starting_cycle = "light"
    )
    Quick_Test.report(out_folder = "./lpr_report_delete")
