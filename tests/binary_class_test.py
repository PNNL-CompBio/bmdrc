import numpy as np
import pandas as pd
import pytest
import bmdrc
from bmdrc import BinaryClass

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=SyntaxWarning)


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
    assert (LongTest.df.columns == ['chemical.id', 'concentration', 'plate.id', 'well', 'endpoint', 'value']).all()

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
    assert (WideTest.df.columns == ['chemical.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']).all()

    # Run a quick check for plate groups
    WideTest.make_plate_groups()
    
# Test wrong inputs for data.frame 
def test_df():

    # The df must be a pandas DataFrame, no exceptions
    with pytest.raises(Exception, match = "df must be a pandas DataFrame."):
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
            df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1),
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
            df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1),
            chemical = "cantelope", 
            plate = "plate.id", 
            well = "well", 
            concentration = "concentration", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )
    
    # The chemical name must not be an unacceptable name 
    with pytest.raises(Exception, match = "bmdrc.Well.ID is not a permitted name. Please rename this column."):

        df2 = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1)
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

# Test wrong inputs for plates
def test_plates():

    # The plate must be a string
    with pytest.raises(Exception, match = "plate must be a name of a column in df."):
        BinaryClass.BinaryClass(
            df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1),
            chemical = "chemical.id", 
            plate = 42, 
            well = "well", 
            concentration = "concentration", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )

    # The plate must be a name in the dataframe 
    with pytest.raises(Exception, match = "taco salad is not in the column names of df."):
        BinaryClass.BinaryClass(
            df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1),
            chemical = "chemical.id", 
            plate = "taco salad", 
            well = "well", 
            concentration = "concentration", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )
    
    # The chemical name must not be an unacceptable name 
    with pytest.raises(Exception, match = "bmdrc.num.tot is not a permitted name. Please rename this column."):

        df2 = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1)
        df2 = df2.rename({"plate.id":"bmdrc.num.tot"}, axis = 1)

        BinaryClass.BinaryClass(
            df = df2,
            chemical = "chemical.id", 
            plate = "bmdrc.num.tot", 
            well = "well", 
            concentration = "concentration", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )

# Test wrong inputs for wells
def test_wells():

    # The well must be a string
    with pytest.raises(Exception, match = "well must be a name of a column in df."):
        BinaryClass.BinaryClass(
            df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1),
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = False, 
            concentration = "concentration", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )

    # The well must be a name in the dataframe 
    with pytest.raises(Exception, match = "Tuesday is not in the column names of df."):
        BinaryClass.BinaryClass(
            df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1),
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "Tuesday", 
            concentration = "concentration", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )
    
    # The well must not be an unacceptable name 
    with pytest.raises(Exception, match = "bmdrc.num.nonna is not a permitted name. Please rename this column."):

        df2 = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1)
        df2 = df2.rename({"well":"bmdrc.num.nonna"}, axis = 1)

        BinaryClass.BinaryClass(
            df = df2,
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "bmdrc.num.nonna", 
            concentration = "concentration", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )

# Test wrong inputs for concentrations
def test_concs():

    # The concentration must be a string
    with pytest.raises(Exception, match = "concentration must be a name of a column in df."):
        BinaryClass.BinaryClass(
            df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1),
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = 22, 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )

    # The concentration must be a name in the dataframe 
    with pytest.raises(Exception, match = "Star Wars is not in the column names of df."):
        BinaryClass.BinaryClass(
            df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1),
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = "Star Wars", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )
    
    # The concentration name must not be an unacceptable name 
    with pytest.raises(Exception, match = "bmdrc.frac.affected is not a permitted name. Please rename this column."):

        df2 = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1)
        df2 = df2.rename({"concentration":"bmdrc.frac.affected"}, axis = 1)

        BinaryClass.BinaryClass(
            df = df2,
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = "bmdrc.frac.affected", 
            endpoint = "endpoint",
            value = "value", 
            format = "long"
        )

# Test wrong inputs for format
def test_format():

    # Format must be wide or long
    with pytest.raises(Exception, match = "format must be 'long' or 'wide'."):
        BinaryClass.BinaryClass(
            df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1),
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = "concentration", 
            endpoint = "endpoint",
            value = "value", 
            format = "short"
        )

# Test wrong inputs for endpoints
def test_endpoints():

    # The endpoint must be a string
    with pytest.raises(Exception, match = "endpoint must be a name of a column in df."):
        BinaryClass.BinaryClass(
            df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1),
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = "concentration", 
            endpoint = 55,
            value = "value", 
            format = "long"
        )

    # The endpoint must be a name in the dataframe 
    with pytest.raises(Exception, match = "Python is not in the column names of df."):
        BinaryClass.BinaryClass(
            df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1),
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = "concentration", 
            endpoint = "Python",
            value = "value", 
            format = "long"
        )
    
    # The endpoint must not be an unacceptable name 
    with pytest.raises(Exception, match = "bmdrc.frac.affected is not a permitted name. Please rename this column."):

        df2 = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1)
        df2 = df2.rename({"endpoint":"bmdrc.frac.affected"}, axis = 1)

        BinaryClass.BinaryClass(
            df = df2,
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = "concentration", 
            endpoint = "bmdrc.frac.affected",
            value = "value", 
            format = "long"
        )

# Test wrong inputs for values
def test_values():

    # Ensure non-string value column names trigger the documented validation error
    with pytest.raises(Exception, match = "value must be a name of a column in df."):
        BinaryClass.BinaryClass(
            df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1),
            chemical = "chemical.id",
            plate = "plate.id",
            well = "well",
            concentration = "concentration",
            endpoint = "endpoint",
            value = 22,
            format = "long"
        )

    # Ensure missing value column names trigger the documented lookup error
    with pytest.raises(Exception, match = "Gogurt is not in the column names of df."):
        BinaryClass.BinaryClass(
            df = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1),
            chemical = "chemical.id",
            plate = "plate.id",
            well = "well",
            concentration = "concentration",
            endpoint = "endpoint",
            value = "Gogurt",
            format = "long"
        )
    
    # Ensure reserved value column names are rejected before any value-content checks run
    with pytest.raises(Exception, match = "bmdrc.filter is not a permitted name. Please rename this column."):

        df2 = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1)
        df2 = df2.rename({"value":"bmdrc.filter"}, axis = 1)

        BinaryClass.BinaryClass(
            df = df2,
            chemical = "chemical.id",
            plate = "plate.id",
            well = "well",
            concentration = "concentration",
            endpoint = "endpoint",
            value = "bmdrc.filter",
            format = "long"
        )

    # Ensure value columns must contain binary or NA values (other values will be rejected)
    with pytest.raises(Exception, match = "The value column must be comprised of only zeroes, ones, and NA values."):

        df3 = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1)
        df3.loc[df3.index[0], "value"] = 2

        BinaryClass.BinaryClass(
            df = df3,
            chemical = "chemical.id",
            plate = "plate.id",
            well = "well",
            concentration = "concentration",
            endpoint = "endpoint",
            value = "value",
            format = "long"
        )


# Ensure that embedded spaces get rejected after casting to strings
def test_binary_class_rejects_spaces_in_string_identifier_columns():
    df4 = pd.read_csv("data/Binary_Simplified_Long.csv").drop("Notes", axis = 1)

    # Ensure chemical validation catches embedded spaces in chemical identifiers
    chem_df = df4.copy()
    chem_df["chemical.id"] = chem_df["chemical.id"].astype(str)
    chem_df.loc[chem_df.index[0], "chemical.id"] = " 1"
    with pytest.raises(Exception, match = "spaces are not permitted in chemical names. Check this column in your dataframe."):
        BinaryClass.BinaryClass(
            df = chem_df,
            chemical = "chemical.id",
            plate = "plate.id",
            well = "well",
            concentration = "concentration",
            endpoint = "endpoint",
            value = "value",
            format = "long"
        )

    # Ensure plate validation catches embedded spaces in plate identifiers
    plate_df = df4.copy()
    plate_df["plate.id"] = plate_df["plate.id"].astype(str)
    plate_df.loc[plate_df.index[0], "plate.id"] = " A"
    with pytest.raises(Exception, match = "spaces are not permitted in plate names. Check this column in your dataframe."):
        BinaryClass.BinaryClass(
            df = plate_df,
            chemical = "chemical.id",
            plate = "plate.id",
            well = "well",
            concentration = "concentration",
            endpoint = "endpoint",
            value = "value",
            format = "long"
        )

    # Ensure well validation catches embedded spaces in well identifiers
    well_df = df4.copy()
    well_df["well"] = well_df["well"].astype(str)
    well_df.loc[well_df.index[0], "well"] = " 1"
    with pytest.raises(Exception, match = "spaces are not permitted in well names. Check this column in your dataframe."):
        BinaryClass.BinaryClass(
            df = well_df,
            chemical = "chemical.id",
            plate = "plate.id",
            well = "well",
            concentration = "concentration",
            endpoint = "endpoint",
            value = "value",
            format = "long"
        )

    # Ensure endpoint validation catches embedded spaces in endpoint identifiers
    endpoint_df = df4.copy()
    endpoint_df["endpoint"] = endpoint_df["endpoint"].astype(str)
    endpoint_df.loc[endpoint_df.index[0], "endpoint"] = "DP 24 "
    with pytest.raises(Exception, match = "spaces are not permitted in endpoint names. Check this column in your dataframe."):
        BinaryClass.BinaryClass(
            df = endpoint_df,
            chemical = "chemical.id",
            plate = "plate.id",
            well = "well",
            concentration = "concentration",
            endpoint = "endpoint",
            value = "value",
            format = "long"
        )

## Test for running BinaryClass and fit it using test data 
@pytest.fixture
def fitted():
    """Build a BinaryClass using a written example and fit it"""
    data_test = np.random.default_rng(0)            # seed = recreate data  
    rows = []
    for conc, p in zip([0,1,5,10,50], [0.0,0.1,0.3,0.6,0.9]):
        for w in range(12):
            rows.append({
                "chemical.id" : "1",
                "conc" : conc,
                "plate.id" : "A",
                "well" : f"w{conc}_w",
                "endpoint" :"DP24",
                "value" : int(data_test.random() < p), 
            })
    obj = BinaryClass.BinaryClass(
        df = pd.DataFrame(rows),
        chemical = "chemical.id", plate = "plate.id" , well = "well" ,
        concentration = "conc", endpoint = "endpoint" , value = "value" ,
        format = "long" , 
    )
    obj.make_plate_groups()
    obj.fit_models()
    return obj

# Test bmds table 
def test_fit_bmds_tables(fitted):
    assert isinstance(fitted.bmds, pd.DataFrame)
    assert not fitted.bmds.empty

# Test bmds has expected columns
def test_bmds_expected_columns(fitted):
    expected = {"bmdrc.Endpoint.ID", "Model", "BMD10", "BMDL", "BMD50"}
    assert expected.issubset(set(fitted.bmds.columns))

# Test bmds values are equal 
def test_bmd_equal_values(fitted):
    row = fitted.bmds.iloc[0]
    assert row["BMD10"] > 0
    assert row["BMD10"] <= row["BMD10"] <= row["BMD50"]
    assert 0 < row["BMD10"] <= 50

# Test output methods
def test_output_methods(fitted, tmp_path):
    fitted.output_benchmark_dose(path=str(tmp_path / "bmd.csv"))
    fitted.output_dose_table(path=str(tmp_path / "dose.csv"))
    fitted.output_fits_table(path=str(tmp_path / "fits.csv"))

# Test report runs
def test_report_runs(fitted, tmp_path):
    fitted.report(out_folder=str(tmp_path))

# Test fitted response curve 
def test_response_curve(fitted):
    fitted.response_curve(
        chemical_name="1",
        endpoint_name="DP24",
        model="quantal linear",
    )