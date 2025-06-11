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

    # The value must be a string
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

    # The value must be a name in the dataframe 
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
    
    # The value must not be an unacceptable name 
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