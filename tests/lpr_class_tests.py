import pandas as pd
import bmdrc.LPRClass
import pytest
import bmdrc
from bmdrc import LPRClass

## How to calculate coverage (from within main package directory): 
# coverage run --source=bmdrc -m pytest -x tests/*
# coverage report
# coverage html

## Binary Class Tests ## 

# Test to ensure long data runs without error 
def test_lpr_class():

    # Light cycle is tested in  model_fitting_and_output_tests.py file, since LPR has slightly different reports. Here we will test dark. 
    LPR_Test = LPRClass.LPRClass(
        df = pd.read_csv("data/LPR_Long.csv"),
        chemical = "chemical.id",
        plate = "plate.id",
        well = "well",
        concentration = "conc",
        time = "variable",
        value = "value",
        cycle_length = 20.0,
        cycle_cooldown = 10.0, 
        starting_cycle = "dark"
    )
    assert isinstance(LPR_Test, bmdrc.LPRClass.LPRClass)
    assert (LPR_Test.df.columns == ['chemical.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']).all()
    
# Test wrong inputs for data.frame 
def test_df():

    # The df must be a pandas DataFrame, no exceptions
    with pytest.raises(Exception, match = "df must be a pandas DataFrame."):
        LPRClass.LPRClass(
            df = "celery",
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
        
    
    # The df must be a pandas DataFrame with data in it 
    with pytest.raises(Exception, match = "df cannot be empty. Please provide a pandas DataFrame."):
        LPRClass.LPRClass(
            df = pd.DataFrame(),
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
        

# Test wrong inputs for chemicals
def test_chemical():

    # The chemical must be a string
    with pytest.raises(Exception, match = "chemical must be a name of a column in df."):
        LPRClass.LPRClass(
            df = pd.read_csv("data/LPR_Long.csv"),
            chemical = 3, 
            plate = "plate.id",
            well = "well",
            concentration = "conc",
            time = "variable",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

    # The chemical must be a name in the dataframe 
    with pytest.raises(Exception, match = "cantelope is not in the column names of df."):
        LPRClass.LPRClass(
            df = pd.read_csv("data/LPR_Long.csv"),
            chemical = "cantelope", 
            plate = "plate.id",
            well = "well",
            concentration = "conc",
            time = "variable",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )
    
    # The chemical name must not be an unacceptable name 
    with pytest.raises(Exception, match = "bmdrc.Well.ID is not a permitted name. Please rename this column."):

        df2 = pd.read_csv("data/LPR_Long.csv")
        df2 = df2.rename({"chemical.id":"bmdrc.Well.ID"}, axis = 1)

        LPRClass.LPRClass(
            df = df2,
            chemical = "bmdrc.Well.ID", 
            plate = "plate.id",
            well = "well",
            concentration = "conc",
            time = "variable",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

# Test wrong inputs for plates
def test_plates():

    # The plate must be a string
    with pytest.raises(Exception, match = "plate must be a name of a column in df."):
        LPRClass.LPRClass(
            df =  pd.read_csv("data/LPR_Long.csv"),
            chemical = "chemical.id", 
            plate = 42, 
            well = "well",
            concentration = "conc",
            time = "variable",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

    # The plate must be a name in the dataframe 
    with pytest.raises(Exception, match = "taco salad is not in the column names of df."):
        LPRClass.LPRClass(
            df = pd.read_csv("data/LPR_Long.csv"),
            chemical = "chemical.id", 
            plate = "taco salad", 
            well = "well",
            concentration = "conc",
            time = "variable",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )
    
    # The chemical name must not be an unacceptable name 
    with pytest.raises(Exception, match = "bmdrc.num.tot is not a permitted name. Please rename this column."):

        df2 = pd.read_csv("data/LPR_Long.csv")
        df2 = df2.rename({"plate.id":"bmdrc.num.tot"}, axis = 1)

        LPRClass.LPRClass(
            df = df2,
            chemical = "chemical.id", 
            plate = "bmdrc.num.tot", 
            well = "well",
            concentration = "conc",
            time = "variable",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

# Test wrong inputs for wells
def test_wells():

    # The well must be a string
    with pytest.raises(Exception, match = "well must be a name of a column in df."):
        LPRClass.LPRClass(
            df = pd.read_csv("data/LPR_Long.csv"),
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = False, 
            concentration = "conc",
            time = "variable",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

    # The well must be a name in the dataframe 
    with pytest.raises(Exception, match = "Tuesday is not in the column names of df."):
        LPRClass.LPRClass(
            df = pd.read_csv("data/LPR_Long.csv"),
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "Tuesday", 
            concentration = "conc",
            time = "variable",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )
    
    # The well must not be an unacceptable name 
    with pytest.raises(Exception, match = "bmdrc.num.nonna is not a permitted name. Please rename this column."):

        df2 = pd.read_csv("data/LPR_Long.csv")
        df2 = df2.rename({"well":"bmdrc.num.nonna"}, axis = 1)

        LPRClass.LPRClass(
            df = df2,
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "bmdrc.num.nonna", 
            concentration = "conc",
            time = "variable",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

# Test wrong inputs for concentrations
def test_concs():

    # The concentration must be a string
    with pytest.raises(Exception, match = "concentration must be a name of a column in df."):
        LPRClass.LPRClass(
            df = pd.read_csv("data/LPR_Long.csv"),
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = 22, 
            time = "variable",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

    # The concentration must be a name in the dataframe 
    with pytest.raises(Exception, match = "Star Wars is not in the column names of df."):
        LPRClass.LPRClass(
            df = pd.read_csv("data/LPR_Long.csv"),
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = "Star Wars", 
            time = "variable",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )
    
    # The concentration name must not be an unacceptable name 
    with pytest.raises(Exception, match = "bmdrc.frac.affected is not a permitted name. Please rename this column."):

        df2 = pd.read_csv("data/LPR_Long.csv")
        df2 = df2.rename({"conc":"bmdrc.frac.affected"}, axis = 1)

        LPRClass.LPRClass(
            df = df2,
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = "bmdrc.frac.affected", 
            time = "variable",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

# Test wrong inputs for time 
def test_time():

    # The time column should be a string
    with pytest.raises(Exception, match = "time must be a name of a column in df."):
        LPRClass.LPRClass(
            df = pd.read_csv("data/LPR_Long.csv"),
            chemical = "chemical.id",
            plate = "plate.id",
            well = "well",
            concentration = "conc",
            time = 4,
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

    # The column name for time should be in the data.frame
    with pytest.raises(Exception, match = "4 is not in the column names of df."):
        LPRClass.LPRClass(
            df = pd.read_csv("data/LPR_Long.csv"),
            chemical = "chemical.id",
            plate = "plate.id",
            well = "well",
            concentration = "conc",
            time = "4",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

    # The time name must not be an unacceptable name 
    with pytest.raises(Exception, match = "cycle is not a permitted name. Please rename this column."):

        df2 = pd.read_csv("data/LPR_Long.csv")
        df2 = df2.rename({"variable":"cycle"}, axis = 1)

        LPRClass.LPRClass(
            df = df2,
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = "conc", 
            time = "cycle",
            value = "value",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

# Test wrong inputs for value
def test_value():

    # The value column should be a string
    with pytest.raises(Exception, match = "value must be a name of a column in df."):
        LPRClass.LPRClass(
            df = pd.read_csv("data/LPR_Long.csv"),
            chemical = "chemical.id",
            plate = "plate.id",
            well = "well",
            concentration = "conc",
            time = "variable",
            value = 22,
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

    # The column name for time should be in the data.frame
    with pytest.raises(Exception, match = "Toblerone is not in the column names of df."):
        LPRClass.LPRClass(
            df = pd.read_csv("data/LPR_Long.csv"),
            chemical = "chemical.id",
            plate = "plate.id",
            well = "well",
            concentration = "conc",
            time = "variable",
            value = "Toblerone",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

    # The time name must not be an unacceptable name 
    with pytest.raises(Exception, match = "bmdrc.frac.affected is not a permitted name. Please rename this column."):

        df2 = pd.read_csv("data/LPR_Long.csv")
        df2 = df2.rename({"value":"bmdrc.frac.affected"}, axis = 1)

        LPRClass.LPRClass(
            df = df2,
            chemical = "chemical.id", 
            plate = "plate.id", 
            well = "well", 
            concentration = "conc", 
            time = "variable",
            value = "bmdrc.frac.affected",
            cycle_length = 20.0,
            cycle_cooldown = 10.0, 
            starting_cycle = "light"
        )

# Test wrong inputs for the cycle_length
def test_cycle_length():

    # Cycle length must be a float 
    with pytest.raises(Exception, match = "cycle_length should be a float."):
        LPRClass.LPRClass(
                df = pd.read_csv("data/LPR_Long.csv"),
                chemical = "chemical.id",
                plate = "plate.id",
                well = "well",
                concentration = "conc",
                time = "variable",
                value = "value",
                cycle_length = "Jesse",
                cycle_cooldown = 10.0, 
                starting_cycle = "light"
            )
        
# Test wrong inputs for the cycle_cooldown
def test_cycle_cooldown():

    # Cycle cooldown must be a float 
    with pytest.raises(Exception, match = "cycle_cooldown should be a float."):
        LPRClass.LPRClass(
                df = pd.read_csv("data/LPR_Long.csv"),
                chemical = "chemical.id",
                plate = "plate.id",
                well = "well",
                concentration = "conc",
                time = "variable",
                value = "value",
                cycle_length = 20.0,
                cycle_cooldown = "grapes", 
                starting_cycle = "light"
            )
        
# Test wrong inputs for the starting_cycle
def test_starting_cycle():

    # Cycle length must be a float 
    with pytest.raises(Exception, match = "starting_cycle must be either 'light' or 'dark'."):
        LPRClass.LPRClass(
                df = pd.read_csv("data/LPR_Long.csv"),
                chemical = "chemical.id",
                plate = "plate.id",
                well = "well",
                concentration = "conc",
                time = "variable",
                value = "value",
                cycle_length = 20.0,
                cycle_cooldown = 10.0, 
                starting_cycle = "night"
            )