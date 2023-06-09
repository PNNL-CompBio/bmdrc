#!/usr/bin/python3

from bmdrc.input_data_classes.BinaryClass import BinaryClass
import pandas as pd

####################
## TEST WIDE DATA ##
####################

# Wide Test
morpho_example_wide = pd.read_csv("/Users/degn400/Git_Repos/bmdrc/data/Binary_Morphology_Wide.csv")

Wide = BinaryClass(
    df = morpho_example_wide,
    chemical = "chemical.id",
    plate = "plate.id",
    well = "well",
    concentration = "conc",
    format = "wide"
)

Wide.df
Wide.format
Wide.chemical
Wide.plate
Wide.well
Wide.concentration
Wide.endpoint
Wide.value

####################
## TEST LONG DATA ##
####################

morpho_example_long = pd.read_csv("/Users/degn400/Git_Repos/bmdrc/data/Binary_Morphology_Long.csv")

Long = BinaryClass(
    df = morpho_example_long,
    chemical = "chemical.id",
    plate = "plate.id",
    well = "well",
    concentration = "conc",
    endpoint = "endpoint",
    value = "value"
)

Long.df
Long.format
Long.chemical
Long.plate
Long.well
Long.concentration
Long.endpoint
Long.value