#################
Welcome to bmdrc!
#################

*bmdrc* is a python package for fitting benchmark dose curves to dichotomous, proportional, and light photomotor response data.
This package is a statistics toolkit which can be followed in 5 main steps: `1. Upload Class Options`_, 
`2. Preprocessing Options`_, 

**NOTE!** Backend function names do not match the frontend function names in objects. Pay close attention to the example code and 
vignettes, instead of the backend function descriptions. Thank you.

#######################
1. Upload Class Options
#######################

Upload data using one of the following class options, depending on your input data format:

`Binary Class`_ is for binary data that needs to be converted to propotions. 
`LPR Class`_ is for continuous larval photomotor response data that needs to be converted to proportions. 
`Simplified Class`_ is for data that is already formatted as proportions.

************
Binary Class
************

.. autoclass:: bmdrc.BinaryClass.BinaryClass

.. code-block:: python

   # Long format
   BinaryClass(
      df = pd.read_csv("path/to/longfile.csv"), # Input is a pandas DataFrame
      chemical = "chemical.id", # The name of the chemical column 
      plate = "plate.id", # The name of the plate ID column
      well = "well", # The name of the column with well names
      concentration = "concentration", # The name of the concentration column
      endpoint = "endpoint", # The name of the column with endpoints
      value = "value", # The name of the column with values
      format = "long" # The format of the input data, either 'long' or 'wide' is accepted
   )

   # Wide format
   BinaryClass(
      df = pd.read_csv("path/to/widefile.csv"),
      chemical = "chemical.id",
      plate = "plate.id",
      well = "well",
      concentration = "conc",
      endpoint = "endpoint",
      value = "value",
      format = "wide"
   )

*********
LPR Class
*********

.. autoclass:: bmdrc.LPRClass.LPRClass

.. code-block:: python

   # Convert the continuous data to dichotomous 
   LPRClass(
      df = pd.read_csv("path/to/lpr.csv"),
      chemical = "chemical.id", # Column in file
      plate = "plate.id", # Column in file
      well = "well", # Column in file
      concentration = "conc", # Column in file
      time = "variable", # Column in file
      value = "value", # Column in file
      cycle_length = 20.0, # Length of cycle in 6 second intervals. 20 * 6 = 120 seconds
      cycle_cooldown = 10.0, # Length of cycle in 6 second intervals. 10 * 6 = 60 seconds
      starting_cycle = "light" # Starting cycle
   )

****************
Simplified Class
****************

.. autoclass:: bmdrc.SimplifiedClass.SimplifiedClass

.. code-block:: python

   SimplifiedClass(
      df = pd.read_table("path/to/proportions.csv"), # Input is a pandas DataFrame
      chemical = "chemical.id", # The name of the chemical column 
      endpoint = "endpoint", # The name of the column with endpoints
      concentration = "concentration", # The name of the concentration column
      response = "response" # The name of the column with response value ranging from 0 to 1
   )

#########################
2. Preprocessing Options
#########################

There are currently three preprocessing options including: `Combining Endpoints`_, `Removing Endpoints`_, and `Removing Wells`_.

*******************
Combining Endpoints 
*******************

.. autoclass:: bmdrc.preprocessing.endpoint_combine

.. code-block:: python

   # Dictionary of terms to add
   endpoint_dict = {"ANY24":["NC24", "DP24", "SM24"], "ANY":["NC24", "DP24", "SM24", "JAW"]}

   # Create a bmdrc object and save it as Long. See the vignettes. Add new endpoint
   Long.combine_and_create_new_endpoints(endpoint_dict)

*******************
Removing Endpoints 
*******************

.. autoclass:: bmdrc.preprocessing.remove_endpoints

.. code-block:: python

   # Create a bmdrc object and save it as Long. See the vignettes
   # Remove the endpoint that should not be modeled
   Long.remove_endpoints("DNC")

**************
Removing Wells
**************

.. autoclass:: bmdrc.preprocessing.well_to_na

.. code-block:: python

   # Create a bmdrc object and save it as Long. See the vignettes
   Long.set_well_to_na(endpoint_name = "DNC", endpoint_value = 1, except_endpoint = ["ANY24"])