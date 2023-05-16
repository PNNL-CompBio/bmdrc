# bmdrc

Python library for the calculation of *b*ench*m*ark *d*ose *r*esponse *c*urves (bmdrc)

# General Schematic 

The main bmdrc function was built to calculate benchmark doses (BMDs) for dichotomous (morphological, light photometer response, or movement response) and continuous 
(transcriptomics) datasets. Potential outputted files include a csv file of all final BMDs and their estimation errors, a csv file of model fits (AIC) for each endpoint, 
and an html report containing information on how much data was filtered and why, as well as interactive response curve plots. Users may specify their outputs of interest. 

![General bmdrc inputs and outputs](./bmdrc.png)

# bmdrc Function Parameters 

**Required Parameters**

*filepath:* A path to a .csv or .txt file (see example input data for acceptable choices)

*datatype:* Specify which of the four input data files the data is

**Optional Parameters**

*format:* Specify whether the data is written in "long" or "wide" format. The default is "long". If in long format, the following columns are required: 'chemical.id', 'conc', 
'plate.id', 'well', 'endpoint', and 'value'. If the data is in wide format, the following columns are required: 'chemical.id', 'conc', 'plate.id', and 'well'. The remaining columns 
are assumed to be endpoints with their respective values. 

*morpho_filepath:* If the input data is not morphological, users may still specify the morphological filepath to filter out specific cases. See "endpoints_to_filter" parameter.

*endpoints_to_filter:* A list of endpoints to filter out everyone. Default is "DNC" (do not count), "DNC_", "MORT" (any mortality), and "MO24" (mortality at 24 hours).  

*endpoints_to_combine:* A dictionary of new endpoints that can be calculated from existing endpoints. For example `{'ANY': ['MO24','DP24','SM24']}`. Default is None.

*measurement_minimum:* The minimum number of measurements to calculate a benchmark dose curve. Default is 3. Users cannot specify a value less than 2. 

*report_options:* A list to specify report information. Default is "missingness" and "response curves" written as `["missingness", "response curves"]`

*outputs:* A list to specify outputs. Default is "benchmark", "model fits", and "report" written as `["benchmark", "model fits", "report"]`

# Example Input Data

Example morphological, light photometer response, movement response, and transcriptomics data are included within the package. See >to do< for more details. 
