# bmdrc

Python library for the calculation of **B**ench**M**ark **D**ose **R**esponse **C**urves (bmdrc)

# General Schematic 

The main bmdrc function was built to calculate benchmark dose (BMD) response curves for dichotomous (morphological, light photometer response, or movement response) and continuous 
(transcriptomics) datasets. Potential outputted files include a csv file of all final BMDs and their estimation errors, a csv file of model fits (AIC) for each endpoint, 
and an html report containing information on how much data was filtered and why, as well as interactive response curve plots. Users may specify their outputs of interest. 

![General bmdrc inputs and outputs](./bmdrc.png)

# bmdrc Function Parameters 

**Required Parameters**

*filepath:* A path to a .csv or .txt file (see example input data for acceptable choices)

*datatype:* Specify input data type. Acceptable options are "morphological" (dichotomous), "lpr" (light photometer response), "movement" (behavorial), and "omics" (transcriptomics, metabolomics, lipidomics,
or proteomics)

**Conditionally Required Parameters**

*omics_metadata:* If the datatype is "omics", this file requires meta information about the dosage, and the column names in the file in the filepath that contain the fold change values and the p-values. 
The required column names in omics_metadata are "dosage", "foldchange", and "pvalue". 

*format:* If the datatype is "morphological", specify whether the data is written in "long" or "wide" format. The default is "long". If in long format, the following columns are required: 
'chemical.id', 'conc', 'plate.id', 'well', 'endpoint', and 'value'. If the data is in wide format, the following columns are required: 'chemical.id', 'conc', 'plate.id', and 'well'. The remaining columns 
are assumed to be endpoints with their respective values. 

**Optional Parameters**

*morpho_filepath:* If the input data is not morphological, users may still specify the morphological filepath to filter out specific cases. See "endpoints_to_filter" parameter.

*endpoints_to_filter:* A list of endpoints to filter out everywhere. Default is "DNC" (do not count), "DNC_", "MORT" (any mortality), and "MO24" (mortality at 24 hours) written as 
`["DNC", "DNC_", "MORT", "MO24"]`  

*endpoints_to_combine:* A dictionary of new endpoints that can be calculated from existing endpoints. For example `{'ANY': ['MO24','DP24','SM24']}`. Default is None.

*bmds:* A list of numeric BMD values to calculate. Even if this parameter is set to None, BMD10 will always be calculated. For example, `[30, 50]` will calculate BMD10, BMD30, and BMD50. 

*minimum_measurement_filter:* The minimum number of measurements to calculate a benchmark dose curve. Default is 3. Users cannot specify a value less than 2. For more details on filters, 
[see here](https://www.nature.com/articles/s41597-023-02021-5)

*negative_control_filter:* The maximum percentage response in the controls. Higher percentages are filtered. Default is 50, which is 50%. For more details on filters, 
[see here](https://www.nature.com/articles/s41597-023-02021-5)

*correlation_filter:* Response must increase with dose. A spearman correlation metric (ranges from -1 to 1) is calculated and curves that fall below the specified threshold are filtered out. Default is 0.2 
For more details on filters, [see here](https://www.nature.com/articles/s41597-023-02021-5)

*equivalent_fit_threshold:* The AIC value threshold for an equivalent fit. Default is 2. [See here](https://doi.org/10.1177/0049124104268644) 

*report_options:* A list to specify report information. Default is "missingness" and "response curves" written as `["missingness", "response curves"]`

*outputs:* A list to specify outputs. Default is "benchmark", "model fits", and "report" written as `["benchmark", "model fits", "report"]`

# Fit Selection

Models with the lowest AIC value are generally considered the optimal fit. In cases where there are models with equivalent fits (see the equivalent_fit_threshold parameter), the 
model with the smallest BMD10 estimation error values are then considered the optimal fit.  

# Accessory Functions

There are a few additional functions besides the main `bmdrc()` function that can be used to extract additional information. Brief descriptions of these functions are listed below:

**predict_response:** Users interested in predicting response curves can supply the "benchmark" file, a string with the endpoint to calculate, and a list of dose values (x variables) to predict on.

# Example Input Data

Example morphological, light photometer response, movement response, and transcriptomics data are included within the package. See >to do< for more details. 
