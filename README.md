# bmdrc

Python library for the calculation of **B**ench**M**ark **D**ose **R**esponse **C**urves (bmdrc)

# General Schematic 

The bmdrc library was built to calculate benchmark dose (BMD) response curves for dichotomous (morphological) and light photomotor response datasets, where continuous variables are transformed to dichotomous as described in [Thomas et al 2019](https://www.sciencedirect.com/science/article/pii/S2468111318300732). Potential outputted files include a csv file of all final BMDs and their estimation errors, a csv file of model fits (AIC) for each endpoint, and an html report containing information on how much data was filtered and why, as well as interactive response curve plots. Users may specify their outputs of interest. 

![General bmdrc inputs and outputs](https://github.com/PNNL-CompBio/bmdrc/blob/main/bmdrc.png)

1. *Input Data Module:* Import data into the python library

2. *Pre-Processing Module:* Combine and remove endpoints as needed

3. *Filtering Modules:* Apply the EPA recommendations for filtering 

4. *Model Fitting Modules:* Fit EPA-recommended models to data

5. *Output Modules:* Select tables to output as csvs. View plots in a HTML report.

# How to install the package

First, install the package from pip using:

`pip install bmdrc`

Or from github using:

`pip install git+https://github.com/PNNL-CompBio/bmdrc`

# How to use the package 

### Vignettes
An example vignette for the dichotomous (binary) data (when plate and well information is available) is located [here](https://github.com/PNNL-CompBio/bmdrc/blob/main/vignettes/Binary%20Class%20Example.ipynb) and an example vignette for when only chemical, endpoint, concentration, and response (ranging from 0-1) information is available is located [here](https://github.com/PNNL-CompBio/bmdrc/blob/main/vignettes/Simplified%20Class%20Example.ipynb). For the light photomotor response data (continuous converted to dichotomous), a vignette is located [here](https://github.com/PNNL-CompBio/bmdrc/blob/main/vignettes/LPR%20Class%20Example.ipynb).
 
### Example Data 
Example data for dichotomous and light photomotor response data can be found [here](https://github.com/PNNL-CompBio/bmdrc/tree/main/data).

### Example Report
A sample for the generated report for dichotomous (binary) data can be found [here](https://github.com/PNNL-CompBio/bmdrc/blob/main/example_report/binary_class/Benchmark%20Dose%20Curves.md)

A sample for the generated report for light photomotor response data can be found [here](https://github.com/PNNL-CompBio/bmdrc/blob/main/example_report/lpr_class/Benchmark%20Dose%20Curves.md)

# Test function coverage using coverage.py

| Name                      | Stmts  | Miss | Cover |
|---------------------------|--------|------|-------|
| bmdrc/BinaryClass.py      |  138   |   0  | 100%  |
| bmdrc/LPRClass.py         |  211   |   0  | 100%  |
| bmdrc/filtering.py        |  146   |   0  | 100%  |
| bmdrc/model_fitting.py    |  687   |   5  |  99%  |
| bmdrc/output_modules.py   |  125   |   2  |  98%  |
| bmdrc/preprocessing.py    |   59   |   0  | 100%  |
| **TOTAL**                 | **1366**|**7**|**99%**|

# Data Classes

For datasets with plates and wells, use the `BinaryClass` object as it contains specific pre-processing and filtering functions that pertain to plates and well. Each value should be 0, 1, or NA, and `bmdrc` will calculate the response column. If your data does not contain plates and wells, use the `SimplifiedClass` object which only requires chemical, endpoint, concentration, and response information. The response must be between 0 and 1. If your data is light photomotor response use the `LPRClass` object. The continuous data will be transformed to dichotomous.