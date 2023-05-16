# bmdrc

Python library for the calculation of benchmark dose response curves 

# General Schematic 

The main bmdrc function was built to calculate benchmark doses (BMDs) for dichotomous (morphological, light photometer response, or movement response) and continuous 
(transcriptomics) datasets. Potential outputted files include a csv file of all final BMDs and their estimation errors, a csv file of model fits (AIC) for each endpoint, 
and an html report containing information on how much data was filtered and why, as well as interactive response curve plots. Users may specify their outputs of interest. 

![General bmdrc inputs and outputs](./bmdrc.png)
