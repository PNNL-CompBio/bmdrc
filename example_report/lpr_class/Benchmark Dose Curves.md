# Benchmark Dose Curves

## Input Data

A **lpr class** object was created. The following column names were set:

|Parameter|Column Name|
|---------|-----------|
|Chemical|chemical.id|
|Plate|plate.id|
|Well|well|
|Concentration|conc|
|Time|variable|
|Value|value|
|Cycle Length|20.0|
|Cycle Cooldown|10.0|
|Starting Cycle|light|

## Pre-Processing

#### **Combine & Make New Endpoints**
New endpoints were made using existing endpoints using 'or', which means that if there is any endpoints with a '1', this new endpoint will also have a '1', regardless of how many zeroes there are in the other endpoints. See a summary table of added endpoints below:

|New Endpoint Name|Combined Existing Endpoints|
|---|---|
|ANY_MOV|MOV1, MOV2, MOV3, MOV4|
|ANY_AUC|AUC1, AUC2, AUC3, AUC4|

#### **Set Invalid Wells to NA**

In some cases, like when a sample fish dies, many affected endpoints need to be set to NA. Here, the 'Endpoint Name' column denotes the specific endpoint that sets this rule. In this example, it could be MORT for mortality. Then, the endpoint value needs to be set, which in this case would be a 1 to indicate sample fish that did die. All endpoints would then be set to NA except for cases where the endpoint should not be affected, which are referred to as 'Endpoint Exceptions.'

This step was not conducted.

#### **Remove Invalid Endpoints**

This step was not conducted.

## Filtering

#### **Negative Control Filter**

Plates with unusually high responses in negative control samples were filtered. The response threshold was set to **50.0**. See a summary below:

|Response|Number of Plates|Filter|
|---|---|---|
|0.0|6|Keep|
|8.3333|8|Keep|
|16.6667|5|Keep|
|25.0|4|Keep|
|33.3333|2|Keep|
|41.6667|1|Keep|
|50.0|2|Filter|
|75.0|2|Filter|

And here is the plot:
![Filter Negative Control](./filter_negative_control.png)

#### **Minimum Concentration Filter**

Endpoints with too few concentration measurements (non-NA) to model are removed. The minimum was set to **3**. See a summary below:

|Number of Concentrations|Number of Endpoints|Filter|
|---|---|---|
|8|10|Keep|

And here is the plot:
![Filter Minimum Concentration](./filter_minimum_concentration.png)

#### **Correlation Score Filter**

Endpoints with little to no positive correlation with dose are unexpected and should be removed. The correlation threshold was set to **0.2**. See a summary below:

|Correlation Score Bin|Number of Endpoints|
|---|---|
|-1.0|0.0|
|-0.8|0.0|
|-0.6|2.0|
|-0.4|1.0|
|-0.2|2.0|
|0.0|3.0|
|0.2|2.0|
|0.4|0.0|
|0.6|0.0|
|0.8|0.0|

And here is the plot:
![Filter Correlation Score](./filter_correlation_score.png)

## Model Fitting & Output Modules

#### **Filter Summary**

Overall, 10 endpoint and chemical combinations were considered. 2 were deemed eligible for modeling, and 8 were not based on filtering selections explained in the previous section.

#### **Model Fitting Selections**

The following model fitting parameters were selected.

|Parameter|Value|Parameter Description|
|---|---|---|
|Goodness of Fit Threshold|0.1|Minimum p-value for fitting a model. Default is 0.1|
|Akaike Information Criterion (AIC) Threshold|2.0|Any models with an AIC within this value are considered an equitable fit. Default is 2.
|Model Selection|lowest BMDL|Either return one model with the lowest BMDL, or combine equivalent fits|

#### **Model Quality Summary**

Below is a summary table of the number of endpoints with a high quality fit (a flag of 1, meaning that the BMD10 value is within the range of measured doses) and those that are not high quality (a flag of 0).

|Flag|Count|
|---|---|
|0|0|
|1|2|

#### **Output Modules**

Below, see a table of useful methods for extracting outputs from bmdrc.

|Method|Description|
|---|---|
|.bmds|Table of fitted benchmark dose values|
|.bmds_filtered|Table of filtered models not eligible for benchmark dose calculations|
|.output_res_benchmark_dose|Table of benchmark doses for all models, regardless of whether they were filtered or not|
|.p_value_df|Table of goodness of fit p-values for every eligible endpoint|
|.aic_df|Table of Akaike Information Criterion values for every eligible endpoint|
|.response_curve|Plot a benchmark dose curve for an endpoint|
