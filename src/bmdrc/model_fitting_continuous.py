import operator
from abc import abstractmethod

import warnings
warnings.filterwarnings("ignore")

from scipy.stats import chi2
from scipy import stats
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from uncertainties import correlated_values
import uncertainties.unumpy as unp 

import numpy as np
import pandas as pd

############################
## CLASSES FOR MODEL FITS ##
############################

## General Class Function ##
class Continuous_Model():
        
    ## Methods which do not change per model ##
    
    @abstractmethod
    def calculate_aic(self, y_obs: np.array, y_pred: np.array, params: np.array):
        '''
        Returns an AIC (Akaike Information Criterion) value for a fit dataset
        
        Parameters
        ----------
        y_obs
            the observed y (Response) variables as a numpy array
        y_pred
            the predicted y (Response) variables as a numpy array
        params
            a numpy array of the parameters

        Returns
        -------
        a single AIC value stored in the object
        '''

        # Number of observations
        n = len(y_obs)
        
        # Residuals
        residuals = y_obs - y_pred
        
        # Variance of residuals (sigma^2)
        sigma2 = np.sum(residuals**2) / n  # Mean squared error (MSE)
        
        # Log-likelihood
        log_likelihood = -n / 2 * np.log(2 * np.pi * sigma2) - np.sum(residuals**2) / (2 * sigma2)
        
        # AIC formula
        aic = 2 * len(params) - 2 * log_likelihood
        
        self.aic = aic

    @abstractmethod
    def get_response_level(self, percentage: float):
        '''
        Calculate the y value (Response) at a percentage

        Parameters
        ----------
        percentage
            the percentage of the response, where 10 means 10%

        Returns
        -------
        the actual response value at that percentage
        '''

        # Check for y_pred
        if self.y_pred is None:
            raise TypeError("Please run the .fit() function first to get a response level.")
        
        # Pull out the increments
        increment = (np.max(self.y_pred) - np.min(self.y_pred)) / 100

        # Return the predicted value
        return np.min(self.y_pred) + (increment * percentage)
    
    @abstractmethod
    def calc_conf_interv(self):
        '''
        Calculate a 95% confidence interval

        Returns
        -------
        a data.frame called "confidence.interval" store as "CI"
        '''

        # Build a statistical summary data.frame to calculate confidence intervals. SEM - standard error measurement (sd / sqrt(n))
        stats_sum = self._toModel[[self._concentration, self._response]].groupby(self._concentration).agg(["mean", "sem", "size"]).reset_index()
        stats_sum.columns = [self._concentration, "mean", "sem", "dof"]

        # Subtract one for the DoF
        stats_sum["dof"] = stats_sum["dof"] - 1

        # Store low/high variables 
        Low = []
        High = []

        # Iterate through and calculate confidence intervals
        for x in range(len(stats_sum)):
            CI = stats.t.interval(0.95, stats_sum["dof"][x], loc = stats_sum["mean"][x], scale = stats_sum["sem"][x])
            Low.append(np.round(CI[0], 4))
            High.append(np.round(CI[1], 4))

        # Add the low and high 
        stats_sum["Low"] = Low
        stats_sum["High"] = High

        # Return confidence intervals
        self.CI = stats_sum
            
    @abstractmethod
    def gen_uneven_spacing(self, doses: list[str], int_steps: int = 10):
        '''
        Support function to calculate inbetween concentration measurements
        
        Parameters
        ----------
        doses
            a list of doses, typically ranging from the min to the max dose, inclusive

        int_steps
            define the number of measurements between doses

        Returns
        -------
        a vector of values with `int_step` number of doses between measurements
        '''
        dose_samples = list()
        for dose_index in range(len(doses) - 1):
            dose_samples.extend(np.linspace(doses[dose_index],doses[dose_index + 1], int_steps).tolist())
        return np.unique(dose_samples)

    ## Attributes used by all models ##

    def __init__(self, toModel, concentration, response):
        self.toModel = toModel
        self.concentration = concentration
        self.response = response

    toModel = property(operator.attrgetter('_toModel'))
    concentration = property(operator.attrgetter('_concentration'))
    response = property(operator.attrgetter('_response'))

    @toModel.setter
    def toModel(self, theModelData):
        self._toModel = theModelData

    @concentration.setter
    def concentration(self, concentrationname):
        self._concentration = concentrationname

    @response.setter
    def response(self, responsename):
        self._response = responsename

## Linear Regression ##
class LinReg_Cont(Continuous_Model):

    def fit(self, fixed_intercept: float = 0):
        '''
        Fit a linear regression model, calculate GOF and AIC

        Parameters
        ----------
        fixed_intercept
            a value that the the modeled line must run through
        
        Returns
        ------
        the fitted model, parameters, y_pred, GOF p-value, and AIC. All return in the object.
        '''

        # Save the selected fixed_intercept
        self.fixed_intercept = fixed_intercept

        # Extract out the values
        x = self._toModel[self._concentration].to_numpy()
        y = self._toModel[self._response].to_numpy()

        # Adjust by intercept
        y = [val - fixed_intercept for val in y]

        # Fit linear model without intercept
        model = LinearRegression(fit_intercept=False)
        model.fit(x.reshape(-1, 1), y)
        self.model = model
        self.y_pred = model.predict(x.reshape(-1, 1)) + fixed_intercept
        self.params = model.coef_

        # Get the AIC value
        self.calculate_aic(y, self.y_pred, self.params)

    def predict_x(self, y):
        '''
        Predict the value of x (Dose) to achieve a specific Y (Response)

        Parameters
        ----------
        y
            the continuous response variable

        Returns
        -------
        an estimate of the x (Dose) at that response        
        '''

        # The params attributes are needed
        if self.params is None:
            raise TypeError("Please run the .fit() function first to get a response level.")
        
        # Return the estimate 
        x = (y - self.fixed_intercept) / self.params[0]
        return x if x != 0 else np.nan
    
    def response_curve(self, steps = 10):
        '''
        Create a response curve
        
        Parameters
        ----------
        steps
            the number of modeled points inbetween existing points
        
        Returns
        -------
        a curve pandas DataFrame with the Dose and Response
        '''

        # Determine x values to use
        dose_x_vals = np.round(self.gen_uneven_spacing(self._toModel[self._concentration].tolist(), int_steps = steps), 4)

        # Define curve and its columns
        curve = pd.DataFrame([dose_x_vals, self.model.predict(dose_x_vals.reshape(-1, 1))]).T
        curve.columns = ["Dose in uM", "Response"]
        curve["Response"] = curve["Response"] + self.fixed_intercept
        
        # Save the curve
        self.curve = curve

    # Calculate the bmdl
    def calc_bmdl(self):
        '''Not currently supported'''
        return np.nan

## Quadratic, Cubic, and Quartic Regression ##
class PolyReg_Cont(Continuous_Model):
    
    def fit(self, fixed_intercept: float = 0):
        '''
        Fit a polynomial regression model, calculate GOF and AIC

        Parameters
        ----------
        fixed_intercept
            a value that the the modeled line must run through
        
        Returns
        ------
        the fitted model, parameters, y_pred, GOF p-value, and AIC. All return in the object.
        '''

        # Save the selected fixed_intercept
        self.fixed_intercept = fixed_intercept

        # Extract out the values
        x = self._toModel[self._concentration].to_numpy()
        y = self._toModel[self._response].to_numpy()

        # Adjust by intercept
        y = [val - fixed_intercept for val in y]

        # Fit linear model without intercept and without the additional fits (include_bias - no need for intercept term)
        model = make_pipeline(PolynomialFeatures(self.degree, include_bias = False), LinearRegression(fit_intercept = False))
        model.fit(x.reshape(-1, 1), y)
        self.model = model
        self.y_pred = model.predict(x.reshape(-1, 1)) + fixed_intercept
        self.params = model.named_steps["linearregression"].coef_

        # Get the AIC value
        self.calculate_aic(y, self.y_pred, self.params)

    def predict_x(self, y):
        '''
        Predict the value of x (Dose) to achieve a specific Y (Response)

        Parameters
        ----------
        y
            the continuous response variable

        Returns
        -------
        an estimate of the x (Dose) at that response        
        '''

        # The params attributes are needed
        if self.params is None:
            raise TypeError("Please run the .fit() function first to get a response level.")
        
        # Calculate the roots with numpy - step 1: format
        ordered_params = np.append(np.flip(self.params), y - self.fixed_intercept)

        # Calculate the roots with numpy - step 2: use np roots
        roots = np.roots(ordered_params)

        # Calculate the roots with numpy - step 3: find the smallest non-negative/non-zero root
        smallest_root = np.min([root for root in roots if root > 0])
        return smallest_root
    
    def response_curve(self, steps = 10):
        '''
        Create a response curve
        
        Parameters
        ----------
        steps
            the number of modeled points inbetween existing points
        
        Returns
        -------
        a curve pandas DataFrame with the Dose and Response
        '''

        # Determine x values to use
        dose_x_vals = np.round(self.gen_uneven_spacing(self._toModel[self._concentration].tolist(), int_steps = steps), 4)

        # Define curve and its columns
        curve = pd.DataFrame([dose_x_vals, self.model.predict(dose_x_vals.reshape(-1, 1))]).T
        curve.columns = ["Dose in uM", "Response"]
        curve["Response"] = curve["Response"] + self.fixed_intercept
        
        # Save the curve
        self.curve = curve
    
    # Calculate the bmdl
    def calc_bmdl(self):
        '''Not currently supported'''
        return np.nan
    
    # Expand the init function definition to include the degree of the polynomial
    def __init__(self, toModel, concentration, response, degree):
        
        # toModel, concentration, and response have already been calculated in another function
        super().__init__(toModel, concentration, response)
        self.degree = degree

    degree = property(operator.attrgetter('_degree'))

    @degree.setter
    def degree(self, degree_val):
        self._degree = degree_val

## General Non-Linear Regression ## 
class GenReg_Cont(Continuous_Model):

    def fit(self, fixed_intercept: float = 0):
        '''
        Fit a regression model, calculate GOF and AIC

        Parameters
        ----------
        fixed_intercept
            a value that the the modeled line must run through
        
        Returns
        ------
        the fitted model, parameters, y_pred, GOF p-value, and AIC. All return in the object.
        '''

        # Save the selected fixed_intercept
        self.fixed_intercept = fixed_intercept

        # Extract out the values
        x = self._toModel[self._concentration].to_numpy()
        y = self._toModel[self._response].to_numpy()
        
        # Fit curve - there is no model object to keep. Keep predictions, params, and covariance matrix
        params, cov = curve_fit(self.model_equation, x, y)
        self.y_pred = self.model_equation(x, *params)
        self.params = params
        self.cov = cov

        # Get the AIC value
        self.calculate_aic(y, self.y_pred, self.params)
    
    def response_curve(self, steps = 10):
        '''
        Create a response curve
        
        Parameters
        ----------
        steps
            the number of modeled points inbetween existing points
        
        Returns
        -------
        a curve pandas DataFrame with the Dose and Response
        '''

        # Determine x values to use
        dose_x_vals = np.round(self.gen_uneven_spacing(self._toModel[self._concentration].tolist(), int_steps = steps), 4)

        # Define curve and its columns
        curve = pd.DataFrame({"Dose in uM": dose_x_vals, "Response": [self.model_equation(the_x, *self.params) for the_x in dose_x_vals]})
        
        # Save the curve
        self.curve = curve

## Asymptotic Regression ##
class AsyReg_Cont(GenReg_Cont):

    # Define the model fitting function based on the fixed intercept
    def model_equation(self, x, a, b):
        '''Internal function to fit the curve'''
        return self.fixed_intercept + a * (1 - np.exp(-b * x))
    
    def predict_x(self, y):
        '''
        Predict the value of x (Dose) to achieve a specific Y (Response)

        Parameters
        ----------
        y
            the continuous response variable

        Returns
        -------
        an estimate of the x (Dose) at that response        
        '''
        # return -1 / b * np.log(1 - ((y-100) / a))
        x = -1 / self.params[1] * np.log(1 - ((y-self.fixed_intercept) / self.params[0]))
        return x if x != 0 else np.nan

    # Calculate the bmdl
    def calc_bmdl(self):
        '''Calculate the lowest effect benchmark dose. In this case, it's the lower bound of BMD10'''
        try:
            a, b = correlated_values(self.params, self.cov)
            response = self.fixed_intercept + a * (1 - unp.exp(-b * 0.1))
            BMDL = self.predict_x(self.get_response_level(10) - unp.std_devs(response))
            return BMDL if BMDL >= 0 else np.nan
        except:
            return np.nan
    
## Exponential Regression ## 
class ExpReg_Cont(GenReg_Cont):

    # Update the model equation
    def model_equation(self, x, a, b):
        '''Internal function to fit the curve'''
        return self.fixed_intercept + (a * (np.exp(b * x) - 1))
    
    # Update the x prediction
    def predict_x(self, y):
        '''
        Predict the value of x (Dose) to achieve a specific Y (Response)

        Parameters
        ----------
        y
            the continuous response variable

        Returns
        -------
        an estimate of the x (Dose) at that response        
        '''
        # return ln((y - 100 / a) + 1) / b
        x = np.log(((y - self.fixed_intercept) / self.params[0]) + 1) / self.params[1]
        return x if x != 0 else np.nan
    
    # Calculate the bmdl
    def calc_bmdl(self):
        '''Calculate the lowest effect benchmark dose. In this case, it's the lower bound of BMD10'''
        try:
            a, b = correlated_values(self.params, self.cov)
            response = self.fixed_intercept + (a * (unp.exp(b * 0.1) - 1))
            BMDL = self.predict_x(self.get_response_level(10) - unp.std_devs(response))
            return BMDL if BMDL >= 0 else np.nan
        except:
            return np.nan

## Gompertz Regression ##
class GomReg_Cont(GenReg_Cont):

    # Update the model equation
    def model_equation(self, x, a, b, c):
        '''Internal function to fit the curve'''
        return self.fixed_intercept + a * np.exp(-1 * b * np.exp(-1 * c * x))
    
    # Update the x prediction
    def predict_x(self, y):
        '''
        Predict the value of x (Dose) to achieve a specific Y (Response)

        Parameters
        ----------
        y
            the continuous response variable

        Returns
        -------
        an estimate of the x (Dose) at that response        
        '''
        # return -1/c * ln((ln((y-100)/a)) / -1*b)
        x = -1 / self.params[2] * np.log(np.log((y - 100) / self.params[0]) / (-1 * self.params[1]))
        return x if x != 0 else np.nan
    
    # Calculate the bmdl
    def calc_bmdl(self):
        '''Calculate the lowest effect benchmark dose. In this case, it's the lower bound of BMD10'''
        try:
            a, b, c = correlated_values(self.params, self.cov)
            response = self.fixed_intercept + a * unp.exp(-1 * b * unp.exp(-1 * c * 0.1))
            BMDL = self.predict_x(self.get_response_level(10) - unp.std_devs(response))
            return BMDL if BMDL >= 0 else np.nan
        except:
            return np.nan

## Hill Regression ## 
class HillReg_Cont(GenReg_Cont):

    # Update the model equation
    def model_equation(self, x, a, b):
        '''Internal function to fit the curve'''
        return self.fixed_intercept + ((self.vmax * (x**a)) / (b**a + x**a))
    
    # Update the x prediction
    def predict_x(self, y):
        '''
        Predict the value of x (Dose) to achieve a specific Y (Response)

        Parameters
        ----------
        y
            the continuous response variable

        Returns
        -------
        an estimate of the x (Dose) at that response        
        '''

        # Calculate the y adjustment
        yadj = (y - self.fixed_intercept) / self.vmax
        a = self.params[0]
        b = self.params[1]

        # Return the value ((yadj * b^a) / (1 - yadj))^(1/a)
        x = ((yadj * b**a) / (1 - yadj))**(1/a)
        return x if x != 0 else np.nan
    
    # Calculate the bmdl
    def calc_bmdl(self):
        '''Calculate the lowest effect benchmark dose. In this case, it's the lower bound of BMD10'''
        try:
            a, b = correlated_values(self.params, self.cov)
            response = self.fixed_intercept + ((self.vmax * (0.1**a)) / (b**a + 0.1**a))
            BMDL = self.predict_x(self.get_response_level(10) - unp.std_devs(response))
            return BMDL if BMDL >= 0 else np.nan
        except:
            return np.nan
    
    # Define the vmax value - the highest y value
    def __init__(self, toModel, concentration, response):
        
        # toModel, concentration, and response have already been calculated in another function
        super().__init__(toModel, concentration, response)
        self.vmax = np.max(toModel[response])

## Michaelis-Mentin Regression ## 
class MMReg_Cont(GenReg_Cont):

    # Update the model equation
    def model_equation(self, x, a):
        '''Internal function to fit the curve'''
        return self.fixed_intercept + ((self.vmax * x) / (a + x))
    
    # Update the x prediction
    def predict_x(self, y):
        '''
        Predict the value of x (Dose) to achieve a specific Y (Response)

        Parameters
        ----------
        y
            the continuous response variable

        Returns
        -------
        an estimate of the x (Dose) at that response        
        '''

        # Return -(y-i)a / (y-i)-v
        x = (-1 * (y - self.fixed_intercept) * self.params[0]) / ((y - self.fixed_intercept) - self.vmax)
        return x if x != 0 else np.nan
    
    # Calculate the bmdl
    def calc_bmdl(self):
        '''Calculate the lowest effect benchmark dose. In this case, it's the lower bound of BMD10'''
        try:
            a = correlated_values(self.params, self.cov)
            response = self.fixed_intercept + ((self.vmax * 0.1) / (a + 0.1))
            BMDL = self.predict_x(self.get_response_level(10) - unp.std_devs(response))
            return BMDL if BMDL >= 0 else np.nan
        except:
            return np.nan
    
    # Define the vmax value - the highest y value
    def __init__(self, toModel, concentration, response):
        
        # toModel, concentration, and response have already been calculated in another function
        super().__init__(toModel, concentration, response)
        self.vmax = np.max(toModel[response])

## Power Regression ##
class PowReg_Cont(GenReg_Cont):

    # Update the model equation
    def model_equation(self, x, a, b):
        '''Internal function to fit the curve'''
        return self.fixed_intercept +  a * x**b
    
    # Update the x prediction
    def predict_x(self, y):
        '''
        Predict the value of x (Dose) to achieve a specific Y (Response)

        Parameters
        ----------
        y
            the continuous response variable

        Returns
        -------
        an estimate of the x (Dose) at that response        
        '''
        # return ((y-100)/a)^1/b
        x = ((y-self.fixed_intercept)/self.params[0])**(1/self.params[1])
        return x if x != 0 else np.nan
    
    # Calculate the bmdl
    def calc_bmdl(self):
        '''Calculate the lowest effect benchmark dose. In this case, it's the lower bound of BMD10'''
        try:
            a, b = correlated_values(self.params, self.cov)
            response = self.fixed_intercept + a * 0.1**b
            BMDL = self.predict_x(self.get_response_level(10) - unp.std_devs(response))
            return BMDL if BMDL >= 0 else np.nan
        except:
            return np.nan
    
## Weibull Regression
class WeiReg_Cont(GenReg_Cont):

    # Update the model equation
    def model_equation(self, x, a, b, c):
        '''Internal function to fit the curve'''
        return self.fixed_intercept + (a * (1 - np.exp(-(x / b)**c)))
    
    # Update the x prediction
    def predict_x(self, y):
        '''
        Predict the value of x (Dose) to achieve a specific Y (Response)

        Parameters
        ----------
        y
            the continuous response variable

        Returns
        -------
        an estimate of the x (Dose) at that response        
        '''

        a = self.params[0]
        b = self.params[1]
        c = self.params[2]

        # return b(-ln(1-(y-100)/a))^1/c
        p1 = 1 - ((y - self.fixed_intercept) / a)
        x = b * (-1 * np.log(p1))**(1/c)
        return x if x != 0 else np.nan
    
    # Calculate the bmdl
    def calc_bmdl(self):
        '''Calculate the lowest effect benchmark dose. In this case, it's the lower bound of BMD10'''
        try:
            a, b, c = correlated_values(self.params, self.cov)
            response = self.fixed_intercept + (a * (1 - unp.exp(-(0.1 / b)**c)))
            BMDL = self.predict_x(self.get_response_level(10) - unp.std_devs(response))
            return BMDL if BMDL >= 0 else np.nan
        except:
            return np.nan

#############################
## MODEL FITTING FUNCTIONS ##
#############################

def _removed_endpoints_stats(self):
    '''
    Accessory function to fit_the_models. 
    As the first in the pipeline, this function calculates summary
    statistics for the endpoints that are filtered out. No models are
    fit in these values. 
    '''

    if any(self.plate_groups["bmdrc.filter"] == "Remove"):

        # Make a data frame with all filtered endpoints called low quality and group it 
        low_quality = self.plate_groups[self.plate_groups["bmdrc.filter"] == "Remove"].groupby("bmdrc.Endpoint.ID")

        # Calculate the area under the curve (AUC) and the min and max dose. Model, BMD10, BMDL, and BMD50 are all NA. 
        bmds_filtered = low_quality.apply(lambda df: np.trapezoid(df[self._response], x = df[self._concentration])).reset_index().rename(columns = {0: "AUC"})
        bmds_filtered[["Model", "BMD10", "BMDL", "BMD50"]] = np.nan
        bmds_filtered["Min_Dose"] = low_quality[["bmdrc.Endpoint.ID", self._concentration]].min(self._concentration).reset_index()[self._concentration]
        bmds_filtered["Max_Dose"] = low_quality[["bmdrc.Endpoint.ID", self._concentration]].max(self._concentration).reset_index()[self._concentration]

        # Calculate the total area
        bmds_filtered["Max_Response"] = low_quality[["bmdrc.Endpoint.ID", self._response]].max(self._response).reset_index()[self._response]
        bmds_filtered["Area"] = (bmds_filtered["Max_Dose"] - bmds_filtered["Min_Dose"]) * bmds_filtered["Max_Response"]

        # Normalize the AUC by the area
        bmds_filtered["AUC_Norm"] = round(bmds_filtered["AUC"] / bmds_filtered["Area"], 8)
        bmds_filtered["AUC"] = round(bmds_filtered["AUC"], 4)

        # Order columns
        self.bmds_filtered = bmds_filtered[["bmdrc.Endpoint.ID", "Model", "BMD10", "BMDL", "BMD50", "AUC", "Min_Dose", "Max_Dose", "AUC_Norm"]]

    else:

        self.bmds_filtered = None

def fit_continuous_models(self, aic_threshold: float, model_selection: str, diagnostic_mode: bool):
    '''
    Fit continuous monotonic models to your dataset. 

    Parameters
    ----------    
    aic_threshold
        A float for the Akaike Information Criterion (AIC) threshold. The default is 2.

    model_selection
        A string for the model_selection model. Currently, only "lowest BMDL" is supported.

    diagnostic_mode
        A boolean to indicate whether diagnostic messages should be printed. Default is False
    '''

    # Save parameters 
    self.model_fitting_aic_threshold = aic_threshold
    self.model_fitting_model_selection = model_selection

    # Pull dose_response
    dose_response = self.plate_groups[self.plate_groups["bmdrc.filter"] == "Keep"]

    # Determine list of endpoints
    endpoints = dose_response["bmdrc.Endpoint.ID"].unique().tolist()

    # Calculate stats for endpoints that will not be fit to models
    _removed_endpoints_stats(self)

    ####################
    ## FIT ALL MODELS ##
    ####################

    # Store AICs, BMDLs, BMD10s, and BMD50s
    AICs = []
    BMDLs = []
    BMD10s = []
    BMD50s = []

    # Iterate through each endpoint
    for endpoint in endpoints:

        # Add a message if diagnostic mode
        if diagnostic_mode:
            print(".....fitting models for", endpoint)

        # Pull a specific dataset
        sub_data = dose_response[dose_response["bmdrc.Endpoint.ID"] == endpoint]

        # Calculate the AIC, BMDL, BMD10, and BMD50
        def fit_stats(mod_type):
            '''Wrapper function for fitting models and calculating statistics'''
            try:
                fit_mod = mod_type(sub_data, self._concentration, self._response)
                fit_mod.fit(fixed_intercept = 100)
                return [fit_mod.aic, fit_mod.calc_bmdl(), fit_mod.predict_x(fit_mod.get_response_level(10)), fit_mod.predict_x(fit_mod.get_response_level(50))]
            except:
                return [np.nan, np.nan, np.nan, np.nan]

        # Calculate the model metrics
        fitted_stats = [fit_stats(LinReg_Cont), fit_stats(AsyReg_Cont), fit_stats(ExpReg_Cont), 
                        fit_stats(GomReg_Cont), fit_stats(HillReg_Cont), fit_stats(MMReg_Cont),
                        fit_stats(PowReg_Cont), fit_stats(WeiReg_Cont)]

        # Unpack each of the statistics and append
        AIC_vals, BMDL_vals, BMD10_vals, BMD50_vals = [endpoint], [endpoint], [endpoint], [endpoint]
        AIC_vals.extend([metrics[0] for metrics in fitted_stats])
        AICs.append(AIC_vals)
        BMDL_vals.extend([metrics[1] for metrics in fitted_stats])
        BMDLs.append(BMDL_vals)
        BMD10_vals.extend([metrics[2] for metrics in fitted_stats])
        BMD10s.append(BMD10_vals)
        BMD50_vals.extend([metrics[3] for metrics in fitted_stats])
        BMD50s.append(BMD50_vals)

        # Set the row order
        column_names = ["bmdrc.Endpoint.ID", "linear", "asymptotic", "exponential", "gompertz", "hill", "michaelis-mentin", "power", "weibull"]

        # Save data.frame of calculated values
        self.AIC_df = pd.DataFrame(AICs, columns = column_names)
        self.BMDLs_df = pd.DataFrame(BMDLs, columns = column_names)
        self.BMD10s_df = pd.DataFrame(BMD10s, columns = column_names)
        self.BMD50s_df = pd.DataFrame(BMD50s, columns = column_names)
    
    #######################
    ## SELECT BEST MODEL ##
    #######################

    # Find the best model per dataset
    best_model = {item: "" for item in self.AIC_df["bmdrc.Endpoint.ID"].tolist()}

    # Keep a list of possible model selections
    poss_models = ["asymptotic", "exponential", "gompertz", "hill", "michaelis-mentin", "power", "weibull"]

    # Find the best model
    for row in range(len(self.AIC_df)):
        
        # Extract the endpoint
        endpoint = self.AIC_df["bmdrc.Endpoint.ID"][row]

        # Step 1: Extract the best models based on AIC
        values = np.array(self.AIC_df.iloc[row, 2:].to_list())
        min_value = np.min([val for val in values if not np.isnan(val)])
        model_choices = [poss_models[x] for x in range(len(values)) if not np.isnan(x) and (values[x] - min_value <= 2)] 
        if len(model_choices) == 1:
            best_model[endpoint] = model_choices[0]

        # Step 2: Select the smallest BMDL
        else:

            # Make a dictionary of BMDLs
            remaining = self.BMDLs_df[self.BMDLs_df["bmdrc.Endpoint.ID"] == endpoint][model_choices].reset_index(drop = True).loc[0].to_dict()

            # Select smallest BMDL if it is possible
            if not all(np.isnan(value) for value in remaining.values()):
                best_model[endpoint] = min(remaining, key = remaining.get)

            # Step 3: Otherwise, select smallest BMD10
            else:
                remaining = self.BMD10s_df[self.BMD10s_df["bmdrc.Endpoint.ID"] == endpoint][model_choices].reset_index(drop = True).loc[0].to_dict()
                best_model[endpoint] = min(remaining, key = remaining.get)

    #####################
    ## BUILD BMD TABLE ##
    #####################

    ## Collect fitting statistics for best model ##

    BMDS_model = []

    # Collect statistics for the best model
    for endpoint in best_model:
    
        # Identify the model
        model = best_model[endpoint]
    
        # Build a dictionary that holds the information per row
        rowDict = {
            "bmdrc.Endpoint.ID": endpoint,
            "Model": model,
            "BMD10": self.BMD10s_df[self.BMD10s_df["bmdrc.Endpoint.ID"] == endpoint][model].values[0],
            "BMDL": self.BMDLs_df[self.BMDLs_df["bmdrc.Endpoint.ID"] == endpoint][model].values[0],
            "BMD50": self.BMD50s_df[self.BMD50s_df["bmdrc.Endpoint.ID"] == endpoint][model].values[0]
        }
    
        # Append the list of dictionaries
        BMDS_model.append(rowDict)
    
    # Save half of the table
    bmds_stats = pd.DataFrame(BMDS_model)

    ## Collect other BMD metrics and merge ## 

    dose_response_groups = self.plate_groups[self.plate_groups["bmdrc.filter"] == "Keep"].groupby("bmdrc.Endpoint.ID")

    # Calculate Min_Dose, Max_Dose, AUC, and AUC_Norm
    bmds = dose_response_groups.apply(lambda df: np.trapezoid(df[self._response], x = df[self._concentration])).reset_index().rename(columns = {0: "AUC"})
    bmds["Min_Dose"] = dose_response_groups[["bmdrc.Endpoint.ID", self._concentration]].min(self._concentration).reset_index()[self._concentration]
    bmds["Max_Dose"] = dose_response_groups[["bmdrc.Endpoint.ID", self._concentration]].max(self._concentration).reset_index()[self._concentration]
    
    # Calculate the total area
    bmds["Max_Response"] = dose_response_groups[["bmdrc.Endpoint.ID", self._response]].max(self._response).reset_index()[self._response]
    bmds["Area"] = (bmds["Max_Dose"] - bmds["Min_Dose"]) * bmds["Max_Response"]
    
    # Normalize the AUC by the area
    bmds["AUC_Norm"] = round(bmds["AUC"] / bmds["Area"], 8)
    bmds["AUC"] = round(bmds["AUC"], 4)
    
    # Save resulting dataframe
    self.bmds = bmds_stats.merge(bmds)[["bmdrc.Endpoint.ID", "Model", "BMD10", "BMDL", "BMD50", "AUC", "Min_Dose", "Max_Dose", "AUC_Norm"]]

    # Set a flag that model fitting has been completed
    self.report_model_fits = True 





    