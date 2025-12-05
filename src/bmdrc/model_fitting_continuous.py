import operator
from abc import abstractmethod

from scipy.stats import chi2
from sklearn.linear_model import LinearRegression
import numpy as np


############################
## CLASSES FOR MODEL FITS ##
############################

## General Functions ##
class Continuous_Model():
        
    ## Methods which do not change per model ##

    @abstractmethod
    def gof_p_value(self, y_obs: np.array[np.float64], y_pred: np.array[np.float64], params: np.array[np.float64]):
        '''
        Return a p-value of model fit (Goodness of Fit) for a fit dataset
        
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
        a single p-value from the goodness of fit test stored in the object
        
        '''
    
        # Calculate the chi squared value
        chi_squared = (((y_obs - y_pred)**2)).sum()
        
        # Degrees of freedom: Number of observations minus number of model parameters
        dof = len(y_obs) - len(params)  
        
        # Calculate p-value
        p_value = 1 - chi2.cdf(chi_squared, dof)
        
        self.p_value = p_value
    
    @abstractmethod
    def calculate_aic(self, y_obs: np.array[np.float64], y_pred: np.array[np.float64], params: np.array[np.float64]):
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

    ## Attributes used by all models ##

    def __init__(self, toModel, concentration, response):
        self.toModel = toModel
        self.concentration = concentration
        self.response = response

    toModel = property(operator.attrgetter('_toModel'))
    concentration = property(operator.attrgetter('_concentration'))
    response = property(operator.attrgetter('_response'))


## Linear Regression ##
class LinReg_Cont(Continuous_Model):

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
        return y / self.params[0]

    def fit(self, fixed_intercept = 0):
        '''
        Fit a linear regression model, calculate GOF and AIC

        Parameters
        ----------
        fixed_intercept
            a value that the the model must run through
        
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

        # Get the GOF value
        self.gof_p_value(y, self.y_pred, self.params)

        # Get the AIC value
        self.calculate_aic(y, self.y_pred, self.params)

