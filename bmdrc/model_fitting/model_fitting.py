import numpy as np
import pandas as pd
import scipy.stats as stats
from astropy import stats as astrostats
from statsmodels.base.model import GenericLikelihoodModel

import warnings
warnings.filterwarnings('ignore')

__author__ = ["Paritosh Pande" , "David Degnan"]

############################
## CLASSES FOR MODEL FITS ##
############################

BMR = 0.1
BMR_50 = 0.5

## LOGISTIC CLASSES & FUNCTIONS ##

def logistic_fun(dose, params):
    alpha_ = params[0].astype('float')
    beta_ = params[1].astype('float')
    dose = dose.astype('float')
    prob_dose = 1/(1 + np.exp(-alpha_ - beta_*dose))
    return prob_dose

class Logistic(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Logistic, self).__init__(endog, exog, **kwds)


    def nloglikeobs(self, params):
        alpha_ = params[0]
        beta_ = params[1]
        dose = self.endog[:,0].flatten()

        params = [alpha_, beta_]
        probs = logistic_fun(dose, params)
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            mu_0 = self.endog[:,0].flatten().mean()
            s_0  = np.sqrt(3)*np.std(self.endog[:,0].flatten())/np.pi
            alpha_0 = -mu_0/s_0
            beta_0 = 1/s_0
            start_params = np.array([alpha_0, beta_0])

        return super(Logistic, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[None, None],[1e-5,None]], disp=0, **kwds)


class Logistic_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Logistic_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        alpha_ = params[0]
        bmdl_ = params[1]

        p_0 = 1/(1 + np.exp(-alpha_))

        chi_ = (1 - p_0) * BMR + p_0
        xi_ = np.log((1 - chi_)/chi_)
        beta_reparam = -(alpha_ + xi_)/bmdl_
        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = logistic_fun(dose, [alpha_, beta_reparam])
        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Logistic_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[None,None],[start_params[1],start_params[1]]],  disp = 0, **kwds)


## GAMMA CLASSES & FUNCTIONS ##

def gamma_fun(dose, params):
    g_ = params[0].astype('float')
    alpha_ = params[1].astype('float')
    beta_ = params[2].astype('float')
    dose = dose.astype('float')
    prob_dose = g_ + (1 - g_) * stats.gamma.cdf(dose, a = alpha_, scale = 1/beta_)
    return prob_dose

class Gamma(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Gamma, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0].astype('float')
        alpha_ = params[1].astype('float')
        beta_ = params[2].astype('float')
        dose = self.endog[:,0].flatten()
        params = [g_, alpha_, beta_]
        probs = gamma_fun(dose, params)

        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            g_0 = 0.1
            beta_0 = self.endog[:,0].flatten().mean()/self.endog[:,0].flatten().var()
            alpha_0 = self.endog[:,0].flatten().mean() * beta_0
            start_params = np.array([g_0, alpha_0, beta_0])

        return super(Gamma, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[0.2, 18],[1e-5, None]],  disp = 0, **kwds)

class Gamma_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Gamma_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        g_ = params[0].astype('float')
        alpha_ = params[1].astype('float')
        bmdl_ = params[2].astype('float')

        beta_reparam = stats.gamma.ppf(BMR, alpha_)/bmdl_

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = gamma_fun(dose, [g_, alpha_, beta_reparam])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Gamma_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[0.2, 18],[start_params[2],start_params[2]]], disp = 0, **kwds)


## WEIBULL CLASSES & FUNCTIONS ##

def weibull_fun(dose, params):
    g_ = params[0].astype('float')
    alpha_ = params[1].astype('float')
    beta_ = params[2].astype('float')
    dose = dose.astype('float')
    prob_dose = g_ + (1 - g_) * (1 - np.exp(-beta_ * (dose.astype('float') ** alpha_)))
    return prob_dose

class Weibull(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Weibull, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0].astype('float')
        alpha_ = params[1].astype('float')
        beta_ = params[2].astype('float')
        dose = self.endog[:,0].flatten()
        params = [g_, alpha_, beta_]
        probs = weibull_fun(dose, params)

        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            g_0 = 0.1

            dose = self.endog[:,0].flatten()
            num_affected = self.endog[:,1].flatten()
            num_total = self.endog[:,2].flatten()
            frac_affected = num_affected/num_total

            X = np.append(np.ones([len(dose[1:]),1]),np.log(np.reshape(dose[1:],(len(dose[1:]),1))),1)
            Y = np.array(np.reshape(np.log(-np.log(1 - frac_affected[1:])),(len(dose[1:]),1)))

            betas = np.linalg.inv((X.T).dot(X)).dot(X.T).dot(Y)

            alpha_0 = betas[1]
            beta_0 = np.exp(betas[0])

            start_params = np.array([g_0, alpha_0, beta_0])

        return super(Weibull, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[1e-5,None],[1e-9,None]], disp = 0, **kwds)

class Weibull_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Weibull_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0].astype('float')
        alpha_ = params[1].astype('float')
        bmdl_ = params[2].astype('float')

        beta_reparam = -np.log(1-BMR)/(bmdl_**alpha_)

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = weibull_fun(dose, [g_, alpha_, beta_reparam])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Weibull_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[1e-5,None],[start_params[2],start_params[2]]], disp = 0,**kwds)


## LOG-LOGISTIC CLASSES & FUNCTIONS ##

def log_logistic_fun(dose, params):
    g_ = params[0].astype('float')
    alpha_ = params[1].astype('float')
    beta_ = params[2].astype('float')
    dose = dose.astype('float')
    dose_nonzero = dose.copy()
    dose_nonzero[dose_nonzero == 0] = 1e-9
    prob_dose = g_ + (1 - g_)/(1 + np.exp(-alpha_ - beta_*np.log(dose_nonzero.astype('float'))))
    return prob_dose

class Log_Logistic(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Log_Logistic, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0].astype('float')
        alpha_ = params[1].astype('float')
        beta_ = params[2].astype('float')
        dose = self.endog[:,0].flatten()

        params = [g_,alpha_, beta_]
        probs = log_logistic_fun(dose, params)
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            doses = self.endog[:,0].copy().flatten()
            nonzero_doses = doses[1:]
            g_0 = 0.1
            mu_0 = np.log(nonzero_doses).mean()
            s_0  = np.sqrt(3)*np.std(np.log(nonzero_doses))/np.pi
            alpha_0 = -mu_0/s_0
            beta_0 = 1/s_0
            start_params = np.array([g_0, alpha_0, beta_0])

        return super(Log_Logistic, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[None, None],[None, None]],  disp = 0, **kwds)

class Log_Logistic_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Log_Logistic_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        g_ = params[0].astype('float')
        beta_ = params[1].astype('float')
        bmdl_ = params[2].astype('float')

        alpha_reparam = np.log(BMR/(1-BMR)) - beta_*np.log(bmdl_)

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = log_logistic_fun(dose, [g_, alpha_reparam, beta_])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Log_Logistic_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[start_params[1]/2,start_params[1]*2],[start_params[2],start_params[2]]], disp = 0,**kwds)


## PROBIT CLASSES & FUNCTIONS ##

def probit_fun(dose, params):
    alpha_ = params[0].astype('float')
    beta_ = params[1].astype('float')
    dose = dose.astype('float')
    prob_dose = stats.norm.cdf((alpha_ + beta_ * dose), loc=0, scale=1)
    return prob_dose

class Probit(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Probit, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        alpha_ = params[0].astype('float')
        beta_ = params[1].astype('float')
        dose = self.endog[:,0].flatten()

        probs = probit_fun(dose, [alpha_, beta_])

        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            dose = self.endog[:,0].flatten()
            num_affected = self.endog[:,1].flatten()
            num_total = self.endog[:,2].flatten()
            alpha_0 = stats.norm.ppf(num_affected[0]/num_total[0])
            beta_0 = (stats.norm.ppf(num_affected[-1]/num_total[-1]) - alpha_0)/dose[-1]
            start_params = np.array([alpha_0, beta_0])

        return super(Probit, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[None, None],[1e-5,None]],  disp = 0, **kwds)

class Probit_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Probit_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        alpha_ = params[0].astype('float')
        bmdl_ = params[1].astype('float')

        p_0 = stats.norm.cdf(alpha_)
        xi_ = (1 - p_0) * BMR + p_0

        beta_reparam = (stats.norm.ppf(xi_) - alpha_)/bmdl_

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = probit_fun(dose, [alpha_, beta_reparam])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Probit_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[None,None],[start_params[1],start_params[1]]], disp = 0,**kwds)


## LOG PROBIT CLASSES & FUNCTIONS ##

def log_probit_fun(dose, params):
    g_ = params[0].astype('float')
    alpha_ = params[1].astype('float')
    beta_ = params[2].astype('float')
    dose_nonzero = dose.copy().astype('float')
    dose_nonzero[dose_nonzero == 0] = 1e-9
    prob_dose = g_ + (1 - g_) * stats.norm.cdf((alpha_ + beta_ * np.log(dose_nonzero)), loc=0, scale=1)
    return prob_dose

class Log_Probit(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Log_Probit, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0]
        alpha_ = params[1]
        beta_ = params[2]
        dose = self.endog[:,0].flatten()

        probs = log_probit_fun(dose, [g_, alpha_, beta_])
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            g_0 = 0.1

            dose = self.endog[:,0].flatten()
            num_affected = self.endog[:,1].flatten()
            num_total = self.endog[:,2].flatten()
            frac_affected = num_affected/num_total
            X = np.array([[1, dose[1]],[1, dose[-1]]])
            Y = (np.array([stats.norm.ppf(1 - frac_affected[1]), \
                           stats.norm.ppf(1 - frac_affected[-1])])).T

            betas = np.linalg.inv((X.T).dot(X)).dot(X.T).dot(Y)

            alpha_0 = betas[0]
            beta_0 = max(1e-5,betas[1])

            dose = self.endog[:,0].flatten()
            start_params = np.array([g_0, alpha_0, beta_0])

        return super(Log_Probit, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[None,None],[1e-5,None]],  disp = 0, **kwds)

class Log_Probit_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Log_Probit_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        g_ = params[0].astype('float')
        alpha_ = params[1].astype('float')
        bmdl_ = params[2].astype('float')

        beta_reparam = (stats.norm.ppf(BMR) - alpha_)/np.log(bmdl_)

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = log_probit_fun(dose, [g_, alpha_, beta_reparam])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Log_Probit_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-9,0.99],[None,None],[start_params[2],start_params[2]]],  disp = 0, **kwds)


## MULTISTAGE 2 CLASSES & FUNCTIONS ##

def multistage_2_fun(dose, params):
    g_ = params[0]
    beta1_ = params[1]
    beta2_ = params[2]
    prob_dose = g_ + (1 - g_) * (1 - np.exp(-(beta1_ * dose) \
                                            -(beta2_ * (dose ** 2))))
    return prob_dose

class Multistage_2(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Multistage_2, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0].astype('float')
        beta1_ = params[1].astype('float')
        beta2_ = params[2].astype('float')
        dose = self.endog[:,0].flatten().astype('float')

        probs = multistage_2_fun(dose, [g_, beta1_, beta2_])
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            g_0 = 0.05

            dose = self.endog[:,0].flatten()
            num_affected = self.endog[:,1].flatten()
            num_total = self.endog[:,2].flatten()
            frac_affected = num_affected/num_total

            X = np.append(np.reshape(dose[1:],(len(dose[1:]),1)),(np.reshape(dose[1:],(len(dose[1:]),1)))**2,1)
            Y = np.array(-np.log(1 - frac_affected[1:]))

            betas = np.linalg.inv((X.T).dot(X)).dot(X.T).dot(Y)

            beta1_0 = betas[0]
            beta2_0 = betas[1]

            start_params = np.array([g_0, beta1_0, beta2_0])

        return super(Multistage_2, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-9,0.99],[1e-9,None],[1e-9,None]],  disp = 0, **kwds)

class Multistage_2_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Multistage_2_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        g_ = params[0].astype('float')
        beta1_ = params[1].astype('float')
        bmdl_ = params[2].astype('float')

        beta2_reparam = -(np.log(1-BMR) + beta1_*bmdl_)/(bmdl_**2)

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = multistage_2_fun(dose, [g_, beta1_, beta2_reparam])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Multistage_2_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-9,0.99],[1e-9,None],[start_params[2],start_params[2]]], disp = 0, **kwds)


## QUANTAL LINEAR CLASSES & FUNCTIONS ##

def quantal_linear_fun(dose, params):
    g_ = params[0].astype('float')
    beta_ = params[1].astype('float')
    dose = dose.astype('float')
    prob_dose = g_ + (1 - g_) * (1 - np.exp(-(beta_ * dose)))
    return prob_dose

class Quantal_Linear(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Quantal_Linear, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0].astype('float')
        beta_ = params[1].astype('float')

        dose = self.endog[:,0].flatten()

        probs = quantal_linear_fun(dose, [g_, beta_])
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            g_0 = 0.1
            beta_0 = 1/((self.endog[:,0].flatten().mean()))/np.log(2)
            start_params = np.array([g_0, beta_0])

        return super(Quantal_Linear, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[1e-5,None]],  disp = 0, **kwds)

class Quantal_Linear_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Quantal_Linear_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        g_ = params[0].astype('float')
        bmdl_ = params[1].astype('float')
        beta_reparam = -np.log(1-BMR)/bmdl_

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = quantal_linear_fun(dose, [g_, beta_reparam])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Quantal_Linear_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[start_params[1],start_params[1]]], disp = 0, **kwds)


#############################
## MODEL FITTING FUNCTIONS ##
#############################

def removed_endpoints_stats(self):
    '''
    Accessory function to fit the models. 
    As the first in the pipeline, this function calculates summary
    statistics for the endpoints that are filtered out. No models are
    fit in these values. 
    '''

    if any(self.plate_groups["bmdrc.filter"] == "Remove"):

        # Make a data frame with all filtered endpoints called low quality
        low_quality = self.plate_groups[self.plate_groups["bmdrc.filter"] == "Remove"]

        # Calculate the fraction affected 
        low_quality["frac.affected"] = low_quality["bmdrc.num.affected"] / low_quality["bmdrc.num.nonna"]

        # Group values by endpoint ID
        low_quality = low_quality.groupby("bmdrc.Endpoint.ID")

        # Calculate values 
        bmds_filtered = low_quality.apply(lambda df: np.trapz(df["frac.affected"], x = df["conc"])).reset_index().rename(columns = {0: "AUC"})
        bmds_filtered[["Model", "BMD10", "BMDL", "BMD50"]] = np.nan
        bmds_filtered["Min_Dose"] = round(low_quality[["bmdrc.Endpoint.ID", "conc"]].min("conc").reset_index()["conc"], 4)
        bmds_filtered["Max_Dose"] = round(low_quality[["bmdrc.Endpoint.ID", "conc"]].max("conc").reset_index()["conc"], 4)
        bmds_filtered["AUC_Norm"] = bmds_filtered["AUC"] / (bmds_filtered["Max_Dose"] - bmds_filtered["Min_Dose"])

        # Order columns
        self.bmds_filtered = bmds_filtered[["bmdrc.Endpoint.ID", "Model", "BMD10", "BMDL", "BMD50", "AUC", "Min_Dose", "Max_Dose", "AUC_Norm"]]

    else:

        self.bmds_filtered = None


def select_and_run_models(self):
    '''Iterate through each acceptable ID and fit models.'''

    # Pull dose_response
    dose_response = self.plate_groups[self.plate_groups["bmdrc.filter"] == "Keep"]

    # Calculate the fraction affected 
    dose_response["frac.affected"] = dose_response["bmdrc.num.affected"] / dose_response["bmdrc.num.nonna"]

    # Pull all values to fit
    to_fit = dose_response["bmdrc.Endpoint.ID"].unique().tolist()

    # Create a dictionary to hold all model results
    model_results = {}

    for endpoint in to_fit:

        # Subset to endpoint
        sub_data = dose_response[dose_response["bmdrc.Endpoint.ID"] == endpoint]

        # Calculate P-Value Function
        def calc_p_value(PredictedValues, Params):
            '''Return a p-value of model fit for each unique ID and Model dataframe pairing'''

            # Get the experimental values 
            ExperimentalValues = sub_data["frac.affected"].tolist()

            # Get count of non-na values
            NonNATotals = sub_data["bmdrc.num.nonna"].tolist() 

            # Now, calculate the chi squared value
            ChiSquared = ((NonNATotals / (PredictedValues * (1 - PredictedValues))) * (ExperimentalValues - PredictedValues)**2).sum()

            # Calculate a p-value of fit 
            return(1 - stats.chi2.sf(ChiSquared, len(NonNATotals) - len(Params)))

        # Regression model function
        def run_regression_model(sub_data, modelfun, fittedfun):
            '''Fit the regression model and return the parameters, fitted_values, and the p_value'''

            # Run the model
            model = modelfun(sub_data[[self.concentration, "bmdrc.num.affected", "bmdrc.num.nonna"]].astype('float').copy())

            # Get the model parameters
            model_params = model.fit().params

            # Get the model's fitted values
            model_fittedvals = fittedfun(sub_data[self.concentration], model_params)

            # Get the p_value
            model_pval = calc_p_value(model_fittedvals, model_params)

            # Get the AIC
            AIC = -2*model.fit().llf + (2 * len(model_params))

            # Return a list
            return([model, model_params, model_fittedvals, model_pval, AIC])

        # Run regression models in a dictionary
        models = {

            ## Logistic ##
            "Logistic": run_regression_model(sub_data, Logistic, logistic_fun),

            ## Gamma ## 
            "Gamma": run_regression_model(sub_data, Gamma, gamma_fun),

            ## Weibull ##
            "Weibull": run_regression_model(sub_data, Weibull, weibull_fun),

            ## Log-Logistic ##
            "Log Logistic": run_regression_model(sub_data, Log_Logistic, log_logistic_fun),

            ## Probit ##
            "Probit": run_regression_model(sub_data, Probit, probit_fun),

            ## Log-Probit ##
            "Log Probit": run_regression_model(sub_data, Log_Probit, log_probit_fun),

            ## Multistage ##
            "Multistage": run_regression_model(sub_data, Multistage_2, multistage_2_fun),

            ## Quantal Linear ##
            "Quantal Linear": run_regression_model(sub_data, Quantal_Linear, quantal_linear_fun),

        }

        # Iterate through all p-values 
        p_values = {}
        for key in models.keys():
            p_values[key] = models[key][3]

        # Iterate through all AICs 
        aics = {}
        for key in models.keys():
            aics[key] = models[key][4]

        # Determine the best model
        BestModel = min(aics, key=lambda k: aics[k]) 

        # Return results 
        model_results[endpoint] = [p_values, models[BestModel], BestModel, aics]

    self.model_fits = model_results


def fit_the_models(self, models, fit_threshold, BMD_Measurements):
    '''
    Fit the EPA recommended models to your dataset. 
    '''

    ##############################
    ## MAKE GROUPS IF NECESSARY ##
    ##############################

    try:
        self.plate_groups
    except AttributeError:
        print("No filters have been applied to this dataset, which is unusual. Proceeding with analysis.")
        self.make_plate_groups(self)

    ################
    ## FIT MODELS ##
    ################

    # Calculate statistics for endpoints that are filtered out
    removed_endpoints_stats(self)

    # Fit models and calculate statistics for endpoints that are not filtered out
    select_and_run_models(self)




