import numpy as np
from scipy.stats import norm


class NullConfint:
    """Class to return Nones as confint"""

    def __init__(self):
        pass

    def confint(self, X, coefs):

        return {
            "lower": coefs * np.nan,
            "upper": coefs * np.nan,
        }


class WaldConfint:
    """class to compute wald se of the proportions deconvolved through a linear prob. model"""

    def __init__(self, level=0.95, scale="linear", pseudofrac=0.001):
        self.level = level
        self.scale = scale
        self.pseudofrac = pseudofrac

    def standard_error(self, X, coefs):
        """compute standard errors on the linear scale"""

        pseudocoefs = coefs / np.sum(coefs)
        # compute fitted values of mutation props.
        fitted_mut = X.dot(pseudocoefs)
        # add pseudofrac to avoid zero fitted mutations
        fitted_mut = fitted_mut + self.pseudofrac
        fitted_mut = fitted_mut / (fitted_mut + (1 - fitted_mut + 2 * self.pseudofrac))
        # compute Fisher information of in terms of binomial model
        binom_info = np.diag(1 / (fitted_mut * (1 - fitted_mut)))
        # compute Fisher information of in terms of linear prob. model
        linprob_info = X.T.dot(binom_info).dot(X)
        # try to invert Fisher information matrix
        try:
            se = np.sqrt(np.diag(np.linalg.inv(linprob_info)))
        except np.linalg.LinAlgError:
            se = coefs * np.nan

        return se

    def logit_standard_error(self, X, coefs):
        """compute the standard error on the logit scale using the Delta method"""

        pseudocoefs = coefs / np.sum(coefs)
        # compute fitted values of mutation props.
        fitted_mut = X.dot(pseudocoefs)
        # add pseudofrac to avoid zero fitted mutations
        fitted_mut = fitted_mut + self.pseudofrac
        fitted_mut = fitted_mut / (fitted_mut + (1 - fitted_mut + 2 * self.pseudofrac))
        # compute Fisher information of in terms of binomial model
        binom_info = np.diag(1 / (fitted_mut * (1 - fitted_mut)))
        # compute Fisher information of in terms of linear prob. model
        linprob_info = X.T.dot(binom_info).dot(X)
        # try to invert Fisher information matrix
        try:
            inv_info = np.linalg.inv(linprob_info)
        except np.linalg.LinAlgError:
            inv_info = linprob_info * np.nan

        # construct the jacobian
        pseudocoefs = (pseudocoefs + self.pseudofrac) / (
            np.sum(pseudocoefs) + 2 * self.pseudofrac
        )
        jacobian = -1 * np.tile(
            np.expand_dims(1 / (pseudocoefs * (1 - pseudocoefs)), 1), pseudocoefs.size
        )
        jacobian = jacobian + 2 * np.diag(-1 * np.diag(jacobian))

        # project inverse fisher inf. on logit scale
        logit_inv_info = jacobian.dot(inv_info).dot(jacobian.T)

        return np.sqrt(np.diag(logit_inv_info))

    def confint(self, X, coefs):
        """compute confidence intervals on the linear scale"""

        if self.scale == "linear":
            se = self.standard_error(X, coefs)
            return {
                "lower": coefs - norm.ppf(1 - (1 - self.level) / 2) * se,
                "upper": coefs + norm.ppf(1 - (1 - self.level) / 2) * se,
            }

        elif self.scale == "logit":
            se = self.logit_standard_error(X, coefs)
