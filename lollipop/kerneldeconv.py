import pandas as pd
import numpy as np
# from scipy.optimize import nnls, least_squares
from .kernels import GaussianKernel, BoxKernel
from .regressors import NnlsReg, RobustReg
from .confints import NullConfint, WaldConfint


class KernelDeconv:
    """
    Compute kernel deconvolution, using specified kernel function and regressor.
    """

    def __init__(
        self,
        X,
        y,
        dates,
        weights=None,
        kernel=GaussianKernel(),
        reg=NnlsReg(),
        confint=WaldConfint(),
    ):
        """
        X (pd.DataFrame): dataframe of variant definition (design matrix)
        y (pd.Series): series of observed mutation frequencies
        dates (pd.Series): series of observation dates
        weights (pd.Series): series of weights for the observations
        kernel (kernel object): object with methods to compute kernel weighting
        reg (regressor object): object with methods to compute the regression
        confint (confint object): object with method to compute confidence bands
        """
        self.X = X
        self.y = y
        self.dates = dates
        if weights is not None:
            self.weights = weights
        else:
            self.weights = np.ones_like(self.y)
        self.kernel = kernel
        self.reg = reg
        self.confint = confint
        self.variant_names = X.columns

    def deconv(self, date, min_tol=1e-10, renormalize=True):
        """
        compute kernel deconvolution centered on specific date, returns fitted regression object
        """
        # compute kernel values
        kvals = (
            self.kernel.values(0, (date - self.dates) / pd.to_timedelta(1, unit="D"))
            * self.weights
        )
        # compute and return fitted coefs
        regfit = self.reg.fit(
            self.X.values[kvals.values >= min_tol, :],
            self.y.values.flatten()[kvals.values >= min_tol],
            kvals.values[kvals.values >= min_tol],
        )

        # renormalize
        if renormalize:
            regfit.fitted = regfit.fitted / np.sum(regfit.fitted)

        # compute and return confint
        regfit.conf_band = self.confint.confint(
            X=self.X.values[kvals.values >= min_tol, :]
            * np.expand_dims(kvals.values[kvals.values >= min_tol], 1),
            coefs=regfit.fitted,
            y=self.y.values.flatten()[kvals.values >= min_tol],
            kvals=kvals.values[kvals.values >= min_tol]
        )

        return regfit

    def deconv_all(self, min_tol=1e-10, renormalize=True):
        """
        compute kernel deconvolution for all dates

        self.fitted (pd.DataFrame):
        """
        #         deconvolved = [self.deconv(date).__dict__ for date in self.dates.unique()]
        #         self.fitted = pd.DataFrame(
        #             np.array([dec["fitted"] for dec in deconvolved]),
        #             columns = self.variant_names
        #         )
        #         self.loss = np.array([dec["loss"] for dec in deconvolved])

        #         return self
        fitted = []
        loss = []
        lower = []
        upper = []
        for date in self.dates.unique():
            deconv = self.deconv(date, min_tol, renormalize)
            fitted.append(deconv.fitted)
            loss.append(deconv.loss)
            lower.append(deconv.conf_band["lower"])
            upper.append(deconv.conf_band["upper"])

        self.fitted = pd.DataFrame(
            np.array(fitted), columns=self.variant_names, index=self.dates.unique()
        )
        self.loss = np.array(loss)
        self.conf_bands = {
            "lower": pd.DataFrame(
                np.array(lower), columns=self.variant_names, index=self.dates.unique()
            ),
            "upper": pd.DataFrame(
                np.array(upper), columns=self.variant_names, index=self.dates.unique()
            ),
        }

        return self

    def renormalize(self):
        """renormalize variants proportion so that they sum to 1"""
        self.fitted = self.fitted.divide(self.fitted.sum(axis=1), axis=0)

        return self
