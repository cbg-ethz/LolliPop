import pandas as pd
import numpy as np
from scipy.optimize import nnls, least_squares


class NnlsReg:
    """
    wrapper around nnls to feed to kernel_deconv
    """

    def __init__(self):
        pass

    def fit(self, X, y, k):
        """
        fit nnls to f(X) ~= y, weighted by values of k

        X (np.array): array of variant definition (design matrix)
        y (np.array): array of observed mutation frequencies
        k (np.array): kernel weighting values
        """
        self.fitted, self.loss = nnls(np.expand_dims(k, 1) * X, k * y)
        return self


class RobustReg:
    """
    wrapper around sp.optimize.least_squares to feed to kernel_deconv
    """

    def __init__(self, loss_type="soft_l1", f_scale=0.1):
        self.loss_type = loss_type
        self.f_scale = f_scale

    def fit(self, X, y, k, b0=None):
        """
        fit robust reg to f(X) ~= y, weighted by values of k

        X (np.array): array of variant definition (design matrix)
        y (np.array): array of observed mutation frequencies
        k (np.array): kernel weighting values
        b0 (np.array): starting values for optimization
        """
        # make starting values
        if b0 is None:
            b0 = np.ones(X.shape[1]) / X.shape[1]
        # regress
        ls = least_squares(
            lambda beta: (np.expand_dims(k, 1) * X).dot(beta) - (k * y),
            b0,
            bounds=(0, 1),
            loss=self.loss_type,
            f_scale=self.f_scale,
        )
        self.fitted, self.loss = ls.x, ls.cost
        return self
