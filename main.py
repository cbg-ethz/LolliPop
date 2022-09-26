import pandas as pd
import numpy as np
from scipy.optimize import nnls, least_squares
from .preprocessors import *
from .kernels import *
from .regressors import *
from .confints import *


# import matplotlib.pyplot as plt
# import seaborn as sns
from tqdm.notebook import tqdm, trange


# temporary, globals
tally_data = "./tallymut_line.tsv"
out_dir = (
    "./out"
)
variants_list = [
    "B.1.1.7",
    "B.1.351",
    "P.1",
    "B.1.617.2",
    "B.1.617.1",
    "BA.1",
    "BA.2",
    "BA.4",
    "BA.5"
]
variants_pangolin = {
    "al": "B.1.1.7",
    "be": "B.1.351",
    "ga": "P.1",
    "C36": "C.36.3",
    "ka": "B.1.617.1",
    "de": "B.1.617.2",
    "AY42": "AY.4.2",
    "B16173": "B.1.617.3",
    "om1": "BA.1",
    "om2": "BA.2",
    "BA1": "custom-BA.1",
    "BA2": "custom-BA.2",
}
variants_not_reported = [
    "custom-BA.1",
    "custom-BA.2",
    "C.36.3",
    "B.1.617.3",
    "AY.4.2",
    "mu",
    "d614g",
]
start_date = "2020-12-08"
to_drop = ["subset", "shared"]
cities_list = [
    "Lugano (TI)",
    "Zürich (ZH)",
    "Chur (GR)",
    "Altenrhein (SG)",
    "Laupen (BE)",
    "Genève (GE)",
    "Lausanne (VD)",
    "Basel (catchment area ARA Basel)",
    "Kanton Zürich",
]


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


def main():
    # load data
    print("load data")
    df_tally = pd.read_csv(tally_data, sep="\t")

    print("preprocess data")
    preproc = DataPreprocesser(df_tally)
    preproc = preproc.general_preprocess(
        variants_list=variants_list,
        variants_pangolin=variants_pangolin,
        variants_not_reported=variants_not_reported,
        to_drop=to_drop,
        start_date=start_date,
        remove_deletions=True,
    )
    preproc = preproc.filter_mutations()

    print("deconvolve all")
    linear_deconv = []
    for city in tqdm(cities_list):
        temp_df = preproc.df_tally[preproc.df_tally["plantname"] == city]
        t_kdec = KernelDeconv(
            temp_df[variants_list + ["undetermined"]],
            temp_df["frac"],
            temp_df["date"],
            kernel=GaussianKernel(10),
            reg=NnlsReg(),
        )
        t_kdec = t_kdec.deconv_all()
        res = t_kdec.renormalize().fitted
        res["city"] = city
        linear_deconv.append(res)
    linear_deconv_df = pd.concat(linear_deconv)

    print("output data")
    linear_deconv_df_flat = linear_deconv_df.melt(
        id_vars="city",
        value_vars=variants_list + ["undetermined"],
        var_name="variant",
        value_name="frac",
        ignore_index=False,
    )
    # linear_deconv_df_flat
    linear_deconv_df_flat.to_csv(out_dir + "/deconvolved.csv", index_label="date")


if __name__ == "__main__":
    main()
