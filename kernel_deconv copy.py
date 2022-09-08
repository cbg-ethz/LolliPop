import pandas as pd
import numpy as np

# temporary, globals
tally_data = "data/.tsv"
variants_list = [
    "B.1.1.7",
    "B.1.351",
    "P.1",
    "B.1.617.2",
    "B.1.617.1",
    "BA.1",
    "BA.2",
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

# data preprocessing


class DataPreprocesser:
    """
    General class to preprocess tallymut data before deconvolution.
    """

    def __init__(self, df_tally):
        self.df_tally = df_tally

    def make_complement(self, df_tally, variants_list):
        """return a dataframe with the complement of mutations signatures and mutations fracs"""
        t_data = df_tally.copy()
        t_data["mutations"] = "-" + t_data["mutations"]
        t_data["frac"] = 1 - t_data["frac"]
        t_data[variants_list] = 1 - t_data[variants_list]
        t_data["undetermined"] = 1

        return t_data

    def general_preprocess(
        self,
        variants_list,
        variants_pangolin,
        variants_not_reported,
        to_drop,
        start_date=None,
        end_date=None,
        remove_deletions=True,
    ):
        """General preprocessing steps"""
        # rename columns
        self.df_tally = self.df_tally.rename(columns=variants_pangolin)
        # drop non reported variants
        self.df_tally = self.df_tally.drop(variants_not_reported, axis=1)
        # drop rows without estimated frac or date
        self.df_tally.dropna(subset=["frac", "date"], inplace=True)
        # create column with mutation signature
        self.df_tally["mutations"] = (
            self.df_tally["pos"].astype(str) + self.df_tally["base"]
        )
        # convert date string to date object
        self.df_tally["date"] = pd.to_datetime(self.df_tally["date"])
        # filter by minimum and maximum dates
        if start_date is not None:
            self.df_tally = self.df_tally[(self.df_tally["date"] >= start_date)]
        if end_date is not None:
            self.df_tally = self.df_tally[(self.df_tally["date"] < end_date)]
        # remove deletions
        if remove_deletions:
            self.df_tally = self.df_tally[~(self.df_tally["base"] == "-")]

        # df_data = df_data[df_data.columns.difference(['pos', 'gene', 'base'], sort=False)]
        for v in variants_list:
            self.df_tally = self.df_tally[~self.df_tally[v].isin(to_drop)]
        # drop index
        self.df_tally = self.df_tally.reset_index(drop=True)

        # this should be done very differently: create 0-1 matrix of definitions
        self.df_tally = self.df_tally.replace(np.nan, 0)
        self.df_tally = self.df_tally.replace(["extra", "mut", "shared"], 1)

        # make complement of mutation signatures for undetermined cases
        self.df_tally.insert(self.df_tally.columns.size - 1, "undetermined", 0)
        self.df_tally = pd.concat(
            [self.df_tally, self.make_complement(self.df_tally, variants_list)]
        )

        return self

    def filter_mutations(self):
        """very temporary function, to filter out hardcoded problematic mutations"""

        self.df_tally = self.df_tally[
            ~self.df_tally["mutations"].isin(["28461G", "11201G", "26801C"])
        ]

        self.df_tally = self.df_tally[
            ~(
                (pd.to_datetime(self.df_tally.date) > np.datetime64("2021-11-20"))
                & (self.df_tally.pos >= 22428)
                & (self.df_tally.pos <= 22785)
            )
        ]  # amplicon75
        self.df_tally = self.df_tally[
            ~(
                (pd.to_datetime(self.df_tally.date) > np.datetime64("2021-11-20"))
                & (self.df_tally.pos >= 22677)
                & (self.df_tally.pos <= 23028)
            )
        ]  # amplicon76
        self.df_tally = self.df_tally[
            ~(
                (pd.to_datetime(self.df_tally.date) > np.datetime64("2021-11-20"))
                & (self.df_tally.pos >= 22974)
                & (self.df_tally.pos <= 23327)
            )
        ]  # amplicon77
        self.df_tally = self.df_tally[
            ~(
                (pd.to_datetime(self.df_tally.date) > np.datetime64("2021-11-20"))
                & (self.df_tally.pos >= 26277)
                & (self.df_tally.pos <= 26635)
            )
        ]  # amplicon88
        self.df_tally = self.df_tally[
            ~(
                (pd.to_datetime(self.df_tally.date) > np.datetime64("2021-11-20"))
                & (self.df_tally.pos >= 26895)
                & (self.df_tally.pos <= 27256)
            )
        ]  # amplicon90
        self.df_tally = self.df_tally[
            ~(
                (pd.to_datetime(self.df_tally.date) > np.datetime64("2021-11-20"))
                & (self.df_tally.pos == 26709)
            )
        ]  # other
        self.df_tally = self.df_tally[
            ~(
                (pd.to_datetime(self.df_tally.date) > np.datetime64("2021-11-20"))
                & (self.df_tally.pos == 27807)
            )
        ]  # other
        self.df_tally = self.df_tally[
            ~(
                (pd.to_datetime(self.df_tally.date) > np.datetime64("2021-11-20"))
                & (self.df_tally.pos == 2832)
            )
        ]  # other
        self.df_tally = self.df_tally[
            ~(
                (pd.to_datetime(self.df_tally.date) > np.datetime64("2021-11-20"))
                & (self.df_tally.pos == 10449)
            )
        ]  # other

        return self


# Kernels


class GaussianKernel:
    """compute gaussian kernel weighting between y1 and y2"""

    def __init__(self, bandwidth=1.0):
        """
        bandwith (float): bandwith argument to kernel, default = 1.0
        """
        self.bandwidth = bandwidth

    def values(self, y1, y2):
        """
        compute gaussian kernel between y1 and y2
        """
        return np.exp(-((y1 - y2) ** 2) / 2 / self.bandwidth)


class BoxKernel:
    """compute box kernel between y1 and y2"""

    def __init__(self, bandwidth=1.0):
        """
        bandwith (float): bandwith argument to kernel, default = 1.0
        """
        self.bandwidth = bandwidth

    def values(self, y1, y2):
        """
        compute box kernel between y1 and y2
        """
        return 1.0 * (np.abs(y1 - y2) <= self.bandwidth / 2)


# Regressors


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


# Kernel Deconvolution


class KernelDeconv:
    """
    Compute kernel deconvolution, using specified kernel function and regressor.
    """

    def __init__(self, X, y, dates, kernel=GaussianKernel(), reg=NnlsReg()):
        """
        X (pd.DataFrame): dataframe of variant definition (design matrix)
        y (pd.Series): series of observed mutation frequencies
        dates (pd.Series): series of observation dates
        kernel (kernel object): object with methods to compute kernel weighting
        """
        self.X = X
        self.y = y
        self.dates = dates
        self.kernel = kernel
        self.reg = reg
        self.variant_names = X.columns

    def deconv(self, date, min_tol=1e-10):
        """
        compute kernel deconvolution centered on specific date, returns fitted regression object
        """
        # compute kernel values
        kvals = self.kernel.values(
            0, (date - self.dates) / pd.to_timedelta(1, unit="D")
        )
        # compute and return fitted coefs
        return self.reg.fit(
            self.X.values[kvals.values >= min_tol, :],
            self.y.values.flatten()[kvals.values >= min_tol],
            kvals.values[kvals.values >= min_tol],
        )

    def deconv_all(self, min_tol=1e-10):
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
        for date in self.dates.unique():
            deconv = self.deconv(date, min_tol)
            fitted.append(deconv.fitted)
            loss.append(deconv.loss)

        self.fitted = pd.DataFrame(
            np.array(fitted), columns=self.variant_names, index=self.dates.unique()
        )
        self.loss = np.array(loss)

        return self

    def renormalize(self):
        """renormalize variants proportion so that they sum to 1"""
        self.fitted = self.fitted.divide(self.fitted.sum(axis=1), axis=0)

        return self


def main():
    # load data

    df_tally = pd.read_csv(tally_data, sep="\t")
    df_proc = preprocess_data(
        df_tally,
        variants_list,
        variants_pangolin,
        variants_not_reported,
        start_date,
        remove_deletions=True,
    )


if __name__ == "__main__":
    main()
