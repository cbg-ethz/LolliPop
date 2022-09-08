import pandas as pd
import numpy as np


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
            ~self.df_tally["mutations"].isin(
                ["28461G", "11201G", "26801C"] + ["-28461G", "-11201G", "-26801C"]
            )
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
