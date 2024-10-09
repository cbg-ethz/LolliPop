import pandas as pd
import numpy as np
from functools import reduce
import re
import sys
from pandas.api.types import is_numeric_dtype


class DataPreprocesser:
    """
    General class to preprocess tallymut data before deconvolution.
    """

    def __init__(self, df_tally, namefield="mutations"):
        self.df_tally = df_tally

    def make_complement(self, df_tally, variants_list, namefield):
        """return a dataframe with the complement of mutations signatures and mutations fracs"""
        t_data = df_tally.copy()
        if namefield in t_data.columns:
            t_data[namefield] = "-" + t_data[namefield]
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
        no_date=False,
        remove_deletions=True,
        make_complement=True,
        namefield="mutations",
    ):
        """General preprocessing steps"""
        # rename columns
        assert len(variants_pangolin.values()) == len(
            set(variants_pangolin.values())
        ), f"duplicate values in:\n{variants_pangolin}"
        self.df_tally = self.df_tally.rename(columns=variants_pangolin)
        # drop non reported variants
        self.df_tally = self.df_tally.drop(
            variants_not_reported, axis=1, errors="ignore"
        )
        # drop rows without estimated frac or date
        self.df_tally.dropna(
            subset=["frac", "date"] if not no_date else ["frac"], inplace=True
        )
        # create column with mutation signature
        if ("base" in self.df_tally.columns) and ("pos" in self.df_tally.columns):
            # NOTE if cojac-based instead of SNV-bsed deconvolution: there is no single mutation cojac tabmut has an option to generate mutations strings
            self.df_tally["mutations"] = (
                self.df_tally["pos"].astype(str) + self.df_tally["base"]
            )
        # convert date string to date object
        # (also convert any dummy date of 'no_date'
        self.df_tally["date"] = pd.to_datetime(self.df_tally["date"])
        # filter by minimum and maximum dates
        if not no_date:
            if start_date is not None:
                self.df_tally = self.df_tally[(self.df_tally["date"] >= start_date)]
            if end_date is not None:
                self.df_tally = self.df_tally[(self.df_tally["date"] < end_date)]
        # remove deletions
        if remove_deletions:
            if "base" in self.df_tally.columns:
                self.df_tally = self.df_tally[~(self.df_tally["base"] == "-")]
            else:
                # NOTE if cojac-based instead of SNV-bsed deconvolution: there is no single mutation
                print(
                    f"Warning, 'remove_deletions' is set in configuration, but no 'base' column is present in columns {self.df_tally.columns}",
                    file=sys.stderr,
                )

        # df_data = df_data[df_data.columns.difference(['pos', 'gene', 'base'], sort=False)]

        # delete lines with mutation of the type that we want to delete (to_drop)
        # e.g.: remove all 'subset' mutations
        absentcol = set(variants_list) - set(self.df_tally.columns)
        if len(absentcol):
            # check for missing
            print(
                f"Warning, variants_list's {absentcol} is not present in columns {self.df_tally.columns}",
                file=sys.stderr,
            )
        for v in variants_list:
            if v in self.df_tally.columns:
                drop_mask = self.df_tally[v].isin(to_drop)
                if any(drop_mask):
                    self.df_tally = self.df_tally[~drop_mask]
        # drop index
        self.df_tally = self.df_tally.reset_index(drop=True)

        # this should be done very differently: create 0-1 matrix of definitions
        self.df_tally = self.df_tally.replace(np.nan, 0)
        # self.df_tally = self
        self.df_tally = self.df_tally.replace(
            ["extra", "mut", "shared", "revert", "subset"], 1
        ).infer_objects()

        # remove uninformative mutations
        variants_columns = list(set(variants_list) & set(self.df_tally.columns))
        self.df_tally = self.df_tally[
            ~self.df_tally[variants_columns]
            .sum(axis=1)
            .isin([0, len(variants_columns)])
        ]

        # make complement of mutation signatures for undetermined cases
        self.df_tally.insert(self.df_tally.columns.size - 1, "undetermined", 0)
        if make_complement:
            self.df_tally = pd.concat(
                [
                    self.df_tally,
                    self.make_complement(self.df_tally, variants_columns, namefield),
                ]
            )

        return self

    def filter_mutations(self, filters=None):
        """filter out hardcoded problematic mutations"""
        if not filters:
            return self

        types = self.df_tally.dtypes

        rxprser = re.compile(
            r"^ *(?P<qc>['\"]?)(?:(?P<col>"
            + r"|".join(self.df_tally.columns)
            + r")|(?P<bad>\w+))(?P=qc) *(?P<op>in|[<>=~!]*) *(?P<qv>['\"]?)(?P<val>.+)(?P=qv) *$"
        )

        def apply_filter_statement(name, fs):
            """parse a single statement from a filter and apply it, returning a boolean series"""
            m = rxprser.search(fs)
            assert m, f"Cannot parse statement <{fs}> in filter {name}"
            m = m.groupdict()

            assert m[
                "col"
            ], f"bad column name {m['bad']}, not in list: {self.df_tally.columns}, while parsing statement <{fs}> in filter {name}"

            # HACK handle 'date' column differently, to force datatypes
            col = (
                pd.to_datetime(self.df_tally["date"])
                if "date" == m["col"]
                else self.df_tally[m["col"]]
            )
            val = (
                np.datetime64(m["val"])
                if "date" == m["col"]
                else (
                    pd.to_numeric(m["val"])
                    if is_numeric_dtype(types[m["col"]])
                    else m["val"]
                )
            )

            # apply operator
            match m["op"]:
                case "=" | "==" | "" as e:
                    if e == "":
                        assert (
                            " " not in val
                        ), "Do not use values with space <{val}> when using no operator (implicit 'equals'). (while parsing statement <{fs}> in filter {name})"
                    return col == val
                case "!=" | "!":
                    return col != val
                case "<":
                    return col < val
                case "<=" | "=<":
                    return col <= val
                case ">=" | ">=":
                    return col >= val
                case ">":
                    return col > val
                case "in":
                    # unpack list
                    return col.isin(
                        [
                            v.strip("\"' ")
                            for v in val.lstrip("[ ").rstrip(" ]").split(",")
                        ]
                    )
                case "~" | "=~" | "~=":
                    return col.str.contains(
                        val[1:-1] if not m["qv"] and val[0] == val[-1] in "/@" else val
                    )
                case "!~" | "~!":
                    return ~(
                        col.str.contains(
                            val[1:-1]
                            if not m["qv"] and val[0] == val[-1] in "/@"
                            else val
                        )
                    )
                case _ as o:
                    raise ValueError(
                        f"unknown operator {o}, while parsing statement <{fs}> in filter {name}"
                    )

        for name, fl in filters.items():
            print(f"filter {name}")

            self.df_tally = self.df_tally[
                ~reduce(
                    (lambda x, y: x & y),
                    [apply_filter_statement(name, fstatmt) for fstatmt in fl],
                )
            ]

        # HACK completely disable filters
        return self
