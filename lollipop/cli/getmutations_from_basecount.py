#!/usr/bin/env python3
"""
Script retrieving the coverage of mutations from a single sample
"""
import pandas as pd
import numpy as np
import os
import click
from click_option_group import optgroup
import sys

__author__ = "Matteo Carrara"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


#####


def scan_basecnt(basecnt, tsvbase, mut):
    # warning that table is *tsvbase*-based
    basecount = (
        pd.read_csv(
            basecnt,
            sep="\t",
            header=[0, 1],
            index_col=[0, 1],
        )
        .droplevel("ref")
        .T.droplevel("sample")
        .T
    )
    # total coverage
    basecount["cov"] = basecount.apply(sum, axis=1)
    # look mutations per position
    return pd.DataFrame(
        data=mut.apply(
            lambda x: pd.concat(
                [
                    pd.Series(
                        [
                            x.gene,
                            x.position,
                            x.variant,
                            # -1 : 1-based to 0-based
                            basecount.loc[x.position - (1 - tsvbase)]["cov"],
                            basecount.loc[x.position - (1 - tsvbase)][x.variant],
                            (
                                basecount.loc[x.position - (1 - tsvbase)][x.variant]
                                / basecount.loc[x.position - (1 - tsvbase)]["cov"]
                                if basecount.loc[x.position - (1 - tsvbase)]["cov"]
                                else np.nan
                            ),
                        ],
                        index=[
                            "gene",
                            "pos",
                            "base",
                            "cov",
                            "var",
                            "frac",
                        ],
                    ),
                    pd.Series(x[4:]),
                ]
            ),
            axis=1,
        )
    )


###
"""
Helper functions
"""


def build_outname(outname, sample, batch):
    return f"{sample}_{batch}_mutations.txt" if outname is None else outname


###
@click.command(
    help="Search mutations and retrieve frequency from a TSV table produced by V-pipe",
)
@click.option(
    "--outname",
    "--output",
    "-o",
    required=False,
    type=click.Path(),
    help="Filename of the final output table. If not provided, it defaults to <samplename>_mutations.txt",
)
@click.option(
    "--muttable",
    "--mutationtable",
    "-m",
    required=False,
    default="mutlist.txt",
    type=click.Path(exists=True),
    help="Mutations helper table",
)
@click.option(
    "--based",
    "-a",
    "base",
    required=False,
    default=1,
    type=int,
    help="Are the positions in the tsv 0-based or 1-based?",
)
@click.argument(
    "basecnt",
    metavar="BASECOUNT",
    nargs=1,
    type=click.Path(exists=True),
)
@optgroup.group(
    "Argument used for simple concatenation",
    help="These options allows subsequently building simply by concatenation (using `xsv`, or even `tail` & `head`)",
)
@optgroup.option(
    "--location",
    "-l",
    required=False,
    type=str,
    default=None,
    help="Location of this sample",
)
@optgroup.option(
    "--date",
    "-d",
    required=False,
    type=str,
    default=None,
    help="Date of this sample",
)
@optgroup.group(
    "Argument use for V-pipe integration",
    help="These options help tracking output to the 2-level samples structure used by V-pipe",
)
@optgroup.option(
    "-s",
    "--sample",
    "--samplename",
    required=False,
    type=str,
    default=None,
    help="'sample_name' as found in the first column of the V-pipe samples.tsv",
)
@optgroup.option(
    "-b",
    "--batch",
    required=False,
    type=str,
    default=None,
    help="'batch'/'date' as in the second column of the V-pipe samples.tsv",
)
def from_basecount(outname, muttable, base, basecnt, location, date, sample, batch):
    outname = build_outname(outname, sample, batch)

    # list of mutations to search
    mut = pd.read_csv(muttable, sep="\t").astype({"position": "int"})

    # seach them!
    table = scan_basecnt(basecnt=basecnt, tsvbase=base, mut=mut)
    assert table.shape[0] > 0, "Generated an empty mutation table!"

    idx = []
    # add extra columns
    if sample:
        table["sample"] = sample
        idx += ["sample"]
    if batch:
        table["batch"] = batch
        idx += ["batch"]
    if location:
        table["location"] = location
        idx += ["location"]
    if date:
        table["date"] = date
        idx += ["date"]

    # set index
    idx += ["pos"]
    table.set_index(idx, inplace=True)

    print(outname)
    table.to_csv(outname, sep="\t", compression={"method": "infer"})


if __name__ == "__main__":
    from_basecount()
