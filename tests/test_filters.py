import pandas as pd
import numpy as np
import lollipop as ll
import ruamel.yaml
from pandas.testing import assert_frame_equal

yaml = ruamel.yaml.YAML(typ="rt")


def assert_frame_NOT_equal(*args, **kwargs):
    try:
        assert_frame_equal(*args, **kwargs)
    except AssertionError:
        # frames are not equal
        pass
    else:
        # frames are equal
        raise AssertionError


def hardcoded_filter(df_tally):
    """this is the old hard-coded filter"""
    df_tally = df_tally[
        ~df_tally["mutations"].isin(
            ["28461G", "11201G", "26801C"] + ["-28461G", "-11201G", "-26801C"]
        )
    ]

    df_tally = df_tally[
        ~(
            (pd.to_datetime(df_tally.date) > np.datetime64("2021-11-20"))
            & (df_tally.pos >= 22428)
            & (df_tally.pos <= 22785)
        )
    ]  # amplicon75
    df_tally = df_tally[
        ~(
            (pd.to_datetime(df_tally.date) > np.datetime64("2021-11-20"))
            & (df_tally.pos >= 22677)
            & (df_tally.pos <= 23028)
        )
    ]  # amplicon76
    df_tally = df_tally[
        ~(
            (pd.to_datetime(df_tally.date) > np.datetime64("2021-11-20"))
            & (df_tally.pos >= 22974)
            & (df_tally.pos <= 23327)
        )
    ]  # amplicon77
    df_tally = df_tally[
        ~(
            (pd.to_datetime(df_tally.date) > np.datetime64("2021-11-20"))
            & (df_tally.pos >= 26277)
            & (df_tally.pos <= 26635)
        )
    ]  # amplicon88
    df_tally = df_tally[
        ~(
            (pd.to_datetime(df_tally.date) > np.datetime64("2021-11-20"))
            & (df_tally.pos >= 26895)
            & (df_tally.pos <= 27256)
        )
    ]  # amplicon90
    df_tally = df_tally[
        ~(
            (pd.to_datetime(df_tally.date) > np.datetime64("2021-11-20"))
            & (df_tally.pos == 26709)
        )
    ]  # other
    df_tally = df_tally[
        ~(
            (pd.to_datetime(df_tally.date) > np.datetime64("2021-11-20"))
            & (df_tally.pos == 27807)
        )
    ]  # other
    df_tally = df_tally[
        ~(
            (pd.to_datetime(df_tally.date) > np.datetime64("2021-11-20"))
            & (df_tally.pos == 2832)
        )
    ]  # other
    df_tally = df_tally[
        ~(
            (pd.to_datetime(df_tally.date) > np.datetime64("2021-11-20"))
            & (df_tally.pos == 10449)
        )
    ]  # other

    return df_tally


def test_filters():
    tally_file = "preprint/data/tallymut_line_full.tsv.zst"
    varconf_file = "config_preprint.yaml"
    filter_file = "filters_preprint.yaml"

    # filter to test
    with open(filter_file, "r") as f:
        filters = yaml.load(f)

    # load config
    with open(varconf_file, "r") as f:
        conf_yaml = yaml.load(f)
    variants_list = conf_yaml["variants_list"]
    variants_pangolin = conf_yaml["variants_pangolin"]
    variants_not_reported = conf_yaml.get("variants_not_reported", [])
    to_drop = conf_yaml.get("to_drop", [])
    start_date = conf_yaml.get("start_date", None)
    end_date = conf_yaml.get("end_date", None)
    no_date = conf_yaml.get("no_date", False)
    remove_deletions = conf_yaml.get("remove_deletions", True)

    # load data
    df_tally = pd.read_csv(
        tally_file, sep="\t", parse_dates=["date"], dtype={"location_code": "str"}
    )

    # pre-process data
    preproc = ll.DataPreprocesser(df_tally)
    preproc = preproc.general_preprocess(
        variants_list=variants_list,
        variants_pangolin=variants_pangolin,
        variants_not_reported=variants_not_reported,
        to_drop=to_drop,
        start_date=start_date,
        end_date=end_date,
        no_date=no_date,
        remove_deletions=remove_deletions,
    )

    # keep copy of original
    df_unfil = preproc.df_tally

    # old filtering
    df_hc = hardcoded_filter(df_unfil)

    assert_frame_NOT_equal(df_unfil, df_hc)

    # modern filtering
    preproc = preproc.filter_mutations(filters=filters)

    df_filt = preproc.df_tally

    assert_frame_NOT_equal(df_unfil, df_filt)

    # compare
    assert_frame_equal(df_filt[df_filt["proto"] == "v3"], df_hc[df_hc["proto"] == "v3"])
    assert_frame_equal(
        df_filt[df_filt["proto"] != "v3"], df_unfil[df_unfil["proto"] != "v3"]
    )
