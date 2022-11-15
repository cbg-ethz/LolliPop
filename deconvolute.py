#!/usr/bin/env python3
import pandas as pd
import numpy as np
import lollipop as ll
from scipy.optimize import nnls, least_squares
from tqdm import tqdm, trange

import click
import ruamel.yaml
import os
import sys


kernels = {
    "gaussian": ll.GaussianKernel,
    "box": ll.BoxKernel,
}
confints = {
    "null": ll.NullConfint,
    "wald": ll.WaldConfint,
}
regressors = {
    "nnls": ll.NnlsReg,
    "robust": ll.RobustReg,
}


@click.command(
    help="Deconvolution for Wastewater Genomics",
    # epilog="",
)
@click.option(
    "--output",
    "-o",
    metavar="CSV",
    required=False,
    default="deconvolved.csv",
    type=str,
    help="Write results to this output CSV instead of 'deconvolved.csv'",
)
@click.option(
    "--variants-config",
    "--var",
    "-c",
    metavar="YAML",
    required=True,
    type=str,
    help="Variant configuration used during deconvolution",
)
@click.option(
    "--variants-dates",
    "--vd",
    metavar="YAML",
    required=False,
    default=None,
    type=str,
    help="Variants to scan per periods (as determined with cojac)",
)
@click.option(
    "--deconv-config",
    "--dec",
    "-k",
    metavar="YAML",
    required=True,
    type=str,
    help="configuration of parameters for kernel deconvolution",
)
@click.option(
    "--loc",
    "--location",
    "--wwtp",
    "--catchment",
    "-l",
    metavar="NAME",
    required=False,
    multiple=True,
    default=None,
    help="Name(s) of location/wastewater treatment plant/catchment area to process",
)
@click.option(
    "--seed",
    "-s",
    metavar="SEED",
    required=False,
    default=None,
    type=int,
    help="Seed the random generator",
)
@click.argument("tally_data", metavar="TALLY_TSV", nargs=1)
def deconvolute(
    variants_config, variants_dates, deconv_config, loc, seed, output, tally_data
):
    # load data
    print("load data")
    with open(variants_config, "r") as file:
        conf_yaml = ruamel.yaml.load(file, Loader=ruamel.yaml.Loader)
    variants_list = conf_yaml["variants_list"]
    variants_pangolin = conf_yaml["variants_pangolin"]
    variants_not_reported = conf_yaml.get("variants_not_reported", [])
    to_drop = conf_yaml.get("to_drop", [])
    start_date = conf_yaml.get("start_date", None)
    end_date = conf_yaml.get("end_date", None)
    remove_deletions = conf_yaml.get("remove_deletions", True)
    locations_list = loc if loc and len(loc) else conf_yaml.get("locations_list", None)

    # TODO support date-less dataset
    # dates intervals for which to apply different variants as discovered using cojac
    if variants_dates:
        with open(variants_dates, "r") as file:
            var_dates = ruamel.yaml.load(file, Loader=ruamel.yaml.Loader)
    else:
        # search for all, always
        var_dates = {
            "var_dates": {conf_yaml.get("start_date", "2020-01-01"): variants_list}
        }
        print(
            "Warning: deconvoluting for all variants on all dates. Consider writing a var_dates YAML based on cojac detections",
            file=sys.stderr,
        )
    # build the intervals pairs
    d = list(var_dates["var_dates"].keys())
    date_intervals = list(zip(d, d[1:] + [None]))
    for mindate, maxdate in date_intervals:
        if maxdate:
            assert (
                mindate < maxdate
            ), f"out of order dates: {mindate} >= {maxdate}. Please fix the content of {variants_date}"
            print(f"from {mindate} to {maxdate}: {var_dates['var_dates'][mindate]}")
        else:
            print(f"from {mindate} onward: {var_dates['var_dates'][mindate]}")

    # kernel deconvolution params
    with open(deconv_config, "r") as file:
        deconv = ruamel.yaml.load(file, Loader=ruamel.yaml.Loader)

    # data
    df_tally = pd.read_csv(
        tally_data, sep="\t", parse_dates=["date"], dtype={"location_code": "str"}
    )
    if locations_list is None:
        # remember to remove empty cells: nan or empty cells
        locations_list = list(set(df_tally["location"].unique()) - {"", np.nan})
        print(locations_list)
    else:
        bad_locations = set(locations_list) - set(df_tally["location"].unique())
        assert 0 == len(
            bad_locations
        ), f"Bad locations in list: {bad_locations}, please fix {variants_config}."
        # locations_list = list(set(locations_list) - bad_locations)

    print("preprocess data")
    preproc = ll.DataPreprocesser(df_tally)
    preproc = preproc.general_preprocess(
        variants_list=variants_list,
        variants_pangolin=variants_pangolin,
        variants_not_reported=variants_not_reported,
        to_drop=to_drop,
        start_date=start_date,
        end_date=end_date,
        remove_deletions=remove_deletions,
    )
    preproc = preproc.filter_mutations()

    print("deconvolve all")
    np.random.seed(seed)
    all_deconv = []
    # TODO parameters sanitation (e.g.: JSON schema, check in list)
    # bootstrap
    bootstrap = deconv.get("bootstrap", 0)
    # kernel
    kernel = kernels.get(deconv.get("kernel"), ll.GaussianKernel)
    kernel_params = deconv.get("kernel_params", {})
    # confint
    confint = confints.get(deconv.get("confint"), ll.NullConfint)
    have_confint = confint != ll.NullConfint
    assert not (
        have_confint and bootstrap > 1
    ), f"either use bootstrapping or a confint class, not both at the same time.\nbootstrap: {bootstrap}, confint: {confint}"
    confint_name = deconv["confint"].capitalize() if have_confint else None
    confint_params = deconv.get("confint_params", {})
    # regressor
    regressor = regressors.get(deconv.get("regressor"), ll.NnlsReg)
    regressor_params = deconv.get("regressor_params", {})
    # deconv
    deconv_params = deconv.get("deconv_params", {})
    print(
        f""" parameters:
  bootstrap: {bootstrap}
  kernel: {kernel}
   params: {kernel_params}
  confint: {confint}
   params: {confint_params}
   name: {confint_name}
   non-dummy: {have_confint}
  regressor: {regressor}
   params: {regressor_params}
  deconv:
   params: {deconv_params}"""
    )

    # do it
    for location in tqdm(locations_list) if len(locations_list) > 1 else locations_list:
        if bootstrap <= 1 and len(date_intervals) <= 1:
            tqdm.write(location)
        # select the current location
        loc_df = preproc.df_tally[preproc.df_tally["location"] == location]
        for b in (
            trange(bootstrap, desc=location, leave=(len(locations_list) > 1))
            if bootstrap > 1
            else [0]
        ):
            if bootstrap > 1:
                # resample if we're doing bootstrapping
                temp_dfb = ll.resample_mutations(loc_df, loc_df.mutations.unique())[0]
                weights = {"weights": temp_dfb["resample_value"]}
            else:
                # just run one on everything
                temp_dfb = loc_df
                weights = {}

            for mindate, maxdate in (
                tqdm(date_intervals, desc=location)
                if bootstrap <= 1 and len(date_intervals) > 1
                else date_intervals
            ):
                if maxdate is not None:
                    temp_df2 = temp_dfb[
                        temp_dfb.date.between(mindate, maxdate, inclusive="left")
                    ]
                else:
                    temp_df2 = temp_dfb[temp_dfb.date >= mindate]
                if temp_df2.size == 0:
                    continue

                # deconvolution
                t_kdec = ll.KernelDeconv(
                    temp_df2[var_dates["var_dates"][mindate] + ["undetermined"]],
                    temp_df2["frac"],
                    temp_df2["date"],
                    kernel=kernel(**kernel_params),
                    reg=regressor(**regressor_params),
                    confint=confint(**confint_params),
                    **weights,
                )
                t_kdec = t_kdec.deconv_all(**deconv_params)
                if have_confint:
                    # with conf int
                    res = t_kdec.fitted.copy()
                    res["location"] = location
                    res["estimate"] = "MSE"
                    all_deconv.append(res)

                    res_lower = t_kdec.conf_bands["lower"].copy()
                    res_lower["location"] = location
                    res_lower["estimate"] = f"{confint_name}_lower"
                    all_deconv.append(res_lower)

                    res_upper = t_kdec.conf_bands["upper"].copy()
                    res_upper["location"] = location
                    res_upper["estimate"] = f"{confint_name}_upper"
                    all_deconv.append(res_upper)
                else:
                    # without conf int
                    res = t_kdec.fitted
                    res["location"] = location
                    all_deconv.append(res)

    deconv_df = pd.concat(all_deconv)
    if not have_confint:
        deconv_df = deconv_df.fillna(0)

    print("output data")
    id_vars = ["location"]
    if have_confint:
        id_vars += ["estimate"]

    # variants actually in dataframe
    found_var = list(set(variants_list) & set(deconv_df.columns))
    if len(found_var) < len(variants_list):
        print(
            f"some variants never found in dataset {set(variants_list) - set(found_var)}. Check the dates in {variants_dates}",
            file=sys.stderr,
        )

    deconv_df_flat = deconv_df.melt(
        id_vars=id_vars,
        value_vars=found_var + ["undetermined"],
        var_name="variant",
        value_name="frac",
        ignore_index=False,
    )
    # linear_deconv_df_flat
    deconv_df_flat.to_csv(output, sep="\t", index_label="date")


if __name__ == "__main__":
    deconvolute()
