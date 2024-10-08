#!/usr/bin/env python3
import pandas as pd
import numpy as np
import lollipop as ll
from scipy.optimize import nnls, least_squares
from tqdm import tqdm, trange

import click
import ruamel.yaml
import json

import sys

import logging
import time

from typing import List, Tuple, Union 

import dask

from dask.distributed import Client, LocalCluster
from dask_jobqueue import SLURMCluster

import os

def get_dask_client():
    if 'SLURM_JOB_ID' in os.environ:
        # We're on a SLURM cluster
        from dask_jobqueue import SLURMCluster
        
        # Use default values that can be overridden by environment variables
        queue = os.environ.get('SLURM_QUEUE', 'default')
        project = os.environ.get('SLURM_PROJECT', 'default')
        cores = int(os.environ.get('SLURM_CORES', 1))
        memory = os.environ.get('SLURM_MEMORY', '4GB')
        jobs = int(os.environ.get('SLURM_JOBS', 1))
        
        cluster = SLURMCluster(
            queue=queue,
            project=project,
            cores=cores,
            memory=memory
        )
        cluster.scale(jobs=jobs)
        return Client(cluster)
    else:
        # We're running locally
        return Client(LocalCluster())

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


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


@dask.delayed
def _deconvolute_bootstrap(
        location: str, 
        preproc: ll.DataPreprocesser, 
        bootstrap : int, 
        locations_list: List, 
        no_loc: bool, 
        no_date: bool, 
        date_intervals: List[Tuple],
        var_dates: dict,
        kernel: Union[ll.GaussianKernel, ll.BoxKernel],
        kernel_params: dict,
        regressor: Union[ll.NnlsReg, ll.RobustReg],
        regressor_params: dict,
        confint: Union[ll.confints.NullConfint, ll.confints.WaldConfint],
        confint_params: dict,
        deconv_params: dict,
        have_confint: bool,
        confint_name: str
        ) -> List[pd.DataFrame]:
    """
    Deconvolute the data for a given location and bootstrap iteration.
    
    Parameters
    ----------
    location : str
        The location to deconvolute.
    preproc : ll.DataPreprocesser
        The preprocessed data object.
    bootstrap : int
        The number of bootstrap iterations to perform.
    locations_list : List
        The list of locations to deconvolute in total. Used for progress bar. 
        TODO: consider removing this parameter.
    no_loc : bool
        Whether to ignore location information.
    no_date : bool
        Whether to ignore date information.
    date_intervals : List[Tuple]
        The date intervals to deconvolute.
    var_dates : dict
        The variants to scan per periods.
    kernel : ll.Kernel
        The kernel to use for deconvolution.
    kernel_params : dict
        The parameters for the kernel.
    regressor : ll.Regressor
        The regressor to use for deconvolution.
    regressor_params : dict
        The parameters for the regressor.
    confint : ll.Confint
        The confidence interval to use for deconvolution.
    confint_params : dict
        The parameters for the confidence interval.
    deconv_params : dict
        The parameters for the deconvolution.
    have_confint : bool
        Whether to use a confidence interval.
    confint_name : str
        The name of the confidence interval.
    
    Returns
    -------
    List[pd.DataFrame]
        The deconvolution results for the location and bootstrap iterations.
    """

    # deconvolution results
    deconv = []

    if bootstrap <= 1 and len(date_intervals) <= 1:
        tqdm.write(location)

    # select the current location
    loc_df = (
        preproc.df_tally[preproc.df_tally["location"] == location]
        if not no_loc
        else preproc.df_tally
    )
    # print the memory usage of the current location
    logging.info(f"memory usage: {loc_df.memory_usage().sum() / 1024**2} MB")

    for b in (
        trange(bootstrap, desc=location, leave=(len(locations_list) > 1))
        if bootstrap > 1
        else [0]
    ):
        logging.info(f"bootstrap: {b}")
        start_time_b = time.time()
        if bootstrap > 1:
            # resample if we're doing bootstrapping
            temp_dfb = ll.resample_mutations(loc_df, loc_df.mutations.unique())[0]
        else:
            # just run one on everything
            temp_dfb = loc_df

        # print the memory usage of the current bootstrap
        logging.info(f"memory usage: {temp_dfb.memory_usage().sum() / 1024**2} MB")

        for mindate, maxdate in (
            tqdm(date_intervals, desc=location)
            if bootstrap <= 1 and len(date_intervals) > 1
            else date_intervals
        ):
            logging.info(f"date: {mindate} - {maxdate}")
            start_time_d = time.time()
            if not no_date:
                # filter by time period for period-specific variants list
                if maxdate is not None:
                    temp_df2 = temp_dfb[
                        temp_dfb.date.between(mindate, maxdate, inclusive="left")
                    ]
                else:
                    temp_df2 = temp_dfb[temp_dfb.date >= mindate]
            else:
                # no date => no filtering
                temp_df2 = temp_dfb
            if temp_df2.size == 0:
                continue
            # print the memory usage of the current date
            logging.info(f"memory usage: {temp_df2.memory_usage().sum() / 1024**2} MB")

            # remove uninformative mutations (present either always or never)
            variants_columns = list(
                set(var_dates["var_dates"][mindate]) & set(temp_df2.columns)
            )
            temp_df2 = temp_df2[
                ~temp_df2[variants_columns]
                .sum(axis=1)
                .isin([0, len(variants_columns)])
            ]
            if temp_df2.size == 0:
                continue

            # resampling weights
            if bootstrap > 1:
                weights = {"weights": temp_df2["resample_value"]}
            else:
                # just run one on everything
                weights = {}

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
                deconv.append(res)

                res_lower = t_kdec.conf_bands["lower"].copy()
                res_lower["location"] = location
                res_lower["estimate"] = f"{confint_name}_lower"
                deconv.append(res_lower)

                res_upper = t_kdec.conf_bands["upper"].copy()
                res_upper["location"] = location
                res_upper["estimate"] = f"{confint_name}_upper"
                deconv.append(res_upper)
            else:
                # without conf int
                res = t_kdec.fitted
                res["location"] = location
                deconv.append(res)
            logging.info(f"date took {time.time() - start_time_d} seconds")
        logging.info(f"bootstrap took {time.time() - start_time_b} seconds")

    return deconv


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
    "--fmt-columns",
    "-C",
    is_flag=True,
    default=False,
    help="Change output CSV format to one column per variant (normally, variants are each on a separate line)",
)
@click.option(
    "--out-json",
    "--oj",
    metavar="JSON",
    required=False,
    default=None,
    type=click.Path(),
    help="Also write a JSON results for upload to Cov-spectrum, etc.",
)
@click.option(
    "--variants-config",
    "--var",
    "-c",
    metavar="YAML",
    required=True,
    type=str,
    help="Variants configuration used during deconvolution",
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
    help="Configuration of parameters for kernel deconvolution",
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
    "--filters",
    "-fl",
    metavar="YAML",
    required=False,
    default=None,
    type=str,
    help="List of filters for removing problematic mutations from tally",
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
    variants_config,
    variants_dates,
    deconv_config,
    loc,
    filters,
    seed,
    output,
    fmt_columns,
    out_json,
    tally_data,
):
    # load data
    yaml = ruamel.yaml.YAML(typ="rt")
    print("load data")
    with open(variants_config, "r") as file:
        conf_yaml = yaml.load(file)
    variants_pangolin = conf_yaml["variants_pangolin"]
    variants_list = conf_yaml.get("variants_list", None)
    variants_not_reported = conf_yaml.get("variants_not_reported", [])
    to_drop = conf_yaml.get("to_drop", [])
    no_date = conf_yaml.get("no_date", False)
    no_loc = conf_yaml.get("no_loc", False)
    start_date = conf_yaml.get("start_date", None)
    end_date = conf_yaml.get("end_date", None)
    remove_deletions = conf_yaml.get("remove_deletions", True)
    locations_list = loc if loc and len(loc) else conf_yaml.get("locations_list", None)

    # kernel deconvolution params
    with open(deconv_config, "r") as file:
        deconv = yaml.load(file)

    # problematic mutation filters
    if filters:
        with open(filters, "r") as file:
            filters = yaml.load(file)

        print(f"{len(filters)} filter{ '' if len(filters) == 1 else 's' } loaded")
    else:
        filters = None

    # data
    try:
        df_tally = pd.read_csv(
            tally_data, sep="\t", parse_dates=["date"], dtype={"location_code": "str"}
        )
    except ValueError:
        df_tally = pd.read_csv(tally_data, sep="\t", dtype={"location_code": "str"})

    # handle location
    if not no_loc and "location" not in df_tally.columns:
        if "location_code" in df_tally.columns:
            print("NOTE: No location fullnames, using codes instead")
            df_tally["location"] = df_tally["location_code"]
        elif len(locations_list) == 1:
            print(
                "WARNING: No location in input data, assuming everything is {locations_list[0]}"
            )
            df_tally["location"] = locations_list[0]
        elif locations_list is None:
            print(
                f"WARNING: No location in input data. Either pass one with `--loc`/`locations_list` parameter  or set true the `no_loc` parameter {variants_config}"
            )
            no_loc = True
        else:
            print(
                f"ERROR: No location in input data. Either pass exactly one with `--loc`/`locations_list` parameter or set true the `no_loc` parameter {variants_config}"
            )
            sys.exit(1)

    if no_loc:
        if "location" in df_tally:
            locations_list = list(set(df_tally["location"].unique()) - {"", np.nan})
            if len(locations_list):
                print(
                    f"WARNING: no_loc is set, but there are still locations in input: {locations_list}"
                )
        else:
            print(
                "no_loc: ignoring location information and treating all input as a single location"
            )

        df_tally["location"] = "location"
        locations_list = ["location"]

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

    # check if dates are present
    if "date" not in df_tally.columns or all(df_tally["date"].isna()):
        if not no_date:
            no_date = True
            print(
                f"WARNING: No dates found in input data, automatically switching `no_date` !!!\n\tPlease either check input if this is not expected or add set true the `no_date` parameter to your {variants_config}"
            )
    # no date!!!
    if no_date:
        print("no_date: special mode for deconvoluting without time component")
        # HACK dummy date to keep the deconvolution kernel happy
        # add dummy date
        date_dict = dict(
            zip(
                df_tally["sample"].unique(),
                [
                    str(np.datetime64("1999-12-01") + np.timedelta64(i, "D"))
                    for i in range(len(df_tally["sample"].unique()))
                ],
            )
        )
        df_tally["date"] = pd.to_datetime(
            np.array([date_dict[i] for i in df_tally["sample"]])
        )

    # dates intervals for which to apply different variants as discovered using cojac
    if variants_dates:
        if no_date:
            print(
                f"WARNING: running in `no_date` mode, still var_dates specified in {variants_dates}"
            )
        with open(variants_dates, "r") as file:
            var_dates = yaml.load(file)

        all_var_dates = set(
            [var for lst in var_dates["var_dates"].values() for var in lst]
        )

        if variants_list is None:
            # build list of all variants from var_dates (if we did lack one)
            variants_list = list(all_var_dates)
        else:
            # have list => double - check it against var_dates
            not_on_date = list(set(variants_list) - all_var_dates)
            if len(not_on_date):
                print(
                    f"NOTE: {not_on_date} never used in {variants_dates}, despite being in variants_list"
                )
            not_on_list = list(all_var_dates - set(variants_list))
            if len(not_on_list):
                print(
                    f"WARNING: {variants_dates} lists variants: {not_on_list}, but they are not in variants_list"
                )
                variants_list += not_on_list
    else:
        if variants_list is None:
            # build list of all variants from lineage map (if we did lack one)
            variants_list = list(set(variants_pangolin.values()))

        if no_date:
            # dummy date
            var_dates = {"var_dates": {"1999-12-01": variants_list}}
        else:
            # search for all, always
            var_dates = {
                "var_dates": {conf_yaml.get("start_date", "2020-01-01"): variants_list}
            }
            print(
                "NOTE: deconvoluting for all variants on all dates. Consider writing a var_dates YAML based on cojac detections",
                file=sys.stderr,
            )

    # build the intervals pairs
    d = list(var_dates["var_dates"].keys())
    date_intervals = list(zip(d, d[1:] + [None]))
    if not no_date:
        for mindate, maxdate in date_intervals:
            if maxdate:
                assert (
                    mindate < maxdate
                ), f"out of order dates: {mindate} >= {maxdate}. Please fix the content of {variants_date}"
                print(f"from {mindate} to {maxdate}: {var_dates['var_dates'][mindate]}")
            else:
                print(f"from {mindate} onward: {var_dates['var_dates'][mindate]}")

    print("preprocess data")
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
    preproc = preproc.filter_mutations(filters=filters)

    print("deconvolve all")
    np.random.seed(seed)
    all_deconv = []
    # TODO parameters sanitation (e.g.: JSON schema, check in list)
    # bootstrap
    bootstrap = deconv.get("bootstrap", 0)
    # kernel
    kernel = kernels.get(deconv.get("kernel"), ll.GaussianKernel)
    kernel_params = deconv.get("kernel_params", {})
    if no_date:
        print("no_date: overriding kernel bandwidth")
        kernel_params["bandwidth"] = 1e-17
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

    # print out the iterations that will be done
    logging.info(f"locations: {locations_list}")
    logging.info(f"bootstrap: {bootstrap}")
    logging.info(f"date_intervals: {date_intervals}")
    
    # starting computation
    logging.info("starting computation")
    # start a timer
    start_time = time.time()
    # print the memory usage of the dataframe
    logging.info(f"memory usage: {df_tally.memory_usage().sum() / 1024**2} MB")

    # CORE DECONVOLUTION
    # iterate over locations
    # Create delayed objects for each location
    delayed_results = [_deconvolute_bootstrap(
        location=location,
        preproc=preproc,
        bootstrap=bootstrap,
        locations_list=locations_list,
        no_loc=no_loc,
        no_date=no_date,
        date_intervals=date_intervals,
        var_dates=var_dates,
        kernel=kernel,
        kernel_params=kernel_params,
        regressor=regressor,
        regressor_params=regressor_params,
        confint=confint,
        confint_params=confint_params,
        deconv_params=deconv_params,
        have_confint=have_confint,
        confint_name=confint_name
    ) for location in locations_list]

    # Compute all results in parallel
    # logging.info("Dask Multiprocessing Dashboard link: %s", client.dashboard_link)
    all_deconv = dask.compute(*delayed_results)

    # Flatten the results if necessary
    all_deconv = [item for sublist in all_deconv for item in sublist]

    logging.info(f"all locations took {time.time() - start_time} seconds")
    print("post-process data")
    deconv_df = pd.concat(all_deconv) # what does this do.
    if not have_confint:
        deconv_df = deconv_df.fillna(0)

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

    # deconv output
    deconv_df_flat = deconv_df.melt(
        id_vars=id_vars,
        value_vars=found_var + ["undetermined"],
        var_name="variant",
        value_name="frac",
        ignore_index=False,
    )
    # deconv_df_flat.to_csv(out_flat, sep="\t", index_label="date")

    # aggregation
    agg_columns = ["location", "variant", "index"]
    if bootstrap > 1:
        # bootstrap => mean + quantiles
        deconv_df_agg = (
            deconv_df_flat.reset_index()
            .groupby(agg_columns)
            .agg(
                [
                    "mean",
                    lambda x: np.quantile(x, q=0.025),
                    lambda x: np.quantile(x, q=0.975),
                ]
            )
            .reset_index()
        )

        export_columns = {
            ("index", ""): "date",
            ("frac", "mean"): "proportion",
            ("frac", "<lambda_0>"): "proportionLower",
            ("frac", "<lambda_1>"): "proportionUpper",
        }
    elif have_confint:
        # wald => pivot
        deconv_df_agg = (
            deconv_df_flat.reset_index()
            .pivot(index=agg_columns, columns="estimate")
            .reset_index()
        )

        export_columns = {
            ("index", ""): "date",
            ("frac", "MSE"): "proportion",
            ("frac", f"{confint_name}_lower"): "proportionLower",
            ("frac", f"{confint_name}_upper"): "proportionUpper",
        }
    else:
        # no conf => as-is
        deconv_df_agg = deconv_df_flat.reset_index()[agg_columns + ["frac"]]
        export_columns = {
            "index": "date",
            "frac": "proportion",
        }
    deconv_df_agg.columns = [
        export_columns.get(col, "".join(col) if type(col) is tuple else col)
        for col in deconv_df_agg.columns.values
    ]
    deconv_df_agg = deconv_df_agg.sort_values(by=["location", "variant", "date"])
    # reverse logit scale
    if have_confint and confint_params["scale"] == "logit":
        deconv_df_agg[["proportionLower", "proportionUpper"]] = deconv_df_agg[
            ["proportionLower", "proportionUpper"]
        ].applymap(
            lambda x: np.exp(np.clip(x, -100, 100))
            / (1 + np.exp(np.clip(x, -100, 100)))
        )

    ### CSV output
    if fmt_columns:
        output_df = (
            deconv_df_agg.reset_index()
            .pivot(
                index=["location", "date"],
                columns="variant",
                values=list(set(export_columns.values()) - {"date"}),
            )
            .reset_index()
        )
        output_df.columns = [
            (
                col
                if type(col) is not tuple
                else (
                    col[1]
                    if col[0] == "proportion"
                    else (
                        f"{col[1]}_{col[0][len('proportion'):]}"
                        if len(col[0]) > len("proportion")
                        else "".join(col)
                    )
                )
            )
            for col in output_df.columns.values
        ]
    else:
        output_df = deconv_df_agg
    print("output data")
    output_df.drop(
        (["location"] if no_loc else []) + (["date"] if no_date else []),
        axis=1,
        errors="ignore",
    ).to_csv(output, sep="\t", index=None)

    ### JSON
    print("output json")
    if out_json:
        update_data = {}

        loc_uniq = deconv_df_agg["location"].unique()
        var_uniq = deconv_df_agg["variant"].unique()

        json_columns = export_columns.values()
        if no_date:
            json_columns = list(set(json_columns) - {"date"})
        for loc in tqdm(loc_uniq, desc="Location", position=0):
            update_data[loc] = {}
            for var in tqdm(var_uniq, desc=loc, position=1, leave=False):
                tt_df = deconv_df_agg.loc[
                    (deconv_df_agg["variant"] == var)
                    & (deconv_df_agg["location"] == loc),
                    json_columns,
                ].copy()
                if not no_date:
                    tt_df["date"] = tt_df["date"].astype("str")

                update_data[loc][var] = {
                    "timeseriesSummary": [
                        dict(tt_df.iloc[i,]) for i in range(tt_df.shape[0])
                    ]
                }

        with open(out_json, "w") as file:
            file.write(
                json.dumps(update_data).replace("NaN", "null")
            )  # syntactically standard compliant JSON vs. python numpy's output.


if __name__ == "__main__":
    client = get_dask_client()
    deconvolute()
