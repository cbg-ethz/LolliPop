import pandas as pd
import numpy as np
import lollipop as ll
from scipy.optimize import nnls, least_squares
from tqdm import tqdm, trange

import click
import ruamel.yaml
import os


kernels = {
    "gaussian" : ll.GaussianKernel,
    "box" : ll.BoxKernel,
}
confints = {
    "null": ll.NullConfint,
    "wald": ll.WaldConfint,
}
regressors = {
    "nnls" : ll.NnlsReg,
    "robust" : ll.RobustReg,
}


@click.command(
    help="Deconvolution for Wastewater Genomics",
    #epilog="",
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
    "--deconv-config",
    "--dec",
    "-d",
    metavar="YAML",
    required=True,
    type=str,
    help="configuration of parameters for deconvolution",
)
@click.option(
    "--plant",
    "--catchment",
    "--wwtp",
    "-t",
    metavar="NAME",
    required=False,
    multiple=True,
    default=None,
    help="Name(s) of wastewater treatment plant/catchment area to process",
)
@click.option(
    "--seed",
    "-s",
    metavar="SEED",
    required=False,
    default=None,
    type=int,
    help="Seed the random generator"
)
@click.argument("tally_data", metavar="TALLY_TSV", nargs=1)
def deconvolute(variants_config, deconv_config, plant, seed, output, tally_data):
    # load data
    print("load data")
    with open(variants_config, 'r') as file:
        conf_yaml = ruamel.yaml.load(file, Loader=ruamel.yaml.Loader)
    variants_list = conf_yaml['variants_list']
    variants_pangolin = conf_yaml['variants_pangolin']
    variants_not_reported = conf_yaml.get('variants_not_reported', [ ])
    to_drop=conf_yaml.get('to_drop', [])
    start_date=conf_yaml.get('start_date')
    end_date=conf_yaml.get('end_date')
    remove_deletions=conf_yaml.get('remove_deletions', True)
    cities_list = plant if plant and len(plant) else conf_yaml.get('cities_list', None)

    with open(deconv_config, 'r') as file:
        deconv = ruamel.yaml.load(file, Loader=ruamel.yaml.Loader)

    df_tally = pd.read_csv(tally_data, sep="\t")
    if cities_list is None:
        cities_list = df_tally['plantname'].unique()

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
    linear_deconv = []
    # TODO parameters sanitation (e.g.: JSON schema, check in list)
    # kernel
    kernel = kernels.get(deconv.get("kernel"),ll.GaussianKernel)
    kernel_params = deconv.get("kernel_params", { })
    # confint
    confint = confints.get(deconv.get("confint"),ll.NullConfint)
    have_confint = confint != ll.NullConfint
    confint_name = deconv["confint"].capitalize() if have_confint else None
    confint_params = deconv.get("confint_params", { })
    # regressor
    regressor = regressors.get(deconv.get("regressor"),ll.NnlsReg)
    # deconv
    deconv_params = deconv.get("deconv_params", { })
    print(f""" parameters:
  kernel: {kernel}
   params: {kernel_params}
  confint: {confint}
   params: {confint_params}
   name: {confint_name}
   does actually generate confints: {have_confint}
  regressor: {regressor}""")
    for city in tqdm(cities_list) if len(cities_list) > 1 else cities_list:
        tqdm.write(city)
        temp_df = preproc.df_tally[preproc.df_tally["plantname"] == city]
        t_kdec = ll.KernelDeconv(
            temp_df[variants_list + ["undetermined"]],
            temp_df["frac"],
            temp_df["date"],
#         weights=temp_df["resample_value"],
            kernel=kernel(**kernel_params),
            reg=regressor(),
            confint=confint(**confint_params)
        )
        t_kdec = t_kdec.deconv_all(**deconv_params)
        if have_confint:
            # with conf int
            res = t_kdec.fitted.copy()
            res["city"] = city
            res["estimate"] = "MSE"
            linear_deconv.append(res)

            res_lower = t_kdec.conf_bands["lower"].copy()
            res_lower["city"] = city
            res_lower["estimate"] = f"{confint_name}_lower"
            linear_deconv.append(res_lower)

            res_upper = t_kdec.conf_bands["upper"].copy()
            res_upper["city"] = city
            res_upper["estimate"] = f"{confint_name}_upper"
            linear_deconv.append(res_upper)
        else:
            # without conf int
            res = t_kdec.fitted
            res["city"] = city
            linear_deconv.append(res)

    linear_deconv_df = pd.concat(linear_deconv)
    if not have_confint:
        linear_deconv_df = linear_deconv_df.fillna(0)

    print("output data")
    id_vars=["city"]
    if have_confint:
        id_vars += [ "estimate" ]

    linear_deconv_df_flat = linear_deconv_df.melt(
        id_vars=id_vars,
        value_vars=variants_list + ["undetermined"],
        var_name="variant",
        value_name="frac",
        ignore_index=False,
    )
    # linear_deconv_df_flat
    linear_deconv_df_flat.to_csv(output) #, index_label="date")


if __name__ == "__main__":
    deconvolute()
