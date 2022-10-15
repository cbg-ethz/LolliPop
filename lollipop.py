import pandas as pd
import numpy as np
import lollipop as ll
from scipy.optimize import nnls, least_squares
from tqdm import tqdm, trange

import click
import yaml
import os


@click.command(
    help="Deconvolution for Wastewater Genomics",
    #epilog="",
)
@click.option(
    "--out-dir",
    "--output-dir",
    "-o",
    metavar="DIR",
    required=False,
    default="./out",
    type=str,
    help="Output directory where to write results",
)
@click.option(
    "--config-file",
    "--config",
    "--conf",
    "-c",
    metavar="YAML",
    required=True,
    type=str,
    help="Variant configuration used during deconvolution",
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
def deconvolute(config_file, plant, seed, out_dir, tally_data):
    # load data
    print("load data")
    df_tally = pd.read_csv(tally_data, sep="\t")

    with open(config_file, 'r') as file:
        conf_yaml = yaml.load(file,  Loader=yaml.FullLoader)
    variants_list = conf_yaml['variants_list']
    variants_pangolin = conf_yaml['variants_pangolin']
    variants_not_reported = conf_yaml.get('variants_not_reported', [ ])
    to_drop=conf_yaml.get('to_drop', [])
    start_date=conf_yaml.get('start_date')
    end_date=conf_yaml.get('end_date')
    remove_deletions=conf_yaml.get('remove_deletions', True)
    cities_list = plant if plant and len(plant) else conf_yaml.get('cities_list', df_tally['plantname'].unique())

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
    for city in tqdm(cities_list):
        tqdm.write(city)
        temp_df = preproc.df_tally[preproc.df_tally["plantname"] == city]
        t_kdec = ll.KernelDeconv(
            temp_df[variants_list + ["undetermined"]],
            temp_df["frac"],
            temp_df["date"],
#         weights=temp_df["resample_value"],
            kernel=ll.GaussianKernel(10),
            reg=ll.NnlsReg(),
            confint=ll.NullConfint()
        )
        t_kdec = t_kdec.deconv_all()
        res = t_kdec.fitted
        res["city"] = city
        linear_deconv.append(res)
    linear_deconv_df = pd.concat(linear_deconv)
    linear_deconv_df = linear_deconv_df.fillna(0)

    print("output data")
    linear_deconv_df_flat = linear_deconv_df.melt(
        id_vars="city",
        value_vars=variants_list + ["undetermined"],
        var_name="variant",
        value_name="frac",
        ignore_index=False,
    )
    # linear_deconv_df_flat
    linear_deconv_df_flat.to_csv(os.path.join(out_dir, "deconvolved.csv"), index_label="date")


if __name__ == "__main__":
    deconvolute()
