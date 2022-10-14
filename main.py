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
