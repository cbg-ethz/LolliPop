# LolliPop

[![Bioconda package](https://img.shields.io/conda/dn/bioconda/lollipop.svg?label=Bioconda)](https://bioconda.github.io/recipes/lollipop/README.html)
[![Docker container](https://quay.io/repository/biocontainers/lollipop/status)](https://quay.io/repository/biocontainers/lollipop)
[![bio.tools](https://img.shields.io/badge/bio-tools-blue.svg)](https://bio.tools/LolliPop)
[![Tests](https://github.com/cbg-ethz/LolliPop/actions/workflows/main.yaml/badge.svg)](https://github.com/cbg-ethz/LolliPop/actions/workflows/main.yaml)

LolliPop - a tool for Deconvolution for Wastewater Genomics

The LolliPop tool is part of the [V-pipe workflow for analysing NGS data of short viral genomes](https://github.com/cbg-ethz/v-pipe).

## Description

Wastewater-based monitoring has become an increasingly important source of
information on the spread of SARS-CoV-2 variants since clinical tests are
declining and may eventually disappear.

LolliPop has been developed to improve wastewater-based genomic surveillance
as the number of variants of concern increased and to account for shared
mutations among variants.
It relies on a kernel-based deconvolution, and leverages the time series
nature of the samples.
This approach enables to generate higher confidence relative abundance
curves despite the very high noise and overdispersion present in wastewater
samples.

It has been integrated in conjunction with
[COJAC](https://github.com/cbg-ethz/cojac) into
[V-pipe](https://github.com/cbg-ethz/v-pipe), a workflow designed for the
analysis of next generation sequencing (NGS) data from viral pathogens.
These tools now form the basis of the SARS-CoV-2 wastewater genomic
surveillance commissioned by the Swiss Federal Office of Public Health, a
cornerstone of the COVID-19 pandemic surveillance in Switzerland.
This [surveillance covers daily samples at ten wastewater treatment plants
](https://cov-spectrum.org/story/wastewater-in-switzerland) across Switzerland
from February 2021 onward, and delivers weekly updates of the variants relative
abundance curves.

## Usage

### Notebooks

LolliPop provides several classes that can be used imported in Jupyter notebooks
```python
from lollipop import *
```

See [notebook WwSmoothingKernel.ipynb](preprint/WwSmoothingKernel.ipynb)
in directory [preprint/](preprint/)

### Command line

Here are the available command-line tools:

| command                      | purpose |
| :--------------------------- | :------ |
| `lollipop  generate-mutlist` | Generate the mutlist used when looking for variant using variant signatures |
| `lollipop  getmutations from-basecount` | Search a single sample for mutations and retrieve frequency from a TSV table of per-position base counts produced by V-pipe |
| `lollipop deconvolute`       | Run the deconvolution on a timeline of mutations |

Use option `-h` / `--help` to see available command-line options:

```console
$ lollipop  generate-mutlist --help
Usage: lollipop generate-mutlist [OPTIONS] VOC_YAML

  Generate the mutlist used when looking for variant using variant signatures

Options:
  -o, --output, --out TSV         Write results to this output TSV instead of
                                  'mutlist.tsv'
  -p, --out-pangovars, --output-variants-pangolin YAML
                                  Write a YAML mapping shortnames/columnnames
                                  to the Pangolineages (useful to make the
                                  'variants_pangolin' section of deconvolute's
                                  input configuration)
  -g, --genes GFF                 Add 'gene' column to table
  -d, --voc-dir PATH              Scan directory for additional voc YAML files
  -v, --verbose / -V, --no-verbose
                                  Verbose (dumps table on terminal)
  -h, --help                      Show this message and exit.
```

```console
$ lollipop getmutations from-basecount --help  
Usage: lollipop getmutations from-basecount [OPTIONS] BASECOUNT  

  Search mutations and retrieve frequency from a TSV table produced by V-pipe

Options:
  -o, --outname, --output PATH    Filename of the final output table. If not
                                  provided, it defaults to
                                  <samplename>_mutations.txt
  -m, --muttable, --mutationtable PATH
                                  Mutations helper table
  -a, --based INTEGER             Are the positions in the tsv 0-based or
                                  1-based?
  Argument used for simple concatenation: 
                                  These options allows subsequently building
                                  simply by concatenation (using `xsv`, or
                                  even `tail` & `head`)
    -l, --location TEXT           Location of this sample
    -d, --date TEXT               Date of this sample
  Argument use for V-pipe integration: 
                                  These options help tracking output to the
                                  2-level samples structure used by V-pipe
    -s, --sample, --samplename TEXT
                                  'sample_name' as found in the first column
                                  of the V-pipe samples.tsv
    -b, --batch TEXT              'batch'/'date' as in the second column of
                                  the V-pipe samples.tsv
  -h, --help                      Show this message and exit.
```

```console
$ lollipop deconvolute --help
Usage: lollipop deconvolute [OPTIONS] TALLY_TSV

  Deconvolution for Wastewater Genomics

Options:
  -o, --output CSV                Write results to this output CSV instead of
                                  'deconvolved.csv'
  -C, --fmt-columns               Change output CSV format to one column per
                                  variant (normally, variants are each on a
                                  separate line)
  --out-json, --oj JSON           Also write a JSON results for upload to Cov-
                                  spectrum, etc.
  -c, --variants-config, --var YAML
                                  Variants configuration used during
                                  deconvolution  [required]
  --variants-dates, --vd YAML     Variants to scan per periods (as determined
                                  with cojac)
  -k, --deconv-config, --dec YAML
                                  Configuration of parameters for kernel
                                  deconvolution  [required]
  -l, --loc, --location, --wwtp, --catchment NAME
                                  Name(s) of location/wastewater treatment
                                  plant/catchment area to process
  -fl, --filters YAML             List of filters for removing problematic
                                  mutations from tally
  -s, --seed SEED                 Seed the random generator
  -n, --n-cores  N                Cores for parallel processing for multiple locations,
                                  defaults to 1 for sequential processing
  -nf, --namefield COLUMN         column to use as 'names' for the entries in
                                  tally table. By default, if 'pos' and 'base'
                                  exist a column 'mutations' will be created
                                  and used as name.
  -h, --help                      Show this message and exit.
```

## Howto

### Input data requirements

Analysis can be performed on virus samples sequenced with most tiled
multiplexed PCRs amplification protocols. Having coverage across the whole
genome of the virus increases the chance of some variant-specific mutations
being picked up and increasing the confidence, even if dropouts are experienced
on some other regions of the genome (e.g.: dropouts on the fragment carrying
the binding domain).

Sampling dates are important information to keep track of because LolliPop
leverages time series.

### Mutations lists

Analysis will use variants description YAML that lists mutations to be searched -- 
the same YAMLs as used by [COJAC](https://github.com/cbg-ethz/cojac). You can
refer to COJAC's commands `cojac sig-generate` to help generate exhaustive
lists from requests on [Cov-Spectrum](https://cov-spectrum.org/) or 
[TSV files of Covariants.org](https://github.com/hodcroftlab/covariants/blob/master/defining_mutations/),
or `cojac phe2cojac` to import ready-made manually-curated lists from YMLs available at
[PHE Genomic's _Standardised Variant Definitions_](https://github.com/phe-genomics/variant_definitions).

Generate a list of mutation to be searched:
```bash
lollipop generate-mutlist --output mutlist.tsv --out-pangovars variants_pangolin.yaml --genes Genes_NC_045512.2.GFF3 -- vocs/delta_mutations_full.yaml vocs/omicron_ba1_mutations_full.yaml vocs/omicron_ba2_mutations_full.yaml
```
- Annotating the list with a GFF file is optional: Lollipop's deconvolution
  does not use genes information, but it could be useful for downstream
  visualizations.
- `--out-pangovars` writes a table mapping back short names to full
  Pangolineages. It can be useful to help write (or be used in lieu of) a
  variants' config.

### Search mutations in a single sample

#### basecount table

By default, LolliPop searches the mutations into a basecount TSV, a table that
gives per position coverage of each A, T, C, G bases and deletion. 
[V-pipe](https://github.com/cbg-ethz/v-pipe) generates such a TSV using 
[smallgenomeutilities's command `aln2basecnt`](https://github.com/cbg-ethz/smallgenomeutilities/#aln2basecnt),
you can use it in your workflow when starting from alignments:

```bash
aln2basecnt --first 1 --basecnt sample1.basecnt.tsv.gz --coverage sample1.coverage.tsv.gz --name "sample1" sample1.bam
```
- `--first` is used to specify if the positions in the TSV are 1-based
  (like samtools) or 0-based (like pysam).

Then, search this TSV files for the mutations from the list generated above:
```bash
lollipop getmutations from-basecount --based 1 --output sample1.mut.tsv --location "main plant" --date "2023-02-27" -m mutlist.tsv -- sample1.basecnt.tsv.gz
```
- options `--location` and `--date` are a straightforward way to add the
  time series information for each sample

#### VCF and coverage

> (a future version of LolliPop will be extended to support VCFs and
> coverage TSV as a more standard input)

### Combine the time series

Once the above step has been run on every single sample of the cohort, combine
all individual samples into a single heatmap-like object tracking the mutation
overtime across all samples. This can be done by concatenating all the per-sample
mutations TSVs with a tool such as [xsv](https://github.com/BurntSushi/xsv):

```bash
xsv cat rows --output tallymut.tsv sample*.mut.tsv
```
> - If you have not tagged each individual sample with `--location` and `--date`, now it would be a good time to add extra columns to _tallymut.tsv_, e.g., with a _join_ operation.
> - Note that this file can get quite huge. It is possible to compress it on the fly: `… | xsv fmt --out-delimiter '\t' | gzip -o tallymut.tsv.gz`

### Run the deconvolution

The deconvolution can now be run on this table

#### Kernel deconvolution config

Various aspects of the kernel-based deconvolution can be set with a YAML file:
type of kernel (box vs Gaussian) and its parameters (such as bandwidth),
regressor used, using bootstrapping to generate confidence value, estimating
confidence intervals with Wald, computing the estimates on a logit scale, etc.

Various presets are available in the [presets/](presets/) subdirectory.

For example:
```yaml
kernel: 'gaussian'
kernel_params:
  bandwidth: 10

regressor: 'robust'

deconv_params:
  min_tol: 1e-3
```

#### Variants configuration

This file controls the data set that the deconvolution runs on. At minimum, it
should have a section mapping the short names back to full Pangolineages. This
can be copied by the file generated with `--out-pangovars` on the first step
(or that file reused as-is).

But this can also be used to optionally specify time limits (`start_date`
and/or ` end_date`), the subset of variants (`variants_list`) or locations
(`locations_list`) to run deconvolution onto, variants column to delete
(`variants_not_reported`) before processing any further, not considering the
deletions (`remove_deletions`), etc. 
see [example in config_preprint.yaml](config_preprint.yaml).

#### Variants dates

The deconvolution performs much better if only the variants known to be present
in the mixture are considered. For longer-running experiment, it is therefore
possible to specify, for different time periods, the list of variants to
consider for deconvolution, based on their previous detection with a sensitive
tool, e.g, such as determined running COJAC and looking for amplicons carrying
mutations combinations which are exclusive for certain variants.

For example:
```yaml
var_dates:
  '2022-06-15':
  - BA.1
  - BA.2
  - BA.4
  - BA.5
  - BA.2.75
  '2022-08-15':
  - BA.4
  - BA.5
  - BA.2.75
  - BQ.1.1
  '2022-11-01':
  - BA.4
  - BA.5
  - BA.2.75
  - BQ.1.1
  - XBB
```
see [variants_dates_example.yaml](variants_dates_example.yaml).

#### Filters (optional)

Some mutations might be problematic and need to be taken out -- e.g. 
due to drop-outs in the multiplex PCR amplification, they do not show up in the data
and this could be misinterpreted by LolliPop as proof of absence of a variant.
This optional file contains a collection of filters.
Each filter has a list of statements with the following syntax:
```text
- <column> <op> <value>
```
Valid _op_ are:
- `==` on that line, the value in _column_ is exactly _value_
  - for simple strings this can be omitted: `- proto v3` is synonymous with `- proto == v3`
- `<=` the value is less than or equal to _value_
- `>=` the value is greater than or equal to _value_
- `<` the value is less than _value_
- `>` the value is greater than _value_
- `!=` the value is **not** _value_
- `in` the value is found in the list specidied in _value_
- `~` the value matches the regular expression in _value_
  - regex can be quoted using `/` or `@`
- `!~` the vlue does **not** matche the regular expression in _value_

Any arbitrary column found in the input file can be used.

All statements of a filter are combined with a logical `and` and matching lines are removed from the tally table.

Filters are processed in the order found in the YAML file.

For example:
```yaml
# filter to remove test samples
remove_test:
- sample ~ /^Test/

# filter to remove an amplicon that has drop-outs
amplicon75:
  - proto v3
  - date > 2021-11-20
  - pos >= 22428
  - pos <= 22785
```
see [example in filters_preprint.yaml](filters_preprint.yaml).

#### Running it

```bash
lollipop deconvolution --output=deconvoluted.tsv --out-json=deconvoluted_upload.json --var=variants_conf.yaml --vd=variants_dates.yaml --dec=deconv_linear.yaml --seed=42  --n-cores=8 -- tallymut.tsv
```

### Output

The output is tabular:

| location   |       date | variant | proportion |
| :--------- | :--------- | :------ | ---------: |
| main plant | 2023-02-27 | BA.4    |      0.000 |

Optionally, LolliPop can also package the results in a JSON structure, e.g.,
to be sent to online dashboards:

```json
{
  "mainplant": {
    "BA.4": {
      "timeseriesSummary": [
        {
          "date": "2023-02-27",
          "proportion": 0.000
        },
        {
          "date":  "
          … etc …
          "
        }
      ]
    }
  }
}
```

The repository [cowwid](https://github.com/cbg-ethz/cowwid) contains real-world examples
of downstream analysis of the output of LolliPop.

## Installation

We recommend using [bioconda software repositories](https://bioconda.github.io/index.html)
for easy installation.
You can find instructions to setup your bioconda environment at the following
address:

 - https://bioconda.github.io/index.html#usage

### Prebuilt package

_LolliPop_ and its dependencies are all available in the bioconda repository.
We strongly advise you to install
[this pre-built package](https://bioconda.github.io/recipes/lollipop/README.html)
for a hassle-free experience.

You can install _lollipop_ in its own environment and activate it:

```bash
conda create -n lollipop lollipop
conda activate lollipop
# test it
lollipop --help
```

And to update it to the latest version, run:

```bash
# activate the environment if not already active:
conda activate lollipop
conda update lollipop
```

Or you can add it to the current environment (e.g.: in environment _base_):

```bash
conda install lollipop
```

### Building and deploying yourself

#### within conda environment

If you want to install the software yourself, you can see the list of
dependencies in [`conda_lollipop_env.yaml`](conda_lollipop_env.yaml).

We recommend using conda to install them:

```bash
conda env create -f conda_lollipop_env.yaml
conda activate lollipop
```

Install _lollipop_ using pip:
```bash
# install both the python module and the cli
pip install '.[cli]'
# (this will autodetect dependencies already installed by conda)
```

The command `lollipop` should now be accessible from your PATH

```bash
# activate the environment if not already active:
conda activate lollipop
lollipop --help
```

### Remove conda environment

You can remove the conda environment if you do not need it any more:

```bash
# exit the lollipop environment first:
conda deactivate
conda env remove -n lollipop
```

### Python poetry

LolliPop has its dependencies in a [pyproject.toml](pyproject.toml) managed
with poetry and can be installed with it.

```bash
# If not installed system-wide: manually run poetry-dynamic-versioning
poetry-dynamic-versioning
# (this sets the version string from the git currently cloned and checked out)

poetry install --extras "cli"
```

## Upcoming features

- [ ] Support VCFs and coverage TSV as alternative to basecount TSV

Long term goal:

~~- [x] Inputs other than SNVs: can deconvolute  COJAC's output tables~~

## Contributions

#### Package developers:

- [David Dreifuss ![orcid]](https://orcid.org/0000-0002-5827-5387), [![github]](https://github.com/dr-david)
- [Ivan Topolsky ![orcid]](https://orcid.org/0000-0002-7561-0810), [![github]](https://github.com/dryak)
- [Gordon Koehn ![orcid]](https://orcid.org/0000-0003-3397-7769), [![github]](https://github.com/gordonkoehn)
- [Pelin Icer Baykal ![orcid]](https://orcid.org/0000-0002-9542-5292), [![github]](https://github.com/picerbaykal)
- [Mateo Carrara ![orcid]](	https://orcid.org/0000-0002-8559-8296), [![github]](https://github.com/mcarrara-bioinfo)

#### Corresponding author:

 - [Niko Beerenwinkel ![orcid]](https://orcid.org/0000-0002-0573-6119)

[github]: images/mark-github.svg
[orcid]: images/ORCIDiD_iconvector.svg


## Citation

If you use this software in your research, please cite:

- David Dreifuss, Ivan Topolsky, Pelin Icer Baykal & Niko Beerenwinkel

  "*Tracking SARS-CoV-2 genomic variants in wastewater sequencing data with* LolliPop."

  medRxiv; [doi:10.1101/2022.11.02.22281825](https://doi.org/10.1101/2022.11.02.22281825)

## Contacts

If you experience problems running the software:

- We encourage using the
  [issue tracker on GitHub](https://github.com/cbg-ethz/lollipop/issues)
- For further enquiries, you can also contact the
  [V-pipe Dev Team](https://cbg-ethz.github.io/V-pipe/contact/)
- You can contact the publication’s corresponding author
