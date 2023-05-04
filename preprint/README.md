# Notebooks

This directory contains the Notebooks that have been used to generate the plots for the preprint.

## Jupyter kernel

To install the `lollipop` kernel necessary for [notebook WwSmoothingKernel.ipynb](WwSmoothingKernel.ipynb),
you can also use conda:

```bash
conda env create -f conda_ipylollipop_env.yaml
conda activate ipylollipop
python -m ipykernel install --user --name=lollipop
```
