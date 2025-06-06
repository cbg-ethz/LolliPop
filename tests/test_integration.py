import subprocess
from pathlib import PurePath


def test_dummy_run():
    """
    Does it even load and run?
    """

    # dummy run
    subprocess.check_call(["lollipop", "--version"])


def test_workflow_multiple_locs():
    """
    Tests workflow with multiple locations,
    and with multiple cores.
    """
    # data: its handled with LFS

    # do fast test from preprint
    subprocess.check_call(
        [
            "lollipop",
            "deconvolute",
            "--n-cores=2",
            "--par-bar=1",
            "--output=test_results.csv",
            "--out-json=test_results.json",
            "--fmt-columns",
            "--variants-config=config_preprint.yaml",
            ## no var dates: do on all.
            "--deconv-config=presets/deconv_linear.yaml",
            "--filters=filters_preprint.yaml",
            "--location=Zürich (ZH)",
            "--location=Genève (GE)",
            "--seed=42",
            "preprint/data/tallymut_line_full.tsv.zst",
        ]
    )


def test_workflow_one_loc_two_cores():
    """
    Tests workflow with only one location,
    and with only two cores.
    """
    # data: its handled with LFS

    # do fast test from preprint
    subprocess.check_call(
        [
            "lollipop",
            "deconvolute",
            "--n-cores=2",
            "--par-bar=1",
            "--output=test_results_oneloc2cpu.csv",
            "--out-json=test_results_oneloc2cpu.json",
            "--fmt-columns",
            "--variants-config=config_preprint.yaml",
            "--variants-dates=tests/test_var_dates.yaml",  # only the last month
            "--deconv-config=presets/deconv_linear.yaml",
            "--filters=filters_preprint.yaml",
            "--location=Zürich (ZH)",
            "--seed=42",
            "preprint/data/tallymut_line_full.tsv.zst",
        ]
    )


def test_workflow_one_loc_one_core():
    """
    Tests workflow with only one location,
    and with only one core.
    """
    # data: its handled with LFS

    # do fast test from preprint
    subprocess.check_call(
        [
            "lollipop",
            "deconvolute",
            "--n-cores=1",
            "--par-bar=1",
            "--output=test_results_oneloc.csv",
            "--out-json=test_results_oneloc.json",
            "--fmt-columns",
            "--variants-config=config_preprint.yaml",
            "--variants-dates=tests/test_var_dates.yaml",  # only the last month
            "--deconv-config=presets/deconv_linear.yaml",
            "--filters=filters_preprint.yaml",
            "--location=Zürich (ZH)",
            "--seed=42",
            "preprint/data/tallymut_line_full.tsv.zst",
        ]
    )


def test_workflow_no_loc():
    """
    Tests workflow with no location
    given.
    """
    # data: its handled with LFS

    # do fast test from preprint
    subprocess.check_call(
        [
            "lollipop",
            "deconvolute",
            "--n-cores=1",
            "--output=test_results_oneloc.csv",
            "--out-json=test_results_oneloc.json",
            "--variants-config=config_preprint.yaml",
            "--variants-dates=tests/test_var_dates.yaml",  # only the last month
            "--deconv-config=presets/deconv_linear.yaml",
            "--seed=42",
            "preprint/data/tallymut_line_full.tsv.zst",
        ]
    )


def test_workflow_auto_no_loc():
    """
    Tests workflow to autodetect no_loc.
    """

    datapath = PurePath("tests/test_auto_no_loc")

    # do fast test from preprint
    subprocess.check_call(
        [
            "lollipop",
            "deconvolute",
            "--n-cores=1",
            "--output=test_results_auto_no_loc.csv",
            "--out-json=test_results_auto_no_loc.json",
            "--variants-config",
            datapath / "config.yaml",
            "--namefield=mutation",
            "--deconv-config=presets/deconv_linear.yaml",
            "--seed=42",
            datapath / "tallymut.tsv",
        ]
    )
