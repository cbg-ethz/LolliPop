import subprocess


def test_workflow_multiple_locs():
    """
    Tests workflow with multiple locations,
    and with multiple cores.
    """
    # data: its handled with LFS

    # dummy run
    subprocess.check_call(["lollipop", "--version"])

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
            "--deconv-config=presets/deconv_linear.yaml",
            "--filters=filters_preprint.yaml",
            "--location=Zürich (ZH)",
            "--location=Genève (GE)",
            "--seed=42",
            "preprint/data/tallymut_line_full.tsv.zst",
        ]
    )


def test_workflow_one_loc():
    """
    Tests workflow with only one location,
    and with only one core.
    """
    # data: its handled with LFS

    # dummy run
    subprocess.check_call(["lollipop", "--version"])

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
            "--deconv-config=presets/deconv_linear.yaml",
            "--filters=filters_preprint.yaml",
            "--location=Zürich (ZH)",
            "--seed=42",
            "preprint/data/tallymut_line_full.tsv.zst",
        ]
    )


def test_workflow_one_loc():
    """
    Tests workflow with only one location,
    and with only one core.
    """
    # data: its handled with LFS

    # dummy run
    subprocess.check_call(["lollipop", "--version"])

    # do fast test from preprint
    subprocess.check_call(
        [
            "lollipop",
            "deconvolute",
            "--n-cores=1",
            "--par-bar=1",
            "--output=test_results.csv",
            "--out-json=test_results.json",
            "--fmt-columns",
            "--variants-config=config_preprint.yaml",
            "--deconv-config=presets/deconv_linear.yaml",
            "--filters=filters_preprint.yaml",
            "--location=Zürich (ZH)",
            "--seed=42",
            "preprint/data/tallymut_line_full.tsv.zst",
        ]
    )
