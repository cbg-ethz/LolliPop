import subprocess


def test_workflow():
    # data: its handled with LFS

    # dummy run
    subprocess.check_call(["lollipop", "--version"])

    # do fast test from preprint
    subprocess.check_call(
        [
            "lollipop",
            "deconvolute",
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
