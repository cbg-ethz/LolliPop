"""Test the reproducibility of deconvolution results with lollipop.deconvolute using seeds."""

import subprocess
import pandas as pd
import numpy as np
import tempfile
import os
import json
import pytest
from pathlib import Path
import shutil

def test_basic_reproducibility_quick():
    """Quick test that can run without heavy test data."""
    # This test uses the minimal test data in test_auto_no_loc
    temp_dir = tempfile.mkdtemp()
    try:
        test_data_path = "tests/test_auto_no_loc/tallymut.tsv"
        config_path = "tests/test_auto_no_loc/config.yaml"
        deconv_config = "presets/deconv_linear.yaml"
        
        if not os.path.exists(test_data_path):
            pytest.skip("Test data not available")
        
        # Run 15 times with same seed
        outputs = []
        for i in range(15):
            output_csv = os.path.join(temp_dir, f"quick_test_{i}.csv")
            cmd = [
                "lollipop", "deconvolute",
                "--n-cores", "1",
                "--output", output_csv,
                "--variants-config", config_path,
                "--namefield", "mutation",
                "--deconv-config", deconv_config,
                "--seed", "42",
                test_data_path
            ]
            subprocess.check_call(cmd)
            outputs.append(output_csv)
        
        # Compare outputs
        df_first = pd.read_csv(outputs[0], sep='\t')

        for i in range(1, 15):
            df_next = pd.read_csv(outputs[i], sep='\t')
            try:
                pd.testing.assert_frame_equal(df_first, df_next, check_exact=True)
            except AssertionError as e:
                print(f"Difference found between run 0 and run {i}:")
                print(df_first.compare(df_next))
                df_first.to_csv(os.path.join(temp_dir, "df_first_debug.csv"), index=False)
                df_next.to_csv(os.path.join(temp_dir, f"df_next_{i}_debug.csv"), index=False)
                raise e
    finally:
        print(f"Cleaning up temporary directory: {temp_dir}")
    #    shutil.rmtree(temp_dir, ignore_errors=True)
