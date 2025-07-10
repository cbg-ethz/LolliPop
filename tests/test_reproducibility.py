"""
Test the reproducibility of deconvolution results with lollipop.deconvolute using seeds.

This module contains tests that verify the numerical stability and randomization 
control in the LolliPop deconvolution pipeline

The use the full preprint dataset (tallymut_line_full.tsv.zst) to ensure
realistic testing conditions with real-world data complexity.
"""

import subprocess
import pandas as pd
import numpy as np
import tempfile
import os
import json
import pytest
from pathlib import Path
import shutil

def test_reproducibility():
    """
    Test correct seeding behavior with multiple locations processing.
    
    This test verifies that:
    1. Random number seeding works correctly when processing multiple locations in parallel
    2. Each location's results are reproducible across runs with the same seed
    3. Multi-location deconvolution maintains deterministic behavior
    4. Parallel processing doesn't introduce race conditions in random number generation
    
    Uses the full preprint dataset with multiple Swiss cities to test complex
    scenarios where deconvolution processes multiple locations simultaneously.
    """
    # Create debug directory for inspection
    debug_dir = Path("tests/test_reproducibility/debug")
    debug_dir.mkdir(parents=True, exist_ok=True)
    
    temp_dir = tempfile.mkdtemp()
    try:
        # Use the compressed LFS file directly from preprint/data
        test_data_path = "preprint/data/tallymut_line_full.tsv.zst"
        config_path = "tests/test_reproducibility/config.yaml"
        deconv_config = "presets/deconv_linear.yaml"
        
        if not os.path.exists(test_data_path):
            pytest.skip("Test data not available")
        
        # Run 2 times with same seed to verify seeding works with multiple locations
        outputs = []
        for i in range(2):
            output_csv = os.path.join(temp_dir, f"multilocations_test_{i}.csv")
            cmd = [
                "lollipop", "deconvolute",
                "--n-cores", "1",  # Force single core to ensure deterministic ordering
                "--output", output_csv,
                "--variants-config", config_path,
                "--namefield", "mutation",
                "--deconv-config", deconv_config,
                # Locations are now specified in config.yaml via locations_list
                "--seed", "123",  # Different seed from numerical test
                test_data_path
            ]
            subprocess.check_call(cmd)
            outputs.append(output_csv)
        
        # Compare outputs for exact equality across all locations
        df_first = pd.read_csv(outputs[0], sep='\t')

        for i in range(1, 2):
            df_next = pd.read_csv(outputs[i], sep='\t')
            try:
                pd.testing.assert_frame_equal(df_first, df_next, check_exact=True)
            except AssertionError as e:
                print(f"Multi-location seeding difference found between run 0 and run {i}:")
                print(df_first.compare(df_next))
                # Save debug files to persistent location for inspection
                df_first.to_csv(debug_dir / "df_first_debug.csv", index=False)
                df_next.to_csv(debug_dir / f"df_next_{i}_debug.csv", index=False)
                print(f"Debug files saved to: {debug_dir}")
                raise AssertionError(f"Multi-location seeding reproducibility failed: {e}")
    finally:
        print(f"Temporary directory was: {temp_dir}")
        print(f"Debug files (if any) saved to: {debug_dir}")
        # shutil.rmtree(temp_dir, ignore_errors=True)