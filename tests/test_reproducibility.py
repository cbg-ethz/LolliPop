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

def test_reproducibility_cowwid():
    """
    Test reproducibility of deconvolution results under cowwid production settings.
    
    This test verifies that:
    1. Random number seeding works correctly when processing multiple locations in parallel
    2. Each location's results are reproducible across runs with the same seed
    3. Multi-location deconvolution maintains deterministic behavior
    4. Parallel processing does works with random number generation
    
    Uses the full preprint dataset with multiple Swiss cities to test complex
    scenarios where deconvolution processes multiple locations simultaneously.
    """
    temp_dir = tempfile.mkdtemp()
    try:
        # Use the compressed LFS file directly from preprint/data
        test_data_path = "preprint/data/tallymut_line_full.tsv.zst"
        config_path = "tests/test_reproducibility/config.yaml"
        deconv_config = "presets/deconv_bootstrap_cowwid.yaml"
        
        if not os.path.exists(test_data_path):
            pytest.skip("Test data not available")
        
        # Run x times with same seed to verify seeding works with multiple locations
        n_runs = 2
        outputs = []
        for i in range(n_runs):
            output_csv = os.path.join(temp_dir, f"multilocations_test_{i}.csv")
            cmd = [
                "lollipop", "deconvolute",
                "--n-cores", "3",  
                "--output", output_csv,
                "--variants-config", config_path,
                "--deconv-config", deconv_config,
                # 3 Locations are specified in config.yaml via locations_list
                "--seed", "123", 
                test_data_path
            ]
            subprocess.check_call(cmd)
            outputs.append(output_csv)
        
        # Compare outputs for exact equality across all locations
        df_first = pd.read_csv(outputs[0], sep='\t')
    

        for i in range(1, n_runs):
            df_next = pd.read_csv(outputs[i], sep='\t')
            try:
                pd.set_option('display.float_format', lambda x: '%.20f' % x)
                pd.testing.assert_frame_equal(df_first, df_next, check_exact=True)
            except AssertionError as e:
                print(f"Multi-location seeding difference found between run 0 and run {i}:")
                
                # Create a more readable comparison by showing the actual rows that differ
                comparison = df_first.compare(df_next, align_axis=1)
                if not comparison.empty:
                    print("\nDetailed differences (showing location, variant, date context):")
                    print("=" * 80)
                    
                    # Get indices of differing rows
                    differing_indices = comparison.index.tolist()
                    
                    for idx in differing_indices[:10]:  # Show first 10 differences
                        print(f"\nRow {idx}:")
                        print(f"Location: {df_first.loc[idx, 'location']}")
                        print(f"Variant:  {df_first.loc[idx, 'variant']}")
                        print(f"Date:     {df_first.loc[idx, 'date']}")
                        
                        # Show the differing columns for this row
                        row_comparison = comparison.loc[idx]
                        for col in row_comparison.index.levels[0]:  # Get column names
                            if pd.notna(row_comparison[col, 'self']) or pd.notna(row_comparison[col, 'other']):
                                val_first = row_comparison[col, 'self'] if pd.notna(row_comparison[col, 'self']) else df_first.loc[idx, col]
                                val_next = row_comparison[col, 'other'] if pd.notna(row_comparison[col, 'other']) else df_next.loc[idx, col]
                                print(f"  {col}:")
                                print(f"    Run 0: {val_first}")
                                print(f"    Run {i}: {val_next}")
                        print("-" * 40)
                    
                    if len(differing_indices) > 10:
                        print(f"\n... and {len(differing_indices) - 10} more differences")
                
                raise AssertionError(f"Multi-location seeding reproducibility failed: {e}")
    except subprocess.CalledProcessError as e:
        pytest.fail(f"Deconvolution command failed: {e}")
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)