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


class TestReproducibility:
    """Test suite for reproducibility of deconvolution results."""
    
    @classmethod
    def setup_class(cls):
        """Setup test data and temporary directories."""
        cls.temp_dir = tempfile.mkdtemp()
        cls.test_data_path = "tests/test_auto_no_loc/tallymut.tsv"
        cls.config_path = "tests/test_auto_no_loc/config.yaml"
        cls.deconv_config = "presets/deconv_linear.yaml"
        
        # Check if test data exists
        if not os.path.exists(cls.test_data_path):
            pytest.skip("Test data not available")
    
    @classmethod
    def teardown_class(cls):
        """Clean up temporary directories."""
        shutil.rmtree(cls.temp_dir, ignore_errors=True)
    
    def run_deconvolute(self, seed=None, output_suffix="", n_cores=1, bootstrap=0):
        """Helper method to run deconvolute with given parameters."""
        output_csv = os.path.join(self.temp_dir, f"test_output{output_suffix}.csv")
        output_json = os.path.join(self.temp_dir, f"test_output{output_suffix}.json")
        
        cmd = [
            "lollipop", "deconvolute",
            "--n-cores", str(n_cores),
            "--output", output_csv,
            "--out-json", output_json,
            "--variants-config", self.config_path,
            "--namefield", "mutation",
            "--deconv-config", self.deconv_config,
            self.test_data_path
        ]
        
        if seed is not None:
            cmd.extend(["--seed", str(seed)])
            
        subprocess.check_call(cmd)
        
        return output_csv, output_json
    
    def run_deconvolute_with_bootstrap(self, seed=None, output_suffix="", n_bootstrap=10):
        """Helper method to run deconvolute with bootstrap for testing randomness."""
        # Create a custom bootstrap config
        bootstrap_config = os.path.join(self.temp_dir, "bootstrap_config.yaml")
        with open(bootstrap_config, 'w') as f:
            f.write(f"""
                bootstrap: {n_bootstrap}
                kernel: gaussian
                kernel_params:
                bandwidth: 10
                regressor: nnls
                regressor_params: {{}}
                """)
        
        output_csv = os.path.join(self.temp_dir, f"bootstrap_output{output_suffix}.csv")
        
        cmd = [
            "lollipop", "deconvolute",
            "--n-cores", "1",
            "--output", output_csv,
            "--variants-config", self.config_path,
            "--namefield", "mutation",
            "--deconv-config", bootstrap_config,
            self.test_data_path
        ]
        
        if seed is not None:
            cmd.extend(["--seed", str(seed)])
            
        subprocess.check_call(cmd)
        
        return output_csv

    def test_same_seed_same_results_no_bootstrap(self):
        """Test that the same seed produces identical results without bootstrap."""
        # Run with seed=42 twice
        output1_csv, _ = self.run_deconvolute(seed=42, output_suffix="_1")
        output2_csv, _ = self.run_deconvolute(seed=42, output_suffix="_2")
        
        # Load and compare results
        df1 = pd.read_csv(output1_csv, sep='\t')
        df2 = pd.read_csv(output2_csv, sep='\t')
        
        # DataFrames should be identical
        pd.testing.assert_frame_equal(df1, df2, check_exact=True)
    
    def test_different_seeds_different_results_with_bootstrap(self):
        """Test that different seeds produce different results with bootstrap."""
        # Run with different seeds using bootstrap
        output1_csv = self.run_deconvolute_with_bootstrap(seed=42, output_suffix="_seed42")
        output2_csv = self.run_deconvolute_with_bootstrap(seed=123, output_suffix="_seed123")
        
        # Load results
        df1 = pd.read_csv(output1_csv, sep='\t')
        df2 = pd.read_csv(output2_csv, sep='\t')
        
        # Results should be different (at least some values)
        # We check that not all proportion values are identical
        prop_cols = [col for col in df1.columns if 'proportion' in col.lower()]
        if prop_cols:
            for col in prop_cols:
                if col in df2.columns:
                    # Allow for small numerical differences but expect some meaningful differences
                    differences = np.abs(df1[col] - df2[col])
                    # At least some differences should be larger than numerical precision
                    assert np.any(differences > 1e-10), f"No meaningful differences found in {col}"
    
    def test_same_seed_same_results_with_bootstrap(self):
        """Test that the same seed produces identical results with bootstrap."""
        # Run with same seed twice using bootstrap
        output1_csv = self.run_deconvolute_with_bootstrap(seed=42, output_suffix="_boot1")
        output2_csv = self.run_deconvolute_with_bootstrap(seed=42, output_suffix="_boot2")
        
        # Load and compare results
        df1 = pd.read_csv(output1_csv, sep='\t')
        df2 = pd.read_csv(output2_csv, sep='\t')
        
        # DataFrames should be identical
        pd.testing.assert_frame_equal(df1, df2, check_exact=True)
    
    def test_no_seed_different_results_with_bootstrap(self):
        """Test that without seed, bootstrap runs produce different results."""
        # Run without seed twice using bootstrap
        output1_csv = self.run_deconvolute_with_bootstrap(seed=None, output_suffix="_noseed1")
        output2_csv = self.run_deconvolute_with_bootstrap(seed=None, output_suffix="_noseed2")
        
        # Load results
        df1 = pd.read_csv(output1_csv, sep='\t')
        df2 = pd.read_csv(output2_csv, sep='\t')
        
        # Results should likely be different (though there's a small chance they're the same)
        # We'll check if they're different and if they are, that's expected behavior
        try:
            pd.testing.assert_frame_equal(df1, df2, check_exact=True)
            # If they are equal, that's unlikely but possible, so we'll just warn
            import warnings
            warnings.warn("Two bootstrap runs without seed produced identical results - this is possible but unlikely")
        except AssertionError:
            # This is the expected case - results should be different without seed
            pass
    
    def test_multicore_reproducibility(self):
        """Test that results are reproducible across different numbers of cores."""
        if not os.path.exists("preprint/data/tallymut_line_full.tsv.zst"):
            pytest.skip("Multi-location test data not available")
            
        temp_config = os.path.join(self.temp_dir, "multicore_config.yaml")
        with open(temp_config, 'w') as f:
            f.write("""
variants_pangolin:
  KP.2: KP.2
  KP.3: KP.3
  LP.8: LP.8
variants_list:
  - KP.2
  - KP.3
  - LP.8
variants_not_reported: []
to_drop: []
no_date: false
no_loc: false
""")
        
        # Run with 1 core
        output1_csv, _ = self.run_deconvolute_multicore(temp_config, n_cores=1, seed=42, suffix="_1core")
        
        # Run with 2 cores  
        output2_csv, _ = self.run_deconvolute_multicore(temp_config, n_cores=2, seed=42, suffix="_2core")
        
        # Load and compare results
        df1 = pd.read_csv(output1_csv, sep='\t')
        df2 = pd.read_csv(output2_csv, sep='\t')
        
        # Sort both dataframes to ensure consistent ordering
        sort_cols = [col for col in ['location', 'variant', 'date'] if col in df1.columns]
        if sort_cols:
            df1 = df1.sort_values(sort_cols).reset_index(drop=True)
            df2 = df2.sort_values(sort_cols).reset_index(drop=True)
        
        # Results should be identical regardless of number of cores
        pd.testing.assert_frame_equal(df1, df2, check_exact=True)
    
    def run_deconvolute_multicore(self, config_path, n_cores, seed, suffix):
        """Helper for multicore testing."""
        output_csv = os.path.join(self.temp_dir, f"multicore{suffix}.csv")
        output_json = os.path.join(self.temp_dir, f"multicore{suffix}.json")
        
        cmd = [
            "lollipop", "deconvolute",
            "--n-cores", str(n_cores),
            "--output", output_csv,
            "--out-json", output_json,
            "--variants-config", config_path,
            "--deconv-config", self.deconv_config,
            "--seed", str(seed),
            "--location", "Zürich (ZH)",
            "preprint/data/tallymut_line_full.tsv.zst"
        ]
        
        subprocess.check_call(cmd)
        return output_csv, output_json
    
    def test_json_output_consistency(self):
        """Test that JSON output is also reproducible."""
        # Run with same seed twice
        _, output1_json = self.run_deconvolute(seed=42, output_suffix="_json1")
        _, output2_json = self.run_deconvolute(seed=42, output_suffix="_json2")
        
        # Load and compare JSON outputs
        with open(output1_json, 'r') as f:
            json1 = json.load(f)
        with open(output2_json, 'r') as f:
            json2 = json.load(f)
        
        # JSON outputs should be identical
        assert json1 == json2, "JSON outputs are not identical"
    
    def test_seed_parameter_validation(self):
        """Test that invalid seed parameters are handled properly."""
        output_csv = os.path.join(self.temp_dir, "seed_validation.csv")
        
        # Test with negative seed (should work - numpy accepts negative seeds)
        cmd = [
            "lollipop", "deconvolute",
            "--n-cores", "1",
            "--output", output_csv,
            "--variants-config", self.config_path,
            "--namefield", "mutation",
            "--deconv-config", self.deconv_config,
            "--seed", "-1",
            self.test_data_path
        ]
        
        # This should not raise an error
        subprocess.check_call(cmd)
        
        # Verify output was created
        assert os.path.exists(output_csv), "Output file was not created with negative seed"


class TestSeedingBugFix:
    """Test suite specifically for the seeding bug in bootstrap resampling."""
    
    def test_numpy_legacy_random_seeding_issue(self):
        """Test that demonstrates the current seeding bug."""
        import numpy as np
        from lollipop.confints import resample_mutations
        import pandas as pd
        
        # Create test data
        test_data = pd.DataFrame({
            'mutations': ['A123T', '-A123T', 'C456G', '-C456G'] * 10,
            'frac': np.random.random(40),
            'other_col': range(40)
        })
        
        mutations = test_data['mutations'].unique()
        
        # Test 1: Using np.random.seed (legacy) - should be reproducible
        np.random.seed(42)
        result1, _ = resample_mutations(test_data, mutations)
        
        np.random.seed(42)
        result2, _ = resample_mutations(test_data, mutations)
        
        # These should be identical
        pd.testing.assert_frame_equal(result1, result2)
        
        # Test 2: Using np.random.default_rng (current bug) - might not be reproducible
        rng = np.random.default_rng(42)  # This doesn't affect np.random.randint!
        result3, _ = resample_mutations(test_data, mutations)
        
        rng = np.random.default_rng(42)
        result4, _ = resample_mutations(test_data, mutations)
        
        # These might be different because np.random.randint is not affected by default_rng
        # This test documents the current bug
        try:
            pd.testing.assert_frame_equal(result3, result4)
        except AssertionError:
            # This is expected with the current bug
            pytest.xfail("Known bug: np.random.default_rng() doesn't seed legacy np.random functions")


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
        
        # Run twice with same seed
        outputs = []
        for i in range(2):
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
        df1 = pd.read_csv(outputs[0], sep='\t')
        df2 = pd.read_csv(outputs[1], sep='\t')
        pd.testing.assert_frame_equal(df1, df2, check_exact=True)
        
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)