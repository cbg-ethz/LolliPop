"""
Test to demonstrate and fix the seeding bug in lollipop deconvolute.

The issue is in the _deconvolute_bootstrap_wrapper function where
np.random.default_rng(child_seed) is called but the resample_mutations
function uses the legacy np.random.randint() which is not affected by
the new RNG.
"""

import numpy as np
import pandas as pd
import sys
import os

# Add the lollipop package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

try:
    from lollipop.confints import resample_mutations
except ImportError:
    # Create a minimal version for testing
    def resample_mutations(df_city1, mutations, namefield="mutations"):
        """Simplified version for testing."""
        # This is the problematic line that uses legacy numpy random
        rand_idcs = np.random.randint(0, high=int(len(mutations) / 2), size=int(len(mutations) / 2))
        resamples_counts = np.bincount(rand_idcs, minlength=int(len(mutations) / 2))
        resample_coeff_dict = dict(zip(mutations, np.concatenate([resamples_counts, resamples_counts])))
        df_sampled = df_city1.copy()
        df_sampled.loc[:, "resample_value"] = df_sampled[namefield].map(resample_coeff_dict)
        return df_sampled, rand_idcs


def test_current_seeding_bug():
    """Demonstrate the current seeding bug."""
    # Create test data
    test_data = pd.DataFrame({
        'mutations': ['A123T', '-A123T', 'C456G', '-C456G'] * 5,
        'frac': np.random.random(20),
        'other_col': range(20)
    })
    
    mutations = test_data['mutations'].unique()
    
    print("Testing current seeding approach (buggy):")
    
    # This is what the current code does (incorrectly)
    child_seed = np.random.SeedSequence(42)
    np.random.default_rng(child_seed)  # This doesn't affect np.random.randint!
    result1, indices1 = resample_mutations(test_data, mutations)
    
    child_seed = np.random.SeedSequence(42)
    np.random.default_rng(child_seed)  # This doesn't affect np.random.randint!
    result2, indices2 = resample_mutations(test_data, mutations)
    
    print(f"Indices from run 1: {indices1}")
    print(f"Indices from run 2: {indices2}")
    print(f"Are indices identical? {np.array_equal(indices1, indices2)}")
    
    return np.array_equal(indices1, indices2)


def test_correct_seeding():
    """Demonstrate the correct seeding approach."""
    # Create test data
    test_data = pd.DataFrame({
        'mutations': ['A123T', '-A123T', 'C456G', '-C456G'] * 5,
        'frac': np.random.random(20),
        'other_col': range(20)
    })
    
    mutations = test_data['mutations'].unique()
    
    print("\nTesting correct seeding approach:")
    
    # Correct approach: use legacy np.random.seed
    child_seed = np.random.SeedSequence(42)
    np.random.seed(child_seed.entropy)  # Use entropy from SeedSequence
    result1, indices1 = resample_mutations(test_data, mutations)
    
    child_seed = np.random.SeedSequence(42)
    np.random.seed(child_seed.entropy)  # Use entropy from SeedSequence
    result2, indices2 = resample_mutations(test_data, mutations)
    
    print(f"Indices from run 1: {indices1}")
    print(f"Indices from run 2: {indices2}")
    print(f"Are indices identical? {np.array_equal(indices1, indices2)}")
    
    return np.array_equal(indices1, indices2)


def test_alternative_correct_seeding():
    """Alternative correct seeding using generated integers."""
    # Create test data  
    test_data = pd.DataFrame({
        'mutations': ['A123T', '-A123T', 'C456G', '-C456G'] * 5,
        'frac': np.random.random(20),
        'other_col': range(20)
    })
    
    mutations = test_data['mutations'].unique()
    
    print("\nTesting alternative correct seeding (generate seed integer):")
    
    # Alternative correct approach: generate a seed integer from SeedSequence
    child_seed = np.random.SeedSequence(42)
    seed_int = child_seed.generate_state(1)[0]  # Generate a single integer
    np.random.seed(seed_int)
    result1, indices1 = resample_mutations(test_data, mutations)
    
    child_seed = np.random.SeedSequence(42)
    seed_int = child_seed.generate_state(1)[0]  # Generate a single integer
    np.random.seed(seed_int)
    result2, indices2 = resample_mutations(test_data, mutations)
    
    print(f"Indices from run 1: {indices1}")
    print(f"Indices from run 2: {indices2}")
    print(f"Are indices identical? {np.array_equal(indices1, indices2)}")
    
    return np.array_equal(indices1, indices2)


if __name__ == "__main__":
    print("=" * 60)
    print("TESTING LOLLIPOP SEEDING BUG")
    print("=" * 60)
    
    # Test current (buggy) behavior
    is_reproducible_current = test_current_seeding_bug()
    
    # Test correct approaches
    is_reproducible_correct1 = test_correct_seeding()
    is_reproducible_correct2 = test_alternative_correct_seeding()
    
    print("\n" + "=" * 60)
    print("SUMMARY:")
    print(f"Current approach (buggy): {'PASS' if is_reproducible_current else 'FAIL - Not reproducible'}")
    print(f"Correct approach 1: {'PASS' if is_reproducible_correct1 else 'FAIL'}")
    print(f"Correct approach 2: {'PASS' if is_reproducible_correct2 else 'FAIL'}")
    print("=" * 60)
    
    if not is_reproducible_current:
        print("\n🐛 BUG CONFIRMED: Current seeding does not work with bootstrap!")
        print("   The issue is that np.random.default_rng() doesn't affect np.random.randint()")
    
    if is_reproducible_correct1 and is_reproducible_correct2:
        print("\n✅ SOLUTION VERIFIED: Correct seeding approaches work!")
    
    print("\nTo fix the bug, modify the _deconvolute_bootstrap_wrapper function in")
    print("lollipop/cli/deconvolute.py to use np.random.seed() instead of np.random.default_rng()")
