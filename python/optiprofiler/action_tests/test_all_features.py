#!/usr/bin/env python3
"""Action tests for all features.

This script runs benchmark tests for all 11 features on different OS and Python versions.
These tests are run by GitHub Actions workflows, not by pytest.
"""
import sys
import os
import platform
import time
import matplotlib
matplotlib.use('Agg')

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

from optiprofiler import benchmark
from optiprofiler.action_tests.solvers import SOLVERS, SOLVER_NAMES, UNCONSTRAINED_SOLVERS, UNCONSTRAINED_SOLVER_NAMES


def get_plibs():
    """Get problem libraries based on OS.
    
    Note: pycutest is only available on Linux in CI (macOS ARM64 has installation issues).
    """
    if sys.platform.startswith('linux'):
        return ['s2mpj', 'pycutest', 'custom']
    else:
        return ['s2mpj', 'custom']


def get_benchmark_id():
    """Generate a benchmark ID based on OS and Python version."""
    os_name = platform.system().lower()
    py_version = f"py{sys.version_info.major}{sys.version_info.minor}"
    return f"os_{os_name}_{py_version}"


def run_feature_test(feature_name, benchmark_id=None):
    """Run a benchmark test for a specific feature."""
    if benchmark_id is None:
        benchmark_id = get_benchmark_id()
    
    options = {
        'solver_names': SOLVER_NAMES,
        'solver_verbose': 2,
        'feature_name': feature_name,
        'max_eval_factor': 100,
        'ptype': 'ubln',
        'mindim': 11,
        'maxdim': 11,
        'maxb': 15,
        'maxlcon': 3,
        'maxnlcon': 3,
        'maxcon': 5,
        'benchmark_id': benchmark_id,
        'max_tol_order': 3,
        'plibs': get_plibs(),
    }
    
    print(f"\n{'='*60}")
    print(f"Testing feature: {feature_name}")
    print(f"{'='*60}\n")
    
    start_time = time.time()
    benchmark(SOLVERS, **options)
    elapsed = time.time() - start_time
    print(f"\nTime for {feature_name} test: {elapsed:.2f} seconds\n")


def test_plain(benchmark_id=None):
    """Test plain feature."""
    run_feature_test('plain', benchmark_id)


def test_perturbed_x0(benchmark_id=None):
    """Test perturbed_x0 feature."""
    run_feature_test('perturbed_x0', benchmark_id)


def test_noisy(benchmark_id=None):
    """Test noisy feature."""
    run_feature_test('noisy', benchmark_id)


def test_truncated(benchmark_id=None):
    """Test truncated feature."""
    run_feature_test('truncated', benchmark_id)


def test_permuted(benchmark_id=None):
    """Test permuted feature."""
    run_feature_test('permuted', benchmark_id)


def test_linearly_transformed(benchmark_id=None):
    """Test linearly_transformed feature."""
    run_feature_test('linearly_transformed', benchmark_id)


def test_random_nan(benchmark_id=None):
    """Test random_nan feature."""
    run_feature_test('random_nan', benchmark_id)


def test_unrelaxable_constraints(benchmark_id=None):
    """Test unrelaxable_constraints feature."""
    run_feature_test('unrelaxable_constraints', benchmark_id)


def test_nonquantifiable_constraints(benchmark_id=None):
    """Test nonquantifiable_constraints feature."""
    run_feature_test('nonquantifiable_constraints', benchmark_id)


def test_quantized(benchmark_id=None):
    """Test quantized feature."""
    run_feature_test('quantized', benchmark_id)


def test_custom(benchmark_id=None):
    """Test custom feature."""
    run_feature_test('custom', benchmark_id)


def run_all_feature_tests(benchmark_id=None):
    """Run all feature tests."""
    if benchmark_id is None:
        benchmark_id = get_benchmark_id()
    
    features = [
        'plain',
        'perturbed_x0',
        'noisy',
        'truncated',
        'permuted',
        'linearly_transformed',
        'random_nan',
        'unrelaxable_constraints',
        'nonquantifiable_constraints',
        'quantized',
        'custom',
    ]
    
    for feature in features:
        run_feature_test(feature, benchmark_id)


if __name__ == '__main__':
    benchmark_id = sys.argv[1] if len(sys.argv) > 1 else None
    run_all_feature_tests(benchmark_id)
