#!/usr/bin/env python3
"""Random test for OptiProfiler.

Tests OptiProfiler with randomly generated options.
This test is run by GitHub Actions workflows, not by pytest.
"""
import sys
import os
import platform
import numpy as np
from datetime import datetime

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

from optiprofiler import benchmark
from optiprofiler.action_tests.solvers import SOLVERS, SOLVER_NAMES


def get_plibs():
    """Get problem libraries based on OS."""
    if sys.platform.startswith('linux') or sys.platform == 'darwin':
        return ['s2mpj', 'pycutest']
    else:
        return ['s2mpj']


def get_benchmark_id():
    """Generate a benchmark ID based on OS and Python version."""
    os_name = platform.system().lower()
    py_version = f"py{sys.version_info.major}{sys.version_info.minor}"
    return f"os_{os_name}_{py_version}"


def random_test(benchmark_id=None):
    """Run random test with randomly generated options."""
    if benchmark_id is None:
        benchmark_id = get_benchmark_id()
    
    # Generate seed based on current date
    now = datetime.now()
    seed = 100 * (now.year % 100) + now.isocalendar()[1] + now.day
    print(f"Seed: {seed}\n")
    
    rng = np.random.default_rng(seed)
    
    # Build random options
    options = {
        'n_jobs': rng.integers(1, 6),
        'seed': seed,
        'benchmark_id': benchmark_id,
    }
    
    # Random solver selection (2 or 3 solvers)
    num_solvers = rng.integers(2, 4)
    solver_indices = rng.permutation(len(SOLVERS))[:num_solvers].tolist()
    solvers = [SOLVERS[i] for i in solver_indices]
    options['solver_names'] = [SOLVER_NAMES[i] for i in solver_indices]
    
    # Random profile options
    options['errorbar_type'] = rng.choice(['minmax', 'meanstd'])
    options['hist_aggregation'] = rng.choice(['min', 'mean', 'max'])
    options['max_tol_order'] = rng.integers(1, 17)
    options['max_eval_factor'] = rng.random() * 30
    options['project_x0'] = rng.random() < 0.5
    options['run_plain'] = rng.random() < 0.5
    options['score_only'] = rng.random() < 0.5
    options['summarize_performance_profiles'] = rng.random() < 0.5
    options['summarize_data_profiles'] = rng.random() < 0.5
    options['summarize_log_ratio_profiles'] = len(solvers) == 2 and rng.random() < 0.5
    options['summarize_output_based_profiles'] = rng.random() < 0.5
    options['solver_verbose'] = rng.integers(0, 3)
    options['semilogx'] = rng.random() < 0.5
    options['normalized_scores'] = rng.random() < 0.5
    
    # Random line styles and colors
    colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    options['line_colors'] = list(rng.choice(colors, len(solvers), replace=False))
    options['bar_colors'] = list(rng.choice(colors, len(solvers), replace=False))
    
    styles = ['-', '--', '-.', ':']
    options['line_styles'] = list(rng.choice(styles, len(solvers), replace=True))
    options['line_widths'] = (rng.random(len(solvers)) * 2 + 0.5).tolist()
    
    # Random hist plot mode
    hist_choices = ['none', 'sequential', 'parallel']
    options['draw_hist_plots'] = rng.choice(hist_choices)
    
    # Random feature
    feature_choices = [
        'plain', 'perturbed_x0', 'noisy', 'truncated', 'permuted',
        'linearly_transformed', 'random_nan', 'unrelaxable_constraints',
        'nonquantifiable_constraints', 'quantized'
    ]
    feature_name = rng.choice(feature_choices)
    options['feature_name'] = feature_name
    
    # Feature-specific options
    if feature_name == 'plain':
        options['n_runs'] = 1
    elif feature_name == 'perturbed_x0':
        options['n_runs'] = rng.integers(1, 4)
        options['perturbation_level'] = rng.random()
        options['distribution'] = rng.choice(['gaussian', 'spherical'])
    elif feature_name == 'noisy':
        options['n_runs'] = rng.integers(1, 4)
        options['noise_level'] = rng.random()
        options['distribution'] = rng.choice(['gaussian', 'uniform'])
        options['noise_type'] = rng.choice(['absolute', 'relative', 'mixed'])
    elif feature_name == 'truncated':
        options['n_runs'] = rng.integers(1, 4)
        options['significant_digits'] = rng.integers(1, 11)
        options['perturbed_trailing_digits'] = rng.random() < 0.5
    elif feature_name == 'linearly_transformed':
        options['n_runs'] = rng.integers(1, 4)
        options['rotated'] = rng.random() < 0.5
        options['condition_factor'] = rng.random() * 100
    elif feature_name == 'random_nan':
        options['n_runs'] = rng.integers(1, 4)
        options['nan_rate'] = rng.random() * 0.5
    elif feature_name == 'unrelaxable_constraints':
        options['n_runs'] = 1
        options['unrelaxable_bounds'] = rng.random() < 0.5
        options['unrelaxable_linear_constraints'] = rng.random() < 0.5
        options['unrelaxable_nonlinear_constraints'] = rng.random() < 0.5
    elif feature_name == 'nonquantifiable_constraints':
        options['n_runs'] = 1
    elif feature_name == 'quantized':
        options['n_runs'] = 1
        options['mesh_size'] = rng.random()
        options['mesh_type'] = rng.choice(['absolute', 'relative'])
        options['ground_truth'] = rng.random() < 0.5
    
    # Problem library options
    options['plibs'] = get_plibs()
    
    # Random problem type
    ptype_choices = ['u', 'b', 'l', 'n', 'ub', 'ul', 'un', 'bl', 'bn', 'ln', 
                     'ubl', 'ubn', 'uln', 'bln', 'ubln']
    options['ptype'] = rng.choice(ptype_choices)
    
    # Random dimension constraints
    options['mindim'] = rng.integers(1, 4)
    options['maxdim'] = options['mindim'] + rng.integers(1, 4)
    
    # Random constraint counts
    options['minb'] = 0
    options['maxb'] = rng.integers(20, 41)
    options['minlcon'] = 0
    options['maxlcon'] = rng.integers(10, 21)
    options['minnlcon'] = 0
    options['maxnlcon'] = rng.integers(10, 21)
    options['mincon'] = 0
    options['maxcon'] = rng.integers(20, 41)
    
    # Random label options
    if rng.random() < 0.5:
        options['xlabel_data_profile'] = ''
    if rng.random() < 0.5:
        options['ylabel_data_profile'] = ''
    if rng.random() < 0.5:
        options['xlabel_performance_profile'] = ''
    if rng.random() < 0.5:
        options['ylabel_performance_profile'] = ''
    if rng.random() < 0.5:
        options['xlabel_log_ratio_profile'] = ''
    if rng.random() < 0.5:
        options['ylabel_log_ratio_profile'] = ''
    
    # Print options
    print("Running random test with options:")
    for key, value in options.items():
        print(f"  {key}: {value}")
    print()
    
    # Run benchmark
    benchmark(solvers, **options)


if __name__ == '__main__':
    benchmark_id = sys.argv[1] if len(sys.argv) > 1 else None
    random_test(benchmark_id)
