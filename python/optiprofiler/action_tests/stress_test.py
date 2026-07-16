#!/usr/bin/env python3
"""Stress test for OptiProfiler.

Tests randomly selected high-dimensional (unconstrained) problems from the
bundled S2MPJ library. This test is run by GitHub Actions workflows, not by
pytest.
"""
import sys
import os
import platform
import matplotlib
matplotlib.use('Agg')
import numpy as np
from datetime import datetime

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

from optiprofiler import benchmark
from optiprofiler.action_tests.solvers import UNCONSTRAINED_SOLVERS, UNCONSTRAINED_SOLVER_NAMES


DEFAULT_MAX_PROBLEMS_PER_PLIB = 1

STRESS_PROBLEM_NAMES = {
    's2mpj': [
        'SPMSRTLS',
        'WOODS',
    ],
}


def get_plibs():
    """Return the bundled library exercised by core stress tests."""
    return ['s2mpj']


def get_benchmark_id():
    """Generate a benchmark ID based on OS and Python version."""
    os_name = platform.system().lower()
    py_version = f"py{sys.version_info.major}{sys.version_info.minor}"
    return f"os_{os_name}_{py_version}"


def get_max_problems_per_plib():
    """Get the maximum number of stress problems to run from each library."""
    env_value = os.environ.get('OPTIPROFILER_STRESS_MAX_PROBLEMS_PER_PLIB')
    if env_value is None:
        return DEFAULT_MAX_PROBLEMS_PER_PLIB
    try:
        max_problems = int(env_value)
    except ValueError as exc:
        raise ValueError('OPTIPROFILER_STRESS_MAX_PROBLEMS_PER_PLIB must be an integer.') from exc
    if max_problems < 1:
        raise ValueError('OPTIPROFILER_STRESS_MAX_PROBLEMS_PER_PLIB must be positive.')
    return max_problems


def get_stress_problem_names(plib, seed):
    """Select a small deterministic subset of high-dimensional stress problems."""
    problem_names = STRESS_PROBLEM_NAMES[plib]
    max_problems = min(get_max_problems_per_plib(), len(problem_names))
    rng = np.random.default_rng(seed + sum(ord(c) for c in plib))
    idx = sorted(rng.permutation(len(problem_names))[:max_problems])
    return [problem_names[i] for i in idx]


def stress_test(benchmark_id=None):
    """Run stress test on high-dimensional problems.
    
    Tests randomly selected high-dimensional (2000-10000) unconstrained problems.
    """
    if benchmark_id is None:
        benchmark_id = get_benchmark_id()
    
    # Generate seed based on current date (year and week)
    now = datetime.now()
    seed = 100 * (now.year % 100) + now.isocalendar()[1]
    print(f"Seed: {seed}\n")
    
    # Set options
    common_options = {
        'seed': seed,
        'solver_names': UNCONSTRAINED_SOLVER_NAMES,
        'ptype': 'u',
        'mindim': 2000,
        'maxdim': 10000,
        'benchmark_id': benchmark_id,
    }
    
    # Adjust max_eval_factor based on OS
    if sys.platform.startswith('linux'):
        common_options['max_eval_factor'] = 0.8
    else:
        common_options['max_eval_factor'] = 0.2
    
    for plib in get_plibs():
        options = {
            **common_options,
            'plibs': [plib],
            'problem_names': get_stress_problem_names(plib, seed),
            'benchmark_id': f'{benchmark_id}_{plib}',
        }

        print(f"Running stress test with options:")
        for key, value in options.items():
            print(f"  {key}: {value}")
        print()

        # Run benchmark
        benchmark(UNCONSTRAINED_SOLVERS, **options)


if __name__ == '__main__':
    benchmark_id = sys.argv[1] if len(sys.argv) > 1 else None
    stress_test(benchmark_id)
