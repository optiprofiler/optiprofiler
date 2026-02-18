#!/usr/bin/env python3
"""Stress test for OptiProfiler.

Tests randomly high-dimensional (unconstrained) problems from s2mpj or pycutest.
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
from optiprofiler.action_tests.solvers import UNCONSTRAINED_SOLVERS, UNCONSTRAINED_SOLVER_NAMES


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
    options = {
        'seed': seed,
        'solver_names': UNCONSTRAINED_SOLVER_NAMES,
        'ptype': 'u',
        'mindim': 2000,
        'maxdim': 10000,
        'benchmark_id': benchmark_id,
        'plibs': get_plibs(),
    }
    
    # Adjust max_eval_factor based on OS
    if sys.platform.startswith('linux'):
        options['max_eval_factor'] = 0.8
    else:
        options['max_eval_factor'] = 0.2
    
    print(f"Running stress test with options:")
    for key, value in options.items():
        print(f"  {key}: {value}")
    print()
    
    # Run benchmark
    benchmark(UNCONSTRAINED_SOLVERS, **options)


if __name__ == '__main__':
    benchmark_id = sys.argv[1] if len(sys.argv) > 1 else None
    stress_test(benchmark_id)
