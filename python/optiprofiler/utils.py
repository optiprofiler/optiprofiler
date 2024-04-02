import logging
import os
import platform
import sys
from enum import Enum
from importlib.metadata import PackageNotFoundError, version


class FeatureName(str, Enum):
    """
    Feature names.
    """
    CUSTOM = 'custom'
    NOISY = 'noisy'
    PERMUTED = 'permuted'
    PLAIN = 'plain'
    PERTURBED_X0 = 'perturbed_x0'
    RANDOM_NAN = 'random_nan'
    TRUNCATED = 'truncated'
    UNRELAXABLE_CONSTRAINTS = 'unrelaxable_constraints'


class ProfileOption(str, Enum):
    """
    Profile's options.
    """
    N_JOBS = 'n_jobs'
    BENCHMARK_ID = 'benchmark_id'
    PROJECT_X0 = 'project_x0'


class CUTEstProblemOption(str, Enum):
    """
    Options for loading a CUTEst problem.
    """
    N_MIN = 'n_min'
    N_MAX = 'n_max'
    M_BOUND_MIN = 'm_bound_min'
    M_BOUND_MAX = 'm_bound_max'
    M_LINEAR_MIN = 'm_linear_min'
    M_LINEAR_MAX = 'm_linear_max'
    M_LINEAR_INEQUALITY_MIN = 'm_linear_inequality_min'
    M_LINEAR_INEQUALITY_MAX = 'm_linear_inequality_max'
    M_LINEAR_EQUALITY_MIN = 'm_linear_equality_min'
    M_LINEAR_EQUALITY_MAX = 'm_linear_equality_max'
    M_NONLINEAR_MIN = 'm_nonlinear_min'
    M_NONLINEAR_MAX = 'm_nonlinear_max'
    M_NONLINEAR_INEQUALITY_MIN = 'm_nonlinear_inequality_min'
    M_NONLINEAR_INEQUALITY_MAX = 'm_nonlinear_inequality_max'
    M_NONLINEAR_EQUALITY_MIN = 'm_nonlinear_equality_min'
    M_NONLINEAR_EQUALITY_MAX = 'm_nonlinear_equality_max'
    ALL_VARIABLES_CONTINUOUS = 'all_variables_continuous'


class FeatureOption(str, Enum):
    """
    Feature's options.
    """
    DISTRIBUTION = 'distribution'
    MODIFIER = 'modifier'
    N_RUNS = 'n_runs'
    PERTURBED_TRAILING_ZEROS = 'perturbed_trailing_zeros'
    RATE_NAN = 'rate_nan'
    SIGNIFICANT_DIGITS = 'significant_digits'
    TYPE = 'type'
    UNRELAXABLE_BOUNDS = 'unrelaxable_bounds'
    UNRELAXABLE_LINEAR_CONSTRAINTS = 'unrelaxable_linear_constraints'
    UNRELAXABLE_NONLINEAR_CONSTRAINTS = 'unrelaxable_nonlinear_constraints'


class NoiseType(str, Enum):
    """
    Noise types.
    """
    ABSOLUTE = 'absolute'
    RELATIVE = 'relative'


class ProblemError(Exception):
    """
    Exception raised when a problem cannot be loaded.
    """


def show_versions():
    """
    Display useful system and dependency information.

    When reporting issues, please include this information.
    """
    # Get system and dependency information.
    sys_info = _get_sys_info()
    deps_info = _get_deps_info()

    # Display information.
    print('System settings')
    print('---------------')
    print('\n'.join(
        f'{k:>{max(map(len, sys_info.keys())) + 1}}: {v}'
        for k, v in sys_info.items()
    ))
    print()
    print('Python dependencies')
    print('-------------------')
    print('\n'.join(
        f'{k:>{max(map(len, deps_info.keys())) + 1}}: {v}'
        for k, v in deps_info.items()
    ))


def get_logger(name, level=logging.INFO):
    """
    Get a logger.

    Parameters
    ----------
    name : str
        Name of the logger. If ``None``, the root logger is returned.
    level : int, optional
        Logging level.

    Returns
    -------
    `logging.Logger`
        Logger with the given name. If a logger with the given name has already
        been created, it is returned instead of creating a new one.
    """
    logger = logging.getLogger(name)
    if len(logger.handlers) == 0:
        logger.setLevel(level)
        formatter = logging.Formatter('[%(levelname)-8s] %(message)s')

        # Attach a console handler.
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    return logger


def _get_sys_info():
    """
    Get system information.

    Returns
    -------
    dict
        System information.
    """
    return {
        'python': sys.version.replace(os.linesep, ' '),
        'executable': sys.executable,
        'machine': platform.platform(),
    }


def _get_deps_info():
    """
    Get dependency information.

    Returns
    -------
    dict
        Dependency information.
    """
    deps = [
        'optiprofiler',
        'matplotlib',
        'numpy',
        'pycutest',
        'pypdf',
        'scipy',
        'setuptools',
        'pip',
    ]
    deps_info = {}
    for module in deps:
        try:
            deps_info[module] = version(module)
        except PackageNotFoundError:
            deps_info[module] = None
    return deps_info
