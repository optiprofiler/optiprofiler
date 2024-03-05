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
    PLAIN = 'plain'
    RANDOMIZE_X0 = 'randomize_x0'
    REGULARIZE = 'regularize'
    TOUGH = 'tough'
    TRUNCATE = 'truncate'
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
    M_MIN = 'm_min'
    M_MAX = 'm_max'
    ALL_VARIABLES_CONTINUOUS = 'all_variables_continuous'
    AT_LEAST_ONE_BOUND_CONSTRAINT = 'at_least_one_bound_constraint'
    AT_LEAST_ONE_LINEAR_CONSTRAINT = 'at_least_one_linear_constraint'
    AT_LEAST_ONE_LINEAR_INEQUALITY_CONSTRAINT = 'at_least_one_linear_inequality_constraint'
    AT_LEAST_ONE_LINEAR_EQUALITY_CONSTRAINT = 'at_least_one_linear_equality_constraint'
    AT_LEAST_ONE_NONLINEAR_CONSTRAINT = 'at_least_one_nonlinear_constraint'
    AT_LEAST_ONE_NONLINEAR_INEQUALITY_CONSTRAINT = 'at_least_one_nonlinear_inequality_constraint'
    AT_LEAST_ONE_NONLINEAR_EQUALITY_CONSTRAINT = 'at_least_one_nonlinear_equality_constraint'


class FeatureOption(str, Enum):
    """
    Feature's options.
    """
    DISTRIBUTION = 'distribution'
    MODIFIER = 'modifier'
    N_RUNS = 'n_runs'
    ORDER = 'order'
    PARAMETER = 'parameter'
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


def get_logger(name=None, level=logging.INFO):
    """
    Get a logger.

    Parameters
    ----------
    name : str
        Name of the logger. If ``None``, the root logger is returned.
    level : int
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

        # Attach a console handler (thread-safe).
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter('[%(levelname)-8s] %(message)s'))
        logger.addHandler(handler)
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
        'joblib',
        'matplotlib',
        'numpy',
        'pycutest',
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
