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
    REGULARIZED = 'regularized'
    TOUGH = 'tough'
    TRUNCATED = 'truncated'


class ProfileOptionKey(str, Enum):
    """
    Profile's options.
    """
    N_JOBS = 'n_jobs'
    BENCHMARK_ID = 'benchmark_id'


class CUTEstProblemOptionKey(str, Enum):
    """
    Options for loading a CUTEst problem.
    """
    N_MIN = 'n_min'
    N_MAX = 'n_max'
    M_MIN = 'm_min'
    M_MAX = 'm_max'


class FeatureOptionKey(str, Enum):
    """
    Feature's options.
    """
    DISTRIBUTION = 'distribution'
    MODIFIER = 'modifier'
    N_RUNS = 'n_runs'
    ORDER = 'order'
    PARAMETER = 'parameter'
    RATE_ERROR = 'rate_error'
    RATE_NAN = 'rate_nan'
    SIGNIFICANT_DIGITS = 'significant_digits'
    TYPE = 'type'


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
    pass


def show_versions():
    """
    Display useful system and dependency information.

    When reporting issues, please include this information.
    """
    print('System settings')
    print('---------------')
    sys_info = _get_sys_info()
    print('\n'.join(f'{k:>{max(map(len, sys_info.keys())) + 1}}: {v}' for k, v in sys_info.items()))

    print()
    print('Python dependencies')
    print('-------------------')
    deps_info = _get_deps_info()
    print('\n'.join(f'{k:>{max(map(len, deps_info.keys())) + 1}}: {v}' for k, v in deps_info.items()))


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
    deps = ['OptiProfiler', 'joblib', 'matplotlib', 'numpy', 'pycutest', 'setuptools', 'pip']
    deps_info = {}
    for module in deps:
        try:
            deps_info[module] = version(module)
        except PackageNotFoundError:
            deps_info[module] = None
    return deps_info
