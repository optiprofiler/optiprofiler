import logging
import os
import platform
import sys
from importlib.metadata import PackageNotFoundError, version


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
