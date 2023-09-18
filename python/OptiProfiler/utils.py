import logging
import os
import platform
import sys
from importlib.metadata import PackageNotFoundError, version


def get_logger(name=None, level=logging.INFO):
    """
    Get a logger with the given name and level.

    Parameters
    ----------
    name : str
        Name of the logger.
    level : int
        Logging level.

    Returns
    -------
    logging.Logger
        Logger with the given name and level. If a logger with the given name
        already exists, it is returned instead of creating a new one.
    """
    logger = logging.getLogger(name)

    # Multiple calls to get_logger with the same name will return a reference
    # to the same logger. We do not create handlers if some already exist.
    if len(logger.handlers) == 0:
        logger.setLevel(level)

        # Attach a console handler (thread-safe).
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter('[%(levelname)-8s] %(message)s'))
        logger.addHandler(handler)
    return logger


def show_versions():
    """
    Provide useful information when reporting issues.
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
    deps = ['OptiProfiler', 'numpy', 'setuptools', 'pip']
    deps_info = {}
    for module in deps:
        try:
            deps_info[module] = version(module)
        except PackageNotFoundError:
            deps_info[module] = None
    return deps_info
