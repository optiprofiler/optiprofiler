import os
from pathlib import Path


def _get_plib_config_path(plib):
    """
    Locate the config.txt file for a problem library.

    Parameters
    ----------
    plib : str
        Name of the problem library.

    Returns
    -------
    Path or None
        Path to the config.txt file, or None if it does not exist.
    """
    builtin_dir = Path(__file__).parent.resolve() / 'problem_libs' / plib
    config_path = builtin_dir / 'config.txt'
    return config_path if config_path.exists() else None


def _parse_config_file(config_path):
    """
    Parse a config.txt file and return a dict of key-value pairs.

    Lines starting with ``#`` are treated as comments and blank lines are
    ignored.  Values that look like integers are converted to ``int``.

    Parameters
    ----------
    config_path : Path
        Path to the config.txt file.

    Returns
    -------
    dict
        Parsed configuration as ``{variable_name: value}``.
    """
    config = {}
    with open(config_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            if '=' in stripped:
                key, _, value = stripped.partition('=')
                key = key.strip()
                value = value.split('#')[0].split('%')[0].strip()
                try:
                    value = int(value)
                except ValueError:
                    pass
                config[key] = value
    return config


def get_plib_config(plib, verbose=False):
    """
    Read the current configuration of a problem library.

    The returned dictionary reflects the effective values: if an
    environment variable ``<PLIB>_<VARIABLE>`` (all upper-case) has been
    set (e.g.  via `set_plib_config`), it takes precedence over the
    value in the library's ``config.txt``.

    Parameters
    ----------
    plib : str
        Name of the problem library (e.g. ``'s2mpj'``, ``'pycutest'``).
    verbose : bool, optional
        If ``True``, the full contents of ``config.txt`` (including
        comments) are printed so that the user can see all available
        options and their descriptions. Default is ``False``.

    Returns
    -------
    dict
        Effective configuration as ``{variable_name: value}``.

    Raises
    ------
    FileNotFoundError
        If the problem library does not have a ``config.txt`` file.

    Examples
    --------
    .. code-block:: python

        from optiprofiler import get_plib_config

        config = get_plib_config('s2mpj')
        print(config)
        # {'variable_size': 'default', 'test_feasibility_problems': 0}

        get_plib_config('s2mpj', verbose=True)
        # Prints the full config.txt with comments, then returns the dict.
    """
    config_path = _get_plib_config_path(plib)
    if config_path is None:
        raise FileNotFoundError(
            f"No config.txt found for problem library '{plib}'."
        )

    if verbose:
        with open(config_path, 'r') as f:
            print(f.read())

    config = _parse_config_file(config_path)

    for key in config:
        env_key = f'{plib.upper()}_{key.upper()}'
        if env_key in os.environ:
            env_val = os.environ[env_key]
            try:
                env_val = int(env_val)
            except ValueError:
                pass
            config[key] = env_val

    return config


def set_plib_config(plib, **kwargs):
    """
    Override configuration variables for a problem library.

    Each keyword argument is translated to the environment variable
    ``<PLIB>_<VARIABLE>`` (all upper-case) so that subsequent calls to
    the library's ``select`` function (directly or through `benchmark`)
    will pick up the new value.  The override persists for the lifetime
    of the current Python process.

    Parameters
    ----------
    plib : str
        Name of the problem library (e.g. ``'s2mpj'``, ``'pycutest'``).
    **kwargs
        Configuration variables to set.  The variable names must exist
        in the library's ``config.txt``; otherwise a ``ValueError`` is
        raised to prevent silent typos.

    Raises
    ------
    FileNotFoundError
        If the problem library does not have a ``config.txt`` file.
    ValueError
        If a variable name is not recognised.

    Examples
    --------
    .. code-block:: python

        from optiprofiler import set_plib_config, get_plib_config

        set_plib_config('s2mpj', variable_size='min')
        print(get_plib_config('s2mpj'))
        # {'variable_size': 'min', 'test_feasibility_problems': 0}

        set_plib_config('s2mpj',
                        variable_size='all',
                        test_feasibility_problems=2)
    """
    config_path = _get_plib_config_path(plib)
    if config_path is None:
        raise FileNotFoundError(
            f"No config.txt found for problem library '{plib}'."
        )

    known_keys = set(_parse_config_file(config_path).keys())

    for key, value in kwargs.items():
        if key not in known_keys:
            raise ValueError(
                f"Unknown config variable '{key}' for problem library "
                f"'{plib}'. Available variables: {sorted(known_keys)}"
            )
        env_key = f'{plib.upper()}_{key.upper()}'
        os.environ[env_key] = str(value)
