"""Problem-library configuration compatibility and per-run resolution."""

import os
from pathlib import Path

from .problem_libraries import (
    _copy_problem_library_options,
    _resolve_problem_library_options,
    load_problem_library,
    resolve_problem_library,
)


_PLIB_CONFIG_OVERRIDES = {}


def _reference_key(reference):
    """Return a provider-specific key for process-local overrides."""
    return (
        reference.name,
        reference.source,
        reference.locator,
        reference.distribution,
    )


def _config_path_from_reference(reference):
    """Return a legacy adapter's package-owned config file, if visible."""
    if reference.source not in {'builtin', 'custom'}:
        return None
    config_path = Path(reference.locator).parent / 'config.txt'
    return config_path if config_path.is_file() else None


def _get_plib_config_path(plib, custom_problem_libs_path=None):
    """Locate a visible legacy ``config.txt`` file for compatibility."""
    reference = resolve_problem_library(plib, custom_problem_libs_path)
    return _config_path_from_reference(reference)


def _parse_config_file(config_path):
    """Parse a legacy ``config.txt`` file into a key-value dictionary."""
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


def _legacy_file_config(reference):
    """Read the effective config of an old file-based adapter."""
    config_path = _config_path_from_reference(reference)
    if config_path is None:
        raise FileNotFoundError(
            f'Problem library "{reference.name}" does not expose configurable options.'
        )
    config = _parse_config_file(config_path)
    for key in config:
        env_key = f'{reference.name.upper()}_{key.upper()}'
        if env_key in os.environ:
            value = os.environ[env_key]
            try:
                value = int(value)
            except ValueError:
                pass
            config[key] = value
    return config


def _load_library_for_config(plib, custom_problem_libs_path=None):
    reference = resolve_problem_library(plib, custom_problem_libs_path)
    return reference, load_problem_library(reference)


def _resolve_plib_options(
    plib,
    run_options=None,
    custom_problem_libs_path=None,
):
    """Resolve one library's effective options for a benchmark call."""
    reference, plugin = _load_library_for_config(
        plib,
        custom_problem_libs_path,
    )
    session_options = _PLIB_CONFIG_OVERRIDES.get(_reference_key(reference), {})
    run_options = {} if run_options is None else dict(run_options)

    if plugin.get_default_options is None and plugin.validate_options is None:
        if run_options:
            raise ValueError(
                f'Problem library "{plib}" does not support explicit '
                '`plib_options`. Its legacy selector remains usable with its '
                'existing package-level configuration.'
            )
        # Old adapters continue to read their own files/environment. Passing
        # inferred values would falsely imply that the core validated them.
        return {}

    overrides = dict(session_options)
    overrides.update(run_options)
    return _resolve_problem_library_options(plugin, overrides)


def get_plib_config(
    plib,
    verbose=False,
    custom_problem_libs_path=None,
):
    """Read the effective process-level configuration of a problem library.

    Library-owned defaults are combined with overrides previously set by
    :func:`set_plib_config`. Per-run ``benchmark(..., plib_options=...)`` values
    are intentionally not included because they apply only to that benchmark
    call.

    Parameters
    ----------
    plib : str
        Public problem-library name.
    verbose : bool, optional
        Print the package-owned ``config.txt`` when one is visible. For an
        installed plugin without such a file, print the effective mapping.
    custom_problem_libs_path : str or pathlib.Path, optional
        Explicit local provider path, using the same resolution rules as
        :func:`optiprofiler.benchmark`.

    Returns
    -------
    dict
        Effective process-level configuration. Modern plugin values are
        validated by the provider; legacy built-in file values retain their
        historical parsing behavior and are validated when the adapter runs.

    Examples
    --------
    .. code-block:: python

        from optiprofiler import get_plib_config

        print(get_plib_config('s2mpj'))
        # {'variable_size': 'default', 'test_feasibility_problems': 0}
    """
    reference = resolve_problem_library(plib, custom_problem_libs_path)
    config_path = _config_path_from_reference(reference)

    # Preserve the original built-in configuration API without importing the
    # adapter. Reading a package-owned config file should remain import-safe
    # even when a bundled adapter has optional runtime dependencies.
    if reference.source == 'builtin' and config_path is not None:
        config = _legacy_file_config(reference)
        config.update(_copy_problem_library_options(
            _PLIB_CONFIG_OVERRIDES.get(_reference_key(reference), {})
        ))
    else:
        plugin = load_problem_library(reference)
        if plugin.get_default_options is None and plugin.validate_options is None:
            if reference.source == 'entry_point':
                # A modern plugin can explicitly declare that it has no
                # library-specific options by omitting both callbacks.
                config = {}
            else:
                config = _legacy_file_config(reference)
        else:
            overrides = _PLIB_CONFIG_OVERRIDES.get(_reference_key(reference), {})
            config = _resolve_problem_library_options(plugin, overrides)

    if verbose:
        if config_path is not None:
            with open(config_path, 'r') as f:
                print(f.read())
        else:
            print(config)
    return config


def set_plib_config(
    plib,
    *,
    custom_problem_libs_path=None,
    **kwargs,
):
    """Set process-level default overrides for one problem library.

    This backward-compatible API affects subsequent calls in the current
    Python process and scopes each override to the resolved provider. For a
    reproducible experiment, prefer the higher-priority ``plib_options``
    argument of :func:`optiprofiler.benchmark`; its resolved values are stored
    with the experiment.

    Parameters
    ----------
    plib : str
        Public problem-library name.
    custom_problem_libs_path : str or pathlib.Path, optional
        Explicit local provider path, using the same resolution rules as
        :func:`optiprofiler.benchmark`.
    **kwargs
        Library-specific option overrides.

    Examples
    --------
    .. code-block:: python

        from optiprofiler import set_plib_config

        set_plib_config(
            's2mpj',
            variable_size='all',
            test_feasibility_problems=2,
        )
    """
    reference = resolve_problem_library(plib, custom_problem_libs_path)
    reference_key = _reference_key(reference)
    config_path = _config_path_from_reference(reference)

    # Keep package-owned built-in configuration usable without importing an
    # optional runtime. Values retain the historical key-only validation here;
    # a modern adapter validates them when benchmark actually loads the library.
    if reference.source == 'builtin' and config_path is not None:
        config = _legacy_file_config(reference)
        unknown = sorted(set(kwargs) - set(config))
        if unknown:
            raise ValueError(
                f'Unknown config variables for problem library "{plib}": '
                f'{unknown}. Available variables: {sorted(config)}'
            )
        overrides = dict(_PLIB_CONFIG_OVERRIDES.get(reference_key, {}))
        overrides.update(kwargs)
        _PLIB_CONFIG_OVERRIDES[reference_key] = (
            _copy_problem_library_options(overrides)
        )
        config.update(overrides)
    else:
        plugin = load_problem_library(reference)
        if plugin.get_default_options is None and plugin.validate_options is None:
            if reference.source == 'entry_point':
                if kwargs:
                    raise ValueError(
                        f'Problem library "{plib}" does not declare configurable '
                        'options.'
                    )
                return {}
            if reference.source != 'builtin':
                raise ValueError(
                    f'Problem library "{plib}" uses a legacy local adapter without '
                    'the configuration callbacks required for provider-scoped '
                    'process overrides. Define both get_default_options and '
                    'validate_options, or configure the adapter outside OptiProfiler.'
                )
            config = _legacy_file_config(reference)
            unknown = sorted(set(kwargs) - set(config))
            if unknown:
                raise ValueError(
                    f'Unknown config variables for problem library "{plib}": '
                    f'{unknown}. Available variables: {sorted(config)}'
                )
            config.update(kwargs)
        else:
            overrides = dict(_PLIB_CONFIG_OVERRIDES.get(reference_key, {}))
            overrides.update(kwargs)
            effective = _resolve_problem_library_options(plugin, overrides)
            _PLIB_CONFIG_OVERRIDES[reference_key] = (
                _copy_problem_library_options(overrides)
            )
            config = effective

    # Preserve the old built-in adapter behavior for direct selector calls.
    # Custom and installed providers stay isolated in the provider-scoped
    # override table; mirroring them into a name-only environment variable
    # could accidentally alter a different provider with the same public name.
    legacy_keys = set(_parse_config_file(config_path)) if config_path else set()
    for key in kwargs:
        if reference.source == 'builtin' and key in legacy_keys:
            os.environ[f'{plib.upper()}_{key.upper()}'] = str(config[key])
