"""Problem-library discovery and loading for OptiProfiler.

Problem libraries expose a small, versioned contract. Bundled and legacy
path-based libraries are adapted from their ``<name>_tools.py`` modules, while
separately installed libraries are discovered through Python entry points.
Discovery is intentionally separate from loading so that listing or validating
libraries does not import optional dependencies such as PyCUTEst.
"""

from dataclasses import dataclass
from collections.abc import Mapping as MappingABC
from collections.abc import Sequence as SequenceABC
import importlib.util
from importlib import metadata
from pathlib import Path
import pickle
import re
from typing import Any, Callable, Mapping, Optional, Sequence


PROBLEM_LIBRARY_API_VERSION = 1
PROBLEM_LIBRARY_ENTRY_POINT_GROUP = 'optiprofiler.problem_libraries'

_BUILTIN_SOURCE = 'builtin'
_CUSTOM_SOURCE = 'custom'
_ENTRY_POINT_SOURCE = 'entry_point'
_PROBLEM_LIBRARY_NAME_PATTERN = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*$')


@dataclass(frozen=True)
class ProblemLibraryPlugin:
    """Versioned interface implemented by a problem-library adapter.

    Parameters
    ----------
    name : str
        Public name used in the ``plibs`` benchmark option.
    api_version : int
        Major version of the OptiProfiler problem-library protocol.
    select : callable
        ``select(problem_options, library_options)`` returns the names of
        problems matching the common and library-specific options.
    load : callable
        ``load(problem_name, library_options)`` returns an OptiProfiler
        ``Problem`` instance.
    collect_info : callable, optional
        Library-specific metadata collection helper.
    check_available : callable, optional
        Zero-argument check that returns normally when the library can be used
        and raises an informative exception otherwise. ``benchmark`` calls it
        once before selecting problems.
    get_default_options : callable, optional
        Zero-argument callable returning the library-owned raw default option
        mapping. It may read package data or environment-level defaults. The
        core overlays higher-priority process and per-run values before calling
        ``validate_options`` once, so lower-priority invalid values may be
        replaced explicitly by a benchmark run.
    validate_options : callable, optional
        Callable accepting the merged library option mapping and returning a
        validated, normalized mapping. A library with nonempty options must
        provide this callable.

    Notes
    -----
    The entry-point factory may be called independently in the main process
    and in worker processes. Consequently, ``select`` and ``load`` must not
    rely on mutable state shared through one plugin instance. OptiProfiler
    resolves library options in the main process and passes isolated copies of
    the same pickleable effective mapping to the selector and each loader.
    Callbacks should treat the mapping as read-only.
    """

    name: str
    api_version: int
    select: Callable[[Mapping[str, Any], Mapping[str, Any]], Sequence[str]]
    load: Callable[[str, Mapping[str, Any]], Any]
    collect_info: Optional[Callable[..., Any]] = None
    check_available: Optional[Callable[[], None]] = None
    get_default_options: Optional[Callable[[], Mapping[str, Any]]] = None
    validate_options: Optional[
        Callable[[Mapping[str, Any]], Mapping[str, Any]]
    ] = None


@dataclass(frozen=True)
class ProblemLibraryRef:
    """Serializable identity of a resolved problem library."""

    name: str
    source: str
    locator: str
    distribution: Optional[str] = None

    @property
    def is_builtin(self):
        """Whether this reference points to a bundled library."""
        return self.source == _BUILTIN_SOURCE


def _problem_libs_dir():
    return Path(__file__).parent.resolve() / 'problem_libs'


def _is_valid_problem_library_name(name):
    return isinstance(name, str) and bool(_PROBLEM_LIBRARY_NAME_PATTERN.fullmatch(name))


def _validate_problem_library_name(name):
    """Reject names that are unsafe or ambiguous in paths and entry points."""
    if not _is_valid_problem_library_name(name):
        raise ValueError(
            'A problem-library name must be an ASCII Python identifier: it must '
            'start with a letter or underscore and contain only letters, digits, '
            'and underscores.'
        )
    return name


def _normalize_custom_problem_libs_path(custom_libs_path):
    """Return a validated absolute path for local problem libraries."""
    custom_libs_path = Path(custom_libs_path).expanduser().resolve()
    if not custom_libs_path.exists():
        raise ValueError(
            f'Custom problem-library path does not exist: {custom_libs_path}'
        )
    if not custom_libs_path.is_dir():
        raise ValueError(
            f'Custom problem-library path must be a directory: {custom_libs_path}'
        )
    return custom_libs_path


def _custom_problem_library_names(custom_libs_path):
    """Return strict custom-library names available under a path."""
    custom_libs_path = _normalize_custom_problem_libs_path(custom_libs_path)
    if (custom_libs_path / f'{custom_libs_path.name}_tools.py').is_file():
        return [custom_libs_path.name]

    names = []
    for path in custom_libs_path.iterdir():
        if path.is_dir() and not path.name.startswith('.') and path.name != '__pycache__':
            if (path / f'{path.name}_tools.py').is_file():
                names.append(path.name)
    return sorted(names)


def _custom_problem_library_dir(custom_libs_path, name):
    """Return the unambiguous custom directory for ``name``, if present."""
    custom_libs_path = _normalize_custom_problem_libs_path(custom_libs_path)
    candidates = []
    if custom_libs_path.name == name and custom_libs_path.is_dir():
        candidates.append(custom_libs_path)
    parent_candidate = custom_libs_path / name
    if parent_candidate.is_dir():
        candidates.append(parent_candidate)

    resolved = []
    for candidate in candidates:
        candidate = candidate.resolve()
        if candidate not in resolved:
            resolved.append(candidate)
    if len(resolved) > 1:
        raise ValueError(
            f'The custom problem-library path is ambiguous for "{name}". It can '
            f'point to either the parent directory containing "{name}" or the '
            f'"{name}" directory itself, but not both.'
        )
    return resolved[0] if resolved else None


def _builtin_problem_library_names():
    names = []
    problem_dir = _problem_libs_dir()
    if not problem_dir.is_dir():
        return names
    for path in problem_dir.iterdir():
        if path.is_dir() and not path.name.startswith('.') and path.name != '__pycache__':
            if (path / f'{path.name}_tools.py').is_file():
                names.append(path.name)
    return sorted(names)


def _problem_library_entry_points():
    """Return entry points for all installed OptiProfiler problem libraries."""
    entry_points = metadata.entry_points()
    if hasattr(entry_points, 'select'):
        return list(entry_points.select(group=PROBLEM_LIBRARY_ENTRY_POINT_GROUP))
    return list(entry_points.get(PROBLEM_LIBRARY_ENTRY_POINT_GROUP, ()))


def _entry_point_distribution(entry_point):
    distribution = getattr(entry_point, 'dist', None)
    if distribution is None:
        return None
    name = getattr(distribution, 'name', None)
    if name:
        return str(name)
    try:
        return str(distribution.metadata['Name'])
    except (AttributeError, KeyError, TypeError):
        return None


def list_problem_libraries(custom_problem_libs_path=None):
    """List discoverable libraries without importing optional dependencies.

    Parameters
    ----------
    custom_problem_libs_path : str or pathlib.Path, optional
        Parent directory containing local libraries, or one local library
        directory. Local adapters use the strict ``<name>_tools.py`` layout.

    Returns
    -------
    list of str
        Sorted public library names visible to OptiProfiler.

    Notes
    -----
    Discovery reads package metadata and filenames only. A listed library may
    still be unusable when selected if an optional runtime is missing. Provider
    conflicts are checked by ``benchmark`` when a name is resolved.
    """
    names = {
        name for name in _builtin_problem_library_names()
        if _is_valid_problem_library_name(name)
    }
    names.update(
        entry_point.name for entry_point in _problem_library_entry_points()
        if _is_valid_problem_library_name(entry_point.name)
    )
    if custom_problem_libs_path is not None:
        names.update(
            name for name in _custom_problem_library_names(custom_problem_libs_path)
            if _is_valid_problem_library_name(name)
        )
    return sorted(names)


def _normalize_selected_problem_names(problem_names, library_name):
    """Validate and normalize a library selector's ordered result."""
    if isinstance(problem_names, (str, bytes)) or not isinstance(
        problem_names, SequenceABC
    ):
        raise TypeError(
            f'Problem library "{library_name}" select(options) must return an '
            'ordered sequence of problem-name strings, not a single string.'
        )
    problem_names = list(problem_names)
    if any(not isinstance(name, str) for name in problem_names):
        raise TypeError(
            f'Problem library "{library_name}" select(options) returned a '
            'non-string problem name.'
        )
    return problem_names


def _normalize_problem_library_options(options, library_name, source):
    """Materialize and validate one problem library's option mapping."""
    if not isinstance(options, MappingABC):
        raise TypeError(
            f'Problem library "{library_name}" {source} must be a mapping.'
        )
    options = dict(options)
    if any(not isinstance(key, str) for key in options):
        raise TypeError(
            f'Problem library "{library_name}" {source} must use string keys.'
        )
    return options


def _copy_problem_library_options(options):
    """Return an isolated, type-preserving copy of resolved options.

    Resolved options are required to be pickleable. A pickle round trip makes
    sequential callbacks observe the same process-isolation semantics as
    callbacks executed by spawned workers, including for nested mutable values.
    """
    return pickle.loads(pickle.dumps(options))


def _resolve_problem_library_options(plugin, overrides=None):
    """Resolve one plugin's defaults and per-run overrides deterministically."""
    if plugin.get_default_options is None:
        defaults = {}
    else:
        defaults = _normalize_problem_library_options(
            plugin.get_default_options(),
            plugin.name,
            'get_default_options() result',
        )

    if overrides is None:
        overrides = {}
    overrides = _normalize_problem_library_options(
        overrides,
        plugin.name,
        'option overrides',
    )
    try:
        # Defaults and user overrides may contain nested mutable values. Make
        # independent copies before merging so a validator cannot mutate the
        # caller's original ``plib_options`` or package-owned defaults.
        defaults = _copy_problem_library_options(defaults)
        overrides = _copy_problem_library_options(overrides)
    except Exception as exc:
        raise TypeError(
            f'Problem library "{plugin.name}" options must be pickleable so that '
            'they can be stored and passed to worker processes.'
        ) from exc
    merged = defaults.copy()
    merged.update(overrides)

    if plugin.validate_options is None:
        if merged:
            raise ValueError(
                f'Problem library "{plugin.name}" does not declare configurable '
                'options, but nonempty options were provided or returned as defaults.'
            )
        normalized = {}
    else:
        normalized = _normalize_problem_library_options(
            plugin.validate_options(_copy_problem_library_options(merged)),
            plugin.name,
            'validate_options() result',
        )

    try:
        pickle.dumps(normalized)
    except Exception as exc:
        raise TypeError(
            f'Problem library "{plugin.name}" options must be pickleable so that '
            'they can be stored and passed to worker processes.'
        ) from exc
    return normalized


def resolve_problem_library(name, custom_problem_libs_path=None):
    """Resolve a library name to a serializable reference without importing it.

    An explicitly supplied local library takes precedence. Otherwise, a name
    must have exactly one provider: either a bundled adapter or one installed
    entry point. Ambiguous installed providers are rejected rather than chosen
    according to environment-dependent entry-point ordering.
    """
    name = _validate_problem_library_name(name)

    if custom_problem_libs_path is not None:
        custom_dir = _custom_problem_library_dir(custom_problem_libs_path, name)
        if custom_dir is not None:
            tools_file = custom_dir / f'{name}_tools.py'
            if not tools_file.is_file():
                raise ValueError(
                    f'The custom problem library "{name}" at {custom_dir} must '
                    f'contain "{name}_tools.py".'
                )
            return ProblemLibraryRef(name, _CUSTOM_SOURCE, str(tools_file.resolve()))

    builtin_tools_file = _problem_libs_dir() / name / f'{name}_tools.py'
    entry_points = [
        entry_point for entry_point in _problem_library_entry_points()
        if entry_point.name == name
    ]

    if builtin_tools_file.is_file() and entry_points:
        distributions = sorted({
            _entry_point_distribution(entry_point) or entry_point.value
            for entry_point in entry_points
        })
        raise ValueError(
            f'Problem library "{name}" is ambiguous because it is both bundled '
            f'and installed through {distributions}. Remove the duplicate provider.'
        )
    if len(entry_points) > 1:
        distributions = sorted({
            _entry_point_distribution(entry_point) or entry_point.value
            for entry_point in entry_points
        })
        raise ValueError(
            f'Problem library "{name}" is provided by multiple installed '
            f'distributions: {distributions}. Keep exactly one provider.'
        )
    if builtin_tools_file.is_file():
        return ProblemLibraryRef(
            name,
            _BUILTIN_SOURCE,
            str(builtin_tools_file.resolve()),
            distribution='optiprofiler',
        )
    if entry_points:
        entry_point = entry_points[0]
        return ProblemLibraryRef(
            name,
            _ENTRY_POINT_SOURCE,
            entry_point.value,
            distribution=_entry_point_distribution(entry_point),
        )

    available = list_problem_libraries(custom_problem_libs_path)
    raise ValueError(
        f'Problem library "{name}" is not available. Available libraries: {available}.'
    )


def _validate_plugin(plugin, expected_name):
    if not isinstance(plugin, ProblemLibraryPlugin):
        raise TypeError(
            f'Problem library "{expected_name}" must return a '
            'ProblemLibraryPlugin from its entry-point factory.'
        )
    if plugin.name != expected_name:
        raise ValueError(
            f'Problem-library entry point "{expected_name}" returned a plugin '
            f'named "{plugin.name}".'
        )
    if type(plugin.api_version) is not int or plugin.api_version != PROBLEM_LIBRARY_API_VERSION:
        raise ValueError(
            f'Problem library "{expected_name}" uses API version '
            f'{plugin.api_version}; OptiProfiler requires '
            f'{PROBLEM_LIBRARY_API_VERSION}.'
        )
    if not callable(plugin.select) or not callable(plugin.load):
        raise TypeError(
            f'Problem library "{expected_name}" must provide callable select '
            'and load functions.'
        )
    if plugin.collect_info is not None and not callable(plugin.collect_info):
        raise TypeError(
            f'Problem library "{expected_name}" has a non-callable collect_info.'
        )
    if plugin.check_available is not None and not callable(plugin.check_available):
        raise TypeError(
            f'Problem library "{expected_name}" has a non-callable check_available.'
        )
    if plugin.get_default_options is not None and not callable(plugin.get_default_options):
        raise TypeError(
            f'Problem library "{expected_name}" has a non-callable get_default_options.'
        )
    if plugin.validate_options is not None and not callable(plugin.validate_options):
        raise TypeError(
            f'Problem library "{expected_name}" has a non-callable validate_options.'
        )
    if (plugin.get_default_options is None) != (plugin.validate_options is None):
        raise TypeError(
            f'Problem library "{expected_name}" must provide both '
            'get_default_options and validate_options, or neither.'
        )
    return plugin


def _load_tools_module(reference):
    module_path = Path(reference.locator)
    if reference.is_builtin:
        module_name = (
            f'optiprofiler.problem_libs.{reference.name}.'
            f'{reference.name}_tools'
        )
    else:
        module_name = f'{reference.name}.{reference.name}_tools'

    spec = importlib.util.spec_from_file_location(module_name, str(module_path))
    if spec is None or spec.loader is None:
        raise ImportError(
            f'Cannot create an import specification for problem library '
            f'"{reference.name}" at {module_path}.'
        )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_legacy_plugin(reference):
    module = _load_tools_module(reference)
    legacy_select = getattr(module, f'{reference.name}_select')
    legacy_load = getattr(module, f'{reference.name}_load')
    collect_info = getattr(module, f'{reference.name}_collect_info', None)
    if collect_info is None:
        collect_info = getattr(module, f'{reference.name}_get_info', None)
    check_available = getattr(
        module,
        f'{reference.name}_check_available',
        None,
    )
    get_default_options = getattr(
        module,
        f'{reference.name}_get_default_options',
        None,
    )
    validate_options = getattr(
        module,
        f'{reference.name}_validate_options',
        None,
    )
    if (get_default_options is None) != (validate_options is None):
        raise TypeError(
            f'Legacy problem library "{reference.name}" must define both '
            f'"{reference.name}_get_default_options" and '
            f'"{reference.name}_validate_options", or neither.'
        )

    if get_default_options is None:
        def select(problem_options, library_options):
            if library_options:
                raise ValueError(
                    f'Legacy problem library "{reference.name}" does not support '
                    'explicit library options.'
                )
            return legacy_select(problem_options)

        def load(problem_name, library_options):
            if library_options:
                raise ValueError(
                    f'Legacy problem library "{reference.name}" does not support '
                    'explicit library options.'
                )
            return legacy_load(problem_name)
    else:
        def select(problem_options, library_options):
            return legacy_select(problem_options, library_options)

        def load(problem_name, library_options):
            return legacy_load(problem_name, library_options=library_options)

    return ProblemLibraryPlugin(
        name=reference.name,
        api_version=PROBLEM_LIBRARY_API_VERSION,
        select=select,
        load=load,
        collect_info=collect_info,
        check_available=check_available,
        get_default_options=get_default_options,
        validate_options=validate_options,
    )


def _entry_point_from_reference(reference):
    """Rebuild an already resolved entry point without rescanning packages."""
    return metadata.EntryPoint(
        reference.name,
        reference.locator,
        PROBLEM_LIBRARY_ENTRY_POINT_GROUP,
    )


def load_problem_library(reference):
    """Load and validate a previously resolved problem-library reference."""
    if not isinstance(reference, ProblemLibraryRef):
        raise TypeError('reference must be a ProblemLibraryRef.')
    if reference.source in {_BUILTIN_SOURCE, _CUSTOM_SOURCE}:
        plugin = _load_legacy_plugin(reference)
    elif reference.source == _ENTRY_POINT_SOURCE:
        factory = _entry_point_from_reference(reference).load()
        if not callable(factory):
            raise TypeError(
                f'Entry point for problem library "{reference.name}" must '
                'resolve to a zero-argument factory.'
            )
        plugin = factory()
    else:
        raise ValueError(
            f'Unknown problem-library source "{reference.source}" for '
            f'"{reference.name}".'
        )
    return _validate_plugin(plugin, reference.name)


__all__ = [
    'PROBLEM_LIBRARY_API_VERSION',
    'PROBLEM_LIBRARY_ENTRY_POINT_GROUP',
    'ProblemLibraryPlugin',
    'ProblemLibraryRef',
    'list_problem_libraries',
    'load_problem_library',
    'resolve_problem_library',
]
