"""
Trusted compatibility boundary for configurations serialized by earlier versions.

Result archives, ``options_*.pkl`` files and configuration pickles written by
OptiProfiler 1.x reference class paths and enumeration members that no longer
exist in the 2.0 core: ``optiprofiler.utils.FeatureOption`` had an ``n_runs``
member (the run count now belongs to the experiment), the 1.x
``optiprofiler.opclasses.Feature`` stored mixed stage/experiment options, and
the 1.x bridge stored ``optiprofiler.composition.ComposedFeature`` and
``Stage`` objects. Those layouts are decoded here, at the boundary, into small
compatibility values before anything is normalized; fresh core objects never
restore historical state.

* Historical enumeration values decode to the active member when it still
  exists and to a string-valued :class:`LegacyEnumValue` otherwise, so old
  dictionaries keep their keys and values without reintroducing the member.
* Historical Feature layouts decode to :class:`LegacyFeature`,
  :class:`LegacyComposedFeature` and :class:`LegacyStage` containers that hold
  the raw state. :func:`import_legacy_feature` converts a container into a
  canonical :class:`~optiprofiler.Feature` (stage-local options only) plus the
  retained run count, which is passed to ``benchmark(..., n_runs=...)``.
* :func:`load_options` reads ``options_user.pkl``/``options_refined.pkl`` and
  :func:`replay_arguments` maps the supported layouts to ``benchmark`` inputs.

Trust statement: these readers unpickle. Unpickling executes the reduction
hooks of the objects in the file, so they are for archives and configuration
files produced by OptiProfiler on machines you trust. They offer no safety for
untrusted pickle payloads, and they do not resume live ``FeaturedProblem``
objects across versions.
"""

import copy
import io
import pickle
from collections import namedtuple
from collections.abc import Mapping
from enum import Enum
from pathlib import Path

from . import utils as _utils
from .feature_definitions import EXPERIMENT_OPTIONS, STAGE_NAMES, validated_local_options
from .utils import FeatureOption, ProblemOption

REFINED_SCHEMA = 'options_refined-v2'
USER_SCHEMA = 'options_user-v2'


class LegacyConfigurationError(ValueError):
    """A historical configuration cannot be represented as a 2.0 specification."""


def historical_effective_options(name, options):
    """Preserve the grid actually used by old native/refined configurations."""
    options = dict(options)
    mesh = options.get('mesh_type')
    if (name == 'quantized' and isinstance(mesh, str)
            and mesh.lower() in ('absolute', 'relative') and mesh != mesh.lower()):
        # Before the validator fix, any accepted non-lowercase spelling took
        # the absolute runtime branch. Reinterpreting it as fresh input would
        # silently change a saved experiment. Only effective options change;
        # the historical declaration/provenance remains untouched.
        options['mesh_type'] = 'absolute'
    return options


def _historical_specification(specification):
    """Copy saved structured input without modifying the caller's archive data."""
    if isinstance(specification, Mapping):
        entry = dict(specification)
        name, options = entry.get('name'), entry.get('options', {})
        if isinstance(name, str) and isinstance(options, Mapping) and 'options' in entry:
            entry['options'] = historical_effective_options(name.strip().lower(), options)
        return entry
    if isinstance(specification, (list, tuple)):
        entries = [_historical_specification(entry) for entry in specification]
        return tuple(entries) if isinstance(specification, tuple) else entries
    return specification


class LegacyEnumValue(str):
    """
    A historical enumeration value whose member no longer exists (for example
    ``FeatureOption('n_runs')``). It compares and hashes as its string value,
    so dictionaries keyed by it are looked up by the string spelling, and it
    keeps the enumeration name and the value for inspection.
    """

    def __new__(cls, value, enum_name='FeatureOption'):
        self = super().__new__(cls, value)
        self.enum_name = str(enum_name)
        return self

    def __getnewargs__(self):
        return (str(self), self.enum_name)

    @property
    def value(self):
        return str(self)

    @property
    def name(self):
        return str(self).upper()

    def __repr__(self):
        return f'<LegacyEnumValue {self.enum_name}.{self.name}: {str(self)!r}>'


def _tolerant_enum(enum_class):
    """A factory standing in for ``enum_class`` in a pickle: member if it exists, marker otherwise."""

    def decode(value):
        try:
            return enum_class(value)
        except ValueError:
            return LegacyEnumValue(value, enum_class.__name__)

    decode.__name__ = f'legacy_{enum_class.__name__}'
    return decode


class LegacyObject:
    """A historical object decoded by the trusted unpickler: its class path and raw state."""

    class_path = 'optiprofiler.legacy_compat.LegacyObject'

    def __init__(self, state=None):
        self.state = dict(state or {})

    def __setstate__(self, state):
        if isinstance(state, tuple) and len(state) == 2:
            # Objects with __slots__ pickle as (dict state, slot state).
            merged = dict(state[0] or {})
            merged.update(state[1] or {})
        elif isinstance(state, dict):
            merged = dict(state)
        else:
            raise LegacyConfigurationError(f'Unsupported serialized state for {self.class_path}: {type(state).__name__}.')
        self.state = merged

    def __getstate__(self):
        return dict(getattr(self, 'state', {}))

    def __repr__(self):
        return f'{type(self).__name__}({self.class_path!r}, {sorted(map(str, getattr(self, "state", {})))})'


class LegacyFeature(LegacyObject):
    """The 1.x ``optiprofiler.opclasses.Feature`` (``_name`` plus mixed ``_options``)."""

    class_path = 'optiprofiler.opclasses.Feature'


class LegacyComposedFeature(LegacyObject):
    """The 1.x bridge ``optiprofiler.composition.ComposedFeature`` (root ``_options`` plus ``_stages``)."""

    class_path = 'optiprofiler.composition.ComposedFeature'


class LegacyStage(LegacyObject):
    """The 1.x bridge ``optiprofiler.composition.Stage`` (``position``, ``name``, ``code``, ``occurrence``, ``feature``)."""

    class_path = 'optiprofiler.composition.Stage'


_LEGACY_CLASSES = {
    ('optiprofiler.opclasses', 'Feature'): LegacyFeature,
    ('optiprofiler.composition', 'ComposedFeature'): LegacyComposedFeature,
    ('optiprofiler.composition', 'Stage'): LegacyStage,
}


class TrustedUnpickler(pickle.Unpickler):
    """Unpickler mapping historical OptiProfiler class paths to compatibility values."""

    def find_class(self, module, name):
        key = (module, name)
        if key in _LEGACY_CLASSES:
            return _LEGACY_CLASSES[key]
        if module == 'optiprofiler.utils':
            target = getattr(_utils, name, None)
            if isinstance(target, type) and issubclass(target, Enum):
                return _tolerant_enum(target)
        return super().find_class(module, name)


def loads_trusted(data):
    """Unpickle trusted bytes written by any OptiProfiler version (see the module trust statement)."""
    return TrustedUnpickler(io.BytesIO(bytes(data))).load()


def load_trusted(path):
    """Unpickle a trusted file written by any OptiProfiler version."""
    with open(Path(path).expanduser(), 'rb') as stream:
        return TrustedUnpickler(stream).load()


def load_options(path):
    """
    Read a trusted ``options_user.pkl`` or ``options_refined.pkl`` file as a
    dictionary. A file that is not a readable OptiProfiler pickle (truncated,
    empty, foreign bytes, missing classes) raises
    :class:`LegacyConfigurationError` naming the underlying error; a missing
    file raises the usual ``OSError``.
    """
    try:
        options = load_trusted(path)
    except LegacyConfigurationError:
        raise
    except (pickle.UnpicklingError, EOFError, AttributeError, ImportError, IndexError, TypeError, ValueError) as err:
        raise LegacyConfigurationError(f'{path} is not a readable OptiProfiler options pickle '
                                       f'({type(err).__name__}: {err}).') from err
    if not isinstance(options, dict):
        raise LegacyConfigurationError(f'{path} does not contain an options dictionary ({type(options).__name__}).')
    return options


def _option_key(key):
    """The string spelling of an option key (active member, legacy marker or plain string)."""
    value = getattr(key, 'value', None)
    if isinstance(key, str) and isinstance(value, str):
        return value
    return key


def _split_count(options):
    """Split a historical option mapping into stage-local options and the retained run count."""
    local, count = {}, None
    for key, value in dict(options).items():
        key = _option_key(key)
        if not isinstance(key, str):
            raise LegacyConfigurationError(f'Historical option names must be strings, got {type(key).__name__}.')
        if key in EXPERIMENT_OPTIONS:
            count = value
        else:
            local[key] = value
    return local, count


LegacyFeatureImport = namedtuple('LegacyFeatureImport', 'feature n_runs declared_name declared_spec source')
LegacyFeatureImport.__doc__ = """
The result of :func:`import_legacy_feature`.

``feature`` is the canonical specification built from the effective stages and
their stage-local options; its declaration is the one the historical object
recorded (1.x bridge objects) and otherwise unknown (``feature.declared.route``
is ``None``, ``feature.declared_name`` is ``None``); ``n_runs`` is the run count
retained by the historical object (a resolved value, not evidence of an
explicit request; ``None`` when none was stored); ``declared_name`` and
``declared_spec`` are the historical declaration verbatim when the object
recorded one and ``None`` otherwise; ``source`` is the historical class path.
"""


def import_legacy_feature(container):
    """
    Convert a decoded historical Feature container into a canonical ``Feature``
    and the run count it retained.

    Stage-local options are validated by the 2.0 definitions; the root run count
    of a composite wins over any residue stored in its child stages; a recorded
    declaration is restored, and nothing is guessed for a declaration the
    object did not record. Raises
    :class:`LegacyConfigurationError` when the historical options are not a valid
    2.0 specification and ``TypeError`` for anything that is not a legacy
    container.
    """
    if not isinstance(container, LegacyObject):
        raise TypeError('import_legacy_feature expects a legacy Feature container decoded by loads_trusted or '
                        f'load_trusted, not {type(container).__name__}.')
    state = getattr(container, 'state', None) or {}
    if isinstance(container, LegacyComposedFeature):
        entries = []
        for stage in state.get('_stages') or ():
            if not isinstance(stage, LegacyStage):
                raise LegacyConfigurationError(f'Unsupported stage object in {container.class_path}: {type(stage).__name__}.')
            child = stage.state.get('feature')
            child_state = getattr(child, 'state', None) or {} if isinstance(child, LegacyObject) else {}
            local, _ = _split_count(child_state.get('_options') or {})
            name = stage.state.get('name', child_state.get('_name'))
            entries.append({'name': name, 'options': local})
        _, count = _split_count(state.get('_options') or {})
        if not entries:
            entries = [{'name': 'plain', 'options': {}}]
    elif isinstance(container, LegacyFeature):
        if '_name' not in state or '_options' not in state:
            raise LegacyConfigurationError(f'Unsupported {container.class_path} layout: expected `_name` and `_options`, '
                                           f'found {sorted(map(str, state))}.')
        local, count = _split_count(state['_options'])
        name = state['_name']
        if not isinstance(name, str):
            raise LegacyConfigurationError(f'{container.class_path}: `_name` must be a string, got {type(name).__name__}.')
        name = name.strip().lower()
        if name == 'plain':
            # The 1.x identity object: no effective stage. Its local options must
            # still be valid for `plain` (none are accepted), never ignored.
            _validate_recorded_options('plain', local, container.class_path)
            entries = []
        else:
            entries = [{'name': name, 'options': local}]
    else:
        raise TypeError(f'{type(container).__name__} ({container.class_path}) is not a historical Feature container.')
    # The imported Feature is built in its native form: the effective stages
    # are validated by the 2.0 definitions. The declaration is restored only
    # when the historical object recorded one (the 1.x bridge stored
    # `_declared_spec`, and composites their `_route`); objects that recorded
    # none get an unknown declaration (route None, no entries). Nothing is
    # fabricated from the effective stages.
    from .opclasses import FEATURE_NATIVE_VERSION, _rebuild_feature
    declared_spec = state.get('_declared_spec')
    declared_spec = copy.deepcopy(declared_spec) if declared_spec is not None else None
    route, declared_entries = None, ()
    if declared_spec is not None:
        route = state.get('_route', 'feature_name')
        declared_entries = _recorded_declaration(declared_spec, entries, container.class_path)
    try:
        feature = _rebuild_feature(FEATURE_NATIVE_VERSION, route, declared_entries,
                                   tuple((entry['name'], entry['options']) for entry in entries))
    except (TypeError, ValueError) as err:
        raise LegacyConfigurationError(f'The historical feature configuration ({container.class_path}) is not a valid '
                                       f'2.0 specification: {err}') from err
    return LegacyFeatureImport(feature, count, state.get('_declared_name'), declared_spec, container.class_path)


def _validate_recorded_options(name, options, class_path):
    """Validate recorded stage options with the 2.0 definitions of that stage (``plain`` accepts none)."""
    if name not in STAGE_NAMES:
        raise LegacyConfigurationError(f'{class_path}: unknown recorded stage {name!r}.')
    try:
        validated_local_options(name, options)
    except (TypeError, ValueError) as err:
        raise LegacyConfigurationError(f'{class_path}: recorded options of stage {name!r} are invalid: {err}') from err


def _recorded_declaration(declared_spec, entries, class_path):
    """
    The declaration entries recorded by a 1.x bridge object, checked against
    the effective stages it carried: a list of ``{'name', 'options'}`` mappings,
    ``plain`` entries included and kept, every entry's options valid for its
    stage, and the non-plain names, in order, equal to the effective stage
    names. A recorded declaration that does not describe the object is an
    error, never repaired or ignored.
    """
    if not isinstance(declared_spec, (list, tuple)):
        raise LegacyConfigurationError(f'{class_path}: the recorded declaration is not a list of stage entries.')
    recorded = []
    for entry in declared_spec:
        if not isinstance(entry, Mapping) or not isinstance(entry.get('name'), str):
            raise LegacyConfigurationError(f'{class_path}: malformed recorded declaration entry {entry!r}.')
        options = entry.get('options', {})
        if not isinstance(options, Mapping) or any(not isinstance(key, str) for key in options):
            raise LegacyConfigurationError(f'{class_path}: malformed recorded declaration options for {entry["name"]!r}.')
        name = entry['name'].strip().lower()
        options = {key.lower(): value for key, value in options.items()}
        _validate_recorded_options(name, options, class_path)
        recorded.append((name, options))
    effective_names = [name for name, _ in recorded if name != 'plain']
    if effective_names != [entry['name'] for entry in entries]:
        raise LegacyConfigurationError(f'{class_path}: the recorded declaration {effective_names} does not describe the '
                                       f'effective stages {[entry["name"] for entry in entries]}.')
    return tuple(recorded)


def _is_load_invocation(options):
    """Whether an options mapping was written by a ``benchmark(load=...)`` invocation."""
    load_request = options.get('load')
    return options.get('operation') == 'load' or (isinstance(load_request, str) and load_request != '')


def _with_run_count(result, count):
    """Add ``n_runs`` to a replay result only when the file records one."""
    if count is not None:
        result['n_runs'] = count
    return result


def _load_invocation_replay(plain):
    """A file written by a load (re-plot) invocation replays only its recovered archived recipe."""
    archived = plain.get('archived_experiment')
    if isinstance(archived, Mapping) and archived.get('replayable') is True:
        specification = archived.get('feature_specification')
        problem_options = archived.get('problem_options')
        if not isinstance(specification, (list, tuple, Mapping)) or not isinstance(problem_options, Mapping):
            raise LegacyConfigurationError('The archived recipe of this load-written options file is incomplete '
                                           '(feature specification or problem options missing).')
        return _with_run_count({'feature': _historical_specification(specification),
                                'problem_options': dict(problem_options)}, archived.get('n_runs'))
    reason = plain.get('replay_reason') if isinstance(plain.get('replay_reason'), str) else 'no_archived_recipe_recorded'
    raise LegacyConfigurationError(
        'This options file was written by a load (re-plot) invocation and records no replayable archived '
        f'experiment ({reason}); it describes the load, not the archived execution. Replay the source '
        'experiment from its own test_log/options_refined.pkl instead.')


def replay_arguments(options, feature_name=None):
    """
    Map a trusted options dictionary to ``benchmark`` inputs.

    Supported layouts: ``options_user-v2`` (raw current input),
    ``options_refined-v2`` (native ``feature_specification``
    and top-level ``n_runs``); the 1.x bridge layout (``feature_specification``
    without a schema, run count among the flat keys); and the flat 1.x merge of
    profile, feature and problem options, which does not record the feature
    identity and is replayed only with an explicit ``feature_name`` (the flat
    stage options are then broadcast to that name). Returns ``{'feature',
    'problem_options'}`` plus ``'n_runs'`` when the file records a run count
    (a user file written without ``n_runs`` records none, and ``None`` is not a
    valid ``benchmark`` value); the caller supplies solvers and any profile
    options. Callback values stay the native objects stored in the file; nothing
    is reconstructed from JSON descriptions.

    A file written by a ``load`` (re-plot) invocation describes that load, not
    the archived execution. It replays only when it carries the archived recipe
    recovered by :func:`archived_replay_recipe` (``archived_experiment`` with
    ``replayable`` true); otherwise :class:`LegacyConfigurationError` names the
    reason. Case-variant duplicates of one option name in a flat 1.x file are
    rejected the same way (the configuration is ambiguous).
    """
    if not isinstance(options, Mapping):
        raise TypeError(f'replay_arguments expects an options mapping, not {type(options).__name__}.')
    plain = {}
    for key, value in options.items():
        key = _option_key(key)
        if not isinstance(key, str):
            raise LegacyConfigurationError(f'Option names must be strings, got {type(key).__name__}.')
        plain[key] = value
    if _is_load_invocation(plain):
        return _load_invocation_replay(plain)
    problem_keys = {member.value for member in ProblemOption}
    problem_options = {key: value for key, value in plain.items() if key in problem_keys}
    schema = plain.get('schema')
    if schema == USER_SCHEMA:
        # Raw user input is not a historical effective configuration. The
        # writer marker distinguishes fresh case-insensitive mesh input from
        # old unmarked files whose mixed-case grid actually ran as absolute.
        from .opclasses import Feature
        from .provenance import effective_specification
        stage_options = {key: value for key, value in plain.items()
                         if key in FeatureOption.__members__.values() and key not in EXPERIMENT_OPTIONS}
        try:
            if 'feature' in plain:
                if 'feature_name' in plain or stage_options:
                    raise ValueError('Structured feature input cannot be mixed with flat feature options.')
                feature = Feature(plain['feature'])
            else:
                name = feature_name if feature_name is not None else plain.get('feature_name', 'plain')
                feature = Feature(name, **stage_options)
        except (TypeError, ValueError) as err:
            raise LegacyConfigurationError(f'The saved user options are not valid: {err}') from err
        return _with_run_count({'feature': effective_specification(feature), 'problem_options': problem_options},
                               plain.get('n_runs'))
    if schema == REFINED_SCHEMA or (schema is None and 'feature_specification' in plain):
        if 'feature_specification' not in plain:
            raise LegacyConfigurationError(f'An {REFINED_SCHEMA} options file must record `feature_specification`; '
                                           'this one does not.')
        specification = plain['feature_specification']
        if not isinstance(specification, (list, tuple, Mapping)):
            raise LegacyConfigurationError('`feature_specification` must be a mapping or a list/tuple of stage entries.')
        return _with_run_count({'feature': _historical_specification(specification), 'problem_options': problem_options},
                               plain.get('n_runs'))
    if schema is not None:
        raise LegacyConfigurationError(f'Unknown refined options schema {schema!r}.')
    name = feature_name if feature_name is not None else plain.get('feature_name')
    if name is None:
        raise LegacyConfigurationError('This flat 1.x options file does not record the feature identity; pass '
                                       'feature_name=... to replay it (the flat stage options are broadcast to it).')
    stage_options = {key: value for key, value in plain.items()
                     if key in FeatureOption.__members__.values() and key not in EXPERIMENT_OPTIONS}
    if isinstance(name, str) and 'quantized' in [part.strip().lower() for part in name.split('+')]:
        stage_options = historical_effective_options('quantized', stage_options)
    from .opclasses import Feature
    from .provenance import effective_specification
    try:
        feature = Feature(name, **stage_options)
    except (TypeError, ValueError) as err:
        raise LegacyConfigurationError(f'The flat options are not valid for feature {name!r}: {err}') from err
    return _with_run_count({'feature': effective_specification(feature), 'problem_options': problem_options},
                           plain.get('n_runs'))


def _closed_recipe(reason, **details):
    """The recipe of a load whose archived experiment could not be recovered exactly."""
    return {'operation': 'load', 'replayable': False, 'replay_reason': reason,
            'feature_route': None, 'feature_name': None, 'feature_specification': None, 'n_runs': None,
            'archived_experiment': {'replayable': False, 'reason': reason, **details}}


def _recipe_from_archived(archived, results_plibs):
    """Cross-check a recovered archived experiment against the loaded archive and produce the recipe."""
    from .opclasses import Feature
    from .provenance import read_feature_pipeline
    try:
        feature = Feature(archived['feature_specification'])
    except (TypeError, ValueError) as err:
        return _closed_recipe('source_feature_specification_invalid', detail=str(err)[:200])
    n_runs = archived['n_runs']
    if archived.get('seed') is None:
        return _closed_recipe('source_seed_not_recorded')
    stage_names = [stage.name for stage in feature.stages]
    for result in results_plibs or []:
        if not isinstance(result, Mapping):
            continue
        histories = result.get('fun_histories')
        if getattr(histories, 'ndim', 0) == 4 and int(histories.shape[2]) != n_runs:
            return _closed_recipe('archive_run_axis_disagrees_with_source_options',
                                  archive_runs=int(histories.shape[2]), source_n_runs=n_runs)
        stamp = result.get('feature_stamp')
        if isinstance(stamp, str) and isinstance(archived.get('feature_stamp'), str) and stamp != archived['feature_stamp']:
            return _closed_recipe('archive_feature_stamp_disagrees_with_source_options')
        payload, schema = read_feature_pipeline(result.get('feature_pipeline'))
        if isinstance(payload, Mapping) and schema == 'feature_pipeline-v3':
            block = payload.get('feature') if isinstance(payload.get('feature'), Mapping) else {}
            experiment = payload.get('experiment') if isinstance(payload.get('experiment'), Mapping) else {}
            recorded = [stage.get('name') for stage in block.get('stages', []) if isinstance(stage, Mapping)]
            if recorded != stage_names or experiment.get('n_runs') not in (None, n_runs):
                return _closed_recipe('archive_feature_pipeline_disagrees_with_source_options')
        if isinstance(result.get('results_plib_plain'), Mapping) != bool(archived['run_plain']):
            return _closed_recipe('archive_plain_reference_disagrees_with_source_options')
    return {'operation': 'load', 'replayable': True, 'replay_reason': None, 'feature_route': None,
            'feature_name': feature.declared_name, 'feature_specification': archived['feature_specification'],
            'n_runs': n_runs, 'archived_experiment': dict(archived)}


def archived_replay_recipe(source_options_path, results_plibs):
    """
    The replay recipe of a loaded experiment, recovered from the source
    experiment's own native ``test_log/options_refined.pkl`` (exact native
    values, callables included) and cross-checked against the archive that was
    loaded. The report's JSON callback descriptions are never used.

    Returns the fields that ``benchmark(load=...)`` writes into its own
    ``options_refined.pkl``: ``operation='load'``, ``replayable``,
    ``replay_reason`` (``None`` when replayable), the top-level recipe
    (``feature_specification``, ``n_runs`` of the primary role, ``feature_name``,
    ``feature_route=None``) and ``archived_experiment`` (the recovered
    specification, primary ``n_runs``, ``seed``, ``run_plain`` with the fixed
    ``plain_reference_n_runs`` of the reference role, the archived problem
    options, the feature stamp and the source options schema). When the source
    file is missing, unreadable, a flat 1.x file without feature identity, itself
    a load without a recipe, or disagrees with the archive (run axis, feature
    stamp, ``feature_pipeline-v3`` stages or run count, plain reference), the
    recipe fails closed: ``feature_specification`` and ``n_runs`` are ``None``
    and ``replay_reason`` says why. Nothing is inferred from the load label.
    """
    from .experiment import validate_n_runs
    path = Path(source_options_path) if source_options_path is not None else None
    if path is None or not path.is_file():
        return _closed_recipe('source_options_refined_missing')
    try:
        source = load_options(path)
    except LegacyConfigurationError as err:
        cause = err.__cause__ if err.__cause__ is not None else err
        return _closed_recipe('source_options_unreadable', error_type=type(cause).__name__)
    except OSError as err:
        return _closed_recipe('source_options_unreadable', error_type=type(err).__name__)
    if _is_load_invocation(source):
        nested = source.get('archived_experiment')
        if isinstance(nested, Mapping) and nested.get('replayable') is True:
            return _recipe_from_archived(dict(nested), results_plibs)
        return _closed_recipe('source_is_a_load_invocation_without_recipe')
    try:
        replay = replay_arguments(source)
    except LegacyConfigurationError as err:
        return _closed_recipe('source_options_not_replayable', detail=str(err)[:200])
    try:
        n_runs = validate_n_runs(replay['n_runs']) if 'n_runs' in replay else None
    except (TypeError, ValueError):
        n_runs = None
    if n_runs is None:
        return _closed_recipe('source_run_count_not_recorded')
    run_plain = bool(source.get('run_plain', False))
    archived = {'replayable': True, 'reason': None,
                'source_options_schema': source.get('schema') if isinstance(source.get('schema'), str) else 'unversioned',
                'feature_specification': replay['feature'], 'n_runs': n_runs, 'seed': source.get('seed'),
                'run_plain': run_plain, 'plain_reference_n_runs': 1 if run_plain else None,
                'problem_options': replay['problem_options'],
                'feature_stamp': source.get('feature_stamp') if isinstance(source.get('feature_stamp'), str) else None}
    return _recipe_from_archived(archived, results_plibs)
