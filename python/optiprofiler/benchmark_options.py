"""Trusted native settings for fresh benchmark replay and archived recipes.

Replay runs new solvers. It neither resumes runtime state nor reconstructs
callbacks from report metadata. Historical decoding stays in legacy_compat.
"""

import os
import pickle
import warnings
from collections.abc import Mapping
from pathlib import Path

from .legacy_compat import (LegacyConfigurationError, load_trusted, _option_key,
                            _historical_specification, _legacy_feature_specification)
from .feature_definitions import EXPERIMENT_OPTIONS
from .utils import FeatureOption, ProblemOption
from .provenance import FEATURE_PIPELINE_SCHEMA

REFINED_SCHEMA = 'options_refined-v2'
USER_SCHEMA = 'options_user-v2'
_PROBLEM_OPTION_KEYS = frozenset(member.value for member in ProblemOption)


def load_benchmark_options(source, feature_name=None):
    """Read trusted native settings into keywords for a fresh benchmark.

    ``source`` is an options pickle path or an already-loaded mapping. Supports
    ``options_user-v2``, ``options_refined-v2``, the unversioned structured
    bridge and flat 1.x settings. Flat records without feature identity require
    ``feature_name``; stage options use the shorthand's broadcast rules.
    Prefer refined settings to retain resolved defaults and run counts.

    Returns ``feature``, problem-selection keywords (or the original
    ``problem``), and ``n_runs`` only when recorded. Supply solvers and profile
    options, including seed/merit/silent, explicitly. Native callbacks are kept
    as saved objects, never reconstructed from JSON descriptions.

    A direct-Problem refined record needs the original ``problem`` supplied in
    the mapping; otherwise LegacyConfigurationError is raised rather than
    selecting a library. Records predating ``problem_route`` warn when their
    problem source cannot be established. A load-written record replays only
    through its recovered ``archived_experiment`` recipe, not its load label.

    Pickle files must be trusted: decoding can execute user code. Unreadable
    configurations raise LegacyConfigurationError; filesystem errors remain
    OSError; a source that is neither a mapping nor a path raises TypeError.
    Replay runs a new experiment, not a live-state resume or replot.
    """
    if not isinstance(source, (Mapping, str, os.PathLike)):
        raise TypeError('load_benchmark_options expects an options mapping or the path of an options pickle, '
                        f'not {type(source).__name__}.')
    options = source if isinstance(source, Mapping) else _read_options(source)
    return _replay_options(options, feature_name)


def _read_options(path):
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


def _is_load_invocation(options):
    """Whether an options mapping was written by a ``benchmark(load=...)`` invocation."""
    load_request = options.get('load')
    return options.get('operation') == 'load' or (isinstance(load_request, str) and load_request != '')


def _with_run_count(result, count):
    """Add ``n_runs`` to a replay result only when the file records one."""
    if count is not None:
        result['n_runs'] = count
    return result


def _validated_archived(archived, require_flag=True):
    """
    A recovered ``archived_experiment`` record with every field a replay needs,
    or ``(None, reason)`` when it is missing or malformed (a stray key, a
    missing count, a non-boolean ``run_plain``): malformed records fail closed
    instead of surfacing bare exceptions. A *stored* record (a load-written
    file read back) must also say ``replayable`` is true; a candidate handed
    to the cross-check is judged by the cross-check itself.
    """
    from .experiment import validate_n_runs
    if not isinstance(archived, Mapping):
        return None, 'source_recipe_malformed:not_a_mapping'
    if require_flag and archived.get('replayable') is not True:
        return None, 'source_recipe_not_replayable'
    if not isinstance(archived.get('feature_specification'), (list, tuple, Mapping)):
        return None, 'source_recipe_malformed:feature_specification'
    try:
        n_runs = validate_n_runs(archived.get('n_runs'))
    except (TypeError, ValueError):
        return None, 'source_recipe_malformed:n_runs'
    seed = archived.get('seed')
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        return None, 'source_recipe_malformed:seed'
    if not isinstance(archived.get('run_plain'), bool):
        return None, 'source_recipe_malformed:run_plain'
    if not isinstance(archived.get('problem_options'), Mapping):
        return None, 'source_recipe_malformed:problem_options'
    validated = dict(archived)
    validated['n_runs'] = n_runs
    validated['problem_options'] = dict(archived['problem_options'])
    return validated, None


def _load_invocation_replay(plain):
    """A file written by a load (re-plot) invocation replays only its recovered archived recipe."""
    archived, problem = _validated_archived(plain.get('archived_experiment'))
    if archived is not None:
        return _with_run_count({**archived['problem_options'],
                                'feature': _historical_specification(archived['feature_specification'])}, archived['n_runs'])
    if problem != 'source_recipe_not_replayable' and problem != 'source_recipe_malformed:not_a_mapping':
        raise LegacyConfigurationError(f'The archived recipe of this load-written options file is malformed ({problem}); '
                                       'replay the source experiment from its own test_log/options_refined.pkl instead.')
    reason = plain.get('replay_reason') if isinstance(plain.get('replay_reason'), str) else 'no_archived_recipe_recorded'
    raise LegacyConfigurationError(
        'This options file was written by a load (re-plot) invocation and records no replayable archived '
        f'experiment ({reason}); it describes the load, not the archived execution. Replay the source '
        'experiment from its own test_log/options_refined.pkl instead.')


def _replay_options(options, feature_name=None, *, _warn_route=True):
    """Map a decoded native configuration to fresh benchmark keywords."""
    plain = {}
    for key, value in options.items():
        key = _option_key(key)
        if not isinstance(key, str):
            raise LegacyConfigurationError(f'Option names must be strings, got {type(key).__name__}.')
        plain[key] = value
    if _is_load_invocation(plain):
        return _load_invocation_replay(plain)
    route = plain.get('problem_route')
    if route not in (None, 'direct_problem', 'library_selection'):
        raise LegacyConfigurationError(f"`problem_route` must be 'direct_problem' or 'library_selection', not {route!r}.")
    if plain.get('problem') is not None:
        from .opclasses import Problem
        if not isinstance(plain['problem'], Problem):
            raise LegacyConfigurationError(f'`problem` must be a Problem, not {type(plain["problem"]).__name__}.')
        # The same Problem, never a library selection alongside it.
        problem_options = {'problem': plain['problem']}
    elif route == 'direct_problem':
        raise LegacyConfigurationError(
            'These options come from a benchmark of a directly supplied Problem, which options_refined.pkl '
            'does not store; a replay without it would benchmark a problem library instead. Set the original '
            "Problem first (options['problem'] = problem) and call load_benchmark_options again.")
    else:
        problem_options = {key: value for key, value in plain.items() if key in _PROBLEM_OPTION_KEYS}
        # options_user files hold the Problem of a direct run themselves.
        if _warn_route and route is None and plain.get('schema') != USER_SCHEMA:
            warnings.warn(
                'These options do not record whether the experiment selected its problems from a library or '
                'used a directly supplied Problem. If it used a Problem, this replay would benchmark the recorded '
                "library selectors instead: set the original Problem first (options['problem'] = problem).",
                stacklevel=3)
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
        return _with_run_count({**problem_options, 'feature': effective_specification(feature)},
                               plain.get('n_runs'))
    if schema == REFINED_SCHEMA or (schema is None and 'feature_specification' in plain):
        if 'feature_specification' not in plain:
            raise LegacyConfigurationError(f'An {REFINED_SCHEMA} options file must record `feature_specification`; '
                                           'this one does not.')
        specification = plain['feature_specification']
        if not isinstance(specification, (list, tuple, Mapping)):
            raise LegacyConfigurationError('`feature_specification` must be a mapping or a list/tuple of stage entries.')
        return _with_run_count({**problem_options, 'feature': _historical_specification(specification)},
                               plain.get('n_runs'))
    if schema is not None:
        raise LegacyConfigurationError(f'Unknown refined options schema {schema!r}.')
    specification = _legacy_feature_specification(plain, feature_name)
    return _with_run_count({**problem_options, 'feature': specification}, plain.get('n_runs'))


def _closed_recipe(reason, **details):
    """The recipe of a load whose archived experiment could not be recovered exactly."""
    return {'operation': 'load', 'replayable': False, 'replay_reason': reason,
            'feature_route': None, 'feature_name': None, 'feature_specification': None, 'n_runs': None,
            'archived_experiment': {'replayable': False, 'reason': reason, **details}}


def _is_callback_descriptor(value):
    return isinstance(value, Mapping) and value.get('kind') == 'callback'


def _projection_equal(archived, current):
    """
    Equality of two safe-metadata projections of a stage option. Numbers
    compare by value, containers element-wise. Two callback descriptors compare
    by module and name only: the archive keeps no executable callback, so the
    descriptor is the only fact available and it is not a semantic identity.
    Returns ``(equal, callbacks_compared)``.
    """
    if _is_callback_descriptor(archived) or _is_callback_descriptor(current):
        if not (_is_callback_descriptor(archived) and _is_callback_descriptor(current)):
            return False, True
        same = archived.get('module') == current.get('module') and archived.get('name') == current.get('name')
        return same, True
    if isinstance(archived, bool) or isinstance(current, bool):
        return type(archived) is type(current) and archived == current, False
    if isinstance(archived, (int, float)) and isinstance(current, (int, float)):
        return archived == current, False
    if isinstance(archived, list) and isinstance(current, list):
        if len(archived) != len(current):
            return False, False
        callbacks = False
        for left, right in zip(archived, current):
            same, seen = _projection_equal(left, right)
            callbacks = callbacks or seen
            if not same:
                return False, callbacks
        return True, callbacks
    if isinstance(archived, Mapping) and isinstance(current, Mapping):
        if set(archived) != set(current):
            return False, False
        callbacks = False
        for key in archived:
            same, seen = _projection_equal(archived[key], current[key])
            callbacks = callbacks or seen
            if not same:
                return False, callbacks
        return True, callbacks
    return type(archived) is type(current) and archived == current, False


def _stages_disagreement(recorded, feature):
    """
    Compare the stage records of an archived ``feature_pipeline-v3`` block with
    the recovered specification, projected by the same metadata encoder:
    stage count, name, occurrence, identity, code, position and every
    stage-local option. Returns ``(detail, callbacks_compared)`` where
    ``detail`` names the first disagreement or is ``None``.
    """
    from .provenance import describe_feature
    described = describe_feature(feature)['stages']
    if not isinstance(recorded, list) or len(recorded) != len(described):
        return 'stage_count', False
    callbacks = False
    for archived, current in zip(recorded, described):
        if not isinstance(archived, Mapping):
            return 'stage_record', callbacks
        for key in ('name', 'occurrence', 'identity', 'code', 'position'):
            if key in archived and archived[key] != current[key]:
                return f'stage_{key}', callbacks
        options = archived.get('options')
        if not isinstance(options, Mapping):
            return 'stage_options', callbacks
        if set(options) != set(current['options']):
            return 'stage_option_keys', callbacks
        for name in options:
            same, seen = _projection_equal(options[name], current['options'][name])
            callbacks = callbacks or seen
            if not same:
                return f'stage_option:{name}', callbacks
    return None, callbacks


def _recipe_from_archived(archived, results_plibs):
    """
    Cross-check a recovered archived experiment against the loaded archive and
    produce the recipe. The primary question comes first: does the source
    describe the archived experiment? Agreement is required on every fact the
    archive retains: run axis, feature stamp, and for a ``feature_pipeline-v3``
    payload the stage records including their options (callbacks by descriptor
    module and name only, which is recorded as such) and the experiment run
    count and effective name; the plain reference must match ``run_plain``.
    Only a candidate that agrees is then required to be complete enough to
    replay (``source_recipe_malformed:<field>`` otherwise); a malformed field
    the cross-check itself needs is skipped there and reported here.
    """
    from .experiment import validate_n_runs
    from .opclasses import Feature
    from .provenance import read_feature_pipeline
    if not isinstance(archived, Mapping):
        return _closed_recipe('source_recipe_malformed:not_a_mapping')
    if not isinstance(archived.get('feature_specification'), (list, tuple, Mapping)):
        return _closed_recipe('source_recipe_malformed:feature_specification')
    try:
        feature = Feature(archived['feature_specification'])
    except (TypeError, ValueError) as err:
        return _closed_recipe('source_feature_specification_invalid', detail=str(err)[:200])
    try:
        n_runs = validate_n_runs(archived.get('n_runs'))
    except (TypeError, ValueError):
        n_runs = None
    run_plain = archived.get('run_plain') if isinstance(archived.get('run_plain'), bool) else None
    checked = set()
    comparison = {'callback_comparison': 'not_applicable'}
    for result in results_plibs or []:
        if not isinstance(result, Mapping):
            continue
        histories = result.get('fun_histories')
        if n_runs is not None and getattr(histories, 'ndim', 0) == 4:
            if int(histories.shape[2]) != n_runs:
                return _closed_recipe('archive_run_axis_disagrees_with_source_options',
                                      archive_runs=int(histories.shape[2]), source_n_runs=n_runs)
            checked.add('run_axis')
        stamp = result.get('feature_stamp')
        if isinstance(stamp, str) and isinstance(archived.get('feature_stamp'), str):
            if stamp != archived['feature_stamp']:
                return _closed_recipe('archive_feature_stamp_disagrees_with_source_options')
            checked.add('feature_stamp')
        payload, schema = read_feature_pipeline(result.get('feature_pipeline'))
        if isinstance(payload, Mapping) and schema == FEATURE_PIPELINE_SCHEMA:
            block = payload.get('feature') if isinstance(payload.get('feature'), Mapping) else {}
            experiment = payload.get('experiment') if isinstance(payload.get('experiment'), Mapping) else {}
            detail, callbacks = _stages_disagreement(block.get('stages'), feature)
            if detail is None and n_runs is not None and experiment.get('n_runs') not in (None, n_runs):
                detail = 'experiment_n_runs'
            if detail is None and block.get('effective_name') not in (None, feature.name):
                detail = 'effective_name'
            if detail is not None:
                return _closed_recipe('archive_feature_pipeline_disagrees_with_source_options', detail=detail)
            checked.add('feature_pipeline_v3_stages_and_options')
            if callbacks:
                comparison['callback_comparison'] = 'descriptor_module_and_name_only_not_identity_or_semantics'
        if run_plain is not None:
            if isinstance(result.get('results_plib_plain'), Mapping) != run_plain:
                return _closed_recipe('archive_plain_reference_disagrees_with_source_options')
            if 'results_plib_plain' in result:
                checked.add('plain_reference')
    comparison['archive_comparison'] = ';'.join(
        field for field in ('run_axis', 'feature_stamp', 'plain_reference',
                            'feature_pipeline_v3_stages_and_options') if field in checked
    ) or 'not_performed'
    archived, problem = _validated_archived(archived, require_flag=False)
    if archived is None:
        return _closed_recipe(problem)
    return {'operation': 'load', 'replayable': True, 'replay_reason': None, 'feature_route': None,
            'feature_name': feature.declared_name, 'feature_specification': archived['feature_specification'],
            'n_runs': archived['n_runs'], 'archived_experiment': {**dict(archived), **comparison}}


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
        source = _read_options(path)
    except LegacyConfigurationError as err:
        cause = err.__cause__ if err.__cause__ is not None else err
        return _closed_recipe('source_options_unreadable', error_type=type(cause).__name__)
    except OSError as err:
        return _closed_recipe('source_options_unreadable', error_type=type(err).__name__)
    if _is_load_invocation(source):
        nested, problem = _validated_archived(source.get('archived_experiment'))
        if nested is not None:
            return _recipe_from_archived(nested, results_plibs)
        if problem in ('source_recipe_not_replayable', 'source_recipe_malformed:not_a_mapping'):
            return _closed_recipe('source_is_a_load_invocation_without_recipe')
        return _closed_recipe(problem)
    try:
        # A load replays nothing, and a direct Problem leaves no archive to load.
        replay = _replay_options(source, _warn_route=False)
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
                'problem_options': {key: value for key, value in replay.items()
                                    if key in _PROBLEM_OPTION_KEYS},
                'feature_stamp': source.get('feature_stamp') if isinstance(source.get('feature_stamp'), str) else None}
    return _recipe_from_archived(archived, results_plibs)
