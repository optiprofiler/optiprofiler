"""
Pure definitions of the feature stages: names, stage-local options, validators,
local defaults, replicate hints, stochasticity, stamps and the frozen seed codes,
plus the normalization of user input into the canonical stage records.

This module imports nothing from the problem, view, benchmark or report code so
that every layer (specification, experiment planning, execution, provenance,
documentation and the MATLAB mapping) reads one table. Nothing here executes a
user callback: callables are checked for callability only and stored by
reference.

Experiment options (``n_runs``) are not stage options; they are rejected inside a
stage and belong to ``optiprofiler.experiment``.
"""

from collections.abc import Mapping
from types import MappingProxyType

import numpy as np

from .utils import FeatureName, FeatureOption

#: Options that belong to the experiment, never to a stage.
EXPERIMENT_OPTIONS = frozenset({'n_runs'})

#: Frozen literal seed codes of the stage kinds (policy ``seedsequence-v2``).
#: Never derived from enumeration or table order; changing a code changes experiments.
STAGE_CODES = {
    FeatureName.PERTURBED_X0.value: 1,
    FeatureName.NOISY.value: 2,
    FeatureName.TRUNCATED.value: 3,
    FeatureName.PERMUTED.value: 4,
    FeatureName.LINEARLY_TRANSFORMED.value: 5,
    FeatureName.RANDOM_NAN.value: 6,
    FeatureName.UNRELAXABLE_CONSTRAINTS.value: 7,
    FeatureName.NONQUANTIFIABLE_CONSTRAINTS.value: 8,
    FeatureName.QUANTIZED.value: 9,
    FeatureName.CUSTOM.value: 10,
}

#: Seed policies recorded in provenance (values are frozen identifiers).
SEED_POLICY_SINGLE = 'legacy-run-seed'
SEED_POLICY_COMPOSED = 'seedsequence-v2'

_CUSTOM_CALLBACKS = (FeatureOption.MOD_X0, FeatureOption.MOD_BOUNDS, FeatureOption.MOD_LINEAR_UB,
                     FeatureOption.MOD_LINEAR_EQ, FeatureOption.MOD_AFFINE, FeatureOption.MOD_FUN,
                     FeatureOption.MOD_CUB, FeatureOption.MOD_CEQ)

#: Stage-local options accepted by each kind. ``noisy.distribution`` and
#: ``perturbed_x0.distribution`` are distinct definitions that share a spelling.
LOCAL_OPTIONS = {
    FeatureName.PLAIN.value: (),
    FeatureName.CUSTOM.value: tuple(option.value for option in _CUSTOM_CALLBACKS),
    FeatureName.NOISY.value: (FeatureOption.DISTRIBUTION.value, FeatureOption.NOISE_LEVEL.value,
                              FeatureOption.NOISE_TYPE.value, FeatureOption.NOISE_MODE.value,
                              FeatureOption.NOISE_MAP.value),
    FeatureName.PERTURBED_X0.value: (FeatureOption.DISTRIBUTION.value, FeatureOption.PERTURBATION_LEVEL.value),
    FeatureName.RANDOM_NAN.value: (FeatureOption.NAN_RATE.value,),
    FeatureName.TRUNCATED.value: (FeatureOption.PERTURBED_TRAILING_DIGITS.value, FeatureOption.SIGNIFICANT_DIGITS.value),
    FeatureName.UNRELAXABLE_CONSTRAINTS.value: (FeatureOption.UNRELAXABLE_BOUNDS.value,
                                                FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS.value,
                                                FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value),
    FeatureName.LINEARLY_TRANSFORMED.value: (FeatureOption.ROTATED.value, FeatureOption.CONDITION_FACTOR.value),
    FeatureName.QUANTIZED.value: (FeatureOption.MESH_SIZE.value, FeatureOption.MESH_TYPE.value,
                                  FeatureOption.GROUND_TRUTH.value),
    FeatureName.PERMUTED.value: (),
    FeatureName.NONQUANTIFIABLE_CONSTRAINTS.value: (),
}

STAGE_NAMES = tuple(LOCAL_OPTIONS)


def local_options(name):
    """The stage-local options accepted by the stage kind ``name``."""
    try:
        return LOCAL_OPTIONS[name]
    except KeyError:
        raise ValueError(f'Unknown feature: {name}.') from None


def validate_option(name, key, value):
    """Validate one stage-local option of kind ``name`` and return its normalized value."""
    if key in EXPERIMENT_OPTIONS:
        raise ValueError(f'Option `{key}` is experiment-wide and is not a stage option; '
                         f'give it at the top level, not inside a stage.')
    if key not in FeatureOption.__members__.values():
        raise ValueError(f'Unknown option for feature: {key}.')
    if key not in local_options(name):
        raise ValueError(f"Option `{key}` is not valid for feature '{name}'.")
    if key == FeatureOption.DISTRIBUTION:
        if isinstance(value, str):
            if name == FeatureName.NOISY and value not in ['gaussian', 'uniform']:
                raise ValueError(f'Option `{key}` for feature `{name}` must be either "gaussian" or "uniform" when specified as a string.')
            elif name == FeatureName.PERTURBED_X0 and value not in ['gaussian', 'spherical']:
                raise ValueError(f'Option `{key}` for feature `{name}` must be either "gaussian" or "spherical" when specified as a string.')
        elif not callable(value):
            raise TypeError(f'Option `{key}` must be a string or it must be callable.')
    elif key == FeatureOption.NAN_RATE:
        if not isinstance(value, (int, float)):
            raise TypeError(f'Option `{key}` must be a number.')
        if not (0.0 <= value <= 1.0):
            raise ValueError(f'Option `{key}` must be between 0 and 1.')
    elif key == FeatureOption.SIGNIFICANT_DIGITS:
        if isinstance(value, (float, np.floating)) and float(value).is_integer():
            value = int(value)
        if isinstance(value, np.integer):
            value = int(value)
        if not isinstance(value, int):
            raise TypeError(f'Option `{key}` must be an integer.')
        if value <= 0:
            raise ValueError(f'Option `{key}` must be positive.')
    elif key in [FeatureOption.NOISE_LEVEL, FeatureOption.CONDITION_FACTOR]:
        if not isinstance(value, (int, float)):
            raise TypeError(f'Option `{key}` must be a number.')
        if value < 0.0:
            raise ValueError(f'Option `{key}` must be nonnegative.')
    elif key == FeatureOption.NOISE_TYPE:
        if not isinstance(value, str):
            raise TypeError(f'Option {key} must be a string.')
        if value.lower() not in ['absolute', 'relative', 'mixed']:
            raise ValueError(f"Option `{key}` must be one of 'absolute', 'relative', or 'mixed'.")
        value = value.lower()
    elif key == FeatureOption.NOISE_MODE:
        if not isinstance(value, str):
            raise TypeError(f'Option {key} must be a string.')
        if value.lower() not in ['random', 'deterministic']:
            raise ValueError(f"Option `{key}` must be either 'random' or 'deterministic'.")
        value = value.lower()
    elif key == FeatureOption.NOISE_MAP:
        if isinstance(value, str):
            if value.lower() != 'chebyshev':
                raise ValueError(f'Option `{key}` must be "chebyshev" when specified as a string.')
            value = value.lower()
        elif not callable(value):
            raise TypeError(f'Option `{key}` must be a string or it must be callable.')
    elif key in [FeatureOption.PERTURBED_TRAILING_DIGITS, FeatureOption.ROTATED, FeatureOption.UNRELAXABLE_BOUNDS,
                 FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS, FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS,
                 FeatureOption.GROUND_TRUTH]:
        if not isinstance(value, bool):
            raise TypeError(f'Option `{key}` must be a boolean.')
    elif key == FeatureOption.MESH_SIZE:
        if not isinstance(value, (int, float)):
            raise TypeError(f'Option `{key}` must be a number.')
        if value <= 0.0:
            raise ValueError(f'Option `{key}` must be positive.')
    elif key == FeatureOption.MESH_TYPE:
        if not isinstance(value, str):
            raise TypeError(f'Option `{key}` must be a string.')
        if value.lower() not in ['absolute', 'relative']:
            raise ValueError(f"Option `{key}` must be 'absolute' or 'relative'.")
    elif key in _CUSTOM_CALLBACKS:
        if not callable(value):
            raise TypeError(f'Option `{key}` must be callable.')
    return value


def apply_local_defaults(name, options):
    """Fill the unspecified stage-local options of kind ``name`` in place (values are plain data)."""
    if name in [FeatureName.PLAIN, FeatureName.CUSTOM, FeatureName.NONQUANTIFIABLE_CONSTRAINTS, FeatureName.PERMUTED]:
        pass
    elif name == FeatureName.NOISY:
        options.setdefault(FeatureOption.NOISE_MODE.value, 'random')
        options.setdefault(FeatureOption.DISTRIBUTION.value, 'gaussian')
        options.setdefault(FeatureOption.NOISE_MAP.value, 'chebyshev')
        options.setdefault(FeatureOption.NOISE_LEVEL.value, 1e-3)
        options.setdefault(FeatureOption.NOISE_TYPE.value, 'mixed')
    elif name == FeatureName.LINEARLY_TRANSFORMED:
        options.setdefault(FeatureOption.ROTATED.value, True)
        options.setdefault(FeatureOption.CONDITION_FACTOR.value, 0)
    elif name == FeatureName.PERTURBED_X0:
        options.setdefault(FeatureOption.DISTRIBUTION.value, 'spherical')
        options.setdefault(FeatureOption.PERTURBATION_LEVEL.value, 1e-3)
    elif name == FeatureName.RANDOM_NAN:
        options.setdefault(FeatureOption.NAN_RATE.value, 0.05)
    elif name == FeatureName.TRUNCATED:
        options.setdefault(FeatureOption.PERTURBED_TRAILING_DIGITS.value, False)
        options.setdefault(FeatureOption.SIGNIFICANT_DIGITS.value, 6)
    elif name == FeatureName.UNRELAXABLE_CONSTRAINTS:
        options.setdefault(FeatureOption.UNRELAXABLE_BOUNDS.value, True)
        options.setdefault(FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS.value, False)
        options.setdefault(FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value, False)
    elif name == FeatureName.QUANTIZED:
        options.setdefault(FeatureOption.MESH_SIZE.value, 1e-3)
        options.setdefault(FeatureOption.MESH_TYPE.value, 'absolute')
        options.setdefault(FeatureOption.GROUND_TRUTH.value, True)
    else:
        raise NotImplementedError(f'Unknown feature: {name}.')
    return options


def replicate_hint(name, options):
    """
    The established default number of runs of one stage kind, from its validated
    local options. A literal table, not a function of stochasticity: ``custom``
    hints one run although it is stochastic; an unrotated linear transformation
    hints five although it is deterministic. The experiment resolver combines the
    hints of the effective stages; a stage never stores a run count.
    """
    if name in [FeatureName.PLAIN, FeatureName.CUSTOM, FeatureName.NONQUANTIFIABLE_CONSTRAINTS,
                FeatureName.UNRELAXABLE_CONSTRAINTS, FeatureName.QUANTIZED]:
        return 1
    elif name == FeatureName.NOISY:
        return 1 if options[FeatureOption.NOISE_MODE] == 'deterministic' else 5
    elif name == FeatureName.TRUNCATED:
        return 5 if options[FeatureOption.PERTURBED_TRAILING_DIGITS] else 1
    elif name in [FeatureName.PERMUTED, FeatureName.LINEARLY_TRANSFORMED, FeatureName.PERTURBED_X0, FeatureName.RANDOM_NAN]:
        return 5
    raise NotImplementedError(f'Unknown feature: {name}.')


def is_stochastic(name, options):
    """Whether the stage kind ``name`` with these local options draws random numbers."""
    if name == FeatureName.NOISY:
        return options[FeatureOption.NOISE_MODE] == 'random'
    elif name in [FeatureName.PERTURBED_X0, FeatureName.PERMUTED, FeatureName.RANDOM_NAN, FeatureName.CUSTOM]:
        return True
    elif name == FeatureName.TRUNCATED:
        return bool(options[FeatureOption.PERTURBED_TRAILING_DIGITS])
    elif name == FeatureName.LINEARLY_TRANSFORMED:
        return bool(options[FeatureOption.ROTATED])
    return False


def stamp(name, options):
    """The folder-name stamp of one stage kind (the established 1.x conventions)."""
    if name == FeatureName.PERTURBED_X0:
        text = f"{name}_{options[FeatureOption.PERTURBATION_LEVEL]}"
        dist = options.get(FeatureOption.DISTRIBUTION.value)
        if isinstance(dist, str) and dist in ('gaussian', 'spherical'):
            text = f"{text}_{dist}"
    elif name == FeatureName.NOISY:
        text = f"{name}_{options[FeatureOption.NOISE_LEVEL]}_{options[FeatureOption.NOISE_TYPE]}"
        if options[FeatureOption.NOISE_MODE] == 'deterministic':
            text = f"{text}_deterministic"
            noise_map = options.get(FeatureOption.NOISE_MAP.value)
            if isinstance(noise_map, str) and noise_map == 'chebyshev':
                text = f"{text}_{noise_map}"
        else:
            dist = options.get(FeatureOption.DISTRIBUTION.value)
            if isinstance(dist, str) and dist in ('gaussian', 'uniform'):
                text = f"{text}_{dist}"
    elif name == FeatureName.TRUNCATED:
        text = f"{name}_{options[FeatureOption.SIGNIFICANT_DIGITS]}"
        if options[FeatureOption.PERTURBED_TRAILING_DIGITS]:
            text = f"{text}_perturbed_trailing_digits"
    elif name == FeatureName.LINEARLY_TRANSFORMED:
        text = str(name)
        if options[FeatureOption.ROTATED]:
            text = f"{text}_rotated"
        if options[FeatureOption.CONDITION_FACTOR] != 0:
            text = f"{text}_cond{options[FeatureOption.CONDITION_FACTOR]}"
    elif name == FeatureName.RANDOM_NAN:
        text = f"{name}_{options[FeatureOption.NAN_RATE]}"
    elif name == FeatureName.UNRELAXABLE_CONSTRAINTS:
        text = str(name)
        if options[FeatureOption.UNRELAXABLE_BOUNDS]:
            text = f"{text}_bounds"
        if options[FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS]:
            text = f"{text}_linear"
        if options[FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS]:
            text = f"{text}_nonlinear"
    elif name == FeatureName.QUANTIZED:
        text = f"{name}_{options[FeatureOption.MESH_SIZE]}"
        if options[FeatureOption.GROUND_TRUTH]:
            text = f"{text}_ground_truth"
    else:
        text = str(name)
    return text


def own_value(value):
    """
    Defensive copy of an option value at ingress: exact built-in containers
    (``list``, ``tuple``, ``dict``) are copied recursively and exact
    ``numpy.ndarray`` values become read-only owned copies (an object-dtype
    array is copied, its elements are not), so a specification never aliases
    user-owned ordinary data. Callables are checked first and always keep
    their identity (a callable container subclass stays the opaque callback it
    is); scalars, strings, subclasses of the containers and any other object
    keep their identity too. No copy, iteration or reduction hook of a user
    object is invoked.
    """
    if callable(value):
        return value
    kind = type(value)
    if kind is np.ndarray:
        copied = np.array(value, copy=True, subok=False)
        copied.setflags(write=False)
        return copied
    if kind is list:
        return [own_value(item) for item in value]
    if kind is tuple:
        return tuple(own_value(item) for item in value)
    if kind is dict:
        return {key: own_value(item) for key, item in value.items()}
    return value


def view_value(value):
    """
    Isolated read-only view of a stored option value: exact ``numpy.ndarray``
    values are returned as fresh read-only copies (the stored array is never
    exposed, so re-enabling writes on the returned object cannot reach the
    specification), exact lists and tuples as fresh copies, exact dictionaries
    as read-only mappings over fresh copies; callables and everything else as
    they are.
    """
    if callable(value):
        return value
    kind = type(value)
    if kind is np.ndarray:
        copied = np.array(value, copy=True, subok=False)
        copied.setflags(write=False)
        return copied
    if kind is list:
        return [view_value(item) for item in value]
    if kind is tuple:
        return tuple(view_value(item) for item in value)
    if kind is dict:
        return MappingProxyType({key: view_value(item) for key, item in value.items()})
    return value


def own_options(options):
    """A private, owned copy of an option mapping (see ``own_value``)."""
    return {key: own_value(value) for key, value in dict(options).items()}


def view_options(options):
    """A read-only view of an owned option mapping (see ``view_value``)."""
    return MappingProxyType({key: view_value(value) for key, value in options.items()})


def validated_local_options(name, options):
    """Validate and default an owned copy of the stage-local ``options`` of kind ``name``."""
    if any(not isinstance(key, str) for key in options):
        raise TypeError('option names must be strings.')
    validated = {}
    for key, value in options.items():
        key = key.lower()
        validated[key] = own_value(validate_option(name, key, value))
    return apply_local_defaults(name, validated)


class StageRecord:
    """
    One effective stage of a pipeline specification: its kind, occurrence index
    among stages of the same kind, identity, frozen seed code and validated
    stage-local options. Immutable; the options are exposed through a read-only
    view of a private copy, so the record pickles and is never shared state.
    """

    __slots__ = ('_name', '_occurrence', '_options', '_position')

    def __init__(self, name, occurrence, options, position=0):
        object.__setattr__(self, '_name', name)
        object.__setattr__(self, '_occurrence', occurrence)
        object.__setattr__(self, '_options', own_options(options))
        object.__setattr__(self, '_position', position)

    def __setattr__(self, key, value):
        raise AttributeError('StageRecord is immutable.')

    def __delattr__(self, key):
        raise AttributeError('StageRecord is immutable.')

    def __getstate__(self):
        return (self._name, self._occurrence, self._options, self._position)

    def __setstate__(self, state):
        name, occurrence, options, position = state
        object.__setattr__(self, '_name', name)
        object.__setattr__(self, '_occurrence', occurrence)
        object.__setattr__(self, '_options', own_options(options))
        object.__setattr__(self, '_position', position)

    @property
    def name(self):
        return self._name

    @property
    def position(self):
        """Index of the stage among the effective stages of its pipeline."""
        return self._position

    @property
    def occurrence(self):
        return self._occurrence

    @property
    def identity(self):
        return f'{self._name}#{self._occurrence}'

    @property
    def code(self):
        return STAGE_CODES[self._name]

    @property
    def options(self):
        """Read-only view of the validated stage-local options (arrays read-only, containers copied)."""
        return view_options(self._options)

    def native_options(self):
        """An owned plain copy of the validated stage-local options (native values, for transport and replay)."""
        return own_options(self._options)

    @property
    def is_stochastic(self):
        return is_stochastic(self._name, self._options)

    @property
    def replicate_hint(self):
        return replicate_hint(self._name, self._options)

    @property
    def stamp(self):
        return stamp(self._name, self._options)

    def __repr__(self):
        return f'StageRecord({self.identity!r}, {sorted(self._options)})'


class Declaration:
    """
    The specification as the user declared it: the input route (``'feature_name'``
    for the shorthand string, ``'feature'`` for structured entries) and the
    entries as given. A specification imported from a historical object that
    recorded no declaration has ``route`` ``None`` and no entries: its effective
    stages are known, its declaration is unknown and never fabricated.
    """

    __slots__ = ('_route', '_entries')

    def __init__(self, route, entries):
        if route not in ('feature_name', 'feature', None):
            raise ValueError(f'Unknown declaration route {route!r}.')
        entries = tuple((name, own_options(options)) for name, options in entries)
        if route is None and entries:
            raise ValueError('An unknown declaration carries no entries.')
        object.__setattr__(self, '_route', route)
        object.__setattr__(self, '_entries', entries)

    def __setattr__(self, key, value):
        raise AttributeError('Declaration is immutable.')

    def __getstate__(self):
        return (self._route, self._entries)

    def __setstate__(self, state):
        route, entries = state
        object.__setattr__(self, '_route', route)
        object.__setattr__(self, '_entries', tuple((name, own_options(options)) for name, options in entries))

    @property
    def route(self):
        return self._route

    @property
    def entries(self):
        """The entries as given, as read-only views (arrays read-only, containers copied)."""
        return tuple((name, view_options(options)) for name, options in self._entries)

    def native_entries(self):
        """Owned plain copies of the entries (native values, for transport)."""
        return tuple((name, own_options(options)) for name, options in self._entries)

    @property
    def name(self):
        """The declared name (``plain`` tokens included), or ``None`` when the declaration is unknown."""
        if self._route is None:
            return None
        return '+'.join(name for name, _ in self._entries)

    @property
    def is_known(self):
        return self._route is not None

    def __repr__(self):
        return f'Declaration({self._route!r}, {self.name!r})'


_SPEC_TYPE_MESSAGE = ('The first input argument for `Feature` must be a feature name string, a structured '
                      'specification (a mapping or a list/tuple of stage entries) or a Feature.')


def parse_feature_name(name):
    """
    Parse a feature name into its declared tokens and its effective tokens.

    Tokens are separated by ``+``, lowercased and stripped of surrounding
    whitespace. Empty or unknown tokens are rejected. ``plain`` tokens are kept
    in the declared list but removed from the effective list.
    """
    if not isinstance(name, str):
        raise TypeError('The first input argument for `Feature` must be a string.')
    tokens = [token.strip().lower() for token in name.split('+')]
    if any(token == '' for token in tokens):
        raise ValueError(f'Invalid feature name {name!r}: empty stage token in a "+"-separated composition.')
    for token in tokens:
        if token not in FeatureName.__members__.values():
            raise ValueError(f'Unknown feature: {token}.')
    return tokens, [token for token in tokens if token != FeatureName.PLAIN.value]


def reject_experiment_options(options, where='a feature'):
    """Reject experiment options supplied where only stage-local options belong."""
    for key in options:
        if key in EXPERIMENT_OPTIONS:
            example = options[key]
            raise ValueError(f'Option `{key}` is an experiment option, not {where} option; pass it as '
                             f'benchmark(..., {key}={example!r}) instead of Feature(..., {key}=...).')


def reject_flat_stage_options(options):
    """With a structured specification (or a Feature object), no stage option may be a keyword."""
    extra = sorted(key for key in options if key not in EXPERIMENT_OPTIONS)
    if extra:
        raise ValueError('With a structured feature specification, stage options belong inside each entry\'s '
                         f'"options"; only experiment options such as `n_runs` may be given as keywords. '
                         f'Unexpected keyword(s): {extra}.')


def _stage_records(effective_entries):
    records = []
    occurrences = {}
    for position, (name, options) in enumerate(effective_entries):
        occurrence = occurrences.get(name, 0)
        occurrences[name] = occurrence + 1
        records.append(StageRecord(name, occurrence, options, position))
    return tuple(records)


def normalize_shorthand(name, broadcast):
    """
    Normalize ``'a+b+c'`` with flat broadcast options: every supplied option goes
    to every declared token that accepts it, each token validates and defaults
    its own copy, and options accepted by no token are rejected.

    Returns the declaration and the tuple of effective stage records.
    """
    declared_tokens, effective_tokens = parse_feature_name(name)
    supplied = {}
    for key, value in broadcast.items():
        key = key.lower()
        if key not in FeatureOption.__members__.values() and key not in EXPERIMENT_OPTIONS:
            raise ValueError(f'Unknown option for feature: {key}.')
        supplied[key] = value
    reject_experiment_options(supplied)
    entries = []
    effective = []
    for position, token in enumerate(declared_tokens):
        routed = {key: value for key, value in supplied.items() if key in local_options(token)}
        entries.append((token, routed))
        if token == FeatureName.PLAIN.value:
            continue
        try:
            validated = validated_local_options(token, routed)
        except (TypeError, ValueError) as err:
            raise type(err)(f"Invalid options for stage {position + 1} '{token}' of feature "
                            f"'{'+'.join(effective_tokens)}': {err}") from err
        effective.append((token, validated))
    for key in supplied:
        if not any(key in local_options(token) for token in effective_tokens):
            raise ValueError(f"Option `{key}` is not valid for feature '{'+'.join(effective_tokens) or name}'.")
    return Declaration('feature_name', entries), _stage_records(effective)


def normalize_entries(spec):
    """
    Normalize a structured specification: one stage mapping
    ``{'name': ..., 'options': {...}}`` or an ordered list/tuple of such mappings
    and bare names. Every entry, ``plain`` included, is validated before ``plain``
    entries are removed from the effective stages. Entry keys and option names
    are checked for type before any comparison or formatting, so hostile key
    objects never have their hooks run.

    Returns the declaration and the tuple of effective stage records.
    """
    if isinstance(spec, Mapping):
        items = [spec]
    elif isinstance(spec, (list, tuple)):
        items = list(spec)
    else:
        raise TypeError(_SPEC_TYPE_MESSAGE)
    if not items:
        raise ValueError('A structured feature specification must contain at least one stage entry.')
    entries = []
    effective = []
    for position, entry in enumerate(items, start=1):
        label = f'entry {position} of the feature specification'
        if isinstance(entry, str):
            name, options = entry, {}
        elif isinstance(entry, Mapping):
            if any(not isinstance(key, str) for key in entry):
                raise TypeError(f'{label}: entry keys must be strings ("name" and "options").')
            unknown = sorted(key for key in entry if key not in ('name', 'options'))
            if unknown:
                raise ValueError(f'{label}: unknown key(s) {unknown}; a stage entry accepts only "name" and "options".')
            if 'name' not in entry:
                raise ValueError(f'{label}: missing "name".')
            name = entry['name']
            options = entry.get('options', {})
            if not isinstance(options, Mapping):
                raise TypeError(f'{label}: "options" must be a mapping of stage options.')
        else:
            raise TypeError(f'{label}: a stage entry must be a mapping with "name" and optional "options", '
                            f'or a feature name string.')
        if not isinstance(name, str):
            raise TypeError(f'{label}: "name" must be a string.')
        token = name.strip().lower()
        if '+' in token:
            raise ValueError(f'{label}: stage names are atomic; {name!r} declares a composition, '
                             f'give one entry per stage.')
        if token not in FeatureName.__members__.values():
            raise ValueError(f'{label}: unknown feature {name!r}.')
        if any(not isinstance(key, str) for key in options):
            raise TypeError(f"{label} (stage '{token}'): option names must be strings.")
        supplied = {key.lower(): value for key, value in options.items()}
        try:
            validated = validated_local_options(token, supplied)
        except (TypeError, ValueError) as err:
            raise type(err)(f"{label} (stage '{token}'): {err}") from err
        entries.append((token, supplied))
        if token != FeatureName.PLAIN.value:
            effective.append((token, validated))
    return Declaration('feature', entries), _stage_records(effective)
