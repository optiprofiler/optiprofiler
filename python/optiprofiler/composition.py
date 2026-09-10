"""
Ordered composition of features as lightweight problem views.

A composed feature is written as ``feature_name='a+b+c'``. Stage ``a`` is
applied to the original problem first, then ``b``, then ``c``, so the problem
handed to a solver is ``c(b(a(P)))``. Each stage is a lazy *view* of its
predecessor: it exposes the structure of a problem (initial point, bounds,
linear and nonlinear constraints) together with two evaluation channels.

- The *observed* channel is what a solver sees. A query made to the final
  view descends outer-to-inner: an outer value transform first asks its
  predecessor for the value it needs and then modifies it; a point map such
  as ``quantized`` first moves the point and then queries its predecessor
  once at the moved point. No stage evaluates the original problem itself,
  so the number of original callbacks does not grow with the chain length.
- The *reference* channel is the scoring truth. Value transforms (``noisy``,
  ``truncated``, ``random_nan``, ``nonquantifiable_constraints``) leave the
  inherited reference untouched; coordinate transforms transport it; and
  ``quantized`` with ``ground_truth=True`` snaps the inherited reference
  queries as well. Reference reads never advance any random stream.

Exactly one recorder, :class:`ComposedFeaturedProblem`, owns the evaluation
budget, the histories, the cached last values and the termination rule. The
views own no accounting at all. A feature name whose effective pipeline has
at most one stage (``plain`` tokens are removed) does not use this module:
it runs through the established single-feature path of
:class:`optiprofiler.opclasses.FeaturedProblem`, including its seeds, its
callback pattern and its stamp.

Random streams. Every stage of a composition receives its own seed derived
from the run seed with :class:`numpy.random.SeedSequence` and a spawn key
made of the stage's frozen numeric code (:data:`STAGE_CODES`) and its
occurrence index among stages of the same name. Inserting or removing
``plain`` therefore never shifts another stage's stream, and repeated
stages have distinct streams. Within a stage, the served-query counter of
each channel (``fun``, ``cub``, ``ceq``) plays the role that the outer
history length plays for a single feature: it advances on every observed
query the stage serves, including probes made by a later ``custom`` or
``unrelaxable_constraints`` stage, and it is separate from the solver
budget kept by the recorder.
"""

import copy

import numpy as np

from .metadata import safe_metadata
from .opclasses import (Feature, FeatureName, FeatureOption, FeaturedProblem, Problem,
                        _process_1d_array, _validate_max_eval, _validate_seed)
from .utils import get_logger, shorten_log_message

# Frozen numeric stage codes used in seed derivation. Never renumber; append
# new features with new codes. ``plain`` has no code because it never forms a
# stage of an effective pipeline.
STAGE_CODES = {
    'perturbed_x0': 1,
    'noisy': 2,
    'truncated': 3,
    'permuted': 4,
    'linearly_transformed': 5,
    'random_nan': 6,
    'unrelaxable_constraints': 7,
    'nonquantifiable_constraints': 8,
    'quantized': 9,
    'custom': 10,
}

# Identifier of the stage seed derivation recorded in archives and reports.
SEED_POLICY = 'seedsequence-v1'

CHANNELS = ('fun', 'cub', 'ceq')

_DERIVATIVE_MESSAGE = ('Derivatives of a composed feature are not provided: the composed problem '
                       'is not the original problem and its derivatives are not those of the '
                       'original callbacks.')


def parse_feature_name(name):
    """
    Parse a feature name into its declared form and its effective stages.

    Tokens are separated by ``+``, lowercased and stripped of surrounding
    whitespace. Empty or unknown tokens are rejected. ``plain`` tokens are
    kept in the declared name but removed from the effective stage list.

    Returns
    -------
    str
        The normalized declared name, e.g. ``'noisy+plain+truncated'``.
    list of str
        The effective stage names in order, e.g. ``['noisy', 'truncated']``.
    """
    if not isinstance(name, str):
        raise TypeError('The first input argument for `Feature` must be a string.')
    tokens = [token.strip().lower() for token in name.split('+')]
    if any(token == '' for token in tokens):
        raise ValueError(f'Invalid feature name {name!r}: empty stage token in a "+"-separated composition.')
    for token in tokens:
        if token not in FeatureName.__members__.values():
            raise ValueError(f'Unknown feature: {token}.')
    declared = '+'.join(tokens)
    effective = [token for token in tokens if token != FeatureName.PLAIN.value]
    return declared, effective


class Stage:
    """One effective stage of a composition: identity plus its validated child feature."""

    __slots__ = ('position', 'name', 'code', 'occurrence', 'feature')

    def __init__(self, position, name, occurrence, feature):
        self.position = position
        self.name = name
        self.code = STAGE_CODES[name]
        self.occurrence = occurrence
        self.feature = feature

    @property
    def identity(self):
        return f'{self.name}#{self.occurrence}'


def stage_seed(run_seed, stage):
    """Derive the 32-bit seed of a stage from the run seed and the stage identity."""
    entropy = 0 if run_seed is None else int(run_seed)
    sequence = np.random.SeedSequence(entropy, spawn_key=(stage.code, stage.occurrence))
    return int(sequence.generate_state(1, dtype=np.uint32)[0])


class StageContext:
    """Private per-run context of one stage: its seed and served-query counters."""

    __slots__ = ('seed', 'served')

    def __init__(self, seed):
        self.seed = seed
        self.served = {channel: 0 for channel in CHANNELS}

    def next_index(self, channel):
        index = self.served[channel]
        self.served[channel] = index + 1
        return index


def _describe_options(options):
    """Encode effective options without executing user code; never raise."""
    try:
        return safe_metadata(dict(options))
    except Exception as exc:  # defensive: the encoder itself must not abort a benchmark
        return {'value': None, 'reason': 'options_not_described', 'error_type': type(exc).__name__}


def _stored_state(obj):
    """The instance dictionary of ``obj`` without invoking any property or descriptor."""
    try:
        return object.__getattribute__(obj, '__dict__')
    except (AttributeError, TypeError):
        return {}


def describe_pipeline(feature, feature_stamp=None, full_feature_stamp=None):
    """
    Describe the ordered pipeline of a feature as plain, JSON-serializable data.

    Only stored state is read (no property, modifier or callback is invoked),
    and option values are encoded with the shared metadata encoder, which
    describes callables by name only and never invokes user representations.
    A single feature is described as a one-stage pipeline with the legacy
    seed policy; a composition lists its effective stages with their
    identities and effective options. An option set that cannot be encoded is
    recorded as an explicit reason record rather than aborting the benchmark.
    """
    state = _stored_state(feature)
    stages = state.get('_stages')
    if stages:
        entries = [{
            'position': stage.position,
            'name': stage.name,
            'code': stage.code,
            'occurrence': stage.occurrence,
            'identity': stage.identity,
            'options': _describe_options(_stored_state(stage.feature).get('_options', {})),
        } for stage in stages]
        seed_policy = SEED_POLICY
    else:
        name = state.get('_name')
        entries = [{
            'position': 0,
            'name': name,
            'code': STAGE_CODES.get(name),
            'occurrence': 0,
            'identity': f'{name}#0',
            'options': _describe_options(state.get('_options', {})),
        }]
        seed_policy = 'legacy-run-seed'
    return {
        'schema': 'feature_pipeline-v1',
        'declared_name': state.get('_declared_name', state.get('_name')),
        'effective_name': state.get('_name'),
        'seed_policy': seed_policy,
        'feature_stamp': feature_stamp,
        'full_feature_stamp': full_feature_stamp,
        'stages': entries,
    }


class ComposedFeature(Feature):
    """
    A feature with at least two effective stages.

    Instances are created by ``Feature('a+b+c', **options)``; the base class
    dispatches here from ``Feature.__new__``. Supplied options are routed to
    every stage that owns the key, and each stage validates and defaults its
    options independently. ``n_runs`` is global.
    """

    def __init__(self, name, **feature_options):
        # The legacy single-feature constructor is deliberately not called.
        declared, effective = parse_feature_name(name)
        self._declared_name = declared
        self._name = '+'.join(effective)
        options = {key.lower(): value for key, value in feature_options.items()}
        for key in options:
            if key not in FeatureOption.__members__.values():
                raise ValueError(f'Unknown option for feature: {key}.')
        n_runs = options.pop(FeatureOption.N_RUNS.value, None)
        stages = []
        occurrences = {}
        for position, stage_name in enumerate(effective):
            occurrence = occurrences.get(stage_name, 0)
            occurrences[stage_name] = occurrence + 1
            owned = Feature._known_options(stage_name)
            routed = {key: value for key, value in options.items() if key in owned}
            try:
                child = Feature(stage_name, **routed)
            except (TypeError, ValueError) as err:
                raise type(err)(f"Invalid options for stage {position + 1} '{stage_name}' (occurrence "
                                f"{occurrence + 1}) of feature '{self._name}': {err}") from err
            stages.append(Stage(position, stage_name, occurrence, child))
        for key in options:
            if not any(key in Feature._known_options(stage.name) for stage in stages):
                raise ValueError(f"Option `{key}` is not valid for feature '{self._name}'.")
        if n_runs is None:
            n_runs = max(stage.feature.options[FeatureOption.N_RUNS] for stage in stages)
        else:
            n_runs = Feature._validate_n_runs(n_runs)
        self._options = dict(options)
        self._options[FeatureOption.N_RUNS.value] = n_runs
        self._stages = tuple(stages)

    @property
    def is_stochastic(self):
        return any(stage.feature.is_stochastic for stage in self._stages)

    def _unsupported_modifier(self, *args, **kwargs):
        raise NotImplementedError('The modifier methods of a composed feature are not available; '
                                  'apply the composition through `FeaturedProblem(problem, feature, max_eval, seed)`.')

    modifier_x0 = modifier_affine = modifier_bounds = modifier_linear_ub = modifier_linear_eq = _unsupported_modifier
    modifier_fun = modifier_cub = modifier_ceq = _unsupported_modifier


def _structural_violation(view, x):
    """Bound and linear violations of ``view``'s own structure at ``x`` (legacy formulas)."""
    xl, xu = view.xl, view.xu
    cv_bounds = 0.0
    if np.any(np.isfinite(xl)):
        cv_bounds = np.max(xl - x, initial=0.0)
    if np.any(np.isfinite(xu)):
        cv_bounds = np.max(x - xu, initial=cv_bounds)
    aub, bub, aeq, beq = view.aub, view.bub, view.aeq, view.beq
    cv_linear = 0.0
    if aub.size > 0:
        cv_linear = np.max(aub @ x - bub, initial=0.0)
    if aeq.size > 0:
        cv_linear = np.max(np.abs(aeq @ x - beq), initial=cv_linear)
    return cv_bounds, cv_linear


def _nonlinear_violation(cub, ceq, m_ub, m_eq, x):
    """Nonlinear violation from the given constraint evaluators (legacy formulas)."""
    cv_nonlinear = 0.0
    if m_ub > 0:
        cv_nonlinear = np.max(cub(x), initial=0.0)
    if m_eq > 0:
        cv_nonlinear = np.max(np.abs(ceq(x)), initial=cv_nonlinear)
    return cv_nonlinear


class ProblemView(Problem):
    """
    Lazy view of a predecessor problem.

    The base view is the identity: structure and both channels are delegated
    to the predecessor. Subclasses override only what their stage changes.
    ``Problem.__init__`` is deliberately not called: it probes the nonlinear
    constraints at ``x0``, which would be a hidden query. The Problem-facing
    methods ``fun``, ``cub``, ``ceq`` and ``maxcv`` expose the observed
    channel, which is what a custom callback handed this view should see.
    """

    def __init__(self, predecessor, stage=None, context=None):
        self._predecessor = predecessor
        self._stage = stage
        self._feature = stage.feature if stage is not None else None
        self._context = context
        self._name = predecessor._name
        self._x0 = predecessor.x0
        self._xl = predecessor.xl
        self._xu = predecessor.xu
        self._aub = predecessor.aub
        self._bub = predecessor.bub
        self._aeq = predecessor.aeq
        self._beq = predecessor.beq
        self._m_nonlinear_ub = predecessor.m_nonlinear_ub
        self._m_nonlinear_eq = predecessor.m_nonlinear_eq
        self._fun = self._cub = self._ceq = None
        self._grad = self._hess = self._jcub = self._jceq = self._hcub = self._hceq = None

    # Problem-facing API: the observed channel.
    def fun(self, x):
        return self.observed_fun(x)

    def cub(self, x):
        return self.observed_cub(x)

    def ceq(self, x):
        return self.observed_ceq(x)

    def maxcv(self, x):
        return self.observed_maxcv_detailed(x)[0]

    def _maxcv(self, x):
        return self.observed_maxcv_detailed(x)

    def _no_derivatives(self, x):
        raise NotImplementedError(_DERIVATIVE_MESSAGE)

    grad = hess = jcub = jceq = hcub = hceq = _no_derivatives

    # Observed channel (identity by default).
    def observed_fun(self, x):
        return self._predecessor.observed_fun(x)

    def observed_cub(self, x):
        return self._predecessor.observed_cub(x)

    def observed_ceq(self, x):
        return self._predecessor.observed_ceq(x)

    def observed_violation_structural(self, x):
        return _structural_violation(self, x)

    def observed_violation_nonlinear(self, x):
        return _nonlinear_violation(self.observed_cub, self.observed_ceq, self.m_nonlinear_ub, self.m_nonlinear_eq, x)

    def observed_maxcv_detailed(self, x):
        cv_bounds, cv_linear = self.observed_violation_structural(x)
        cv_nonlinear = self.observed_violation_nonlinear(x)
        return np.max([cv_bounds, cv_linear, cv_nonlinear]), cv_bounds, cv_linear, cv_nonlinear

    # Reference channel (inherited by default).
    def reference_fun(self, x):
        return self._predecessor.reference_fun(x)

    def reference_cub(self, x):
        return self._predecessor.reference_cub(x)

    def reference_ceq(self, x):
        return self._predecessor.reference_ceq(x)

    def reference_violation_structural(self, x):
        return self._predecessor.reference_violation_structural(x)

    def reference_violation_nonlinear(self, x):
        return self._predecessor.reference_violation_nonlinear(x)

    def reference_maxcv_detailed(self, x):
        cv_bounds, cv_linear = self.reference_violation_structural(x)
        cv_nonlinear = self.reference_violation_nonlinear(x)
        return np.max([cv_bounds, cv_linear, cv_nonlinear]), cv_bounds, cv_linear, cv_nonlinear

    def reference_maxcv(self, x):
        return self.reference_maxcv_detailed(x)[0]


class RootView(ProblemView):
    """The original problem seen as a view: both channels are its own callbacks."""

    def __init__(self, problem):
        self._problem = problem
        self._predecessor = None
        self._stage = None
        self._feature = None
        self._context = None
        self._name = problem._name
        self._x0 = problem.x0
        self._xl = problem.xl
        self._xu = problem.xu
        self._aub = problem.aub
        self._bub = problem.bub
        self._aeq = problem.aeq
        self._beq = problem.beq
        self._m_nonlinear_ub = problem.m_nonlinear_ub
        self._m_nonlinear_eq = problem.m_nonlinear_eq
        self._fun = self._cub = self._ceq = None
        self._grad = self._hess = self._jcub = self._jceq = self._hcub = self._hceq = None

    def observed_fun(self, x):
        return self._problem.fun(x)

    def observed_cub(self, x):
        return self._problem.cub(x)

    def observed_ceq(self, x):
        return self._problem.ceq(x)

    reference_fun = observed_fun
    reference_cub = observed_cub
    reference_ceq = observed_ceq

    def reference_violation_structural(self, x):
        return _structural_violation(self, x)

    def reference_violation_nonlinear(self, x):
        return _nonlinear_violation(self._problem.cub, self._problem.ceq, self.m_nonlinear_ub, self.m_nonlinear_eq, x)


class NoisyView(ProblemView):
    """Additive/relative/mixed noise on the observed values; reference inherited."""

    def observed_fun(self, x):
        f = self._predecessor.observed_fun(x)
        index = self._context.next_index('fun')
        noise = self._feature._compute_noise(x, self._context.seed, index, f)
        return self._feature._apply_noise(f, noise)

    def _observed_vector(self, channel, values, x):
        index = self._context.next_index(channel)
        if values.size == 0:
            return values
        noise = self._feature._compute_noise(x, self._context.seed, index, values, values.size)
        return self._feature._apply_noise(values, noise)

    def observed_cub(self, x):
        return self._observed_vector('cub', self._predecessor.observed_cub(x), x)

    def observed_ceq(self, x):
        return self._observed_vector('ceq', self._predecessor.observed_ceq(x), x)


class TruncatedView(ProblemView):
    """Rounding of the observed values to significant digits; reference inherited."""

    def observed_fun(self, x):
        f = self._predecessor.observed_fun(x)
        index = self._context.next_index('fun')
        return self._feature._truncate_scalar(f, x, self._context.seed, index)

    def _observed_vector(self, channel, values, x):
        index = self._context.next_index(channel)
        if values.size == 0:
            return values
        return self._feature._truncate_vector(np.array(values, dtype=float), x, self._context.seed, index)

    def observed_cub(self, x):
        return self._observed_vector('cub', self._predecessor.observed_cub(x), x)

    def observed_ceq(self, x):
        return self._observed_vector('ceq', self._predecessor.observed_ceq(x), x)


class RandomNanView(ProblemView):
    """Observed values replaced by NaN at the configured rate; reference inherited."""

    def observed_fun(self, x):
        f = self._predecessor.observed_fun(x)
        index = self._context.next_index('fun')
        return self._feature._random_nan_scalar(f, x, self._context.seed, index)

    def _observed_vector(self, channel, values, x):
        index = self._context.next_index(channel)
        if values.size == 0:
            return values
        return self._feature._random_nan_vector(np.array(values, dtype=float), x, self._context.seed, index)

    def observed_cub(self, x):
        return self._observed_vector('cub', self._predecessor.observed_cub(x), x)

    def observed_ceq(self, x):
        return self._observed_vector('ceq', self._predecessor.observed_ceq(x), x)


class NonquantifiableView(ProblemView):
    """Observed nonlinear constraints reduced to violated/satisfied flags; reference inherited."""

    def observed_cub(self, x):
        values = self._predecessor.observed_cub(x)
        self._context.next_index('cub')
        if values.size == 0:
            return values
        return self._feature._nonquantifiable_cub(np.array(values, dtype=float))

    def observed_ceq(self, x):
        values = self._predecessor.observed_ceq(x)
        self._context.next_index('ceq')
        if values.size == 0:
            return values
        return self._feature._nonquantifiable_ceq(np.array(values, dtype=float))


class PerturbedX0View(ProblemView):
    """A perturbed initial point; objective, constraints and reference are inherited."""

    def __init__(self, predecessor, stage, context):
        super().__init__(predecessor, stage, context)
        self._x0 = self._feature.modifier_x0(context.seed, predecessor)


class AffineView(ProblemView):
    """
    A change of variables ``x_predecessor = A @ x + b`` (``permuted`` and
    ``linearly_transformed``).

    Both channels are evaluated at the mapped point, and the structure is
    transported with the feature's own modifiers: the initial point is pulled
    back through the inverse, finite bounds become linear constraints when
    ``A`` is not diagonal, and linear constraints are composed with ``A``.
    The observed structure is the transported one, which is what a solver
    is handed; the reference violation is measured by the predecessor at the
    mapped point, so scoring stays in the predecessor's coordinates.
    """

    def __init__(self, predecessor, stage, context):
        super().__init__(predecessor, stage, context)
        feature, seed = self._feature, context.seed
        self._A, self._b, self._inv = feature.modifier_affine(seed, predecessor)
        self._x0 = feature.modifier_x0(seed, predecessor)
        self._xl, self._xu = feature.modifier_bounds(seed, predecessor)
        self._aub, self._bub = feature.modifier_linear_ub(seed, predecessor)
        self._aeq, self._beq = feature.modifier_linear_eq(seed, predecessor)

    def _map(self, x):
        return self._A @ x + self._b

    def observed_fun(self, x):
        return self._predecessor.observed_fun(self._map(x))

    def observed_cub(self, x):
        return self._predecessor.observed_cub(self._map(x))

    def observed_ceq(self, x):
        return self._predecessor.observed_ceq(self._map(x))

    def reference_fun(self, x):
        return self._predecessor.reference_fun(self._map(x))

    def reference_cub(self, x):
        return self._predecessor.reference_cub(self._map(x))

    def reference_ceq(self, x):
        return self._predecessor.reference_ceq(self._map(x))

    def reference_violation_structural(self, x):
        return self._predecessor.reference_violation_structural(self._map(x))

    def reference_violation_nonlinear(self, x):
        return self._predecessor.reference_violation_nonlinear(self._map(x))


class QuantizedView(ProblemView):
    """
    Evaluation on a mesh: every observed query is served by the predecessor
    at the snapped point, once per channel. With ``ground_truth=False`` the
    inherited reference is unchanged; with ``ground_truth=True`` the
    reference objective and nonlinear constraints are also read at the
    snapped point. Bound and linear violations are always measured at the
    unsnapped point, and the solver's coordinates are never snapped.
    """

    def __init__(self, predecessor, stage, context):
        super().__init__(predecessor, stage, context)
        self._ground_truth = bool(self._feature.options[FeatureOption.GROUND_TRUTH])

    def _snap(self, x):
        return self._feature._quantize_point(x)

    def _reference_point(self, x):
        return self._snap(x) if self._ground_truth else x

    def observed_fun(self, x):
        return self._predecessor.observed_fun(self._snap(x))

    def observed_cub(self, x):
        return self._predecessor.observed_cub(self._snap(x))

    def observed_ceq(self, x):
        return self._predecessor.observed_ceq(self._snap(x))

    def reference_fun(self, x):
        return self._predecessor.reference_fun(self._reference_point(x))

    def reference_cub(self, x):
        return self._predecessor.reference_cub(self._reference_point(x))

    def reference_ceq(self, x):
        return self._predecessor.reference_ceq(self._reference_point(x))

    def reference_violation_nonlinear(self, x):
        return self._predecessor.reference_violation_nonlinear(self._reference_point(x))


class UnrelaxableView(ProblemView):
    """
    The objective becomes infinite where the immediate predecessor's observed
    constraints of an enabled category are violated. Categories refer to the
    predecessor's own representation: after a rotation, former bounds are
    linear constraints. Violations containing NaN compare as not violated,
    exactly as in the single-feature implementation. Constraint samples drawn
    by the gate are genuine observed queries of the predecessor; they never
    touch the recorder's constraint budget or histories.
    """

    def observed_fun(self, x):
        f = self._predecessor.observed_fun(x)
        options = self._feature.options
        cv_bounds, cv_linear = self._predecessor.observed_violation_structural(x)
        if options[FeatureOption.UNRELAXABLE_BOUNDS] and cv_bounds > 0.0:
            return np.inf
        elif options[FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS] and cv_linear > 0.0:
            return np.inf
        elif options[FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS] \
                and self._predecessor.observed_violation_nonlinear(x) > 0.0:
            return np.inf
        return f


def _stage_label(stage):
    return f"stage {stage.position + 1} '{stage.name}' (occurrence {stage.occurrence + 1})"


def _custom_scalar(value, stage, key):
    """
    ``Problem.fun`` scalar policy for a custom objective output: a real
    scalar or a one-element real array converts to ``float``; anything else,
    including complex values, is logged and recorded as NaN.
    """
    try:
        if isinstance(value, np.ndarray):
            if value.size != 1:
                raise TypeError(f'an array of shape {value.shape} is not a scalar')
            value = value.reshape(-1)[0]
        if isinstance(value, (str, bytes)) or np.iscomplexobj(value):
            raise TypeError(f'{type(value).__name__} is not a real scalar')
        return float(value)
    except Exception as exc:
        get_logger(__name__).warning(
            f'The callback `{key}` of {_stage_label(stage)} returned a value that is not a real scalar '
            f'({shorten_log_message(exc)}); the observed objective value is recorded as NaN.')
        return np.nan


def _custom_vector(value, size, stage, key, what='an array'):
    """A custom callback output as a fresh real one-dimensional float array of the required size."""
    prefix = f'The callback `{key}` of {_stage_label(stage)} returned'
    requirement = f'a real one-dimensional array of size {size} is required'
    try:
        array = np.empty(0) if value is None else np.asarray(value)
    except Exception as exc:
        raise ValueError(f'{prefix} {what} that is not array-like; {requirement}.') from exc
    if np.iscomplexobj(array):
        raise ValueError(f'{prefix} complex values; {requirement}.')
    try:
        array = np.array(array, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{prefix} values that cannot be converted to real numbers; {requirement}.') from exc
    array = np.atleast_1d(np.squeeze(array))
    if array.ndim != 1:
        raise ValueError(f'{prefix} {what} of shape {tuple(np.shape(value))}; {requirement}.')
    if array.size != size:
        raise ValueError(f'{prefix} {what} of size {array.size}; {requirement}.')
    return array


def _custom_matrix(value, shape, stage, key, what='a matrix'):
    """A custom callback output as a fresh real two-dimensional float array of the required shape."""
    prefix = f'The callback `{key}` of {_stage_label(stage)} returned'
    requirement = f'a real matrix of shape {shape} is required'
    try:
        array = np.asarray(value)
    except Exception as exc:
        raise ValueError(f'{prefix} {what} that is not array-like; {requirement}.') from exc
    if np.iscomplexobj(array):
        raise ValueError(f'{prefix} complex values; {requirement}.')
    try:
        array = np.atleast_2d(np.array(array, dtype=float))
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{prefix} values that cannot be converted to real numbers; {requirement}.') from exc
    if array.shape != tuple(shape):
        raise ValueError(f'{prefix} {what} of shape {array.shape}; {requirement}.')
    return array


def _normalized_custom_feature(feature, predecessor, stage):
    """
    A copy of a custom child feature whose user callbacks are wrapped so that
    every output crosses the composition boundary as validated numeric data.
    The legacy modifiers then run unchanged on this copy; the original child
    feature (and its provenance) is untouched. The wrappers are small
    picklable callables, so a composed problem with top-level callbacks can
    still be pickled.
    """
    n = predecessor.n
    options = dict(feature._options)
    wrappers = {
        FeatureOption.MOD_FUN: lambda user: _ScalarCallback(user, stage),
        FeatureOption.MOD_CUB: lambda user: _VectorCallback(user, predecessor.m_nonlinear_ub, stage, 'mod_cub'),
        FeatureOption.MOD_CEQ: lambda user: _VectorCallback(user, predecessor.m_nonlinear_eq, stage, 'mod_ceq'),
        FeatureOption.MOD_X0: lambda user: _InitialPointCallback(user, n, stage),
        FeatureOption.MOD_BOUNDS: lambda user: _BoundsCallback(user, n, stage),
        FeatureOption.MOD_LINEAR_UB: lambda user: _LinearCallback(user, n, stage, 'mod_linear_ub'),
        FeatureOption.MOD_LINEAR_EQ: lambda user: _LinearCallback(user, n, stage, 'mod_linear_eq'),
        FeatureOption.MOD_AFFINE: lambda user: _AffineCallback(user, n, stage),
    }
    for key, wrap in wrappers.items():
        if key in options:
            options[key.value] = wrap(options[key])
    normalized = copy.copy(feature)
    normalized._options = options
    return normalized


def _pair(result, stage, key):
    if not isinstance(result, (tuple, list)) or len(result) != 2:
        raise ValueError(f'The callback `{key}` of {_stage_label(stage)} must return a pair of arrays.')
    return result


class _ScalarCallback:
    """``mod_fun`` wrapper applying the scalar policy."""

    def __init__(self, user, stage):
        self.user, self.stage = user, stage

    def __call__(self, x, rng, problem):
        return _custom_scalar(self.user(x, rng, problem), self.stage, 'mod_fun')


class _VectorCallback:
    """``mod_cub``/``mod_ceq`` wrapper requiring the predecessor's channel size."""

    def __init__(self, user, size, stage, key):
        self.user, self.size, self.stage, self.key = user, size, stage, key

    def __call__(self, x, rng, problem):
        return _custom_vector(self.user(x, rng, problem), self.size, self.stage, self.key)


class _InitialPointCallback:

    def __init__(self, user, n, stage):
        self.user, self.n, self.stage = user, n, stage

    def __call__(self, rng, problem):
        return _custom_vector(self.user(rng, problem), self.n, self.stage, 'mod_x0', 'an initial point')


class _BoundsCallback:

    def __init__(self, user, n, stage):
        self.user, self.n, self.stage = user, n, stage

    def __call__(self, rng, problem):
        xl, xu = _pair(self.user(rng, problem), self.stage, 'mod_bounds')
        return (_custom_vector(xl, self.n, self.stage, 'mod_bounds', 'lower bounds'),
                _custom_vector(xu, self.n, self.stage, 'mod_bounds', 'upper bounds'))


class _LinearCallback:
    """``mod_linear_ub``/``mod_linear_eq`` wrapper: a real (m, n) matrix and a size-m right-hand side."""

    def __init__(self, user, n, stage, key):
        self.user, self.n, self.stage, self.key = user, n, stage, key

    def __call__(self, rng, problem):
        matrix, rhs = _pair(self.user(rng, problem), self.stage, self.key)
        prefix = f'The callback `{self.key}` of {_stage_label(self.stage)} returned'
        try:
            probe = np.asarray(matrix, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(f'{prefix} a coefficient matrix whose values cannot be converted to real numbers; '
                             f'a real matrix with {self.n} columns is required.') from exc
        rows = 0 if probe.size == 0 else np.atleast_2d(probe).shape[0]
        if rows == 0:
            matrix = np.empty((0, self.n))
        else:
            matrix = _custom_matrix(matrix, (rows, self.n), self.stage, self.key, 'a coefficient matrix')
        rhs = _custom_vector(rhs, rows, self.stage, self.key, 'a right-hand side')
        return matrix, rhs


class _AffineCallback:

    def __init__(self, user, n, stage):
        self.user, self.n, self.stage = user, n, stage

    def __call__(self, rng, problem):
        result = self.user(rng, problem)
        if not isinstance(result, (tuple, list)) or len(result) != 3:
            raise ValueError(f'The callback `mod_affine` of {_stage_label(self.stage)} must return a matrix, '
                             f'a vector and an inverse matrix.')
        return (_custom_matrix(result[0], (self.n, self.n), self.stage, 'mod_affine'),
                _custom_vector(result[1], self.n, self.stage, 'mod_affine', 'a translation vector'),
                _custom_matrix(result[2], (self.n, self.n), self.stage, 'mod_affine', 'an inverse matrix'))


class CustomView(AffineView):
    """
    User-supplied modifiers at any position of a composition. The callbacks
    receive the immediate predecessor view as their ``problem`` argument, so
    ``problem.fun(x)`` is a genuine observed query of that predecessor, and a
    stochastic predecessor draws a fresh sample for each probe. As in the
    single-feature implementation, the objective and constraint callbacks
    only change observations, ``mod_affine`` transports both channels, and
    the stream handed to a value callback depends on the value read first.

    Callback outputs are the only untrusted values entering a composition,
    so they are normalized here, once: the objective follows the
    ``Problem.fun`` scalar policy (NaN with a logged warning when the output
    is not a real scalar), nonlinear constraint outputs must be real
    one-dimensional arrays of the predecessor's channel size, and
    construction outputs must have the shapes the problem structure requires.
    Violations raise ``ValueError`` naming the stage, occurrence and callback.
    """

    def __init__(self, predecessor, stage, context):
        stage_with_normalized_feature = copy.copy(stage)
        stage_with_normalized_feature.feature = _normalized_custom_feature(stage.feature, predecessor, stage)
        super().__init__(predecessor, stage_with_normalized_feature, context)
        self._stage = stage

    def observed_fun(self, x):
        xm = self._map(x)
        f = self._predecessor.observed_fun(xm)
        index = self._context.next_index('fun')
        options = self._feature.options
        if FeatureOption.MOD_FUN in options:
            rng_custom = Feature.get_default_rng(self._context.seed, f, *xm, index)
            return options[FeatureOption.MOD_FUN](xm, rng_custom, self._predecessor)
        return f

    def _observed_vector(self, channel, key, x):
        xm = self._map(x)
        values = getattr(self._predecessor, f'observed_{channel}')(xm)
        index = self._context.next_index(channel)
        options = self._feature.options
        if values.size == 0 or key not in options:
            return values
        rng_custom = Feature.get_default_rng(self._context.seed, *values, *xm, index)
        return options[key](xm, rng_custom, self._predecessor)

    def observed_cub(self, x):
        return self._observed_vector('cub', FeatureOption.MOD_CUB, x)

    def observed_ceq(self, x):
        return self._observed_vector('ceq', FeatureOption.MOD_CEQ, x)


_VIEW_CLASSES = {
    'perturbed_x0': PerturbedX0View,
    'noisy': NoisyView,
    'truncated': TruncatedView,
    'permuted': AffineView,
    'linearly_transformed': AffineView,
    'random_nan': RandomNanView,
    'unrelaxable_constraints': UnrelaxableView,
    'nonquantifiable_constraints': NonquantifiableView,
    'quantized': QuantizedView,
    'custom': CustomView,
}


def build_view(predecessor, stage, context):
    """Create the view of ``stage`` over ``predecessor``."""
    try:
        view_class = _VIEW_CLASSES[stage.name]
    except KeyError:
        raise NotImplementedError(f"Stage '{stage.name}' is not yet supported in a composition.") from None
    return view_class(predecessor, stage, context)


class ComposedFeaturedProblem(FeaturedProblem):
    """
    The single recorder around the final view of a composition.

    Created by ``FeaturedProblem(problem, feature, max_eval, seed)`` when the
    feature has at least two effective stages. Budget, histories, cached last
    values and termination follow the single-feature wrapper exactly; the
    observed value comes from the final view's observed channel and the
    recorded history from its reference channel. No original callback is
    stored on the recorder, so nothing can bypass an earlier stage.
    """

    def __init__(self, problem, feature, max_eval, seed=None):
        # The legacy wrapper constructor is deliberately not called.
        self._problem = problem
        self._feature = feature
        self._max_eval = _validate_max_eval(max_eval)
        self._seed = _validate_seed(seed)
        self._real_n_eval_fun = 0
        self._real_n_eval_cub = 0
        self._real_n_eval_ceq = 0
        self._fun_hist = []
        self._cub_hist = []
        self._ceq_hist = []
        self._maxcv_hist = []
        self._last_fun = np.nan
        self._last_cub = np.nan
        self._last_ceq = np.nan

        self._contexts = tuple(StageContext(stage_seed(self._seed, stage)) for stage in feature._stages)
        view = RootView(problem)
        views = []
        for stage, context in zip(feature._stages, self._contexts):
            view = build_view(view, stage, context)
            views.append(view)
        self._views = tuple(views)
        self._final = view

        # Structure of the final view, in solver coordinates.
        self._name = view._name
        self._x0 = view.x0
        self._xl = view.xl
        self._xu = view.xu
        self._aub = view.aub
        self._bub = view.bub
        self._aeq = view.aeq
        self._beq = view.beq
        self._m_nonlinear_ub = view.m_nonlinear_ub
        self._m_nonlinear_eq = view.m_nonlinear_eq
        self._fun = self._cub = self._ceq = None
        self._grad = self._hess = self._jcub = self._jceq = self._hcub = self._hceq = None

        self._fun_init, self._maxcv_init = self._evaluate_truth(self._x0)

    def __getnewargs__(self):
        # Pickling support: ``FeaturedProblem.__new__`` requires the constructor
        # arguments, so hand them to the unpickler; the instance state (views,
        # contexts, histories) is then restored from the instance dictionary.
        return self._problem, self._feature, self._max_eval, self._seed

    def _point(self, x, method):
        x = _process_1d_array(x, f'The argument `x` for method `{method}` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `{method}` in problem must have size {self.n}.')
        return x

    def _evaluate_truth(self, x):
        """Scoring reference at solver coordinates, without an observed query."""
        x = self._point(x, 'fun')
        return self._final.reference_fun(x), self._final.reference_maxcv(x)

    def fun(self, x):
        if self._real_n_eval_fun >= 2 * self._max_eval:
            raise StopIteration(f'The number of the objective function evaluations has reached {2 * self._max_eval} (two times the maximum function evaluations).')
        self._real_n_eval_fun += 1
        if self.n_eval_fun >= self._max_eval:
            return self._last_fun
        x = self._point(x, 'fun')
        f = self._final.observed_fun(x)
        self._last_fun = f
        self._fun_hist.append(self._final.reference_fun(x))
        try:
            self._maxcv_hist.append(self._final.reference_maxcv(x))
        except Exception:
            self._maxcv_hist.append(np.nan)
        return f

    def cub(self, x, record_hist=True):
        if self._real_n_eval_cub >= 2 * self._max_eval:
            raise StopIteration(f'The number of the nonlinear inequality constraint evaluations has reached {2 * self._max_eval} (two times the maximum function evaluations).')
        self._real_n_eval_cub += 1
        if self.n_eval_cub >= self._max_eval:
            return self._last_cub
        x = self._point(x, 'cub')
        c = self._final.observed_cub(x)
        self._last_cub = c
        if record_hist:
            self._cub_hist.append(self._final.reference_cub(x))
        return c

    def ceq(self, x, record_hist=True):
        if self._real_n_eval_ceq >= 2 * self._max_eval:
            raise StopIteration(f'The number of the nonlinear equality constraint evaluations has reached {2 * self._max_eval} (two times the maximum function evaluations).')
        self._real_n_eval_ceq += 1
        if self.n_eval_ceq >= self._max_eval:
            return self._last_ceq
        x = self._point(x, 'ceq')
        c = self._final.observed_ceq(x)
        self._last_ceq = c
        if record_hist:
            self._ceq_hist.append(self._final.reference_ceq(x))
        return c

    def maxcv(self, x):
        return self._final.reference_maxcv(self._point(x, 'maxcv'))

    def _no_derivatives(self, x):
        raise NotImplementedError(_DERIVATIVE_MESSAGE)

    grad = hess = jcub = jceq = hcub = hceq = _no_derivatives
