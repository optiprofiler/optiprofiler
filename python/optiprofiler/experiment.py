"""
Experiment planning: how many times a featured experiment is repeated, why,
and which execution strategy runs it.

A feature specification never stores a run count. The experiment resolves one
plan per role (``primary`` for the benchmarked feature, ``plain_reference`` for
the optional plain baseline) once the solver metadata is known, and the plan
travels explicitly to aggregation, to every worker and into provenance.

Stored run slots, actual solver executions and deterministic copies stay
distinct: ``n_runs`` is the size of the retained run axis; ``actual_runs``
gives the executions per solver, and the remaining slots are copies of the
single deterministic run (reported as repeated, never as actual).
"""

import numpy as np

from .feature_definitions import EXPERIMENT_OPTIONS  # noqa: F401  (re-exported for the option partition)

#: Identifier of the literal 1.x default-run policy implemented by ``resolve_plan``.
RUN_POLICY = 'legacy-hints-v1'
# Language-local identifier of the runtime that executes the plans of this
# producer (MATLAB records its own); stated once per role, never per run.
RUNTIME_POLICY = 'python-featured-problem-v1'

#: Roles of an experiment execution.
PRIMARY = 'primary'
PLAIN_REFERENCE = 'plain_reference'

#: Execution strategies of the recorder (explicit, recorded in provenance).
STRATEGY_IDENTITY = 'identity'
STRATEGY_SINGLE = 'legacy-single'
STRATEGY_COMPOSED = 'composed-views'

_ABSENT = object()
#: Sentinel meaning "the caller supplied no run count" (distinct from an explicit ``None``).
ABSENT = _ABSENT


def validate_n_runs(value):
    """Validate an explicitly supplied run count (a present ``None`` is a value, not an omission)."""
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        value = int(value)
    if isinstance(value, np.integer):
        value = int(value)
    if not isinstance(value, int):
        raise TypeError('Option `n_runs` must be an integer.')
    if value <= 0:
        raise ValueError('Option `n_runs` must be positive.')
    return value


def select_execution_strategy(feature):
    """The recorder strategy for a specification, chosen by its effective stage count."""
    count = len(feature.stages)
    if count == 0:
        return STRATEGY_IDENTITY
    if count == 1:
        return STRATEGY_SINGLE
    return STRATEGY_COMPOSED


class ExperimentPlan:
    """Resolved repetitions of one execution role; immutable and picklable."""

    __slots__ = ('_role', '_n_runs', '_origin', '_execution_strategy')

    def __init__(self, role, n_runs, origin, execution_strategy):
        object.__setattr__(self, '_role', role)
        object.__setattr__(self, '_n_runs', int(n_runs))
        object.__setattr__(self, '_origin', origin)
        object.__setattr__(self, '_execution_strategy', execution_strategy)

    def __setattr__(self, key, value):
        raise AttributeError('ExperimentPlan is immutable.')

    def __getstate__(self):
        return (self._role, self._n_runs, self._origin, self._execution_strategy)

    def __setstate__(self, state):
        role, n_runs, origin, strategy = state
        object.__setattr__(self, '_role', role)
        object.__setattr__(self, '_n_runs', int(n_runs))
        object.__setattr__(self, '_origin', origin)
        object.__setattr__(self, '_execution_strategy', strategy)

    role = property(lambda self: self._role)
    n_runs = property(lambda self: self._n_runs)
    origin = property(lambda self: self._origin)
    run_policy = property(lambda self: RUN_POLICY)
    execution_strategy = property(lambda self: self._execution_strategy)

    def actual_runs(self, feature, solver_isrand, n_solvers):
        """Actual executions per solver: the stored count for stochastic features or randomized solvers, else one."""
        return np.array([self._n_runs if feature.is_stochastic or (solver_isrand is not None and solver_isrand[i]) else 1
                         for i in range(n_solvers)], dtype=int)

    def describe(self):
        """Plain data for provenance."""
        return {'role': self._role, 'n_runs': self._n_runs, 'origin': self._origin,
                'run_policy': RUN_POLICY, 'execution_strategy': self._execution_strategy,
                'runtime_policy': RUNTIME_POLICY}

    def __repr__(self):
        return f'ExperimentPlan({self._role!r}, n_runs={self._n_runs}, origin={self._origin!r})'


def resolve_plan(feature, role=PRIMARY, requested=_ABSENT, solver_isrand=None, is_load=False):
    """
    Resolve the run count of one role.

    ``primary``: an explicitly supplied value wins (validated, ``None`` is
    invalid); otherwise five runs when any solver is declared randomized and the
    experiment is not a load; otherwise the largest literal replicate hint of
    the effective stages (one for an identity pipeline). ``plain_reference``:
    always one stored and one actual run per solver (the established reference
    policy); an explicit count is not accepted for that role.
    """
    strategy = select_execution_strategy(feature)
    if role == PLAIN_REFERENCE:
        if requested is not _ABSENT:
            raise ValueError('The plain reference has a fixed run count of 1; it does not accept `n_runs`.')
        return ExperimentPlan(role, 1, 'reference_policy', strategy)
    if role != PRIMARY:
        raise ValueError(f'Unknown experiment role: {role!r}.')
    if requested is not _ABSENT:
        return ExperimentPlan(role, validate_n_runs(requested), 'explicit', strategy)
    if not is_load and solver_isrand is not None and any(solver_isrand):
        return ExperimentPlan(role, 5, 'randomized_solvers', strategy)
    hints = [stage.replicate_hint for stage in feature.stages]
    return ExperimentPlan(role, max(hints) if hints else 1, 'stage_hints', strategy)
