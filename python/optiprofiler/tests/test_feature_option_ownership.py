"""
Ownership of the run count against the options each stage owns.

``n_runs`` belongs to the experiment: the experiment layer resolves one plan per
role from the feature specification and the solver metadata, and the stages of
a specification own only their validated local options. Every expectation is a
worked value read through public construction, ``optiprofiler.experiment`` and
the provenance description.
"""

import pickle

import numpy as np
import pytest

from optiprofiler.experiment import PLAIN_REFERENCE, PRIMARY, resolve_plan
from optiprofiler.opclasses import Feature, FeaturedProblem, Problem
from optiprofiler.provenance import describe_feature

NOISY_DEFAULTS = {'noise_mode': 'random', 'distribution': 'gaussian', 'noise_map': 'chebyshev',
                  'noise_level': 0.001, 'noise_type': 'mixed'}
TRUNCATED_DEFAULTS = {'perturbed_trailing_digits': False, 'significant_digits': 6}


def mod_fun_plus_one(x, rng, problem):
    return problem.fun(x) + 1.0


def sphere(x):
    return float(np.dot(x, x))


def stage_options(feature):
    return [stage['options'] for stage in describe_feature(feature)['stages']]


class TestStageLocalOwnership:

    def test_composite_stage_records_hold_local_options_only(self):
        feature = Feature('noisy+truncated')
        assert [stage.identity for stage in feature.stages] == ['noisy#0', 'truncated#0']
        assert [dict(stage.options) for stage in feature.stages] == [NOISY_DEFAULTS, TRUNCATED_DEFAULTS]
        description = describe_feature(feature)
        assert 'common_options' not in description and 'options' not in description
        assert stage_options(feature) == [NOISY_DEFAULTS, TRUNCATED_DEFAULTS]

    def test_single_feature_is_described_like_any_pipeline(self):
        feature = Feature('noisy', noise_level=0.5)
        description = describe_feature(feature)
        assert description['seed_policy'] == 'legacy-run-seed'
        assert description['effective_name'] == 'noisy' and description['declaration_route'] == 'feature_name'
        assert description['declared'] == [{'name': 'noisy', 'options': {'noise_level': 0.5}}]
        assert stage_options(feature) == [{**NOISY_DEFAULTS, 'noise_level': 0.5}]

    def test_identity_has_no_stage_to_describe(self):
        description = describe_feature(Feature('plain+plain'))
        assert description['stages'] == [] and description['effective_name'] == 'plain'
        assert description['declared'] == [{'name': 'plain', 'options': {}}, {'name': 'plain', 'options': {}}]

    def test_returned_mappings_are_copies(self):
        feature = Feature('noisy+truncated')
        copied = dict(feature.stages[0].options)
        copied['noise_level'] = 99.0
        described = describe_feature(feature)
        described['stages'][0]['options']['noise_level'] = 99.0
        assert feature.stages[0].options['noise_level'] == 0.001
        assert describe_feature(feature)['stages'][0]['options']['noise_level'] == 0.001

    def test_fresh_pickle_round_trip_keeps_ownership(self):
        clone = pickle.loads(pickle.dumps(Feature('noisy+truncated')))
        assert stage_options(clone) == [NOISY_DEFAULTS, TRUNCATED_DEFAULTS]
        assert clone.declared.entries == (('noisy', {}), ('truncated', {}))


class TestExperimentPlans:

    @pytest.mark.parametrize('name, options, expected_runs, stochastic', [
        # custom keeps its established single run although it is stochastic.
        ('custom+truncated', {'mod_fun': mod_fun_plus_one}, 1, True),
        # An unrotated linear transformation keeps five runs although it is deterministic.
        ('linearly_transformed+quantized', {'rotated': False}, 5, False),
        ('noisy+truncated', {}, 5, True),
        ('truncated+quantized', {}, 1, False),
        ('truncated+quantized', {'perturbed_trailing_digits': True}, 5, True),
        ('noisy+truncated', {'noise_mode': 'deterministic'}, 1, False),
        ('noisy+perturbed_x0', {'noise_mode': 'deterministic'}, 5, True),
        ('plain+noisy+plain', {}, 5, True),
        ('unrelaxable_constraints+nonquantifiable_constraints', {}, 1, False),
        ('plain', {}, 1, False),
        ('custom', {'mod_fun': mod_fun_plus_one}, 1, True),
        ('linearly_transformed', {'rotated': False}, 5, False),
    ])
    def test_default_count_is_the_largest_established_stage_hint(self, name, options, expected_runs, stochastic):
        feature = Feature(name, **options)
        plan = resolve_plan(feature)
        assert (plan.role, plan.n_runs, plan.origin, plan.run_policy) == (PRIMARY, expected_runs, 'stage_hints', 'legacy-hints-v1')
        assert feature.is_stochastic is stochastic
        assert all('n_runs' not in options_of_stage for options_of_stage in stage_options(feature))

    def test_explicit_count_and_randomized_solver_rule(self):
        feature = Feature('noisy+perturbed_x0+random_nan')
        explicit = resolve_plan(feature, requested=2)
        assert (explicit.n_runs, explicit.origin) == (2, 'explicit')
        assert resolve_plan(feature, requested=2.0).n_runs == 2
        randomized = resolve_plan(Feature('truncated+quantized'), solver_isrand=[False, True])
        assert (randomized.n_runs, randomized.origin) == (5, 'randomized_solvers')
        # An explicit value beats the randomized-solver rule; a load never applies that rule.
        assert resolve_plan(Feature('truncated'), requested=3, solver_isrand=[True]).n_runs == 3
        loaded = resolve_plan(Feature('truncated'), solver_isrand=[True], is_load=True)
        assert (loaded.n_runs, loaded.origin) == (1, 'stage_hints')
        with pytest.raises(TypeError, match='must be an integer'):
            resolve_plan(feature, requested=None)
        with pytest.raises(ValueError, match='must be positive'):
            resolve_plan(feature, requested=0)

    def test_plain_reference_has_its_own_fixed_plan(self):
        reference = resolve_plan(Feature('plain'), PLAIN_REFERENCE)
        assert (reference.role, reference.n_runs, reference.origin) == (PLAIN_REFERENCE, 1, 'reference_policy')
        assert reference.execution_strategy == 'identity'
        with pytest.raises(ValueError, match='fixed run count'):
            resolve_plan(Feature('plain'), PLAIN_REFERENCE, requested=3)
        with pytest.raises(ValueError, match='Unknown experiment role'):
            resolve_plan(Feature('plain'), 'secondary')

    def test_execution_strategy_is_explicit(self):
        assert resolve_plan(Feature('plain')).execution_strategy == 'identity'
        assert resolve_plan(Feature('noisy')).execution_strategy == 'legacy-single'
        assert resolve_plan(Feature('plain+noisy')).execution_strategy == 'legacy-single'
        assert resolve_plan(Feature('noisy+truncated')).execution_strategy == 'composed-views'
        problem = Problem(sphere, np.array([1.0, 2.0]))
        assert FeaturedProblem(problem, Feature('noisy'), 5, 1).execution_strategy == 'legacy-single'
        assert FeaturedProblem(problem, Feature('noisy+truncated'), 5, 1).execution_strategy == 'composed-views'

    def test_actual_runs_distinguish_stored_slots_from_executions(self):
        deterministic = Feature('truncated+quantized')
        plan = resolve_plan(deterministic, requested=3)
        np.testing.assert_array_equal(plan.actual_runs(deterministic, None, 2), [1, 1])
        np.testing.assert_array_equal(plan.actual_runs(deterministic, [True, False], 2), [3, 1])
        stochastic = Feature('noisy+truncated')
        np.testing.assert_array_equal(resolve_plan(stochastic, requested=3).actual_runs(stochastic, None, 2), [3, 3])

    def test_plans_pickle_for_worker_transport(self):
        plan = resolve_plan(Feature('noisy'), requested=4, solver_isrand=[True])
        clone = pickle.loads(pickle.dumps(plan))
        assert clone.describe() == plan.describe() == {'role': 'primary', 'n_runs': 4, 'origin': 'explicit',
                                                        'run_policy': 'legacy-hints-v1',
                                                        'execution_strategy': 'legacy-single',
                                                        'runtime_policy': 'python-featured-problem-v1'}
        with pytest.raises(AttributeError):
            clone.n_runs = 2
