"""
The structured feature specification: ``Feature(spec, n_runs=...)`` where
``spec`` is a stage mapping ``{'name': ..., 'options': {...}}`` or an ordered
list/tuple of such mappings and bare names.

Stage options live inside each entry; the run count stays experiment-wide.
A specification and the ``'a+b+c'`` shorthand with equal settings build the
same feature and produce the same numbers. Every expectation is a worked
value or a literal comparison with the established shorthand.
"""

import json
import pickle

import numpy as np
import pytest

from optiprofiler.composition import ComposedFeature, describe_pipeline
from optiprofiler.opclasses import Feature, FeaturedProblem, Problem

ALL_TEN = ['plain', 'perturbed_x0', 'noisy', 'truncated', 'permuted', 'linearly_transformed', 'random_nan',
           'unrelaxable_constraints', 'nonquantifiable_constraints', 'quantized']


def sphere(x):
    return float(np.dot(x, x))


def constant_five(x):
    return 5.0


def cub_offset(x):
    return np.asarray(x, dtype=float) - 1.0


def ceq_sum(x):
    return np.array([float(np.sum(x)) - 1.0])


def ones_distribution(rng, size):
    return np.ones(size)


def mod_fun_plus_one(x, rng, problem):
    return problem.fun(x) + 1.0


def constrained_problem():
    return Problem(sphere, np.array([0.5, -1.5, 2.0]), xl=np.array([-3.0, -3.0, -3.0]), xu=np.array([3.0, 3.0, 3.0]),
                   cub=cub_offset, ceq=ceq_sum)


def stage_options(feature):
    return [stage['options'] for stage in describe_pipeline(feature)['stages']]


class TestConstruction:

    def test_single_stage_mapping_builds_the_legacy_feature(self):
        structured = Feature({'name': 'noisy', 'options': {'noise_level': 1e-2}})
        shorthand = Feature('noisy', noise_level=1e-2)
        assert type(structured) is Feature and structured.name == 'noisy'
        assert structured.options == shorthand.options
        assert structured.is_stochastic is True
        pipeline = describe_pipeline(structured)
        assert pipeline['route'] == 'feature'
        assert pipeline['seed_policy'] == 'legacy-run-seed'
        assert pipeline['declared_name'] == 'noisy'
        assert pipeline['declared_spec'] == [{'name': 'noisy', 'options': {'noise_level': 0.01}}]
        assert pipeline['common_options'] == {'n_runs': 5}

    def test_sequence_builds_the_same_composition_as_the_shorthand(self):
        structured = Feature(['noisy', 'truncated'])
        shorthand = Feature('noisy+truncated')
        assert isinstance(structured, ComposedFeature)
        assert structured.name == shorthand.name == 'noisy+truncated'
        assert structured.options == shorthand.options == {'n_runs': 5}
        assert stage_options(structured) == stage_options(shorthand)
        pipeline = describe_pipeline(structured)
        assert pipeline['route'] == 'feature' and pipeline['seed_policy'] == 'seedsequence-v2'
        assert pipeline['declared_spec'] == [{'name': 'noisy', 'options': {}}, {'name': 'truncated', 'options': {}}]
        # The shorthand's declared specification is projected from its tokens
        # and the supplied broadcast options (none here).
        assert describe_pipeline(shorthand)['route'] == 'feature_name'
        assert describe_pipeline(shorthand)['declared_spec'] == pipeline['declared_spec']

    def test_shorthand_declared_specification_projects_supplied_options_to_owning_tokens(self):
        composite = Feature('noisy+plain+perturbed_x0', distribution='gaussian', noise_level=0.5, n_runs=2)
        assert describe_pipeline(composite)['declared_spec'] == [
            {'name': 'noisy', 'options': {'distribution': 'gaussian', 'noise_level': 0.5}},
            {'name': 'plain', 'options': {}},
            {'name': 'perturbed_x0', 'options': {'distribution': 'gaussian'}}]
        # Declared options are the supplied ones, before defaults; the
        # effective options of the stage carry the defaults.
        single = Feature('plain+truncated', significant_digits=4)
        pipeline = describe_pipeline(single)
        assert pipeline['declared_spec'] == [{'name': 'plain', 'options': {}},
                                             {'name': 'truncated', 'options': {'significant_digits': 4}}]
        assert pipeline['stages'][0]['options'] == {'perturbed_trailing_digits': False, 'significant_digits': 4}
        assert describe_pipeline(Feature('noisy'))['declared_spec'] == [{'name': 'noisy', 'options': {}}]
        assert describe_pipeline(Feature('noisy', n_runs=2))['declared_spec'] == [{'name': 'noisy', 'options': {}}]
        assert describe_pipeline(Feature('plain'))['declared_spec'] == [{'name': 'plain', 'options': {}}]

    def test_names_are_normalized_and_tuples_are_accepted(self):
        feature = Feature(({'name': ' Noisy '}, 'TRUNCATED'))
        assert feature.name == 'noisy+truncated'
        assert [stage['identity'] for stage in describe_pipeline(feature)['stages']] == ['noisy#0', 'truncated#0']

    def test_plain_entries_are_validated_then_removed(self):
        feature = Feature(['plain', {'name': 'truncated', 'options': {'significant_digits': 4}}, 'plain'])
        assert type(feature) is Feature and feature.name == 'truncated'
        assert feature.options == Feature('truncated', significant_digits=4).options
        pipeline = describe_pipeline(feature)
        assert pipeline['declared_name'] == 'plain+truncated+plain'
        assert pipeline['effective_name'] == 'truncated'
        assert pipeline['seed_policy'] == 'legacy-run-seed'
        assert [entry['name'] for entry in pipeline['declared_spec']] == ['plain', 'truncated', 'plain']
        assert Feature(['plain', 'plain']).name == 'plain'
        assert Feature({'name': 'plain'}).options == {'n_runs': 1}

    def test_repeated_stages_own_independent_options(self):
        feature = Feature([{'name': 'noisy', 'options': {'noise_level': 1e-3, 'noise_type': 'absolute', 'distribution': ones_distribution}},
                           {'name': 'noisy', 'options': {'noise_level': 1e-1, 'noise_type': 'absolute', 'distribution': ones_distribution}}])
        assert [stage['identity'] for stage in describe_pipeline(feature)['stages']] == ['noisy#0', 'noisy#1']
        assert [options['noise_level'] for options in stage_options(feature)] == [1e-3, 1e-1]
        problem = Problem(constant_five, np.array([1.0, 2.0]), cub=lambda x: np.array([1.0, 2.0]))
        featured = FeaturedProblem(problem, feature, 10, seed=0)
        # Stage noisy#0 adds 1e-3 * 1, then noisy#1 adds 1e-1 * 1, in this order.
        assert featured.fun(problem.x0) == (5.0 + 1e-3 * 1.0) + 1e-1 * 1.0
        np.testing.assert_array_equal(featured.cub(problem.x0), np.array([(1.0 + 1e-3) + 1e-1, (2.0 + 1e-3) + 1e-1]))

    def test_shared_key_can_differ_between_stages(self):
        feature = Feature([{'name': 'noisy', 'options': {'distribution': 'uniform'}},
                           {'name': 'perturbed_x0', 'options': {'distribution': 'gaussian'}}])
        assert [options['distribution'] for options in stage_options(feature)] == ['uniform', 'gaussian']
        assert feature.options == {'n_runs': 5}
        with pytest.raises(ValueError, match='perturbed_x0'):
            Feature('noisy+perturbed_x0', distribution='uniform')

    def test_run_count_is_experiment_wide(self):
        assert Feature(['noisy', 'truncated'], n_runs=3).options == {'n_runs': 3}
        assert Feature(['truncated', 'noisy']).options['n_runs'] == 5
        assert Feature([{'name': 'custom', 'options': {'mod_fun': mod_fun_plus_one}}, 'truncated']).options['n_runs'] == 1
        assert Feature([{'name': 'linearly_transformed', 'options': {'rotated': False}}, 'quantized']).options['n_runs'] == 5
        assert Feature({'name': 'noisy'}, n_runs=2).options['n_runs'] == 2
        assert Feature(['noisy', 'truncated'], n_runs=2.0).options['n_runs'] == 2
        for stage in describe_pipeline(Feature(['noisy', 'truncated'], n_runs=3))['stages']:
            assert 'n_runs' not in stage['options']

    def test_pickle_round_trip(self):
        problem = Problem(sphere, np.array([1.0, -2.0]))
        x = np.array([0.4, 0.6])
        for spec in (['noisy', {'name': 'quantized', 'options': {'mesh_size': 0.5}}],
                     ['plain', {'name': 'noisy', 'options': {'noise_level': 0.25}}]):
            feature = Feature(spec, n_runs=2)
            clone = pickle.loads(pickle.dumps(feature))
            assert clone.name == feature.name and clone.options == feature.options
            assert describe_pipeline(clone) == describe_pipeline(feature)
            assert FeaturedProblem(problem, clone, 5, 3).fun(x) == FeaturedProblem(problem, feature, 5, 3).fun(x)
        # A composed recorder built from a specification pickles like one built from the shorthand.
        featured = FeaturedProblem(problem, Feature(['noisy', {'name': 'quantized', 'options': {'mesh_size': 0.5}}]), 5, 3)
        restored = pickle.loads(pickle.dumps(featured))
        assert restored.fun(x) == FeaturedProblem(problem, Feature('noisy+quantized', mesh_size=0.5), 5, 3).fun(x)

    def test_provenance_is_json_serializable_with_callables_described(self):
        feature = Feature([{'name': 'noisy', 'options': {'distribution': ones_distribution}}, 'truncated'])
        pipeline = describe_pipeline(feature)
        text = json.dumps(pipeline, sort_keys=True)
        assert 'ones_distribution' in text
        assert pipeline['declared_spec'][0]['options']['distribution'] != ones_distribution


class TestEquivalenceWithShorthand:

    @pytest.mark.parametrize('spec, name, flat', [
        (['noisy', 'truncated'], 'noisy+truncated', {}),
        ([{'name': 'noisy', 'options': {'noise_level': 1e-2}}, 'quantized'], 'noisy+quantized', {'noise_level': 1e-2}),
        (['plain', {'name': 'truncated', 'options': {'significant_digits': 4}}, 'plain'], 'plain+truncated+plain',
         {'significant_digits': 4}),
        ([{'name': 'linearly_transformed', 'options': {'rotated': False}}, 'permuted'], 'linearly_transformed+permuted',
         {'rotated': False}),
        (['perturbed_x0', {'name': 'unrelaxable_constraints', 'options': {'unrelaxable_nonlinear_constraints': True}},
          'random_nan'], 'perturbed_x0+unrelaxable_constraints+random_nan', {'unrelaxable_nonlinear_constraints': True}),
        (list(ALL_TEN), '+'.join(ALL_TEN), {}),
        (list(reversed(ALL_TEN)), '+'.join(reversed(ALL_TEN)), {}),
    ])
    @pytest.mark.parametrize('seed', [0, 3])
    def test_same_settings_give_the_same_numbers(self, spec, name, flat, seed):
        structured = Feature(spec)
        shorthand = Feature(name, **flat)
        assert type(structured) is type(shorthand)
        assert structured.name == shorthand.name
        # ``.options`` of a composition is the shorthand's broadcast projection
        # plus the common count; the structured route has nothing to broadcast,
        # so only the common count and the per-stage options are compared.
        assert structured.options['n_runs'] == shorthand.options['n_runs']
        assert structured.is_stochastic == shorthand.is_stochastic
        assert stage_options(structured) == stage_options(shorthand)
        problems = [FeaturedProblem(constrained_problem(), feature, 12, seed) for feature in (structured, shorthand)]
        np.testing.assert_array_equal(problems[0].x0, problems[1].x0)
        for x in (problems[0].x0, problems[0].x0 + 0.3, np.zeros(3)):
            assert_same_observation(problems[0], problems[1], x)
        assert problems[0].n_eval_fun == problems[1].n_eval_fun
        np.testing.assert_array_equal(problems[0].fun_hist, problems[1].fun_hist)
        np.testing.assert_array_equal(problems[0].maxcv_hist, problems[1].maxcv_hist)


def assert_same_observation(first, second, x):
    np.testing.assert_array_equal(first.fun(x), second.fun(x))
    np.testing.assert_array_equal(first.cub(x), second.cub(x))
    np.testing.assert_array_equal(first.ceq(x), second.ceq(x))
    np.testing.assert_array_equal(first.maxcv(x), second.maxcv(x))


class TestErrors:

    @pytest.mark.parametrize('spec, error, message', [
        ([], ValueError, 'at least one stage entry'),
        ((), ValueError, 'at least one stage entry'),
        ({}, ValueError, r'entry 1 of the feature specification: missing "name"'),
        ([{'name': 'noisy', 'level': 1e-2}], ValueError, r"entry 1 of the feature specification: unknown key\(s\) \['level'\]"),
        ([{'name': 'noisy+truncated'}], ValueError, 'entry 1 of the feature specification: stage names are atomic'),
        ([{'name': 'unknown'}], ValueError, "entry 1 of the feature specification: unknown feature 'unknown'"),
        (['noisy', ''], ValueError, "entry 2 of the feature specification: unknown feature ''"),
        ([{'name': 5}], TypeError, r'entry 1 of the feature specification: "name" must be a string'),
        ([{'name': 'noisy', 'options': 5}], TypeError, r'entry 1 of the feature specification: "options" must be a mapping'),
        ([{'name': 'noisy', 'options': None}], TypeError, r'entry 1 of the feature specification: "options" must be a mapping'),
        ([5], TypeError, 'entry 1 of the feature specification: a stage entry must be a mapping'),
        ([{'name': 'noisy', 'options': {1: 2}}], TypeError,
         r"entry 1 of the feature specification \(stage 'noisy'\): option names must be strings"),
        (['noisy', {'name': 'perturbed_x0', 'options': {'noise_level': 1e-2}}], ValueError,
         r"entry 2 of the feature specification \(stage 'perturbed_x0'\): Option `noise_level` is not valid for feature 'perturbed_x0'"),
        ([{'name': 'noisy', 'options': {'n_runs': 3}}], ValueError,
         r"entry 1 of the feature specification \(stage 'noisy'\): Option `n_runs` is experiment-wide"),
        ([{'name': 'plain', 'options': {'noise_level': 1.0}}, 'noisy'], ValueError,
         r"entry 1 of the feature specification \(stage 'plain'\): Option `noise_level` is not valid for feature 'plain'"),
        ([{'name': 'plain', 'options': {'n_runs': 1}}], ValueError,
         r"entry 1 of the feature specification \(stage 'plain'\): Option `n_runs` is experiment-wide"),
        ([{'name': 'noisy', 'options': {'noise_type': 'loud'}}], ValueError,
         r"entry 1 of the feature specification \(stage 'noisy'\): Option `noise_type` must be one of"),
        ([{'name': 'truncated', 'options': {'significant_digits': 2.5}}], TypeError,
         r"entry 1 of the feature specification \(stage 'truncated'\): Option `significant_digits` must be an integer"),
    ])
    def test_invalid_specifications(self, spec, error, message):
        with pytest.raises(error, match=message):
            Feature(spec)

    def test_hostile_entry_and_option_keys_are_rejected_without_running_their_hooks(self):
        class HostileKey:
            hooks = 0

            def __hash__(self):
                return 1

            def __eq__(self, other):
                HostileKey.hooks += 1
                raise RuntimeError('parser_called_eq')

            def __str__(self):
                HostileKey.hooks += 1
                raise RuntimeError('parser_called_str')

            def __repr__(self):
                HostileKey.hooks += 1
                raise RuntimeError('parser_called_repr')

        with pytest.raises(TypeError, match='entry 1 of the feature specification: entry keys must be strings'):
            Feature([{'name': 'noisy', HostileKey(): 1e-2}])
        with pytest.raises(TypeError, match=r"entry 1 of the feature specification \(stage 'noisy'\): option names must be strings"):
            Feature([{'name': 'noisy', 'options': {HostileKey(): 1e-2}}])
        assert HostileKey.hooks == 0

    def test_flat_stage_options_are_rejected_with_a_specification(self):
        with pytest.raises(ValueError, match=r"only `n_runs` may be given as a keyword.*\['noise_level'\]"):
            Feature(['noisy', 'truncated'], noise_level=1e-2)
        with pytest.raises(ValueError, match=r"only `n_runs` may be given as a keyword"):
            Feature({'name': 'noisy'}, noise_level=1e-2)

    def test_invalid_run_count_with_a_specification(self):
        with pytest.raises(TypeError):
            Feature(['noisy', 'truncated'], n_runs=1.5)
        with pytest.raises(ValueError):
            Feature(['noisy', 'truncated'], n_runs=0)

    @pytest.mark.parametrize('build', [
        lambda **common: Feature('truncated', **common),
        lambda **common: Feature({'name': 'truncated'}, **common),
        lambda **common: Feature('truncated+quantized', **common),
        lambda **common: Feature(['truncated', 'quantized'], **common),
    ])
    def test_explicit_none_run_count_is_rejected_on_every_route(self, build):
        # An explicitly supplied ``n_runs=None`` is a value, not an omission: it
        # must fail the same validation everywhere instead of silently taking
        # the default on the composed routes.
        with pytest.raises(TypeError, match='must be an integer'):
            build(n_runs=None)
        with pytest.raises(ValueError, match='must be positive'):
            build(n_runs=0)
        assert build(n_runs=3).options['n_runs'] == 3
        assert build().options['n_runs'] == 1

    def test_none_and_other_types_are_rejected(self):
        with pytest.raises(TypeError):
            Feature(None)
        with pytest.raises(TypeError):
            Feature(5)
        with pytest.raises(TypeError):
            Feature(b'noisy')
