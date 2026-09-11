"""
Feature 2.0 ownership: a ``Feature`` is a pipeline specification that owns
stage-local settings only. Repetitions, budgets, seeds and histories belong to
the experiment and to the per-run recorder, never to the specification.

Every expectation is public behaviour (construction, inspection, one benchmark
call) or an explicitly labelled structural check of the raw specification state.
"""

import json
import pickle
import warnings

import matplotlib
import numpy as np
import pytest

matplotlib.use('Agg')

from optiprofiler import benchmark
from optiprofiler.loader import load_results_from_h5
from optiprofiler.opclasses import Feature, FeaturedProblem, Problem

NOISY_DEFAULTS = {'noise_mode': 'random', 'distribution': 'gaussian', 'noise_map': 'chebyshev',
                  'noise_level': 0.001, 'noise_type': 'mixed'}
TRUNCATED_DEFAULTS = {'perturbed_trailing_digits': False, 'significant_digits': 6}
ARRAYS = ('fun_histories', 'maxcv_histories', 'fun_outs', 'maxcv_outs', 'fun_inits', 'maxcv_inits', 'n_evals')


def sphere(x):
    return float(np.dot(x, x))


def ones_distribution(rng, size):
    return np.ones(size)


def solver_stay(fun, x0):
    fun(x0)
    return x0


def solver_step(fun, x0):
    best_x, best_f = np.array(x0, dtype=float), fun(x0)
    for i in range(len(x0)):
        for step in (0.5, -0.5):
            trial = best_x.copy()
            trial[i] += step
            value = fun(trial)
            if value < best_f:
                best_x, best_f = trial, value
    return best_x


def common_kwargs(savepath, benchmark_id):
    return dict(plibs=['s2mpj'], ptype='u', mindim=2, maxdim=2, max_eval_factor=10, benchmark_id=benchmark_id,
                savepath=str(savepath), n_jobs=1, silent=True, draw_hist_plots='none', problem_names=['ROSENBR'])


def archive_of(directory):
    archives = list(directory.rglob('data_for_loading.h5'))
    assert len(archives) == 1
    return str(archives[0])


def stage_views(feature):
    return [(stage.name, stage.occurrence, stage.identity, stage.code, dict(stage.options)) for stage in feature.stages]


class TestCanonicalRepresentation:

    def test_single_feature_is_a_one_stage_pipeline(self):
        feature = Feature('noisy', noise_level=0.5)
        assert stage_views(feature) == [('noisy', 0, 'noisy#0', 2, {**NOISY_DEFAULTS, 'noise_level': 0.5})]
        assert feature.name == 'noisy' and feature.declared_name == 'noisy'
        assert feature.declared.route == 'feature_name'
        assert feature.declared.entries == (('noisy', {'noise_level': 0.5}),)
        assert feature.is_stochastic is True and feature.is_identity is False

    @pytest.mark.parametrize('spec', ['plain', 'plain+plain', ['plain'], ['plain', 'plain'], {'name': 'plain'}])
    def test_identity_input_has_no_effective_stage(self, spec):
        feature = Feature(spec)
        assert feature.stages == ()
        assert feature.name == 'plain'
        assert feature.is_identity is True and feature.is_stochastic is False
        assert all(name == 'plain' and options == {} for name, options in feature.declared.entries)

    def test_identity_declaration_is_retained(self):
        assert Feature('plain+plain').declared_name == 'plain+plain'
        assert Feature(['plain', 'plain']).declared.entries == (('plain', {}), ('plain', {}))
        with pytest.raises(ValueError):
            Feature([{'name': 'plain', 'options': {'noise_level': 1.0}}])

    def test_effective_single_and_multi_stage_share_the_representation(self):
        single = Feature('plain+truncated', significant_digits=4)
        assert stage_views(single) == [('truncated', 0, 'truncated#0', 3, {**TRUNCATED_DEFAULTS, 'significant_digits': 4})]
        assert single.declared_name == 'plain+truncated' and single.name == 'truncated'
        assert single.declared.entries == (('plain', {}), ('truncated', {'significant_digits': 4}))
        multi = Feature([{'name': 'noisy', 'options': {'noise_level': 1e-3}}, 'plain',
                         {'name': 'noisy', 'options': {'noise_level': 1e-1, 'noise_type': 'absolute'}}])
        assert [stage.identity for stage in multi.stages] == ['noisy#0', 'noisy#1']
        assert [stage.options['noise_level'] for stage in multi.stages] == [1e-3, 1e-1]
        assert [stage.options['noise_type'] for stage in multi.stages] == ['mixed', 'absolute']
        assert multi.name == 'noisy+noisy' and multi.declared_name == 'noisy+plain+noisy'
        assert multi.declared.route == 'feature'
        assert multi.declared.entries == (('noisy', {'noise_level': 1e-3}), ('plain', {}),
                                          ('noisy', {'noise_level': 1e-1, 'noise_type': 'absolute'}))

    def test_stage_options_are_owned_and_read_only(self):
        options = {'noise_level': 0.5}
        entries = [{'name': 'noisy', 'options': options}]
        feature = Feature(entries)
        options['noise_level'] = 0.9
        entries.append('truncated')
        assert feature.stages[0].options['noise_level'] == 0.5 and len(feature.stages) == 1
        with pytest.raises(TypeError):
            feature.stages[0].options['noise_level'] = 1.0
        with pytest.raises((AttributeError, TypeError)):
            feature.stages = ()
        # STRUCTURAL: the raw specification state carries no experiment or runtime keys.
        assert 'n_runs' not in dict(feature.stages[0].options)
        assert not hasattr(feature, '_options')

    def test_experiment_options_are_rejected_on_the_specification(self):
        with pytest.raises(ValueError, match=r'benchmark\(.*n_runs=3'):
            Feature('noisy', n_runs=3)
        with pytest.raises(ValueError, match=r'benchmark\(.*n_runs='):
            Feature(['noisy', 'truncated'], n_runs=3)
        with pytest.raises(ValueError, match='experiment'):
            Feature([{'name': 'noisy', 'options': {'n_runs': 3}}])

    def test_feature_object_input_is_reused_without_reparsing(self):
        base = Feature([{'name': 'noisy', 'options': {'noise_level': 0.5}}, 'truncated'])
        again = Feature(base)
        assert stage_views(again) == stage_views(base)
        assert again.name == base.name and again.declared_name == base.declared_name
        assert again.declared.route == base.declared.route and again.declared.entries == base.declared.entries
        with pytest.raises(ValueError, match='benchmark'):
            Feature(base, n_runs=2)
        with pytest.raises(ValueError):
            Feature(base, noise_level=1.0)

    def test_pickle_round_trip_keeps_records_and_native_callables(self):
        feature = Feature([{'name': 'noisy', 'options': {'distribution': ones_distribution}}, 'quantized'])
        clone = pickle.loads(pickle.dumps(feature))
        assert stage_views(clone) == stage_views(feature)
        assert clone.stages[0].options['distribution'] is ones_distribution
        assert clone.name == feature.name and clone.declared.route == feature.declared.route

    def test_deprecated_single_stage_conveniences(self):
        with pytest.warns(DeprecationWarning):
            assert Feature('noisy', noise_level=0.5).options == {**NOISY_DEFAULTS, 'noise_level': 0.5}
        with pytest.warns(DeprecationWarning):
            assert Feature('plain').options == {}
        with pytest.raises(ValueError, match='stages'):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                Feature('noisy+truncated').options
        problem = Problem(sphere, np.array([1.0, 2.0]))
        with pytest.warns(DeprecationWarning):
            x0 = Feature('perturbed_x0').modifier_x0(3, problem)
        assert x0.shape == (2,) and not np.array_equal(x0, problem.x0)
        with pytest.raises(ValueError, match='stages'):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                Feature('noisy+truncated').modifier_x0(3, problem)


class TestFreshRuntimePerTrial:

    def test_featured_problems_share_a_specification_without_sharing_state(self):
        feature = Feature([{'name': 'noisy', 'options': {'noise_type': 'absolute', 'distribution': ones_distribution}},
                           'truncated'])
        problem = Problem(sphere, np.array([1.0, 2.0]))
        first = FeaturedProblem(problem, feature, 5, 7)
        second = FeaturedProblem(problem, feature, 5, 7)
        x = np.array([0.3, -0.7])
        assert first.fun(x) == second.fun(x)
        assert first.n_eval_fun == 1 and second.n_eval_fun == 1
        assert first.fun(x) == second.fun(x)
        assert stage_views(feature) == stage_views(Feature(feature))
        # STRUCTURAL: nothing about the trial leaked into the specification.
        assert not any(name.startswith('_fun') or name.startswith('_maxcv') for name in vars(feature))


class TestBenchmarkWithFeatureObjects:

    def test_feature_object_runs_and_is_reusable(self, tmp_path):
        spec = Feature([{'name': 'noisy', 'options': {'noise_level': 1e-2}}, 'truncated'])
        first = benchmark([solver_stay, solver_step], feature=spec, n_runs=2, **common_kwargs(tmp_path, 'a'))
        second = benchmark([solver_stay, solver_step], feature=spec, n_runs=2, **common_kwargs(tmp_path, 'b'))
        third = benchmark([solver_stay, solver_step], feature=[{'name': 'noisy', 'options': {'noise_level': 1e-2}}, 'truncated'],
                          n_runs=2, **common_kwargs(tmp_path, 'c'))
        np.testing.assert_array_equal(first[0], second[0])
        np.testing.assert_array_equal(first[0], third[0])
        results = [load_results_from_h5(archive_of(tmp_path / name))[0] for name in ('a', 'b', 'c')]
        for key in ARRAYS:
            np.testing.assert_array_equal(results[0][key], results[1][key], err_msg=key)
            np.testing.assert_array_equal(results[0][key], results[2][key], err_msg=key)
        assert results[0]['fun_histories'].shape[2] == 2
        assert len({r['feature_stamp'] for r in results}) == 1
        payload = json.loads(results[0]['feature_pipeline'])
        assert payload['schema'] == 'feature_pipeline-v3'
        assert payload['feature']['route'] == 'feature'
        assert [stage['identity'] for stage in payload['feature']['stages']] == ['noisy#0', 'truncated#0']
        assert all('n_runs' not in stage['options'] for stage in payload['feature']['stages'])
        assert payload['feature']['stages'][0]['options']['noise_level'] == 0.01
        assert payload['experiment']['role'] == 'primary'
        assert payload['experiment']['n_runs'] == 2 and payload['experiment']['origin'] == 'explicit'
        assert stage_views(spec) == stage_views(Feature([{'name': 'noisy', 'options': {'noise_level': 1e-2}}, 'truncated']))

    def test_feature_object_obeys_the_structured_route_rules(self, tmp_path):
        spec = Feature('noisy')
        with pytest.raises(ValueError, match='cannot both be given'):
            benchmark([solver_stay, solver_step], feature=spec, feature_name='noisy', **common_kwargs(tmp_path, 'rejected'))
        with pytest.raises(ValueError, match='only'):
            benchmark([solver_stay, solver_step], feature=spec, noise_level=1e-2, **common_kwargs(tmp_path, 'rejected'))
        with pytest.raises(TypeError, match='n_runs'):
            benchmark([solver_stay, solver_step], feature=spec, n_runs=None, **common_kwargs(tmp_path, 'rejected'))
        assert not (tmp_path / 'rejected').exists()

    def test_default_run_count_is_resolved_by_the_experiment(self, tmp_path):
        benchmark([solver_stay, solver_step], feature=Feature('truncated+quantized'), **common_kwargs(tmp_path, 'hint'))
        results = load_results_from_h5(archive_of(tmp_path / 'hint'))[0]
        assert results['fun_histories'].shape[2] == 1
        payload = json.loads(results['feature_pipeline'])
        assert payload['experiment']['n_runs'] == 1 and payload['experiment']['origin'] == 'stage_hints'
        benchmark([solver_stay, solver_step], feature=Feature('truncated+quantized'), solver_isrand=[True, False],
                  **common_kwargs(tmp_path, 'isrand'))
        payload = json.loads(load_results_from_h5(archive_of(tmp_path / 'isrand'))[0]['feature_pipeline'])
        assert payload['experiment']['n_runs'] == 5 and payload['experiment']['origin'] == 'randomized_solvers'
