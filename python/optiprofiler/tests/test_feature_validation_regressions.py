"""Configuration errors must fail before an experiment builds runtime state."""

import pickle

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem, benchmark
from optiprofiler.provenance import describe_feature
from optiprofiler.opclasses import _rebuild_feature
from optiprofiler.legacy_compat import LegacyFeature, import_legacy_feature, load_options, replay_arguments


def identity_objective(x):
    return float(x[0])


def stay_at_initial_point(fun, x0):
    fun(x0)
    return x0


@pytest.mark.parametrize('composed', [False, True])
@pytest.mark.parametrize('mesh_type', ['relative', 'RELATIVE', 'Relative'])
def test_relative_mesh_case_survives_native_transport(composed, mesh_type):
    entries = [{'name': 'quantized', 'options': {
        'mesh_type': mesh_type, 'mesh_size': 0.5, 'ground_truth': False}}]
    if composed:
        entries.append({'name': 'noisy', 'options': {'noise_level': 0.0}})
    feature = Feature(entries)
    for candidate in (feature, pickle.loads(pickle.dumps(feature))):
        assert candidate.stages[0].options['mesh_type'] == 'relative'
        assert describe_feature(candidate)['stages'][0]['options']['mesh_type'] == 'relative'
        problem = FeaturedProblem(Problem(identity_objective, [2.6]), candidate, 10, 17)
        assert problem.fun(np.array([2.6])) == pytest.approx(2.6)


MAGNITUDES = [('perturbed_x0', 'perturbation_level'), ('noisy', 'noise_level'),
              ('linearly_transformed', 'condition_factor'), ('quantized', 'mesh_size'),
              ('random_nan', 'nan_rate')]


@pytest.mark.parametrize('name, option', MAGNITUDES)
@pytest.mark.parametrize('structured', [False, True])
@pytest.mark.parametrize('value', [float('nan'), float('inf'), -float('inf'),
                                   'bad', None, 1j, True])
def test_nonfinite_or_nonscalar_magnitude_rejected(name, option, structured, value):
    with pytest.raises((ValueError, TypeError), match=option):
        if structured:
            Feature([{'name': name, 'options': {option: value}}])
        else:
            Feature(name, **{option: value})


@pytest.mark.parametrize('name, option', MAGNITUDES)
@pytest.mark.parametrize('value', [0.1, np.float32(0.1), np.float64(0.1)])
def test_finite_real_scalar_magnitude_supported(name, option, value):
    feature = Feature(name, **{option: value})
    assert feature.stages[0].options[option] == value


@pytest.mark.parametrize('name, option', MAGNITUDES)
def test_magnitude_sign_and_zero(name, option):
    with pytest.raises(ValueError, match=option):
        Feature(name, **{option: -0.1})
    if option == 'mesh_size':
        with pytest.raises(ValueError, match=option):
            Feature(name, **{option: 0.0})
    else:
        assert Feature(name, **{option: 0.0}).stages[0].options[option] == 0.0


@pytest.mark.parametrize('level', [[0.1, 0.2], (0.1, 0.2), np.array([0.1, 0.2])])
def test_coordinatewise_perturbation_keeps_array_behavior(level):
    base = Problem(identity_objective, [1.0, 2.0])
    expected = FeaturedProblem(base, Feature('perturbed_x0', perturbation_level=np.array([0.1, 0.2])), 5, 17).x0
    actual = FeaturedProblem(base, Feature('perturbed_x0', perturbation_level=level), 5, 17).x0
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize('level', [[], [[0.1]], [0.1, float('nan')], [0.1, -0.2],
                                   [0.1, 'bad'], np.array([object()], dtype=object)])
def test_invalid_perturbation_vectors_fail_before_runtime(level):
    with pytest.raises((TypeError, ValueError), match='perturbation_level'):
        Feature('perturbed_x0', perturbation_level=level)


def test_perturbation_dimension_is_checked_when_problem_is_known():
    spec = Feature('perturbed_x0', perturbation_level=[0.1, 0.2])
    with pytest.raises(ValueError, match='perturbation_level'):
        FeaturedProblem(Problem(identity_objective, [1.0, 2.0, 3.0]), spec, 5, 17)


@pytest.mark.parametrize('value', [True, np.bool_(True), float('nan'), float('inf')])
def test_significant_digits_is_a_count_not_a_boolean(value):
    with pytest.raises((TypeError, ValueError), match='significant_digits'):
        Feature('truncated', significant_digits=value)


@pytest.mark.parametrize('legacy', [False, True])
def test_historical_uppercase_mesh_replays_original_oracle(legacy):
    # e7 and 1.x accepted this spelling but executed an absolute grid. The
    # fresh-input normalization fix must not rewrite saved experiment meaning.
    options = {'mesh_type': 'RELATIVE', 'mesh_size': 0.5, 'ground_truth': False}
    if legacy:
        feature = import_legacy_feature(LegacyFeature({
            '_name': 'quantized', '_options': options})).feature
    else:
        feature = _rebuild_feature(1, 'feature', (('quantized', options),),
                                   (('quantized', options),))
        assert feature.declared.entries[0][1]['mesh_type'] == 'RELATIVE'
    problem = FeaturedProblem(Problem(identity_objective, [2.6]), feature, 10, 17)
    assert problem.fun(np.array([2.6])) == pytest.approx(2.5)
    assert feature.stages[0].options['mesh_type'] == 'absolute'


@pytest.mark.parametrize('layout', ['refined', 'bridge', 'flat'])
def test_saved_options_replay_preserves_mesh_without_mutating_source(layout):
    options = {'mesh_type': 'RELATIVE', 'mesh_size': 0.5, 'ground_truth': False}
    if layout == 'flat':
        saved = {**options, 'feature_name': 'quantized', 'n_runs': 1}
    else:
        saved = {'feature_specification': [{'name': 'quantized', 'options': options}], 'n_runs': 1}
        if layout == 'refined':
            saved['schema'] = 'options_refined-v2'
    original = pickle.dumps(saved)
    replay = replay_arguments(saved)
    feature = Feature(replay['feature'])
    assert pickle.dumps(saved) == original
    assert options['mesh_type'] == 'RELATIVE'
    problem = FeaturedProblem(Problem(identity_objective, [2.6]), feature, 10, 17)
    assert problem.fun(np.array([2.6])) == pytest.approx(2.5)


@pytest.mark.parametrize('route', ['feature_name', 'feature', 'instance'])
def test_current_user_and_refined_archives_replay_the_current_grid(route, tmp_path):
    options = {'mesh_type': 'RELATIVE', 'mesh_size': 0.5, 'ground_truth': False}
    entries = [{'name': 'quantized', 'options': options}]
    request = ({'feature_name': 'quantized', **options} if route == 'feature_name'
               else {'feature': Feature(entries) if route == 'instance' else entries})
    benchmark([stay_at_initial_point, stay_at_initial_point],
              problem=Problem(identity_objective, [2.6]), **request,
              n_runs=1, n_jobs=1, silent=True, draw_hist_plots='none',
              savepath=str(tmp_path), benchmark_id='current-grid')
    for filename in ('options_user.pkl', 'options_refined.pkl'):
        paths = list(tmp_path.rglob(filename))
        assert len(paths) == 1
        saved = load_options(paths[0])
        before = pickle.dumps(saved)
        replay = replay_arguments(saved)
        feature = Feature(replay['feature'])
        fp = FeaturedProblem(Problem(identity_objective, [2.6]), feature, 5, 17)
        assert fp.fun(np.array([2.6])) == pytest.approx(2.6)
        assert pickle.dumps(saved) == before
        if filename == 'options_user.pkl' and route != 'instance':
            raw = saved if route == 'feature_name' else saved['feature'][0]['options']
            assert raw['mesh_type'] == 'RELATIVE'


def test_marked_user_options_keep_stage_local_repeated_options():
    saved = {'schema': 'options_user-v2', 'feature': [
        {'name': 'noisy', 'options': {'noise_level': 0.1, 'distribution': 'uniform'}},
        {'name': 'noisy', 'options': {'noise_level': 0.2, 'distribution': 'gaussian'}}],
        'n_runs': 3, 'mindim': 2, 'maxdim': 8}
    replay = replay_arguments(saved)
    stages = Feature(replay['feature']).stages
    assert [s.options['noise_level'] for s in stages] == [0.1, 0.2]
    assert [s.options['distribution'] for s in stages] == ['uniform', 'gaussian']
    assert replay['n_runs'] == 3
    assert replay['problem_options'] == {'mindim': 2, 'maxdim': 8}


def test_marked_user_options_without_explicit_feature_replay_as_plain():
    assert Feature(replay_arguments({'schema': 'options_user-v2'})['feature']).is_identity


@pytest.mark.parametrize('extra', [{'feature_name': 'noisy'}, {'noise_level': 0.1}])
def test_marked_user_options_reject_mixed_routes(extra):
    with pytest.raises(ValueError, match='cannot be mixed'):
        replay_arguments({'schema': 'options_user-v2', 'feature': 'noisy', **extra})
