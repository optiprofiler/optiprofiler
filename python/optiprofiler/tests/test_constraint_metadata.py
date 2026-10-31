"""Static constraint data: what construction refuses and what the legal sentinels mean."""

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem


def sphere(x):
    return float(np.dot(x, x))


MALFORMED = [
    ('xl', dict(xl=[np.nan])),
    ('xl', dict(xl=[np.inf])),
    ('xl', dict(xl=[0.0, np.nan])),
    ('xu', dict(xu=[np.nan])),
    ('xu', dict(xu=[-np.inf])),
    ('xu', dict(xu=[1.0, -np.inf])),
    ('bub', dict(aub=[[1.0]], bub=[-np.inf])),
    ('bub', dict(aub=[[1.0]], bub=[np.nan])),
    ('bub', dict(aub=[[1.0], [1.0]], bub=[np.inf, -np.inf])),
    ('beq', dict(aeq=[[1.0]], beq=[np.nan])),
    ('beq', dict(aeq=[[1.0]], beq=[np.inf])),
    ('beq', dict(aeq=[[1.0]], beq=[-np.inf])),
    ('aub', dict(aub=[[np.nan]], bub=[1.0])),
    ('aub', dict(aub=[[-np.inf]], bub=[1.0])),
    ('aub', dict(aub=[[np.nan]], bub=[np.inf])),
    ('aeq', dict(aeq=[[np.inf]], beq=[1.0])),
    ('aeq', dict(aeq=[[np.nan]], beq=[1.0])),
]


@pytest.mark.parametrize('name, options', MALFORMED)
def test_malformed_data_is_refused_before_any_callback(name, options):
    calls = []

    def cub(x):
        calls.append(np.copy(x))
        return np.array([-1.0])

    size = len(options.get('xl', options.get('xu', [0.0])))
    with pytest.raises(ValueError, match=f'`{name}` for problem has the invalid entry'):
        Problem(sphere, np.zeros(size), cub=cub, **options)
    assert calls == []


@pytest.mark.parametrize('name, value', [('_xl', [np.inf, np.inf]), ('_xl', [np.nan, np.nan]), ('_xu', [-np.inf, -np.inf])])
@pytest.mark.parametrize('stages', [['plain'], ['noisy'], ['linearly_transformed'], ['permuted'], ['noisy', 'permuted'],
                                    ['perturbed_x0', 'linearly_transformed']])
def test_no_feature_hides_malformed_data_of_an_old_object(name, value, stages):
    # A rotation poses only finite bounds, so the posed problem would look fine.
    problem = Problem(sphere, [0.0, 0.0], xl=[-1.0, -1.0], xu=[1.0, 1.0])
    setattr(problem, name, np.array(value))
    with pytest.raises(ValueError, match=f'`{name[1:]}` of the problem has the invalid entry'):
        FeaturedProblem(problem, Feature(stages), 10, seed=7)


def test_one_sided_and_absent_bounds_are_legal():
    problem = Problem(sphere, [0.0, 0.0], xl=[-np.inf, 0.0], xu=[1.0, np.inf])
    assert problem.ptype == 'b'
    assert problem.mb == 2
    assert problem.maxcv([2.0, -3.0]) == 3.0
    assert Problem(sphere, [0.0], xl=[-np.inf], xu=[np.inf]).ptype == 'u'


def test_row_with_infinite_right_hand_side_imposes_no_constraint():
    only = Problem(sphere, [5.0], aub=[[1.0]], bub=[np.inf])
    assert only.ptype == 'u'
    assert only.m_linear_ub == 0
    assert only.maxcv([5.0]) == 0.0
    mixed = Problem(sphere, [5.0], aub=[[1.0], [2.0]], bub=[np.inf, 1.0])
    assert mixed.ptype == 'l'
    assert mixed.m_linear_ub == 1
    assert mixed.aub.shape == (2, 1)
    np.testing.assert_array_equal(mixed.bub, [np.inf, 1.0])
    assert mixed.maxcv([5.0]) == 9.0


def test_inconsistent_shapes_name_the_required_shape():
    with pytest.raises(ValueError, match=r'`aub` for problem must have shape \(1, 2\)\.$'):
        Problem(sphere, [0.0, 0.0], aub=[[1.0]], bub=[np.inf])
    with pytest.raises(ValueError, match=r'`aub` for problem must have shape \(2, 1\)\.$'):
        Problem(sphere, [0.0], aub=[[1.0]], bub=[1.0, np.inf])
    with pytest.raises(ValueError, match=r'`aeq` for problem must have shape \(1, 2\)\.$'):
        Problem(sphere, [0.0, 0.0], aeq=[[1.0]], beq=[1.0])


def test_reversed_finite_bounds_are_measured_not_refused():
    problem = Problem(sphere, [0.0], xl=[1.0], xu=[0.0])
    assert problem.ptype == 'b'
    assert problem.maxcv([0.0]) == 1.0
    featured = FeaturedProblem(problem, Feature('plain'), 2, seed=7)
    assert featured.maxcv_init == 1.0
    featured.fun([0.0])
    np.testing.assert_array_equal(featured.maxcv_hist, [1.0])


def with_sentinels():
    return Problem(sphere, [0.5, 0.5], xl=[0.0, -np.inf], xu=[np.inf, 1.0],
                   aub=[[1.0, 1.0], [1.0, 0.0]], bub=[np.inf, 1.0],
                   aeq=[[1.0, -1.0]], beq=[0.0])


SENTINEL_FEATURES = [
    ('plain', {}),
    ('noisy', {}),
    ('perturbed_x0', {}),
    ('permuted', {}),
    ('linearly_transformed', {'rotated': True}),
    ('linearly_transformed', {'rotated': False}),
    ('quantized', {'ground_truth': False}),
    ('perturbed_x0+linearly_transformed', {}),
    ('permuted+noisy', {}),
]


@pytest.mark.parametrize('name, options', SENTINEL_FEATURES)
def test_sentinels_keep_initial_and_history_truth_under_features(name, options):
    problem = with_sentinels()
    featured = FeaturedProblem(problem, Feature(name, **options), 10, seed=7)
    assert problem.m_linear_ub == 1
    assert np.isfinite(featured.maxcv_init)
    for point in (featured.x0, featured.x0 + 3.0):
        featured.fun(point)
        assert np.isfinite(featured.maxcv(point))
    history = featured.maxcv_hist
    assert history.shape == (2,) and np.all(np.isfinite(history))
    # The finite rows and bounds stay in force at the far point.
    assert history[1] > 0.0
    assert history[0] == pytest.approx(featured.maxcv_init)


def test_sentinel_row_survives_transport_as_a_row_without_force():
    problem = with_sentinels()
    feature = Feature('linearly_transformed', rotated=True)
    featured = FeaturedProblem(problem, feature, 10, seed=7)
    assert np.isposinf(featured.bub).sum() == 1
    assert featured.aub.shape[0] == featured.bub.size
    # Two finite bounds become rows, as does the finite original row.
    assert featured.m_linear_ub == 3
    assert featured.ptype == 'l'


def test_quantized_ground_truth_measures_only_rows_with_force():
    problem = Problem(sphere, [5.0], aub=[[1.0], [2.0]], bub=[np.inf, 1.0])
    featured = FeaturedProblem(problem, Feature('quantized', ground_truth=True), 10, seed=7)
    assert featured.maxcv_init == 9.0
    only = Problem(sphere, [5.0], aub=[[1.0]], bub=[np.inf])
    featured = FeaturedProblem(only, Feature('quantized', ground_truth=True), 10, seed=7)
    assert featured.ptype == 'u'
    assert featured.maxcv_init == 0.0


BAD_OUTPUTS = [
    ('mod_bounds', lambda rng, p: (np.array([np.nan]), p.xu), 'xl'),
    ('mod_bounds', lambda rng, p: (np.array([np.inf]), p.xu), 'xl'),
    ('mod_bounds', lambda rng, p: (p.xl, np.array([np.nan])), 'xu'),
    ('mod_bounds', lambda rng, p: (p.xl, np.array([-np.inf])), 'xu'),
    ('mod_linear_ub', lambda rng, p: (np.array([[1.0]]), np.array([-np.inf])), 'bub'),
    ('mod_linear_ub', lambda rng, p: (np.array([[1.0]]), np.array([np.nan])), 'bub'),
    ('mod_linear_ub', lambda rng, p: (np.array([[-np.inf]]), np.array([1.0])), 'aub'),
    ('mod_linear_eq', lambda rng, p: (np.array([[1.0]]), np.array([np.inf])), 'beq'),
    ('mod_linear_eq', lambda rng, p: (np.array([[np.nan]]), np.array([0.0])), 'aeq'),
]


@pytest.mark.parametrize('key, callback, name', BAD_OUTPUTS)
def test_single_custom_stage_cannot_pose_malformed_data(key, callback, name):
    problem = Problem(sphere, [0.0])
    with pytest.raises(ValueError, match=f'`{name}` posed by the feature `custom` has the invalid entry'):
        FeaturedProblem(problem, Feature('custom', **{key: callback}), 10, seed=7)


@pytest.mark.parametrize('key, callback, name', BAD_OUTPUTS)
@pytest.mark.parametrize('first', [True, False])
def test_composed_custom_stage_cannot_pose_malformed_data(key, callback, name, first):
    problem = Problem(sphere, [0.0])
    custom = {'name': 'custom', 'options': {key: callback}}
    stages = [custom, 'noisy'] if first else ['noisy', custom]
    position = 1 if first else 2
    with pytest.raises(ValueError, match=f"`{name}` returned by the callback `{key}` of stage {position} 'custom'"):
        FeaturedProblem(problem, Feature(stages), 10, seed=7)


@pytest.mark.parametrize('stages', [
    [{'name': 'custom', 'options': {'mod_linear_ub': lambda rng, p: (np.array([[1.0]]), np.array([np.inf]))}}],
    [{'name': 'custom', 'options': {'mod_linear_ub': lambda rng, p: (np.array([[1.0]]), np.array([np.inf]))}}, 'noisy'],
])
def test_custom_stage_may_pose_the_legal_sentinels(stages):
    def bounds(rng, problem):
        return np.array([-np.inf]), np.array([np.inf])

    stages = [dict(stage, options=dict(stage['options'], mod_bounds=bounds)) if isinstance(stage, dict) else stage
              for stage in stages]
    featured = FeaturedProblem(Problem(sphere, [0.0]), Feature(stages), 10, seed=7)
    assert featured.m_linear_ub == 0
    assert featured.mb == 0
    assert featured.ptype == 'u'


@pytest.mark.parametrize('matrix', [np.diag([2.0, -1.0]), np.array([[1.0, 1.0], [0.0, 1.0]])])
@pytest.mark.parametrize('composed', [False, True])
def test_sentinel_row_keeps_no_force_under_a_custom_affine_map(matrix, composed):
    shift = np.array([0.25, -0.5])
    inverse = np.linalg.inv(matrix)
    stage = {'name': 'custom', 'options': {'mod_affine': lambda rng, p: (matrix, shift, inverse)}}
    problem = with_sentinels()
    featured = FeaturedProblem(problem, Feature([stage, 'noisy'] if composed else [stage]), 10, seed=7)
    assert np.isposinf(featured.bub).sum() == 1
    assert featured.aub.shape[0] == featured.bub.size
    assert featured.maxcv_init == pytest.approx(problem.maxcv(problem.x0))
    point = featured.x0 + 3.0
    featured.fun(point)
    assert featured.maxcv_hist[0] == pytest.approx(problem.maxcv(matrix @ point + shift))
