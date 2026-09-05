"""Regression tests for constraint validity and affine equivalence."""

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem


@pytest.mark.parametrize('kind', ['cub', 'ceq', 'mixed', 'linear'])
def test_nan_constraint_is_not_feasible(kind):
    """An undefined required constraint cannot be hidden by a finite one."""
    options = dict(xl=[-1.0], xu=[1.0])
    if kind == 'linear':
        options.update(aub=[[1.0]], bub=[np.nan])
    else:
        options[kind if kind != 'mixed' else 'cub'] = lambda x: np.array([np.nan, -2.0])
        if kind == 'mixed':
            options['ceq'] = lambda x: np.array([3.0])
    problem = Problem(lambda x: 7.0, [0.0], **options)
    assert np.isnan(problem.maxcv(problem.x0))


@pytest.mark.parametrize('feature_name', ['plain', 'quantized'])
def test_nan_constraint_reaches_initial_value_and_history(feature_name):
    problem = Problem(lambda x: 7.0, [0.0], cub=lambda x: np.array([np.nan]))
    featured = FeaturedProblem(problem, Feature(feature_name), 10, seed=0)
    assert np.isnan(featured.maxcv_init)
    assert np.isnan(featured.maxcv(featured.x0))
    assert featured.fun(featured.x0) == 7.0
    np.testing.assert_array_equal(featured.maxcv_hist, [np.nan])
    assert featured.n_eval_fun == 1


def test_failed_later_constraint_does_not_erase_valid_points():
    problem = Problem(lambda x: x[0] ** 2, [0.0],
                      cub=lambda x: np.array([-1.0 if x[0] == 0 else np.nan]))
    featured = FeaturedProblem(problem, Feature('plain'), 10, seed=0)
    featured.fun([0.0])
    featured.fun([1.0])
    featured.fun([0.0])
    assert featured.maxcv_init == 0.0
    np.testing.assert_array_equal(featured.maxcv_hist, [0.0, np.nan, 0.0])


def test_infinite_bounds_and_satisfied_negative_infinite_constraint_remain_valid():
    problem = Problem(lambda x: 0.0, [0.0], xl=[-np.inf], xu=[np.inf],
                      cub=lambda x: np.array([-np.inf]))
    assert problem.maxcv(problem.x0) == 0.0


@pytest.mark.parametrize('kind', ['cub', 'ceq'])
@pytest.mark.parametrize('container', [list, tuple, np.array])
def test_array_like_constraint_count_is_inferred_once(kind, container):
    calls = []

    def constraint(x):
        calls.append(x.copy())
        return container([x[0] - 1.0])

    problem = Problem(lambda x: x[0] ** 2, [2.0], **{kind: constraint})
    assert len(calls) == 1
    assert getattr(problem, 'm_nonlinear_' + ('ub' if kind == 'cub' else 'eq')) == 1
    assert problem.ptype == 'n'
    assert problem.maxcv([2.0]) == 1.0


@pytest.mark.parametrize('kind', ['cub', 'ceq'])
@pytest.mark.parametrize('empty', [None, [], (), np.empty(0)])
def test_empty_constraint_return_still_means_no_constraints(kind, empty):
    problem = Problem(lambda x: 0.0, [0.0], **{kind: lambda x: empty})
    assert problem.ptype == 'u'
    assert problem.maxcv(problem.x0) == 0.0


@pytest.mark.parametrize('kind', ['cub', 'ceq'])
def test_failed_constraint_count_does_not_construct_unconstrained_problem(kind):
    def failing_constraint(x):
        raise RuntimeError('no constraint values at this point')

    with pytest.raises(ValueError, match='number of nonlinear'):
        Problem(lambda x: 0.0, [0.0], **{kind: failing_constraint})


@pytest.mark.parametrize('kind', ['cub', 'ceq'])
def test_invalid_constraint_shape_is_rejected_during_construction(kind):
    with pytest.raises(ValueError, match='number of nonlinear'):
        Problem(lambda x: 0.0, [0.0], **{kind: lambda x: np.ones((2, 2))})


@pytest.mark.parametrize('x0, aeq, beq, expected', [
    ([1.0], [[1.0]], [3.0], [3.0]),
    ([3.0], [[1.0]], [3.0], [3.0]),
    ([1.0, 4.0], [[1.0, 1.0]], [7.0], [2.0, 5.0]),
    ([1.0, 4.0], [[1.0, 1.0], [2.0, 2.0]], [7.0, 14.0], [2.0, 5.0]),
])
def test_equality_projection_applies_displacement(x0, aeq, beq, expected):
    problem = Problem(lambda x: x @ x, x0, aeq=aeq, beq=beq)
    problem.project_x0()
    np.testing.assert_allclose(problem.x0, expected, atol=1e-12)
    assert problem.maxcv(problem.x0) < 1e-12


@pytest.mark.parametrize('matrix', [
    np.diag([-2.0, 1.0, 0.5]),
    np.array([[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [0.0, 0.0, 1.0]]),
    np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]),
])
@pytest.mark.parametrize('fixed_and_linear', [False, True])
def test_custom_affine_preserves_solver_visible_feasible_region(matrix, fixed_and_linear):
    options = dict(xl=[0.0, 0.0, -np.inf], xu=[1.0, 1.0, np.inf])
    if fixed_and_linear:
        options.update(xl=[0.0, 0.5, -np.inf], xu=[1.0, 0.5, np.inf],
                       aub=[[1.0, 0.0, 1.0]], bub=[1.5],
                       aeq=[[0.0, 0.0, 1.0]], beq=[0.0])
    original = Problem(lambda x: x @ x, [0.5, 0.5, 0.0], **options)
    shift = np.array([0.25, -0.5, 1.0])
    inverse = np.linalg.inv(matrix)
    feature = Feature('custom', mod_affine=lambda rng, p: (matrix, shift, inverse))
    featured = FeaturedProblem(original, feature, 20, seed=0)
    solver_view = Problem(lambda y: y @ y, featured.x0,
                          xl=featured.xl, xu=featured.xu,
                          aub=featured.aub, bub=featured.bub,
                          aeq=featured.aeq, beq=featured.beq)
    diagonal = np.count_nonzero(matrix - np.diag(np.diag(matrix))) == 0
    expected_inequalities = original.m_linear_ub + (0 if diagonal else (2 if fixed_and_linear else 4))
    expected_equalities = original.m_linear_eq + (int(fixed_and_linear) if not diagonal else 0)
    assert solver_view.m_linear_ub == expected_inequalities
    assert solver_view.m_linear_eq == expected_equalities
    for x1 in [-1.0, 0.0, 0.5, 1.0, 2.0]:
        for x2 in [0.0, 0.5, 1.0]:
            for x3 in [-1.0, 0.0, 1.0]:
                x = np.array([x1, x2, x3])
                y = inverse @ (x - shift)
                assert (solver_view.maxcv(y) < 1e-12) == (original.maxcv(x) < 1e-12)


@pytest.mark.parametrize('name, options', [
    ('permuted', {}),
    ('linearly_transformed', {'rotated': True}),
    ('linearly_transformed', {'rotated': False}),
    ('plain', {}), ('noisy', {}), ('perturbed_x0', {}),
])
def test_builtin_feature_solver_visible_affine_constraints(name, options):
    original = Problem(lambda x: x @ x, [0.5, 0.5, 0.0],
                       xl=[0.0, 0.5, -np.inf], xu=[1.0, 0.5, np.inf],
                       aub=[[1.0, 0.0, 1.0]], bub=[1.5],
                       aeq=[[0.0, 0.0, 1.0]], beq=[0.0])
    feature = Feature(name, **options)
    featured = FeaturedProblem(original, feature, 20, seed=7)
    matrix, shift, inverse = feature.modifier_affine(7, original)
    solver_view = Problem(lambda y: y @ y, featured.x0,
                          xl=featured.xl, xu=featured.xu,
                          aub=featured.aub, bub=featured.bub,
                          aeq=featured.aeq, beq=featured.beq)
    for x in [[0.5, 0.5, 0.0], [-1.0, 0.5, 0.0], [0.5, 1.0, 0.0], [0.5, 0.5, 1.0]]:
        x = np.asarray(x)
        y = inverse @ (x - shift)
        assert (solver_view.maxcv(y) < 1e-12) == (original.maxcv(x) < 1e-12)
