"""Quantized truth must agree at initialization, history, and solver output."""

import copy

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem
from optiprofiler.profile_utils import get_default_profile_options
from optiprofiler.profiles import _solve_one_problem


def make_problem(kind='cub', x0=0.49):
    options = {}
    if kind in ('cub', 'mixed'):
        options['cub'] = lambda x: np.array([x[0] ** 2 - 0.25 ** 2])
    if kind in ('ceq', 'mixed'):
        options['ceq'] = lambda x: np.array([x[0] ** 2])
    if kind == 'linear':
        options.update(xl=[0.2], xu=[0.4], aub=[[1.0]], bub=[0.3])
    return Problem(lambda x: x[0], [x0], name='QUANTIZED_TRUTH', **options)


def expected_truth(problem, kind, truth, mesh_type):
    x = problem.x0[0]
    mesh = max(1.0, abs(x)) if mesh_type == 'relative' else 1.0
    q = mesh * np.round(x / mesh) if truth else x
    if kind == 'linear':
        # The feature quantizes nonlinear oracles, not bounds/linear residuals.
        cv = problem.maxcv(problem.x0)
    else:
        cv = max([0.0] + ([q * q - 0.25 ** 2] if kind in ('cub', 'mixed') else [])
                 + ([abs(q * q)] if kind in ('ceq', 'mixed') else []))
    return q, cv


@pytest.mark.parametrize('truth', [False, True])
@pytest.mark.parametrize('kind', ['u', 'cub', 'ceq', 'mixed', 'linear'])
@pytest.mark.parametrize('mesh_type,x0', [('absolute', 0.49), ('absolute', -1.49), ('relative', 1.49)])
def test_initial_history_and_actual_solver_output_agree(truth, kind, mesh_type, x0):
    problem = make_problem(kind, x0)
    feature = Feature('quantized', mesh_size=1.0, mesh_type=mesh_type,
                      ground_truth=truth, n_runs=1)
    expected_f, expected_cv = expected_truth(problem, kind, truth, mesh_type)

    def solver(fun, x, *constraints):
        original = x.copy()
        fun(x)
        if kind in ('cub', 'ceq', 'mixed'):
            constraints[-2](x)
            constraints[-1](x)
        np.testing.assert_array_equal(x, original)
        return x

    options = get_default_profile_options([solver], feature, dict(
        solver_names=['stay'], solver_isrand=[False], project_x0=False,
        max_eval_factor=1, silent=True, solver_verbose=2, score_only=True,
        draw_hist_plots='none', seed=17))
    result = _solve_one_problem([solver], problem, feature, problem.name,
                                len(problem.name), options, False, None)
    for field in ('fun_init', 'fun_history', 'fun_out'):
        np.testing.assert_allclose(result[field], expected_f)
    for field in ('maxcv_init', 'maxcv_history', 'maxcv_out'):
        np.testing.assert_allclose(result[field], expected_cv)
    np.testing.assert_array_equal(result['n_eval'], [[1]])
    assert not result['solver_abnormal_termination'].any()
    assert not result['solver_output_fallback'].any()
    np.testing.assert_array_equal(problem.x0, [x0])


def test_quantized_constraint_truth_does_not_consume_budget_or_change_oracle_cache():
    fp = FeaturedProblem(make_problem('mixed'), Feature('quantized', mesh_size=1.0), 1, seed=17)
    assert [fp._real_n_eval_fun, fp._real_n_eval_cub, fp._real_n_eval_ceq] == [0, 0, 0]
    fp.fun(fp.x0)
    fp.cub(fp.x0)
    fp.ceq(fp.x0)
    state_names = ['_real_n_eval_fun', '_real_n_eval_cub', '_real_n_eval_ceq',
                   '_last_fun', '_last_cub', '_last_ceq', '_fun_hist', '_cub_hist',
                   '_ceq_hist', '_maxcv_hist']
    before = {name: copy.deepcopy(getattr(fp, name)) for name in state_names}
    rng_before = np.random.get_state()
    # This point is not the last solver evaluation. Truth must not reuse the
    # cached last oracle value after the solver has exhausted its budget.
    for _ in range(8):
        assert fp.maxcv(np.array([1.49])) == 1.0
    for name in state_names:
        np.testing.assert_equal(getattr(fp, name), before[name])
    for actual, expected in zip(np.random.get_state(), rng_before):
        np.testing.assert_equal(actual, expected)


def test_nonquantifiable_nan_stays_undefined():
    p = Problem(lambda x: 0.0, [0.0],
                cub=lambda x: np.array([np.nan, -np.inf, 0.0, np.inf]),
                ceq=lambda x: np.array([np.nan, -np.inf, 0.0, np.inf]))
    feature = Feature('nonquantifiable_constraints')
    np.testing.assert_array_equal(feature.modifier_cub(p.x0, 0, p, 0), [np.nan, 0, 0, 1])
    np.testing.assert_array_equal(feature.modifier_ceq(p.x0, 0, p, 0), [np.nan, 1, 0, 1])


@pytest.mark.parametrize('n', [1, 2, 4])
@pytest.mark.parametrize('rotated', [False, True])
def test_documented_condition_factor_matches_matrix(n, rotated):
    p = Problem(lambda x: x @ x, np.ones(n))
    feature = Feature('linearly_transformed', condition_factor=2, rotated=rotated)
    a, _, inv_a = feature.modifier_affine(17, p)
    assert np.linalg.cond(a) == pytest.approx(2 ** np.sqrt(n) if n > 1 else 1)
    np.testing.assert_allclose(a @ inv_a, np.eye(n), atol=1e-14)
