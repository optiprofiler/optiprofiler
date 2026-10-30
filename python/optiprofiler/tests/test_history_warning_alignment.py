"""A history diagnostic must not break the pairing of scoring observations."""

import warnings

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem, benchmark


def malformed_constraint_at_best_point(x):
    return [-1., -1.] if np.all(x == 1.) else [-1.]


@pytest.mark.parametrize('name', ['plain', 'truncated+truncated'])
def test_escalated_history_warning_preserves_paired_histories(name):
    problem = Problem(lambda x: float(np.sum((x - 1.)**2)), [0., 0.],
                      cub=malformed_constraint_at_best_point)
    featured = FeaturedProblem(problem, Feature(name), 2, 17)

    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        with pytest.raises(RuntimeWarning, match='Reference constraint violation'):
            featured.fun([1., 1.])

    # Raising the diagnostic is the caller's policy. The evaluated point
    # nevertheless has an unavailable constraint value, not a missing slot.
    assert featured.n_eval_fun == 1
    np.testing.assert_array_equal(featured.fun_hist, [0.])
    np.testing.assert_array_equal(featured.maxcv_hist, [np.nan])
    assert featured.n_eval_cub == featured.n_eval_ceq == 0

    assert featured.fun([2., 2.]) == 2.
    assert featured.n_eval_fun == 2
    np.testing.assert_array_equal(featured.fun_hist, [0., 2.])
    np.testing.assert_array_equal(featured.maxcv_hist, [np.nan, 0.])

    # Both the cached tail and hard stop still count the escalated query.
    assert featured.fun([3., 3.]) == 2.
    assert featured.fun([4., 4.]) == 2.
    with pytest.raises(StopIteration):
        featured.fun([5., 5.])
    np.testing.assert_array_equal(featured.maxcv_hist, [np.nan, 0.])


@pytest.mark.parametrize('name', ['plain', 'truncated+truncated'])
def test_caught_history_warning_cannot_make_an_unavailable_point_feasible(name):
    problem = Problem(lambda x: float(np.sum((x - 1.)**2)), [0., 0.],
                      cub=malformed_constraint_at_best_point)
    caught = []

    def solver(fun, x0, *args):
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            try:
                fun(np.ones(2))
            except RuntimeWarning as warning:
                caught.append(str(warning))
        fun(np.full(2, 2.))
        return x0

    scores, _, _ = benchmark(
        [solver, solver], problem=problem, feature=Feature(name), n_jobs=1,
        n_runs=1, score_only=True, silent=True, run_plain=False, project_x0=False,
        draw_hist_plots='none', max_eval_factor=1, max_tol_order=1,
    )
    # f=0 has undefined feasibility; the only assessable evaluation is f=2,
    # the initial value. The invalid best objective cannot earn improvement.
    np.testing.assert_array_equal(scores, [0., 0.])
    assert len(caught) == 2
    assert all('Reference constraint violation' in warning for warning in caught)
