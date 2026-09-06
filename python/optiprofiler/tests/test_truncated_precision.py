"""Regression tests for the significant-digit oracle shared by all channels."""

import numpy as np
import pytest
import warnings
from decimal import DefaultContext, Inexact, localcontext

from optiprofiler.opclasses import Feature, FeaturedProblem, Problem


def evaluate_truncated(channel, values, significant_digits, perturb=False, seed=17):
    values = np.asarray(values)
    problem = Problem(
        lambda x: values[0], np.zeros(1),
        cub=lambda x: values.copy(), ceq=lambda x: values.copy(),
    )
    feature = Feature(
        'truncated', significant_digits=significant_digits,
        perturbed_trailing_digits=perturb,
    )
    return np.atleast_1d(
        getattr(feature, 'modifier_' + channel)(problem.x0, seed, problem, 1)
    )


@pytest.mark.parametrize('channel', ['fun', 'cub', 'ceq'])
@pytest.mark.parametrize('value,digits,expected', [
    (0.0123456, 2, 0.012),
    (-0.012345, 3, -0.0123),
])
def test_small_values_keep_requested_significant_digits(channel, value, digits, expected):
    np.testing.assert_array_equal(
        evaluate_truncated(channel, [value], digits), [expected]
    )


@pytest.mark.parametrize('channel', ['fun', 'cub', 'ceq'])
@pytest.mark.parametrize('value,digits,expected', [
    (1.25, 2, 1.3), (-1.25, 2, -1.3),
    (0.0125, 2, 0.013), (-0.0125, 2, -0.013),
    (1.005, 3, 1.01), (2.675, 3, 2.68),
    (1250., 2, 1300.), (-1250., 2, -1300.),
    (2.05, 2, 2.1), (np.nextafter(2.05, 0.), 2, 2.),
    (1.25 - 1e-10, 2, 1.2), (1.25 + 1e-10, 2, 1.3),
])
def test_decimal_rounding_matches_matlab(channel, value, digits, expected):
    np.testing.assert_array_equal(
        evaluate_truncated(channel, [value], digits), [expected]
    )


@pytest.mark.parametrize('channel', ['fun', 'cub', 'ceq'])
@pytest.mark.parametrize('dtype', [np.float16, np.float32, np.float64])
@pytest.mark.parametrize('value,expected', [(1.25, 1.3), (-1.25, -1.3)])
def test_numpy_scalar_callbacks_use_same_decimal_tie_policy(channel, dtype, value, expected):
    np.testing.assert_array_equal(
        evaluate_truncated(channel, np.array([value], dtype=dtype), 2), [expected]
    )


@pytest.mark.parametrize('channel', ['fun', 'cub', 'ceq'])
@pytest.mark.parametrize('perturb', [False, True])
@pytest.mark.parametrize('value', [np.nan, np.inf, -np.inf])
def test_nonfinite_values_are_preserved_without_rounding_warnings(channel, perturb, value):
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        actual = evaluate_truncated(channel, [value], 2, perturb)
    np.testing.assert_array_equal(actual, [value])


@pytest.mark.parametrize('channel', ['cub', 'ceq'])
def test_mixed_constraint_values_keep_nonfinite_entries_and_round_finite_entries(channel):
    values = [0., np.nan, np.inf, -np.inf, 0.0123456, -1.25]
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        actual = evaluate_truncated(channel, values, 2)
    np.testing.assert_array_equal(actual, [0., np.nan, np.inf, -np.inf, 0.012, -1.3])


@pytest.mark.parametrize('channel', ['fun', 'cub', 'ceq'])
@pytest.mark.parametrize('value,digits,expected', [
    (0., 1, 0.), (-0., 2, -0.), (0.01, 2, 0.01),
    (10., 2, 10.), (99.99, 3, 100.), (-99.99, 3, -100.),
    (1.2345, 1000, 1.2345), (1e-310, 2, 1e-310),
    (np.nextafter(0., 1.), 2, np.nextafter(0., 1.)),
    (np.finfo(float).tiny, 2, 2.2e-308),
    (np.finfo(float).max, 2, np.inf),
])
def test_zero_scale_extremes_and_integer_boundaries(channel, value, digits, expected):
    np.testing.assert_array_equal(
        evaluate_truncated(channel, [value], digits), [expected]
    )


@pytest.mark.parametrize('channel', ['fun', 'cub', 'ceq'])
def test_caller_decimal_context_cannot_change_oracle(channel):
    with localcontext() as context:
        context.prec = 1
        context.Emax = 1
        np.testing.assert_array_equal(evaluate_truncated(channel, [1250.], 2), [1300.])


@pytest.mark.parametrize('channel', ['fun', 'cub', 'ceq'])
@pytest.mark.parametrize('setting', ['bounds', 'traps'])
def test_modified_default_decimal_context_cannot_change_oracle(channel, setting, monkeypatch):
    if setting == 'bounds':
        monkeypatch.setattr(DefaultContext, 'Emax', 1)
        monkeypatch.setattr(DefaultContext, 'Emin', -1)
    else:
        monkeypatch.setitem(DefaultContext.traps, Inexact, True)
    np.testing.assert_array_equal(evaluate_truncated(channel, [1250.], 2), [1300.])


@pytest.mark.parametrize('channel', ['fun', 'cub', 'ceq'])
@pytest.mark.parametrize('sign', [-1., 1.])
@pytest.mark.parametrize('seed', [0, 1, 17])
def test_trailing_perturbation_uses_last_significant_decimal_place(channel, sign, seed):
    coarse = evaluate_truncated(channel, [sign * .0123456], 2, True, seed)[0]
    fine = evaluate_truncated(channel, [sign * .0123456], 3, True, seed)[0]
    coarse_draw = (sign * coarse - .012) / .001
    fine_draw = (sign * fine - .0123) / .0001
    assert 0 <= coarse_draw < 1
    assert 0 <= fine_draw < 1
    assert coarse_draw == pytest.approx(fine_draw, abs=2e-14)
    assert coarse == evaluate_truncated(channel, [sign * .0123456], 2, True, seed)[0]


@pytest.mark.parametrize('channel', ['cub', 'ceq'])
def test_mixed_nonfinite_constraints_can_perturb_finite_values(channel):
    values = [0., np.nan, np.inf, -np.inf, .0123456, -.0123456]
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        actual = evaluate_truncated(channel, values, 2, True)
    np.testing.assert_array_equal(actual[1:4], [np.nan, np.inf, -np.inf])
    assert 0 <= actual[0] < .1
    assert .012 <= actual[4] < .013
    assert -.013 < actual[5] <= -.012
    np.testing.assert_array_equal(actual, evaluate_truncated(channel, values, 2, True))


def test_rounding_preserves_featured_problem_budget_and_true_history():
    problem = Problem(lambda x: .0123456, [0.], cub=lambda x: [.0123456], ceq=lambda x: [-.0123456])
    featured = FeaturedProblem(problem, Feature('truncated', significant_digits=2), 3, seed=17)
    assert (featured.n_eval_fun, featured.n_eval_cub, featured.n_eval_ceq) == (0, 0, 0)
    assert featured.fun(featured.x0) == .012
    np.testing.assert_array_equal(featured.cub(featured.x0), [.012])
    np.testing.assert_array_equal(featured.ceq(featured.x0), [-.012])
    assert (featured.n_eval_fun, featured.n_eval_cub, featured.n_eval_ceq) == (1, 1, 1)
    np.testing.assert_array_equal(featured.fun_hist, [.0123456])
    np.testing.assert_array_equal(featured.cub_hist, [[.0123456]])
    np.testing.assert_array_equal(featured.ceq_hist, [[-.0123456]])


@pytest.mark.parametrize('channel', ['fun', 'cub', 'ceq'])
def test_existing_finite_rng_stream_is_unchanged(channel):
    # Recorded through the public modifiers at pre-fix commit 259c131,
    # using values whose rounding is unchanged by this correction.
    expected = [1.234114681786892] if channel == 'fun' else [
        1.2347436532835494, -1.2397549500007972,
        0.005784930951980017, 123.00900001609428,
    ]
    np.testing.assert_array_equal(
        evaluate_truncated(channel, [1.23456, -1.23456, 0., 123.456], 3, True), expected
    )
