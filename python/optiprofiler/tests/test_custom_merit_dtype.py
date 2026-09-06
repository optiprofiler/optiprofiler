"""Custom merit callbacks retain real values regardless of their first result."""

import numpy as np
import pytest

from optiprofiler.profile_utils import _default_merit, compute_merit_values


def test_first_integer_merit_does_not_truncate_later_fraction():
    def merit(fun, maxcv, maxcv_init):
        return 0 if fun == 0 else 1.5

    result = compute_merit_values(merit, [0.0, 1.0], [0.0, 0.0], 0.0)

    np.testing.assert_array_equal(result, [0.0, 1.5])


@pytest.mark.parametrize('first, expected', [
    (0, [0.0, 1.5]),
    (np.int64(0), [0.0, 1.5]),
    (False, [0.0, 1.5]),
    (np.bool_(False), [0.0, 1.5]),
    (True, [1.0, 1.5]),
    (np.array(0), [0.0, 1.5]),
])
def test_scalar_result_types_do_not_make_merit_depend_on_input_order(first, expected):
    def merit(fun, maxcv, maxcv_init):
        return first if fun == 0 else 1.5

    forward = compute_merit_values(merit, [0.0, 1.0], [0.0, 0.0], 0.0)
    reverse = compute_merit_values(merit, [1.0, 0.0], [0.0, 0.0], 0.0)

    np.testing.assert_array_equal(forward, expected)
    np.testing.assert_array_equal(reverse, expected[::-1])
    assert np.issubdtype(forward.dtype, np.floating)
    assert forward.shape == reverse.shape == (2,)


@pytest.mark.parametrize('value', [0, np.int64(0), False, np.float64(1.5)])
def test_scalar_inputs_retain_scalar_shape_and_real_dtype(value):
    result = compute_merit_values(lambda f, cv, cv_init: value, 1.0, 0.0, 0.0)

    assert result.shape == ()
    assert np.issubdtype(result.dtype, np.floating)
    assert result == value


def test_fractional_merit_preserves_history_axes_and_initial_value_broadcasting():
    fun = np.arange(16).reshape(2, 2, 2, 2)
    inits = np.array([[0.0, 0.25], [0.5, 0.75]])

    def merit(f, cv, cv_init):
        return 0 if f == 0 else f / 4 + cv + cv_init

    result = compute_merit_values(merit, fun, np.zeros_like(fun), inits)

    expected = np.array([
        [[[0.0, 0.25], [0.75, 1.0]], [[1.0, 1.25], [1.75, 2.0]]],
        [[[2.5, 2.75], [3.25, 3.5]], [[3.5, 3.75], [4.25, 4.5]]],
    ])
    np.testing.assert_array_equal(result, expected)
    assert result.shape == fun.shape


def test_integer_first_result_does_not_discard_nonfinite_callback_results():
    values = [0, np.nan, np.inf, -np.inf]
    result = compute_merit_values(lambda f, cv, cv_init: values[int(f)],
                                  [0.0, 1.0, 2.0, 3.0], np.zeros(4), 0.0)

    np.testing.assert_array_equal(result, [0.0, np.nan, np.inf, -np.inf])


def test_invalid_raw_evaluations_still_override_finite_custom_merit():
    def merit(f, cv, cv_init):
        return 0 if f == 0 else 1.5

    result = compute_merit_values(merit, [0.0, np.nan, 2.0, 3.0],
                                  [0.0, 0.0, np.nan, 0.0], 0.0)

    np.testing.assert_array_equal(result, [0.0, np.inf, np.inf, 1.5])


def test_valid_infinite_inputs_keep_custom_merit_policy():
    def merit(f, cv, cv_init):
        return 0 if f == 0 else 1.5

    result = compute_merit_values(merit, [0.0, np.inf, -np.inf, 1.0],
                                  [0.0, 0.0, 0.0, np.inf], 0.0)

    np.testing.assert_array_equal(result, [0.0, 1.5, 1.5, 1.5])


def test_default_merit_accepts_integer_objectives_followed_by_infeasibility():
    result = compute_merit_values(_default_merit, [0, 1], [0, 1], 0)

    np.testing.assert_array_equal(result, [0.0, np.inf])
