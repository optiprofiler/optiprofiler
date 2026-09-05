"""Invalid raw evaluations must not become finite custom-merit references."""

import numpy as np
import pytest

from optiprofiler.profile_utils import compute_merit_values


@pytest.mark.parametrize('shape', [(3,), (1, 3), (1, 1, 3), (1, 1, 1, 3)])
def test_custom_merit_cannot_hide_nan_objective_or_constraint(shape):
    fun = np.array([1.0, np.nan, 3.0]).reshape(shape)
    maxcv = np.array([0.0, 0.0, np.nan]).reshape(shape)
    merits = compute_merit_values(lambda f, cv, cv_init: 0, fun, maxcv, 0.0)
    np.testing.assert_array_equal(merits, np.array([0.0, np.inf, np.inf]).reshape(shape))


def test_valid_custom_merit_outputs_are_unchanged():
    merits = compute_merit_values(lambda f, cv, cv_init: f + cv,
                                  [1.5, 2.25], [0.0, 0.5], 0.0)
    np.testing.assert_array_equal(merits, [1.5, 2.75])


def test_infinite_values_are_not_confused_with_nan():
    merits = compute_merit_values(lambda f, cv, cv_init: 7.0,
                                  [np.inf, -np.inf], [0.0, np.inf], 0.0)
    np.testing.assert_array_equal(merits, [7.0, 7.0])
