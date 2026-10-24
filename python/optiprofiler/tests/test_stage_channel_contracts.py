"""Explicit observed/reference expectations for every built-in stage kind."""
import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem
from optiprofiler.feature_definitions import STAGE_NAMES, STAGE_CODES
from optiprofiler.composition import _VIEW_CLASSES

# Constants make coordinate transformations irrelevant, so value policies are
# independently checkable without reproducing the implementation's algorithms.
CASES = {
    'plain': ({}, (1., 2., -3.)),
    'custom': ({}, (1., 2., -3.)),
    'perturbed_x0': ({}, (1., 2., -3.)),
    'permuted': ({}, (1., 2., -3.)),
    'linearly_transformed': ({'rotated': False, 'condition_factor': 0.}, (1., 2., -3.)),
    'truncated': ({'significant_digits': 15}, (1., 2., -3.)),
    'quantized': ({}, (1., 2., -3.)),
    'noisy': ({'noise_level': 0.}, (1., 2., -3.)),
    'random_nan': ({'nan_rate': 1.}, (np.nan, np.nan, np.nan)),
    'unrelaxable_constraints': ({'unrelaxable_nonlinear_constraints': True}, (np.inf, 2., -3.)),
    'nonquantifiable_constraints': ({}, (1., 1., 1.)),
}


def test_new_stage_requires_explicit_channel_classification():
    assert set(STAGE_NAMES) == set(CASES)
    assert set(_VIEW_CLASSES) == set(STAGE_CODES)


@pytest.mark.parametrize('name', CASES)
@pytest.mark.parametrize('composed', [False, True])
def test_observed_and_reference_slots(name, composed):
    options, expected = CASES[name]
    stage = {'name': name, 'options': options}
    if composed and name != 'plain':
        feature = Feature([stage, {'name': 'truncated', 'options': {'significant_digits': 15}}])
    else:
        feature = Feature(name, **options)
    root = Problem(lambda x: 1., [1., 2.], cub=lambda x: np.array([2.]),
                   ceq=lambda x: np.array([-3.]))
    fp = FeaturedProblem(root, feature, 3, 17)
    actual = [fp.fun(fp.x0), fp.cub(fp.x0)[0], fp.ceq(fp.x0)[0]]
    np.testing.assert_allclose(actual, expected, equal_nan=True)
    np.testing.assert_array_equal(fp.fun_hist, [1.])
    np.testing.assert_array_equal(fp.cub_hist, [[2.]])
    np.testing.assert_array_equal(fp.ceq_hist, [[-3.]])
