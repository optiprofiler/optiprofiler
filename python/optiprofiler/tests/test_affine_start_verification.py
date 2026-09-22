"""An inflated inverse must not certify its own inaccurate initial point."""

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem
from optiprofiler.opclasses import _pulled_back


def make_featured(A, inverse, x0, composed, bounds=None, shift=None):
    b = np.zeros(len(x0)) if shift is None else shift
    stage = {'name': 'custom', 'options': {
        'mod_affine': lambda rng, problem: (A, b, inverse)}}
    feature = Feature([stage, 'noisy'] if composed else [stage])
    options = {} if bounds is None else {'xl': bounds[0], 'xu': bounds[1]}
    problem = Problem(lambda x: np.sum((x - x0) ** 2), x0, **options)
    return FeaturedProblem(problem, feature, 10, 7)


@pytest.mark.parametrize('composed', [False, True])
@pytest.mark.parametrize('delta', [3e-14, 1e-14])
@pytest.mark.parametrize('factor', [1., 10.])
@pytest.mark.parametrize('unrelated_coordinate', [False, True])
def test_ill_conditioned_inverse_cannot_move_a_feasible_start(composed, delta, factor, unrelated_coordinate):
    A = np.array([[1., 1.], [1., 1. + delta]])
    # This is the reviewer's explicit fabricated inverse, not factor times
    # the exact inverse. Even factor=1 omits a unit that matters at x0.
    inverse = (factor / delta) * np.array([[1., -1.], [-1., 1.]])
    x0 = np.array([1., 2.])
    lower, upper = [.5, 1.5], [1.5, 2.5]
    if unrelated_coordinate:
        # A global normwise accuracy cap would let this unrelated coordinate
        # conceal the same wrong start in the first two coordinates.
        A = np.block([[A, np.zeros((2, 1))], [np.zeros((1, 2)), np.ones((1, 1))]])
        inverse = np.block([[inverse, np.zeros((2, 1))], [np.zeros((1, 2)), np.ones((1, 1))]])
        x0 = np.append(x0, 1e14)
        lower, upper = lower + [-np.inf], upper + [np.inf]
    try:
        featured = make_featured(A, inverse, x0, composed, bounds=(lower, upper))
    except ValueError as exc:
        assert 'initial point' in str(exc) or 'affine transformation' in str(exc)
    else:
        # A solve can recover this particular structured example accurately.
        # Rejecting an unverifiable map is safe; accepting its old image
        # (0, 1) or (0, 10) with an inflated rounding bound is not.
        np.testing.assert_allclose(A @ featured.x0, x0, rtol=0., atol=1e-8)
        assert featured.maxcv_init == 0.


@pytest.mark.parametrize('composed', [False, True])
def test_well_conditioned_rotation_retains_its_rounding_policy(composed):
    # Like S2MPJ STREG/STREGNE: a perfectly conditioned rotation mixes small
    # and large entries. An unconditional local forward cap would reject
    # its ordinary evaluation rounding, although inverse checking is sound.
    A = .5 * np.array([[1., 1., 1., 1.], [1., -1., 1., -1.],
                      [1., 1., -1., -1.], [1., -1., -1., 1.]])
    x0 = np.array([-1.2, 1., 1e10, 1e10])
    candidate = A.T @ x0
    assert np.any(np.abs(A @ candidate - x0) > np.sqrt(np.finfo(float).eps) * np.maximum(1., np.abs(x0)))
    featured = make_featured(A, A.T, x0, composed)
    # The local threshold requests an independent solve, not a refusal of a
    # representable problem merely because its coordinates have mixed scales.
    np.testing.assert_array_equal(featured.x0, np.linalg.solve(A, x0))


@pytest.mark.parametrize('condition_factor', [0, 10, 6000])
def test_generated_rotation_does_not_gain_a_supplied_inverse_solve_trigger(condition_factor):
    problem = Problem(lambda x: np.sum((x / 1e10) ** 2), [-1.2, 1., 1e10, 1e10])
    feature = Feature('linearly_transformed', rotated=True, condition_factor=condition_factor)
    featured = FeaturedProblem(problem, feature, 10, 7)
    A, b, inverse = featured._runtime.modifier_affine(featured._seed, problem)
    # This also covers the baseline fallback if its finite/rounding guard
    # fails; a custom-candidate threshold must not silently change built-ins.
    np.testing.assert_array_equal(featured.x0, _pulled_back(A, inverse, problem.x0, supplied_inverse=False))


def test_empty_affine_start_remains_supported_by_the_helper():
    np.testing.assert_array_equal(_pulled_back(np.empty((0, 0)), np.empty((0, 0)), np.empty(0)), [])


@pytest.mark.parametrize('composed', [False, True])
@pytest.mark.parametrize('scale', [1e-9, 1e-100])
@pytest.mark.parametrize('delta', [3e-14, 1e-14])
def test_changing_units_cannot_hide_a_bad_inverse(composed, scale, delta):
    A = np.array([[1., 1.], [1., 1. + delta]])
    inverse = (10. / delta) * np.array([[1., -1.], [-1., 1.]])
    x0 = scale * np.array([1., 2.])
    featured = make_featured(A, inverse, x0, composed)
    # A fixed absolute floor of 1 would retain the many-times-wrong candidate
    # after a change of units. Compare the actual independent solve, not an
    # absolute assertion tolerance that would hide the same error in the test.
    np.testing.assert_array_equal(featured.x0, np.linalg.solve(A, x0))
    assert not np.array_equal(featured.x0, inverse @ x0)


@pytest.mark.parametrize('composed', [False, True])
def test_zero_coordinates_cannot_borrow_an_unrelated_coordinate_scale(composed):
    delta = 1e-14
    A = np.array([[1., 1., 0.], [1., 1. + delta, 0.], [0., 0., 1.]])
    inverse = np.eye(3)
    inverse[:2, :2] = (10. / delta) * np.array([[1., -1.], [-1., 1.]])
    x0, shift = np.array([0., 0., 1e14]), np.array([1., 2., 0.])
    featured = make_featured(A, inverse, x0, composed, shift=shift)
    # Even a floor based on min(nonzero(abs(x0))) is unsafe: the independent
    # third coordinate would lend 1e14 to desired zeros in the bad block.
    np.testing.assert_array_equal(featured.x0, np.linalg.solve(A, x0 - shift))
    np.testing.assert_array_equal(A @ featured.x0 + shift, x0)


@pytest.mark.parametrize('composed', [False, True])
def test_centering_cannot_hide_an_error_in_the_pulled_back_right_hand_side(composed):
    delta = 1e-14
    A = np.array([[1., 1.], [1., 1. + delta]])
    inverse = (10. / delta) * np.array([[1., -1.], [-1., 1.]])
    x0 = np.full(2, 1e14)
    shift = x0 + [1., 2.]
    # x0 is large but the right-hand side being pulled back is small. A cap
    # relative only to x0 would accept an error of 8 in the centered block.
    featured = make_featured(A, inverse, x0, composed, bounds=(x0 - .5, x0 + .5), shift=shift)
    np.testing.assert_array_equal(featured.x0, np.linalg.solve(A, x0 - shift))
    np.testing.assert_array_equal(A @ featured.x0 + shift, x0)
    assert featured.maxcv_init == 0.


@pytest.mark.parametrize('composed', [False, True])
def test_large_shift_cannot_relax_the_accuracy_of_a_zero_target(composed):
    A, inverse = np.array([[3.]]), np.array([[1. / 3.]])
    x0, shift = np.zeros(1), np.array([2. ** 46 + 2. ** -6])
    # The target is zero but the right-hand side is large. Measuring accuracy
    # only against x0-b would retain the multiply-by-inverse candidate whose
    # image is 1/64, although direct division recovers zero on this example.
    featured = make_featured(A, inverse, x0, composed, shift=shift)
    np.testing.assert_array_equal(featured.x0, np.linalg.solve(A, x0 - shift))
    assert not np.array_equal(featured.x0, inverse @ (x0 - shift))
    np.testing.assert_array_equal(A @ featured.x0 + shift, x0)
    assert featured.fun_init == 0.


@pytest.mark.parametrize('composed', [False, True])
@pytest.mark.parametrize('inverse_error', [0., 1e-9])
def test_nonfinite_allowance_alone_cannot_certify_a_finite_candidate(composed, inverse_error):
    A = np.array([[1., -1.], [0., 1.]])
    inverse = np.array([[1., 1. + inverse_error], [0., 1.]])
    # Candidate and signed mapped point are finite. Only the absolute product
    # overflows; neither of the other finite checks pins this premise.
    with pytest.raises(ValueError, match='initial point'):
        make_featured(A, inverse, np.array([0., 1e308]), composed)


@pytest.mark.parametrize('composed', [False, True])
@pytest.mark.parametrize('delta', [3e-14, 1e-14])
def test_honest_ill_conditioned_map_keeps_a_recoverable_start(composed, delta):
    A = np.array([[1., 1.], [1., 1. + delta]])
    x0 = np.array([1., 2.])
    featured = make_featured(A, np.linalg.inv(A), x0, composed,
                             bounds=([.5, 1.5], [1.5, 2.5]))
    np.testing.assert_allclose(A @ featured.x0, x0, rtol=0., atol=1e-8)
    assert featured.maxcv_init == 0.


@pytest.mark.parametrize('composed', [False, True])
def test_independent_solve_retains_the_evaluation_rounding_policy(composed):
    A = np.array([[1., 1e15], [0., 1.]])
    inverse = np.array([[1., -1e15], [0., 1.]])
    # This valid map has exactly recoverable points, but at (1/3, 1/7) the
    # double-precision coordinate grid and cancellation are too coarse.
    featured = make_featured(A, inverse, np.array([1., 0.]), composed)
    np.testing.assert_array_equal(A @ featured.x0, [1., 0.])
    x0 = np.array([1. / 3., 1. / 7.])
    featured = make_featured(A, inverse, x0, composed)
    np.testing.assert_array_equal(featured.x0, np.linalg.solve(A, x0))
    error = np.abs(A @ featured.x0 - x0)
    assert np.any(error > np.sqrt(np.finfo(float).eps) * np.maximum(1., np.abs(x0)))
    assert np.all(error <= 64 * 2 * np.finfo(float).eps * (np.abs(A) @ np.abs(featured.x0) + np.abs(x0)))


@pytest.mark.parametrize('composed', [False, True])
@pytest.mark.parametrize('delta,distortion', [(5e-13, .1), (1e-11, .001)])
def test_smaller_inverse_errors_also_require_an_independent_start(composed, delta, distortion):
    A = np.array([[1., 1.], [1., 1. + delta]])
    inverse = np.linalg.inv(A) + (distortion / delta) * np.array([[1., -1.], [-1., 1.]])
    # Checking only the band where this norm-product budget reaches 1 misses
    # accepted inverse errors that already move the start by 0.1 or 0.001.
    assert 64 * 2 * np.finfo(float).eps * np.linalg.norm(A, np.inf) * np.linalg.norm(inverse, np.inf) < 1.
    x0 = np.array([1., 2.])
    featured = make_featured(A, inverse, x0, composed, bounds=([.5, 1.5], [1.5, 2.5]))
    np.testing.assert_array_equal(featured.x0, np.linalg.solve(A, x0))
    np.testing.assert_array_equal(A @ featured.x0, x0)
    assert featured.maxcv_init == 0.
