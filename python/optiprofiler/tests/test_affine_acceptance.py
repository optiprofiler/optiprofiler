"""Public-constructor regressions from the independent affine acceptance review."""

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem


def build(problem, A, b, inverse, composed):
    stage = {'name': 'custom', 'options': {'mod_affine': lambda rng, predecessor: (A, b, inverse)}}
    feature = Feature([stage, 'noisy'] if composed else [stage])
    return FeaturedProblem(problem, feature, 10, 7)


@pytest.mark.parametrize('composed', [False, True])
@pytest.mark.parametrize('scale', [1e14, 1e15])
def test_singular_map_cannot_retain_an_optimum(composed, scale):
    # The proposed inverse annihilates A on both sides. A termwise rounding
    # allowance greater than one must not certify this as an identity.
    A = np.ones((2, 2))
    inverse = scale * np.array([[1., -1.], [-1., 1.]])
    problem = Problem(lambda x: (x[0] - 1.) ** 2 + (x[1] + 1.) ** 2, [0., 0.],
                      reference={'merit': 0., 'kind': 'optimum', 'source': 'analytic',
                                 'mapping': 'feasible_objective/1'})
    with pytest.raises(ValueError, match='affine transformation'):
        build(problem, A, np.zeros(2), inverse, composed)


@pytest.mark.parametrize('composed', [False, True])
@pytest.mark.parametrize('A,inverse,x0,b', [
    (0.5, 2., 1e308, 0.),
    (1e-300, 1e300, 1e100, 0.),
    (2., 0.5, 1e308, -1e308),
])
def test_finite_start_cannot_escape_through_infinite_comparison(composed, A, inverse, x0, b):
    # The constant objective isolates the representation guard: it does not
    # itself reject Inf or hide the constructor's failure behind an overflow.
    problem = Problem(lambda x: 1., [x0])
    with pytest.raises(ValueError, match='initial point'):
        build(problem, np.array([[A]]), np.array([b]), np.array([[inverse]]), composed)


@pytest.mark.parametrize('composed', [False, True])
@pytest.mark.parametrize('x0,A,inverse,expected', [(1e308, 1., 1., 1e308), (1., .5, 2., 2.)])
def test_representable_start_remains_valid(composed, x0, A, inverse, expected):
    # Summing unscaled absolute terms overflows for the valid huge identity.
    featured = build(Problem(lambda x: 1., [x0]), np.array([[A]]), np.zeros(1),
                     np.array([[inverse]]), composed)
    np.testing.assert_array_equal(featured.x0, [expected])


@pytest.mark.parametrize('composed', [False, True])
def test_large_finite_allowance_does_not_accept_an_inaccurate_inverse(composed):
    # The inaccurate candidate is finite. Only the old allowance overflows,
    # hiding its error; the corrected check must use the exact solve instead.
    featured = build(Problem(lambda x: 1., [1e308]), np.eye(1), np.zeros(1),
                     np.array([[1. + 1e-9]]), composed)
    np.testing.assert_array_equal(featured.x0, [1e308])


@pytest.mark.parametrize('composed', [False, True])
def test_overflowing_forward_check_does_not_certify_a_finite_candidate(composed):
    A = np.array([[1., 1., -1.], [0., 1., 0.], [0., 0., 1.]])
    inverse = np.array([[1., -1., 1.], [0., 1., 0.], [0., 0., 1.]])
    # The absolute sum in the verification overflows, even on BLAS kernels
    # that happen to compute a finite signed sum. Refuse an unverifiable map.
    with pytest.raises(ValueError, match='initial point'):
        build(Problem(lambda x: 1., np.full(3, 1e308)), A, np.zeros(3), inverse, composed)


@pytest.mark.parametrize('composed', [False, True])
@pytest.mark.parametrize('dense', [False, True])
def test_translation_cannot_collapse_a_strict_interval(composed, dense):
    A = np.array([[1., 1.], [0., 1.]]) if dense else np.eye(2)
    inverse = np.array([[1., -1.], [0., 1.]]) if dense else np.eye(2)
    problem = Problem(lambda x: 1., [.5, 0.], xl=[0., -np.inf], xu=[1., np.inf])
    with pytest.raises(ValueError, match='bound.*interval'):
        build(problem, A, np.array([1e16, 0.]), inverse, composed)


@pytest.mark.parametrize('composed', [False, True])
def test_scaling_cannot_collapse_a_strict_interval(composed):
    spacing = np.spacing(1.)
    lower, upper = 1.5 + 2 * spacing, 1.5 + 3 * spacing
    assert lower < upper and .75 * lower == .75 * upper
    problem = Problem(lambda x: 1., [lower], xl=[lower], xu=[upper])
    with pytest.raises(ValueError, match='bound.*interval'):
        build(problem, np.array([[4. / 3.]]), np.zeros(1), np.array([[.75]]), composed)


@pytest.mark.parametrize('composed', [False, True])
@pytest.mark.parametrize('lower,upper,shift', [(0., 0., 1e16), (0., 4., 1e16), (0., 1e-200, 0.)])
def test_fixed_and_representable_intervals_are_not_refused(composed, lower, upper, shift):
    featured = build(Problem(lambda x: 1., [lower], xl=[lower], xu=[upper]),
                     np.eye(1), np.array([shift]), np.eye(1), composed)
    np.testing.assert_array_equal(featured.xl, [lower - shift])
    np.testing.assert_array_equal(featured.xu, [upper - shift])


@pytest.mark.parametrize('composed', [False, True])
def test_replaced_bounds_are_not_checked_as_original_bounds(composed):
    problem = Problem(lambda x: 1., [.5], xl=[0.], xu=[1.])
    stage = {'name': 'custom', 'options': {
        'mod_affine': lambda rng, p: (np.eye(1), np.array([1e16]), np.eye(1)),
        'mod_bounds': lambda rng, p: (np.array([-np.inf]), np.array([np.inf]))}}
    featured = FeaturedProblem(problem, Feature([stage, 'noisy'] if composed else [stage]), 10, 7)
    assert np.isneginf(featured.xl[0]) and np.isposinf(featured.xu[0])
