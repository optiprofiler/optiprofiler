"""
Geometry and randomness of composed features: initialization, affine and
permutation transport, stage streams and served-query counters.
"""

import numpy as np
import pytest

from optiprofiler.composition import STAGE_CODES, Stage, stage_seed
from optiprofiler.opclasses import Feature, FeaturedProblem, Problem


def quadratic(z):
    z = np.asarray(z, dtype=float)
    return float(np.sum(np.array([1.0, 2.0, 0.5]) * (z - np.array([0.5, -1.0, 2.0])) ** 2) + 0.5)


def quadratic_grad(z):
    z = np.asarray(z, dtype=float)
    return 2.0 * np.array([1.0, 2.0, 0.5]) * (z - np.array([0.5, -1.0, 2.0]))


def cub_two(z):
    return np.array([z[0] ** 2 + z[1] - 1.5, z[2] - 3.0])


def ceq_one(z):
    return np.array([z[0] * z[1] - 0.2])


X0 = np.array([0.3, -1.2, 2.5])
XL = np.array([-2.0, -3.0, -1.0])
XU = np.array([2.0, 3.0, 4.0])
AUB = np.array([[1.0, 1.0, 0.0]])
BUB = np.array([2.5])
AEQ = np.array([[0.0, 1.0, -1.0]])
BEQ = np.array([-0.5])


def constrained_problem():
    return Problem(quadratic, X0, xl=XL, xu=XU, aub=AUB, bub=BUB, aeq=AEQ, beq=BEQ, cub=cub_two, ceq=ceq_one,
                   name='geometry')


def bounded_problem():
    return Problem(quadratic, X0, xl=XL, xu=XU, grad=quadratic_grad)


class TestPerturbedX0:

    def test_changes_initialization_only(self):
        featured = FeaturedProblem(constrained_problem(), Feature('perturbed_x0+truncated', perturbation_level=0.5), 10, 3)
        root = constrained_problem()
        assert not np.allclose(featured.x0, X0)
        # Spherical perturbation: the step has the exact configured length.
        assert np.linalg.norm(featured.x0 - X0) == pytest.approx(0.5 * max(1.0, np.linalg.norm(X0)))
        # Structure and reference are untouched; the reference at x0 is the root value there.
        np.testing.assert_array_equal(featured.xl, XL)
        np.testing.assert_array_equal(featured.aub, AUB)
        assert featured.fun_init == root.fun(featured.x0)
        assert featured.maxcv_init == root.maxcv(featured.x0)
        x = featured.x0
        featured.fun(x)
        assert featured.fun_hist.tolist() == [root.fun(x)]

    def test_commutes_with_a_value_stage(self):
        a = FeaturedProblem(constrained_problem(), Feature('noisy+perturbed_x0'), 10, 11)
        b = FeaturedProblem(constrained_problem(), Feature('perturbed_x0+noisy'), 10, 11)
        np.testing.assert_array_equal(a.x0, b.x0)
        for x in (a.x0, a.x0 + 0.1, a.x0):
            assert a.fun(x) == b.fun(x)
            np.testing.assert_array_equal(a.cub(x), b.cub(x))
        np.testing.assert_array_equal(a.fun_hist, b.fun_hist)


class TestAffineTransport:

    def test_scaling_then_scaling_composes_and_transports_x0(self):
        # rotated=False gives A = diag(2 ** power) with power spread by sqrt(condition_factor * n / 2).
        feature = Feature('linearly_transformed+linearly_transformed', rotated=False, condition_factor=2.0)
        featured = FeaturedProblem(bounded_problem(), feature, 10, 5)
        spread = np.sqrt(2.0 * 3 / 2)
        d = 2.0 ** np.linspace(-spread / 2, spread / 2, 3)
        # Two identical diagonal scalings: the solver point x maps to d * d * x.
        np.testing.assert_allclose(featured.x0, X0 / (d * d), rtol=1e-12)
        x = np.array([0.4, -0.7, 1.1])
        assert featured.fun(x) == pytest.approx(quadratic(d * d * x), rel=1e-12)
        assert featured.fun_hist[0] == pytest.approx(quadratic(d * d * x), rel=1e-12)
        # Bounds are scaled twice as well and the initial point stays feasible.
        np.testing.assert_allclose(featured.xl, XL / (d * d), rtol=1e-12)
        np.testing.assert_allclose(featured.xu, XU / (d * d), rtol=1e-12)
        assert featured.maxcv_init == 0.0

    def test_observed_value_at_transported_x0_equals_root_value_at_x0(self):
        for name in ('permuted+linearly_transformed', 'linearly_transformed+permuted+truncated',
                     'permuted+permuted+noisy'):
            options = {}
            if 'noisy' in name:
                options['noise_level'] = 0.0
            if 'truncated' in name:
                options['significant_digits'] = 12
            if 'linearly_transformed' in name:
                options['condition_factor'] = 1.0
            feature = Feature(name, **options)
            featured = FeaturedProblem(constrained_problem(), feature, 10, 9)
            assert featured.fun_init == pytest.approx(quadratic(X0), rel=1e-9)
            assert featured.fun(featured.x0) == pytest.approx(quadratic(X0), rel=1e-9)
            assert featured.maxcv_init == pytest.approx(constrained_problem().maxcv(X0), abs=1e-12)

    def test_rotation_turns_bounds_into_linear_constraints(self):
        feature = Feature('linearly_transformed+noisy', condition_factor=0.0)
        featured = FeaturedProblem(bounded_problem(), feature, 10, 2)
        assert np.all(np.isinf(featured.xl)) and np.all(np.isinf(featured.xu))
        # Six finite bounds become six linear inequality rows; the rotation is orthogonal.
        assert featured.aub.shape == (6, 3)
        np.testing.assert_allclose(featured.bub, np.concatenate([XU, -XL]))
        q = featured.aub[:3]
        np.testing.assert_allclose(q @ q.T, np.eye(3), atol=1e-12)
        np.testing.assert_allclose(featured.aub[3:], -q, atol=1e-12)
        # The reference violation is measured on the original problem at the mapped point.
        x = featured.x0 + np.array([3.0, 0.0, 0.0])
        assert featured.maxcv(x) == pytest.approx(bounded_problem().maxcv(q @ x), abs=1e-12)

    def test_permutation_keeps_nonlinear_shapes_and_reference(self):
        featured = FeaturedProblem(constrained_problem(), Feature('permuted+nonquantifiable_constraints'), 10, 4)
        assert featured.m_nonlinear_ub == 2 and featured.m_nonlinear_eq == 1 and featured.ptype == 'n'
        x = featured.x0
        observed = featured.cub(x)
        assert observed.shape == (2,) and set(observed.tolist()) <= {0.0, 1.0}
        np.testing.assert_allclose(featured.cub_hist[0], cub_two(X0))


class TestStreams:

    def test_stage_seeds_are_frozen_literals(self):
        noisy_0 = Stage(0, 'noisy', 0, Feature('noisy'))
        noisy_1 = Stage(1, 'noisy', 1, Feature('noisy'))
        truncated_0 = Stage(2, 'truncated', 0, Feature('truncated'))
        assert STAGE_CODES == {'perturbed_x0': 1, 'noisy': 2, 'truncated': 3, 'permuted': 4,
                               'linearly_transformed': 5, 'random_nan': 6, 'unrelaxable_constraints': 7,
                               'nonquantifiable_constraints': 8, 'quantized': 9, 'custom': 10}
        seeds = {stage_seed(0, noisy_0), stage_seed(0, noisy_1), stage_seed(0, truncated_0), stage_seed(1, noisy_0)}
        assert len(seeds) == 4
        assert all(0 <= seed < 2 ** 32 for seed in seeds)
        assert stage_seed(0, noisy_0) == stage_seed(None, noisy_0)
        # Same identity, same seed, independent of pipeline position.
        assert stage_seed(7, Stage(5, 'noisy', 1, Feature('noisy'))) == stage_seed(7, noisy_1)

    def test_repeated_stages_draw_from_distinct_streams(self):
        options = dict(noise_type='absolute', noise_level=1.0)
        double = FeaturedProblem(bounded_problem(), Feature('noisy+noisy', **options), 10, 0)
        single = FeaturedProblem(bounded_problem(), Feature('noisy+plain', **options), 10, 0)
        base = quadratic(X0)
        two_draws = double.fun(X0) - base
        one_draw = single.fun(X0) - base
        assert two_draws != one_draw and two_draws != 2.0 * one_draw

    def test_reference_reads_do_not_advance_observed_streams(self):
        feature = Feature('noisy+truncated', significant_digits=12)
        quiet = FeaturedProblem(bounded_problem(), feature, 10, 21)
        busy = FeaturedProblem(bounded_problem(), feature, 10, 21)
        outputs_quiet, outputs_busy = [], []
        for x in (X0, X0, X0 + 0.2):
            outputs_quiet.append(quiet.fun(x))
            busy.maxcv(x)
            busy.maxcv(x + 1.0)
            outputs_busy.append(busy.fun(x))
            busy.maxcv(x)
        assert outputs_quiet == outputs_busy

    def test_unrecorded_constraint_queries_still_advance_streams(self):
        featured = FeaturedProblem(constrained_problem(), Feature('noisy+truncated', significant_digits=12), 10, 8)
        first = featured.cub(X0, record_hist=False)
        second = featured.cub(X0, record_hist=False)
        assert not np.array_equal(first, second)
        assert featured.n_eval_cub == 0


class TestCapabilities:

    def test_composite_derivatives_are_explicitly_unavailable(self):
        featured = FeaturedProblem(bounded_problem(), Feature('noisy+truncated'), 10, 0)
        with pytest.raises(NotImplementedError, match='composed feature'):
            featured.grad(X0)

    def test_single_stage_keeps_legacy_gradient_passthrough(self):
        featured = FeaturedProblem(bounded_problem(), Feature('noisy+plain'), 10, 0)
        np.testing.assert_array_equal(featured.grad(X0), quadratic_grad(X0))
