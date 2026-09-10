"""
Quantized point maps, the observed unrelaxable gate and custom stages inside
compositions, plus callback evidence for long chains.
"""

import time

import numpy as np
import pytest

from optiprofiler.opclasses import Feature, FeaturedProblem, Problem


class Counting:
    """Root callbacks that record every invocation per channel."""

    def __init__(self):
        self.calls = {'fun': [], 'cub': [], 'ceq': []}

    def wrap(self, channel, function):
        def wrapped(x):
            self.calls[channel].append(np.array(x, dtype=float))
            return function(x)
        return wrapped

    def count(self, channel):
        return len(self.calls[channel])

    def reset(self):
        for channel in self.calls:
            self.calls[channel] = []


def quadratic(z):
    z = np.asarray(z, dtype=float)
    return float(np.sum(np.array([1.0, 2.0, 0.5]) * (z - np.array([0.5, -1.0, 2.0])) ** 2) + 0.5)


def cub_two(z):
    return np.array([z[0] ** 2 + z[1] - 1.5, z[2] - 3.0])


def ceq_one(z):
    return np.array([z[0] * z[1] - 0.2])


X0 = np.array([0.3, -1.2, 2.5])
XL = np.array([-2.0, -3.0, -1.0])
XU = np.array([2.0, 3.0, 4.0])


def snap(x, mesh):
    return mesh * np.round(np.asarray(x, dtype=float) / mesh)


def counting_problem(counting, kind='n'):
    fun = counting.wrap('fun', quadratic)
    if kind == 'u':
        return Problem(fun, X0)
    if kind == 'b':
        return Problem(fun, X0, xl=XL, xu=XU)
    return Problem(fun, X0, xl=XL, xu=XU, aub=np.array([[1.0, 1.0, 0.0]]), bub=np.array([2.5]),
                   cub=counting.wrap('cub', cub_two), ceq=counting.wrap('ceq', ceq_one))


class TestQuantizedPointMap:

    @pytest.mark.parametrize('ground_truth', [False, True])
    def test_one_observed_and_one_reference_read_per_query(self, ground_truth):
        counting = Counting()
        # A zero-level noise stage is an exact identity on values.
        featured = FeaturedProblem(counting_problem(counting), Feature('quantized+noisy', mesh_size=0.3,
                                                                       ground_truth=ground_truth, noise_level=0.0), 10, 0)
        counting.reset()
        x = np.array([0.44, -1.31, 2.66])
        observed = featured.fun(x)
        # One observed and one reference objective read; the scoring violation
        # recorded with the history reads each constraint callback once.
        assert counting.count('fun') == 2 and counting.count('cub') == 1 and counting.count('ceq') == 1
        assert observed == quadratic(snap(x, 0.3))
        expected_reference = quadratic(snap(x, 0.3)) if ground_truth else quadratic(x)
        assert featured.fun_hist[0] == expected_reference
        # The scoring violation keeps bounds unsnapped; nonlinear reference follows the mode.
        cub_point = snap(x, 0.3) if ground_truth else x
        expected_cv = max(0.0, float(np.max(cub_two(cub_point))), abs(float(ceq_one(cub_point)[0])))
        assert featured.maxcv_hist[0] == pytest.approx(expected_cv)

    def test_relative_mesh_and_initial_point(self):
        featured = FeaturedProblem(counting_problem(Counting()), Feature('quantized+noisy', mesh_size=0.3,
                                                                         mesh_type='relative', ground_truth=True,
                                                                         noise_level=0.0), 10, 0)
        np.testing.assert_array_equal(featured.x0, X0)
        mesh = 0.3 * np.maximum(1.0, np.abs(X0))
        assert featured.fun_init == quadratic(mesh * np.round(X0 / mesh))
        assert featured.fun(X0) == quadratic(mesh * np.round(X0 / mesh))

    def test_snapping_composes_in_declared_order(self):
        # a = quantized(0.3), b = linearly_transformed(diagonal scaling d), c = quantized(0.3):
        # the root is evaluated at snap(d * snap(x)).
        feature = Feature('quantized+linearly_transformed+quantized', mesh_size=0.3, rotated=False,
                          condition_factor=2.0, ground_truth=False)
        featured = FeaturedProblem(counting_problem(Counting(), 'u'), feature, 10, 0)
        spread = np.sqrt(2.0 * 3 / 2)
        d = 2.0 ** np.linspace(-spread / 2, spread / 2, 3)
        x = np.array([0.44, -1.31, 2.66])
        assert featured.fun(x) == pytest.approx(quadratic(snap(d * snap(x, 0.3), 0.3)), rel=1e-12)
        assert featured.fun_hist[0] == pytest.approx(quadratic(d * x), rel=1e-12)

    def test_constraint_channels_snap_lazily(self):
        counting = Counting()
        featured = FeaturedProblem(counting_problem(counting), Feature('quantized+quantized', mesh_size=0.3,
                                                                       ground_truth=False), 10, 0)
        counting.reset()
        x = np.array([0.44, -1.31, 2.66])
        np.testing.assert_array_equal(featured.cub(x), cub_two(snap(x, 0.3)))
        np.testing.assert_array_equal(featured.ceq(x), ceq_one(snap(x, 0.3)))
        # One observed read at the snapped point and one unsnapped reference read per channel.
        assert counting.count('cub') == 2 and counting.count('ceq') == 2
        np.testing.assert_array_equal(featured.cub_hist[0], cub_two(x))


class TestLongChains:

    @pytest.mark.parametrize('length', [8, 32, 64])
    @pytest.mark.parametrize('ground_truth', [False, True])
    def test_repeated_quantized_reads_do_not_grow_with_length(self, length, ground_truth):
        counting = Counting()
        feature = Feature('+'.join(['quantized'] * length), mesh_size=0.3, ground_truth=ground_truth)
        featured = FeaturedProblem(counting_problem(counting), feature, 1000, 0)
        counting.reset()
        x = np.array([0.44, -1.31, 2.66])
        featured.fun(x)
        # Objective: one observed and one reference read; the recorded scoring
        # violation reads each constraint callback once, whatever the length.
        assert (counting.count('fun'), counting.count('cub'), counting.count('ceq')) == (2, 1, 1)
        featured.cub(x)
        assert (counting.count('fun'), counting.count('cub'), counting.count('ceq')) == (2, 3, 1)
        featured.ceq(x)
        assert (counting.count('fun'), counting.count('cub'), counting.count('ceq')) == (2, 3, 3)

    def test_mixed_long_chain_runtime_is_linear(self):
        counting = Counting()
        feature = Feature('+'.join(['quantized', 'noisy', 'truncated', 'random_nan'] * 16), mesh_size=0.3)
        featured = FeaturedProblem(counting_problem(counting, 'u'), feature, 10000, 0)
        counting.reset()
        x = np.array([0.44, -1.31, 2.66])
        start = time.perf_counter()
        for _ in range(200):
            featured.fun(x)
        elapsed = time.perf_counter() - start
        assert counting.count('fun') == 400
        # Generous bound: exponential work would take far longer than this.
        assert elapsed < 60.0


class TestUnrelaxableGate:

    @staticmethod
    def line(x):
        return float(x[0])

    @staticmethod
    def c_line(x):
        return np.array([x[0] - 0.4])

    def make(self, name, **options):
        problem = Problem(self.line, np.array([0.0]), cub=self.c_line)
        return FeaturedProblem(problem, Feature(name, **options), 10, 0)

    def test_gate_uses_predecessor_observed_constraints(self):
        # c(x) = x - 0.4, quantized(mesh=1, ground_truth=False), then the nonlinear gate at x = 0.45:
        # the predecessor observes c(0) = -0.4 and permits the objective; the reference violation is 0.05.
        featured = self.make('quantized+unrelaxable_constraints', mesh_size=1.0, ground_truth=False,
                             unrelaxable_bounds=False, unrelaxable_nonlinear_constraints=True)
        assert featured.fun(np.array([0.45])) == 0.0
        assert featured.maxcv_hist[0] == pytest.approx(0.05)
        assert featured.fun_hist[0] == pytest.approx(0.45)

    def test_order_decides_which_constraints_gate(self):
        options = dict(noise_mode='deterministic', noise_map=lambda x: 1.0, noise_type='absolute', noise_level=1.0,
                       unrelaxable_bounds=False, unrelaxable_nonlinear_constraints=True)
        x = np.array([0.2])
        # noisy then gate: the gate sees c(x) + 1 = 0.8 > 0.
        assert self.make('noisy+unrelaxable_constraints', **options).fun(x) == np.inf
        # gate then noisy: the gate sees c(x) = -0.2, then the objective is noised.
        assert self.make('unrelaxable_constraints+noisy', **options).fun(x) == 1.2

    def test_nan_mixed_with_finite_violations_keeps_legacy_comparison(self):
        def c_nan(x):
            return np.array([np.nan, 0.5])

        problem = Problem(self.line, np.array([0.0]), xl=np.array([-1.0]), xu=np.array([1.0]), cub=c_nan)
        options = dict(unrelaxable_bounds=True, unrelaxable_nonlinear_constraints=True)
        featured = FeaturedProblem(problem, Feature('truncated+unrelaxable_constraints', significant_digits=12,
                                                    **options), 10, 0)
        legacy = FeaturedProblem(problem, Feature('unrelaxable_constraints', **options), 10, 0)
        # A NaN entry makes the nonlinear violation NaN, and NaN > 0 is false: the finite violation
        # 0.5 does not close the gate. This is the established single-feature comparison behaviour.
        assert featured.fun(np.array([0.3])) == 0.3
        assert legacy.fun(np.array([0.3])) == 0.3
        # A violated bound is checked first and still closes the gate.
        assert featured.fun(np.array([1.5])) == np.inf
        assert np.isnan(featured.maxcv_hist[0])

    def test_gate_samples_advance_predecessor_streams_not_the_budget(self):
        counting = Counting()
        featured = FeaturedProblem(counting_problem(counting), Feature('noisy+unrelaxable_constraints',
                                                                       noise_type='absolute', noise_level=1.0,
                                                                       unrelaxable_nonlinear_constraints=True), 10, 0)
        counting.reset()
        featured.fun(X0)
        # The gate samples the noisy constraints once (observed), and the
        # recorded scoring violation reads them once more (reference).
        assert counting.count('cub') == 2 and counting.count('ceq') == 2
        assert featured.n_eval_cub == 0 and featured.n_eval_ceq == 0
        first = featured.cub(X0)
        second = featured.cub(X0)
        assert not np.array_equal(first, second)

    def test_categories_follow_the_predecessor_representation(self):
        problem = Problem(quadratic, X0, xl=XL, xu=XU)
        outside = np.array([9.0, 0.0, 0.0])
        # After a rotation the finite bounds are linear constraints of the predecessor.
        bounds_gate = FeaturedProblem(problem, Feature('linearly_transformed+unrelaxable_constraints',
                                                       condition_factor=0.0, unrelaxable_bounds=True), 10, 1)
        linear_gate = FeaturedProblem(problem, Feature('linearly_transformed+unrelaxable_constraints',
                                                       condition_factor=0.0, unrelaxable_bounds=False,
                                                       unrelaxable_linear_constraints=True), 10, 1)
        assert np.isfinite(bounds_gate.fun(outside))
        assert linear_gate.fun(outside) == np.inf


class TestCustomStage:

    def test_custom_probes_see_distinct_predecessor_samples(self):
        seen = []

        def mod_fun(x, rng, problem):
            seen.append((problem.fun(x), problem.fun(x)))
            assert isinstance(problem, Problem) and problem.n == 3
            return seen[-1][0]

        counting = Counting()
        feature = Feature('noisy+custom', noise_type='absolute', noise_level=1.0, mod_fun=mod_fun)
        featured = FeaturedProblem(counting_problem(counting, 'u'), feature, 10, 0)
        counting.reset()
        featured.fun(X0)
        first, second = seen[0]
        assert first != second
        # One read for the stream derivation, two probes, one reference read.
        assert counting.count('fun') == 4
        assert featured.n_eval_fun == 1
        assert featured.fun_hist[0] == quadratic(X0)

    def test_custom_affine_and_initial_point_then_noise(self):
        d = np.array([2.0, 0.5, 1.0])

        def mod_affine(rng, problem):
            return np.diag(d), np.zeros(3), np.diag(1.0 / d)

        def mod_x0(rng, problem):
            return problem.x0 + 1.0

        feature = Feature('custom+noisy', mod_affine=mod_affine, mod_x0=mod_x0, noise_level=0.0)
        featured = FeaturedProblem(counting_problem(Counting(), 'b'), feature, 10, 0)
        np.testing.assert_array_equal(featured.x0, X0 + 1.0)
        x = np.array([0.4, -0.7, 1.1])
        assert featured.fun(x) == quadratic(d * x)
        assert featured.fun_hist[0] == quadratic(d * x)
        np.testing.assert_allclose(featured.xl, XL / d)

    def test_custom_in_the_middle(self):
        def mod_fun(x, rng, problem):
            return problem.fun(x) + 100.0

        feature = Feature('truncated+custom+noisy', mod_fun=mod_fun, noise_level=0.0, significant_digits=3)
        featured = FeaturedProblem(counting_problem(Counting(), 'u'), feature, 10, 0)
        x = np.array([0.4, -0.7, 1.2])
        assert quadratic(x) == pytest.approx(1.01)
        assert featured.fun(x) == 1.01 + 100.0
        assert featured.fun_hist[0] == quadratic(x)
