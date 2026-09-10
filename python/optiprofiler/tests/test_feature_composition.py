"""
Public behaviour of composed features written as ``feature_name='a+b+c'``.

Stage ``a`` is applied first, then ``b``, then ``c``; the final problem is
``c(b(a(P)))``. Every expectation below is a worked value computed by hand
or an independent literal oracle, never a recreation of the implementation.
"""

import numpy as np
import pytest

from optiprofiler.opclasses import Feature, FeaturedProblem, Problem


class CountingProblem:
    """A root problem whose callbacks record every invocation per channel."""

    def __init__(self, fun, x0, **kwargs):
        self.calls = {'fun': [], 'cub': [], 'ceq': []}
        self._fun = fun
        cub = kwargs.pop('cub', None)
        ceq = kwargs.pop('ceq', None)
        if cub is not None:
            kwargs['cub'] = self._wrap('cub', cub)
        if ceq is not None:
            kwargs['ceq'] = self._wrap('ceq', ceq)
        self.problem = Problem(self._wrap('fun', fun), x0, **kwargs)
        self.reset()

    def _wrap(self, channel, function):
        def wrapped(x):
            self.calls[channel].append(np.array(x, dtype=float))
            return function(x)
        return wrapped

    def reset(self):
        for channel in self.calls:
            self.calls[channel] = []


def rosenbrock(x):
    return float(np.sum(100.0 * (x[1:] - x[:-1] ** 2) ** 2 + (1.0 - x[:-1]) ** 2))


def constant_1_23456789(x):
    return 1.23456789


def noise_map_0_123456(x):
    return 0.123456


class TestCompositeConstruction:

    def test_name_normalization_and_plain_elimination(self):
        assert Feature(' Noisy + Truncated ').name == 'noisy+truncated'
        assert Feature('noisy+plain').name == 'noisy'
        assert Feature('plain+plain').name == 'plain'
        assert Feature('PLAIN+noisy+plain+truncated').name == 'noisy+truncated'
        assert isinstance(Feature('noisy+truncated'), Feature)

    @pytest.mark.parametrize('name', ['noisy+', '+noisy', 'noisy++truncated', '+', ' + '])
    def test_empty_tokens_are_rejected(self, name):
        with pytest.raises(ValueError, match='empty'):
            Feature(name)

    def test_unknown_tokens_are_rejected(self):
        with pytest.raises(ValueError, match='Unknown feature: unknown.'):
            Feature('noisy+unknown')

    def test_effective_single_stage_is_the_legacy_feature(self):
        x0 = np.array([0.3, -1.2, 2.5])
        legacy = FeaturedProblem(Problem(rosenbrock, x0), Feature('noisy'), 10, 7)
        declared = FeaturedProblem(Problem(rosenbrock, x0), Feature('plain+noisy+plain'), 10, 7)
        for _ in range(3):
            assert declared.fun(x0) == legacy.fun(x0)
        assert declared.fun(x0 + 0.5) == legacy.fun(x0 + 0.5)
        np.testing.assert_array_equal(declared.fun_hist, legacy.fun_hist)


class TestOptionRouting:

    def test_supplied_options_are_broadcast_only_to_owners(self):
        assert Feature('noisy+perturbed_x0').options == {'n_runs': 5}
        feature = Feature('noisy+perturbed_x0', distribution='gaussian', noise_level=0.2)
        assert feature.options == {'distribution': 'gaussian', 'noise_level': 0.2, 'n_runs': 5}
        assert Feature('noisy+noisy', noise_level=0.2).options == {'noise_level': 0.2, 'n_runs': 5}

    def test_stage_specific_validation_error(self):
        with pytest.raises(ValueError) as excinfo:
            Feature('noisy+perturbed_x0', distribution='uniform')
        message = str(excinfo.value)
        assert "perturbed_x0" in message and 'stage 2' in message and '"gaussian" or "spherical"' in message

    def test_option_owned_by_no_stage_is_rejected(self):
        with pytest.raises(ValueError, match="Option `nan_rate` is not valid for feature 'noisy\\+truncated'"):
            Feature('noisy+truncated', nan_rate=0.1)
        with pytest.raises(ValueError, match='Unknown option for feature: foo.'):
            Feature('noisy+truncated', foo=1)


class TestGlobalRuns:

    def test_default_is_maximum_of_child_legacy_defaults(self):
        assert Feature('noisy+truncated').options['n_runs'] == 5
        assert Feature('truncated+quantized').options['n_runs'] == 1
        assert Feature('noisy+truncated', noise_mode='deterministic').options['n_runs'] == 1
        assert Feature('noisy+perturbed_x0', noise_mode='deterministic').options['n_runs'] == 5

    def test_legacy_default_not_actual_stochasticity(self):
        feature = Feature('linearly_transformed+truncated', rotated=False)
        assert feature.options['n_runs'] == 5
        assert feature.is_stochastic is False

    def test_stochasticity_is_or_of_children(self):
        assert Feature('truncated+quantized').is_stochastic is False
        assert Feature('truncated+quantized', perturbed_trailing_digits=True).is_stochastic is True
        assert Feature('noisy+truncated', noise_mode='deterministic').is_stochastic is False

    def test_explicit_n_runs_is_global_and_validated(self):
        assert Feature('noisy+perturbed_x0', n_runs=3).options['n_runs'] == 3
        assert Feature('noisy+perturbed_x0', n_runs=2.0).options['n_runs'] == 2
        with pytest.raises(TypeError):
            Feature('noisy+perturbed_x0', n_runs=1.5)
        with pytest.raises(ValueError):
            Feature('noisy+perturbed_x0', n_runs=0)


class TestOrderedValueComposition:
    """noisy (deterministic map 0.123456, absolute) and truncated (6 digits)."""

    OPTIONS = dict(noise_mode='deterministic', noise_map=noise_map_0_123456, noise_type='absolute',
                   noise_level=1.0, significant_digits=6)

    def make(self, name):
        counting = CountingProblem(constant_1_23456789, np.array([0.7]))
        featured = FeaturedProblem(counting.problem, Feature(name, **self.OPTIONS), 5, 0)
        counting.reset()
        return counting, featured

    def test_noisy_then_truncated(self):
        # 1.23456789 + 0.123456 = 1.35802389, then six significant digits.
        counting, featured = self.make('noisy+truncated')
        assert featured.fun(np.array([0.7])) == 1.35802
        assert featured.fun_hist.tolist() == [1.23456789]

    def test_truncated_then_noisy(self):
        # 1.23457 after six significant digits, then + 0.123456.
        counting, featured = self.make('truncated+noisy')
        assert featured.fun(np.array([0.7])) == pytest.approx(1.358026, abs=1e-15)
        assert featured.fun_hist.tolist() == [1.23456789]

    def test_plain_insertion_is_neutral(self):
        _, a = self.make('noisy+truncated')
        _, b = self.make('plain+noisy+plain+truncated+plain')
        assert a.fun(np.array([0.7])) == b.fun(np.array([0.7]))

    def test_one_observed_and_one_reference_root_read_per_query(self):
        counting, featured = self.make('noisy+truncated')
        featured.fun(np.array([0.7]))
        featured.fun(np.array([0.9]))
        assert len(counting.calls['fun']) == 4
        assert featured.n_eval_fun == 2

    def test_budget_and_termination_do_not_depend_on_composition(self):
        counting, featured = self.make('noisy+truncated')
        x = np.array([0.7])
        outputs = [featured.fun(x) for _ in range(10)]
        assert outputs == [1.35802] * 10
        assert featured.n_eval_fun == 5
        with pytest.raises(StopIteration):
            featured.fun(x)
