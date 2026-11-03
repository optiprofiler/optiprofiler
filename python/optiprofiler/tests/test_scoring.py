"""Numerical merit, integration and historical default-callback contracts."""

import pickle

import numpy as np
import pytest

from optiprofiler.legacy_compat import loads_trusted
from optiprofiler.scoring import (_default_merit, _default_score_fun,
                                 _default_score_weight_fun, compute_merit_values, integrate)
from optiprofiler.utils import ProfileOption


class TestDefaultMerit:

    def test_feasible(self):
        assert _default_merit(1.0, 0.0, 0.0) == 1.0

    def test_nan_fun(self):
        assert np.isinf(_default_merit(np.nan, 0.0, 0.0))

    def test_nan_maxcv(self):
        assert np.isinf(_default_merit(1.0, np.nan, 0.0))

    def test_small_violation(self):
        result = _default_merit(1.0, 0.005, 0.0)
        assert result > 1.0

    def test_large_violation(self):
        result = _default_merit(1.0, 1.0, 0.0)
        assert np.isinf(result)

    def test_with_nonzero_maxcv_init(self):
        result = _default_merit(1.0, 0.0, 1.0)
        assert result == 1.0


class TestDefaultScoreFunctions:

    def test_default_score_weight_fun(self):
        assert _default_score_weight_fun(0) == 1
        assert _default_score_weight_fun(100) == 1

    def test_default_score_fun(self):
        x = np.random.rand(3, 5, 2, 2)
        result = _default_score_fun(x)
        expected = np.mean(x[:, :, 0, 0], axis=1)
        np.testing.assert_array_equal(result, expected)


class TestComputeMeritValues:

    def test_custom_merit(self):
        def my_merit(f, v, v0):
            return f + v
        fun_values = np.array([1.0, 2.0])
        maxcv_values = np.array([0.1, 0.2])
        maxcv_inits = np.array([0.0, 0.0])
        result = compute_merit_values(my_merit, fun_values, maxcv_values, maxcv_inits)
        np.testing.assert_array_almost_equal(result, [1.1, 2.2])

    def test_scalar_maxcv_init(self):
        fun_values = np.array([1.0, 2.0, 3.0])
        maxcv_values = np.array([0.0, 0.0, 0.0])
        result = compute_merit_values(_default_merit, fun_values, maxcv_values, 0.0)
        np.testing.assert_array_equal(result, fun_values)

    def test_array_maxcv_init(self):
        fun_values = np.array([[1.0, 2.0], [3.0, 4.0]])
        maxcv_values = np.array([[0.0, 0.0], [0.0, 0.0]])
        maxcv_init = np.array([0.0, 0.0])
        result = compute_merit_values(_default_merit, fun_values, maxcv_values, maxcv_init)
        np.testing.assert_array_equal(result, fun_values)

    def test_nan_returns_inf(self):
        fun_values = np.array([np.nan])
        maxcv_values = np.array([0.0])
        result = compute_merit_values(_default_merit, fun_values, maxcv_values, 0.0)
        assert np.isinf(result[0])


class TestIntegrate:

    def test_none_curve(self):
        assert integrate(None, 'perf', {}) == 0.0

    def test_empty_tuple(self):
        assert integrate((), 'perf', {}) == 0.0

    def test_empty_arrays(self):
        assert integrate((np.array([]), np.array([])), 'perf', {}) == 0.0

    def test_perf_step_function(self):
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([0.5, 0.8, 1.0])
        result = integrate((x, y), 'perf', {})
        expected = 0.5 * 1.0 + 0.8 * 1.0
        np.testing.assert_allclose(result, expected)

    def test_data_step_function(self):
        x = np.array([0.0, 1.0, 2.0])
        y = np.array([0.0, 0.5, 1.0])
        result = integrate((x, y), 'data', {})
        expected = 0.0 * 1.0 + 0.5 * 1.0
        np.testing.assert_allclose(result, expected)

    def test_log_ratio(self):
        x = np.array([1, 2, 3])
        y = np.array([0.5, -0.3, 0.2])
        result = integrate((x, y), 'log_ratio', {})
        expected = np.sum(np.abs(y))
        np.testing.assert_allclose(result, expected)

    def test_nonfinite_segments_are_ignored(self):
        result = integrate((np.array([1.0, np.inf]), np.array([0.0, 1.0])), 'perf', {})
        assert result == 0.0
        result = integrate((np.array([1.0]), np.array([np.inf])), 'log_ratio', {})
        assert result == 0.0

    def test_ndarray_input(self):
        curve = np.array([[1.0, 2.0, 3.0], [0.5, 0.8, 1.0]])
        result = integrate(curve, 'perf', {})
        expected = 0.5 * 1.0 + 0.8 * 1.0
        np.testing.assert_allclose(result, expected)

    def test_single_point(self):
        x = np.array([1.0])
        y = np.array([0.5])
        result = integrate((x, y), 'perf', {})
        assert result == 0.0

    def test_custom_weight_fun(self):
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([0.5, 0.8, 1.0])
        opts = {ProfileOption.SCORE_WEIGHT_FUN: lambda t: 2.0}
        result = integrate((x, y), 'perf', opts)
        expected = 2.0 * (0.5 * 1.0 + 0.8 * 1.0)
        np.testing.assert_allclose(result, expected)

    def test_unknown_type_returns_zero(self):
        x = np.array([1.0, 2.0])
        y = np.array([0.5, 0.8])
        result = integrate((x, y), 'unknown', {})
        assert result == 0.0


class TestComputeMeritValuesShapeParity:
    """Cross-language parity: the single-problem loaded-history case takes
    (n_solvers, n_runs, n_evals) values with an (n_runs,) init vector, which is
    the shape the load path slices out of the saved per-problem arrays."""

    def test_single_problem_history_with_run_vector(self):
        fun_values = np.zeros((2, 5, 6))
        maxcv_values = np.zeros_like(fun_values)
        inits = np.arange(1.0, 6.0)
        merits = compute_merit_values(lambda f, cv, cv_init: f + cv + cv_init, fun_values, maxcv_values, inits)
        assert merits.shape == (2, 5, 6)
        np.testing.assert_array_equal(merits[0, :, 0], inits)

    def test_single_problem_out_with_run_vector(self):
        fun_values = np.zeros((2, 5))
        maxcv_values = np.zeros_like(fun_values)
        inits = np.arange(1.0, 6.0)
        merits = compute_merit_values(lambda f, cv, cv_init: f + cv + cv_init, fun_values, maxcv_values, inits)
        assert merits.shape == (2, 5)
        np.testing.assert_array_equal(merits[1, :], inits)


@pytest.mark.parametrize('name', ['_default_merit', '_default_score_fun', '_default_score_weight_fun'])
def test_saved_default_callbacks_resolve_to_the_numerical_implementation(name):
    # Protocol-0 GLOBAL is the exact historical module/name reference, without
    # relying on a current callable's new __module__ when making the fixture.
    historical = b'coptiprofiler.profile_utils\n' + name.encode('ascii') + b'\n.'
    expected = globals()[name]
    assert loads_trusted(historical) is expected
    assert pickle.loads(pickle.dumps(expected)) is expected


@pytest.mark.parametrize('name', ['_default_merit', '_default_score_fun', '_default_score_weight_fun'])
def test_plain_pickle_of_an_old_options_file_finds_the_default_callbacks(name):
    # 1.x and early 2.0 option pickles name the callbacks under profile_utils,
    # and users read them with plain pickle as well as through OptiProfiler.
    # Protocol 0: MARK and DICT, the key, GLOBAL module and name, SETITEM and STOP.
    historical = b''.join([b'(d', b'V' b'merit_fun\n', b'c' b'optiprofiler.profile_utils\n',
                           name.encode('ascii') + b'\n', b's.'])
    assert pickle.loads(historical)['merit_fun'] is globals()[name]
