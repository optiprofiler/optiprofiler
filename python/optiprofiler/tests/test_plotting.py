"""Tests for the plotting module."""
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pytest
from matplotlib import pyplot as plt

from optiprofiler.plotting import (
    _data_formatter,
    _draw_perf_detail,
    _draw_data_detail,
    _get_extended_performances_data_profile_axes,
    _get_performance_data_profile_axes,
    _perf_formatter,
    data_ticks,
    draw_profiles,
    format_float_scientific_latex,
    perf_ticks,
    set_profile_context,
    compute_y_shift,
    process_hist_y_axes,
)
from optiprofiler.utils import ProfileOption


class TestPerfFormatter:

    def test_integer_power(self):
        result = _perf_formatter(0.0, None)
        assert result == '1'

    def test_integer_power_2(self):
        result = _perf_formatter(3.0, None)
        assert result == '8'

    def test_non_integer_power(self):
        result = _perf_formatter(1.5, None)
        assert '$2^{' in result
        assert '}$' in result


class TestDataFormatter:

    def test_integer_power(self):
        result = _data_formatter(0.0, None)
        assert result == '0'

    def test_integer_power_2(self):
        result = _data_formatter(1.0, None)
        assert result == '1'

    def test_integer_power_3(self):
        result = _data_formatter(3.0, None)
        assert result == '7'

    def test_non_integer_power(self):
        result = _data_formatter(1.5, None)
        assert '$2^{' in result
        assert '}-1$' in result


class TestPerfTicks:

    def test_semilogx_large(self):
        ticks, labels = perf_ticks(10.0, True)
        assert len(ticks) > 0
        assert len(labels) == len(ticks)
        assert ticks[0] == 0
        assert labels[0] == '1'

    def test_semilogx_medium(self):
        ticks, labels = perf_ticks(3.0, True)
        assert len(ticks) > 0
        assert ticks[0] == 0

    def test_semilogx_small(self):
        ticks, labels = perf_ticks(0.5, True)
        assert len(ticks) == 2
        assert ticks[0] == 0

    def test_semilogx_very_small(self):
        ticks, labels = perf_ticks(1e-4, True)
        assert len(ticks) == 1
        assert ticks[0] == 0

    def test_linear_large(self):
        ticks, labels = perf_ticks(10.0, False)
        assert len(ticks) > 0
        assert ticks[0] == 1

    def test_linear_medium(self):
        ticks, labels = perf_ticks(3.0, False)
        assert len(ticks) > 0

    def test_linear_small(self):
        ticks, labels = perf_ticks(1.5, False)
        assert len(ticks) == 2

    def test_linear_very_small(self):
        ticks, labels = perf_ticks(1.0, False)
        assert len(ticks) == 1


class TestDataTicks:

    def test_semilogx_large(self):
        ticks, labels = data_ticks(10.0, True)
        assert len(ticks) > 0
        assert ticks[0] == 0
        assert labels[0] == '0'

    def test_semilogx_medium(self):
        ticks, labels = data_ticks(3.0, True)
        assert len(ticks) > 0

    def test_semilogx_small(self):
        ticks, labels = data_ticks(0.5, True)
        assert len(ticks) == 2

    def test_semilogx_very_small(self):
        ticks, labels = data_ticks(0.05, True)
        assert len(ticks) == 1

    def test_linear_large(self):
        ticks, labels = data_ticks(10.0, False)
        assert len(ticks) > 0

    def test_linear_medium(self):
        ticks, labels = data_ticks(3.0, False)
        assert len(ticks) > 0

    def test_linear_small(self):
        ticks, labels = data_ticks(0.5, False)
        assert len(ticks) == 2

    def test_linear_very_small(self):
        ticks, labels = data_ticks(0.05, False)
        assert len(ticks) == 1


class TestSetProfileContext:

    def test_returns_dict(self):
        ctx = set_profile_context({})
        assert isinstance(ctx, dict)
        assert 'font.family' in ctx
        assert ctx['font.family'] == 'serif'
        assert ctx['font.size'] == 16
        assert 'text.usetex' in ctx


class TestFormatFloatScientificLatex:

    def test_zero(self):
        raw, formatted = format_float_scientific_latex(0)
        assert raw == '0'
        assert formatted == '0'

    def test_one(self):
        raw, formatted = format_float_scientific_latex(1.0)
        assert '1' in raw

    def test_large_number(self):
        raw, formatted = format_float_scientific_latex(1e10)
        assert '10' in formatted

    def test_small_number(self):
        raw, formatted = format_float_scientific_latex(1e-5)
        assert '10' in formatted
        assert '-5' in formatted

    def test_negative(self):
        raw, formatted = format_float_scientific_latex(-3.14e5)
        assert '-' in raw
        assert '-' in formatted

    def test_custom_digits(self):
        raw, formatted = format_float_scientific_latex(1.23456e3, digits=2)
        assert isinstance(raw, str)
        assert isinstance(formatted, str)

    def test_coefficient_one(self):
        raw, formatted = format_float_scientific_latex(1e3)
        assert '10^{3}' in formatted

    def test_exponent_zero(self):
        raw, formatted = format_float_scientific_latex(5.0, digits=2)
        assert '10' not in formatted


class TestComputeYShift:

    def test_basic(self):
        history = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        opts = {ProfileOption.ERRORBAR_TYPE: 'minmax'}
        shift = compute_y_shift(history, opts)
        assert isinstance(shift, (int, float, np.floating))

    def test_constant_history(self):
        history = np.array([3.0, 3.0, 3.0])
        opts = {ProfileOption.ERRORBAR_TYPE: 'minmax'}
        shift = compute_y_shift(history, opts)
        assert isinstance(shift, (int, float, np.floating))
        assert shift == 0

    def test_with_negative_values(self):
        history = np.array([-1.0, 0.0, 1.0])
        opts = {ProfileOption.ERRORBAR_TYPE: 'minmax'}
        shift = compute_y_shift(history, opts)
        assert shift > 0

    def test_meanstd_mode(self):
        history = np.array([[1.0, 2.0], [3.0, 4.0]])
        opts = {ProfileOption.ERRORBAR_TYPE: 'meanstd'}
        shift = compute_y_shift(history, opts)
        assert isinstance(shift, (int, float, np.floating))


class TestProcessHistYAxes:

    def test_basic(self):
        histories = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        init = 0.5
        result = process_hist_y_axes(histories, init)
        assert isinstance(result, np.ndarray)

    def test_with_nan(self):
        histories = np.array([[1.0, np.nan, 3.0], [4.0, 5.0, 6.0]])
        init = 0.5
        result = process_hist_y_axes(histories, init)
        assert isinstance(result, np.ndarray)


class TestGetPerformanceDataProfileAxes:

    def test_basic(self):
        work = np.array([[[10.0], [20.0]], [[5.0], [15.0]]])
        denominator = lambda i_problem, i_run: np.nanmin(work[i_problem, :, i_run])
        x, y, ratio_max = _get_performance_data_profile_axes(work, denominator)
        assert np.isfinite(ratio_max)
        assert ratio_max > 0

    def test_all_nan(self):
        work = np.full((3, 2, 1), np.nan)
        denominator = lambda i_problem, i_run: np.nanmin(work[i_problem, :, i_run], initial=np.inf)
        x, y, ratio_max = _get_performance_data_profile_axes(work, denominator)
        assert ratio_max == np.finfo(float).eps

    def test_single_problem(self):
        work = np.array([[[5.0], [10.0]]])
        denominator = lambda i_problem, i_run: np.nanmin(work[i_problem, :, i_run])
        x, y, ratio_max = _get_performance_data_profile_axes(work, denominator)
        assert ratio_max == 2.0


class TestGetExtendedPerformancesDataProfileAxes:

    def _make_profile_options(self):
        return {ProfileOption.SEMILOGX: True}

    def _make_curves(self, n_solvers, n_runs):
        return {
            'perf': [[None] * (n_runs + 1) for _ in range(n_solvers)],
            'data': [[None] * (n_runs + 1) for _ in range(n_solvers)],
            'log_ratio': [None] * n_solvers,
        }

    def test_basic(self):
        work = np.array([[[10.0], [20.0]], [[5.0], [15.0]]])
        problem_dims = np.array([2, 3])
        curves = self._make_curves(2, 1)
        x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data, curves = \
            _get_extended_performances_data_profile_axes(work, problem_dims, self._make_profile_options(), curves)
        assert np.isfinite(ratio_max_perf)
        assert np.isfinite(ratio_max_data)
        assert x_perf.shape[1] == 2
        assert x_data.shape[1] == 2

    def test_semilogx_false(self):
        work = np.array([[[10.0], [20.0]]])
        problem_dims = np.array([2])
        opts = {ProfileOption.SEMILOGX: False}
        curves = self._make_curves(2, 1)
        x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data, curves = \
            _get_extended_performances_data_profile_axes(work, problem_dims, opts, curves)
        assert np.isfinite(ratio_max_perf)


class TestDrawProfiles:

    def _make_profile_options(self):
        return {
            ProfileOption.SEMILOGX: True,
            ProfileOption.SCORE_ONLY: True,
            ProfileOption.LINE_COLORS: ['#1f77b4', '#ff7f0e'],
            ProfileOption.LINE_STYLES: ['-', '--'],
            ProfileOption.LINE_WIDTHS: [1.5],
            ProfileOption.ERRORBAR_TYPE: 'minmax',
            ProfileOption.XLABEL_PERFORMANCE_PROFILE: '',
            ProfileOption.YLABEL_PERFORMANCE_PROFILE: '',
            ProfileOption.XLABEL_DATA_PROFILE: '',
            ProfileOption.YLABEL_DATA_PROFILE: '',
            ProfileOption.XLABEL_LOG_RATIO_PROFILE: '',
            ProfileOption.YLABEL_LOG_RATIO_PROFILE: '',
        }

    def test_score_only(self):
        n_problems, n_solvers, n_runs = 3, 2, 1
        work = np.random.rand(n_problems, n_solvers, n_runs) * 100 + 1
        problem_dims = np.array([2, 3, 5])
        curves = {
            'perf': [[None] * (n_runs + 1) for _ in range(n_solvers)],
            'data': [[None] * (n_runs + 1) for _ in range(n_solvers)],
            'log_ratio': [None] * n_solvers,
        }
        fig_perf, fig_data, fig_log_ratio, curves = draw_profiles(
            work, problem_dims, ['s1', 's2'], '0.1', 0,
            None, None, None, False, True, True, False,
            self._make_profile_options(), curves,
        )
        assert fig_perf is not None
        assert fig_data is not None
        assert curves['perf'][0][0] is not None
        assert curves['data'][0][0] is not None
        plt.close('all')

    def test_with_nan_work(self):
        n_problems, n_solvers, n_runs = 2, 2, 1
        work = np.full((n_problems, n_solvers, n_runs), np.nan)
        work[0, 0, 0] = 5.0
        problem_dims = np.array([2, 3])
        curves = {
            'perf': [[None] * (n_runs + 1) for _ in range(n_solvers)],
            'data': [[None] * (n_runs + 1) for _ in range(n_solvers)],
            'log_ratio': [None] * n_solvers,
        }
        fig_perf, fig_data, fig_log_ratio, curves = draw_profiles(
            work, problem_dims, ['s1', 's2'], '0.1', 0,
            None, None, None, False, True, True, False,
            self._make_profile_options(), curves,
        )
        assert fig_perf is not None
        plt.close('all')


class TestDrawPerfDataDetail:

    def _make_profile_options(self):
        return {
            ProfileOption.SEMILOGX: True,
            ProfileOption.LINE_COLORS: ['#1f77b4', '#ff7f0e'],
            ProfileOption.LINE_STYLES: ['-', '--'],
            ProfileOption.LINE_WIDTHS: [1.5],
            ProfileOption.ERRORBAR_TYPE: 'minmax',
            ProfileOption.XLABEL_PERFORMANCE_PROFILE: 'Perf ratio',
            ProfileOption.YLABEL_PERFORMANCE_PROFILE: 'Perf profiles (tol = %s)',
            ProfileOption.XLABEL_DATA_PROFILE: 'Simplex gradients',
            ProfileOption.YLABEL_DATA_PROFILE: 'Data profiles (tol = %s)',
        }

    def test_draw_perf_detail_finite(self):
        fig, ax = plt.subplots()
        x = np.array([[0, 1, 2, 3], [0, 1, 2, 3]], dtype=float).T
        y = np.array([[[0, 0.5, 0.8, 1.0]], [[0, 0.3, 0.6, 0.9]]], dtype=float).T
        y = y.reshape(4, 2, 1)
        _draw_perf_detail(ax, x, y, 3.0, ['s1', 's2'], self._make_profile_options(), '0.1')
        plt.close(fig)

    def test_draw_perf_detail_inf_ratio(self):
        fig, ax = plt.subplots()
        x = np.array([[0, np.inf], [0, np.inf]], dtype=float).T
        y = np.array([[[0, 0]], [[0, 0]]], dtype=float).T
        y = y.reshape(2, 2, 1)
        _draw_perf_detail(ax, x, y, np.inf, ['s1', 's2'], self._make_profile_options(), '0.1')
        plt.close(fig)

    def test_draw_data_detail_finite(self):
        fig, ax = plt.subplots()
        x = np.array([[0, 1, 2], [0, 1, 2]], dtype=float).T
        y = np.array([[[0, 0.5, 1.0]], [[0, 0.3, 0.8]]], dtype=float).T
        y = y.reshape(3, 2, 1)
        _draw_data_detail(ax, x, y, 2.0, ['s1', 's2'], self._make_profile_options(), '0.1')
        plt.close(fig)

    def test_draw_data_detail_inf_ratio(self):
        fig, ax = plt.subplots()
        x = np.array([[0, np.inf], [0, np.inf]], dtype=float).T
        y = np.array([[[0, 0]], [[0, 0]]], dtype=float).T
        y = y.reshape(2, 2, 1)
        _draw_data_detail(ax, x, y, np.inf, ['s1', 's2'], self._make_profile_options(), '0.1')
        plt.close(fig)
