"""Tests for the plotting module."""
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pytest

from optiprofiler.plotting import (
    _data_formatter,
    _perf_formatter,
    data_ticks,
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
