"""Display-only regressions: extreme histories must export without losing curves."""
import io

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pytest
from matplotlib.figure import Figure

from optiprofiler import plotting
from optiprofiler.utils import ProfileOption


def options(errorbar='minmax'):
    return {
        ProfileOption.ERRORBAR_TYPE: errorbar,
        ProfileOption.HIST_AGGREGATION: 'min',
        ProfileOption.LINE_COLORS: ['tab:blue', 'tab:orange'],
        ProfileOption.LINE_STYLES: ['-', '--'],
        ProfileOption.LINE_WIDTHS: [1, 1],
        ProfileOption.XLABEL_DATA_PROFILE: 'Evaluations / (n + 1)',
    }


def test_nonfinite_replacement_is_run_local():
    original = np.array([[[1, 2, 3], [100, 200, np.inf]]] * 2)
    actual = plotting.process_hist_y_axes(original, [1, 100])
    expected = original.copy()
    expected[:, 1, 2] = 250
    np.testing.assert_array_equal(actual, expected)
    assert np.isinf(original[0, 1, 2])


@pytest.mark.parametrize('values,initial', [
    ([-1e308, np.inf, 1e308], 0),
    ([1, np.inf, 3], np.inf),
    ([np.nan, np.inf, -np.inf], np.nan),
])
def test_nonfinite_replacement_stays_finite(values, initial):
    original = np.tile(values, (2, 2, 1))
    before = original.copy()
    actual = plotting.process_hist_y_axes(original, [initial, initial])
    assert np.isfinite(actual).all()
    np.testing.assert_array_equal(original, before)


@pytest.mark.parametrize('values', [
    [1, 1e290, 2], [1, 1e300, 2], [-1e308, 0, 1e308],
    [np.nan, np.inf, -np.inf], [1e-308, 1e-200, 1e-150],
    [1e308] * 3, [-1e308] * 3, [0] * 3, [1e-308] * 3,
])
@pytest.mark.parametrize('cum', [False, True])
@pytest.mark.parametrize('errorbar', ['minmax', 'meanstd'])
def test_extreme_histories_export_and_remain_in_view(values, cum, errorbar, monkeypatch):
    # Render the real history path, not a stand-alone substitute axis. Disable
    # optional external LaTeX only, so this regression needs no TeX installation.
    monkeypatch.setattr(plotting, 'set_profile_context', lambda _: {'text.usetex': False})
    history = np.tile(values, (2, 2, 1)).astype(float)
    before = history.copy()
    figure = Figure(figsize=(6, 4))
    ax = figure.subplots()
    plotting.draw_hist(history, np.zeros_like(history), history,
                       [values[0]] * 2, [0, 0], [values[0]] * 2,
                       ['a', 'b'], [ax], cum, 'u', 2,
                       np.full((2, 2), 3), options(errorbar), 4)
    output = io.BytesIO()
    figure.savefig(output, format='pdf')
    assert output.getvalue().startswith(b'%PDF')
    low, high = ax.get_ylim()
    assert np.isfinite([low, high]).all() and low < high
    assert np.isfinite(ax.get_yticks()).all()
    for line in ax.lines:
        shown = np.asarray(line.get_ydata())
        assert np.isfinite(shown).all()
        assert np.all((shown >= low) & (shown <= high)), (shown, (low, high))
    if np.any(np.isfinite(history) & (np.abs(history) > 1e100)):
        assert any('display' in text.get_text().lower() for text in ax.texts)
    np.testing.assert_array_equal(history, before)


@pytest.mark.parametrize('n_solvers,n_evals', [(1, 3), (2, 1), (1, 1)])
def test_singleton_history_dimensions(n_solvers, n_evals, monkeypatch):
    monkeypatch.setattr(plotting, 'set_profile_context', lambda _: {'text.usetex': False})
    history = np.ones((n_solvers, 2, n_evals))
    figure = Figure(figsize=(6, 4))
    ax = figure.subplots()
    plotting.draw_hist(history, history, history, [1, 1], [1, 1], [1, 1],
                       ['a', 'b'][:n_solvers], [ax], True, 'u', 2,
                       np.full((n_solvers, 2), n_evals), options(), 4)
    figure.savefig(io.BytesIO(), format='pdf')


@pytest.mark.parametrize('runs,mean,population_std', [
    ([[1, 4, 2]], [1, 4, 2], [0, 0, 0]),
    ([[0, 0, 0], [4, 8, 2]], [2, 4, 1], [2, 4, 1]),
    ([[1, 4, 2], [3, 8, 6], [5, 6, 10]], [3, 6, 6],
     [np.sqrt(8 / 3), np.sqrt(8 / 3), np.sqrt(32 / 3)]),
])
@pytest.mark.parametrize('cum', [False, True])
def test_meanstd_retains_population_normalization(runs, mean, population_std, cum, monkeypatch):
    # Unequal runs distinguish Python's long-standing N normalization from
    # MATLAB's N-1 normalization; duplicated runs cannot detect this difference.
    monkeypatch.setattr(plotting, 'set_profile_context', lambda _: {'text.usetex': False})
    history = np.array([runs], dtype=float)
    before = history.copy()
    mean = np.asarray(mean, dtype=float)
    deviation = np.asarray(population_std, dtype=float)
    expected_shift = np.finfo(float).eps if len(runs) == 2 else 0
    shift = plotting.compute_y_shift(history, options('meanstd'))
    assert shift == expected_shift
    expected = [mean + shift, mean - deviation + shift, mean + deviation + shift]
    if cum:
        expected = [np.minimum.accumulate(values) for values in expected]
    figure = Figure(figsize=(6, 4))
    ax = figure.subplots()
    bands = []
    real_fill_between = ax.fill_between

    def record_band(x, lower, upper, *args, **kwargs):
        bands.append((np.asarray(x), np.asarray(lower), np.asarray(upper)))
        return real_fill_between(x, lower, upper, *args, **kwargs)

    monkeypatch.setattr(ax, 'fill_between', record_band)
    plotting.draw_fun_maxcv_merit_hist(ax, history, ['solver'], cum, 2, shift,
                                     np.full((1, len(runs)), 3), options('meanstd'))
    np.testing.assert_allclose(ax.lines[0].get_ydata(), expected[0])
    if len(runs) > 1:
        assert len(bands) == 1
        np.testing.assert_allclose(bands[0][0], [1 / 3, 2 / 3, 1])
        np.testing.assert_allclose(bands[0][1], expected[1])
        np.testing.assert_allclose(bands[0][2], expected[2])
    else:
        assert not bands
    np.testing.assert_array_equal(history, before)
    figure.savefig(io.BytesIO(), format='pdf')


@pytest.mark.parametrize('aggregation', ['min', 'max', 'mean'])
@pytest.mark.parametrize('cum', [False, True])
@pytest.mark.parametrize('count', [3, 1002, 1003, 2002])
def test_aggregated_evaluation_coordinates_are_one_based(aggregation, cum, count, monkeypatch):
    monkeypatch.setattr(plotting, 'set_profile_context', lambda _: {'text.usetex': False})
    history = np.tile(np.arange(1., count + 1), (2, 2, 1))
    figure = Figure(figsize=(6, 4))
    ax = figure.subplots()
    opts = options()
    opts[ProfileOption.HIST_AGGREGATION] = aggregation
    plotting.draw_hist(history, history, history, [1, 1], [1, 1], [1, 1],
                       ['a', 'b'], [ax], cum, 'u', 2, np.full((2, 2), count), opts, 4)
    for line in ax.lines:
        evaluations = np.rint(line.get_xdata() * 3).astype(int)
        assert evaluations[0] == 1 and evaluations[-1] == count
        assert np.all(np.diff(evaluations) > 0)
        if count <= 1002:
            np.testing.assert_array_equal(evaluations, np.arange(1, count + 1))
        # No-aggregation monotone histories must put their exact value at the
        # matching evaluation; cumulative minimum is identically one here.
        if count <= 1002:
            np.testing.assert_array_equal(line.get_ydata(), np.ones(count) if cum else evaluations)
