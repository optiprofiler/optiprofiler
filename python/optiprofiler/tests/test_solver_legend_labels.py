"""Solver identifiers are literal labels, not Matplotlib legend directives."""
import io

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pytest
from matplotlib.figure import Figure

from optiprofiler import plotting
from optiprofiler.utils import ProfileOption


def options():
    return {
        ProfileOption.ERRORBAR_TYPE: 'minmax',
        ProfileOption.HIST_AGGREGATION: 'min',
        ProfileOption.LINE_COLORS: ['tab:blue', 'tab:orange'],
        ProfileOption.LINE_STYLES: ['-', '--'],
        ProfileOption.LINE_WIDTHS: [1],
        ProfileOption.SEMILOGX: True,
        ProfileOption.XLABEL_PERFORMANCE_PROFILE: 'Performance ratio',
        ProfileOption.YLABEL_PERFORMANCE_PROFILE: 'Solved fraction',
        ProfileOption.XLABEL_DATA_PROFILE: 'Evaluations / (n + 1)',
        ProfileOption.YLABEL_DATA_PROFILE: 'Solved fraction',
    }


@pytest.mark.parametrize('kind', ['raw', 'cummin', 'performance', 'data'])
@pytest.mark.parametrize('n_solvers', [2, 11])
def test_literal_solver_names_survive_final_pdf(kind, n_solvers, monkeypatch):
    monkeypatch.setattr(plotting, 'set_profile_context', lambda _: {'text.usetex': False})
    # This is a valid Python function name (and an allowed solver_names value),
    # not a request to hide a plotted solver from the legend.
    names = ['_hidden_solver'] + [f'solver_{i}' for i in range(1, n_solvers)]
    figure = Figure(figsize=(8, 5))
    ax = figure.subplots()
    if kind in ('raw', 'cummin'):
        history = np.tile([3., 2., 1.], (n_solvers, 2, 1))
        plotting.draw_hist(history, history, history, [3, 3], [3, 3], [3, 3],
                           names, [ax], kind == 'cummin', 'u', 2,
                           np.full((n_solvers, 2), 3), options(), 4)
    else:
        x = np.tile([0., 1., 2.], (n_solvers, 1)).T
        y = np.tile(np.array([0., .5, 1.])[:, None, None], (1, n_solvers, 2))
        draw = plotting._draw_perf_detail if kind == 'performance' else plotting._draw_data_detail
        draw(ax, x, y, 2., names, options(), '0.1')
    output = io.BytesIO()
    figure.savefig(output, format='pdf', bbox_inches='tight')
    assert output.getvalue().startswith(b'%PDF')
    assert len(ax.lines) == n_solvers
    assert [text.get_text() for text in ax.get_legend().get_texts()] == names


def test_solver_without_evaluations_does_not_create_fake_legend_entry(monkeypatch):
    monkeypatch.setattr(plotting, 'set_profile_context', lambda _: {'text.usetex': False})
    figure = Figure()
    ax = figure.subplots()
    history = np.ones((2, 1, 1))
    plotting.draw_hist(history, history, history, [1], [1], [1],
                       ['_not_run', '_one_evaluation'], [ax], False, 'u', 2,
                       np.array([[0], [1]]), options(), 4)
    figure.savefig(io.BytesIO(), format='pdf')
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ['_one_evaluation']
