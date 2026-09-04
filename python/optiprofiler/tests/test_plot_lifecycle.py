"""Public benchmark regressions for file-only plotting and saved results."""

import matplotlib
matplotlib.use('Agg')

from pathlib import Path
import hashlib
import os
import subprocess
import sys

import numpy as np
import pytest
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from pypdf import PdfReader

import optiprofiler
from optiprofiler import benchmark


def initial_point_solver(fun, x0):
    """Record a value without introducing a solver dependency."""
    fun(x0)
    return x0


def zero_point_solver(fun, x0):
    """Record the initial value and the value at the origin."""
    fun(x0)
    x = np.zeros_like(x0)
    fun(x)
    return x


def benchmark_options(tmp_path):
    return dict(
        plibs=['custom'],
        problem_names=['custom1'],
        ptype='u',
        mindim=1,
        maxdim=2,
        feature_name='plain',
        n_runs=1,
        n_jobs=1,
        max_eval_factor=3,
        max_tol_order=1,
        draw_hist_plots='none',
        savepath=str(tmp_path),
        benchmark_id='plot-lifecycle',
        solver_names=['initial', 'zero'],
        silent=True,
    )


@pytest.mark.parametrize('line_styles', [None, ['-', '--'], ['solid', 'dashdot'], [(0, (1, 2)), 'dotted']])
def test_score_only_needs_no_figures_or_files(tmp_path, monkeypatch, caplog, line_styles):
    caller_figure, caller_axes = plt.subplots()
    caller_axes.plot([0, 1], [1, 0])
    backend = matplotlib.get_backend()

    def reject_figure_creation(*args, **kwargs):
        raise AssertionError('score_only must not create a Figure')

    monkeypatch.setattr(Figure, '__init__', reject_figure_creation)
    options = benchmark_options(tmp_path)
    options.update(score_only=True, silent=False)
    if line_styles is not None:
        options['line_styles'] = line_styles
    try:
        for _ in range(2):
            scores, profile_scores, curves = benchmark(
                [initial_point_solver, zero_point_solver], **options,
            )
            assert scores.shape == (2,)
            assert np.all(np.isfinite(scores))
            assert len(profile_scores) > 0
            assert curves
            assert plt.fignum_exists(caller_figure.number)
            assert len(caller_axes.lines) == 1
        assert matplotlib.get_backend() == backend
        assert list(tmp_path.iterdir()) == []
        assert 'Failed to save the data' not in caplog.text
    finally:
        plt.close(caller_figure)


@pytest.mark.parametrize('line_styles', [[':', '--'], ['solid', 'dashdot'], [(0, (1, 2)), 'dotted']])
def test_file_plots_do_not_use_pyplot_managers(tmp_path, monkeypatch, line_styles):
    caller_figure, caller_axes = plt.subplots()
    caller_axes.plot([0, 1], [1, 0])
    existing_numbers = plt.get_fignums()
    backend = matplotlib.get_backend()
    font_size = matplotlib.rcParams['font.size']

    def reject_figure_manager(*args, **kwargs):
        raise AssertionError('File export must not create a pyplot manager')

    monkeypatch.setattr(plt, 'new_figure_manager', reject_figure_manager)
    options = benchmark_options(tmp_path)
    options['draw_hist_plots'] = 'sequential'
    options['line_styles'] = line_styles
    try:
        results = benchmark([initial_point_solver, zero_point_solver], **options)
        scores, _, _ = results
        assert np.all(np.isfinite(scores))
        pdfs = list(tmp_path.rglob('*.pdf'))
        assert any(path.name.startswith('summary_') for path in pdfs)
        assert any(path.name == 'custom1.pdf' for path in pdfs)
        for path in pdfs:
            assert len(PdfReader(path).pages) > 0
        assert plt.get_fignums() == existing_numbers
        assert matplotlib.get_backend() == backend
        assert matplotlib.rcParams['font.size'] == font_size
        assert len(caller_axes.lines) == 1

        original_hashes = {
            path: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in tmp_path.rglob('*') if path.is_file()
        }

        def reject_figure_creation(*args, **kwargs):
            raise AssertionError('score_only reload must not create a Figure')

        monkeypatch.setattr(Figure, '__init__', reject_figure_creation)
        monkeypatch.chdir(tmp_path)
        options.pop('solver_names')
        options.update(load='latest', score_only=True)
        loaded_results = benchmark(None, **options)
        np.testing.assert_equal(loaded_results, results)
        assert {
            path: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in tmp_path.rglob('*') if path.is_file()
        } == original_hashes
        assert plt.get_fignums() == existing_numbers
    finally:
        plt.close(caller_figure)


@pytest.mark.parametrize('line_styles', [[], ['bad-style'], ['r--'], [(0, (1, 'bad'))]])
def test_invalid_line_styles_do_not_touch_caller_figure(tmp_path, line_styles):
    caller_figure, caller_axes = plt.subplots()
    caller_axes.plot([0, 1], [1, 0])
    try:
        options = benchmark_options(tmp_path)
        options.update(score_only=True, line_styles=line_styles)
        with pytest.raises(ValueError, match='line styles|empty list'):
            benchmark([initial_point_solver, zero_point_solver], **options)
        assert plt.fignum_exists(caller_figure.number)
        assert len(caller_axes.lines) == 1
        assert list(tmp_path.iterdir()) == []
    finally:
        plt.close(caller_figure)


@pytest.mark.parametrize('load_saved', [False, True])
def test_results_survive_sequential_renderer_exit(tmp_path, load_saved, monkeypatch):
    monkeypatch.chdir(tmp_path)
    options = benchmark_options(tmp_path)
    original_hashes = {}
    if load_saved:
        benchmark([initial_point_solver, zero_point_solver], **options)
        original_hashes = {
            path: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in tmp_path.rglob('data_for_loading.h5')
        }
        assert len(original_hashes) == 1

    # Use a separate process so a native-renderer-style exit cannot kill pytest
    # or leave the interrupted benchmark's logging resources in this process.
    package_file = Path(optiprofiler.__file__).resolve()
    env = os.environ.copy()
    env['PYTHONPATH'] = str(package_file.parent.parent) + os.pathsep + env.get('PYTHONPATH', '')
    completed = subprocess.run(
        [sys.executable, __file__, str(tmp_path), str(int(load_saved)), str(package_file)],
        capture_output=True, text=True, timeout=120, cwd=tmp_path, env=env,
    )
    assert completed.returncode == 73, completed.stdout + completed.stderr
    saved_files = set(tmp_path.rglob('data_for_loading.h5'))
    new_files = saved_files - original_hashes.keys()
    assert len(new_files) == 1
    for path, digest in original_hashes.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest

    # Verify the interrupted experiment through the public load API, with no
    # solver execution and no dependence on the HDF5 implementation details.
    monkeypatch.chdir(next(iter(new_files)).parent.parent)
    options.pop('solver_names')
    options.update(load='latest', score_only=True, benchmark_id='.')
    scores, profile_scores, curves = benchmark(None, **options)
    np.testing.assert_array_equal(scores, [1.0, 0.0] if load_saved else [0.0, 1.0])
    assert len(profile_scores) > 0
    assert curves
    assert set(tmp_path.rglob('data_for_loading.h5')) == saved_files


def _exit_during_history_render(output_path, load_saved):
    original_savefig = Figure.savefig

    def interrupted_savefig(figure, path, *args, **kwargs):
        if isinstance(path, (str, Path)) and 'history_plots' in Path(path).parts:
            os._exit(73)
        return original_savefig(figure, path, *args, **kwargs)

    Figure.savefig = interrupted_savefig
    options = benchmark_options(Path(output_path))
    options['draw_hist_plots'] = 'sequential'
    if load_saved:
        options.pop('solver_names')
        options.update(load='latest', solvers_to_load=[1, 0])
        benchmark(None, **options)
    else:
        benchmark([initial_point_solver, zero_point_solver], **options)


if __name__ == '__main__':
    assert Path(optiprofiler.__file__).resolve() == Path(sys.argv[3])
    _exit_during_history_render(sys.argv[1], bool(int(sys.argv[2])))
