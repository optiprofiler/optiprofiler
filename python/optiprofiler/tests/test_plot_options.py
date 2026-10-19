"""Native plotting properties are observable in public benchmark artifacts."""
import hashlib

import matplotlib
import numpy as np
import pytest
from pypdf import PdfReader

from optiprofiler import Problem, benchmark


def stay(fun, x0):
    fun(x0)
    return x0


def zero(fun, x0):
    point = np.zeros_like(x0)
    fun(point)
    return point


@pytest.mark.parametrize('jobs', [1, 2])
def test_public_profile_pdf_uses_user_font_sizes_without_changing_scores(tmp_path, jobs):
    library = tmp_path / 'style_library'
    library.mkdir()
    (library / 'style_library_tools.py').write_text(
        'from optiprofiler import Problem\n'
        'def style_library_select(options): return ["P"]\n'
        'def style_library_load(name):\n'
        '    return Problem(lambda x: float(x@x), [2.], name=name)\n')
    baseline = benchmark([stay, zero], plibs=['style_library'], custom_problem_libs_path=library,
                         n_jobs=1, n_runs=1, score_only=True, silent=True,
                         max_tol_order=1, run_plain=False)
    old_family = list(matplotlib.rcParams['font.family'])
    styled = benchmark([stay, zero], plibs=['style_library'], custom_problem_libs_path=library,
                       n_jobs=jobs, n_runs=1,
                       draw_hist_plots='parallel' if jobs == 2 else 'sequential', silent=True, max_tol_order=1,
                       run_plain=False, savepath=str(tmp_path),
                       xlabel_performance_profile='X_STYLE_PROBE',
                       ylabel_performance_profile='Y_STYLE_PROBE',
                       plot_options={'text': {'fontsize': 13},
                                     'xlabel': {'fontsize': 21}, 'ylabel': {'fontsize': 22},
                                     'title': {'fontsize': 23}, 'legend': {'fontsize': 18},
                                     'tick_params': {'labelsize': 17},
                                     'rc': {'text.usetex': False, 'font.family': 'monospace'}})
    np.testing.assert_array_equal(styled[0], baseline[0])
    np.testing.assert_array_equal(styled[1], baseline[1])
    pdf = next(tmp_path.rglob('perf_hist_1.pdf'))
    observations = pdf_text(pdf)
    assert ('X_STYLE_PROBE', 21) in observations
    assert ('Y_STYLE_PROBE', 22) in observations
    assert ('stay', 18) in observations
    fonts = []
    PdfReader(pdf).pages[0].extract_text(visitor_text=lambda text, cm, tm, font, size:
                                       fonts.append(str((font or {}).get('/BaseFont', ''))))
    assert any('Mono' in font for font in fonts), fonts
    assert list(matplotlib.rcParams['font.family']) == old_family
    history = next(path for path in tmp_path.rglob('P.pdf') if path.parent.name not in {'raw', 'cummin'})
    history_text = pdf_text(history)
    assert ('stay', 18) in history_text
    assert any('Solving' in text and size == 23 for text, size in history_text)
    summary = next(tmp_path.rglob('summary_*.pdf'))
    assert any('Profiles with' in text and size == 23 for text, size in pdf_text(summary))


def pdf_text(path):
    observations = []
    for page in PdfReader(path).pages:
        page.extract_text(visitor_text=lambda text, cm, tm, font, size:
                          observations.append((text.strip(), size)))
    return observations


def test_load_can_restyle_without_solver_calls_or_archive_changes(tmp_path, monkeypatch):
    library = tmp_path / 'reload_style'
    library.mkdir()
    (library / 'reload_style_tools.py').write_text(
        'from optiprofiler import Problem\n'
        'def reload_style_select(options): return ["P"]\n'
        'def reload_style_load(name): return Problem(lambda x: float(x@x), [2.], name=name)\n')
    options = dict(plibs=['reload_style'], custom_problem_libs_path=library,
                   n_jobs=1, n_runs=1, max_tol_order=1, max_eval_factor=2, run_plain=False,
                   silent=True, savepath=str(tmp_path), benchmark_id='saved',
                   draw_hist_plots='none', plot_options={'rc': {'text.usetex': False}})
    original = benchmark([stay, zero], **options)
    files = [path for path in (tmp_path / 'saved').rglob('*') if path.is_file()]
    hashes = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in files}
    monkeypatch.chdir(tmp_path)
    output = tmp_path / 'restyled'
    loaded = benchmark(None, load='latest', benchmark_id='saved', savepath=str(output),
                       max_tol_order=1, n_jobs=1, silent=True, run_plain=False,
                       xlabel_performance_profile='RELOADED_STYLE',
                       plot_options={'xlabel': {'fontsize': 25}, 'rc': {'text.usetex': False}})
    np.testing.assert_array_equal(loaded[0], original[0])
    np.testing.assert_array_equal(loaded[1], original[1])
    assert all(hashlib.sha256(path.read_bytes()).hexdigest() == digest for path, digest in hashes.items())
    assert ('RELOADED_STYLE', 25) in pdf_text(next(output.rglob('perf_hist_1.pdf')))


@pytest.mark.parametrize('options', [
    {'axess': {}}, {'axes': {'unknown_property': 1}}, {'xlabel': {'fontsize': 'invalid'}},
    {'legend': {'unknown_property': 1}}, {'rc': {'backend': 'MacOSX'}},
    {'legend': {'handles': []}}, {'legend': {'labels': []}},
    {'figure': {'canvas': None}}, {'axes': []},
])
def test_invalid_plot_options_are_rejected_before_solver_execution(options):
    calls = []

    def solver(fun, x0):
        calls.append(1)
        return x0

    with pytest.raises((ValueError, TypeError), match='plot_options'):
        benchmark([solver, solver], problem=Problem(lambda x: 0., [1.]),
                  n_jobs=1, score_only=True, silent=True, plot_options=options)
    assert calls == []
