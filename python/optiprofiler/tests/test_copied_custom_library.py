"""The bundled custom example remains usable after copying it outside the package."""

from pathlib import Path
import shutil

import numpy as np
import pytest

from optiprofiler import benchmark
from optiprofiler.problem_libs.custom import custom_tools
from optiprofiler.problem_libraries import (
    _load_tools_module, load_problem_library, resolve_problem_library,
)


def copy_custom_library(root, *, initial=17, offset=700, init_file=True):
    folder = root / 'custom'
    shutil.copytree(Path(custom_tools.__file__).parent, folder,
                    ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    if not init_file:
        (folder / '__init__.py').unlink()
    (folder / 'python_problems' / 'custom1.py').write_text(
        'def custom1():\n'
        f'    return dict(name="custom1", x0=[{initial}., 19.], fun=fun)\n'
        f'def fun(x): return {offset}. + sum(x*x)\n', encoding='utf-8')
    return folder


@pytest.mark.parametrize('init_file', [True, False])
def test_copied_template_loads_and_benchmarks_its_external_problem(tmp_path, init_file):
    folder = copy_custom_library(tmp_path, init_file=init_file)
    reference = resolve_problem_library('custom', folder)
    library = load_problem_library(reference)
    assert reference.source == 'custom'
    assert library.select({'ptype': 'u'}, {}) == ['custom1']
    problem = library.load('custom1', {})
    np.testing.assert_array_equal(problem.x0, [17., 19.])
    assert problem.fun(problem.x0) == 1350.

    observed = []

    def solver(fun, x0):
        observed.append((x0.tolist(), float(fun(x0))))
        return x0

    benchmark([solver, solver], plibs=['custom'], custom_problem_libs_path=folder,
              problem_names=['custom1'], n_jobs=1, score_only=True, silent=True,
              draw_hist_plots='none', run_plain=False, max_eval_factor=1,
              max_tol_order=1)
    assert observed == [([17., 19.], 1350.), ([17., 19.], 1350.)]


def test_copied_template_collects_new_external_problem_and_keeps_copies_isolated(tmp_path):
    first = copy_custom_library(tmp_path / 'first')
    second = copy_custom_library(tmp_path / 'second', initial=23, offset=900)
    (first / 'python_problems' / 'external_problem.py').write_text(
        'def external_problem():\n'
        '    return dict(name="external_problem", x0=[31.], fun=fun)\n'
        'def fun(x): return float(sum(x*x))\n', encoding='utf-8')
    first_ref = resolve_problem_library('custom', first)
    info = _load_tools_module(first_ref).custom_get_info()
    assert 'external_problem' in info['name'].tolist()
    first_library = load_problem_library(first_ref)
    assert 'external_problem' in first_library.select({'ptype': 'u'}, {})
    assert first_library.load('external_problem', {}).fun([31.]) == 961.
    second_library = load_problem_library(resolve_problem_library('custom', second))
    np.testing.assert_array_equal(first_library.load('custom1', {}).x0, [17., 19.])
    np.testing.assert_array_equal(second_library.load('custom1', {}).x0, [23., 19.])
    np.testing.assert_array_equal(custom_tools.custom_load('custom1').x0, [1., 1.])
