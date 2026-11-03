"""A profile branch in which no solver meets the convergence test scores 0, as in MATLAB."""

import json

import numpy as np
import pytest

from optiprofiler import benchmark
from optiprofiler.tests.eval_report_contract import assert_valid

HISTORY, OUTPUT = 0, 1
TOOLS = """
from optiprofiler import Problem


def objective(x):
    return float(x @ x)


def undrawable_library_select(options):
    return ['A', 'B']


def undrawable_library_load(name):
    x0 = [2.0] if name == 'A' else [3.0]
    return Problem(objective, x0, name=name, reference={'fun': 0.0, 'maxcv': 0.0, 'kind': 'feasible', 'source': 'analytic'})
"""


def to_three_tenths(fun, x0):
    fun(x0)
    fun(0.3 * x0)
    return 0.3 * x0


def to_a_quarter(fun, x0):
    fun(x0)
    fun(x0 / 4)
    return x0 / 4


def crash_before_evaluating(fun, x0):
    raise RuntimeError('crashes before its first evaluation')


@pytest.fixture
def options(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    root = tmp_path / 'libraries' / 'undrawable_library'
    root.mkdir(parents=True)
    (root / 'undrawable_library_tools.py').write_text(TOOLS, encoding='utf-8')
    return dict(plibs=['undrawable_library'], custom_problem_libs_path=str(root.parent), ptype='u', mindim=1,
                maxdim=1, n_jobs=1, max_eval_factor=2, max_tol_order=3, score_only=True, silent=True,
                solver_verbose=0, report_path=str(tmp_path / 'report.json'))


def test_tolerances_no_solver_reaches_score_zero_in_every_profile(options):
    # Both solvers reach 9 % and 6 % of f(x0) against the reference 0, so
    # only the loosest tolerance, 1e-1, has a solver that meets the test.
    _, profiles, curves = benchmark([to_three_tenths, to_a_quarter], **options)
    assert np.all(profiles[:, 0, HISTORY, :2] > 0)
    np.testing.assert_array_equal(profiles[:, 1:], 0.0)
    for curve in curves[1:]:
        for branch in ('hist', 'out'):
            assert curve[branch]['log_ratio'] == [None, None]
    report = json.loads(open(options['report_path'], encoding='utf-8').read())
    assert_valid(report, 'eval_report-v2.schema.json')
    assert report['profiles']['plot_refs']
    assert all(ref.startswith('profile-1-') for ref in report['profiles']['plot_refs'])


def test_solvers_that_all_crash_score_zero_in_every_profile(options):
    _, profiles, curves = benchmark([crash_before_evaluating, crash_before_evaluating], **options)
    np.testing.assert_array_equal(profiles, 0.0)
    report = json.loads(open(options['report_path'], encoding='utf-8').read())
    assert_valid(report, 'eval_report-v2.schema.json')
    assert report['profiles']['plot_refs'] == []
