"""A run whose solver raised or returned nothing has no output to pass the output-based test."""

import json

import h5py
import numpy as np
import pytest

from optiprofiler import Problem, benchmark

HISTORY, OUTPUT = 0, 1
TOOLS = """
import numpy as np
from optiprofiler import Problem


def objective(x):
    return {objective}


def fallback_library_select(options):
    return ['A', 'B']


def fallback_library_load(name):
    x0 = [2.0] if name == 'A' else [3.0]
    return Problem(objective, x0, name=name, reference={reference})
"""
KINDS = {
    # The initial point meets the target of a reference at f(x0) ...
    'loose_reference': {'objective': 'float(x @ x)',
                        'reference': "{'fun': float(x0[0] ** 2), 'maxcv': 0.0, 'kind': 'merit', 'source': 'f(x0)'}"},
    # ... and the convention for merit_init = +Inf passes its merit of +Inf.
    'infinite_initial_merit': {'objective': "float('inf') if x[0] > 1.0 else float(x @ x)", 'reference': 'None'},
}


def reaches_zero(fun, x0):
    fun(x0)
    fun(np.zeros_like(x0))
    return np.zeros_like(x0)


def crash_before_evaluating(fun, x0):
    raise RuntimeError('crashes before its first evaluation')


def crash_after_evaluating(fun, x0):
    fun(x0)
    fun(np.zeros_like(x0))
    raise RuntimeError('crashes after reaching the best point')


def returns_none(fun, x0):
    fun(x0)


def returns_empty(fun, x0):
    fun(x0)
    return np.empty(0)


def returns_x0(fun, x0):
    fun(x0)
    return x0


def returns_x0_without_evaluating(fun, x0):
    return x0


def crashes_on_problem_b(fun, x0):
    if x0[0] == 3.0:
        raise RuntimeError('crashes on one problem')
    return reaches_zero(fun, x0)


@pytest.fixture
def library(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def options(kind):
        root = tmp_path / 'libraries' / 'fallback_library'
        root.mkdir(parents=True, exist_ok=True)
        (root / 'fallback_library_tools.py').write_text(TOOLS.format(**KINDS[kind]), encoding='utf-8')
        return dict(plibs=['fallback_library'], custom_problem_libs_path=str(root.parent), ptype='u', mindim=1,
                    maxdim=1, n_jobs=1, max_eval_factor=2, max_tol_order=1, silent=True, solver_verbose=0)

    return options


LOAD = dict(load='latest', benchmark_id='saved', score_only=True, silent=True, max_eval_factor=2, max_tol_order=1)
NO_OUTPUT = [crash_before_evaluating, crash_after_evaluating, returns_none, returns_empty]


@pytest.mark.parametrize('kind', sorted(KINDS))
@pytest.mark.parametrize('no_output', NO_OUTPUT)
def test_run_without_an_output_never_passes_the_output_test(library, kind, no_output):
    profiles = benchmark([reaches_zero, no_output], score_only=True, **library(kind))[1]
    np.testing.assert_array_equal(profiles[0, 0, OUTPUT], [1.0, 1.0, 1.0])
    np.testing.assert_array_equal(profiles[1, 0, OUTPUT], [0.0, 0.0, 0.0])


@pytest.mark.parametrize('kind', sorted(KINDS))
def test_solver_that_returns_x0_keeps_passing(library, kind):
    profiles = benchmark([reaches_zero, returns_x0], score_only=True, **library(kind))[1]
    np.testing.assert_array_equal(profiles[1, 0, OUTPUT], [1.0, 1.0, 1.0])


@pytest.mark.parametrize('kind', sorted(KINDS))
def test_output_without_evaluations_is_still_charged_one_unit(library, kind):
    profiles = benchmark([reaches_zero, returns_x0_without_evaluating], score_only=True, **library(kind))[1]
    # One evaluation against two: the genuine output is the fastest.
    np.testing.assert_array_equal(profiles[1, 0, OUTPUT], [1.0, 1.0, 1.0])
    assert np.all(profiles[0, 0, OUTPUT] < 1.0)


def test_only_the_flagged_run_loses_its_output(library, tmp_path):
    seen = set()

    def crashes_on_each_first_run(fun, x0):
        if float(x0[0]) not in seen:
            seen.add(float(x0[0]))
            raise RuntimeError('crashes in the first run of a problem')
        return reaches_zero(fun, x0)

    report_path = tmp_path / 'report.json'
    benchmark([reaches_zero, crashes_on_each_first_run], score_only=True, feature_name='noisy', n_runs=2,
              report_path=str(report_path), **library('loose_reference'))
    report = json.loads(report_path.read_text(encoding='utf-8'))
    convergence = report['profiles']['convergence'][0]
    assert [(entry['hits'], entry['total']) for entry in convergence['output']] == [(4, 4), (2, 4)]
    assert [entry['hits'] for entry in convergence['history']] == [4, 2]
    runs = [run for problem in report['problems'] for run in problem['runs'] if run['solver_index'] == 2]
    assert [run['output_fallback'] for run in runs] == [True, False, True, False]


def test_solver_free_load_reproduces_the_scores(library):
    options = library('loose_reference')
    scores, profiles, _ = benchmark([reaches_zero, crash_after_evaluating], benchmark_id='saved',
                                    draw_hist_plots='none', **options)
    np.testing.assert_array_equal(profiles[1, 0, OUTPUT], [0.0, 0.0, 0.0])
    loaded_scores, loaded_profiles, _ = benchmark(None, **LOAD)
    np.testing.assert_array_equal(loaded_scores, scores)
    np.testing.assert_array_equal(loaded_profiles, profiles)


def test_load_slices_the_flags_with_the_solvers_and_problems(library):
    options = library('loose_reference')
    benchmark([reaches_zero, crash_after_evaluating, crashes_on_problem_b], benchmark_id='saved',
              draw_hist_plots='none', **options)
    # Solver 2 has lost its flags if a stale slice gives it those of solver 1.
    selected = benchmark(None, solvers_to_load=[0, 1], **LOAD)[1]
    np.testing.assert_array_equal(selected[1, 0, OUTPUT], [0.0, 0.0, 0.0])
    selected = benchmark(None, solvers_to_load=[0, 2], **LOAD)[1]
    assert np.all(selected[1, 0, OUTPUT, :2] > 0.0) and np.all(selected[1, 0, OUTPUT, :2] < 1.0)
    # Problem A has no fallback of this solver, whatever the order of the problems.
    selected = benchmark(None, solvers_to_load=[0, 2], problem_names=['A'], **LOAD)[1]
    np.testing.assert_array_equal(selected[1, 0, OUTPUT, :2], [1.0, 1.0])
    selected = benchmark(None, solvers_to_load=[0, 2], problem_names=['B'], **LOAD)[1]
    np.testing.assert_array_equal(selected[1, 0, OUTPUT, :2], [0.0, 0.0])


def test_archive_without_the_flags_keeps_its_behaviour(library, tmp_path):
    options = library('loose_reference')
    profiles = benchmark([reaches_zero, returns_none], benchmark_id='saved', draw_hist_plots='none', **options)[1]
    archives = list((tmp_path / 'saved').rglob('data_for_loading.h5'))
    assert len(archives) == 1
    removed = []

    def remove_flags(name, node):
        if name.endswith('/solver_output_fallbacks'):
            removed.append(name)

    with h5py.File(archives[0], 'r+') as archive:
        archive.visititems(remove_flags)
        for name in removed:
            del archive[name]
    assert removed
    legacy = benchmark(None, **LOAD)[1]
    # Without the flags the recorded value decides, as in an old archive ...
    np.testing.assert_array_equal(legacy[1, 0, OUTPUT], [1.0, 1.0, 1.0])
    # ... and nothing that is history-based depends on them.
    np.testing.assert_array_equal(legacy[:, :, HISTORY], profiles[:, :, HISTORY])


def test_history_based_results_are_the_same_with_and_without_an_output(library):
    options = library('loose_reference')
    without = benchmark([reaches_zero, crash_after_evaluating], score_only=True, **options)
    with_output = benchmark([reaches_zero, reaches_zero], score_only=True, **options)
    np.testing.assert_array_equal(without[0], with_output[0])
    np.testing.assert_array_equal(without[1][:, :, HISTORY], with_output[1][:, :, HISTORY])


def test_report_text_states_the_rule(library, tmp_path):
    benchmark([reaches_zero, crash_before_evaluating], benchmark_id='saved', draw_hist_plots='none',
              **library('loose_reference'))
    reports = list((tmp_path / 'saved').rglob('report.txt'))
    assert len(reports) == 1
    assert 'output fallback (never passes the output-based test)' in reports[0].read_text(encoding='utf-8')


def test_single_problem_score_reads_histories_only():
    problem = Problem(lambda x: float(x @ x), [2.0],
                      reference={'fun': 4.0, 'maxcv': 0.0, 'kind': 'merit', 'source': 'f(x0)'})
    scores = benchmark([reaches_zero, crash_after_evaluating], problem=problem, score_only=True, silent=True,
                       solver_verbose=0, n_jobs=1, max_eval_factor=2)[0]
    np.testing.assert_array_equal(scores[0], scores[1])
