"""Action tests must not count an unavailable solver as exercised coverage."""

import json

import numpy as np
import pytest

from optiprofiler.action_tests import solvers


def initial_solver(fun, x0):
    fun(x0)
    return x0


def test_old_scipy_omits_cobyqa_but_retains_two_supported_solvers(monkeypatch):
    def unavailable(**kwargs):
        raise ValueError('Unknown method cobyqa')

    monkeypatch.setattr(solvers, 'show_options', unavailable)
    for unconstrained in (False, True):
        selected, names, skipped = solvers.get_solver_configuration(unconstrained)
        assert len(selected) == len(names) == 2
        assert 'COBYQA' not in names
        assert 'COBYLA' in names and 'Nelder-Mead' in names
        assert skipped == ['COBYQA']


def test_supported_scipy_keeps_cobyqa(monkeypatch):
    monkeypatch.setattr(solvers, 'show_options', lambda **kwargs: 'COBYQA options')
    selected, names, skipped = solvers.get_solver_configuration()
    assert solvers.scipy_cobyqa in selected
    assert names == ['COBYLA', 'COBYQA', 'Nelder-Mead']
    assert not skipped


@pytest.mark.parametrize('case', ['unknown_method', 'no_evaluations'])
def test_preflight_rejects_solver_before_benchmark_can_swallow_failure(case):
    def broken(fun, x0):
        if case == 'unknown_method':
            raise ValueError('Unknown solver COBYQA')
        return x0

    with pytest.raises(RuntimeError, match='COBYQA'):
        solvers.verify_solver_support([broken], ['COBYQA'])


def test_coverage_is_required_for_each_solver_and_preserves_budget_stops():
    report = {'problems': [{'role': 'primary', 'status': 'completed', 'runs': [
        {'solver_index': 1, 'evaluations': 2, 'abnormal_termination': True},
        {'solver_index': 2, 'evaluations': 0, 'abnormal_termination': True},
    ]}]}
    with pytest.raises(RuntimeError, match='COBYQA.*no recorded objective evaluations'):
        solvers.verify_benchmark_coverage(report, ['COBYLA', 'COBYQA'])
    report['problems'][0]['runs'][1]['evaluations'] = 1
    solvers.verify_benchmark_coverage(report, ['COBYLA', 'COBYQA'])


def test_checked_benchmark_reads_real_histories_and_reports_skips(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(solvers, 'SKIPPED_SOLVERS', ['COBYQA'])
    target = tmp_path / 'action_report.json'
    result = solvers.checked_benchmark(
        [initial_solver, initial_solver], solver_names=['first', 'second'],
        plibs=['custom'], problem_names=['custom1'], n_jobs=1, score_only=True,
        silent=True, draw_hist_plots='none', run_plain=False, max_eval_factor=1,
        max_tol_order=1, report_path=target)
    assert np.all(np.isfinite(result[0]))
    report = json.loads(target.read_text(encoding='utf-8'))
    runs = report['problems'][0]['runs']
    assert {run['solver_index'] for run in runs if run['evaluations'] > 0} == {1, 2}
    output = capsys.readouterr().out
    assert 'SKIP COBYQA' in output
    assert 'first: 1 evaluated runs' in output
    assert 'second: 1 evaluated runs' in output
