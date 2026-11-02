"""`silent=True` means that nothing of a solver is printed, so the level is 0."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from optiprofiler import Feature, benchmark
from optiprofiler.benchmark_options import _read_options
from optiprofiler import load_benchmark_options
from optiprofiler.loader import save_options
from optiprofiler.profile_utils import check_validity_profile_options, get_default_profile_options
from optiprofiler.utils import ProfileOption

WARNING = '[WARNING] The option `solver_verbose` will be set to 0 because silent is True.'


def failing(fun, x0):
    fun(x0)
    raise RuntimeError('deliberate solver failure')


def stay(fun, x0):
    fun(x0)
    return x0


def levels_after_checks(given):
    solvers = [stay, stay]
    checked = check_validity_profile_options(solvers, dict(given))
    return get_default_profile_options(solvers, Feature('plain'), checked)[ProfileOption.SOLVER_VERBOSE]


@pytest.mark.parametrize('level', [None, 0, 1, 2])
def test_silent_makes_the_level_zero(level, capsys):
    given = {'silent': True} if level is None else {'silent': True, 'solver_verbose': level}
    assert levels_after_checks(given) == 0
    printed = ' '.join(capsys.readouterr().out.split())
    assert printed == (WARNING if level in (1, 2) else '')


@pytest.mark.parametrize('given, level', [({}, 1), ({'solver_verbose': 0}, 0), ({'solver_verbose': 2}, 2),
                                          ({'silent': False, 'solver_verbose': 1}, 1)])
def test_without_silent_the_level_is_kept(given, level, capsys):
    assert levels_after_checks(given) == level
    assert capsys.readouterr().out == ''


FAILING_RUN = '''import multiprocessing
import os
from pathlib import Path
import sys
from optiprofiler import benchmark

def failing(fun, x0):
    fun(x0)
    Path(os.environ['SOLVER_PROCESS_FILE']).write_text(multiprocessing.current_process().name)
    print('SOLVER STDOUT')
    print('SOLVER STDERR', file=sys.stderr)
    raise RuntimeError('deliberate solver failure')

def stay(fun, x0):
    fun(x0)
    return x0

if __name__ == '__main__':
    level = {} if sys.argv[1] == 'omitted' else {'solver_verbose': int(sys.argv[1])}
    benchmark([failing, stay], plibs=['custom'], problem_names=['custom1'], n_jobs=int(sys.argv[2]), n_runs=1,
              max_eval_factor=1, max_tol_order=1, run_plain=False, draw_hist_plots='none', silent=True,
              savepath=sys.argv[3], benchmark_id='run', report_path=sys.argv[4], **level)
'''


@pytest.mark.parametrize('level, n_jobs', [('omitted', 1), ('0', 1), ('1', 1), ('2', 1), ('omitted', 2), ('2', 2)])
def test_a_failing_solver_prints_nothing_under_silent_and_is_recorded(level, n_jobs, tmp_path):
    script = tmp_path / 'failing_run.py'
    script.write_text(FAILING_RUN)
    report_path = tmp_path / 'report.json'
    process_file = tmp_path / 'solver_process.txt'
    result = subprocess.run([sys.executable, str(script), level, str(n_jobs), str(tmp_path), str(report_path)],
                            cwd=tmp_path, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120,
                            env=dict(os.environ, MPLBACKEND='Agg', SOLVER_PROCESS_FILE=str(process_file)))
    assert result.returncode == 0, result.stdout
    printed = ' '.join(result.stdout.split())
    assert printed == (WARNING if level in ('1', '2') else ''), result.stdout
    # A run falls back to the main process without a word under `silent`.
    in_a_worker = n_jobs > 1 and os.cpu_count() > 1
    assert (process_file.read_text() != 'MainProcess') == in_a_worker

    report = json.loads(report_path.read_text(encoding='utf-8'))
    assert report['configuration']['effective']['profile_options']['solver_verbose'] == 0
    abnormal = [(problem['name'], run['solver_index']) for problem in report['problems'] for run in problem['runs']
                if run['abnormal_termination']]
    assert abnormal == [('custom1', 1)]
    (text,) = tmp_path.glob('run/*/test_log/report.txt')
    section = text.read_text(encoding='utf-8').split('## Solver runs that terminated abnormally')[1].split('## ')[0]
    assert 'custom1' in section


def test_a_silent_record_replays_and_loads_without_printing(tmp_path, monkeypatch, capfd):
    # Earlier builds recorded silent=True with level 1. The caller of a replay
    # supplies the profile options, and a load runs no solver.
    monkeypatch.chdir(tmp_path)
    options = dict(n_jobs=1, max_eval_factor=1, max_tol_order=1, run_plain=False, draw_hist_plots='none', silent=True)
    benchmark([failing, stay], plibs=['custom'], problem_names=['custom1'], n_runs=1, benchmark_id='old', **options)
    path = next(Path(tmp_path).rglob('options_refined.pkl'))
    record = _read_options(path)
    assert record['silent'] is True and record['solver_verbose'] == 0
    record['solver_verbose'] = 1
    save_options(record, path)
    capfd.readouterr()

    replay = load_benchmark_options(_read_options(path))
    assert {'feature', 'n_runs', 'plibs', 'problem_names'} <= set(replay)
    assert 'problem_options' not in replay and 'silent' not in replay
    benchmark([failing, stay],  benchmark_id='replay',
              **replay, **options)
    benchmark(None, load='latest', benchmark_id='old', score_only=True, max_tol_order=1, silent=True,
              draw_hist_plots='none')
    captured = capfd.readouterr()
    assert captured.out == '' and captured.err == ''
