"""Saved benchmarks must not fork a parent's active logging threads."""

import json
import multiprocessing as mp
import os
from pathlib import Path
import signal
import subprocess
import sys

import pytest


def context_solver(fun, x0):
    receipt = Path(os.environ['OP_WORKER_RECEIPTS']) / f'{os.getpid()}.json'
    receipt.write_text(json.dumps({'pid': os.getpid(), 'method': mp.get_start_method()}))
    fun(x0)
    return x0


def saved_child(method, folder):
    mp.set_start_method(method, force=True)
    from optiprofiler import benchmark

    receipts = folder / 'workers'
    receipts.mkdir()
    os.environ['OP_WORKER_RECEIPTS'] = str(receipts)
    # Keep stdout buffered: logging starts before workers, which used to make
    # forked workers wait forever for an inherited output lock during exit.
    for _ in range(3):
        benchmark([context_solver, context_solver], plibs=['custom'],
                  problem_names=['custom1'], n_jobs=2, n_runs=1,
                  max_eval_factor=1, max_tol_order=1, run_plain=False,
                  silent=True, draw_hist_plots='parallel', savepath=str(folder))
    workers = [json.loads(path.read_text()) for path in receipts.glob('*.json')]
    assert workers and all(item['pid'] != os.getpid() for item in workers)
    assert all(item['method'] == 'spawn' for item in workers), workers
    assert mp.get_start_method() == method, 'benchmark changed caller process policy'
    assert len(list(folder.rglob('data_for_loading.h5'))) == 3
    assert not mp.active_children(), 'benchmark leaked workers'
    print('SAVED_WORKERS_CLEAN', flush=True)


@pytest.mark.parametrize('method', mp.get_all_start_methods())
def test_saved_workers_use_safe_context_without_changing_caller(method, tmp_path):
    command = [sys.executable, str(Path(__file__).resolve()), method, str(tmp_path)]
    env = dict(os.environ, MPLBACKEND='Agg')
    env.pop('PYTHONUNBUFFERED', None)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, env=env, start_new_session=(os.name != 'nt'))
    try:
        output, _ = proc.communicate(timeout=60)
    except subprocess.TimeoutExpired:
        if os.name != 'nt':
            os.killpg(proc.pid, signal.SIGKILL)
        else:
            subprocess.run(['taskkill', '/PID', str(proc.pid), '/T', '/F'], capture_output=True)
        output, _ = proc.communicate()
        pytest.fail('saved parallel benchmark did not finish:\n' + output)
    assert proc.returncode == 0, output
    assert 'SAVED_WORKERS_CLEAN' in output


@pytest.mark.parametrize('entry', ['stdin', 'package'])
def test_unavailable_entry_point_falls_back_before_executing_solvers(entry, tmp_path):
    source = '''import os
from optiprofiler import benchmark
calls = []
def solver(fun, x0):
    calls.append(os.getpid())
    fun(x0)
    return x0
if __name__ == '__main__':
    benchmark([solver, solver], plibs=['custom'], problem_names=['custom1'],
              n_jobs=2, score_only=True, silent=True, run_plain=False,
              max_eval_factor=1, max_tol_order=1)
    assert calls == [os.getpid(), os.getpid()], calls
    print('ENTRY_POINT_FALLBACK_ONCE')
'''
    if entry == 'package':
        package = tmp_path / 'entry_fixture'
        package.mkdir()
        (package / '__init__.py').write_text('')
        (package / '__main__.py').write_text(source)
        command, input_text = [sys.executable, '-m', 'entry_fixture'], None
    else:
        command, input_text = [sys.executable, '-'], source
    result = subprocess.run(command, input=input_text, cwd=tmp_path,
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'ENTRY_POINT_FALLBACK_ONCE' in result.stdout


if __name__ == '__main__':
    saved_child(sys.argv[1], Path(sys.argv[2]))
