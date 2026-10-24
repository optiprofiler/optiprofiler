"""Real public workers, with external watchdogs for cancellation regressions."""
import json
import multiprocessing as mp
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time

import pytest


def solver(fun, x0):
    fun(x0)
    return x0


def sleep_child():
    time.sleep(60)


def broken_initializer(*args):
    raise RuntimeError('initializer probe')


def check_child(method, case, folder):
    saving = case.startswith('saved_')
    if saving: case = case[6:]
    mp.set_start_method(method, force=True)
    from optiprofiler import benchmark
    import optiprofiler.profiles as profiles
    from concurrent.futures.process import BrokenProcessPool
    borrowed = mp.Process(target=sleep_child)
    marker = folder / 'started'
    provider = folder / 'aborttoy'; provider.mkdir()
    (provider / 'aborttoy_tools.py').write_text('''import os, time
from pathlib import Path
from optiprofiler import Problem

def aborttoy_select(options): return ['slow', 'bad']
def aborttoy_load(name):
    case = os.environ['OP_ABORT_CASE']
    if name == 'slow' and case != 'ordinary':
        Path(os.environ['OP_ABORT_MARKER']).write_text('ready')
        time.sleep(60)
    if name == 'bad':
        if case == 'exit': raise SystemExit(71)
        if case == 'interrupt': raise KeyboardInterrupt('worker cancellation')
        if case == 'death': os._exit(71)
        if case == 'ordinary': raise ValueError('ordinary load failure')
    return Problem(lambda x: float(sum(x*x)), [1.,2.], name=name)
''')
    os.environ['OP_ABORT_CASE'] = case
    os.environ['OP_ABORT_MARKER'] = str(marker)
    borrowed.start()
    import faulthandler
    faulthandler.dump_traceback_later(20)
    if case == 'initializer':
        profiles.setup_worker_logging = broken_initializer
    timer = None
    if case == 'parent_interrupt':
        def cancel():
            for _ in range(200):
                if marker.exists():
                    os.kill(os.getpid(), signal.SIGINT)
                    return
                time.sleep(.025)
        timer = threading.Thread(target=cancel, daemon=True); timer.start()
    expected = {'exit': SystemExit, 'interrupt': KeyboardInterrupt,
                'death': BrokenProcessPool, 'initializer': BrokenProcessPool,
                'parent_interrupt': KeyboardInterrupt}.get(case)
    before = time.monotonic()
    try:
        try:
            benchmark([solver, solver], plibs=['aborttoy'], custom_problem_libs_path=str(folder),
                      n_jobs=2, n_runs=1, score_only=not saving, silent=True,
                      savepath=str(folder / "output"), draw_hist_plots="none",
                      run_plain=False, max_eval_factor=1, max_tol_order=1)
        except BaseException as err:
            assert expected is not None and isinstance(err, expected), repr(err)
            if case == 'exit': assert err.code == 71
        else:
            assert expected is None
        assert time.monotonic() - before < 20
        assert borrowed.is_alive(), 'caller-owned process was terminated'
        assert mp.active_children() == [borrowed], 'benchmark leaked workers'
        # A fresh benchmark remains usable; never rerun the aborted tasks.
        profiles.setup_worker_logging = __import__('optiprofiler.utils', fromlist=['setup_worker_logging']).setup_worker_logging
        benchmark([solver, solver], plibs=['custom'], problem_names=['custom1'], n_jobs=2,
                  score_only=True, silent=True, run_plain=False,
                  max_eval_factor=1, max_tol_order=1)
        print(json.dumps({'method': method, 'case': case, 'owned_cleanup': True}), flush=True)
    finally:
        faulthandler.cancel_dump_traceback_later()
        borrowed.terminate(); borrowed.join(5)
        if timer is not None: timer.join(1)


@pytest.mark.parametrize('method', mp.get_all_start_methods())
@pytest.mark.parametrize('case', ['ordinary', 'exit', 'interrupt', 'death', 'initializer', 'parent_interrupt', 'saved_exit', 'saved_parent_interrupt'])
def test_public_parallel_abort(method, case, tmp_path):
    if case in ('parent_interrupt', 'saved_parent_interrupt') and os.name == 'nt':
        pytest.skip('POSIX SIGINT injection; Windows Ctrl-C needs console integration coverage')
    command = [sys.executable, str(Path(__file__).resolve()), '--child', method, case, str(tmp_path)]
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, start_new_session=(os.name != 'nt'))
    try:
        output, _ = proc.communicate(timeout=30)
    except subprocess.TimeoutExpired:
        if os.name != 'nt': os.killpg(proc.pid, signal.SIGKILL)
        else: subprocess.run(['taskkill', '/PID', str(proc.pid), '/T', '/F'], capture_output=True)
        output, _ = proc.communicate()
        pytest.fail('parallel abort watchdog expired:\n' + output)
    assert proc.returncode == 0, output
    assert '"owned_cleanup": true' in output


if __name__ == '__main__':
    check_child(sys.argv[2], sys.argv[3], Path(sys.argv[4]))
