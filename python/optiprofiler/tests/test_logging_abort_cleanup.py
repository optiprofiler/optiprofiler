"""Abort-path logging cleanup must finish after a worker died holding the log queue."""
import faulthandler
import logging
import logging.handlers
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import warnings

import pytest


def solver(fun, x0):
    fun(x0)
    return x0


def hold_writer(queue, ready, seconds):
    # A real feeder holds this lock while it writes a record. Pin that window
    # so that terminating the worker inside it is repeatable.
    queue._wlock.acquire()
    ready.set()
    time.sleep(seconds)
    queue._wlock.release()


def check_child(case, folder):
    from optiprofiler._parallel import worker_context
    from optiprofiler.profiles import _close_logging_resources
    from optiprofiler import benchmark

    faulthandler.dump_traceback_later(25, exit=True)
    root = logging.getLogger()
    before = (list(root.handlers), root.level)
    if case in ('control', 'abort'):
        ctx = worker_context()
        queue = ctx.Queue(-1)
        listener = logging.handlers.QueueListener(queue, logging.NullHandler())
        listener.start()
        ready = ctx.Event()
        child = ctx.Process(target=hold_writer, args=(queue, ready, 30 if case == 'abort' else 1))
        child.start()
        assert ready.wait(20), 'the worker never held the write lock'
        if case == 'abort':
            child.terminate()
        child.join(10)
        assert child.exitcode == (-signal.SIGTERM if case == 'abort' else 0), child.exitcode
        started = time.monotonic()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            _close_logging_resources({'listener': listener, 'log_queue': queue}, aborting=True)
        elapsed = time.monotonic() - started
        messages = [str(item.message) for item in caught]
        if case == 'abort':
            assert elapsed < 15, elapsed
            assert any('did not stop within' in message for message in messages), messages
        else:
            # A healthy queue stops promptly: no abandoned listener.
            assert elapsed < 5 and not messages, (elapsed, messages)
    else:
        provider = folder / 'locktoy'
        provider.mkdir()
        (provider / 'locktoy_tools.py').write_text('''import logging, os
from optiprofiler import Problem

def locktoy_select(options): return ['p1', 'p2']
def locktoy_load(name):
    queue = next(h.queue for h in logging.getLogger().handlers if hasattr(h, 'queue'))
    queue._wlock.acquire()
    os._exit(71)
''')
        from concurrent.futures.process import BrokenProcessPool
        started = time.monotonic()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            with pytest.raises(BrokenProcessPool):
                benchmark([solver, solver], plibs=['locktoy'], custom_problem_libs_path=str(folder),
                          n_jobs=2, n_runs=1, silent=True, savepath=str(folder / 'output'),
                          draw_hist_plots='none', max_eval_factor=1, max_tol_order=1)
        assert time.monotonic() - started < 20
        assert any('did not stop within' in str(item.message) for item in caught), \
            [str(item.message) for item in caught]
    assert (list(root.handlers), root.level) == before, 'root logging state was not restored'
    faulthandler.cancel_dump_traceback_later()
    print('LOGGING_CLEANUP_RETURNED', case, flush=True)


@pytest.mark.skipif(os.name == 'nt', reason='Windows queues have no cross-process writer lock')
@pytest.mark.parametrize('case', ['control', 'abort', 'public_worker_death'])
def test_abort_cleanup_is_bounded_after_writer_death(case, tmp_path):
    command = [sys.executable, str(Path(__file__).resolve()), case, str(tmp_path)]
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, start_new_session=True, env=dict(os.environ, MPLBACKEND='Agg'))
    try:
        output, _ = proc.communicate(timeout=40)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        output, _ = proc.communicate()
        pytest.fail('logging cleanup watchdog expired:\n' + output)
    assert proc.returncode == 0, output
    assert f'LOGGING_CLEANUP_RETURNED {case}' in output, output


if __name__ == '__main__':
    check_child(sys.argv[1], Path(sys.argv[2]))
