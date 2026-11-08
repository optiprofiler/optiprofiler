"""A pool whose worker died writing its result cannot hang benchmark after a failure."""

import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pytest

import optiprofiler

PYTHON_ROOT = Path(optiprofiler.__file__).resolve().parents[1]

# A worker writes the header of a large result and part of its body, then
# stops, as a worker killed while writing leaves the pipe. Meanwhile another
# task fails, and its exception is slow to rebuild in the parent, so the
# manager thread is inside recv() when run_tasks terminates the workers.
CHILD = textwrap.dedent('''
    import faulthandler
    import inspect
    import struct
    import sys
    import time

    sys.path.insert(0, sys.argv[1])
    from optiprofiler import _parallel


    class SlowToRebuild(Exception):
        def __reduce__(self):
            return (_rebuild_slowly, (str(self),))


    def _rebuild_slowly(message):
        time.sleep(1.5)
        return SlowToRebuild(message)


    def task(kind, payload=b''):
        if kind == 'fail':
            time.sleep(0.3)
            raise SlowToRebuild('deliberate task failure')
        if kind == 'sleep':
            time.sleep(120)
        time.sleep(0.6)
        result_queue = inspect.stack(0)[1].frame.f_locals['result_queue']
        lock = getattr(result_queue, '_wlock', None)
        if lock is not None:
            lock.acquire()
        result_queue._writer._send(struct.pack('!i', 1 << 20))
        result_queue._writer._send(b'x' * 1000)
        time.sleep(120)


    if __name__ == '__main__':
        faulthandler.dump_traceback_later(float(sys.argv[3]), exit=True)
        pad = b'p' * 200_000 if sys.argv[2] == 'large_arguments' else b''
        args = [('fail', pad), ('partial', pad)] + ([('sleep', pad)] * 4 if pad else [])
        try:
            _parallel.run_tasks(task, args, 2, None, ())
        except SlowToRebuild as err:
            print('RAISED', err, flush=True)
        faulthandler.cancel_dump_traceback_later()
        # Interpreter exit joins registered manager threads and queue feeders.
        faulthandler.dump_traceback_later(float(sys.argv[3]), exit=True)
''')


@pytest.mark.skipif(sys.platform == 'win32', reason='the partial message is written to a POSIX pipe')
@pytest.mark.parametrize('arguments', ['small_arguments', 'large_arguments'])
def test_failure_propagates_and_the_process_exits(tmp_path, arguments):
    script = tmp_path / 'child.py'
    script.write_text(CHILD, encoding='utf-8')
    watchdog = 30
    start = time.monotonic()
    completed = subprocess.run([sys.executable, '-I', str(script), str(PYTHON_ROOT), arguments, str(watchdog)],
                               capture_output=True, text=True, timeout=3 * watchdog)
    elapsed = time.monotonic() - start
    assert completed.returncode == 0, completed.stderr[-4000:]
    assert 'RAISED deliberate task failure' in completed.stdout, completed.stderr[-4000:]
    assert elapsed < watchdog
