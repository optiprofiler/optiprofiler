"""Exceptions after output setup must release the benchmark-owned logger."""

import logging
from pathlib import Path

import pytest

from optiprofiler import benchmark
import optiprofiler.profiles as profiles


def initial_solver(fun, x0):
    fun(x0)
    return x0


def test_public_exception_closes_owned_logging(tmp_path, monkeypatch):
    class Queue:
        closed = 0
        joined = 0

        def close(self):
            self.closed += 1

        def join_thread(self):
            self.joined += 1

    class Listener:
        handlers = []
        stopped = 0

        def stop(self):
            self.stopped += 1

    class Handler(logging.Handler):
        def __init__(self, queue):
            super().__init__()
            self.queue = queue

        def emit(self, record):
            pass

    queue, listener = Queue(), Listener()
    handler = Handler(queue)
    root_logger = logging.getLogger()
    root_logger.addHandler(handler)
    monkeypatch.setattr(profiles, 'setup_main_process_logging',
                        lambda **kwargs: (queue, listener))

    def fail_after_logging(*args, **kwargs):
        raise RuntimeError('post-setup failure')

    monkeypatch.setattr(profiles, '_solve_all_problems', fail_after_logging)
    try:
        with pytest.raises(RuntimeError, match='post-setup failure'):
            benchmark([initial_solver, initial_solver], plibs=['custom'],
                      problem_names=['custom1'], ptype='u', mindim=1, maxdim=2,
                      n_jobs=1, draw_hist_plots='none', savepath=str(tmp_path),
                      benchmark_id='logging-failure', silent=True)
        assert listener.stopped == 1
        assert queue.closed == 1
        assert queue.joined == 1
        assert handler not in root_logger.handlers
    finally:
        root_logger.removeHandler(handler)
        handler.close()


def test_public_benchmark_preserves_caller_script(tmp_path):
    from optiprofiler import Problem
    benchmark([initial_solver, initial_solver],
              problem=Problem(lambda x: float(x[0] ** 2), [1.0]),
              n_jobs=1, draw_hist_plots='none', savepath=str(tmp_path),
              benchmark_id='caller-script', silent=True)
    copies = list(tmp_path.glob('caller-script/*/test_log/' + Path(__file__).name))
    assert len(copies) == 1
    assert copies[0].read_bytes() == Path(__file__).read_bytes()
    assert not list(tmp_path.glob('caller-script/*/test_log/profiles.py'))
