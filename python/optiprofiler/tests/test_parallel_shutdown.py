"""Shutting down a pool after a failure cannot hang on a stuck manager thread."""

import threading
import time
from concurrent.futures import process as futures_process

from optiprofiler import _parallel


class Wakeup:
    def __init__(self, event):
        self.event = event

    def wakeup(self):
        self.event.set()


class Python38Executor:
    """The executor attributes shutdown_after_failure relies on in Python 3.8."""

    def __init__(self, thread, wakeup):
        self._queue_management_thread = thread
        self._queue_management_thread_wakeup = wakeup
        self._shutdown_lock = threading.Lock()
        self._shutdown_thread = False
        self.shutdowns = []

    def shutdown(self, wait=True):
        self.shutdowns.append(wait)


def started(target):
    thread = threading.Thread(target=target, daemon=True)
    thread.start()
    return thread


def test_a_manager_thread_that_cannot_finish_is_abandoned_after_the_bound(monkeypatch):
    monkeypatch.setattr(_parallel, 'MANAGER_STOP_SECONDS', 0.2)
    release = threading.Event()
    thread = started(release.wait)
    futures_process._threads_wakeups[thread] = Wakeup(threading.Event())
    executor = Python38Executor(thread, Wakeup(threading.Event()))
    try:
        start = time.monotonic()
        _parallel.shutdown_after_failure(executor)
        assert time.monotonic() - start < 2
        assert executor._shutdown_thread
        assert thread not in futures_process._threads_wakeups
        assert executor.shutdowns == []
    finally:
        release.set()
        thread.join()


def test_a_manager_thread_that_finishes_is_shut_down_normally(monkeypatch):
    monkeypatch.setattr(_parallel, 'MANAGER_STOP_SECONDS', 5)
    woken = threading.Event()
    thread = started(woken.wait)
    executor = Python38Executor(thread, Wakeup(woken))
    _parallel.shutdown_after_failure(executor)
    assert not thread.is_alive()
    assert executor.shutdowns == [True]


def test_an_executor_without_that_thread_is_shut_down_normally():
    class Executor:
        def __init__(self):
            self.shutdowns = []

        def shutdown(self, wait=True):
            self.shutdowns.append(wait)

    executor = Executor()
    _parallel.shutdown_after_failure(executor)
    assert executor.shutdowns == [True]
