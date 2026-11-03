"""Ordered benchmark tasks with explicit worker-loss and abort cleanup."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures import process as futures_process
import multiprocessing as mp
from multiprocessing.reduction import ForkingPickler
from pathlib import Path
import pickle
import sys
import threading
import time


def worker_context():
    """Use fresh workers without changing the caller's global start method.

    Saved benchmarks already have logging threads when workers start. Forking
    then can inherit locked output buffers and hang during the child's final
    flush, even after its task has returned. The log queue uses this context too.
    """
    return mp.get_context('spawn')


def reject_worker_bootstrap_call():
    """Refuse a benchmark started while a spawned worker imports ``__main__``."""
    # CPython sets this flag while a spawn child prepares and imports the main
    # module; it drives the interpreter's own refusal to start a process then.
    if getattr(mp.current_process(), '_inheriting', False):
        raise RuntimeError(
            'benchmark was called while a worker process was importing the main script. Parallel '
            'workers (n_jobs > 1) are started with the spawn method and import the script again, so '
            'an unprotected call would rerun the whole benchmark in every worker. Protect the call '
            "with `if __name__ == '__main__':`.")


def worker_main_available():
    """Whether fresh workers can start from this entry point (unlike stdin)."""
    main = sys.modules.get('__main__')
    spec = getattr(main, '__spec__', None)
    if spec is not None:
        return True
    filename = getattr(main, '__file__', None)
    return bool(filename and Path(filename).stem != 'ipython' and Path(filename).is_file())


class WorkerArgumentPickler(ForkingPickler):
    """Use executor reducers and reject globals that spawn cannot restore."""

    def persistent_id(self, obj):
        spec = getattr(sys.modules.get('__main__'), '__spec__', None)
        if (spec is not None and (spec.name == '__main__' or spec.name.endswith('.__main__'))
                and getattr(obj, '__module__', None) in {'__main__', '__mp_main__'}):
            # Spawn skips rebuilding package __main__.py. Imported module-level
            # callables still work (including under pytest); only globals from
            # the skipped entry point need the pre-execution serial fallback.
            raise pickle.PicklingError('A callable or type belongs to an unavailable package entry point.')
        return None


MANAGER_STOP_SECONDS = 5


def close_once(connection):
    """Let several threads close one pipe end: only the first closes it."""
    lock = threading.Lock()
    close = connection._close
    closed = []

    def close_first():
        with lock:
            if not closed:
                closed.append(True)
                close()
    connection._close = close_first


def shutdown_after_failure(executor):
    """Shut down an executor whose workers were already terminated."""
    # A worker killed while writing its result leaves part of a message in the
    # result pipe, and the manager thread then waits in recv() for the rest.
    # Once the workers are dead, closing the parent's write end turns that wait
    # into EOF (unless another process still holds a copy). Only this thread
    # closes that end, here and in shutdown().
    result_queue = getattr(executor, '_result_queue', None)
    if result_queue is not None:
        result_queue._writer.close()
    # Closing the parent's read end of the call queue makes a feeder thread
    # blocked on a full call pipe fail; Python 3.9 and 3.10 would otherwise
    # hang in the manager's broken-pool exit (gh-94777). CPython also closes
    # this end from the manager or the feeder thread, so run_tasks made the
    # close single-shot.
    call_queue = getattr(executor, '_call_queue', None)
    if call_queue is not None:
        call_queue._reader.close()
    # Python 3.8's manager is a daemon thread that interpreter exit joins: if
    # it still cannot finish, abandon it rather than leave benchmark hanging.
    # Python 3.9+ cannot abandon its non-daemon manager, so shutdown waits.
    thread = getattr(executor, '_queue_management_thread', None)
    if thread is not None:
        with executor._shutdown_lock:
            executor._shutdown_thread = True
        executor._queue_management_thread_wakeup.wakeup()
        thread.join(MANAGER_STOP_SECONDS)
        if thread.is_alive():
            futures_process._threads_wakeups.pop(thread, None)
            return
    executor.shutdown(wait=True)


def run_tasks(function, args, n_jobs, initializer, initargs):
    if not args:
        return []
    # Eager startup must not create idle interpreters for a small selection.
    executor = ProcessPoolExecutor(max_workers=min(n_jobs, len(args)), mp_context=worker_context(),
                                   initializer=initializer, initargs=initargs)
    # Eager startup lets the manager observe every worker sentinel before a
    # task can kill a newly spawned worker. Python 3.8 already starts eagerly.
    if hasattr(executor, '_safe_to_dynamically_spawn_children'):
        executor._safe_to_dynamically_spawn_children = False
    # Before any thread can close it; Connection.close is not safe against a
    # concurrent close, which could hit a reused descriptor.
    close_once(executor._call_queue._reader)
    futures = {}
    owned = set()
    results = [None] * len(args)
    try:
        for index, arg in enumerate(args):
            future = executor.submit(function, *arg)
            futures[future] = index
            # Python 3.8 has no public executor termination API. Retain only
            # this executor's process handles, even if its manager later exits.
            owned.update((executor._processes or {}).values())
        for future in as_completed(futures):
            results[futures[future]] = future.result()
    except BaseException:
        try:
            owned.update((executor._processes or {}).values())
            for future in futures:
                future.cancel()
            # shutdown(wait=True) alone waits for running tasks after cancellation.
            # Terminate our workers first; never enumerate caller-owned children.
            for process in owned:
                if process.is_alive():
                    process.terminate()
            deadline = time.monotonic() + 5
            for process in owned:
                process.join(max(0, deadline - time.monotonic()))
            for process in owned:
                if process.is_alive():
                    process.kill()
        finally:
            shutdown_after_failure(executor)
        raise
    executor.shutdown(wait=True)
    return results
