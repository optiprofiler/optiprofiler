"""Ordered benchmark tasks with explicit worker-loss and abort cleanup."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from multiprocessing.reduction import ForkingPickler
from pathlib import Path
import pickle
import sys
import time


def worker_context():
    """Use fresh workers without changing the caller's global start method.

    Saved benchmarks already have logging threads when workers start. Forking
    then can inherit locked output buffers and hang during the child's final
    flush, even after its task has returned. The log queue uses this context too.
    """
    return mp.get_context('spawn')


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
        raise
    finally:
        executor.shutdown(wait=True)
    return results
