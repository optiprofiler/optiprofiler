"""Public provider identity and real spawned-worker configuration regressions."""

import importlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import signal
import subprocess
import sys

import pytest


def initial_solver(fun, x0):
    fun(x0)
    return x0


def _write_provider(root, tag):
    folder = root / 'identitytoy'
    folder.mkdir(parents=True)
    (folder / '__init__.py').write_text('TAG = ' + repr(tag) + '\n')
    (folder / 'settings.py').write_text(
        'from dataclasses import dataclass\n'
        '@dataclass\nclass Settings:\n'
        '    scale: float = 3.0\n    owner_pid: int = 0\n'
    )
    events = root / 'events'
    events.mkdir()
    (folder / 'identitytoy_tools.py').write_text(
        'import json, os\nfrom pathlib import Path\n'
        'from .settings import Settings\nfrom . import TAG\n'
        'from optiprofiler import Problem\n'
        'EVENTS = Path(' + repr(str(events)) + ')\n'
        'def record(event, **values):\n'
        # Windows does not guarantee atomic appends from separate processes.
        # One file per process preserves the evidence without weakening checks.
        '    with (EVENTS / (str(os.getpid()) + ".jsonl")).open("a") as stream:\n'
        '        stream.write(json.dumps(dict(event=event, pid=os.getpid(), '
        'tag=TAG, **values)) + "\\n")\n'
        'record("import")\n'
        'def identitytoy_get_default_options():\n'
        '    record("defaults")\n'
        '    return {"settings": Settings(owner_pid=os.getpid())}\n'
        'def identitytoy_validate_options(options):\n'
        '    record("validate")\n'
        '    assert type(options["settings"]) is Settings\n'
        '    return options\n'
        'def identitytoy_select(problem_options, library_options):\n'
        '    assert type(library_options["settings"]) is Settings\n'
        '    return ["one", "two"]\n'
        'def identitytoy_load(name, *, library_options):\n'
        '    settings = library_options["settings"]\n'
        '    assert type(settings) is Settings\n'
        '    record("load", name=name, owner_pid=settings.owner_pid, scale=settings.scale)\n'
        '    return Problem(lambda x: settings.scale * float(x @ x), '
        '[1., 2.], name=name)\n'
    )
    return folder, events


def _read_events(events):
    return [json.loads(line) for path in sorted(events.glob('*.jsonl'))
            for line in path.read_text().splitlines()]


def _worker_select(reference, options):
    from optiprofiler.problem_libraries import load_problem_library
    library = load_problem_library(reference)
    return library.select({}, options), options['settings'].owner_pid, os.getpid()


def _child(case, mode, root):
    from optiprofiler import benchmark, get_plib_config
    from optiprofiler.problem_libraries import resolve_problem_library

    folder, events = _write_provider(root / 'selected', 'selected')
    original_package = None
    if case == 'canonical':
        sys.path.insert(0, str(folder.parent))
    elif case == 'conflict':
        other, _ = _write_provider(root / 'other', 'other')
        sys.path.insert(0, str(other.parent))
        original_package = importlib.import_module('identitytoy')
        # Fresh workers see the selected package first, unlike this process's
        # already-loaded canonical package. The selected namespace must be pinned.
        sys.path.insert(0, str(folder.parent))

    options = get_plib_config('identitytoy', custom_problem_libs_path=folder)
    expected_scale = 3.0
    overrides = {}
    if case == 'canonical':
        settings = importlib.import_module('identitytoy.settings').Settings
        assert type(options['settings']) is settings
        overrides = {'identitytoy': {'settings': settings(4.0, os.getpid())}}
        expected_scale = 4.0
    else:
        assert type(options['settings']).__module__.startswith(
            '_optiprofiler_provider_identitytoy_')

    if mode in {'spawn', 'forkserver'}:
        from concurrent.futures import ProcessPoolExecutor
        from optiprofiler.problem_libraries import (
            _initialize_problem_library_worker, _problem_library_worker_package,
        )
        reference = resolve_problem_library('identitytoy', folder)
        with ProcessPoolExecutor(
            max_workers=1, mp_context=mp.get_context(mode),
            initializer=_initialize_problem_library_worker,
            initargs=(reference, None, (), _problem_library_worker_package(reference)),
        ) as executor:
            names, owner_pid, worker_pid = executor.submit(
                _worker_select, reference, options).result(timeout=15)
        assert names == ['one', 'two']
        assert owner_pid == os.getpid()
        assert worker_pid != os.getpid()
    else:
        import numpy as np
        jobs = int(mode)
        report_path = root / 'report.json'
        scores, _, _ = benchmark(
            [initial_solver, initial_solver], plibs=['identitytoy'],
            custom_problem_libs_path=folder, plib_options=overrides,
            n_jobs=jobs, score_only=True, silent=True, run_plain=False,
            max_eval_factor=1, max_tol_order=1,
            report_path=report_path,
        )
        assert np.all(np.isfinite(scores))
        # Independent runtime evidence: finite aggregate scores alone would not
        # detect a problem silently dropped after a provider load failure.
        report = json.loads(report_path.read_text())
        assert report['coverage']['selected'] == report['coverage']['completed'] == 2, report['coverage']
        assert report['coverage']['load_failed'] == 0, report['coverage']
        assert len(report['problems']) == 2
        for problem in report['problems']:
            assert len(problem['runs']) == 2
            for run in problem['runs']:
                assert run['evaluations'] == 1, run
                assert not run['abnormal_termination'], run
                assert not run['output_fallback'], run
        records = _read_events(events)
        loads = [record for record in records if record['event'] == 'load']
        assert sorted(record['name'] for record in loads) == ['one', 'two'], loads
        assert all(record['owner_pid'] == os.getpid() for record in loads)
        assert all(record['scale'] == expected_scale for record in loads)
        assert all((record['pid'] == os.getpid()) == (jobs == 1) for record in loads)

    records = _read_events(events)
    assert all(record['tag'] == 'selected' for record in records)
    assert all(record['pid'] == os.getpid() for record in records
               if record['event'] in {'defaults', 'validate'})
    imports = [record['pid'] for record in records if record['event'] == 'import']
    assert len(imports) == len(set(imports)), 'tools executed twice in one process'
    if original_package is not None:
        assert sys.modules['identitytoy'] is original_package
        assert original_package.TAG == 'other'
    print(json.dumps({'case': case, 'mode': mode, 'provider_identity': 'passed'}), flush=True)


def _run_child(case, mode, tmp_path):
    command = [sys.executable, str(Path(__file__).resolve()),
               '--provider-child', case, str(mode), str(tmp_path)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, start_new_session=(os.name != 'nt'))
    try:
        output, _ = process.communicate(timeout=40)
    except subprocess.TimeoutExpired:
        if os.name != 'nt':
            os.killpg(process.pid, signal.SIGKILL)
        else:
            subprocess.run(['taskkill', '/PID', str(process.pid), '/T', '/F'],
                           capture_output=True)
        output, _ = process.communicate()
        pytest.fail('provider worker watchdog expired:\n' + output)
    assert process.returncode == 0, output
    assert '"provider_identity": "passed"' in output


@pytest.mark.parametrize('case', ['canonical', 'synthetic', 'conflict'])
@pytest.mark.parametrize('jobs', [1, 2])
@pytest.mark.parametrize('_attempt', range(3) if os.name == 'nt' else [0])
def test_public_benchmark_preserves_provider_option_identity(case, jobs, _attempt, tmp_path):
    # Repeat on Windows, where shared-file receipts previously lost load events.
    _run_child(case, jobs, tmp_path)


@pytest.mark.parametrize('method', [m for m in ['spawn', 'forkserver']
                                  if m in mp.get_all_start_methods()])
@pytest.mark.parametrize('case', ['synthetic', 'conflict'])
def test_worker_bootstrap_pins_namespace_before_unpickling(case, method, tmp_path):
    _run_child(case, method, tmp_path)


def test_canonical_import_failure_rolls_back_new_provider_modules(tmp_path, monkeypatch):
    from optiprofiler.problem_libraries import load_problem_library, resolve_problem_library

    folder = tmp_path / 'rollback_identitytoy'
    folder.mkdir()
    (folder / '__init__.py').write_text('from . import helper\nraise RuntimeError("retry")\n')
    (folder / 'helper.py').write_text('VALUE = 1\n')
    (folder / 'rollback_identitytoy_tools.py').write_text(
        'def rollback_identitytoy_select(options): return []\n'
        'def rollback_identitytoy_load(name): return None\n')
    monkeypatch.syspath_prepend(str(tmp_path))
    reference = resolve_problem_library('rollback_identitytoy', folder)
    before = set(sys.modules)
    with pytest.raises(RuntimeError, match='retry'):
        load_problem_library(reference)
    assert not any(name.startswith(('rollback_identitytoy',
                                   '_optiprofiler_provider_rollback_identitytoy_'))
                   for name in set(sys.modules) - before)
    (folder / '__init__.py').write_text('from . import helper\n')
    assert load_problem_library(reference).select({}, {}) == []


def test_canonical_failure_preserves_preexisting_package_attributes(tmp_path, monkeypatch):
    from optiprofiler.problem_libraries import (
        _problem_library_worker_package, load_problem_library, resolve_problem_library,
    )

    name = 'rollback_existingtoy'
    folder = tmp_path / name
    folder.mkdir()
    (folder / '__init__.py').write_text('helper = "caller-owned"\n')
    (folder / 'helper.py').write_text('VALUE = 1\n')
    tools = folder / (name + '_tools.py')
    tools.write_text('from .helper import VALUE\nraise RuntimeError("retry")\n')
    monkeypatch.syspath_prepend(str(tmp_path))
    package = importlib.import_module(name)
    reference = resolve_problem_library(name, folder)
    with pytest.raises(RuntimeError, match='retry'):
        load_problem_library(reference)
    assert sys.modules[name] is package
    assert package.helper == 'caller-owned'
    assert name + '.helper' not in sys.modules
    assert name + '.' + name + '_tools' not in sys.modules
    tools.write_text(
        'from .helper import VALUE\n'
        'def rollback_existingtoy_select(options): return []\n'
        'def rollback_existingtoy_load(name): return None\n')
    assert load_problem_library(reference).select({}, {}) == []
    assert _problem_library_worker_package(reference) == name


if __name__ == '__main__':
    _child(sys.argv[2], sys.argv[3], Path(sys.argv[4]))
