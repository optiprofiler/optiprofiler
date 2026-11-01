"""Public and provider regressions for runtime contract boundaries."""

import json
import logging
import pickle
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem, benchmark
from optiprofiler.legacy_compat import _recipe_from_archived
from optiprofiler.metadata import describe_callback
from optiprofiler.problem_libraries import load_problem_library, resolve_problem_library
from optiprofiler.profile_utils import check_validity_profile_options
import optiprofiler.profiles as profiles


def initial_solver(fun, x0):
    fun(x0)
    return x0


def provider(tmp_path, name='datatoy'):
    folder = tmp_path / name
    folder.mkdir()
    counter = tmp_path / 'imports.txt'
    (folder / (name + '_tools.py')).write_text(
        'from __future__ import annotations\n'
        'from dataclasses import dataclass\n'
        'from pathlib import Path\n'
        '@dataclass\nclass Settings:\n    size: int = 2\n'
        'Path({!r}).open("a").write("import\\n")\n'.format(str(counter)) +
        'def {}_select(options): return ["toy"]\n'.format(name) +
        'def {}_load(name): return None\n'.format(name)
    )
    return resolve_problem_library(name, folder), counter


def test_file_provider_has_module_identity_and_executes_once(tmp_path):
    ref, counter = provider(tmp_path)
    with ThreadPoolExecutor(max_workers=4) as pool:
        plugins = list(pool.map(lambda _: load_problem_library(ref), range(8)))
    assert all(p.select({}, {}) == ['toy'] for p in plugins)
    assert counter.read_text().splitlines() == ['import']


def test_same_provider_name_has_distinct_directory_identity(tmp_path):
    a = tmp_path / 'a'; a.mkdir()
    b = tmp_path / 'b'; b.mkdir()
    ra, ca = provider(a)
    rb, cb = provider(b)
    load_problem_library(ra)
    load_problem_library(rb)
    assert ca.read_text() == cb.read_text() == 'import\n'


def test_provider_failed_import_can_be_retried(tmp_path):
    ref, counter = provider(tmp_path)
    path = Path(ref.locator)
    original = path.read_text()
    path.write_text('raise RuntimeError("broken provider")\n')
    with pytest.raises(RuntimeError, match='broken provider'):
        load_problem_library(ref)
    path.write_text(original)
    load_problem_library(ref)
    assert counter.read_text() == 'import\n'


def test_custom_provider_ignores_bare_module_collision(monkeypatch):
    collision = ModuleType('custom1')
    collision.custom1 = lambda: dict(fun=lambda x: 0, x0=[99, 99])
    monkeypatch.setitem(sys.modules, 'custom1', collision)
    actual = load_problem_library(resolve_problem_library('custom')).load('custom1', {})
    assert actual.x0.tolist() != [99, 99]
    assert sys.modules['custom1'] is collision


@pytest.mark.parametrize('signal', [KeyboardInterrupt, SystemExit])
@pytest.mark.parametrize('report', [False, True])
def test_problem_loading_propagates_process_signals(tmp_path, monkeypatch, signal, report):
    def cancel(*args):
        raise signal('cancelled')
    monkeypatch.setattr(profiles, 'load_problem_library', lambda ref: SimpleNamespace(load=cancel))
    with pytest.raises(signal, match='cancelled'):
        profiles._solve_one_problem_wrapper([], None, None, 'toy', 3,
            {'silent': True, '_eval_report_enabled': report}, False, tmp_path,
            SimpleNamespace(name='canceltoy'), {})


@pytest.mark.parametrize('value', [float('nan'), float('inf'), -float('inf'), 0, -1])
def test_budget_rejects_nonfinite_and_nonpositive_values(value):
    with pytest.raises(ValueError, match='max_eval_factor'):
        check_validity_profile_options([], {'max_eval_factor': value})


def test_budget_rejects_boolean():
    with pytest.raises(TypeError, match='max_eval_factor'):
        check_validity_profile_options([], {'max_eval_factor': True})


def test_zero_tolerance_rejected_before_solver_calls():
    calls = []
    def solver(fun, x0, *args):
        calls.append(1)
        return x0
    with pytest.raises(ValueError, match='max_tol_order'):
        benchmark([solver, solver], plibs=['custom'], problem_names=['custom1'],
                  score_only=True, silent=True, n_jobs=1, max_tol_order=0)
    assert not calls


@pytest.mark.parametrize('entries', [[], [{}]])
def test_archive_diagnostics_only_claim_available_checks(entries):
    candidate = dict(feature_specification=[{'name': 'plain', 'options': {}}],
                     n_runs=1, seed=0, run_plain=False, problem_options={}, feature_stamp='plain')
    recipe = _recipe_from_archived(candidate, entries)
    assert recipe['replayable'] is True
    assert recipe['archived_experiment']['archive_comparison'] == 'not_performed'


def test_callback_classes_keep_distinct_descriptors():
    class CallbackA:
        pass
    class CallbackB:
        pass
    a, b = describe_callback(CallbackA), describe_callback(CallbackB)
    assert a != b
    assert a['module'] == __name__
    assert a['name'].endswith('CallbackA')


@pytest.mark.parametrize('name', ['plain', 'truncated+truncated'])
def test_bad_reference_shape_warns_when_recording_nan(name):
    p = Problem(lambda x: float(sum(x)), [0, 0],
                cub=lambda x: np.zeros(1 if np.all(x == 0) else 2))
    fp = FeaturedProblem(p, Feature(name), 10, 0)
    with pytest.warns(RuntimeWarning, match='Reference constraint violation.*ValueError'):
        assert fp.fun([1, 1]) == 2
    assert np.isnan(fp.maxcv_hist[-1])
    with pytest.raises(ValueError, match='size 1'):
        fp.maxcv([1, 1])


@pytest.mark.parametrize('name', ['plain', 'truncated+truncated'])
def test_detailed_reference_read_does_not_consume_observed_budget(name):
    fp = FeaturedProblem(Problem(lambda x: 0, [0, 0], cub=lambda x: [1]), Feature(name), 10, 0)
    assert fp._maxcv([1, 1])[0] == 1
    assert fp.n_eval_cub == fp._real_n_eval_cub == len(fp.cub_hist) == 0


@pytest.mark.parametrize('name', ['plain', 'truncated+truncated'])
@pytest.mark.parametrize('channel', ['cub', 'ceq'])
def test_no_history_skips_reference_reads(name, channel):
    calls = []
    def constraint(x):
        calls.append(1)
        return np.ones(1)
    fp = FeaturedProblem(Problem(lambda x: 0, [0, 0], **{channel: constraint}), Feature(name), 10, 0)
    calls.clear()
    getattr(fp, channel)([1, 1], record_hist=False)
    assert calls == [1]
    assert len(getattr(fp, channel + '_hist')) == 0


@pytest.mark.parametrize('kwargs', [{}, {'xl': [0, 1], 'xu': [2, 3]}, {'cub': lambda x: [1]}])
def test_quantized_maxcv_accepts_array_like_and_validates_dimension(kwargs):
    fp = FeaturedProblem(Problem(lambda x: 0, [0, 0], **kwargs), Feature('quantized'), 10, 0)
    assert fp.maxcv([1, 2]) == fp.maxcv(np.array([1, 2]))
    with pytest.raises(ValueError, match='size 2'):
        fp.maxcv([1])


@pytest.mark.parametrize('score_only', [False, True])
@pytest.mark.parametrize('fail', [False, True])
def test_public_saved_benchmark_restores_host_logging(tmp_path, monkeypatch, fail, score_only):
    root = logging.getLogger()
    original_handlers, original_level = list(root.handlers), root.level
    sentinel = logging.NullHandler()
    root.handlers[:] = [sentinel]
    root.setLevel(logging.ERROR)
    if fail:
        def abort(*args, **kwargs):
            raise RuntimeError('failed after logging setup')
        monkeypatch.setattr(profiles, '_solve_all_problems', abort)
    try:
        args = dict(plibs=['custom'], problem_names=['custom1'], n_jobs=1, silent=True, score_only=score_only,
                    savepath=str(tmp_path), draw_hist_plots='none', max_tol_order=1, max_eval_factor=2)
        if fail:
            with pytest.raises(RuntimeError, match='failed after logging setup'):
                benchmark([initial_solver, initial_solver], **args)
        else:
            benchmark([initial_solver, initial_solver], **args)
        assert root.handlers == [sentinel]
        assert root.level == logging.ERROR
    finally:
        root.handlers[:] = original_handlers
        root.setLevel(original_level)
        sentinel.close()


def _infinite_initial_merit(fun, cv, cv_init):
    fun = np.asarray(fun, dtype=float)
    return np.where(fun == 2.0, np.inf, fun)


@pytest.mark.parametrize('silent', [False, True])
def test_score_only_logs_through_the_benchmark_console(tmp_path, monkeypatch, capfd, silent):
    # A direct problem logs from this process; a library run from workers.
    monkeypatch.chdir(tmp_path)
    problem = Problem(lambda x: float(np.sum((x - 1.0) ** 2)), [0.0, 0.0], name='QUAD')
    benchmark([initial_solver, initial_solver], problem=problem, merit_fun=_infinite_initial_merit,
              n_runs=1, max_eval_factor=3, score_only=True, silent=silent)
    library = tmp_path / 'libraries' / 'halfload'
    library.mkdir(parents=True)
    (library / 'halfload_tools.py').write_text(
        'from optiprofiler import Problem\n'
        "def halfload_select(options): return ['GOOD', 'BROKEN']\n"
        'def halfload_load(name):\n'
        "    if name == 'BROKEN': raise ImportError('simulated missing runtime')\n"
        '    return Problem(lambda x: float(x @ x), [1.0, 2.0], name=name)\n', encoding='utf-8')
    benchmark([initial_solver, initial_solver], plibs=['halfload'], custom_problem_libs_path=str(library.parent),
              n_jobs=2, n_runs=1, max_eval_factor=2, max_tol_order=1, score_only=True, silent=silent,
              draw_hist_plots='none')
    out, err = capfd.readouterr()
    assert not list(tmp_path.rglob('log.txt'))
    if silent:
        assert 'cannot be scored' not in out + err and 'Failed to load' not in out + err
        assert 'Scores of the solvers' not in out
    else:
        assert "[WARNING] Problem 'QUAD': runs with an undefined or infinite initial merit" in out
        assert '[WARNING] Failed to load    BROKEN from "halfload"' in out
        assert '[INFO   ] Loading problem   GOOD   from "halfload".' in out
        assert '[INFO   ] Scores of the solvers:' in out
        assert 'cannot be scored' not in err and 'Failed to load' not in err


def test_reference_history_diagnostic_survives_warning_suppression(caplog):
    import warnings
    from optiprofiler.recording import reference_maxcv_or_nan

    def invalid(x):
        raise ValueError('wrong constraint size')

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert np.isnan(reference_maxcv_or_nan(invalid, [1.]))
    assert 'ValueError: wrong constraint size' in caplog.text


@pytest.mark.parametrize('with_init', [False, True])
def test_provider_relative_import_and_dataclass_have_package_context(tmp_path, with_init):
    ref, counter = provider(tmp_path)
    folder = Path(ref.locator).parent
    (folder / 'helper.py').write_text('NAMES = ["relative"]\n')
    if with_init:
        (folder / '__init__.py').write_text('from .helper import NAMES\n')
    tools = Path(ref.locator)
    tools.write_text(tools.read_text() + 'from .helper import NAMES\n'
                     'def datatoy_select(options): return NAMES\n')
    collision = ModuleType('datatoy')
    original = sys.modules.get('datatoy')
    sys.modules['datatoy'] = collision
    try:
        assert load_problem_library(ref).select({}, {}) == ['relative']
        assert load_problem_library(ref).select({}, {}) == ['relative']
        assert counter.read_text() == 'import\n'
        assert sys.modules['datatoy'] is collision
    finally:
        if original is None:
            sys.modules.pop('datatoy', None)
        else:
            sys.modules['datatoy'] = original


def test_failed_package_import_rolls_back_relative_modules(tmp_path):
    ref, counter = provider(tmp_path)
    tools = Path(ref.locator)
    (tools.parent / 'helper.py').write_text('NAMES = ["relative"]\n')
    original = tools.read_text()
    tools.write_text('from .helper import NAMES\nraise RuntimeError("retry")\n')
    before = set(sys.modules)
    with pytest.raises(RuntimeError, match='retry'):
        load_problem_library(ref)
    assert not any(n.startswith('_optiprofiler_provider_datatoy_') for n in set(sys.modules)-before)
    tools.write_text(original + 'from .helper import NAMES\n')
    assert load_problem_library(ref).select({}, {}) == ['toy']


@pytest.fixture
def reexporting_provider(tmp_path):
    ref, counter = provider(tmp_path, name='registry_toy')
    tools = Path(ref.locator)
    tools.write_text(
        'from pathlib import Path\n'
        'from optiprofiler import Problem\n'
        'with Path({!r}).open("a") as stream: stream.write("import\\n")\n'.format(str(counter)) +
        'class Payload: pass\n'
        'def registry_toy_select(options):\n'
        '    from . import NAMES_BY_TYPE\n'
        '    return NAMES_BY_TYPE[Payload]\n'
        'def registry_toy_load(name):\n'
        '    return Problem(fun=lambda x: sum(x**2), x0=[1., 2.], name=name)\n'
    )
    (tools.parent / '__init__.py').write_text(
        'from .registry_toy_tools import Payload\n'
        'NAMES_BY_TYPE = {Payload: ["toy"]}\n'
    )
    return ref, counter


def test_package_reexport_reuses_tools_module_identity(reexporting_provider):
    from optiprofiler.problem_libraries import _load_tools_module
    ref, counter = reexporting_provider
    tools = _load_tools_module(ref)
    package = sys.modules[tools.__package__]
    assert counter.read_text().splitlines() == ['import']
    assert package.registry_toy_tools is tools
    assert package.Payload is tools.Payload
    assert sys.modules[tools.__name__] is tools
    assert package.NAMES_BY_TYPE[tools.Payload] == ['toy']
    assert _load_tools_module(ref) is tools
    assert load_problem_library(ref).select({}, {}) == ['toy']
    assert counter.read_text().splitlines() == ['import']


def test_package_reexport_selects_in_public_benchmark(reexporting_provider):
    ref, counter = reexporting_provider
    scores, _, _ = benchmark(
        [initial_solver, initial_solver], plibs=['registry_toy'],
        custom_problem_libs_path=str(Path(ref.locator).parent),
        n_jobs=1, score_only=True, silent=True, run_plain=False,
        max_eval_factor=1, max_tol_order=1,
    )
    assert np.all(np.isfinite(scores))
    assert load_problem_library(ref).select({}, {}) == ['toy']
    assert counter.read_text().splitlines() == ['import']


@pytest.mark.parametrize('name', ['plain', 'truncated+truncated'])
def test_reference_history_diagnostic_is_bounded_per_trial(name, caplog):
    def malformed(x):
        return np.array([1.]) if np.all(x == 0) else np.array([1., 2.])
    root = Problem(lambda x: float(sum(x)), [0., 0.], cub=malformed)
    fp = FeaturedProblem(root, Feature(name), 3, 17)
    with pytest.warns(RuntimeWarning, match='Reference constraint') as recorded:
        for _ in range(3):
            assert fp.fun(np.array([1., 1.])) == 2.
    assert len(recorded) == 1
    assert len([r for r in caplog.records if r.name == 'optiprofiler.recording']) == 1
    assert np.all(np.isnan(fp.maxcv_hist))
    other = FeaturedProblem(root, Feature(name), 3, 17)
    with pytest.warns(RuntimeWarning, match='Reference constraint'):
        other.fun(np.array([1., 1.]))


def _probe_size_mismatch(fun, x0, *constraints):
    fun(np.ones(2))
    fun(np.ones(2))
    fun(np.full(2, 2.0))
    return x0


def _write_reftoy(tmp_path):
    library = tmp_path / 'reftoy'
    library.mkdir()
    (library / 'reftoy_tools.py').write_text(
        'import numpy as np\n'
        'from optiprofiler import Problem\n'
        'def malformed(x):\n'
        '    return [-1.0, -1.0] if np.all(np.asarray(x) == 1.0) else [-1.0]\n'
        'def reftoy_select(options):\n'
        "    return ['P']\n"
        'def reftoy_load(name):\n'
        '    return Problem(lambda x: float(np.sum((x - 1.0) ** 2)), [0.0, 0.0], cub=malformed, name=name)\n',
        encoding='utf-8')


@pytest.mark.parametrize('n_jobs', [1, 2])
@pytest.mark.parametrize('feature_name', ['plain', 'noisy+noisy'])
def test_reference_history_failure_reaches_the_report(tmp_path, monkeypatch, n_jobs, feature_name):
    monkeypatch.chdir(tmp_path)
    _write_reftoy(tmp_path)
    report_path = tmp_path / 'report.json'
    benchmark([_probe_size_mismatch, _probe_size_mismatch], plibs=['reftoy'], custom_problem_libs_path=str(tmp_path),
              feature_name=feature_name, n_jobs=n_jobs, n_runs=2, max_eval_factor=5, max_tol_order=1, run_plain=False,
              score_only=True, silent=True, report_path=str(report_path))
    report = json.loads(report_path.read_text(encoding='utf-8'))
    failures = [entry for entry in report['diagnostics'] if entry['code'] == 'reference_history_unavailable']
    # A deterministic feature executes one run per solver; later runs are copies.
    runs = (1,) if feature_name == 'plain' else (1, 2)
    assert {(entry['scope']['solver_index'], entry['scope']['run_index']) for entry in failures} == {
        (solver, run) for solver in (1, 2) for run in runs}
    assert all(entry['stage'] == 'numerical' and entry['scope']['exception_type'] == 'ValueError'
               and entry['scope']['library'] == 'reftoy' and entry['scope']['problem'] == 'P' for entry in failures)


def _probe_without_failure(fun, x0, *constraints):
    fun(np.full(2, 2.0))
    return x0


def _reference_failures(report_path):
    report = json.loads(report_path.read_text(encoding='utf-8'))
    return sorted((entry['scope']['role'], entry['scope']['solver_index'], entry['scope']['run_index'],
                   entry['scope']['exception_type'], entry['scope']['library'], entry['scope']['problem'])
                  for entry in report['diagnostics'] if entry['code'] == 'reference_history_unavailable')


@pytest.mark.parametrize('n_jobs', [1, 2])
def test_load_replays_archived_reference_history_failures(tmp_path, monkeypatch, n_jobs):
    import h5py
    monkeypatch.chdir(tmp_path)
    _write_reftoy(tmp_path)
    benchmark([_probe_size_mismatch, _probe_without_failure, _probe_size_mismatch], plibs=['reftoy'],
              custom_problem_libs_path=str(tmp_path), feature_name='noisy', n_jobs=n_jobs, n_runs=2,
              max_eval_factor=5, max_tol_order=1, run_plain=True, silent=True, draw_hist_plots='none',
              savepath=str(tmp_path), benchmark_id='saved', report_path=str(tmp_path / 'saved.json'))
    saved = _reference_failures(tmp_path / 'saved.json')
    # Stochastic primary runs are actual trials; the plain reference runs once.
    assert saved == sorted([('primary', solver, run, 'ValueError', 'reftoy', 'P')
                            for solver in (1, 3) for run in (1, 2)]
                           + [('plain_reference', solver, 1, 'ValueError', 'reftoy', 'P') for solver in (1, 3)])

    def load(name, **options):
        benchmark(None, load='latest', benchmark_id='saved', savepath=str(tmp_path / 'out'), score_only=True,
                  silent=True, draw_hist_plots='none', report_path=str(tmp_path / name), **options)
        return _reference_failures(tmp_path / name)

    assert load('loaded.json') == saved
    # Retained solvers are renumbered: archived solver 3 is now solver 1.
    assert load('subset.json', solvers_to_load=[2, 1]) == sorted(
        [('primary', 1, run, 'ValueError', 'reftoy', 'P') for run in (1, 2)]
        + [('plain_reference', 1, 1, 'ValueError', 'reftoy', 'P')])
    # An archive written before the field existed still loads; its failures are unknown.
    (archive,) = tmp_path.glob('saved/*/test_log/data_for_loading.h5')
    with h5py.File(archive, 'r+') as stream:
        for group in stream.values():
            del group['reference_history_failures']
            del group['results_plib_plain']['reference_history_failures']
    assert load('legacy.json') == []


@pytest.mark.parametrize('damage', ['not_json', 'wrong_keys', 'solver_out_of_range', 'not_one_entry_per_problem',
                                    'deep_nesting', 'ragged_pickled_list', 'zero_dimensional_pickled_array'])
def test_load_tolerates_a_damaged_reference_history_field(tmp_path, monkeypatch, damage):
    import h5py
    from optiprofiler.eval_report import schema_for_document, schema_resource
    from optiprofiler.tests.eval_report_contract import assert_valid
    monkeypatch.chdir(tmp_path)
    _write_reftoy(tmp_path)
    benchmark([_probe_size_mismatch, _probe_size_mismatch], plibs=['reftoy'], custom_problem_libs_path=str(tmp_path),
              n_jobs=1, n_runs=1, max_eval_factor=5, max_tol_order=1, run_plain=True, silent=True,
              draw_hist_plots='none', savepath=str(tmp_path), benchmark_id='saved')
    damaged = {'not_json': ['{broken'], 'wrong_keys': ['[{"solver": 1}]'],
               'solver_out_of_range': ['[{"solver_index": 9, "run_index": 1, "exception_type": "ValueError"}]'],
               'not_one_entry_per_problem': ['[]', '[]'], 'deep_nesting': ['[' * 200000],
               'ragged_pickled_list': [['a', 'b'], ['c']], 'zero_dimensional_pickled_array': np.array(5)}[damage]
    (archive,) = tmp_path.glob('saved/*/test_log/data_for_loading.h5')
    with h5py.File(archive, 'r+') as stream:
        group = stream['plib_0']
        del group['reference_history_failures']
        if 'pickled' in damage:
            group.create_dataset('reference_history_failures_pickled', data=np.void(pickle.dumps(damaged)))
        else:
            group.create_dataset('reference_history_failures', data=np.array(damaged, dtype=h5py.special_dtype(vlen=str)))
    report_path = tmp_path / 'loaded.json'
    scores = benchmark(None, load='latest', benchmark_id='saved', savepath=str(tmp_path / 'out'), score_only=True,
                       silent=True, draw_hist_plots='none', max_tol_order=1, solvers_to_load=[1, 0],
                       report_path=str(report_path))[0]
    assert np.shape(scores) == (2,)
    report = json.loads(report_path.read_text(encoding='utf-8'))
    assert_valid(report, schema_resource(*schema_for_document(report)))
    entries = [entry for entry in report['diagnostics']
               if entry['code'] in ('reference_history_unavailable', 'reference_metadata_unavailable')]
    # The damaged primary field is reported as unavailable, not replayed;
    # the intact plain-reference field is still replayed.
    assert [entry['scope'] for entry in entries if entry['scope']['role'] == 'primary'] == [
        {'library': 'reftoy', 'problem': 'P', 'role': 'primary', 'field': 'reference_history_failures'}]
    assert entries[0]['code'] == 'reference_metadata_unavailable'
    assert {(entry['scope']['solver_index'], entry['scope']['run_index']) for entry in entries
            if entry['scope']['role'] == 'plain_reference'} == {(1, 1), (2, 1)}
