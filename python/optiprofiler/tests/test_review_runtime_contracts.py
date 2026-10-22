"""Public and provider regressions for runtime contract boundaries."""

import logging
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


@pytest.mark.parametrize('fail', [False, True])
def test_public_saved_benchmark_restores_host_logging(tmp_path, monkeypatch, fail):
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
        args = dict(plibs=['custom'], problem_names=['custom1'], n_jobs=1, silent=True,
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
