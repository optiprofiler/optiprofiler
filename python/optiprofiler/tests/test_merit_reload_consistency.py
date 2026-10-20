"""Cached merits, reference targets and raw archives must stay in one merit space.

These public save/load tests pin three contracts:

* Reloading an archive whose problems carry explicit references never compares
  cached merits from one merit function with targets converted by another.
* Archives without explicit references keep the V1 behaviour: a load without
  ``merit_fun`` reuses the cached merits.
* A merit function that fails after solving does not discard the raw results,
  and a reference that the selected merit maps to a non-finite target is
  reported before any solver runs.
"""
import hashlib
import json

import numpy as np
import pytest

from optiprofiler import benchmark
from optiprofiler.eval_report import schema_resource
from optiprofiler.tests.eval_report_contract import assert_valid

SOLVER_CALLS = []


def finish(fun, x0, *args):
    SOLVER_CALLS.append('finish')
    x = np.full_like(x0, 0.5)
    fun(x)
    return x


def stay(fun, x0, *args):
    SOLVER_CALLS.append('stay')
    fun(x0)
    return x0


def affine_merit(fun, maxcv, initial_maxcv):
    # Increasing affine map of the default merit for feasible points.
    return 2.0 * fun + 100.0 + maxcv


def square_merit(fun, maxcv, initial_maxcv):
    # Not affine: rescoring with the default merit must change the profiles.
    return fun * fun + maxcv


def failing_merit(fun, maxcv, initial_maxcv):
    raise RuntimeError('deliberate merit failure')


def _library(tmp_path, reference):
    library = tmp_path / 'libraries' / 'merit_reload'
    library.mkdir(parents=True, exist_ok=True)
    (library / 'merit_reload_tools.py').write_text(
        "from optiprofiler import Problem\n"
        "def objective(x): return float(x @ x)\n"
        "def merit_reload_select(options): return ['P1', 'P2']\n"
        "def merit_reload_load(name):\n"
        f"    reference = {reference!r}\n"
        "    x0 = [2.0] if name == 'P1' else [3.0]\n"
        "    return Problem(objective, x0, name=name, reference=reference)\n",
        encoding='utf-8')
    return library.parent


def _run(libraries, benchmark_id, **options):
    defaults = dict(plibs=['merit_reload'], custom_problem_libs_path=str(libraries),
                    ptype='u', mindim=1, maxdim=1, n_jobs=1, max_eval_factor=2,
                    max_tol_order=3, draw_hist_plots='none', silent=True,
                    benchmark_id=benchmark_id, seed=7)
    defaults.update(options)
    return benchmark([finish, stay], **defaults)


def _load(benchmark_id, **options):
    return benchmark(None, load='latest', benchmark_id=benchmark_id, silent=True,
                     score_only=True, max_tol_order=3, **options)


def _archive_hash(root):
    archives = sorted(root.rglob('data_for_loading.h5'))
    assert len(archives) == 1
    return hashlib.sha256(archives[0].read_bytes()).hexdigest()


FEASIBLE = {'fun': 0.0, 'maxcv': 0.0, 'kind': 'feasible', 'source': 'analytic test'}


@pytest.mark.parametrize('merit_fun', [affine_merit, square_merit])
@pytest.mark.parametrize('feature_name, n_runs', [('plain', 1), ('perturbed_x0+noisy', 2)])
@pytest.mark.parametrize('run_plain', [False, True])
def test_reload_without_merit_fun_rescores_reference_archives_consistently(
        tmp_path, monkeypatch, merit_fun, feature_name, n_runs, run_plain):
    monkeypatch.chdir(tmp_path)
    libraries = _library(tmp_path, FEASIBLE)
    common = dict(feature_name=feature_name, n_runs=n_runs, run_plain=run_plain)
    custom = _run(libraries, 'custom', merit_fun=merit_fun, **common)
    default = _run(libraries, 'default', **common)
    source_hash = _archive_hash(tmp_path / 'custom')

    report_path = tmp_path / 'reload.json'
    reloaded = _load('custom', report_path=str(report_path))
    # The archive carries explicit references, so the omitted merit means the
    # default merit for every quantity: histories, initial values, the
    # optional plain baseline and the converted targets.
    np.testing.assert_array_equal(reloaded[0], default[0])
    np.testing.assert_array_equal(reloaded[1], default[1])
    report = json.loads(report_path.read_text(encoding='utf-8'))
    assert_valid(report, schema_resource('eval_report', 2))
    codes = [diagnostic['code'] for diagnostic in report.get('diagnostics', [])]
    assert 'cached_merits_recomputed' in codes

    # Supplying the original merit reproduces the original experiment.
    same = _load('custom', merit_fun=merit_fun)
    np.testing.assert_array_equal(same[0], custom[0])
    np.testing.assert_array_equal(same[1], custom[1])
    # Loading never changes the source archive.
    assert _archive_hash(tmp_path / 'custom') == source_hash


def test_reload_of_default_merit_reference_archive_is_unchanged_and_silent(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    libraries = _library(tmp_path, FEASIBLE)
    original = _run(libraries, 'default')
    report_path = tmp_path / 'reload.json'
    reloaded = _load('default', report_path=str(report_path))
    np.testing.assert_array_equal(reloaded[0], original[0])
    np.testing.assert_array_equal(reloaded[1], original[1])
    report = json.loads(report_path.read_text(encoding='utf-8'))
    codes = [diagnostic['code'] for diagnostic in report.get('diagnostics', [])]
    # Recomputing with the default merit reproduced the cached values exactly.
    assert 'cached_merits_recomputed' not in codes


def test_reload_without_references_keeps_cached_custom_merits(tmp_path, monkeypatch):
    # V1 behaviour: targets are empirical minima of the cached merits, so
    # reusing those merits on load is consistent and must stay unchanged.
    monkeypatch.chdir(tmp_path)
    libraries = _library(tmp_path, None)
    custom = _run(libraries, 'custom', merit_fun=square_merit, run_plain=True)
    reloaded = _load('custom')
    np.testing.assert_array_equal(reloaded[0], custom[0])
    np.testing.assert_array_equal(reloaded[1], custom[1])


def test_failing_merit_with_references_fails_before_any_solver_call(tmp_path, monkeypatch):
    # The reference preflight already calls the merit, so a broken merit is
    # reported before any expensive solve instead of after all of them.
    monkeypatch.chdir(tmp_path)
    libraries = _library(tmp_path, FEASIBLE)
    SOLVER_CALLS.clear()
    with pytest.raises(RuntimeError, match='deliberate merit failure'):
        _run(libraries, 'failed', merit_fun=failing_merit)
    assert SOLVER_CALLS == []


def test_failing_merit_keeps_a_reloadable_raw_archive(tmp_path, monkeypatch):
    # Without references nothing calls the merit before solving: it fails on
    # the completed results, which must be saved before the error propagates.
    monkeypatch.chdir(tmp_path)
    libraries = _library(tmp_path, None)
    SOLVER_CALLS.clear()
    with pytest.raises(RuntimeError, match='deliberate merit failure'):
        _run(libraries, 'failed', merit_fun=failing_merit, report_path=str(tmp_path / 'failed.json'))
    assert SOLVER_CALLS
    report = json.loads((tmp_path / 'failed.json').read_text(encoding='utf-8'))
    assert_valid(report, schema_resource('eval_report', 2))
    assert report['stages']['persistence']['status'] == 'completed'
    assert report['status'] != 'completed'
    # The completed solver results were persisted before the merit failure
    # propagated, and the archive is publishable as a loadable experiment.
    assert len(list((tmp_path / 'failed').rglob('time_stamp_*.txt'))) == 1
    source_hash = _archive_hash(tmp_path / 'failed')
    fresh = _run(libraries, 'fresh', merit_fun=affine_merit)
    rescued = _load('failed', merit_fun=affine_merit)
    np.testing.assert_array_equal(rescued[0], fresh[0])
    np.testing.assert_array_equal(rescued[1], fresh[1])
    default = _run(libraries, 'fresh_default')
    rescued_default = _load('failed')
    np.testing.assert_array_equal(rescued_default[0], default[0])
    np.testing.assert_array_equal(rescued_default[1], default[1])
    assert _archive_hash(tmp_path / 'failed') == source_hash


INCOMPATIBLE = {'fun': -3.0, 'maxcv': 0.5, 'kind': 'merit', 'source': 'best known infeasible point'}


def _bounded_library(tmp_path, reference, x0):
    library = tmp_path / 'libraries' / 'merit_preflight'
    library.mkdir(parents=True, exist_ok=True)
    (library / 'merit_preflight_tools.py').write_text(
        "from optiprofiler import Problem\n"
        "def objective(x): return float(x @ x)\n"
        "def merit_preflight_select(options): return ['B1']\n"
        "def merit_preflight_load(name):\n"
        f"    return Problem(objective, {x0!r}, name=name, xl=[0.0], xu=[1.0], reference={reference!r})\n",
        encoding='utf-8')
    return library.parent


def _run_bounded(libraries, benchmark_id, **options):
    defaults = dict(plibs=['merit_preflight'], custom_problem_libs_path=str(libraries),
                    ptype='b', mindim=1, maxdim=1, n_jobs=1, max_eval_factor=2,
                    max_tol_order=2, draw_hist_plots='none', silent=True,
                    benchmark_id=benchmark_id, seed=7)
    defaults.update(options)
    return benchmark([finish, stay], **defaults)


@pytest.mark.parametrize('jobs', [1, 2])
def test_reference_incompatible_with_merit_fails_before_any_solver_call(tmp_path, monkeypatch, jobs):
    monkeypatch.chdir(tmp_path)
    libraries = _bounded_library(tmp_path, INCOMPATIBLE, [0.5])
    SOLVER_CALLS.clear()
    with pytest.raises(ValueError, match='incompatible with the selected merit') as error:
        _run_bounded(libraries, 'preflight', n_jobs=jobs)
    # The message distinguishes a valid record from an invalid one and keeps
    # the established wording of the scoring-time check; a spawned worker
    # re-raises the same error in the parent.
    assert 'finite real scalar' in str(error.value)
    assert 'B1' in str(error.value)
    if jobs == 1:
        # Solver calls are only observable in this process.
        assert SOLVER_CALLS == []
    assert not list((tmp_path / 'preflight').rglob('data_for_loading.h5'))


def test_reference_compatible_with_custom_merit_is_scored(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    libraries = _bounded_library(tmp_path, INCOMPATIBLE, [0.5])
    scores = _run_bounded(libraries, 'compatible', merit_fun=affine_merit)[0]
    assert np.all(np.isfinite(scores))


def test_plain_baseline_never_triggers_the_reference_preflight(tmp_path, monkeypatch):
    # The primary runs start from a strongly infeasible perturbed point, so the
    # default merit is finite at the reference there. The plain baseline starts
    # from the feasible x0, where it is not, but plain records are not targets.
    monkeypatch.chdir(tmp_path)
    libraries = _bounded_library(tmp_path, INCOMPATIBLE, [0.5])
    scores = _run_bounded(libraries, 'plain_baseline', feature_name='perturbed_x0',
                          perturbation_level=50.0, n_runs=1, run_plain=True)[0]
    assert np.all(np.isfinite(scores))
