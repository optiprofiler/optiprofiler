"""The opt-in report is a public file contract, not another execution mode."""
import json
import hashlib
import os
import stat
from pathlib import Path

import numpy as np
import pytest

from optiprofiler import Problem, benchmark
from optiprofiler.plotting import prepare_history_plot_data
from optiprofiler.utils import ProfileOption
from optiprofiler.tests.eval_report_contract import assert_valid


def stay(fun, x0):
    fun(x0)
    return np.asarray(x0)


def zero(fun, x0):
    x = np.zeros_like(x0)
    fun(x)
    return x


def crash(fun, x0):
    fun(x0)
    raise RuntimeError('deliberate solver failure')


def long_walk(fun, x0):
    for value in [-1, -2, -3] + list(range(1, 129)):
        fun(np.full_like(x0, value))
    return np.asarray(x0)


def spike_walk(fun, x0):
    for value in range(1, 1001):
        fun(np.full_like(x0, value))
    return np.asarray(x0)


def bound_stay(fun, x0, *constraints):
    return stay(fun, x0)


def bound_zero(fun, x0, *constraints):
    return zero(fun, x0)


def worker_only(fun, x0):
    from multiprocessing import current_process
    assert current_process().name != 'MainProcess', 'solver did not execute on worker'
    fun(x0)
    return np.asarray(x0)


def _read(path):
    """Read a main report; every report read by a test must satisfy the schema.

    The schema is the shared Python/MATLAB contract, so a key spelled
    differently by this emitter fails here (MATLAB's fixture does the same).
    Reports are UTF-8 by contract: decode explicitly instead of trusting the
    platform default (cp1252 on Windows runners turned UTF-8 into mojibake).
    """
    report = json.loads(path.read_text(encoding='utf-8'),
                        parse_constant=lambda token: (_ for _ in ()).throw(ValueError(token)))
    if report.get('schema') == 'optiprofiler.eval_report/1':
        assert_valid(report, 'eval_report.schema.json')
    return report


def _plot_data(path, report):
    companion_path = path.parent / report['plot_data']['path']
    raw = companion_path.read_bytes()
    assert len(raw) == report['plot_data']['bytes']
    assert hashlib.sha256(raw).hexdigest() == report['plot_data']['sha256']
    companion = _read(companion_path)
    assert companion['evaluation_id'] == report['evaluation_id']
    assert companion['schema'] == 'optiprofiler.plot_data/1'
    assert_valid(companion, 'plot_data.schema.json')
    return companion


def test_saved_experiment_load_has_compact_main_and_complete_numeric_detail(tmp_path, monkeypatch):
    fixture = os.environ.get('OPTIPROFILER_REPORT_ARCHIVE')
    if not fixture:
        pytest.skip('set OPTIPROFILER_REPORT_ARCHIVE to an immutable saved experiment directory')
    fixture = Path(fixture).resolve()
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in fixture.rglob('*') if p.is_file()}
    monkeypatch.chdir(fixture)
    options = dict(load='latest', benchmark_id='.', savepath=str(tmp_path),
                   score_only=True, draw_hist_plots='none', silent=True,
                   ptype='u', mindim=1, maxdim=3, max_tol_order=3)
    expected = benchmark(None, **options)
    target = tmp_path / 'compact.json'
    actual = benchmark(None, report_path=target, **options)
    np.testing.assert_equal(actual[:2], expected[:2])
    report = _read(target)
    detail = _plot_data(target, report)
    assert report['artifacts'] == []
    assert 'previews' not in report['profiles']
    ids = {h['id'] for h in detail['histories']}
    assert ids
    for problem in report['problems']:
        for run in problem['runs']:
            assert run['history_ref'] in ids
            assert 'history_preview' not in run
            assert 'budget_note' not in run
            assert 'best_semantics' not in run['objective']
    plots = {plot['id']: plot for plot in detail['plots']}
    assert all(ref in plots for ref in report['profiles']['plot_refs'])
    assert any(plot['kind'] == 'history' for plot in plots.values())
    assert any(plot['kind'] == 'log_ratio' for plot in plots.values())
    assert len(detail['target_work']) == 3
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() == digest for p, digest in before.items())
    assert set(p.name for p in tmp_path.iterdir()) == {'compact.json', 'compact.plot_data.json'}


def _assert_artifacts_verified(report_path, report):
    """Every listed artifact is present with exact bytes and a matching SHA256."""
    assert report['artifacts']
    for artifact in report['artifacts']:
        assert artifact['status'] == 'present' and 'reason' not in artifact
        path = report_path.parent / report['artifact_root'] / artifact['path']
        assert path.stat().st_size == artifact['bytes']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == artifact['sha256']


def _library(root, name, problems):
    """An actual developer-path provider; discovery/loading/workers stay real."""
    directory = root / name
    directory.mkdir(parents=True)
    (directory / f'{name}_tools.py').write_text(
        'from optiprofiler import Problem\n'
        f'def {name}_select(options):\n    return {problems!r}\n'
        f'def {name}_load(name):\n'
        "    if name == 'BROKEN': raise ValueError('deliberate load failure')\n"
        "    return Problem(lambda x: float(x @ x), [1.0, 2.0], name=name)\n")


def _options(tmp_path, names):
    return dict(plibs=names, custom_problem_libs_path=tmp_path / 'libraries',
                ptype='u', mindim=1, maxdim=10, score_only=True,
                draw_hist_plots='none', silent=True, savepath=str(tmp_path),
                solver_names=['stay', 'zero'], n_runs=3, seed=17,
                n_jobs=1, max_eval_factor=5, max_tol_order=2)


def test_report_only_preserves_single_problem_scores_and_creates_no_plots(tmp_path):
    options = dict(problem=Problem(lambda x: float(x @ x), [1.0], name='QUAD'),
                   score_only=True, draw_hist_plots='none', silent=True,
                   savepath=str(tmp_path), n_runs=1)
    expected = benchmark([stay, zero], **options)
    target = tmp_path / 'eval_report.json'
    actual = benchmark([stay, zero], report_path=target, **options)
    np.testing.assert_array_equal(actual[0], expected[0])
    np.testing.assert_array_equal(actual[0], [0.0, 1.0])
    assert actual[1:] == expected[1:] == (None, None)
    report = _read(target)
    assert report['schema'] == 'optiprofiler.eval_report/1'
    assert report['status'] == 'completed'
    assert report['stages']['rendering']['status'] == 'not_requested'
    assert report['coverage']['completed'] == 1
    assert report['scores']['solver_scores'] == [0.0, 1.0]
    assert not list(tmp_path.rglob('*.pdf'))
    assert not list(tmp_path.rglob('*.h5'))
    detail = _plot_data(target, report)
    assert len(detail['histories']) == 2
    assert len(detail['plots']) == 2  # Raw/cummin objective panels, no constraints.


def test_empty_path_disables_report(tmp_path):
    result = benchmark([stay, zero], problem=Problem(lambda x: float(x @ x), [1.0]),
                       score_only=True, silent=True, report_path='', savepath=str(tmp_path))
    np.testing.assert_array_equal(result[0], [0, 1])
    assert not list(tmp_path.iterdir())


def test_whole_libraries_preserve_identity_partial_coverage_and_actual_runs(tmp_path):
    _library(tmp_path / 'libraries', 'evaltoy_a', ['SAME', 'BROKEN'])
    _library(tmp_path / 'libraries', 'evaltoy_b', ['SAME'])
    options = _options(tmp_path, ['evaltoy_a', 'evaltoy_b'])
    expected = benchmark([stay, zero], **options)
    target = tmp_path / 'whole.json'
    actual = benchmark([stay, zero], report_path=target, **options)
    np.testing.assert_array_equal(actual[0], expected[0])
    np.testing.assert_array_equal(actual[1], expected[1])
    report = _read(target)
    assert report['status'] == 'partial'
    assert report['coverage']['selected'] == 3
    assert report['coverage']['completed'] == 2
    assert report['coverage']['load_failed'] == 1
    solved = [p for p in report['problems'] if p['status'] == 'completed']
    detail = _plot_data(target, report)
    histories = {h['id']: h for h in detail['histories']}
    assert len({p['id'] for p in solved}) == 2
    for problem in solved:
        assert problem['provider']['source'] == 'custom'
        assert problem['library'] in ('evaltoy_a', 'evaltoy_b')
        assert problem['budget'] == {'evaluations': 10, 'rule': 'ceil(max_eval_factor*dimension)'}
        assert len(problem['runs']) == 6
        for run in problem['runs']:
            assert run['evaluations'] == 1
            assert 'budget' not in run and run['budget_reached'] is False
            # The integer 0 means every retained evaluation was finite.
            assert run['objective']['invalid_evaluations'] == 0
            assert 'first_invalid_evaluation_index' not in run['objective']
            assert 'availability_reason' not in run['objective']
            assert 'convergence' not in run
            assert histories[run['history_ref']]['channels']['objective']['count'] == 1
            assert run['execution']['kind'] == ('actual' if run['run_index'] == 1 else 'repeated')
            assert run['oracle_seed'] == 23333 * 17
            assert run.get('oracle_seed_reason') == (None if run['run_index'] == 1 else 'copied_from_source_run')
    assert report['scores']['profile_scores'] == actual[1].tolist()
    assert len(report['profiles']['convergence']) == 2
    assert report['profiles']['convergence'][0]['history'][1]['hits'] == 6
    assert report['stages']['scoring']['status'] == 'completed'
    assert any(d['code'] == 'problem_load_failed' for d in report['diagnostics'])


def test_empty_and_all_failed_selection_are_not_success(tmp_path):
    _library(tmp_path / 'libraries', 'evalempty', [])
    _library(tmp_path / 'libraries', 'evalfailed', ['BROKEN'])
    for library, status in [('evalempty', 'empty'), ('evalfailed', 'failed')]:
        target = tmp_path / (library + '.json')
        benchmark([stay, zero], report_path=target, **_options(tmp_path, [library]))
        report = _read(target)
        assert report['status'] == status
        assert report['coverage']['completed'] == 0
        # Same stage vocabulary as MATLAB for these two outcomes.
        expected = ('not_applicable', 'empty_selection') if status == 'empty' else ('failed', 'all_selected_problem_loads_failed')
        assert (report['stages']['numerical']['status'], report['stages']['numerical']['reason']) == expected
        assert report['stages']['scoring'] == {'status': 'not_applicable', 'reason': 'no_loaded_problems'}


def test_collision_refuses_before_calling_solver_and_preserves_bytes(tmp_path):
    target = tmp_path / 'existing.json'
    target.write_text('owned by another evaluation')
    before = target.read_bytes()
    with pytest.raises(FileExistsError):
        benchmark([stay, zero], report_path=target,
                  problem=Problem(lambda x: float(x @ x), [1.0]), score_only=True)
    assert target.read_bytes() == before


def test_companion_collision_refuses_before_execution_and_leaves_main_absent(tmp_path):
    target = tmp_path / 'evaluation.json'
    companion = tmp_path / 'evaluation.plot_data.json'
    companion.write_bytes(b'owned by another evaluation')
    with pytest.raises(FileExistsError):
        # Invalid solver/configuration would fail if preflight did not reject
        # the companion first; no numerical work is needed for this contract.
        benchmark(None, report_path=target, nonexistent_option=True)
    assert companion.read_bytes() == b'owned by another evaluation'
    assert not target.exists()


def test_fatal_configuration_failure_leaves_strict_failed_report(tmp_path):
    target = tmp_path / 'failed.json'
    with pytest.raises(ValueError, match='Unknown option'):
        benchmark([stay, zero], report_path=target, nonexistent_option=True)
    report = _read(target)
    assert report['status'] == 'failed'
    assert report['coverage']['completed'] == 0
    assert report['diagnostics'][-1]['code'] == 'benchmark_exception'


def test_solver_abnormality_is_observed_without_reclassifying_numerical_completion(tmp_path):
    _library(tmp_path / 'libraries', 'evalcrash', ['QUAD'])
    target = tmp_path / 'solver.json'
    benchmark([crash, zero], report_path=target, **_options(tmp_path, ['evalcrash']))
    report = _read(target)
    assert report['stages']['numerical']['status'] == 'completed'
    run = report['problems'][0]['runs'][0]
    assert run['abnormal_termination'] is True
    assert run['output_fallback'] is True
    assert run['evaluations'] == 1
    assert 'convergence' not in run  # never inferred; see semantics.convergence
    assert report['semantics']['convergence'].startswith('Never inferred')


def test_merit_failure_preserves_completed_raw_measurements(tmp_path):
    def bad_merit(f, cv, cv0):
        raise ArithmeticError('deliberate merit failure')
    target = tmp_path / 'merit.json'
    with pytest.raises(ArithmeticError, match='deliberate merit failure'):
        benchmark([stay, zero], problem=Problem(lambda x: float(x @ x), [1.0]),
                  merit_fun=bad_merit, score_only=True, silent=True, report_path=target)
    report = _read(target)
    assert report['stages']['numerical']['status'] == 'completed'
    assert report['stages']['scoring']['status'] == 'failed'
    assert report['coverage']['completed'] == 1
    assert report['problems'][0]['runs'][1]['objective']['output'] == 0


def test_reporting_does_not_call_user_callbacks_or_consume_numpy_rng(tmp_path):
    calls = []
    def merit(f, cv, cv0):
        calls.append((f, cv, cv0))
        return float(f) + float(cv)
    problem = Problem(lambda x: float(x @ x), [1.0], name='QUAD')
    options = dict(problem=problem, score_only=True, silent=True, merit_fun=merit,
                   solver_isrand=[True, True], n_runs=3, seed=42)
    np.random.seed(734)
    without = benchmark([stay, zero], **options)
    count_without = len(calls)
    state_without = np.random.get_state()
    calls.clear()
    np.random.seed(734)
    with_report = benchmark([stay, zero], report_path=tmp_path / 'rng.json', **options)
    assert len(calls) == count_without
    state_after = np.random.get_state()
    assert state_after[0] == state_without[0]
    np.testing.assert_array_equal(state_after[1], state_without[1])
    assert state_after[2:] == state_without[2:]
    np.testing.assert_array_equal(with_report[0], without[0])


def test_save_load_records_source_and_artifacts_without_reexecuting_provider(tmp_path, monkeypatch):
    _library(tmp_path / 'libraries', 'evalsaved', ['QUAD'])
    options = _options(tmp_path, ['evalsaved'])
    options.update(score_only=False, n_jobs=2, n_runs=1, max_tol_order=1,
                   benchmark_id='saved', draw_hist_plots='sequential')
    saved_path = tmp_path / 'saved.json'
    original = benchmark([stay, zero], report_path=saved_path, **options)
    saved_report = _read(saved_path)
    # Every supported platform has a secure hashing path (openat on POSIX,
    # identity-checked opens on Windows), so a saved run is fully verified.
    assert saved_report['status'] == 'completed'
    assert saved_report['stages']['persistence'] == {'status': 'completed'}
    artifacts = saved_report['artifacts']
    assert any(a['path'].endswith('data_for_loading.h5') for a in artifacts)
    assert any('summary_' in a['path'] and a['path'].endswith('.pdf') for a in artifacts)
    _assert_artifacts_verified(saved_path, saved_report)
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (tmp_path / 'saved').rglob('*') if p.is_file()}
    # Public load must work with no solver functions and no provider tree.
    (tmp_path / 'libraries').rename(tmp_path / 'removed_libraries')
    monkeypatch.chdir(tmp_path)
    target = tmp_path / 'loaded.json'
    loaded = benchmark(None, load='latest', benchmark_id='saved', savepath=str(tmp_path),
                       score_only=True, silent=True, max_tol_order=1, report_path=target)
    np.testing.assert_array_equal(loaded[0], original[0])
    np.testing.assert_array_equal(loaded[1], original[1])
    report = _read(target)
    assert report['operation'] == 'load'
    assert report['status'] == 'completed'
    assert report['coverage']['scope'] == 'retained_archive'
    assert report['coverage']['load_failed'] is None
    assert report['source']['sha256'] == next(v for p, v in before.items() if p.name == 'data_for_loading.h5')
    assert report['artifacts'] == []
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() == value for p, value in before.items())
    # Load cannot infer the original budget from today's options.
    assert report['problems'][0]['budget'] == {'evaluations': None, 'reason': 'original_execution_budget_not_retained'}
    assert report['problems'][0]['runs'][0]['budget_reached'] is None
    assert report['problems'][0]['runs'][0]['oracle_seed'] is None
    assert report['problems'][0]['runs'][0]['oracle_seed_reason'] == 'execution_metadata_not_retained'
    assert report['configuration']['scope'] == 'current_load_selection_reanalysis_and_rendering'
    assert report['configuration']['effective']['feature']['scope'] == 'current_load_context_not_original_execution_feature'

    def invalid_merit(f, cv, cv0):
        raise ArithmeticError('invalid reload merit')
    failed_path = tmp_path / 'loaded-failed-merit.json'
    with pytest.raises(ValueError, match='invalid reload merit'):
        benchmark(None, load='latest', benchmark_id='saved', savepath=str(tmp_path),
                  score_only=True, silent=True, merit_fun=invalid_merit, report_path=failed_path)
    failed = _read(failed_path)
    assert failed['stages']['numerical']['status'] == 'completed'
    assert failed['stages']['scoring']['status'] == 'failed'
    assert failed['coverage']['completed'] == 1
    assert failed['problems'][0]['runs'][0]['objective']['initial'] == 5.0
    merit = failed['problems'][0]['runs'][0]['merit']
    assert merit['best'] is None
    assert merit['availability_reason'] == 'history_or_evaluation_count_unavailable'
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() == value for p, value in before.items())


def test_export_failure_is_not_a_numerical_failure(tmp_path, monkeypatch):
    from matplotlib.figure import Figure
    _library(tmp_path / 'libraries', 'evalrender', ['QUAD'])
    def failed_export(self, *args, **kwargs):
        raise OSError('deliberate file export failure')
    monkeypatch.setattr(Figure, 'savefig', failed_export)  # Real renderer's I/O boundary only.
    options = _options(tmp_path, ['evalrender'])
    options.update(score_only=False, n_runs=1, max_tol_order=1)
    target = tmp_path / 'render.json'
    with pytest.raises(OSError, match='deliberate file export failure'):
        benchmark([stay, zero], report_path=target, **options)
    report = _read(target)
    assert report['stages']['numerical']['status'] == 'completed'
    # The archive was saved and hashed before the export failure.
    assert report['stages']['persistence'] == {'status': 'completed'}
    assert report['stages']['rendering']['status'] == 'failed'
    assert report['stages']['scoring']['status'] == 'unknown'
    assert any(d['code'] == 'profile_export_failed' for d in report['diagnostics'])
    assert any(a['path'].endswith('data_for_loading.h5') for a in report['artifacts'])
    _assert_artifacts_verified(target, report)


def _saved_run(tmp_path, name):
    _library(tmp_path / 'libraries', name, ['QUAD'])
    options = _options(tmp_path, [name])
    options.update(score_only=False, n_runs=1, max_tol_order=1, benchmark_id=name)
    target = tmp_path / f'{name}.json'
    scores = benchmark([stay, zero], report_path=target, **options)
    return target, _read(target), scores


def test_hashing_unavailable_keeps_numerical_data_and_reports_unverified_artifacts(tmp_path, monkeypatch):
    from optiprofiler import eval_report as module
    monkeypatch.setattr(module, '_hashing_capability', lambda: 'unavailable')
    target, report, scores = _saved_run(tmp_path, 'evalnohash')
    # Integrity verification is missing; the saved numbers are not.
    assert report['status'] == 'partial'
    assert report['stages']['numerical'] == {'status': 'completed'}
    assert report['stages']['scoring'] == {'status': 'completed'}
    assert report['stages']['persistence'] == {'status': 'partial', 'reason': 'artifact_hash_unavailable'}
    assert report['artifacts']
    for artifact in report['artifacts']:
        assert artifact['status'] == 'unverified'
        assert artifact['bytes'] is None and artifact['sha256'] is None
        assert artifact['reason'] == 'secure_artifact_hashing_unavailable_on_platform'
        assert (target.parent / report['artifact_root'] / artifact['path']).is_file()
    assert [d for d in report['diagnostics'] if d['code'] == 'artifact_hash_unavailable']
    archive = next(a['path'] for a in report['artifacts'] if a['path'].endswith('data_for_loading.h5'))
    assert (target.parent / report['artifact_root'] / archive).stat().st_size > 0
    monkeypatch.chdir(tmp_path)
    loaded = benchmark(None, load='latest', benchmark_id='evalnohash', savepath=str(tmp_path),
                       score_only=True, silent=True, max_tol_order=1)
    np.testing.assert_array_equal(loaded[0], scores[0])


def test_identity_checked_hashing_matches_openat_and_refuses_symlinks(tmp_path, monkeypatch):
    """The Windows hashing path, exercised on every platform.

    It must produce the same receipts as the POSIX openat path and must not
    hash through a symlink placed inside the benchmark-owned output tree.
    """
    from optiprofiler import eval_report as module
    if module._hashing_capability() == 'unavailable':
        pytest.skip('no secure hashing primitive on this platform')
    reference_target, reference, _ = _saved_run(tmp_path, 'evalopenat')
    reference_root = reference_target.parent / reference['artifact_root']
    monkeypatch.setattr(module, '_hashing_capability', lambda: 'identity')
    # Both methods must produce the same receipt for the same files.
    for artifact in reference['artifacts']:
        receipt = module._digest(reference_root / artifact['path'], directory=reference_root)
        assert receipt == {'bytes': artifact['bytes'], 'sha256': artifact['sha256']}
    target, report, _ = _saved_run(tmp_path / 'identity', 'evalidentity')
    assert report['status'] == 'completed'
    _assert_artifacts_verified(target, report)
    root = target.parent / report['artifact_root']
    outside = tmp_path / 'outside.txt'
    outside.write_bytes(b'not an artifact')
    try:
        os.symlink(outside, root / 'planted_link.txt')
    except (OSError, NotImplementedError):
        pytest.skip('symlink creation is not permitted for this user')
    collector = module.EvalReport(tmp_path / 'identity' / 'again.json', {})
    collector.configure({}, {'score_only': False}, object(), output_dir=root)
    collector._initial_artifacts = set()
    collector._harvest()
    names = {a['path'].split('/')[-1] for a in collector.document['artifacts']}
    assert 'planted_link.txt' not in names and 'data_for_loading.h5' in names


def test_schemas_are_package_resources_shared_by_installed_tests_and_docs():
    from optiprofiler import eval_report as module
    for name, identifier in (('eval_report', 'urn:optiprofiler:eval_report:1'), ('plot_data', 'urn:optiprofiler:plot_data:1')):
        schema = module.load_schema(name)
        assert schema['$id'] == identifier
        assert schema['$schema'] == 'https://json-schema.org/draft/2020-12/schema'
        # The text is the resource itself (no doc copy to drift from).
        assert json.loads(module.schema_text(name)) == schema
    with pytest.raises(ValueError, match='Unknown EvalReport schema'):
        module.schema_text('agent_report')


def test_lossy_bins_preserve_extrema_indices_and_all_nonfinite_counts(tmp_path):
    def objective(x):
        if x[0] == -1:
            return np.nan
        if x[0] == -2:
            return np.inf
        if x[0] == -3:
            return -np.inf
        return float(x @ x)
    target = tmp_path / 'preview.json'
    benchmark([long_walk, zero], problem=Problem(objective, [1.0]),
              score_only=True, silent=True, max_eval_factor=200, report_path=target)
    report = _read(target)
    run = report['problems'][0]['runs'][0]
    assert run['evaluations'] == 131
    assert run['objective']['invalid_evaluations'] == dict(
        nan=1, positive_infinity=1, negative_infinity=1, observed_evaluations=131)
    detail = _plot_data(target, report)
    history = next(h for h in detail['histories'] if h['id'] == run['history_ref'])
    channel = history['channels']['objective']
    assert len(channel['bins']) == 32
    assert channel['bins'][0]['start_index'] == 1
    assert channel['bins'][-1]['end_index'] == 131
    # Exact integer edges shared with MATLAB: bin k covers (k-1)*n//32+1..k*n//32.
    assert [(b['start_index'], b['end_index']) for b in channel['bins']] == \
        [(k * 131 // 32 + 1, (k + 1) * 131 // 32) for k in range(32)]
    assert channel['bins'][0]['first'] == {'value': None, 'reason': 'nan'}
    assert channel['bins'][0]['nonfinite'] == {'nan': 1, 'positive_infinity': 1, 'negative_infinity': 1}
    # Finite extrema exclude infinities; the scalar best does not and can be -inf.
    assert channel['bins'][0]['finite_min'] == {'value': 1.0, 'evaluation_index': 4}
    assert channel['bins'][-1]['finite_max'] == {'value': 16384.0, 'evaluation_index': 131}
    assert run['objective']['first_invalid_evaluation_index'] == 1
    assert run['objective']['best'] == {'value': None, 'reason': 'negative_infinity'}
    assert run['objective']['best_evaluation_index'] == 3


def test_compact_bins_keep_spike_missed_by_uniform_preview_and_full_plot_vertices(tmp_path):
    def objective(x):
        return 1e6 if x[0] == 17 else float(x[0])
    target = tmp_path / 'spike.json'
    # Budget 1100 > 1002 on purpose: the padded history (max_eval = 1100 for
    # a 1-D problem) makes the renderer block-aggregate the raw display
    # series, so an array position is not an evaluation number.
    benchmark([spike_walk, zero], problem=Problem(objective, [1.0]),
              max_eval_factor=1100, n_runs=1, score_only=True, silent=True,
              report_path=target)
    report = _read(target)
    detail = _plot_data(target, report)
    run = report['problems'][0]['runs'][0]
    assert report['problems'][0]['budget']['evaluations'] == 1100
    history = next(h for h in detail['histories'] if h['id'] == run['history_ref'])
    bins = history['channels']['objective']['bins']
    assert len(bins) == 32
    assert bins[0]['finite_max'] == {'value': 1e6, 'evaluation_index': 17}
    raw = next(p for p in detail['plots'] if p['kind'] == 'history' and p['mode'] == 'raw')
    assert raw['padded_length'] == 1100
    assert raw['aggregation_trigger'] == 'padded_history_length_above_1002'
    series = raw['series'][0]
    indices = series['evaluation_indices']
    assert indices[0] == 1 and indices[-1] == 1000
    assert 256 < len(indices) < 1000  # Aggregated display copy, not one point per evaluation.
    # Non-decreasing; the renderer appends the last actual evaluation
    # unconditionally, so it may repeat the final interior block's vertex.
    assert indices == sorted(indices) and len(set(indices)) >= len(indices) - 1
    assert 3 not in indices  # The first two-point 'min' block dropped evaluation 3.
    # Locate the spike through its evaluation index, never through position.
    position = indices.index(17)
    assert position != 16
    assert series['mean'][position] == 1e6
    assert len(series['x']) == len(series['mean']) == len(series['lower']) == len(series['upper']) == len(indices)
    assert series['x'][position] == pytest.approx(17 / 2)


def _display_history(values):
    """A (solver=1, run=1, evaluation) display copy with one spike at 17."""
    history = np.asarray(values, dtype=float)[np.newaxis, np.newaxis, :]
    options = {ProfileOption.ERRORBAR_TYPE: 'minmax', ProfileOption.HIST_AGGREGATION: 'min'}
    n_eval = np.array([[history.shape[2]]])
    indices, means, _, _ = prepare_history_plot_data(history, False, 0.0, n_eval, options)
    return list(indices[0]), list(means[0])


def test_history_block_aggregation_boundary_is_the_padded_length():
    # 1002 padded evaluations = 1000 interior blocks of one point: every
    # evaluation keeps its own vertex, so position == evaluation index - 1.
    values = [float(i) for i in range(1, 1003)]
    values[16] = 1e6
    indices, means = _display_history(values)
    assert indices == list(range(1, 1003))
    assert means[16] == 1e6
    # 1003 padded evaluations: one block now holds two points and the 'min'
    # mode drops the larger one, so a spike at 17 is found by index only.
    values = [float(i) for i in range(1, 1004)]
    values[16] = 1e6
    indices, means = _display_history(values)
    assert len(indices) == 1002 and indices[0] == 1 and indices[-1] == 1003
    assert means[indices.index(17)] == 1e6
    # Aggregation is keyed on the padded length, not on actual evaluations:
    # 200 actual evaluations inside a 1003-long padded history are still
    # subject to the block rule (here the 'min' block merges two points).
    n_eval = np.array([[200]])
    history = np.asarray(values, dtype=float)[np.newaxis, np.newaxis, :]
    options = {ProfileOption.ERRORBAR_TYPE: 'minmax', ProfileOption.HIST_AGGREGATION: 'min'}
    short, _, _, _ = prepare_history_plot_data(history, False, 0.0, n_eval, options)
    assert short[0][-1] == 200 and 3 not in short[0] and 2 in short[0]


def test_report_captures_rendering_merit_without_extra_stateful_callback_calls(tmp_path, monkeypatch):
    from matplotlib.figure import Figure
    calls = []
    def stateful_merit(f, cv, cv0):
        calls.append(1)
        return float(f) + float(cv) + len(calls) / 100
    problem = Problem(lambda x: float(x @ x), [1.0], xl=[-2.0], xu=[2.0], name='BOUND')
    options = dict(problem=problem, merit_fun=stateful_merit, n_runs=1,
                   silent=True, savepath=str(tmp_path), solver_names=['stay', 'zero'])
    benchmark([bound_stay, bound_zero], **options)
    expected_calls = len(calls)
    calls.clear()
    captured = {}
    original_savefig = Figure.savefig
    def capture(self, *args, **kwargs):
        for ax in self.axes:
            if 'merit function' in ax.get_ylabel().lower():
                mode = 'cummin' if 'cummin' in ax.get_ylabel().lower() else 'raw'
                captured[mode] = [line.get_ydata().tolist() for line in ax.lines]
        return original_savefig(self, *args, **kwargs)
    monkeypatch.setattr(Figure, 'savefig', capture)  # Export I/O observation only.
    target = tmp_path / 'stateful.json'
    benchmark([bound_stay, bound_zero], report_path=target, **options)
    assert len(calls) == expected_calls
    report = _read(target)
    detail = _plot_data(target, report)
    # The direct-problem run rendered its history PDFs: say so, as MATLAB does.
    assert report['stages']['rendering'] == {'status': 'completed'}
    assert report['problems'][0]['rendering'] == {'status': 'completed'}
    assert any(a['path'].endswith('.pdf') for a in report['artifacts'])
    for plot in detail['plots']:
        if plot.get('channel') == 'merit':
            assert plot['observation_scope'] == 'rendering_inputs'
            assert [s['mean'] for s in plot['series']] == captured[plot['mode']]


def test_exact_large_score_tensor_and_many_solver_axis(tmp_path):
    _library(tmp_path / 'libraries', 'evalscores', ['QUAD'])
    options = _options(tmp_path, ['evalscores'])
    options.update(solver_names=['stay', 'zero', 'zero2', 'zero3', 'zero4'],
                   n_runs=1, max_tol_order=16)
    target = tmp_path / 'scores.json'
    _, profile_scores, _ = benchmark([stay, zero, zero, zero, zero], report_path=target, **options)
    assert profile_scores.size > 256
    report = _read(target)
    np.testing.assert_array_equal(report['scores']['profile_scores'], profile_scores)
    assert report['scores']['axis_values']['profile_type'] == ['performance', 'data']


def test_secondary_full_disk_reporting_error_preserves_original_save_error(tmp_path, monkeypatch):
    import os
    _library(tmp_path / 'libraries', 'evaldisk', ['QUAD'])
    real_replace = os.replace
    archive_failed = False
    def full_disk(source, destination, *args, **kwargs):
        nonlocal archive_failed
        if str(destination).endswith('data_for_loading.h5'):
            archive_failed = True
            raise OSError('original archive full disk')
        if archive_failed:
            raise OSError('secondary report full disk')
        return real_replace(source, destination, *args, **kwargs)
    monkeypatch.setattr(os, 'replace', full_disk)  # Filesystem atomic-rename boundary.
    options = _options(tmp_path, ['evaldisk'])
    options.update(score_only=False, n_runs=1)
    with pytest.raises(RuntimeError, match='Failed to save the experiment') as failure:
        benchmark([stay, zero], report_path=tmp_path / 'disk.json', **options)
    assert str(failure.value.__cause__.__cause__) == 'original archive full disk'


def test_score_callback_failure_does_not_blame_completed_rendering(tmp_path):
    _library(tmp_path / 'libraries', 'evalscorefailure', ['QUAD'])
    options = _options(tmp_path, ['evalscorefailure'])
    options.update(score_only=False, n_runs=1, max_tol_order=1)
    def failed_score(values):
        raise ArithmeticError('deliberate score callback failure')
    target = tmp_path / 'score-failure.json'
    with pytest.raises(ArithmeticError, match='deliberate score callback failure'):
        benchmark([stay, zero], score_fun=failed_score, report_path=target, **options)
    report = _read(target)
    assert report['stages']['scoring']['status'] == 'failed'
    assert report['stages']['rendering']['status'] == 'completed'
    assert report['stages']['numerical']['status'] == 'completed'


def test_failed_requested_plain_reference_keeps_primary_counts_but_is_partial(tmp_path):
    _library(tmp_path / 'libraries', 'evalplainpartial', ['QUAD'])
    tools = tmp_path / 'libraries' / 'evalplainpartial' / 'evalplainpartial_tools.py'
    source = tools.read_text(encoding='utf-8').replace(
        "    if name == 'BROKEN':",
        "    from pathlib import Path\n"
        "    marker = Path(__file__).with_suffix('.loaded')\n"
        "    if marker.exists(): raise ValueError('plain reference load failed')\n"
        "    marker.touch()\n"
        "    if name == 'BROKEN':")
    tools.write_text(source, encoding='utf-8')
    options = _options(tmp_path, ['evalplainpartial'])
    options.update(feature_name='noisy', n_runs=1, run_plain=True)
    target = tmp_path / 'plain-partial.json'
    benchmark([stay, zero], report_path=target, **options)
    report = _read(target)
    assert report['coverage']['selected'] == report['coverage']['completed'] == 1
    assert report['coverage']['load_failed'] == 0
    assert report['status'] == 'partial'
    assert report['stages']['numerical']['reason'] == 'requested_plain_reference_incomplete'
    assert report['stages']['scoring']['status'] == 'completed'


def test_parallel_report_collects_actual_worker_results(tmp_path):
    _library(tmp_path / 'libraries', 'evalworkers', ['FIRST', 'SECOND'])
    options = _options(tmp_path, ['evalworkers'])
    options.update(n_jobs=2, n_runs=1)
    target = tmp_path / 'workers.json'
    benchmark([worker_only, zero], report_path=target, **options)
    report = _read(target)
    assert report['coverage']['completed'] == 2
    for problem in report['problems']:
        run = problem['runs'][0]
        assert run['abnormal_termination'] is False
        assert run['output_fallback'] is False
        assert run['evaluations'] == 1


def test_report_files_are_owner_only_on_posix_and_say_so(tmp_path):
    target = tmp_path / 'private.json'
    benchmark([stay, zero], problem=Problem(lambda x: float(x @ x), [1.0], name='QUAD'),
              score_only=True, silent=True, report_path=target)
    report = _read(target)
    policy = report['report_files']
    assert policy['permission_policy'] == 'owner_read_write_only_best_effort'
    if os.name != 'posix':
        pytest.skip('POSIX file modes are not enforced on this platform; the report says so')
    for path in (target, tmp_path / report['plot_data']['path']):
        assert stat.S_IMODE(path.stat().st_mode) == 0o600
    assert policy['permissions_applied'] is True


def test_single_problem_report_shares_the_cross_language_vocabulary(tmp_path):
    target = tmp_path / 'vocabulary.json'
    benchmark([stay, zero], problem=Problem(lambda x: float(x @ x), [1.0], name='QUAD'),
              score_only=True, silent=True, n_runs=2, report_path=target)
    report = _read(target)
    detail = _plot_data(target, report)
    assert report['producer']['language'] == 'python'
    problem = report['problems'][0]
    assert problem['library'] == 'user'  # MATLAB uses the same label for a direct problem.
    assert problem['id'] == '["user","QUAD","primary"]'
    assert set(report['semantics']) == {'index_base', 'metric_best', 'budget', 'convergence',
                                        'run_defaults', 'configuration', 'paths', 'privacy'}
    assert set(detail['semantics']) == {'index_base', 'history_bins', 'history_plots', 'profile_plots',
                                        'error_bands', 'target_work', 'privacy'}
    assert report['stages']['persistence'] == {'status': 'not_applicable', 'reason': 'single_problem_has_no_reload_archive'}
    assert report['stages']['rendering'] == {'status': 'not_requested', 'reason': 'score_only'}
    for plot in detail['plots']:
        assert plot['x_transform'] == 'evaluation_index/(dimension+1)'
        assert plot['display_limit'] == 1e100
        assert plot['n_runs'] == 2 and plot['std_ddof'] == 0
        assert plot['observation_scope'] == 'retained_scoring_observations'
    run = problem['runs'][0]
    assert set(run) == {'solver_index', 'run_index', 'evaluations', 'budget_reached', 'objective',
                        'constraint', 'merit', 'abnormal_termination', 'output_fallback', 'execution',
                        'history_ref', 'oracle_seed', 'elapsed_seconds'}
    assert run['merit']['availability_reason'] == 'output_or_initial_unavailable'


def test_profile_plots_share_the_cross_language_vocabulary(tmp_path):
    _library(tmp_path / 'libraries', 'evalvocab', ['ONE', 'TWO'])
    options = _options(tmp_path, ['evalvocab'])
    options.update(n_runs=2, max_tol_order=1)
    target = tmp_path / 'profiles.json'
    # Two identical solvers: every problem/run pair is a log-ratio tie.
    benchmark([zero, zero], report_path=target, **options)
    report = _read(target)
    detail = _plot_data(target, report)
    plots = {plot['kind']: plot for plot in detail['plots'] if plot['kind'] != 'history'}
    assert set(plots) == {'performance', 'data', 'log_ratio'}
    assert plots['performance']['x_transform'] == 'log2(work/best_work)'
    assert plots['data']['x_transform'] == 'log2(1+work/(dimension+1))'
    for kind in ('performance', 'data'):
        plot = plots[kind]
        assert plot['observation_scope'] == 'scoring_profile_work'
        assert plot['target_work_ref'] == 'target-work-1'
        assert plot['failure_placeholder'] == pytest.approx(1.1 * plot['ratio_max'])
        assert plot['n_runs'] == 2 and plot['std_ddof'] == 0
        assert all(series['geometry'] == 'step' and series['band_visible'] for series in plot['series'])
    log_ratio = plots['log_ratio']
    assert log_ratio['x_transform'] == 'sorted_bar_position'
    assert log_ratio['problem_mapping'] is None
    assert log_ratio['problem_mapping_reason'] == 'legacy_bar_sort_does_not_retain_identity'
    assert log_ratio['tie_pairs'] == 4 and log_ratio['both_failed_pairs'] == 0
    assert log_ratio['bar_groups'] == [] and not any(log_ratio['series'][0]['visible'])
    assert log_ratio['failure_placeholder'] == pytest.approx(1.1 * log_ratio['ratio_max'])
    assert set(report['profiles']['plot_refs']) == {plot['id'] for plot in detail['plots'] if plot['kind'] != 'history'}
    assert {plot['history_or_output'] for plot in detail['plots'] if plot['kind'] != 'history'} == {'history', 'output'}


def test_unicode_identity_is_preserved_and_sensitive_unknown_option_is_redacted(tmp_path):
    # English label with a non-ASCII mathematical symbol: the repository
    # content is English, but user metadata is arbitrary Unicode and must
    # survive the UTF-8 round trip byte for byte on every platform.
    name = 'quadratic_\u03c6'  # Greek small phi, escaped so the source file stays ASCII
    target = tmp_path / 'unicode.json'
    benchmark([stay, zero], problem=Problem(lambda x: float(x @ x), [1.0], name=name),
              score_only=True, silent=True, report_path=target)
    raw = target.read_bytes()
    assert name.encode('utf-8') in raw and b'\\u03c6' not in raw  # written as UTF-8, not escaped ASCII
    assert _read(target)['problems'][0]['name'] == name
    assert _read(target)['problems'][0]['id'] == f'["user","{name}","primary"]'
    secret_path = tmp_path / 'secret.json'
    with pytest.raises(ValueError, match='Unknown option'):
        benchmark([stay, zero], report_path=secret_path, auth_token='do-not-copy-this-secret')
    assert 'do-not-copy-this-secret' not in secret_path.read_text(encoding='utf-8')
    assert _read(secret_path)['configuration']['request']['auth_token']['reason'] == 'redacted_sensitive_option'
