"""Keep plain references attached to their saved problem-library container."""

from copy import deepcopy
import hashlib

import numpy as np
import pytest

from optiprofiler.loader import load_results, load_results_from_h5, save_results_to_h5
from optiprofiler.profile_utils import process_results, write_report
from optiprofiler.utils import ProfileOption


def _block(plib, names, values, n_runs=2, n_evals=3):
    values = np.asarray(values, dtype=float)
    histories = np.broadcast_to(values[:, None, None, None],
                                (len(names), 3, n_runs, n_evals)).copy()
    inits = np.broadcast_to(values[:, None], (len(names), n_runs)).copy()
    block = {
        'plib': plib, 'plib_options': {'variant': plib},
        'problem_names': names, 'problem_dims': np.full(len(names), 2),
        'problem_types': ['u'] * len(names), 'solver_names': ['a', 'b', 'c'],
        'fun_histories': histories, 'maxcv_histories': np.zeros_like(histories),
        'fun_outs': histories[..., -1], 'maxcv_outs': np.zeros_like(histories[..., -1]),
        'fun_inits': inits, 'maxcv_inits': np.zeros_like(inits),
        'merit_histories': histories, 'merit_outs': histories[..., -1],
        'merit_inits': inits, 'n_evals': np.full(histories.shape[:3], n_evals),
        'computation_times': np.ones(histories.shape[:3]),
        'solvers_successes': np.ones(histories.shape[:3], dtype=bool),
        'feature_stamp': 'noisy', 'ptype': 'u', 'mindim': 1, 'maxdim': 10,
    }
    for key in ['mbs', 'mlcons', 'mnlcons', 'mcons']:
        block['problem_' + key] = np.zeros(len(names), dtype=int)
    return block


def _paired(plib, value, plain_value):
    block = _block(plib, ['DUP'], [value])
    block['results_plib_plain'] = _block(plib, ['DUP'], [plain_value], n_runs=3, n_evals=2)
    return block


def _mins(blocks, enabled=True):
    return process_results(blocks, {ProfileOption.RUN_PLAIN: enabled})[3]


def test_same_name_in_different_libraries_is_not_the_same_problem():
    blocks = [_paired('library_a', 5, 1), _paired('library_b', 20, 15)]
    np.testing.assert_array_equal(_mins(blocks), [[1, 1], [15, 15]])


def test_same_provider_with_different_saved_options_stays_separate():
    blocks = [_paired('provider', 5, 1), _paired('provider', 20, 15)]
    blocks[0]['plib_options'] = {'n': 2}
    blocks[1]['plib_options'] = {'n': 10}
    np.testing.assert_array_equal(_mins(blocks), [[1, 1], [15, 15]])


def test_reordered_and_missing_plain_rows_match_inside_the_library():
    block = _block('library_a', ['A', 'B', 'C'], [10, 20, 30])
    block['results_plib_plain'] = _block('library_a', ['C', 'A'], [3, 1], n_evals=5)
    np.testing.assert_array_equal(_mins([block]), [[1, 1], [20, 20], [3, 3]])


def test_repeated_processing_does_not_pad_caller_result_arrays():
    first = _block('library_a', ['DUP'], [5], n_evals=2)
    first['results_plib_plain'] = _block('library_a', ['DUP'], [1])
    second = _block('library_b', ['DUP'], [20], n_evals=5)
    second['results_plib_plain'] = _block('library_b', ['DUP'], [15])
    blocks = [first, second]
    original = deepcopy(blocks)
    np.testing.assert_array_equal(_mins(blocks), [[1, 1], [15, 15]])
    np.testing.assert_array_equal(_mins(blocks), [[1, 1], [15, 15]])
    for actual, before in zip(blocks, original):
        np.testing.assert_array_equal(actual['merit_histories'], before['merit_histories'])
        np.testing.assert_array_equal(actual['fun_histories'], before['fun_histories'])


@pytest.mark.parametrize('empty_plain', [None, {}])
def test_absent_or_failed_plain_library_does_not_hide_later_libraries(empty_plain):
    first = _block('library_a', ['DUP'], [5])
    if empty_plain is not None:
        first['results_plib_plain'] = empty_plain
    np.testing.assert_array_equal(_mins([first, _paired('library_b', 20, 15)]),
                                  [[5, 5], [15, 15]])


def test_all_nan_plain_row_cannot_erase_a_valid_featured_reference():
    np.testing.assert_array_equal(_mins([_paired('library_a', 5, np.nan)]), [[5, 5]])


def test_plain_reference_uses_all_solvers_and_runs():
    block = _paired('library_a', 20, 15)
    block['results_plib_plain']['merit_histories'][0, 2, 2, 1] = 1
    np.testing.assert_array_equal(_mins([block]), [[1, 1]])


def test_disabled_plain_keeps_featured_reference():
    np.testing.assert_array_equal(_mins([_paired('library_a', 5, 1)], False), [[5, 5]])


def test_no_plain_option_or_data_keeps_featured_reference():
    result = process_results([_block('library_a', ['A'], [5])], {})
    np.testing.assert_array_equal(result[3], [[5, 5]])


def test_report_accepts_partial_plain_results(tmp_path):
    block = _block('library_a', ['A', 'B'], [10, 20])
    block['results_plib_plain'] = _block('library_a', ['B'], [15])
    report = tmp_path / 'report.txt'
    readme = tmp_path / 'README.txt'
    readme.touch()
    write_report({'score_only': False, 'run_plain': True, 'silent': False},
                 [block], report, readme)
    text = report.read_text()
    assert 'Number of problems selected: 2' in text
    assert 'Wall-clock time spent by all the solvers: 18.00 secs' in text


@pytest.mark.parametrize('where', ['featured', 'plain'])
def test_ambiguous_duplicate_names_within_a_library_are_rejected(where):
    block = _block('library_a', ['A'], [10])
    block['results_plib_plain'] = _block('library_a', ['A'], [1])
    duplicate = _block('library_a', ['A', 'A'], [1, 2])
    if where == 'plain':
        block['results_plib_plain'] = duplicate
    else:
        duplicate['results_plib_plain'] = block['results_plib_plain']
        block = duplicate
    with pytest.raises(ValueError, match='duplicate problem names'):
        _mins([block])


def test_existing_h5_layout_reload_solver_subset_and_problem_filter(tmp_path, monkeypatch):
    blocks = [_paired('library_a', 5, 1), _paired('library_b', 20, 15)]
    # A better plain value from an excluded solver must not enter the reference.
    blocks[1]['results_plib_plain']['merit_histories'][:, 0, :, :] = -100
    log_dir = tmp_path / 'bench' / 'saved' / 'log'
    log_dir.mkdir(parents=True)
    path = log_dir / 'data_for_loading.h5'
    save_results_to_h5(blocks, path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    (log_dir / 'time_stamp_20200101_000000.txt').touch()
    monkeypatch.chdir(tmp_path)
    options = {'load': '20200101_000000', 'benchmark_id': 'bench',
               'solvers_to_load': [2, 1], 'run_plain': True}
    loaded, loaded_options = load_results({'problem_names': ['DUP']}, options)
    assert loaded_options['solver_names'] == ['c', 'b']
    np.testing.assert_array_equal(_mins(loaded), [[1, 1], [15, 15]])
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    # Old data with per-problem rather than per-run initial values still loads.
    legacy = deepcopy(blocks)
    for block in legacy:
        block['merit_inits'] = block['merit_inits'][:, 0]
        block['results_plib_plain']['merit_inits'] = block['results_plib_plain']['merit_inits'][:, 0]
    legacy_path = tmp_path / 'legacy.h5'
    save_results_to_h5(legacy, legacy_path)
    np.testing.assert_array_equal(_mins(load_results_from_h5(legacy_path)), [[1, 1], [-100, -100]])


def test_legacy_cached_merits_cannot_hide_raw_nan_after_h5_reload(tmp_path):
    block = _paired('library_a', 5, 1)
    block['fun_histories'] = block['fun_histories'].copy()
    block['fun_histories'][..., 0] = np.nan
    block['merit_histories'][..., 0] = -100
    block['fun_outs'] = np.full_like(block['fun_outs'], np.nan)
    block['merit_outs'] = np.zeros_like(block['merit_outs'])
    # Legacy raw initial values can be 1-D even when cached merits are 2-D.
    block['maxcv_inits'] = np.array([np.nan])
    block['merit_inits'] = np.zeros_like(block['merit_inits'])
    plain = block['results_plib_plain']
    plain['maxcv_histories'][..., 0] = np.nan
    plain['merit_histories'][..., 0] = -200
    path = tmp_path / 'cached-merits.h5'
    save_results_to_h5([block], path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    loaded = load_results_from_h5(path)
    processed = process_results(loaded, {'run_plain': True})
    np.testing.assert_array_equal(processed[3], [[1, 1]])
    assert np.all(np.isposinf(processed[1]))
    assert np.all(np.isposinf(processed[2]))
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest


def test_merit_only_helper_data_remains_supported():
    block = _paired('library_a', 5, 1)
    for record in [block, block['results_plib_plain']]:
        for prefix in ['fun_', 'maxcv_']:
            for suffix in ['histories', 'outs', 'inits']:
                record.pop(prefix + suffix)
    np.testing.assert_array_equal(_mins([block]), [[1, 1]])


def test_legacy_single_init_merit_expands_without_losing_per_run_validity():
    block = _block('library_a', ['A', 'B'], [5, 15])
    block['merit_inits'] = np.array([5.0, 15.0])
    block['maxcv_inits'][1, 1] = np.nan
    processed = process_results([block], {'run_plain': False})
    np.testing.assert_array_equal(processed[2], [[5, 5], [15, np.inf]])


def _initial_solver(fun, x0):
    fun(x0)
    return x0


def _zero_solver(fun, x0):
    fun(x0)
    x = np.zeros_like(x0)
    fun(x)
    return x


def test_public_benchmark_then_reload_keeps_library_references(tmp_path, monkeypatch):
    from optiprofiler import benchmark
    import optiprofiler.profiles as profiles

    monkeypatch.chdir(tmp_path)
    library_root = tmp_path / 'libraries'
    for name, lower_value in [('library_a', 1.0), ('library_b', 15.0)]:
        root = library_root / name
        root.mkdir(parents=True)
        (root / (name + '_tools.py')).write_text(
            'import numpy as np\n'
            'from optiprofiler.opclasses import Problem\n'
            f'def {name}_select(options):\n'
            "    return ['DUP']\n"
            f'def {name}_load(name):\n'
            f'    return Problem(lambda x: {lower_value} + float(x @ x), np.ones(2), name=name)\n',
            encoding='utf-8',
        )
    observed = []
    original_process = profiles.process_results

    def record_reference(results, options):
        processed = original_process(results, options)
        observed.append(processed[3].copy())
        return processed

    monkeypatch.setattr(profiles, 'process_results', record_reference)
    common = dict(savepath=str(tmp_path), benchmark_id='reference-public',
                  run_plain=True, draw_hist_plots='none', silent=True,
                  n_jobs=1, max_tol_order=1)
    benchmark([_initial_solver, _zero_solver],
              plibs=['library_a', 'library_b'],
              custom_problem_libs_path=str(library_root),
              ptype='u', mindim=2, maxdim=2,
              feature_name='perturbed_x0', n_runs=2, max_eval_factor=2,
              solver_names=['initial', 'zero'], **common)
    np.testing.assert_array_equal(observed[-1], [[1, 1], [15, 15]])
    source = next((tmp_path / 'reference-public').rglob('data_for_loading.h5'))
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    benchmark(None, load='latest', solvers_to_load=[1, 0], score_only=True, **common)
    np.testing.assert_array_equal(observed[-1], [[1, 1], [15, 15]])
    assert hashlib.sha256(source.read_bytes()).hexdigest() == digest
