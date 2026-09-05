"""Saving another experiment must never overwrite an earlier invocation."""

from datetime import datetime, timezone
from concurrent.futures import ProcessPoolExecutor
import hashlib
import multiprocessing
from pathlib import Path
import shutil

import numpy as np
import pytest

from optiprofiler import Problem, benchmark
import optiprofiler.profiles as profiles
from optiprofiler.loader import load_results, save_results_to_h5


class FixedDateTime(datetime):
    @classmethod
    def now(cls, tz=None):
        value = cls(2026, 9, 5, 12, 34, 56, tzinfo=timezone.utc)
        return value if tz is None else value.astimezone(tz)


def identity_solver(fun, x0):
    fun(x0)
    return np.asarray(x0)


def zero_solver(fun, x0):
    x = np.zeros_like(x0)
    fun(x)
    return x


def _concurrent_benchmark(root):
    profiles.datetime = FixedDateTime
    benchmark([identity_solver, zero_solver],
              problem=Problem(lambda x: x[0] ** 2, [1.0]),
              savepath=root, benchmark_id='parallel',
              draw_hist_plots='none', silent=True)


def test_concurrent_benchmarks_reserve_distinct_outputs(tmp_path):
    # Separate benchmark processes also isolate logging resources as in real
    # independent callers; all four share the same filesystem and frozen clock.
    with ProcessPoolExecutor(max_workers=4, mp_context=multiprocessing.get_context('spawn')) as pool:
        list(pool.map(_concurrent_benchmark, [str(tmp_path)] * 4))
    assert len(list((tmp_path / 'parallel').iterdir())) == 4


@pytest.fixture
def saved_same_second_runs(tmp_path, monkeypatch):
    monkeypatch.setattr(profiles, 'datetime', FixedDateTime)
    monkeypatch.chdir(tmp_path)
    library = tmp_path / 'identitytoy'
    library.mkdir()
    (library / 'identitytoy_tools.py').write_text(
        'from optiprofiler import Problem\n'
        'def identitytoy_select(options):\n'
        '    return ["TOY"]\n'
        'def identitytoy_load(name):\n'
        '    return Problem(lambda x: x[0] ** 2, [1.0], name=name)\n',
        encoding='utf-8')
    scores = []
    for _ in range(2):
        scores.append(benchmark(
            [identity_solver, zero_solver], plibs=['identitytoy'],
            custom_problem_libs_path=library, benchmark_id='saved',
            draw_hist_plots='none', silent=True, max_tol_order=1,
            summarize_performance_profiles=False, summarize_data_profiles=False,
        )[0])
    markers = sorted((tmp_path / 'saved').glob('*/test_log/time_stamp_*.txt'))
    assert len(markers) == 2
    return tmp_path / 'saved', markers, scores


def test_save_load_roundtrip_preserves_source_and_accepts_collision_suffix(saved_same_second_runs):
    root, markers, scores = saved_same_second_runs
    originals = {path: hashlib.sha256(path.read_bytes()).hexdigest()
                 for path in root.glob('*/test_log/data_for_loading.h5')}
    for marker, expected in zip(markers, scores):
        loaded = benchmark(None, load=marker.stem[11:], benchmark_id='saved',
                           draw_hist_plots='none', silent=True, max_tol_order=1,
                           summarize_performance_profiles=False,
                           summarize_data_profiles=False)[0]
        np.testing.assert_equal(loaded, expected)
    assert len(list(root.glob('*/test_log/data_for_loading.h5'))) == 4
    for path, digest in originals.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest


def test_latest_ignores_incomplete_and_malformed_markers(saved_same_second_runs):
    root, markers, _ = saved_same_second_runs
    incomplete = root / 'incomplete' / 'test_log'
    incomplete.mkdir(parents=True)
    (incomplete / 'time_stamp_20990101_000000.txt').touch()
    malformed = root / 'malformed' / 'test_log'
    malformed.mkdir(parents=True)
    (malformed / 'time_stamp_not_an_id.txt').touch()
    shutil.copyfile(markers[0].parent / 'data_for_loading.h5', malformed / 'data_for_loading.h5')
    results, _ = load_results({}, {'load': 'latest', 'benchmark_id': str(root)})
    assert results[0]['problem_names'] == ['TOY']


def test_latest_orders_collision_suffixes_numerically(saved_same_second_runs):
    root, markers, _ = saved_same_second_runs
    results, _ = load_results({}, {'load': 'latest', 'benchmark_id': str(root)})
    base_id = markers[0].stem[11:26]
    for suffix in ('999', '1000'):
        directory = root / f'ordered_{base_id}_{suffix}' / 'test_log'
        directory.mkdir(parents=True)
        results[0]['feature_stamp'] = suffix
        save_results_to_h5(results, directory / 'data_for_loading.h5')
        (directory / f'time_stamp_{base_id}_{suffix}.txt').touch()
    latest, _ = load_results({}, {'load': 'latest', 'benchmark_id': str(root)})
    assert latest[0]['feature_stamp'] == '1000'


def test_failed_data_write_does_not_publish_load_marker(saved_same_second_runs, monkeypatch):
    import h5py
    root, markers, _ = saved_same_second_runs
    real_file = h5py.File

    def fail_output_write(name, mode='r', *args, **kwargs):
        if mode == 'w' and Path(name).name == 'data_for_loading.h5':
            raise OSError('simulated output filesystem failure')
        return real_file(name, mode, *args, **kwargs)

    monkeypatch.setattr(h5py, 'File', fail_output_write)
    benchmark(None, load='latest', benchmark_id='saved',
              draw_hist_plots='none', silent=True, max_tol_order=1,
              summarize_performance_profiles=False, summarize_data_profiles=False)
    assert len(list(root.glob('*/test_log/time_stamp_*.txt'))) == len(markers)
    assert len([path for path in root.iterdir() if path.is_dir()]) == 3


def test_explicit_ambiguous_legacy_timestamp_requires_narrower_path(saved_same_second_runs, monkeypatch):
    root, markers, _ = saved_same_second_runs
    duplicate = root / 'duplicate' / 'test_log'
    shutil.copytree(markers[0].parent, duplicate)
    with pytest.raises(ValueError, match='multiple saved experiments'):
        load_results({}, {'load': markers[0].stem[11:], 'benchmark_id': str(root)})
    monkeypatch.chdir(markers[0].parent.parent)
    scores = benchmark(None, load=markers[0].stem[11:], benchmark_id='.',
                       score_only=True, silent=True, max_tol_order=1)[0]
    assert scores.shape == (2,)


def test_same_second_benchmarks_preserve_previous_output(tmp_path, monkeypatch):
    monkeypatch.setattr(profiles, 'datetime', FixedDateTime)
    options = dict(problem=Problem(lambda x: x[0] ** 2, [1.0]),
                   savepath=str(tmp_path), benchmark_id='collision',
                   draw_hist_plots='none', silent=True)
    first_score = benchmark([identity_solver, identity_solver], **options)[0]
    root = tmp_path / 'collision'
    first_dir, = [path for path in root.iterdir() if path.is_dir()]
    sentinel = first_dir / 'keep.txt'
    sentinel.write_text('first experiment', encoding='utf-8')

    second_score = benchmark([identity_solver, identity_solver], **options)[0]

    assert sentinel.read_text(encoding='utf-8') == 'first experiment'
    assert len([path for path in root.iterdir() if path.is_dir()]) == 2
    np.testing.assert_equal(first_score, second_score)
