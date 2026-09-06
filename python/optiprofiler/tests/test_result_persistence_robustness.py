"""Failure must not replace a usable archive or advertise an unsaved run."""

import hashlib
import os
import stat
from io import StringIO
from pathlib import Path

import h5py
import numpy as np
import pytest

import optiprofiler.loader as loader
import optiprofiler.profiles as profiles
from optiprofiler import Feature, benchmark
from optiprofiler.profile_utils import process_results, _write_history_display_note
from .test_plain_reference_identity import _block
from .test_experiment_identity import saved_same_second_runs


@pytest.mark.parametrize('raw_fields', [True, False])
@pytest.mark.parametrize('roundtrip', [True, False])
def test_all_nan_featured_merits_keep_valid_initial_reference(tmp_path, raw_fields, roundtrip):
    # Valid raw evaluations may produce NaN in a user merit callback. The
    # initial point remains an available reference, as in MATLAB omitnan.
    record = _block('toy', ['A'], [3.0])
    record['merit_histories'] = np.full_like(record['merit_histories'], np.nan)
    if not raw_fields:
        for key in ('fun_histories', 'maxcv_histories', 'fun_outs', 'maxcv_outs',
                    'fun_inits', 'maxcv_inits'):
            record.pop(key)
    records = [record]
    if roundtrip:
        path = tmp_path / 'data.h5'
        loader.save_results_to_h5(records, path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        records = loader.load_results_from_h5(path)
    result = process_results(records, {})
    np.testing.assert_array_equal(result[3], [[3.0, 3.0]])
    assert np.isnan(records[0]['merit_histories']).all()
    if roundtrip:
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest


def test_all_nan_history_and_initial_stay_undefined():
    record = _block('toy', ['A'], [np.nan])
    # A legacy record without raw validity information cannot manufacture a
    # finite reference if neither the history nor the initial merit has one.
    for key in ('fun_histories', 'maxcv_histories', 'fun_outs', 'maxcv_outs',
                'fun_inits', 'maxcv_inits'):
        record.pop(key)
    assert np.isnan(process_results([record], {})[3]).all()


def test_nan_raw_initial_is_not_rescued_by_custom_merit():
    record = _block('toy', ['A'], [3.0])
    record['fun_inits'] = np.full_like(record['fun_inits'], np.nan)
    record['merit_histories'] = np.full_like(record['merit_histories'], np.nan)
    # A fabricated finite saved initial merit is masked from the raw fields;
    # the finite-init fallback must not undo the NaN constraint/objective repair.
    result = process_results([record], {})
    assert np.isposinf(result[2]).all()
    assert np.isposinf(result[3]).all()


def test_display_report_is_explicit_without_mutating_results():
    record = _block('toy', ['A'], [3.0])
    record['fun_histories'][0, 0, 0, :] = [1e308, np.inf, 2.0]
    before = record['fun_histories'].copy()
    stream = StringIO()
    _write_history_display_note(stream, [record])
    assert '1e+100' in stream.getvalue()
    assert 'fun_histories: 1 clipped, 1 nonfinite' in stream.getvalue()
    assert 'Raw data, oracle values and scores are not clipped' in stream.getvalue()
    np.testing.assert_array_equal(record['fun_histories'], before)


@pytest.mark.parametrize('truth', [True, False])
def test_new_quantized_truth_note_is_specific(tmp_path, truth):
    paths = [tmp_path / 'report.txt', tmp_path / 'README.txt']
    profiles._append_quantized_truth_note(
        Feature('quantized', ground_truth=truth), False, [], *paths)
    for path in paths:
        text = path.read_text()
        assert f'ground_truth={str(truth).lower()}' in text
        assert ('featured' if truth else 'original') in text
        assert 'returned point is unchanged' in text


@pytest.mark.parametrize('stamp', ['quantized_1_ground_truth', 'custom_label', None])
def test_reload_does_not_claim_legacy_quantized_data_was_repaired(tmp_path, stamp):
    paths = [tmp_path / 'report.txt', tmp_path / 'README.txt']
    # The caller's current feature options cannot certify old saved channels.
    profiles._append_quantized_truth_note(
        Feature('plain'), True, [{'feature_stamp': stamp}], *paths)
    for path in paths:
        text = path.read_text()
        assert 'older quantized archives may contain mixed truth channels' in text
        assert 'Replotting does not repair' in text
        assert 'Loaded quantized experiment' not in text


def test_failed_h5_write_preserves_existing_archive(tmp_path, monkeypatch):
    path = tmp_path / 'saved.h5'
    original = [_block('toy', ['A'], [3.0])]
    loader.save_results_to_h5(original, path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    original_create = h5py.Group.create_dataset
    calls = 0

    def fail_after_one_dataset(self, *args, **kwargs):
        nonlocal calls
        calls += 1
        if calls > 1:
            raise OSError('injected storage failure after a partial write')
        return original_create(self, *args, **kwargs)

    monkeypatch.setattr(h5py.Group, 'create_dataset', fail_after_one_dataset)
    with pytest.raises(OSError, match='injected storage failure'):
        loader.save_results_to_h5([_block('new', ['B'], [4.0])], path)
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    assert list(tmp_path.iterdir()) == [path]


@pytest.mark.skipif(os.name == 'nt', reason='POSIX permission modes')
def test_atomic_archive_preserves_normal_creation_and_existing_modes(tmp_path):
    direct = tmp_path / 'direct.h5'
    loader._write_results_to_h5([_block('toy', ['A'], [3.0])], direct)
    path = tmp_path / 'atomic.h5'
    loader.save_results_to_h5([_block('toy', ['A'], [3.0])], path)
    assert stat.S_IMODE(path.stat().st_mode) == stat.S_IMODE(direct.stat().st_mode)
    path.chmod(0o640)
    loader.save_results_to_h5([_block('toy', ['A'], [4.0])], path)
    assert stat.S_IMODE(path.stat().st_mode) == 0o640


@pytest.mark.skipif(os.name == 'nt', reason='POSIX permission modes')
def test_private_archive_staging_does_not_expose_data(tmp_path, monkeypatch):
    path = tmp_path / 'private.h5'
    loader.save_results_to_h5([_block('toy', ['A'], [3.0])], path)
    path.chmod(0o600)
    original_write = loader._write_results_to_h5
    modes = []

    def inspect_write(results, staging):
        modes.append(stat.S_IMODE(Path(staging).stat().st_mode))
        original_write(results, staging)
        modes.append(stat.S_IMODE(Path(staging).stat().st_mode))

    monkeypatch.setattr(loader, '_write_results_to_h5', inspect_write)
    loader.save_results_to_h5([_block('toy', ['A'], [4.0])], path)
    # Both the empty stage and the completed payload must remain private,
    # not only the finally published archive.
    assert len(modes) == 2 and all(mode & 0o077 == 0 for mode in modes)
    assert stat.S_IMODE(path.stat().st_mode) == 0o600


def test_fsync_handle_is_write_capable(tmp_path, monkeypatch):
    original_open = Path.open

    def require_write_handle(self, mode='r', *args, **kwargs):
        if self.suffix == '.tmp' and 'b' in mode:
            # Windows FlushFileBuffers requires GENERIC_WRITE. Model that
            # boundary explicitly even on POSIX, where read-only fsync works.
            assert '+' in mode or 'w' in mode, 'fsync requires a write-capable handle'
        return original_open(self, mode, *args, **kwargs)

    monkeypatch.setattr(Path, 'open', require_write_handle)
    path = tmp_path / 'saved.h5'
    loader.save_results_to_h5([_block('toy', ['A'], [3.0])], path)
    assert loader.load_results_from_h5(path)[0]['problem_names'] == ['A']


def test_failed_h5_publish_preserves_existing_archive(tmp_path, monkeypatch):
    path = tmp_path / 'saved.h5'
    loader.save_results_to_h5([_block('toy', ['A'], [3.0])], path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()

    def fail_replace(source, target):
        raise OSError('injected atomic publication failure')

    monkeypatch.setattr(loader.os, 'replace', fail_replace)
    with pytest.raises(OSError, match='complete archive is preserved') as caught:
        loader.save_results_to_h5([_block('new', ['B'], [4.0])], path)
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    complete, = tmp_path.glob('.saved.h5.*.tmp')
    assert str(complete) in str(caught.value)
    assert loader.load_results_from_h5(complete)[0]['problem_names'] == ['B']
    assert isinstance(caught.value.__cause__, OSError)


def test_benchmark_save_failure_is_explicit_even_when_silent(saved_same_second_runs, monkeypatch):
    root, markers, _ = saved_same_second_runs

    def fail_save(*args, **kwargs):
        raise OSError('injected save failure')

    monkeypatch.setattr(profiles, 'save_results_to_h5', fail_save)
    with pytest.raises(RuntimeError, match='save.*experiment') as caught:
        benchmark(None, load='latest', benchmark_id='saved',
                  draw_hist_plots='none', silent=True, max_tol_order=1,
                  summarize_performance_profiles=False,
                  summarize_data_profiles=False)
    assert isinstance(caught.value.__cause__, OSError)
    assert len(list(root.glob('*/test_log/time_stamp_*.txt'))) == len(markers)


def test_summary_uses_legacy_compatible_figure_keyword(saved_same_second_runs, monkeypatch):
    root, _, _ = saved_same_second_runs
    real_figure = profiles.Figure

    def legacy_figure(*args, **kwargs):
        if 'layout' in kwargs:
            raise TypeError('Matplotlib 3.4 does not accept layout')
        return real_figure(*args, **kwargs)

    monkeypatch.setattr(profiles, 'Figure', legacy_figure)
    scores, _, _ = benchmark(None, load='latest', benchmark_id='saved',
                             draw_hist_plots='none', silent=True, max_tol_order=1,
                             summarize_performance_profiles=True,
                             summarize_data_profiles=True)
    assert scores.shape == (2,)
    assert list(root.glob('*/summary_*.pdf'))
