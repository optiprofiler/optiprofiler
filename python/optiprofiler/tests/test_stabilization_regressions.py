"""
Regressions repaired after the 2026-09 independent review of the Feature and
EvalReport candidate: public option names in validation messages, duplicate
case-variant option names, view input normalization, boolean run counts,
report metadata sanitizing and reservation rollback, physical path resolution
through symlinks, undecodable artifact names, the manifest cap, nested archive
text, deterministic problem order, the trusted-boundary edge cases, and the
replay recipe written by a ``load`` invocation.
"""

import base64
import hashlib
import json
import os
import pickle
import shutil
import sys

import h5py
import matplotlib
import numpy as np
import pytest

matplotlib.use('Agg')

from optiprofiler import Feature, FeaturedProblem, Problem, benchmark
from optiprofiler import eval_report as eval_report_module
from optiprofiler.eval_report import EvalReport
from optiprofiler.experiment import validate_n_runs
from optiprofiler.legacy_compat import (LegacyConfigurationError, archived_replay_recipe, load_options,
                                        replay_arguments)
from optiprofiler.loader import load_results_from_h5
from optiprofiler.tests.eval_report_contract import assert_valid, load_schema, validation_errors
from optiprofiler.tests.test_legacy_compat import AC4_SINGLE

ARRAYS = ('fun_histories', 'maxcv_histories', 'fun_outs', 'maxcv_outs', 'fun_inits', 'maxcv_inits', 'n_evals')


def sphere(x):
    return float(np.dot(x, x))


def stay(fun, x0):
    fun(x0)
    return x0


def zero(fun, x0):
    x = np.zeros_like(x0)
    fun(x)
    return x


def step(fun, x0):
    best_x, best_f = np.array(x0, dtype=float), fun(x0)
    for i in range(len(x0)):
        for delta in (0.5, -0.5):
            trial = best_x.copy()
            trial[i] += delta
            value = fun(trial)
            if value < best_f:
                best_x, best_f = trial, value
    return best_x


def list_probe(x, rng, predecessor):
    return predecessor.fun(list(x)) + 0.0 * predecessor.maxcv(list(x))


def wrong_size_probe(x, rng, predecessor):
    return predecessor.fun([1.0])


def _library(root, name, problems=('QUAD', 'SHIFT')):
    """A developer-path provider on disk (discovery, loading and workers stay real)."""
    directory = root / name
    directory.mkdir(parents=True)
    (directory / f'{name}_tools.py').write_text(
        'import numpy as np\n'
        'from optiprofiler import Problem\n'
        f'def {name}_select(options):\n    return {list(problems)!r}\n'
        f'def {name}_load(name):\n'
        "    if name == 'SHIFT':\n"
        "        return Problem(lambda x: float(np.sum((x - 1.0) ** 2)), [0.0, 2.0], name=name)\n"
        "    return Problem(lambda x: float(x @ x), [1.0, 2.0], name=name)\n")


def _options(tmp_path, name, **extra):
    options = dict(plibs=[name], custom_problem_libs_path=tmp_path / 'libraries', ptype='u', mindim=1, maxdim=10,
                   score_only=False, draw_hist_plots='none', silent=True, savepath=str(tmp_path),
                   solver_names=['stay', 'step'], n_jobs=1, max_eval_factor=4, max_tol_order=1, seed=7)
    options.update(extra)
    return options


def _archive(directory):
    archives = sorted(directory.rglob('data_for_loading.h5'))
    assert archives, f'no archive under {directory}'
    return archives


def _refined(directory):
    documents = []
    for path in sorted(directory.rglob('options_refined.pkl')):
        documents.append((path, load_options(path)))
    return documents


def _same_numbers(first, second):
    for key in ARRAYS:
        np.testing.assert_array_equal(first[key], second[key], err_msg=key)


# ---------------------------------------------------------------- validation and construction

@pytest.mark.parametrize('value', [float('nan'), -1.0, [0.1, 'bad'], [], np.array([[0.1]]), True])
def test_perturbation_level_errors_name_the_public_option(value):
    with pytest.raises((TypeError, ValueError)) as info:
        Feature('perturbed_x0', perturbation_level=value)
    assert '`perturbation_level`' in str(info.value)
    assert 'FeatureOption' not in str(info.value)


def test_case_variant_duplicate_option_names_are_rejected():
    with pytest.raises(ValueError, match='(?i)duplicate'):
        Feature('noisy', noise_level=0.1, NOISE_LEVEL=0.2)
    with pytest.raises(ValueError, match='(?i)duplicate'):
        Feature({'name': 'noisy', 'options': {'noise_level': 0.1, 'Noise_Level': 0.2}})
    with pytest.raises(ValueError, match='(?i)duplicate'):
        Feature(['plain', {'name': 'noisy', 'options': {'noise_level': 0.1, 'NOISE_LEVEL': 0.2}}])
    # One spelling of any case is still accepted and canonicalized.
    assert Feature('noisy', NOISE_LEVEL=0.2).stages[0].options['noise_level'] == 0.2
    with pytest.raises(ValueError, match='(?i)duplicate'):
        benchmark([stay, zero], problem=Problem(sphere, [1.0]), silent=True, score_only=True,
                  feature=[{'name': 'noisy', 'options': {'noise_level': 0.1, 'NOISE_LEVEL': 0.2}}])


def test_view_methods_normalize_input_like_problem():
    problem = Problem(sphere, [0.26, 0.74])
    feature = Feature('quantized+custom', mesh_size=0.5, mesh_type='absolute', mod_fun=list_probe)
    assert FeaturedProblem(problem, feature, 10, 7).fun([0.26, 0.74]) == pytest.approx(0.5)
    assert FeaturedProblem(problem, Feature('noisy+custom', noise_level=0.0, mod_fun=list_probe), 10, 7).fun(
        [0.26, 0.74]) == pytest.approx(sphere([0.26, 0.74]))
    bad = Feature('quantized+custom', mesh_size=0.5, mod_fun=wrong_size_probe)
    with pytest.raises(ValueError, match='must have size 2'):
        FeaturedProblem(problem, bad, 10, 7).fun([0.26, 0.74])


def test_boolean_run_counts_are_rejected():
    for value in (True, False, np.bool_(True)):
        with pytest.raises(TypeError, match='n_runs'):
            validate_n_runs(value)
    assert validate_n_runs(np.int64(3)) == 3 and validate_n_runs(2.0) == 2
    with pytest.raises(TypeError, match='n_runs'):
        benchmark([stay, zero], problem=Problem(sphere, [1.0]), n_runs=True, silent=True, score_only=True)


# ---------------------------------------------------------------- report robustness

def test_surrogate_metadata_is_sanitized_and_the_benchmark_completes(tmp_path):
    report_path = tmp_path / 'report.json'
    scores, _, _ = benchmark([stay, zero], problem=Problem(sphere, [1.0, 2.0], name='p\udcff'),
                             solver_names=['ok', 'bad\udcff'], score_only=True, n_runs=1, silent=True,
                             report_path=report_path)
    text = report_path.read_text(encoding='utf-8')
    document = json.loads(text)
    assert np.isfinite(scores).all()
    assert document['status'] == 'completed'
    assert '\udcff' not in text and '\ufffd' in text
    # A name that had to be sanitized is identified by a hash of the private text.
    assert document['problems'][0]['id'].startswith('sha256:')
    assert document['scores']['solver_names'] == ['ok', 'bad\ufffd']
    assert_valid(document, 'eval_report')
    assert_valid(json.loads(report_path.with_name('report.plot_data.json').read_text(encoding='utf-8')), 'plot_data')


def test_first_write_failure_rolls_back_only_owned_files(tmp_path, monkeypatch):
    original = EvalReport._write_owned

    def fail_main(self, target, document, owned_identity):
        if target.name == 'retry.json':
            raise OSError('injected initial main write failure')
        return original(self, target, document, owned_identity)

    monkeypatch.setattr(EvalReport, '_write_owned', fail_main)
    with pytest.raises(OSError, match='injected'):
        EvalReport(tmp_path / 'retry.json', {})
    assert not (tmp_path / 'retry.json').exists()
    assert not (tmp_path / 'retry.plot_data.json').exists()
    assert not list(tmp_path.glob('.retry*'))
    monkeypatch.setattr(EvalReport, '_write_owned', original)
    EvalReport(tmp_path / 'retry.json', {})
    assert json.loads((tmp_path / 'retry.json').read_text())['status'] == 'running'

    def foreign_then_fail(self, target, document, owned_identity):
        if target.name == 'foreign.json':
            temp = target.with_name('foreign.tmp')
            temp.write_bytes(b'foreign writer')
            os.replace(temp, target)
            raise OSError('injected failure after a foreign replacement')
        return original(self, target, document, owned_identity)

    monkeypatch.setattr(EvalReport, '_write_owned', foreign_then_fail)
    with pytest.raises(OSError, match='foreign'):
        EvalReport(tmp_path / 'foreign.json', {})
    # The main target was replaced by another writer: preserved. The companion
    # was still ours (its identity is the one we published): removed.
    assert (tmp_path / 'foreign.json').read_bytes() == b'foreign writer'
    assert not (tmp_path / 'foreign.plot_data.json').exists()
    companion = tmp_path / 'kept.plot_data.json'
    companion.write_bytes(b'owner-sentinel')
    with pytest.raises(FileExistsError):
        EvalReport(tmp_path / 'kept.json', {})
    assert companion.read_bytes() == b'owner-sentinel' and not (tmp_path / 'kept.json').exists()


def test_foreign_in_place_rewrite_and_inode_reuse_are_detected(tmp_path):
    # A foreign in-place rewrite keeps the inode, and after two foreign
    # replace-over-target operations ext4 hands the recorded inode back; an
    # identity made of device and inode alone would treat both as owned.
    target = tmp_path / 'rewritten.json'
    report = EvalReport(target, {})
    inode = target.lstat().st_ino
    target.write_text('external in-place rewrite', encoding='utf-8')
    assert target.lstat().st_ino == inode
    with pytest.raises(FileExistsError, match='ownership'):
        report._write()
    assert target.read_text(encoding='utf-8') == 'external in-place rewrite'
    companion_owner = tmp_path / 'companion.json'
    companion = tmp_path / 'companion.plot_data.json'
    report = EvalReport(companion_owner, {})
    companion.write_text('external companion rewrite', encoding='utf-8')
    with pytest.raises(FileExistsError, match='ownership'):
        report._write()
    assert companion.read_text(encoding='utf-8') == 'external companion rewrite'
    replaced = tmp_path / 'replaced.json'
    report = EvalReport(replaced, {})
    recorded = replaced.lstat().st_ino
    for k in (1, 2):
        stage = replaced.with_name('replaced.json.replacement')
        stage.write_text(f'external replacement {k}', encoding='utf-8')
        os.replace(stage, replaced)
    reused = replaced.lstat().st_ino == recorded  # true on ext4; the check must not depend on it
    with pytest.raises(FileExistsError, match='ownership'):
        report._write()
    assert replaced.read_text(encoding='utf-8') == 'external replacement 2'
    assert isinstance(reused, bool)
    # A publish of an untouched report still succeeds afterwards.
    fresh = EvalReport(tmp_path / 'fresh.json', {})
    fresh._write()
    assert json.loads((tmp_path / 'fresh.json').read_text(encoding='utf-8'))['status'] == 'running'


@pytest.mark.skipif(os.name == 'nt', reason='symbolic links need privileges on Windows')
def test_symlinked_report_parent_and_output_resolve_physically(tmp_path, monkeypatch):
    physical = tmp_path / 'physical'
    (physical / 'deep' / 'reports').mkdir(parents=True)
    (physical / 'out').mkdir()
    alias = tmp_path / 'alias'
    alias.symlink_to(physical / 'deep' / 'reports', target_is_directory=True)
    out_alias = tmp_path / 'out_alias'
    out_alias.symlink_to(physical / 'out', target_is_directory=True)
    _library(tmp_path / 'libraries', 'symlinked')
    options = _options(tmp_path, 'symlinked', benchmark_id='linked', savepath=str(out_alias), n_runs=1)
    report_path = alias / 'r.json'
    benchmark([stay, step], report_path=report_path, **options)
    report = json.loads(report_path.read_text(encoding='utf-8'))
    assert report['status'] == 'completed' and report['artifacts']
    root = report_path.parent / report['artifact_root']
    assert root.resolve() == (physical / 'out' / 'linked').resolve() / root.resolve().name
    for artifact in report['artifacts']:
        assert '..' not in artifact['path'].split('/')
        assert (root / artifact['path']).is_file(), artifact['path']
    # The physical tree can be copied as a whole; the references still resolve.
    copy = tmp_path / 'copy'
    shutil.copytree(physical, copy)
    copied_report = copy / 'deep' / 'reports' / 'r.json'
    for artifact in report['artifacts']:
        assert (copied_report.parent / report['artifact_root'] / artifact['path']).is_file()
    # A load through the alias records a source path that resolves on disk.
    archive = _archive(physical / 'out' / 'linked')[0]
    monkeypatch.chdir(out_alias)
    target = alias / 'loaded.json'
    benchmark(None, load='latest', benchmark_id='linked', savepath=str(out_alias), score_only=True, silent=True,
              max_tol_order=1, report_path=target)
    loaded = json.loads(target.read_text(encoding='utf-8'))
    source = target.parent / loaded['source']['path']
    assert source.is_file() and source.resolve() == archive.resolve()
    assert loaded['source']['sha256'] == hashlib.sha256(archive.read_bytes()).hexdigest()


@pytest.mark.skipif(os.name == 'nt' or sys.getfilesystemencoding().lower().replace('-', '') != 'utf8',
                    reason='needs a POSIX UTF-8 filesystem encoding to create an undecodable name')
def test_undecodable_artifact_name_is_recorded_with_a_reason(tmp_path):
    output = tmp_path / 'out'
    output.mkdir()
    report = EvalReport(tmp_path / 'r.json', {})
    report.configure({}, {}, Feature('plain'), output_dir=output)
    with open(os.path.join(os.fsencode(str(output)), b'bad\xff.txt'), 'wb') as stream:
        stream.write(b'x')
    (output / 'good.txt').write_text('y')
    report.finish()
    document = json.loads((tmp_path / 'r.json').read_text(encoding='utf-8'))
    bad = [a for a in document['artifacts'] if a['path'] is None]
    good = [a for a in document['artifacts'] if a['path'] == 'good.txt']
    assert len(bad) == 1 and len(good) == 1
    assert bad[0]['path_reason'] == 'artifact_name_not_utf8'
    assert bad[0]['sha256'] == hashlib.sha256(b'x').hexdigest()
    assert any(d['code'] == 'artifact_path_unavailable' for d in document['diagnostics'])
    assert_valid(document, 'eval_report')


def test_diagnostics_are_deduplicated_and_omissions_count_distinct_entries(tmp_path):
    report = EvalReport(tmp_path / 'r.json', {})
    for index in range(140):
        report.add_diagnostic('review_probe', 'runtime', index=index)
        report.add_diagnostic('review_probe', 'runtime', index=index)
    report.finish()
    document = json.loads((tmp_path / 'r.json').read_text(encoding='utf-8'))
    assert len(document['diagnostics']) == 128
    assert document['diagnostics_omitted'] == 12
    assert len({d['scope']['index'] for d in document['diagnostics']}) == 128
    assert_valid(document, 'eval_report')


def test_manifest_cap_counts_only_artifacts(tmp_path, monkeypatch):
    monkeypatch.setattr(eval_report_module, '_MAX_ARTIFACTS', 3)
    output = tmp_path / 'out'
    output.mkdir()
    for i in range(3):
        (output / f'a_pre{i}.txt').write_text('pre')
    report = EvalReport(tmp_path / 'r.json', {})
    report.configure({}, {}, Feature('plain'), output_dir=output)
    for i in range(4):
        (output / f'z_new{i}.txt').write_text('new')
    report.finish()
    document = json.loads((tmp_path / 'r.json').read_text(encoding='utf-8'))
    assert [a['path'] for a in document['artifacts']] == ['z_new0.txt', 'z_new1.txt', 'z_new2.txt']
    assert any(d['code'] == 'artifact_manifest_limit_reached' for d in document['diagnostics'])


def test_schema_requires_a_reason_for_null_paths(tmp_path):
    schema = load_schema('eval_report')
    artifact = {'path': None, 'media_type': 'text/plain', 'bytes': 1, 'sha256': '0' * 64,
                'kind': 'supporting_file', 'status': 'present'}
    partial = {'$ref': '#/$defs/artifact', '$defs': schema['$defs']}
    assert validation_errors(artifact, partial)
    assert not validation_errors({**artifact, 'path_reason': 'artifact_name_not_utf8'}, partial)
    source = {'kind': 'archive', 'name': 'data_for_loading.h5', 'path': None, 'status': 'unknown',
              'bytes': None, 'sha256': None}
    partial = {'$ref': '#/$defs/source', '$defs': schema['$defs']}
    assert validation_errors(source, partial)
    assert not validation_errors({**source, 'path_reason': 'cross_volume_relative_path_unavailable'}, partial)
    report_path = tmp_path / 'r.json'
    benchmark([stay, zero], problem=Problem(sphere, [1.0]), score_only=True, n_runs=1, silent=True,
              report_path=report_path)
    document = json.loads(report_path.read_text(encoding='utf-8'))
    assert document['artifact_root'] is None and document['artifacts'] == []
    assert not validation_errors(document, schema)
    document['artifacts'] = [{**artifact, 'path': 'x.txt'}]
    assert validation_errors(document, schema)
    document['artifact_root_reason'] = 'cross_volume_relative_path_unavailable'
    assert not validation_errors(document, schema)


# ---------------------------------------------------------------- archives and ordering

def test_run_plain_archive_stores_nested_feature_pipeline_as_text(tmp_path):
    _library(tmp_path / 'libraries', 'plain_reference')
    options = _options(tmp_path, 'plain_reference', benchmark_id='rp', n_runs=1, run_plain=True)
    benchmark([stay, step], feature_name='noisy', **options)
    archive = _archive(tmp_path / 'rp')[0]
    with h5py.File(archive, 'r') as handle:
        group = handle['plib_0']['results_plib_plain']
        assert 'feature_pipeline' in group and 'feature_pipeline_pickled' not in group
    loaded = load_results_from_h5(str(archive))[0]
    payload = json.loads(loaded['results_plib_plain']['feature_pipeline'])
    assert payload['schema'] == 'feature_pipeline-v3'
    assert payload['experiment']['role'] == 'plain_reference'


def test_problem_name_intersection_keeps_the_requested_order(tmp_path):
    _library(tmp_path / 'libraries', 'ordered')
    options = _options(tmp_path, 'ordered', benchmark_id='ord', n_runs=1, problem_names=['SHIFT', 'QUAD', 'SHIFT'])
    benchmark([stay, step], **options)
    results = load_results_from_h5(str(_archive(tmp_path / 'ord')[0]))
    assert list(results[0]['problem_names']) == ['SHIFT', 'QUAD']


# ---------------------------------------------------------------- trusted boundary edges

def test_legacy_compat_edge_cases_fail_explicitly(tmp_path):
    with pytest.raises(LegacyConfigurationError, match='feature_specification'):
        replay_arguments({'schema': 'options_refined-v2'})
    for payload in (b'', b'\x80\x04garbage', pickle.dumps({'a': 1})[:5]):
        path = tmp_path / 'bad.pkl'
        path.write_bytes(payload)
        with pytest.raises(LegacyConfigurationError, match='not a readable'):
            load_options(path)
    with pytest.raises(FileNotFoundError):
        load_options(tmp_path / 'missing.pkl')
    with pytest.raises(TypeError, match='legacy_compat'):
        pickle.loads(base64.b64decode(AC4_SINGLE))
    without_count = replay_arguments({'schema': 'options_user-v2', 'feature_name': 'noisy', 'ptype': 'u'})
    assert 'n_runs' not in without_count and without_count['problem_options'] == {'ptype': 'u'}
    assert replay_arguments({'schema': 'options_user-v2', 'feature_name': 'noisy', 'n_runs': 4})['n_runs'] == 4
    # A flat 1.x file written by a load invocation describes the load label.
    with pytest.raises(LegacyConfigurationError, match='load'):
        replay_arguments({'load': '20250101_000000', 'feature_name': 'plain', 'n_runs': 1, 'ptype': 'u'})


# ---------------------------------------------------------------- the recipe written by a load

def _source_experiment(tmp_path, name, benchmark_id='src', **kwargs):
    _library(tmp_path / 'libraries', name)
    options = _options(tmp_path, name, benchmark_id=benchmark_id)
    options.update(kwargs)
    benchmark([stay, step], **options)
    return options


def _load_experiment(tmp_path, benchmark_id='src', **kwargs):
    options = dict(load='latest', benchmark_id=benchmark_id, savepath=str(tmp_path), score_only=False,
                   draw_hist_plots='none', silent=True, max_tol_order=1)
    options.update(kwargs)
    return benchmark(None, **options)


def _load_recipe(tmp_path, benchmark_id='src'):
    documents = [(path, doc) for path, doc in _refined(tmp_path / benchmark_id) if doc.get('operation') == 'load']
    assert len(documents) >= 1
    return documents[-1]


def _source_recipe(tmp_path, benchmark_id='src'):
    documents = [(path, doc) for path, doc in _refined(tmp_path / benchmark_id) if doc.get('operation') != 'load']
    assert len(documents) == 1
    return documents[0]


def _time_stamp_of(source_path):
    markers = list(source_path.parent.glob('time_stamp_*.txt'))
    assert len(markers) == 1
    return markers[0].name[len('time_stamp_'):-len('.txt')]


@pytest.mark.parametrize('label', [{}, {'feature_name': 'truncated'}])
def test_load_recipe_records_the_archived_experiment_and_replays_it(tmp_path, monkeypatch, label):
    _source_experiment(tmp_path, 'recipe', feature_name='noisy', noise_level=0.3, n_runs=2)
    source_path, source = _source_recipe(tmp_path)
    monkeypatch.chdir(tmp_path)
    _load_experiment(tmp_path, **label)
    _, refined = _load_recipe(tmp_path)
    assert refined['schema'] == 'options_refined-v2' and refined['operation'] == 'load'
    assert refined['replayable'] is True and refined['replay_reason'] is None
    # The recipe is the archived experiment, never the load label or its defaults.
    assert refined['feature_specification'] == source['feature_specification']
    assert refined['feature_specification'][0]['options']['noise_level'] == 0.3
    assert refined['n_runs'] == 2 and refined['feature_name'] == 'noisy' and refined['feature_route'] is None
    archived = refined['archived_experiment']
    assert archived['seed'] == 7 and archived['run_plain'] is False and archived['plain_reference_n_runs'] is None
    assert archived['source_options_schema'] == 'options_refined-v2'
    assert archived['problem_options']['plibs'] == ['recipe']
    replay = replay_arguments(refined)
    assert replay['n_runs'] == 2 and replay['feature'] == source['feature_specification']
    options = _options(tmp_path, 'recipe', benchmark_id='replay', seed=archived['seed'], run_plain=archived['run_plain'])
    benchmark([stay, step], feature=replay['feature'], n_runs=replay['n_runs'], **options)
    _same_numbers(load_results_from_h5(str(_archive(tmp_path / 'src')[0]))[0],
                  load_results_from_h5(str(_archive(tmp_path / 'replay')[0]))[0])


def test_load_recipe_covers_compositions_and_the_plain_reference(tmp_path, monkeypatch):
    _source_experiment(tmp_path, 'composed', feature_name='noisy+truncated', noise_level=0.2,
                       significant_digits=4, n_runs=3, run_plain=True)
    source_path, source = _source_recipe(tmp_path)
    monkeypatch.chdir(tmp_path)
    _load_experiment(tmp_path)
    _, refined = _load_recipe(tmp_path)
    assert refined['replayable'] is True
    assert [entry['name'] for entry in refined['feature_specification']] == ['noisy', 'truncated']
    assert refined['n_runs'] == 3
    archived = refined['archived_experiment']
    # Primary and plain-reference counts stay distinct facts.
    assert archived['n_runs'] == 3 and archived['run_plain'] is True and archived['plain_reference_n_runs'] == 1
    replay = replay_arguments(refined)
    options = _options(tmp_path, 'composed', benchmark_id='replay', seed=archived['seed'], run_plain=True)
    benchmark([stay, step], feature=replay['feature'], n_runs=replay['n_runs'], **options)
    original = load_results_from_h5(str(_archive(tmp_path / 'src')[0]))[0]
    again = load_results_from_h5(str(_archive(tmp_path / 'replay')[0]))[0]
    _same_numbers(original, again)
    _same_numbers(original['results_plib_plain'], again['results_plib_plain'])


def test_load_recipe_fails_closed_when_the_source_cannot_be_recovered(tmp_path, monkeypatch):
    _source_experiment(tmp_path, 'closed', feature_name='noisy', noise_level=0.3, n_runs=2)
    source_path, source = _source_recipe(tmp_path)
    monkeypatch.chdir(tmp_path)
    # A flat 1.x source file without feature identity.
    source_path.write_bytes(pickle.dumps({'n_runs': 2, 'noise_level': 0.3, 'ptype': 'u', 'seed': 7}))
    _load_experiment(tmp_path)
    _, refined = _load_recipe(tmp_path)
    assert refined['replayable'] is False
    assert refined['replay_reason'] == 'source_options_not_replayable'
    assert refined['feature_specification'] is None and refined['n_runs'] is None
    assert refined['archived_experiment']['replayable'] is False
    with pytest.raises(LegacyConfigurationError, match='load'):
        replay_arguments(refined)
    with pytest.raises(ValueError, match='cannot be None'):
        benchmark([stay, step], feature=refined['feature_specification'], n_runs=1,
                  **_options(tmp_path, 'closed', benchmark_id='never'))
    # 'latest' now selects the load above, which carries no recipe.
    _load_experiment(tmp_path)
    _, refined = _load_recipe(tmp_path)
    assert refined['replay_reason'] == 'source_is_a_load_invocation_without_recipe'
    # A source that disagrees with the archive (tampered run count).
    stamp = _time_stamp_of(source_path)
    source_path.write_bytes(pickle.dumps({**source, 'n_runs': 3}))
    _load_experiment(tmp_path, load=stamp)
    _, refined = _load_recipe(tmp_path)
    assert refined['replay_reason'] == 'archive_run_axis_disagrees_with_source_options'
    # A missing source file.
    source_path.unlink()
    _load_experiment(tmp_path, load=stamp)
    _, refined = _load_recipe(tmp_path)
    assert refined['replay_reason'] == 'source_options_refined_missing'
    assert archived_replay_recipe(None, [])['replay_reason'] == 'source_options_refined_missing'


def custom_probe_a(x, rng, problem):
    return problem.fun(x)


def custom_probe_b(x, rng, problem):
    return problem.fun(x) + 1.0


def _archive_for(feature, n_runs=2, stamp='same-explicit-label'):
    from optiprofiler.experiment import resolve_plan
    from optiprofiler.provenance import feature_pipeline_text
    return {'fun_histories': np.zeros((1, 2, n_runs, 3)), 'feature_stamp': stamp,
            'feature_pipeline': feature_pipeline_text(feature, resolve_plan(feature, requested=n_runs), feature_stamp=stamp)}


def _native_for(feature, n_runs=2, stamp='same-explicit-label', **extra):
    from optiprofiler.provenance import effective_specification
    record = dict(replayable=True, feature_specification=effective_specification(feature), n_runs=n_runs, seed=7,
                  run_plain=False, feature_stamp=stamp, problem_options={'ptype': 'u'})
    record.update(extra)
    return record


def test_recovered_recipe_requires_matching_stage_options():
    from optiprofiler.legacy_compat import _recipe_from_archived
    actual = Feature('noisy', noise_level=0.1)
    # Same stage name, count and label, different option: not the archived experiment.
    stale = _recipe_from_archived(_native_for(Feature('noisy', noise_level=0.2)), [_archive_for(actual)])
    assert stale['replayable'] is False
    assert stale['replay_reason'] == 'archive_feature_pipeline_disagrees_with_source_options'
    assert stale['archived_experiment']['detail'] == 'stage_option:noise_level'
    assert stale['feature_specification'] is None and stale['n_runs'] is None
    # The cross-check decides whether a candidate is replayable: one without a
    # ``replayable`` flag is judged on the facts, and the stated reason is the
    # option mismatch, not the missing flag.
    candidate = {key: value for key, value in _native_for(Feature('noisy', noise_level=0.2)).items() if key != 'replayable'}
    without_flag = _recipe_from_archived(candidate, [_archive_for(actual)])
    assert without_flag['replay_reason'] == 'archive_feature_pipeline_disagrees_with_source_options'
    assert without_flag['archived_experiment']['detail'] == 'stage_option:noise_level'
    # The archive disagreement is the primary fact: it is reported even when the
    # candidate is also incomplete; completeness is asked of a candidate that agrees.
    partial = {key: value for key, value in candidate.items() if key != 'problem_options'}
    assert _recipe_from_archived(partial, [_archive_for(actual)])['archived_experiment']['detail'] == 'stage_option:noise_level'
    agreeing = {key: value for key, value in _native_for(actual).items() if key != 'problem_options'}
    assert _recipe_from_archived(agreeing, [_archive_for(actual)])['replay_reason'] == 'source_recipe_malformed:problem_options'
    same = _recipe_from_archived(_native_for(actual), [_archive_for(actual)])
    assert same['replayable'] is True
    assert same['archived_experiment']['archive_comparison'].endswith('feature_pipeline_v3_stages_and_options')
    assert same['archived_experiment']['callback_comparison'] == 'not_applicable'
    # A different composition of the same names, and a different count, are rejected too.
    reordered = _recipe_from_archived(_native_for(Feature('noisy+truncated')), [_archive_for(Feature('truncated+noisy'))])
    assert reordered['replay_reason'] == 'archive_feature_pipeline_disagrees_with_source_options'
    count = _recipe_from_archived(_native_for(actual, n_runs=3), [_archive_for(actual, n_runs=3)])
    assert count['replayable'] is True
    # Callbacks are compared by descriptor (module and name) only, and the record says so.
    with_a = Feature('custom', mod_fun=custom_probe_a)
    described = _recipe_from_archived(_native_for(with_a), [_archive_for(with_a)])
    assert described['replayable'] is True
    assert described['archived_experiment']['callback_comparison'] == 'descriptor_module_and_name_only_not_identity_or_semantics'
    other = _recipe_from_archived(_native_for(Feature('custom', mod_fun=custom_probe_b)), [_archive_for(with_a)])
    assert other['replay_reason'] == 'archive_feature_pipeline_disagrees_with_source_options'
    assert other['archived_experiment']['detail'] == 'stage_option:mod_fun'
    # Malformed recovered records fail closed with a reason, never with a bare exception.
    for broken, reason in (({'replayable': True}, 'source_recipe_malformed:feature_specification'),
                           ({**_native_for(actual), 'n_runs': None}, 'source_recipe_malformed:n_runs'),
                           ({**_native_for(actual), 'seed': True}, 'source_recipe_malformed:seed'),
                           ({**_native_for(actual), 'run_plain': 'yes'}, 'source_recipe_malformed:run_plain'),
                           ({**_native_for(actual), 'problem_options': None}, 'source_recipe_malformed:problem_options')):
        assert _recipe_from_archived(broken, [_archive_for(actual)])['replay_reason'] == reason
    with pytest.raises(LegacyConfigurationError, match='malformed'):
        replay_arguments({'operation': 'load', 'archived_experiment': {'replayable': True}})


def test_load_recipe_fails_closed_when_source_options_disagree_with_the_archive(tmp_path, monkeypatch):
    # An explicit label keeps the feature stamp identical, so only the
    # pipeline's stage options can reveal the disagreement.
    _source_experiment(tmp_path, 'stale', feature_name='noisy', noise_level=0.3, n_runs=2, feature_stamp='label')
    source_path, source = _source_recipe(tmp_path)
    monkeypatch.chdir(tmp_path)
    tampered = dict(source)
    tampered['feature_specification'] = [{'name': 'noisy', 'options': {**source['feature_specification'][0]['options'],
                                                                        'noise_level': 0.2}}]
    source_path.write_bytes(pickle.dumps(tampered))
    _load_experiment(tmp_path)
    _, refined = _load_recipe(tmp_path)
    assert refined['replayable'] is False
    assert refined['replay_reason'] == 'archive_feature_pipeline_disagrees_with_source_options'
    assert refined['archived_experiment']['detail'] == 'stage_option:noise_level'
    assert refined['feature_specification'] is None and refined['n_runs'] is None
    with pytest.raises(LegacyConfigurationError, match='load'):
        replay_arguments(refined)
    # A malformed nested recipe in a load-written source also fails closed.
    stamp = _time_stamp_of(source_path)
    source_path.write_bytes(pickle.dumps({'operation': 'load', 'load': stamp,
                                          'archived_experiment': {'replayable': True, 'n_runs': 2}}))
    _load_experiment(tmp_path, load=stamp)
    _, refined = _load_recipe(tmp_path)
    assert refined['replay_reason'] == 'source_recipe_malformed:feature_specification'


def test_load_of_a_load_adopts_the_recovered_recipe(tmp_path, monkeypatch):
    _source_experiment(tmp_path, 'nested', feature_name='noisy', noise_level=0.3, n_runs=2)
    monkeypatch.chdir(tmp_path)
    _load_experiment(tmp_path)
    first_path, first = _load_recipe(tmp_path)
    _load_experiment(tmp_path)
    second_path, second = _load_recipe(tmp_path)
    assert second_path != first_path
    assert second['replayable'] is True
    assert second['feature_specification'] == first['feature_specification'] and second['n_runs'] == 2
