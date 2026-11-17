"""A load names the archived feature in the summary title, the README and the log."""
import json

import h5py
import numpy as np
import pytest
from pypdf import PdfReader

from optiprofiler import Feature, benchmark, provenance
from optiprofiler.profile_utils import _get_default_feature_stamp
from optiprofiler.utils import FeatureName


def stay(fun, x0):
    fun(x0)
    return x0


def zero(fun, x0):
    point = np.zeros_like(x0)
    fun(point)
    return point


def no_solver(*args, **kwargs):
    pytest.fail('loading saved results must not execute solvers')


@pytest.fixture
def library(tmp_path):
    path = tmp_path / 'title_library'
    path.mkdir()
    (path / 'title_library_tools.py').write_text(
        'from optiprofiler import Problem\n'
        'def title_library_select(options): return ["P"]\n'
        'def title_library_load(name):\n'
        '    return Problem(lambda x: float(x@x), [2., 1.], name=name)\n')
    return path


def run(tmp_path, library, **feature):
    """Save one small experiment in ``tmp_path / 'saved'`` and return its directory."""
    benchmark([stay, zero], plibs=['title_library'], custom_problem_libs_path=library, n_jobs=1, n_runs=1,
              max_tol_order=1, max_eval_factor=2, draw_hist_plots='none', silent=True, savepath=str(tmp_path),
              benchmark_id='saved', plot_options={'rc': {'text.usetex': False}}, **feature)
    return next((tmp_path / 'saved').iterdir())


def load(tmp_path, monkeypatch, name, source='saved', **options):
    """Load the experiment below ``tmp_path / source`` into ``tmp_path / name``; return the new directory and the scores."""
    monkeypatch.chdir(tmp_path)
    scores = benchmark(None, load='latest', benchmark_id=source, savepath=str(tmp_path / name), max_tol_order=1, n_jobs=1,
                       draw_hist_plots='none', plot_options={'rc': {'text.usetex': False}}, **options)[0]
    return next((tmp_path / name / source).iterdir()), scores


def named_features(experiment):
    """The feature named by the summary title, the README and the log of an experiment (no log line if silent)."""
    summary = next(experiment.glob('summary_*.pdf'))
    # Feature names have no spaces, so the spacing of the extracted text does not matter.
    title = ''.join(PdfReader(summary).pages[0].extract_text().split())
    # The README is a table; a long file name wraps below the first column, beside the description.
    rows = (experiment / 'README.txt').read_text().splitlines()
    column = next(row for row in rows if row.startswith('File/Folder Name')).index('Description')
    readme = ' '.join(' '.join(row[column:] for row in rows).split())
    log = ' '.join((experiment / 'test_log' / 'log.txt').read_text().split())
    in_title = title.split('with"')[1].split('"feature')[0]
    in_readme = readme.split("for the feature '")[1].split("'.")[0]
    in_log = log.split('Finished creating profiles with the "')[1].split('" feature.')[0] if 'Finished' in log else None
    return in_title, in_readme, in_log


def stamp_of(experiment):
    """The folder name of an experiment without its time stamp."""
    return experiment.name[:-len('_20200101_000000')]


def drop_feature_pipeline(experiment):
    """Turn the archive into the layout of OptiProfiler 1.x, which has the stamp only."""
    with h5py.File(experiment / 'test_log' / 'data_for_loading.h5', 'r+') as archive:
        for group in archive.values():
            del group['feature_pipeline']


@pytest.mark.parametrize('feature, name', [
    (dict(), 'plain'),
    (dict(feature_name='linearly_transformed', rotated=True, condition_factor=0), 'linearly_transformed'),
    (dict(feature_name='linearly_transformed+noisy', rotated=True, noise_level=1e-2), 'linearly_transformed+noisy'),
    (dict(feature=[{'name': 'noisy', 'options': {'noise_level': 1e-2}}, 'truncated']), 'noisy+truncated'),
    (dict(feature_name='plain+quantized'), 'quantized'),
    # The archive names the feature, whatever stamp the user chose for the folder.
    (dict(feature_name='noisy', feature_stamp='my_noise'), 'noisy'),
])
def test_load_without_feature_name_names_the_archived_feature(tmp_path, monkeypatch, library, feature, name):
    saved = run(tmp_path, library, **feature)
    assert named_features(saved) == (name, name, None)
    loaded, scores = load(tmp_path, monkeypatch, 'loaded')
    assert named_features(loaded) == (name, name, name)
    # The output of a load is an archive of the same experiment.
    again, same_scores = load(tmp_path, monkeypatch, 'again', source='loaded', silent=True)
    assert named_features(again) == (name, name, None)
    np.testing.assert_array_equal(same_scores, scores)


MISMATCH = 'given to load differs from the feature of the loaded experiment'


@pytest.mark.parametrize('label', ['noisy', 'plain', 'noisy+truncated'])
def test_feature_name_given_to_a_load_stays_its_label(tmp_path, monkeypatch, capsys, library, label):
    saved = run(tmp_path, library, feature_name='linearly_transformed', rotated=True, condition_factor=0)
    archive = (saved / 'test_log' / 'data_for_loading.h5').read_bytes()
    monkeypatch.setattr('optiprofiler.profiles._solve_all_problems', no_solver)
    _, scores = load(tmp_path, monkeypatch, 'unnamed', silent=True)
    capsys.readouterr()
    loaded, labelled_scores = load(tmp_path, monkeypatch, 'loaded', feature_name=label)
    assert named_features(loaded) == (label, label, label)
    printed = ' '.join(capsys.readouterr().out.split())
    assert f'The feature_name "{label}" {MISMATCH}, "linearly_transformed".' in printed
    # The label changes neither the results nor the source, nor the stamp of the new folder.
    np.testing.assert_array_equal(labelled_scores, scores)
    assert (saved / 'test_log' / 'data_for_loading.h5').read_bytes() == archive
    assert '_linearly_transformed_rotated_' in loaded.name


@pytest.mark.parametrize('feature, label, pipeline, warned', [
    (dict(feature_name='linearly_transformed', rotated=True), 'linearly_transformed', True, False),
    (dict(feature_name='linearly_transformed+noisy'), 'linearly_transformed+noisy', True, False),
    # The declared and the effective name both match.
    (dict(feature_name='plain+quantized'), 'plain+quantized', True, False),
    (dict(feature_name='plain+quantized'), 'quantized', True, False),
    (dict(feature_name='plain+quantized'), 'noisy', True, True),
    # An archive of OptiProfiler 1.x is compared through its stamp.
    (dict(feature_name='noisy'), 'noisy', False, False),
    (dict(feature_name='noisy'), 'plain', False, True),
    # A stamp set by the user names no feature to compare with.
    (dict(feature_name='noisy', feature_stamp='my_noise'), 'plain', False, False),
])
def test_feature_name_given_to_a_load_is_compared_with_the_archive(tmp_path, monkeypatch, capsys, library, feature,
                                                                   label, pipeline, warned):
    saved = run(tmp_path, library, **feature)
    if not pipeline:
        drop_feature_pipeline(saved)
    capsys.readouterr()
    loaded, _ = load(tmp_path, monkeypatch, 'loaded', feature_name=label, silent=True)
    assert (MISMATCH in ' '.join(capsys.readouterr().out.split())) == warned
    assert named_features(loaded)[0] == Feature(label).name


@pytest.mark.parametrize('feature, name', [
    (dict(), 'plain'),
    (dict(feature_name='linearly_transformed', rotated=True, condition_factor=0), 'linearly_transformed'),
    (dict(feature_name='noisy', noise_level=1e-2), 'noisy'),
    (dict(feature_name='truncated', significant_digits=4), 'truncated'),
    # A stamp set by the user is the only name that such an archive has.
    (dict(feature_name='noisy', feature_stamp='my_noise'), 'my_noise'),
])
def test_load_of_an_archive_without_feature_pipeline_reads_the_stamp(tmp_path, monkeypatch, library, feature, name):
    saved = run(tmp_path, library, **feature)
    drop_feature_pipeline(saved)
    loaded, _ = load(tmp_path, monkeypatch, 'loaded')
    assert named_features(loaded) == (name, name, name)
    # The new folder keeps the archived stamp, as before.
    assert stamp_of(loaded) == stamp_of(saved)


def test_load_of_an_archive_without_feature_stamp_keeps_the_default_label(tmp_path, monkeypatch, library):
    saved = run(tmp_path, library, feature_name='noisy')
    with h5py.File(saved / 'test_log' / 'data_for_loading.h5', 'r+') as archive:
        for group in archive.values():
            del group['feature_pipeline'], group['feature_stamp']
    assert named_features(load(tmp_path, monkeypatch, 'loaded')[0]) == ('plain', 'plain', 'plain')


def generated_stamps():
    """Stamps generated for every feature: defaults, other option values, and the conventions of MATLAB."""
    cases = [(Feature(name.value), name.value) for name in FeatureName]
    cases += [
        (Feature('perturbed_x0', perturbation_level=0.5, distribution='gaussian'), 'perturbed_x0'),
        (Feature('perturbed_x0', perturbation_level=[0.1, 0.2]), 'perturbed_x0'),
        (Feature('noisy', noise_level=1e-2, noise_type='absolute', distribution='uniform'), 'noisy'),
        (Feature('noisy', noise_mode='deterministic'), 'noisy'),
        (Feature('truncated', significant_digits=4, perturbed_trailing_digits=True), 'truncated'),
        (Feature('linearly_transformed', rotated=False), 'linearly_transformed'),
        (Feature('linearly_transformed', rotated=False, condition_factor=2), 'linearly_transformed'),
        (Feature('linearly_transformed', condition_factor=1.5), 'linearly_transformed'),
        (Feature('random_nan', nan_rate=0.2), 'random_nan'),
        (Feature('unrelaxable_constraints', unrelaxable_bounds=False), 'unrelaxable_constraints'),
        (Feature('unrelaxable_constraints', unrelaxable_linear_constraints=True,
                 unrelaxable_nonlinear_constraints=True), 'unrelaxable_constraints'),
        (Feature('quantized', mesh_size=0.1, ground_truth=False), 'quantized'),
        (Feature('noisy+truncated'), 'noisy+truncated'),
        (Feature('linearly_transformed+noisy+noisy'), 'linearly_transformed+noisy+noisy'),
    ]
    cases = [(_get_default_feature_stamp(feature, bounded=False), name) for feature, name in cases]
    return cases + [('linearly_transformed_rotated_2', 'linearly_transformed'), ('perturbed_x0_0.001_spherical', 'perturbed_x0'),
                    ('noisy_0.001_mixed_gaussian', 'noisy'), ('quantized_0.001_ground_truth', 'quantized')]


@pytest.mark.parametrize('stamp, name', generated_stamps())
def test_name_of_a_generated_stamp(stamp, name):
    assert provenance.archived_feature_name({'feature_stamp': stamp}) == name


@pytest.mark.parametrize('stamp', [
    'my_label', 'plain_vs_noisy', 'custom_shift', 'permuted2', 'noisy2', 'Noisy_0.001', 'noisy__mine',
    '__noisy_0.001', 'linearly_transformed_rotated__noisy_0.01_mixed_gaussian__trunc_1a2b3c4d',
])
def test_a_stamp_of_another_form_is_the_name(stamp):
    assert provenance.archived_feature_name({'feature_stamp': stamp}) == stamp


def test_feature_pipeline_names_the_feature_before_the_stamp():
    feature = Feature('plain+noisy+truncated')
    record = {'feature_stamp': 'my_label', 'feature_pipeline': provenance.feature_pipeline_text(feature, None)}
    assert provenance.archived_feature_name(record) == 'noisy+truncated'
    # The first schema was flat.
    record['feature_pipeline'] = json.dumps({'schema': 'feature_pipeline-v1', 'effective_name': 'quantized+noisy'})
    assert provenance.archived_feature_name(record) == 'quantized+noisy'
    record['feature_pipeline'] = record['feature_pipeline'].encode()
    assert provenance.archived_feature_name(record) == 'quantized+noisy'
    for damaged in ('{not json', json.dumps([1]), json.dumps({'feature': {'effective_name': 3}}),
                    json.dumps({'feature': {'effective_name': ''}}), json.dumps({'schema': 'feature_pipeline-v3'})):
        record['feature_pipeline'] = damaged
        assert provenance.archived_feature_name(record) == 'my_label'
    assert provenance.archived_feature_name({}) is None
    assert provenance.archived_feature_name({'feature_stamp': ''}) is None
    assert provenance.archived_feature_name({'feature_stamp': b'noisy_0.001_mixed_gaussian'}) == 'noisy'


def test_names_that_an_archive_vouches_for():
    record = {'feature_stamp': 'my_label', 'feature_pipeline': provenance.feature_pipeline_text(Feature('plain+quantized'), None)}
    assert provenance.archived_feature_names(record) == {'quantized', 'plain+quantized'}
    assert provenance.archived_feature_names({'feature_stamp': 'linearly_transformed_rotated__noisy_0.01_mixed_gaussian'}) == \
        {'linearly_transformed+noisy'}
    assert provenance.archived_feature_names({'feature_stamp': b'noisy_0.001_mixed_gaussian'}) == {'noisy'}
    for nothing in ({}, {'feature_stamp': ''}, {'feature_stamp': 'my_label'}, {'feature_stamp': 'my_label', 'feature_pipeline': '{not json'}):
        assert provenance.archived_feature_names(nothing) == set()
