"""A load names the archived feature in the summary title and the log."""
import numpy as np
import pytest

from optiprofiler import Feature, benchmark, profile_utils, profiles
from optiprofiler.utils import FeatureName


def stay(fun, x0):
    fun(x0)
    return x0


def zero(fun, x0):
    point = np.zeros_like(x0)
    fun(point)
    return point


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


@pytest.fixture
def titled(monkeypatch):
    """The feature names that the summary titles of the test are made of."""
    names = []
    format_title = profiles._format_feature_title
    monkeypatch.setattr(profiles, '_format_feature_title',
                        lambda name, context: names.append(name) or format_title(name, context))
    return names


def run(tmp_path, library, **feature):
    """Save one small experiment in ``tmp_path / 'saved'`` and return its directory."""
    benchmark([stay, zero], plibs=['title_library'], custom_problem_libs_path=library, n_jobs=1, n_runs=1,
              max_tol_order=1, max_eval_factor=2, draw_hist_plots='none', silent=True, savepath=str(tmp_path),
              benchmark_id='saved', **feature)
    return next((tmp_path / 'saved').iterdir())


def load(tmp_path, monkeypatch, name, **options):
    """Load the saved experiment into ``tmp_path / name``; return the new directory and the scores."""
    monkeypatch.chdir(tmp_path)
    scores = benchmark(None, load='latest', benchmark_id='saved', savepath=str(tmp_path / name), max_tol_order=1,
                       n_jobs=1, draw_hist_plots='none', **options)[0]
    return next((tmp_path / name / 'saved').iterdir()), scores


def logged_feature(experiment):
    log = ' '.join((experiment / 'test_log' / 'log.txt').read_text().split())
    return log.split('Finished creating profiles with the "')[1].split('" feature.')[0]


@pytest.mark.parametrize('feature, name', [
    (dict(), 'plain'),
    (dict(feature_name='linearly_transformed', rotated=True, condition_factor=0), 'linearly_transformed'),
    (dict(feature_name='noisy', noise_level=1e-2), 'noisy'),
    (dict(feature_name='truncated', significant_digits=4), 'truncated'),
    # A stamp set by the user is the only name that the saved experiment has.
    (dict(feature_name='noisy', feature_stamp='my_noise'), 'my_noise'),
])
def test_load_without_feature_name_names_the_archived_feature(tmp_path, monkeypatch, library, titled, feature, name):
    saved = run(tmp_path, library, **feature)
    loaded, _ = load(tmp_path, monkeypatch, 'loaded')
    assert titled == [feature.get('feature_name', 'plain'), name]
    assert logged_feature(loaded) == name
    # The new folder keeps the archived stamp, as before.
    assert loaded.name[:-len('_20200101_000000')] == saved.name[:-len('_20200101_000000')]


@pytest.mark.parametrize('label', ['noisy', 'plain'])
def test_feature_name_given_to_a_load_stays_its_label(tmp_path, monkeypatch, library, titled, label):
    saved = run(tmp_path, library, feature_name='linearly_transformed', rotated=True, condition_factor=0)
    archive = (saved / 'test_log' / 'data_for_loading.h5').read_bytes()
    _, scores = load(tmp_path, monkeypatch, 'unnamed', silent=True)
    loaded, labelled_scores = load(tmp_path, monkeypatch, 'loaded', feature_name=label)
    assert titled == ['linearly_transformed', 'linearly_transformed', label]
    assert logged_feature(loaded) == label
    # The label changes neither the results nor the source, nor the stamp of the new folder.
    np.testing.assert_array_equal(labelled_scores, scores)
    assert (saved / 'test_log' / 'data_for_loading.h5').read_bytes() == archive
    assert '_linearly_transformed_rotated_' in loaded.name


def default_stamps():
    """Default stamps of every feature: default options, other option values, and the conventions of MATLAB."""
    cases = [Feature(name.value) for name in FeatureName]
    cases += [
        Feature('perturbed_x0', perturbation_level=0.5, distribution='gaussian'),
        Feature('noisy', noise_level=1e-2, noise_type='absolute', distribution='uniform'),
        Feature('noisy', noise_mode='deterministic'),
        Feature('truncated', significant_digits=4, perturbed_trailing_digits=True),
        Feature('linearly_transformed', rotated=False),
        Feature('linearly_transformed', rotated=False, condition_factor=2),
        Feature('linearly_transformed', condition_factor=1.5),
        Feature('random_nan', nan_rate=0.2),
        Feature('unrelaxable_constraints', unrelaxable_linear_constraints=True, unrelaxable_nonlinear_constraints=True),
        Feature('quantized', mesh_size=0.1, ground_truth=False),
    ]
    cases = [(profile_utils._get_default_feature_stamp(feature), feature.name) for feature in cases]
    return cases + [('linearly_transformed_rotated_2', 'linearly_transformed'), ('perturbed_x0_0.001_spherical', 'perturbed_x0'),
                    ('noisy_0.001_mixed_gaussian', 'noisy'), ('quantized_0.001_ground_truth', 'quantized'),
                    # A composition of OptiProfiler 2.
                    ('linearly_transformed_rotated__noisy_0.01_mixed_gaussian', 'linearly_transformed+noisy')]


@pytest.mark.parametrize('stamp, name', default_stamps())
def test_name_of_a_default_stamp(stamp, name):
    assert profile_utils._get_feature_name_from_stamp(stamp) == name


@pytest.mark.parametrize('stamp', ['my_label', 'plain_vs_noisy', 'custom_shift', 'permuted2', 'noisy2', 'Noisy_0.001',
                                   'noisy__mine', '__noisy_0.001'])
def test_a_stamp_of_another_form_is_the_name(stamp):
    assert profile_utils._get_feature_name_from_stamp(stamp) == stamp
