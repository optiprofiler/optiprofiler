"""Reference witnesses are optional metadata, never hidden oracle queries."""
import pickle
import json

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem, ProblemReference, benchmark
from optiprofiler.loader import load_results_from_h5
from optiprofiler.tests.eval_report_contract import assert_valid


def quadratic(x):
    return float(np.dot(x, x))


def test_reference_point_is_owned_optional_and_survives_serialization():
    point = np.array([1.0, 2.0])
    reference = ProblemReference(5.0, 'best_known', 'analytic', x_ref=point)
    problem = Problem(quadratic, [3.0, 4.0], reference=reference)
    point[:] = -1
    np.testing.assert_array_equal(problem.reference.x_ref, [1.0, 2.0])
    exported = problem.reference.x_ref
    exported[:] = 99
    np.testing.assert_array_equal(problem.reference.x_ref, [1.0, 2.0])
    restored = pickle.loads(pickle.dumps(problem))
    np.testing.assert_array_equal(restored.reference.x_ref, [1.0, 2.0])
    assert restored.reference == reference
    assert hash(restored.reference) == hash(reference)
    assert reference.as_dict()['x_ref'] == [1.0, 2.0]
    assert 'x_ref' not in ProblemReference(0.0, 'optimum', 'analytic').as_dict()


def affine(rng, problem):
    return np.diag([2.0, 4.0]), np.array([-1.0, 1.0]), np.diag([0.5, 0.25])


@pytest.mark.parametrize('composed', [False, True])
def test_reference_point_stays_in_original_coordinates(composed):
    reference = ProblemReference(5.0, 'best_known', 'analytic', x_ref=[1.0, 2.0])
    problem = Problem(quadratic, [3.0, 4.0], reference=reference)
    spec = 'linearly_transformed'
    feature = Feature([spec, {'name': 'noisy', 'options': {'noise_level': 0.0}}]) \
        if composed else Feature(spec)
    featured = FeaturedProblem(problem, feature, 10, 17)
    np.testing.assert_array_equal(featured.reference.x_ref, [1.0, 2.0])
    assert featured.reference.fun == 5.0
    assert featured.n_eval_fun == 0
    np.testing.assert_array_equal(problem.reference.x_ref, [1.0, 2.0])


@pytest.mark.parametrize('point', [[], [np.nan], [np.inf], [1 + 1j], [[1.0, 2.0], [3.0, 4.0]]])
def test_reference_point_rejects_invalid_coordinates(point):
    with pytest.raises((TypeError, ValueError)):
        ProblemReference(0.0, 'optimum', 'analytic', x_ref=point)


@pytest.mark.parametrize('point', [[True, False], np.array([True]), ['1', '2.5'], '12', b'1', [None, 1.0]])
def test_reference_point_follows_the_scalar_field_type_rules(point):
    # Scalar fields reject booleans and numeric strings with TypeError; the
    # point must not silently coerce them (MATLAB rejects them as well).
    with pytest.raises(TypeError, match='x_ref'):
        ProblemReference(0.0, 'optimum', 'analytic', x_ref=point)
    with pytest.raises(TypeError):
        ProblemReference(True, 'optimum', 'analytic')


def test_reference_point_accepts_integer_and_float_coordinates():
    assert ProblemReference(0.0, 'optimum', 'analytic', x_ref=np.array([1, 2], dtype=np.int64)).x_ref.tolist() == [1.0, 2.0]
    assert ProblemReference(0.0, 'optimum', 'analytic', x_ref=(0.5, np.float32(1.5))).x_ref.tolist() == [0.5, 1.5]


def test_reference_point_dimension_is_checked_before_constraint_callbacks():
    calls = []
    with pytest.raises(ValueError, match='x_ref'):
        Problem(quadratic, [1.0, 2.0], cub=lambda x: calls.append(x),
                reference=ProblemReference(0.0, 'optimum', 'analytic', x_ref=[0.0]))
    assert calls == []


def test_truth_changing_feature_drops_the_reference_point_with_its_fact():
    problem = Problem(quadratic, [3.0, 4.0],
                      reference=ProblemReference(5.0, 'best_known', 'analytic', x_ref=[1.0, 2.0]))
    featured = FeaturedProblem(problem, Feature('quantized', ground_truth=True), 10, 17)
    assert featured.reference is None


def stay(fun, x0):
    fun(x0)
    return x0


def zero(fun, x0):
    point = np.zeros_like(x0)
    fun(point)
    return point


@pytest.mark.parametrize('dimension', [2, 300])
def test_public_archive_reload_and_report_preserve_the_reference_point(tmp_path, monkeypatch, dimension):
    library = tmp_path / 'xref_library'
    library.mkdir()
    (library / 'xref_library_tools.py').write_text(
        'from optiprofiler import Problem, ProblemReference\n'
        'def xref_library_select(options): return ["P"]\n'
        'def xref_library_load(name):\n'
        f'    return Problem(lambda x: float(x@x), [2.]*{dimension}, name=name,\n'
        f'        reference=ProblemReference(0., "optimum", "analytic", x_ref=[0.]*{dimension}))\n')
    report_path = tmp_path / 'report.json'
    actual = benchmark([stay, zero], plibs=['xref_library'], custom_problem_libs_path=library,
                       n_jobs=1, n_runs=1, run_plain=False, silent=True, draw_hist_plots='none', maxdim=dimension,
                       max_tol_order=1, max_eval_factor=2, savepath=str(tmp_path), benchmark_id='saved',
                       report_path=report_path)
    archive = next(tmp_path.rglob('data_for_loading.h5'))
    stored = load_results_from_h5(str(archive))[0]
    assert json.loads(stored['reference_records'][0])[0]['x_ref'] == [0.] * dimension
    report = json.loads(report_path.read_text())
    assert_valid(report, 'eval_report')
    expected = [0.] * dimension if dimension <= 256 else {
        'dimension': dimension, 'reason': 'reference_point_not_embedded'}
    assert report['problems'][0]['reference_records'][0]['x_ref'] == expected
    # The archive search is cwd/benchmark_id; savepath only controls output.
    monkeypatch.chdir(tmp_path)
    loaded = benchmark(None, load='latest', savepath=str(tmp_path), benchmark_id='saved',
                       score_only=True, silent=True, max_tol_order=1)
    np.testing.assert_array_equal(actual[0], loaded[0])
    np.testing.assert_array_equal(actual[1], loaded[1])
