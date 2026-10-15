"""Reference facts must change scoring, not just survive construction."""
import json

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem, ProblemReference, benchmark
from optiprofiler.loader import load_results_from_h5
from optiprofiler.profile_utils import process_results, _default_merit
from optiprofiler.utils import ProfileOption
from optiprofiler.eval_report import schema_resource
from optiprofiler.tests.eval_report_contract import assert_valid


def raw_reference(fun=0.0, maxcv=0.0, kind='optimum'):
    return dict(fun=fun, maxcv=maxcv, kind=kind, source='analytic test',
                mapping='objective_violation/1')


def results(reference=None, runs=1):
    history = np.full((1, 2, runs, 2), 1.0)
    return dict(problem_names=['P'], problem_dims=np.array([1]),
                merit_histories=history, merit_outs=history[:, :, :, -1],
                merit_inits=np.full((1, runs), 4.0),
                maxcv_inits=np.arange(1, runs + 1, dtype=float)[None, :],
                n_evals=np.full((1, 2, runs), 2),
                reference_records=[json.dumps([reference] * runs)])


@pytest.mark.parametrize('kind', ProblemReference.KINDS)
def test_raw_reference_is_structural_and_has_no_point(kind):
    reference = ProblemReference.from_record(raw_reference(kind=kind))
    assert reference.fun == reference.maxcv == 0.0
    assert reference.as_dict() == raw_reference(kind=kind)
    assert not hasattr(reference, 'point')


@pytest.mark.parametrize('kind', ['lower_bound', 'optimum'])
def test_bound_and_optimum_cannot_claim_nonzero_violation(kind):
    with pytest.raises(ValueError, match='maxcv'):
        ProblemReference.from_record(raw_reference(maxcv=0.1, kind=kind))


@pytest.mark.parametrize('value', [-1.0, np.nan, np.inf])
def test_reference_violation_must_be_finite_nonnegative(value):
    with pytest.raises(ValueError, match='maxcv'):
        ProblemReference.from_record(raw_reference(maxcv=value, kind='target'))


def test_explicit_reference_replaces_cohort_minimum():
    options = {ProfileOption.MERIT_FUN: lambda f, cv, init: f}
    implicit = results()
    legacy = dict(implicit)
    del legacy['reference_records']
    assert process_results([legacy], options)[3].item() == 1.0
    assert process_results([implicit], options)[3].item() == 1.0
    assert process_results([results(raw_reference())], options)[3].item() == 0.0
    # A user-selected target is not lowered to the cohort minimum either.
    assert process_results([results(raw_reference(2.0, kind='target'))], options)[3].item() == 2.0


def test_raw_reference_uses_current_merit_and_each_run_initial_violation():
    record = raw_reference(fun=2.0, maxcv=0.5, kind='best_known')
    options = {ProfileOption.MERIT_FUN: lambda f, cv, init: 2 * f + 10 * cv + init}
    np.testing.assert_array_equal(process_results([results(record, runs=2)], options)[3], [[10.0, 11.0]])


def test_plain_baseline_cannot_override_explicit_reference():
    result = results(raw_reference(2.0, kind='target'))
    result['results_plib_plain'] = dict(problem_names=['P'],
                                      merit_histories=np.full((1, 2, 1, 2), -5.0),
                                      merit_inits=np.array([[4.0]]))
    options = {ProfileOption.MERIT_FUN: lambda f, cv, init: f,
               ProfileOption.RUN_PLAIN: True}
    assert process_results([result], options)[3].item() == 2.0


def test_invalid_converted_reference_is_rejected_not_silently_ignored():
    with pytest.raises(ValueError, match='reference'):
        process_results([results(raw_reference())],
                        {ProfileOption.MERIT_FUN: lambda f, cv, init: np.nan})


def test_default_merit_rejects_a_highly_infeasible_raw_target():
    with pytest.raises(ValueError, match='finite real scalar'):
        process_results([results(raw_reference(maxcv=1e3, kind='target'))],
                        {ProfileOption.MERIT_FUN: _default_merit})


def test_nonzero_violation_is_not_blindly_transported_through_affine():
    reference = ProblemReference.from_record(raw_reference(maxcv=0.5, kind='best_known'))
    problem = Problem(lambda x: float(x @ x), [1.0, 2.0], reference=reference)
    assert FeaturedProblem(problem, Feature('noisy+permuted'), 10, 7).reference == reference
    assert FeaturedProblem(problem, Feature('linearly_transformed'), 10, 7).reference is None
    # Feasibility itself is invariant; a zero-violation bound still survives.
    problem = Problem(lambda x: float(x @ x), [1.0, 2.0], reference=raw_reference())
    assert FeaturedProblem(problem, Feature('linearly_transformed'), 10, 7).reference is not None


def finish(fun, x0):
    x = np.ones_like(x0)
    fun(x)
    return x


def stay(fun, x0):
    fun(x0)
    return x0


def square_merit(fun, maxcv, initial_maxcv):
    return fun * fun + maxcv


@pytest.mark.parametrize('jobs', [1, 2])
def test_public_benchmark_save_load_uses_supplied_targets(tmp_path, monkeypatch, jobs):
    monkeypatch.chdir(tmp_path)
    library = tmp_path / 'libraries' / 'reference_scoring'
    library.mkdir(parents=True)
    (library / 'reference_scoring_tools.py').write_text(
        "import numpy as np\nfrom optiprofiler import Problem\n"
        "def objective(x): return float(x @ x)\n"
        "def reference_scoring_select(options): return ['REF', 'UNKNOWN']\n"
        "def reference_scoring_load(name):\n"
        "    reference = None if name == 'UNKNOWN' else " + repr(raw_reference()) + "\n"
        "    return Problem(objective, [2.0], name=name, reference=reference)\n",
        encoding='utf-8')
    original = benchmark([finish, stay], plibs=['reference_scoring'],
                         custom_problem_libs_path=str(library.parent),
                         ptype='u', mindim=1, maxdim=1, n_jobs=jobs,
                         n_runs=1, max_eval_factor=2, max_tol_order=2,
                         draw_hist_plots='none', silent=True, run_plain=True,
                         benchmark_id='saved', seed=7,
                         report_path=str(tmp_path / 'saved-report.json'))
    report = json.loads((tmp_path / 'saved-report.json').read_text(encoding='utf-8'))
    assert_valid(report, schema_resource('eval_report', 2))
    reported = {problem['name']: problem['reference_records'] for problem in report['problems']
                if problem['role'] == 'primary'}
    assert reported == {'REF': [raw_reference()], 'UNKNOWN': [None]}
    archives = list((tmp_path / 'saved').rglob('data_for_loading.h5'))
    assert len(archives) == 1
    archive = load_results_from_h5(str(archives[0]))[0]
    rows = dict(zip(archive['problem_names'], archive['reference_records']))
    assert json.loads(rows['REF']) == [raw_reference()]
    assert json.loads(rows['UNKNOWN']) == [None]
    targets = process_results([archive], {ProfileOption.MERIT_FUN: lambda f, cv, init: f,
                                         ProfileOption.RUN_PLAIN: True})[3]
    np.testing.assert_array_equal(targets[:, 0], [0.0 if name == 'REF' else 1.0
                                                for name in archive['problem_names']])
    # No solvers supplied and no provider reload: both scoring policies persist.
    loaded = benchmark(None, load='latest', benchmark_id='saved', silent=True,
                       score_only=True, max_tol_order=2)
    np.testing.assert_array_equal(loaded[0], original[0])
    np.testing.assert_array_equal(loaded[1], original[1])
    nonlinear = benchmark(None, load='latest', benchmark_id='saved', silent=True,
                          score_only=True, max_tol_order=2, merit_fun=square_merit)
    assert not np.array_equal(nonlinear[1], original[1])
    # Load filtering must retain the matching raw fact, not the adjacent row.
    benchmark(None, load='latest', benchmark_id='saved', silent=True,
              score_only=True, max_tol_order=2, problem_names=['REF'])


def test_report_observes_reference_without_an_extra_merit_call(tmp_path):
    calls = []

    def merit(fun, cv, initial):
        calls.append((fun, cv, initial))
        return fun + cv

    options = dict(problem=Problem(lambda x: float(x @ x), [2.0], reference=raw_reference()),
                   score_only=True, silent=True, n_jobs=1, merit_fun=merit,
                   max_eval_factor=2)
    expected = benchmark([finish, stay], **options)
    count = len(calls)
    calls.clear()
    actual = benchmark([finish, stay], report_path=str(tmp_path / 'reference.json'), **options)
    assert len(calls) == count
    np.testing.assert_array_equal(actual[0], expected[0])
    report = json.loads((tmp_path / 'reference.json').read_text(encoding='utf-8'))
    assert_valid(report, schema_resource('eval_report', 2))
    assert report['problems'][0]['reference_records'] == [raw_reference()]
