"""The public reference contract has two kinds and an empirical fallback."""

import json

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem, ProblemReference, benchmark
from optiprofiler.tests.eval_report_contract import assert_valid


def test_two_reference_kinds_do_not_require_a_mapping_token():
    feasible = ProblemReference.from_record(dict(fun=0., maxcv=0., kind='feasible', source='analytic'))
    merit = ProblemReference.from_record(dict(fun=-1., maxcv=.25, kind='merit', source='author'))
    assert ProblemReference.KINDS == ('feasible', 'merit')
    assert feasible.as_dict() == dict(fun=0., maxcv=0., kind='feasible', source='analytic')
    assert merit.as_dict() == dict(fun=-1., maxcv=.25, kind='merit', source='author')


@pytest.mark.parametrize('old, canonical', [('lower_bound', 'feasible'), ('optimum', 'feasible'),
                                           ('best_known', 'merit'), ('target', 'merit')])
def test_earlier_provider_records_normalize_to_two_kinds(old, canonical):
    reference = ProblemReference.from_record(dict(fun=1., maxcv=0., kind=old, source='old provider',
                                                 mapping='objective_violation/1'))
    assert reference.as_dict() == dict(fun=1., maxcv=0., kind=canonical, source='old provider')


@pytest.mark.parametrize('name', ['permuted', 'linearly_transformed', 'noisy+linearly_transformed',
                                  'linearly_transformed+noisy'])
@pytest.mark.parametrize('violation', [0., .25])
def test_builtin_references_keep_raw_values_and_original_coordinates(name, violation):
    reference = ProblemReference(5., 'merit', 'author', maxcv=violation, x_ref=[1., 2.])
    problem = Problem(lambda x: float(x @ x), [3., 4.], reference=reference)
    featured = FeaturedProblem(problem, Feature(name), 10, 17)
    assert featured.reference == reference
    np.testing.assert_array_equal(featured.reference.x_ref, [1., 2.])
    assert featured.reference.fun == 5.
    assert featured.reference.maxcv == violation
    assert featured.n_eval_fun == 0


@pytest.mark.parametrize('spec', ['custom', {'name': 'custom', 'options': {'mod_x0': lambda rng, p: p.x0}},
                                 'noisy+custom', 'custom+permuted', 'quantized', 'quantized+noisy',
                                 'noisy+quantized+permuted'])
def test_custom_and_truth_quantization_fall_back_even_with_a_point(spec):
    problem = Problem(lambda x: float(x @ x), [3., 4.],
                      reference=ProblemReference(0., 'feasible', 'analytic', x_ref=[0., 0.]))
    featured = FeaturedProblem(problem, Feature(spec), 10, 17)
    assert featured.reference is None


def stay(fun, x0):
    fun(x0)
    return x0


def zero(fun, x0):
    x = np.zeros_like(x0)
    fun(x)
    return x


@pytest.mark.parametrize('name, reason', [('custom', 'custom_feature'),
                                         ('quantized+custom', 'quantized_ground_truth'),
                                         ('custom+quantized', 'custom_feature')])
def test_public_report_explains_the_first_reference_fallback(tmp_path, name, reason):
    path = tmp_path / 'report.json'
    problem = Problem(lambda x: float(x @ x), [2.],
                      reference=ProblemReference(-1., 'feasible', 'author', x_ref=[0.]))
    benchmark([stay, zero], problem=problem, feature_name=name, score_only=True,
              n_jobs=1, n_runs=1, run_plain=False, silent=True, report_path=path)
    report = json.loads(path.read_text())
    assert_valid(report, 'eval_report')
    assert report['problems'][0]['reference_records'] == [None]
    assert report['problems'][0]['reference_fallback_reasons'] == [reason]
