"""Metamorphic checks across geometry, reuse, channel counts and native transport."""

import itertools
import pickle

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem


def sphere(x):
    return float(np.dot(x, x))


def inequalities(x):
    return np.array([x.sum()-1.0, -x.sum()-2.0])


def equalities(x):
    return np.array([x[0]-0.25, x[-1]+0.5])


def problem(n):
    return Problem(sphere, np.linspace(0.2, 0.6, n), xl=-np.ones(n), xu=np.ones(n),
                   aub=np.ones((1, n)), bub=[0.5], aeq=np.ones((1, n)), beq=[0.2],
                   cub=inequalities, ceq=equalities)


@pytest.mark.parametrize('n', [1, 3, 7])
@pytest.mark.parametrize('order', list(itertools.permutations(['permuted', 'linearly_transformed', 'noisy'])))
def test_geometric_composition_transports_initial_oracle(n, order):
    feature = Feature('+'.join(order), noise_level=0.0, condition_factor=0.2)
    base = problem(n)
    fp = FeaturedProblem(base, feature, 8, 49)
    # No stage here perturbs x0 or reference values. Bounds may become linear
    # rows, but both observed callbacks and the reference must reach base x0.
    assert fp.fun(fp.x0) == pytest.approx(base.fun(base.x0), rel=1e-11, abs=1e-11)
    np.testing.assert_allclose(fp.cub(fp.x0), base.cub(base.x0), rtol=1e-11, atol=1e-11)
    np.testing.assert_allclose(fp.ceq(fp.x0), base.ceq(base.x0), rtol=1e-11, atol=1e-11)
    assert fp.maxcv(fp.x0) == pytest.approx(base.maxcv(base.x0), rel=1e-11, abs=1e-11)


@pytest.mark.parametrize('seed', [0, 17, 255])
@pytest.mark.parametrize('name', ['noisy+truncated+permuted',
                                 'quantized+noisy+quantized',
                                 'random_nan+noisy+truncated',
                                 'linearly_transformed+perturbed_x0+noisy'])
def test_reusing_or_serializing_spec_never_reuses_run_state(seed, name):
    spec = Feature(name)
    a = FeaturedProblem(problem(3), spec, 5, seed)
    probes = [np.array([0.25, -0.5, 0.75]), np.zeros(3)] * 2
    values = [(a.fun(x), a.cub(x).copy(), a.ceq(x).copy()) for x in probes]
    for feature in (spec, pickle.loads(pickle.dumps(spec))):
        b = FeaturedProblem(problem(3), feature, 5, seed)
        for x, expected in zip(probes, values):
            for actual, reference in zip((b.fun(x), b.cub(x), b.ceq(x)), expected):
                np.testing.assert_array_equal(actual, reference)
        for attr in ('fun_hist', 'cub_hist', 'ceq_hist', 'maxcv_hist'):
            np.testing.assert_array_equal(getattr(a, attr), getattr(b, attr))


@pytest.mark.parametrize('feature', ['plain', 'noisy', 'noisy+noisy'])
@pytest.mark.parametrize('budget', [1, 2, 5])
def test_vector_constraint_budget_counts_queries_not_components(feature, budget):
    options = {} if feature == 'plain' else {'noise_level': 0.0}
    fp = FeaturedProblem(problem(3), Feature(feature, **options), budget, 17)
    for k in range(budget):
        x = np.full(3, float(k))
        np.testing.assert_array_equal(fp.cub(x), inequalities(x))
        np.testing.assert_array_equal(fp.ceq(x), equalities(x))
        assert fp.n_eval_cub == fp.n_eval_ceq == k+1
    last_cub, last_ceq = fp.cub_hist[-1].copy(), fp.ceq_hist[-1].copy()
    for _ in range(budget):
        np.testing.assert_array_equal(fp.cub(np.ones(3)*99), last_cub)
        np.testing.assert_array_equal(fp.ceq(np.ones(3)*99), last_ceq)
    with pytest.raises(StopIteration):
        fp.cub(np.zeros(3))
    with pytest.raises(StopIteration):
        fp.ceq(np.zeros(3))
