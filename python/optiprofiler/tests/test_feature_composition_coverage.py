"""
Acceptance coverage of compositions over the ten predefined features: every
ordered pair (repetitions and plain included) on u/b/l/n fixtures, all-ten
sequences in several orders, and long chains without a stage-count cap.

These tests check invariants that hold for every chain: transported
structure keeps its shapes, observed values are floats that may honestly be
NaN or infinite, the scoring reference of these smooth fixtures is finite,
the budget advances by one per query, and termination follows the budget.
"""

import random
import time

import numpy as np
import pytest

from optiprofiler.opclasses import Feature, FeaturedProblem, Problem

BUILTINS = ['plain', 'perturbed_x0', 'noisy', 'truncated', 'permuted', 'linearly_transformed', 'random_nan',
            'unrelaxable_constraints', 'nonquantifiable_constraints', 'quantized']
KINDS = ['u', 'b', 'l', 'n']
MAX_EVAL = 3


def quadratic(z):
    z = np.asarray(z, dtype=float)
    return float(np.sum(np.array([1.0, 2.0, 0.5]) * (z - np.array([0.5, -1.0, 2.0])) ** 2) + 0.5)


def cub_two(z):
    return np.array([z[0] ** 2 + z[1] - 1.5, z[2] - 3.0])


def ceq_one(z):
    return np.array([z[0] * z[1] - 0.2])


def fixture(kind):
    x0 = np.array([0.3, -1.2, 2.5])
    if kind == 'u':
        return Problem(quadratic, x0, name='u')
    xl, xu = np.array([-2.0, -3.0, -1.0]), np.array([2.0, 3.0, 4.0])
    if kind == 'b':
        return Problem(quadratic, x0, xl=xl, xu=xu, name='b')
    linear = dict(aub=np.array([[1.0, 1.0, 0.0]]), bub=np.array([2.5]), aeq=np.array([[0.0, 1.0, -1.0]]),
                  beq=np.array([-0.5]))
    if kind == 'l':
        return Problem(quadratic, x0, xl=xl, xu=xu, name='l', **linear)
    return Problem(quadratic, x0, xl=xl, xu=xu, name='n', cub=cub_two, ceq=ceq_one, **linear)


def exercise(featured, kind):
    n = 3
    assert featured.n == n and featured.x0.shape == (n,)
    assert featured.xl.shape == (n,) and featured.xu.shape == (n,)
    assert featured.aub.shape == (featured.bub.size, n) and featured.aeq.shape == (featured.beq.size, n)
    assert featured.ptype in 'ubln'
    if kind == 'n':
        assert featured.m_nonlinear_ub == 2 and featured.m_nonlinear_eq == 1 and featured.ptype == 'n'
    else:
        assert featured.m_nonlinear_ub == 0 and featured.m_nonlinear_eq == 0
    assert np.isfinite(featured.fun_init)
    assert featured.maxcv_init >= 0.0
    points = [featured.x0, featured.x0, featured.x0 + np.array([0.1, -0.1, 0.1])]
    for x in points:
        value = featured.fun(x)
        assert isinstance(value, float)
    assert featured.n_eval_fun == 3
    assert featured.fun_hist.shape == (3,) and np.all(np.isfinite(featured.fun_hist))
    assert featured.maxcv_hist.shape == (3,) and np.all(featured.maxcv_hist >= 0.0)
    if kind == 'n':
        c = featured.cub(points[0])
        e = featured.ceq(points[0])
        assert c.shape == (2,) and e.shape == (1,)
        assert featured.cub_hist.shape == (1, 2) and featured.ceq_hist.shape == (1, 1)
        assert np.all(np.isfinite(featured.cub_hist)) and np.all(np.isfinite(featured.ceq_hist))
    # Beyond the budget the last value is returned, then termination is raised.
    last = featured.fun(points[0])
    assert featured.n_eval_fun == MAX_EVAL
    with pytest.raises(StopIteration):
        for _ in range(2 * MAX_EVAL):
            featured.fun(points[0])
    assert featured.n_eval_fun == MAX_EVAL
    return last


@pytest.mark.parametrize('kind', KINDS)
@pytest.mark.parametrize('second', BUILTINS)
@pytest.mark.parametrize('first', BUILTINS)
def test_every_ordered_pair(first, second, kind):
    feature = Feature(f'{first}+{second}')
    assert isinstance(feature.is_stochastic, bool)
    assert feature.options['n_runs'] in (1, 5)
    exercise(FeaturedProblem(fixture(kind), feature, MAX_EVAL, 1), kind)


@pytest.mark.parametrize('kind', KINDS)
@pytest.mark.parametrize('order', ['declared', 'reversed', 'shuffled'])
def test_all_ten_in_several_orders(order, kind):
    names = list(BUILTINS)
    if order == 'reversed':
        names.reverse()
    elif order == 'shuffled':
        random.Random(2026).shuffle(names)
    feature = Feature('+'.join(names))
    assert feature.name == '+'.join(name for name in names if name != 'plain')
    assert feature.is_stochastic is True and feature.options['n_runs'] == 5
    exercise(FeaturedProblem(fixture(kind), feature, MAX_EVAL, 2), kind)


def test_repeated_stages_have_no_cap():
    exercise(FeaturedProblem(fixture('n'), Feature('noisy+noisy+noisy'), MAX_EVAL, 3), 'n')
    exercise(FeaturedProblem(fixture('n'), Feature('+'.join(['quantized'] * 12)), MAX_EVAL, 3), 'n')


@pytest.mark.parametrize('length', [32, 64])
def test_long_mixed_chains(length):
    names = [BUILTINS[i % len(BUILTINS)] for i in range(length)]
    feature = Feature('+'.join(names))
    start = time.perf_counter()
    exercise(FeaturedProblem(fixture('n'), feature, MAX_EVAL, 4), 'n')
    assert time.perf_counter() - start < 20.0
