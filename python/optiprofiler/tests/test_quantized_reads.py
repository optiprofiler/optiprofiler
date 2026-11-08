"""Which points a ``quantized`` feature reads from the oracle, and what it records."""

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem

MESH = 0.5
X0 = 0.3
SNAPPED_X0 = 0.5
X = 1.2
SNAPPED_X = 1.0


class Oracle:
    """A problem whose callbacks log, per channel, the points they were called at."""

    def __init__(self):
        self.calls = {'fun': [], 'cub': [], 'ceq': []}

    def _log(self, channel, x):
        self.calls[channel].append(float(x[0]))

    def fun(self, x):
        self._log('fun', x)
        return float(x[0] ** 2)

    def cub(self, x):
        self._log('cub', x)
        return np.array([x[0] - 2.0])

    def ceq(self, x):
        self._log('ceq', x)
        return np.array([x[0] - 3.0])

    def problem(self):
        problem = Problem(self.fun, [X0], cub=self.cub, ceq=self.ceq)
        self.clear()
        return problem

    def clear(self):
        self.calls = {channel: [] for channel in self.calls}


def quantized(ground_truth):
    return Feature('quantized', mesh_size=MESH, mesh_type='absolute', ground_truth=ground_truth)


def violation(x):
    return max(0.0, x - 2.0, abs(x - 3.0))


@pytest.mark.parametrize('ground_truth, truth_x0, truth_x', [(False, X0, X), (True, SNAPPED_X0, SNAPPED_X)])
def test_construction_reads_the_truth_point_once(ground_truth, truth_x0, truth_x):
    oracle = Oracle()
    featured = FeaturedProblem(oracle.problem(), quantized(ground_truth), 10, seed=1)
    assert oracle.calls == {'fun': [truth_x0], 'cub': [truth_x0], 'ceq': [truth_x0]}
    assert featured.fun_init == truth_x0 ** 2
    assert featured.maxcv_init == violation(truth_x0)


@pytest.mark.parametrize('ground_truth, truth_x', [(False, X), (True, SNAPPED_X)])
def test_objective_query_reads_observation_and_truth_only(ground_truth, truth_x):
    oracle = Oracle()
    featured = FeaturedProblem(oracle.problem(), quantized(ground_truth), 10, seed=1)
    oracle.clear()
    assert featured.fun([X]) == SNAPPED_X ** 2
    # The observation at the snapped point; the truth at the point of the
    # solver, unless the snapped problem is the truth, whose value it already is.
    assert oracle.calls['fun'] == ([SNAPPED_X] if ground_truth else [SNAPPED_X, X])
    # The recorded violation is the truth's: no read of the channel at another point.
    assert oracle.calls['cub'] == [truth_x]
    assert oracle.calls['ceq'] == [truth_x]
    np.testing.assert_array_equal(featured.fun_hist, [truth_x ** 2])
    np.testing.assert_array_equal(featured.maxcv_hist, [violation(truth_x)])


@pytest.mark.parametrize('channel, value', [('cub', -2.0), ('ceq', -3.0)])
@pytest.mark.parametrize('ground_truth, truth_x', [(False, X), (True, SNAPPED_X)])
def test_constraint_query_reads_observation_and_truth_only(channel, value, ground_truth, truth_x):
    oracle = Oracle()
    featured = FeaturedProblem(oracle.problem(), quantized(ground_truth), 10, seed=1)
    oracle.clear()
    query = getattr(featured, channel)
    np.testing.assert_array_equal(query([X]), [SNAPPED_X + value])
    assert oracle.calls[channel] == ([SNAPPED_X] if ground_truth else [SNAPPED_X, X])
    assert oracle.calls['fun'] == []
    np.testing.assert_array_equal(getattr(featured, channel + '_hist'), [[truth_x + value]])


@pytest.mark.parametrize('channel', ['cub', 'ceq'])
@pytest.mark.parametrize('ground_truth', [False, True])
def test_unrecorded_constraint_query_reads_the_observation_only(channel, ground_truth):
    oracle = Oracle()
    featured = FeaturedProblem(oracle.problem(), quantized(ground_truth), 10, seed=1)
    oracle.clear()
    getattr(featured, channel)([X], record_hist=False)
    assert oracle.calls[channel] == [SNAPPED_X]


@pytest.mark.parametrize('ground_truth, truth_x', [(False, X), (True, SNAPPED_X)])
def test_violation_reads_each_constraint_at_the_truth_point_once(ground_truth, truth_x):
    oracle = Oracle()
    featured = FeaturedProblem(oracle.problem(), quantized(ground_truth), 10, seed=1)
    oracle.clear()
    assert featured.maxcv([X]) == violation(truth_x)
    assert oracle.calls == {'fun': [], 'cub': [truth_x], 'ceq': [truth_x]}


def test_problem_without_constraints_reads_no_constraint_oracle():
    calls = []

    def fun(x):
        calls.append(float(x[0]))
        return float(x[0] ** 2)

    featured = FeaturedProblem(Problem(fun, [X0]), quantized(True), 10, seed=1)
    assert calls == [SNAPPED_X0]
    assert featured.cub([X]).size == 0 and featured.ceq([X]).size == 0
    assert calls == [SNAPPED_X0]


def test_cached_query_after_the_budget_reads_nothing():
    oracle = Oracle()
    featured = FeaturedProblem(oracle.problem(), quantized(True), 1, seed=1)
    featured.fun([X])
    oracle.clear()
    assert featured.fun([0.0]) == SNAPPED_X ** 2
    assert oracle.calls == {'fun': [], 'cub': [], 'ceq': []}


@pytest.mark.parametrize('ground_truth, truth_x', [(False, X), (True, SNAPPED_X)])
def test_composed_view_keeps_one_read_per_channel(ground_truth, truth_x):
    oracle = Oracle()
    stage = {'name': 'quantized', 'options': {'mesh_size': MESH, 'mesh_type': 'absolute', 'ground_truth': ground_truth}}
    noise = {'name': 'noisy', 'options': {'noise_level': 0.0}}
    featured = FeaturedProblem(oracle.problem(), Feature([stage, noise]), 10, seed=1)
    oracle.clear()
    assert featured.fun([X]) == SNAPPED_X ** 2
    # Each view reads its predecessor once: an observation and a reference.
    assert oracle.calls['fun'] == [SNAPPED_X, truth_x]
    assert oracle.calls['cub'] == [truth_x]
    assert oracle.calls['ceq'] == [truth_x]
    np.testing.assert_array_equal(featured.fun_hist, [truth_x ** 2])
    oracle.clear()
    featured.cub([X])
    assert oracle.calls['cub'] == [SNAPPED_X, truth_x]


def test_failed_snapped_point_is_observed_as_invalid_without_reading_elsewhere():
    calls = []

    def fun(x):
        calls.append(float(x[0]))
        if x[0] == SNAPPED_X:
            raise RuntimeError('no value on this mesh point')
        return float(x[0] ** 2)

    featured = FeaturedProblem(Problem(fun, [X0]), quantized(True), 10, seed=1)
    calls.clear()
    assert np.isnan(featured.fun([X]))
    assert calls == [SNAPPED_X]
    assert np.isnan(featured.fun_hist[0])
