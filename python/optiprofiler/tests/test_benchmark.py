"""Tests for the benchmark function and solve helpers in profiles module."""
import shutil
import tempfile

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pytest

from optiprofiler import benchmark
from optiprofiler.opclasses import Feature, FeaturedProblem, Problem


def rosen(x):
    """Rosenbrock function."""
    return np.sum(1e2 * (x[1:] - x[:-1] ** 2) ** 2 + (1.0 - x[:-1]) ** 2)


def simple_solver_1(fun, x0):
    """A trivial solver that returns x0."""
    return x0


def simple_solver_2(fun, x0):
    """A trivial solver that does one step of coordinate search."""
    n = len(x0)
    best_x = x0.copy()
    best_f = fun(x0)
    for i in range(n):
        x_trial = x0.copy()
        x_trial[i] += 0.1
        f_trial = fun(x_trial)
        if f_trial < best_f:
            best_f = f_trial
            best_x = x_trial.copy()
    return best_x


def problem_solver(problem):
    """A solver that accepts a Problem object."""
    x = problem.x0.copy()
    f0 = problem.fun(x)
    for i in range(problem.n):
        x_trial = x.copy()
        x_trial[i] += 0.01
        f_trial = problem.fun(x_trial)
        if f_trial < f0:
            x = x_trial
            f0 = f_trial
    return x


def _common_kwargs(tmpdir, benchmark_id='test'):
    return dict(
        plibs=['s2mpj'],
        ptype='u',
        mindim=2,
        maxdim=2,
        max_eval_factor=10,
        benchmark_id=benchmark_id,
        savepath=tmpdir,
        n_jobs=1,
        silent=True,
        draw_hist_plots='none',
        problem_names=['ROSENBR'],
    )


class TestBenchmarkBasic:

    def test_two_solvers_plain(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, profile_scores, curves = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                **_common_kwargs(tmpdir, 'test_plain'),
            )
            assert isinstance(scores, np.ndarray)
            assert scores.shape[0] == 2

    def test_problem_interface_solver(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, profile_scores, curves = benchmark(
                [problem_solver, simple_solver_1],
                feature_name='plain',
                **_common_kwargs(tmpdir, 'test_prob_iface'),
            )
            assert isinstance(scores, np.ndarray)
            assert scores.shape[0] == 2

    def test_returns_curves(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, profile_scores, curves = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                **_common_kwargs(tmpdir, 'test_curves'),
            )
            assert curves is not None
            assert isinstance(curves, list)
            assert len(curves) > 0


class TestBenchmarkFeatures:

    @pytest.fixture(autouse=True)
    def setup_tmpdir(self):
        self.tmpdir = tempfile.mkdtemp()
        yield
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _run_benchmark(self, feature_name, **extra_kwargs):
        kwargs = _common_kwargs(self.tmpdir, f'test_{feature_name}')
        kwargs.update(extra_kwargs)
        return benchmark(
            [simple_solver_1, simple_solver_2],
            feature_name=feature_name,
            **kwargs,
        )

    def test_noisy(self):
        scores, _, _ = self._run_benchmark('noisy', n_runs=2)
        assert isinstance(scores, np.ndarray)
        assert scores.shape[0] == 2

    def test_truncated(self):
        scores, _, _ = self._run_benchmark('truncated')
        assert isinstance(scores, np.ndarray)

    def test_perturbed_x0(self):
        scores, _, _ = self._run_benchmark('perturbed_x0', n_runs=2)
        assert isinstance(scores, np.ndarray)

    def test_quantized(self):
        scores, _, _ = self._run_benchmark('quantized')
        assert isinstance(scores, np.ndarray)


class TestBenchmarkErrors:

    def test_no_solvers_no_load(self):
        with pytest.raises((TypeError, ValueError)):
            benchmark(None)

    def test_invalid_feature_name(self):
        with pytest.raises((TypeError, ValueError)):
            benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='invalid_feature',
            )

    def test_single_solver_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(ValueError, match='At least two'):
                benchmark(
                    [simple_solver_1],
                    feature_name='plain',
                    **_common_kwargs(tmpdir, 'test_single'),
                )

    def test_solvers_not_callable(self):
        with pytest.raises(TypeError):
            benchmark(['not_callable', 'also_not_callable'])


class TestBenchmarkOptions:

    def test_custom_solver_names(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                solver_names=['MySolver1', 'MySolver2'],
                **_common_kwargs(tmpdir, 'test_names'),
            )
            assert isinstance(scores, np.ndarray)

    def test_custom_max_tol_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                max_tol_order=3,
                **_common_kwargs(tmpdir, 'test_tol'),
            )
            assert isinstance(scores, np.ndarray)

    def test_custom_seed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                seed=42,
                **_common_kwargs(tmpdir, 'test_seed'),
            )
            assert isinstance(scores, np.ndarray)
