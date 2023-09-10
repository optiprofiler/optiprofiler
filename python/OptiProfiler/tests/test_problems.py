import numpy as np
import pytest

from OptiProfiler import Problem


class TestProblem:

    @staticmethod
    def rosen(x):
        return np.sum(100 * (x[1:] - x[:-1] ** 2) ** 2 + (1 - x[:-1]) ** 2)

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_unconstrained_problem(self, n):
        # Construct an unconstrained problem.
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)

        # Check the problem attributes.
        assert problem.n == n
        assert problem.m_linear_ub == 0
        assert problem.m_linear_eq == 0
        assert problem.m_nonlinear_ub == 0
        assert problem.m_nonlinear_eq == 0
        assert problem.n_fun_eval == 0
        assert problem.xl.shape == (n,)
        assert problem.xu.shape == (n,)
        assert problem.aub.shape == (0, n)
        assert problem.bub.shape == (0,)
        assert problem.aeq.shape == (0, n)
        assert problem.beq.shape == (0,)
        assert problem.fun_hist.shape == (0,)
        assert problem.maxcv_hist.shape == (0,)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.xl, np.full(n, -np.inf))
        np.testing.assert_array_equal(problem.xu, np.full(n, np.inf))

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)
        assert problem.n_fun_eval == 1
        np.testing.assert_array_equal(problem.fun_hist, [self.rosen(x0)])
        np.testing.assert_array_equal(problem.maxcv_hist, [0.0])

        # Perform a second objective function evaluation.
        assert problem.fun(problem.x0 + 1.0) == self.rosen(x0 + 1.0)
        assert problem.n_fun_eval == 2
        np.testing.assert_array_equal(problem.fun_hist, [self.rosen(x0), self.rosen(x0 + 1.0)])
        np.testing.assert_array_equal(problem.maxcv_hist, [0.0, 0.0])

        # Evaluate the nonlinear constraints.
        assert problem.cub(problem.x0).shape == (0,)
        assert problem.ceq(problem.x0).shape == (0,)
        assert problem.n_fun_eval == 2

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_bound_constrained_problem(self, n):
        # Construct a bound-constrained problem.
        x0 = np.zeros(n)
        xl = np.full(n, -2.048)
        xu = np.full(n, 2.048)
        problem = Problem(self.rosen, x0, xl, xu)

        # Check the problem attributes.
        assert problem.n == n
        assert problem.m_linear_ub == 0
        assert problem.m_linear_eq == 0
        assert problem.m_nonlinear_ub == 0
        assert problem.m_nonlinear_eq == 0
        assert problem.n_fun_eval == 0
        assert problem.xl.shape == (n,)
        assert problem.xu.shape == (n,)
        assert problem.aub.shape == (0, n)
        assert problem.bub.shape == (0,)
        assert problem.aeq.shape == (0, n)
        assert problem.beq.shape == (0,)
        assert problem.fun_hist.shape == (0,)
        assert problem.maxcv_hist.shape == (0,)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.xl, xl)
        np.testing.assert_array_equal(problem.xu, xu)

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)
        assert problem.n_fun_eval == 1
        np.testing.assert_array_equal(problem.fun_hist, [self.rosen(x0)])
        np.testing.assert_array_equal(problem.maxcv_hist, [0.0])

        # Perform a second objective function evaluation.
        assert problem.fun(problem.x0 + 3.0) == self.rosen(x0 + 3.0)
        assert problem.n_fun_eval == 2
        np.testing.assert_array_equal(problem.fun_hist, [self.rosen(x0), self.rosen(x0 + 3.0)])
        np.testing.assert_array_almost_equal(problem.maxcv_hist, [0.0, 0.952])

        # Evaluate the nonlinear constraints.
        assert problem.cub(problem.x0).shape == (0,)
        assert problem.ceq(problem.x0).shape == (0,)
        assert problem.n_fun_eval == 2

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_linearly_constrained_problem(self, n):
        # Construct a linearly constrained problem.
        x0 = np.zeros(n)
        aub = np.ones((1, n))
        bub = 1.0
        aeq = np.arange(n)
        beq = 1.0
        problem = Problem(self.rosen, x0, aub=aub, bub=bub, aeq=aeq, beq=beq)

        # Check the problem attributes.
        assert problem.n == n
        assert problem.m_linear_ub == 1
        assert problem.m_linear_eq == 1
        assert problem.m_nonlinear_ub == 0
        assert problem.m_nonlinear_eq == 0
        assert problem.n_fun_eval == 0
        assert problem.xl.shape == (n,)
        assert problem.xu.shape == (n,)
        assert problem.aub.shape == (1, n)
        assert problem.bub.shape == (1,)
        assert problem.aeq.shape == (1, n)
        assert problem.beq.shape == (1,)
        assert problem.fun_hist.shape == (0,)
        assert problem.maxcv_hist.shape == (0,)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.xl, np.full(n, -np.inf))
        np.testing.assert_array_equal(problem.xu, np.full(n, np.inf))

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)
        assert problem.n_fun_eval == 1
        np.testing.assert_array_equal(problem.fun_hist, [self.rosen(x0)])
        np.testing.assert_array_equal(problem.maxcv_hist, [1.0])

        # Perform a second objective function evaluation.
        assert problem.fun(problem.x0 + 1.0) == self.rosen(x0 + 1.0)
        assert problem.n_fun_eval == 2
        np.testing.assert_array_equal(problem.fun_hist, [self.rosen(x0), self.rosen(x0 + 1.0)])
        np.testing.assert_array_almost_equal(problem.maxcv_hist, [1.0, abs(n * (n - 1) / 2 - 1.0)])

        # Evaluate the nonlinear constraints.
        assert problem.cub(problem.x0).shape == (0,)
        assert problem.ceq(problem.x0).shape == (0,)
        assert problem.n_fun_eval == 2
