import numpy as np
import pytest

from OptiProfiler import Problem


class TestProblem:

    @staticmethod
    def rosen(x):
        return np.sum(100 * (x[1:] - x[:-1] ** 2) ** 2 + (1 - x[:-1]) ** 2)

    @staticmethod
    def sum_cos(x):
        return np.sum(np.cos(x))

    @staticmethod
    def sum_sin(x):
        return np.sum(np.sin(x))

    @staticmethod
    def assert_dimensions(problem, n, m_linear_ub, m_linear_eq, m_nonlinear_ub, m_nonlinear_eq):
        assert problem.n == n
        assert problem.m_linear_ub == m_linear_ub
        assert problem.m_linear_eq == m_linear_eq
        assert problem.m_nonlinear_ub == m_nonlinear_ub
        assert problem.m_nonlinear_eq == m_nonlinear_eq
        assert problem.xl.shape == (n,)
        assert problem.xu.shape == (n,)
        assert problem.aub.shape == (m_linear_ub, n)
        assert problem.bub.shape == (m_linear_ub,)
        assert problem.aeq.shape == (m_linear_eq, n)
        assert problem.beq.shape == (m_linear_eq,)

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_unconstrained_problem(self, n):
        # Construct an unconstrained problem.
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)

        # Check the problem attributes.
        self.assert_dimensions(problem, n, 0, 0, 0, 0)
        assert problem.n_fun_eval == 0
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
        self.assert_dimensions(problem, n, 0, 0, 0, 0)
        assert problem.n_fun_eval == 0
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
        self.assert_dimensions(problem, n, 1, 1, 0, 0)
        assert problem.n_fun_eval == 0
        assert problem.fun_hist.shape == (0,)
        assert problem.maxcv_hist.shape == (0,)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.xl, np.full(n, -np.inf))
        np.testing.assert_array_equal(problem.xu, np.full(n, np.inf))
        np.testing.assert_array_equal(problem.aub, aub)
        np.testing.assert_array_equal(problem.bub, bub)
        np.testing.assert_array_equal(problem.aeq, aeq.reshape(1, n))
        np.testing.assert_array_equal(problem.beq, beq)

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)
        assert problem.n_fun_eval == 1
        np.testing.assert_array_equal(problem.fun_hist, [self.rosen(x0)])
        np.testing.assert_array_almost_equal(problem.maxcv_hist, [1.0])

        # Perform a second objective function evaluation.
        assert problem.fun(problem.x0 + 1.0) == self.rosen(x0 + 1.0)
        assert problem.n_fun_eval == 2
        np.testing.assert_array_equal(problem.fun_hist, [self.rosen(x0), self.rosen(x0 + 1.0)])
        np.testing.assert_array_almost_equal(problem.maxcv_hist, [1.0, abs(n * (n - 1) / 2 - 1.0)])

        # Evaluate the nonlinear constraints.
        assert problem.cub(problem.x0).shape == (0,)
        assert problem.ceq(problem.x0).shape == (0,)
        assert problem.n_fun_eval == 2

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_nonlinearly_constrained_problem(self, n):
        # Construct a nonlinearly constrained problem.
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0, cub=self.sum_cos, ceq=self.sum_sin, m_nonlinear_ub=1, m_nonlinear_eq=1)

        # Check the problem attributes.
        self.assert_dimensions(problem, n, 0, 0, 1, 1)
        assert problem.n_fun_eval == 0
        assert problem.fun_hist.shape == (0,)
        assert problem.maxcv_hist.shape == (0,)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.xl, np.full(n, -np.inf))
        np.testing.assert_array_equal(problem.xu, np.full(n, np.inf))

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)
        assert problem.n_fun_eval == 1
        np.testing.assert_array_equal(problem.fun_hist, [self.rosen(x0)])
        np.testing.assert_array_almost_equal(problem.maxcv_hist, [n])

        # Perform a second objective function evaluation.
        assert problem.fun(problem.x0 + np.pi / 2) == self.rosen(x0 + np.pi / 2)
        assert problem.n_fun_eval == 2
        np.testing.assert_array_equal(problem.fun_hist, [self.rosen(x0), self.rosen(x0 + np.pi / 2)])
        np.testing.assert_array_almost_equal(problem.maxcv_hist, [n, n])

    def test_exceptions(self):
        with pytest.raises(ValueError):
            Problem('fun', np.zeros(1))
        with pytest.raises(ValueError):
            Problem(lambda x, y: x + y, np.zeros(2))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), cub='cub')
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), cub=lambda x, y: x + y)
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), ceq='ceq')
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), ceq=lambda x, y: x + y)
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), cub=self.sum_cos, m_nonlinear_ub=1.5)
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), ceq=self.sum_sin, m_nonlinear_eq=1.5)
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), np.zeros(2))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), xu=np.zeros(2))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), aub=np.zeros((1, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), aeq=np.zeros((1, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), m_nonlinear_ub=1)
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), m_nonlinear_eq=1)
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1), cub=self.sum_cos)
            problem.m_nonlinear_ub
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1), ceq=self.sum_sin)
            problem.m_nonlinear_eq
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1))
            problem.fun(np.zeros(2))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1), cub=self.sum_cos)
            problem.cub(np.zeros(2))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1), cub=lambda x: np.zeros(int(np.sum(x))))
            problem.cub(np.zeros(1))
            problem.cub(np.ones(1))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1), ceq=self.sum_sin)
            problem.ceq(np.zeros(2))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1), ceq=lambda x: np.zeros(int(np.sum(x))))
            problem.ceq(np.zeros(1))
            problem.ceq(np.ones(1))
