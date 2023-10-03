import numpy as np
import pytest

from OptiProfiler.features import Feature
from OptiProfiler.problems import FeaturedProblem, Problem, ProblemError, load_cutest, find_cutest_problem_names


class BaseTestProblem:

    @staticmethod
    def rosen(x):
        return np.sum(1e2 * (x[1:] - x[:-1] ** 2) ** 2 + (1.0 - x[:-1]) ** 2)

    @staticmethod
    def sum_cos(x):
        return np.sum(np.cos(x))

    @staticmethod
    def sum_sin(x):
        return np.sum(np.sin(x))

    @staticmethod
    def bad_fun(x):
        raise RuntimeError

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


class TestProblem(BaseTestProblem):

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_unconstrained_problem(self, n):
        # Construct an unconstrained problem.
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)

        # Check the problem attributes.
        self.assert_dimensions(problem, n, 0, 0, 0, 0)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.xl, np.full(n, -np.inf))
        np.testing.assert_array_equal(problem.xu, np.full(n, np.inf))

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)
        assert problem.maxcv(problem.x0) == 0.0

        # Evaluate the nonlinear constraints.
        assert problem.cub(problem.x0).shape == (0,)
        assert problem.ceq(problem.x0).shape == (0,)

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_bound_constrained_problem(self, n):
        # Construct a bound-constrained problem.
        x0 = np.zeros(n)
        xl = np.full(n, -2.048)
        xu = np.full(n, 2.048)
        problem = Problem(self.rosen, x0, xl, xu)

        # Check the problem attributes.
        self.assert_dimensions(problem, n, 0, 0, 0, 0)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.xl, xl)
        np.testing.assert_array_equal(problem.xu, xu)

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)
        assert problem.maxcv(problem.x0) == 0.0

        # Evaluate the nonlinear constraints.
        assert problem.cub(problem.x0).shape == (0,)
        assert problem.ceq(problem.x0).shape == (0,)

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
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.xl, np.full(n, -np.inf))
        np.testing.assert_array_equal(problem.xu, np.full(n, np.inf))
        np.testing.assert_array_equal(problem.aub, aub)
        np.testing.assert_array_equal(problem.bub, bub)
        np.testing.assert_array_equal(problem.aeq, aeq.reshape(1, n))
        np.testing.assert_array_equal(problem.beq, beq)

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)
        np.testing.assert_allclose(problem.maxcv(problem.x0), 1.0)

        # Evaluate the nonlinear constraints.
        assert problem.cub(problem.x0).shape == (0,)
        assert problem.ceq(problem.x0).shape == (0,)

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_nonlinearly_constrained_problem(self, n):
        # Construct a nonlinearly constrained problem.
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0, cub=self.sum_cos, ceq=self.sum_sin, m_nonlinear_ub=1, m_nonlinear_eq=1)

        # Check the problem attributes.
        self.assert_dimensions(problem, n, 0, 0, 1, 1)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.xl, np.full(n, -np.inf))
        np.testing.assert_array_equal(problem.xu, np.full(n, np.inf))

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)
        np.testing.assert_allclose(problem.maxcv(problem.x0), n)

    def test_exceptions(self):
        with pytest.raises(TypeError):
            Problem('fun', np.zeros(1))
        with pytest.raises(TypeError):
            Problem(self.rosen, np.zeros(2), cub='cub')
        with pytest.raises(TypeError):
            Problem(self.rosen, np.zeros(2), ceq='ceq')
        with pytest.raises(TypeError):
            Problem(self.rosen, np.zeros(1), cub=self.sum_cos, m_nonlinear_ub=1.5)
        with pytest.raises(TypeError):
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
            problem.m_nonlinear_ub  # noqa
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1), ceq=self.sum_sin)
            problem.m_nonlinear_eq  # noqa
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
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1))
            problem.maxcv(np.zeros(2))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros((2, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), np.zeros((2, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), xu=np.zeros((2, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), aub=np.zeros((2, 2)), bub=np.zeros((2, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), aeq=np.zeros((2, 2)), beq=np.zeros((2, 2)))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(2))
            problem.fun(np.zeros((2, 2)))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(2))
            problem.cub(np.zeros((2, 2)))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(2), cub=lambda _: np.zeros((2, 2)))
            problem.cub(problem.x0)
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(2))
            problem.ceq(np.zeros((2, 2)))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(2), ceq=lambda _: np.zeros((2, 2)))
            problem.ceq(problem.x0)
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(2))
            problem.maxcv(np.zeros((2, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), aub=np.zeros((2, 2, 2)), bub=np.zeros(2))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), aeq=np.zeros((2, 2, 2)), beq=np.zeros(2))

    def test_catch(self):
        # The value m_nonlinear_ub can be a float.
        problem = Problem(self.rosen, np.zeros(1), cub=self.sum_cos, m_nonlinear_ub=1.0)
        assert problem.m_nonlinear_ub == 1
        assert problem.m_nonlinear_eq == 0

        # The value m_nonlinear_eq can be a float.
        problem = Problem(self.rosen, np.zeros(1), ceq=self.sum_sin, m_nonlinear_eq=1.0)
        assert problem.m_nonlinear_ub == 0
        assert problem.m_nonlinear_eq == 1

        # The objective function can be ill-defined.
        problem = Problem(self.bad_fun, np.zeros(1))
        with pytest.warns(RuntimeWarning):
            assert np.isnan(problem.fun(problem.x0))
        assert problem.maxcv(problem.x0) == 0.0

        # The nonlinear inequality constraint function can be ill-defined.
        problem = Problem(self.rosen, np.zeros(1), cub=self.bad_fun, m_nonlinear_ub=1)
        with pytest.warns(RuntimeWarning):
            assert np.isnan(problem.cub(problem.x0))

        # The nonlinear equality constraint function can be ill-defined.
        problem = Problem(self.rosen, np.zeros(1), ceq=self.bad_fun, m_nonlinear_eq=1)
        with pytest.warns(RuntimeWarning):
            assert np.isnan(problem.ceq(problem.x0))


class TestFeaturedProblem(BaseTestProblem):

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_simple(self, n):
        # Construct a simple problem.
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)

        # Construct a featured problem.
        feature = Feature('plain')
        featured_problem = FeaturedProblem(problem, feature)

        # Check the featured problem attributes.
        assert featured_problem.n_eval == 0
        assert featured_problem.fun_values.shape == (0,)
        assert featured_problem.maxcv_values.shape == (0,)

        # Evaluate the objective function at x0.
        f = featured_problem.fun(x0)
        assert featured_problem.n_eval == 1
        np.testing.assert_array_equal(featured_problem.fun_values, [f])
        np.testing.assert_array_equal(featured_problem.maxcv_values, [0.0])

        # Construct a featured problem with a different feature.
        feature = Feature('custom', modifier=lambda x, f, seed: f + 1)
        featured_problem = FeaturedProblem(problem, feature)

        # Evaluate the objective function at x0.
        f = featured_problem.fun(x0)
        assert featured_problem.n_eval == 1
        np.testing.assert_allclose(featured_problem.fun_values, [f - 1])
        np.testing.assert_array_equal(featured_problem.maxcv_values, [0.0])

    def test_catch(self):
        # Construct a nonlinearly constrained problem.
        x0 = np.zeros(2)
        problem = Problem(self.rosen, x0, cub=self.sum_cos, ceq=self.sum_sin)

        # Construct a featured problem.
        feature = Feature('plain')
        FeaturedProblem(problem, feature)


@pytest.mark.extra
class TestLoadCUTEst:

    @pytest.mark.parametrize('constraint', ['unconstrained', 'fixed', 'bound', 'adjacency', 'linear', 'quadratic', 'other'])
    def test_simple(self, constraint):
        problem_names = find_cutest_problem_names(constraint, n_min=1, n_max=10, m_min=0, m_max=100)
        for i_problem in range(min(10, len(problem_names))):
            try:
                problem = load_cutest(problem_names[i_problem], n_min=1, n_max=10, m_min=0, m_max=100)
                assert isinstance(problem, Problem)
                problem.fun(problem.x0)
                problem.cub(problem.x0)
                problem.ceq(problem.x0)
            except ProblemError:
                pass

    def test_exceptions(self):
        with pytest.raises(TypeError):
            load_cutest(1)
        with pytest.raises(ValueError):
            load_cutest('ARGLINA', N=1)
        with pytest.raises(TypeError):
            load_cutest('ARGLINA', n_min=1.5)
        with pytest.raises(ValueError):
            load_cutest('ARGLINA', n_min=5, n_max=1)
        with pytest.raises(ValueError):
            load_cutest('ARGLINA', m_min=5, m_max=1)

    def test_catch(self):
        problem_names = find_cutest_problem_names('unconstrained', n_min=1, n_max=10)
        load_cutest(problem_names[0], n_min=1.0, n_max=10.0)


@pytest.mark.extra
class TestFindCUTEstProblemNames:

    def test_simple(self):
        for constraint in ['unconstrained', 'fixed', 'bound', 'adjacency', 'linear', 'quadratic', 'other']:
            problem_names = find_cutest_problem_names(constraint)
            assert isinstance(problem_names, list)
            for problem_name in problem_names:
                assert isinstance(problem_name, str)

        # Check that constraint types can be combined.
        unconstrained_problem_names = find_cutest_problem_names('unconstrained')
        bound_problem_names = find_cutest_problem_names('bound')
        both_problem_names = find_cutest_problem_names('unconstrained bound')
        assert set(both_problem_names) == set(unconstrained_problem_names + bound_problem_names)

    def test_parameters(self):
        # Check the n_min and n_max parameters.
        problem_names_1_100 = find_cutest_problem_names('unconstrained', n_min=1, n_max=100)
        problem_names_10_50 = find_cutest_problem_names('unconstrained', n_min=10, n_max=50)
        assert set(problem_names_10_50).issubset(set(problem_names_1_100))

        # Check the m_min and m_max parameters.
        problem_names_1_100 = find_cutest_problem_names('adjacency linear quadratic other', m_min=1, m_max=100)
        problem_names_10_50 = find_cutest_problem_names('adjacency linear quadratic other', m_min=10, m_max=50)
        assert set(problem_names_10_50).issubset(set(problem_names_1_100))

    def test_exceptions(self):
        with pytest.raises(TypeError):
            find_cutest_problem_names(1)
        with pytest.raises(ValueError):
            find_cutest_problem_names('cubic')
        with pytest.raises(TypeError):
            find_cutest_problem_names('unconstrained', n_min=0)
        with pytest.raises(TypeError):
            find_cutest_problem_names('unconstrained', n_min=1.5)
        with pytest.raises(TypeError):
            find_cutest_problem_names('unconstrained', n_max=0)
        with pytest.raises(TypeError):
            find_cutest_problem_names('unconstrained', n_max=1.5)
        with pytest.raises(TypeError):
            find_cutest_problem_names('unconstrained', n_min=3, n_max=2)
        with pytest.raises(TypeError):
            find_cutest_problem_names('unconstrained', m_min=-1)
        with pytest.raises(TypeError):
            find_cutest_problem_names('unconstrained', m_min=1.5)
        with pytest.raises(TypeError):
            find_cutest_problem_names('unconstrained', m_max=-1)
        with pytest.raises(TypeError):
            find_cutest_problem_names('unconstrained', m_max=1.5)
        with pytest.raises(TypeError):
            find_cutest_problem_names('unconstrained', m_min=3, m_max=2)

    def test_catch(self):
        find_cutest_problem_names('quadratic', n_min=1.0, n_max=10.0, m_min=1.0, m_max=100.0)
