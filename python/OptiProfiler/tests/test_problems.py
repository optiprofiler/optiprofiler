import numpy as np
import pytest

from OptiProfiler.features import Feature
from OptiProfiler.problems import FeaturedProblem, Problem, ProblemError, load_cutest_problem, find_cutest_problems


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
    def assert_dimensions(problem, n, num_linear_ub, num_linear_eq, num_nonlinear_ub, num_nonlinear_eq):
        assert problem.dimension == n
        assert problem.num_linear_ub == num_linear_ub
        assert problem.num_linear_eq == num_linear_eq
        assert problem.num_nonlinear_ub == num_nonlinear_ub
        assert problem.num_nonlinear_eq == num_nonlinear_eq
        assert problem.lb.shape == (n,)
        assert problem.ub.shape == (n,)
        assert problem.a_ub.shape == (num_linear_ub, n)
        assert problem.b_ub.shape == (num_linear_ub,)
        assert problem.a_eq.shape == (num_linear_eq, n)
        assert problem.b_eq.shape == (num_linear_eq,)


class TestProblem(BaseTestProblem):

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_unconstrained_problem(self, n):
        # Construct an unconstrained problem.
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)

        # Check the problem attributes.
        self.assert_dimensions(problem, n, 0, 0, 0, 0)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.lb, np.full(n, -np.inf))
        np.testing.assert_array_equal(problem.ub, np.full(n, np.inf))

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)
        assert problem.maxcv(problem.x0) == 0.0

        # Evaluate the nonlinear constraints.
        assert problem.c_ub(problem.x0).shape == (0,)
        assert problem.c_eq(problem.x0).shape == (0,)

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_bound_constrained_problem(self, n):
        # Construct a bound-constrained problem.
        x0 = np.zeros(n)
        lb = np.full(n, -2.048)
        ub = np.full(n, 2.048)
        problem = Problem(self.rosen, x0, lb, ub)

        # Check the problem attributes.
        self.assert_dimensions(problem, n, 0, 0, 0, 0)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.lb, lb)
        np.testing.assert_array_equal(problem.ub, ub)

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)
        assert problem.maxcv(problem.x0) == 0.0

        # Evaluate the nonlinear constraints.
        assert problem.c_ub(problem.x0).shape == (0,)
        assert problem.c_eq(problem.x0).shape == (0,)

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_linearly_constrained_problem(self, n):
        # Construct a linearly constrained problem.
        x0 = np.zeros(n)
        a_ub = np.ones((1, n))
        b_ub = 1.0
        a_eq = np.arange(n)
        b_eq = 1.0
        problem = Problem(self.rosen, x0, a_ub=a_ub, b_ub=b_ub, a_eq=a_eq, b_eq=b_eq)

        # Check the problem attributes.
        self.assert_dimensions(problem, n, 1, 1, 0, 0)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.lb, np.full(n, -np.inf))
        np.testing.assert_array_equal(problem.ub, np.full(n, np.inf))
        np.testing.assert_array_equal(problem.a_ub, a_ub)
        np.testing.assert_array_equal(problem.b_ub, b_ub)
        np.testing.assert_array_equal(problem.a_eq, a_eq.reshape(1, n))
        np.testing.assert_array_equal(problem.b_eq, b_eq)

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)
        np.testing.assert_allclose(problem.maxcv(problem.x0), 1.0)

        # Evaluate the nonlinear constraints.
        assert problem.c_ub(problem.x0).shape == (0,)
        assert problem.c_eq(problem.x0).shape == (0,)

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_nonlinearly_constrained_problem(self, n):
        # Construct a nonlinearly constrained problem.
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0, c_ub=self.sum_cos, c_eq=self.sum_sin, num_nonlinear_ub=1, num_nonlinear_eq=1)

        # Check the problem attributes.
        self.assert_dimensions(problem, n, 0, 0, 1, 1)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.lb, np.full(n, -np.inf))
        np.testing.assert_array_equal(problem.ub, np.full(n, np.inf))

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)
        np.testing.assert_allclose(problem.maxcv(problem.x0), n)

    def test_exceptions(self):
        with pytest.raises(TypeError):
            Problem('fun', np.zeros(1))
        with pytest.raises(TypeError):
            Problem(self.rosen, np.zeros(2), c_ub='c_ub')
        with pytest.raises(TypeError):
            Problem(self.rosen, np.zeros(2), c_eq='c_eq')
        with pytest.raises(TypeError):
            Problem(self.rosen, np.zeros(1), c_ub=self.sum_cos, num_nonlinear_ub=1.5)
        with pytest.raises(TypeError):
            Problem(self.rosen, np.zeros(1), c_eq=self.sum_sin, num_nonlinear_eq=1.5)
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), np.zeros(2))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), ub=np.zeros(2))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), a_ub=np.zeros((1, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), a_eq=np.zeros((1, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), num_nonlinear_ub=1)
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(1), num_nonlinear_eq=1)
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1), c_ub=self.sum_cos)
            problem.num_nonlinear_ub  # noqa
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1), c_eq=self.sum_sin)
            problem.num_nonlinear_eq  # noqa
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1))
            problem.fun(np.zeros(2))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1), c_ub=self.sum_cos)
            problem.c_ub(np.zeros(2))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1), c_ub=lambda x: np.zeros(int(np.sum(x))))
            problem.c_ub(np.zeros(1))
            problem.c_ub(np.ones(1))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1), c_eq=self.sum_sin)
            problem.c_eq(np.zeros(2))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1), c_eq=lambda x: np.zeros(int(np.sum(x))))
            problem.c_eq(np.zeros(1))
            problem.c_eq(np.ones(1))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(1))
            problem.maxcv(np.zeros(2))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros((2, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), np.zeros((2, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), ub=np.zeros((2, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), a_ub=np.zeros((2, 2)), b_ub=np.zeros((2, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), a_eq=np.zeros((2, 2)), b_eq=np.zeros((2, 2)))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(2))
            problem.fun(np.zeros((2, 2)))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(2))
            problem.c_ub(np.zeros((2, 2)))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(2), c_ub=lambda _: np.zeros((2, 2)))
            problem.c_ub(problem.x0)
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(2))
            problem.c_eq(np.zeros((2, 2)))
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(2), c_eq=lambda _: np.zeros((2, 2)))
            problem.c_eq(problem.x0)
        with pytest.raises(ValueError):
            problem = Problem(self.rosen, np.zeros(2))
            problem.maxcv(np.zeros((2, 2)))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), a_ub=np.zeros((2, 2, 2)), b_ub=np.zeros(2))
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), a_eq=np.zeros((2, 2, 2)), b_eq=np.zeros(2))

    def test_catch(self):
        # The value num_nonlinear_ub can be a float.
        problem = Problem(self.rosen, np.zeros(1), c_ub=self.sum_cos, num_nonlinear_ub=1.0)
        assert problem.num_nonlinear_ub == 1
        assert problem.num_nonlinear_eq == 0

        # The value num_nonlinear_eq can be a float.
        problem = Problem(self.rosen, np.zeros(1), c_eq=self.sum_sin, num_nonlinear_eq=1.0)
        assert problem.num_nonlinear_ub == 0
        assert problem.num_nonlinear_eq == 1

        # The objective function can be ill-defined.
        problem = Problem(self.bad_fun, np.zeros(1))
        assert np.isnan(problem.fun(problem.x0))
        assert problem.maxcv(problem.x0) == 0.0

        # The nonlinear inequality constraint function can be ill-defined.
        problem = Problem(self.rosen, np.zeros(1), c_ub=self.bad_fun, num_nonlinear_ub=1)
        assert np.isnan(problem.c_ub(problem.x0))

        # The nonlinear equality constraint function can be ill-defined.
        problem = Problem(self.rosen, np.zeros(1), c_eq=self.bad_fun, num_nonlinear_eq=1)
        assert np.isnan(problem.c_eq(problem.x0))


class TestFeaturedProblem(BaseTestProblem):

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_simple(self, n):
        # Construct a simple problem.
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)

        # Construct a featured problem.
        feature = Feature('plain')
        featured_problem = FeaturedProblem(problem, feature, 500 * n)

        # Check the featured problem attributes.
        assert featured_problem.n_eval == 0
        assert featured_problem.fun_history.shape == (0,)
        assert featured_problem.maxcv_history.shape == (0,)

        # Evaluate the objective function at x0.
        f = featured_problem.fun(x0)
        assert featured_problem.n_eval == 1
        np.testing.assert_array_equal(featured_problem.fun_history, [f])
        np.testing.assert_array_equal(featured_problem.maxcv_history, [0.0])

        # Construct a featured problem with a different feature.
        feature = Feature('custom', modifier=lambda x, f, seed: f + 1)
        featured_problem = FeaturedProblem(problem, feature, 500 * n)

        # Evaluate the objective function at x0.
        f = featured_problem.fun(x0)
        assert featured_problem.n_eval == 1
        np.testing.assert_allclose(featured_problem.fun_history, [f - 1])
        np.testing.assert_array_equal(featured_problem.maxcv_history, [0.0])

    def test_catch(self):
        # Construct a nonlinearly constrained problem.
        x0 = np.zeros(2)
        problem = Problem(self.rosen, x0, c_ub=self.sum_cos, c_eq=self.sum_sin)

        # Construct a featured problem.
        feature = Feature('plain')
        FeaturedProblem(problem, feature, 1000)

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_randomize_x0(self, n):
        # Construct a simple problem.
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)

        # Construct a featured problem.
        feature = Feature('randomize_x0', distribution=lambda rng, n: np.ones(n))
        featured_problem = FeaturedProblem(problem, feature, 500 * n)

        # Evaluate the objective function at x0.
        f = featured_problem.fun(featured_problem.x0)
        np.testing.assert_allclose(f, problem.fun(x0 + 1.0))

    def test_max_eval(self):
        # Construct a simple problem.
        x0 = np.zeros(2)
        problem = Problem(self.rosen, x0)

        # Construct a featured problem.
        feature = Feature('plain')
        featured_problem = FeaturedProblem(problem, feature, 1)

        # Evaluate the objective function at x0 twice.
        featured_problem.fun(featured_problem.x0)
        with pytest.raises(RuntimeError):
            featured_problem.fun(featured_problem.x0)


@pytest.mark.extra
class TestLoadCUTEst:

    @pytest.mark.parametrize('constraint', ['unconstrained', 'fixed', 'bound', 'adjacency', 'linear', 'quadratic', 'other'])
    def test_simple(self, constraint):
        problem_names = find_cutest_problems(constraint, n_min=1, n_max=10, m_min=0, m_max=100)
        for i_problem in range(min(10, len(problem_names))):
            try:
                problem = load_cutest_problem(problem_names[i_problem], n_min=1, n_max=10, m_min=0, m_max=100)
                assert isinstance(problem, Problem)
                problem.fun(problem.x0)
                problem.c_ub(problem.x0)
                problem.c_eq(problem.x0)
            except ProblemError:
                pass

    def test_exceptions(self):
        with pytest.raises(TypeError):
            load_cutest_problem(1)
        with pytest.raises(ValueError):
            load_cutest_problem('ARGLINA', N=1)
        with pytest.raises(TypeError):
            load_cutest_problem('ARGLINA', n_min=1.5)
        with pytest.raises(ValueError):
            load_cutest_problem('ARGLINA', n_min=5, n_max=1)
        with pytest.raises(ValueError):
            load_cutest_problem('ARGLINA', m_min=5, m_max=1)

    def test_catch(self):
        problem_names = find_cutest_problems('unconstrained', n_min=1, n_max=10)
        load_cutest_problem(problem_names[0], n_min=1.0, n_max=10.0)


@pytest.mark.extra
class TestFindCUTEstProblemNames:

    def test_simple(self):
        for constraint in ['unconstrained', 'fixed', 'bound', 'adjacency', 'linear', 'quadratic', 'other']:
            problem_names = find_cutest_problems(constraint)
            assert isinstance(problem_names, list)
            for problem_name in problem_names:
                assert isinstance(problem_name, str)

        # Check that constraint types can be combined.
        unconstrained_problem_names = find_cutest_problems('unconstrained')
        bound_problem_names = find_cutest_problems('bound')
        both_problem_names = find_cutest_problems('unconstrained bound')
        assert set(both_problem_names) == set(unconstrained_problem_names + bound_problem_names)

    def test_parameters(self):
        # Check the n_min and n_max parameters.
        problem_names_1_100 = find_cutest_problems('unconstrained', n_min=1, n_max=100)
        problem_names_10_50 = find_cutest_problems('unconstrained', n_min=10, n_max=50)
        assert set(problem_names_10_50).issubset(set(problem_names_1_100))

        # Check the m_min and m_max parameters.
        problem_names_1_100 = find_cutest_problems('adjacency linear quadratic other', m_min=1, m_max=100)
        problem_names_10_50 = find_cutest_problems('adjacency linear quadratic other', m_min=10, m_max=50)
        assert set(problem_names_10_50).issubset(set(problem_names_1_100))

    def test_exceptions(self):
        with pytest.raises(TypeError):
            find_cutest_problems(1)
        with pytest.raises(ValueError):
            find_cutest_problems('cubic')
        with pytest.raises(TypeError):
            find_cutest_problems('unconstrained', n_min=0)
        with pytest.raises(TypeError):
            find_cutest_problems('unconstrained', n_min=1.5)
        with pytest.raises(TypeError):
            find_cutest_problems('unconstrained', n_max=0)
        with pytest.raises(TypeError):
            find_cutest_problems('unconstrained', n_max=1.5)
        with pytest.raises(TypeError):
            find_cutest_problems('unconstrained', n_min=3, n_max=2)
        with pytest.raises(TypeError):
            find_cutest_problems('unconstrained', m_min=-1)
        with pytest.raises(TypeError):
            find_cutest_problems('unconstrained', m_min=1.5)
        with pytest.raises(TypeError):
            find_cutest_problems('unconstrained', m_max=-1)
        with pytest.raises(TypeError):
            find_cutest_problems('unconstrained', m_max=1.5)
        with pytest.raises(TypeError):
            find_cutest_problems('unconstrained', m_min=3, m_max=2)

    def test_catch(self):
        find_cutest_problems('quadratic', n_min=1.0, n_max=10.0, m_min=1.0, m_max=100.0)
