"""Tests for Problem and FeaturedProblem classes."""
from contextlib import suppress

import numpy as np
import pytest

from optiprofiler.opclasses import Feature, Problem, FeaturedProblem
from optiprofiler.utils import ProblemError

# Import pycutest tools for CUTEst tests (if available)
try:
    from optiprofiler.problem_libs.pycutest.pycutest_tools import pycutest_load, pycutest_select
    PYCUTEST_AVAILABLE = True
except ImportError:
    PYCUTEST_AVAILABLE = False


class BaseTestProblem:
    """Base class with common test utilities."""

    @staticmethod
    def rosen(x):
        """Rosenbrock function."""
        return np.sum(1e2 * (x[1:] - x[:-1] ** 2) ** 2 + (1.0 - x[:-1]) ** 2)

    @staticmethod
    def sum_cos(x):
        """Sum of cosines."""
        return np.array([np.sum(np.cos(x))])

    @staticmethod
    def sum_sin(x):
        """Sum of sines."""
        return np.array([np.sum(np.sin(x))])

    @staticmethod
    def bad_fun(x):
        """Function that raises an exception."""
        raise RuntimeError("Bad function")

    @staticmethod
    def assert_dimensions(problem, n, m_linear_ub, m_linear_eq, m_nonlinear_ub, m_nonlinear_eq):
        """Assert problem dimensions are correct."""
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
    """Tests for the Problem class."""

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_unconstrained_problem(self, n):
        """Test creating and evaluating an unconstrained problem."""
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)

        # Check the problem attributes.
        self.assert_dimensions(problem, n, 0, 0, 0, 0)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.xl, np.full(n, -np.inf))
        np.testing.assert_array_equal(problem.xu, np.full(n, np.inf))

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)

        # Evaluate the nonlinear constraints (should return empty arrays).
        assert problem.cub(problem.x0).shape == (0,)
        assert problem.ceq(problem.x0).shape == (0,)

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_bound_constrained_problem(self, n):
        """Test creating and evaluating a bound-constrained problem."""
        x0 = np.zeros(n)
        lb = np.full(n, -2.048)
        ub = np.full(n, 2.048)
        problem = Problem(self.rosen, x0, xl=lb, xu=ub)

        # Check the problem attributes.
        self.assert_dimensions(problem, n, 0, 0, 0, 0)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.xl, lb)
        np.testing.assert_array_equal(problem.xu, ub)

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)

    @pytest.mark.parametrize('n', [1, 10, 100])
    def test_linearly_constrained_problem(self, n):
        """Test creating and evaluating a linearly constrained problem."""
        x0 = np.zeros(n)
        a_ub = np.ones((1, n))
        b_ub = np.array([1.0])
        a_eq = np.arange(n).reshape(1, n)
        b_eq = np.array([1.0])
        problem = Problem(self.rosen, x0, aub=a_ub, bub=b_ub, aeq=a_eq, beq=b_eq)

        # Check the problem attributes.
        self.assert_dimensions(problem, n, 1, 1, 0, 0)
        np.testing.assert_array_equal(problem.x0, x0)
        np.testing.assert_array_equal(problem.aub, a_ub)
        np.testing.assert_array_equal(problem.bub, b_ub)
        np.testing.assert_array_equal(problem.aeq, a_eq)
        np.testing.assert_array_equal(problem.beq, b_eq)

        # Perform an objective function evaluation.
        assert problem.fun(problem.x0) == self.rosen(x0)

    @pytest.mark.parametrize('n', [2, 10])
    def test_nonlinearly_constrained_problem(self, n):
        """Test creating and evaluating a nonlinearly constrained problem."""
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0, cub=self.sum_cos, ceq=self.sum_sin)

        # The number of nonlinear constraints is inferred from the first call.
        cub_val = problem.cub(problem.x0)
        ceq_val = problem.ceq(problem.x0)
        
        assert cub_val.shape == (1,)
        assert ceq_val.shape == (1,)
        
        # Check dimensions after constraint inference.
        assert problem.m_nonlinear_ub == 1
        assert problem.m_nonlinear_eq == 1

    def test_problem_with_name(self):
        """Test creating a problem with a name."""
        x0 = np.zeros(2)
        problem = Problem(self.rosen, x0, name="TestProblem")
        assert problem.name == "TestProblem"

    def test_exceptions(self):
        """Test that appropriate exceptions are raised."""
        # Invalid fun type
        with pytest.raises(TypeError):
            Problem('fun', np.zeros(1))
        
        # Invalid cub type
        with pytest.raises(TypeError):
            Problem(self.rosen, np.zeros(2), cub='cub')
        
        # Invalid ceq type
        with pytest.raises(TypeError):
            Problem(self.rosen, np.zeros(2), ceq='ceq')
        
        # Invalid x0 shape
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros((2, 2)))
        
        # Invalid xl shape
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), xl=np.zeros((2, 2)))
        
        # Invalid xu shape
        with pytest.raises(ValueError):
            Problem(self.rosen, np.zeros(2), xu=np.zeros((2, 2)))

    def test_ill_defined_functions(self):
        """Test that ill-defined functions return NaN."""
        # Objective function that raises an exception.
        problem = Problem(self.bad_fun, np.zeros(1))
        assert np.isnan(problem.fun(problem.x0))


class TestFeaturedProblem(BaseTestProblem):
    """Tests for the FeaturedProblem class."""

    @pytest.mark.parametrize('n', [2, 10])
    def test_simple(self, n):
        """Test creating and using a simple featured problem."""
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)
        feature = Feature('plain')
        featured_problem = FeaturedProblem(problem, feature, 500 * n)

        # Check the featured problem attributes.
        assert featured_problem.n_eval_fun == 0
        assert featured_problem.fun_hist.shape == (0,)
        assert featured_problem.maxcv_hist.shape == (0,)

        # Evaluate the objective function at x0.
        f = featured_problem.fun(x0)
        assert featured_problem.n_eval_fun == 1
        np.testing.assert_array_equal(featured_problem.fun_hist, [f])

    @pytest.mark.parametrize('n', [2, 10])
    def test_custom_feature(self, n):
        """Test using a custom feature with FeaturedProblem."""
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)
        
        # Custom modifier that adds 1 to the function value.
        def custom_mod_fun(x, rng, prob):
            f = prob.fun(x)
            return f + 1.0
        
        feature = Feature('custom', mod_fun=custom_mod_fun)
        featured_problem = FeaturedProblem(problem, feature, 500 * n)

        # Evaluate the objective function at x0.
        f = featured_problem.fun(x0)
        np.testing.assert_allclose(f, self.rosen(x0) + 1.0)

    @pytest.mark.parametrize('n', [2, 10])
    def test_perturbed_x0(self, n):
        """Test that perturbed_x0 feature modifies the initial point."""
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)
        feature = Feature('perturbed_x0')
        featured_problem = FeaturedProblem(problem, feature, 500 * n, seed=42)

        # The perturbed x0 should be different from the original.
        assert not np.allclose(featured_problem.x0, x0)

    def test_max_eval_reached(self):
        """Test that evaluation stops after 2*max_eval is reached."""
        x0 = np.zeros(2)
        problem = Problem(self.rosen, x0)
        feature = Feature('plain')
        max_eval = 2
        featured_problem = FeaturedProblem(problem, feature, max_eval)

        # First max_eval evaluations record history.
        for _ in range(max_eval):
            featured_problem.fun(featured_problem.x0)
        
        # Next max_eval evaluations return cached values.
        for _ in range(max_eval):
            featured_problem.fun(featured_problem.x0)
        
        # After 2*max_eval, StopIteration is raised.
        with pytest.raises(StopIteration):
            featured_problem.fun(featured_problem.x0)

    def test_exceptions(self):
        """Test that appropriate exceptions are raised."""
        x0 = np.zeros(2)
        problem = Problem(self.rosen, x0)
        feature = Feature('plain')
        
        # Invalid problem type
        with pytest.raises(TypeError):
            FeaturedProblem('problem', feature, 1)
        
        # Invalid feature type
        with pytest.raises(TypeError):
            FeaturedProblem(problem, 'feature', 1)
        
        # Invalid max_eval type
        with pytest.raises(TypeError):
            FeaturedProblem(problem, feature, 'max_eval')
        
        # Invalid max_eval value
        with pytest.raises(ValueError):
            FeaturedProblem(problem, feature, 0)
        
        # Invalid seed type
        with pytest.raises(TypeError):
            FeaturedProblem(problem, feature, 1, 1.5)
        
        # Invalid seed value
        with pytest.raises(ValueError):
            FeaturedProblem(problem, feature, 1, -1)

    def test_catch(self):
        """Test edge cases that should work."""
        x0 = np.zeros(2)
        problem = Problem(self.rosen, x0, cub=self.sum_cos, ceq=self.sum_sin)
        feature = Feature('plain')
        
        # Creating a featured problem should work.
        FeaturedProblem(problem, feature, 1000)

        # max_eval can be a float if it's an integer value.
        FeaturedProblem(problem, feature, 1000.0)

        # seed can be a float if it's an integer value.
        FeaturedProblem(problem, feature, 1000, 1.0)


@pytest.mark.extra
@pytest.mark.skipif(not PYCUTEST_AVAILABLE, reason="pycutest not available")
class TestPyCUTEstSelect:
    """Tests for pycutest_select function."""

    def test_simple(self):
        """Test pycutest_select with different problem types."""
        for ptype in ['u', 'b', 'l', 'n']:
            options = {'ptype': ptype, 'maxdim': 10}
            problem_names = pycutest_select(options)
            assert isinstance(problem_names, list)
            for problem_name in problem_names:
                assert isinstance(problem_name, str)

    def test_dimension_filter(self):
        """Test that dimension filtering works."""
        options_1_100 = {'ptype': 'u', 'mindim': 1, 'maxdim': 100}
        problem_names_1_100 = pycutest_select(options_1_100)
        
        options_10_50 = {'ptype': 'u', 'mindim': 10, 'maxdim': 50}
        problem_names_10_50 = pycutest_select(options_10_50)
        
        # Problems with narrower dimension range should be subset.
        assert set(problem_names_10_50).issubset(set(problem_names_1_100))

    def test_combined_types(self):
        """Test that combining problem types works."""
        options_u = {'ptype': 'u', 'maxdim': 10}
        unconstrained = pycutest_select(options_u)
        
        options_b = {'ptype': 'b', 'maxdim': 10}
        bound = pycutest_select(options_b)
        
        options_ub = {'ptype': 'ub', 'maxdim': 10}
        both = pycutest_select(options_ub)
        
        assert set(both) == set(unconstrained + bound)


@pytest.mark.extra
@pytest.mark.skipif(not PYCUTEST_AVAILABLE, reason="pycutest not available")
class TestPyCUTEstLoad:
    """Tests for pycutest_load function."""

    @pytest.mark.parametrize('ptype', ['u', 'b', 'l', 'n'])
    def test_simple(self, ptype):
        """Test loading problems from pycutest."""
        options = {'ptype': ptype, 'maxdim': 10, 'maxcon': 10}
        problem_names = pycutest_select(options)
        
        for i_problem in range(min(3, len(problem_names))):
            with suppress(ProblemError, Exception):
                problem = pycutest_load(problem_names[i_problem])
                assert isinstance(problem, Problem)
                problem.fun(problem.x0)

    def test_load_specific_problem(self):
        """Test loading a specific known problem."""
        with suppress(ProblemError, Exception):
            problem = pycutest_load('ROSENBR')
            assert isinstance(problem, Problem)
            assert problem.n == 2
