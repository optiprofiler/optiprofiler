import os
import sys
import warnings
from contextlib import redirect_stdout, redirect_stderr

import numpy as np
from scipy.linalg import lstsq, qr
from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint, minimize

from .features import Feature
from .utils import FeatureName, CUTEstProblemOption, FeatureOption, ProblemError, get_logger

# Options for the CUTEst problems.
_cutest_problem_options = {
    CUTEstProblemOption.N_MIN.value: 1,
    CUTEstProblemOption.N_MAX.value: sys.maxsize,
    CUTEstProblemOption.M_BOUND_MIN: 0,
    CUTEstProblemOption.M_BOUND_MAX: sys.maxsize,
    CUTEstProblemOption.M_LINEAR_MIN.value: 0,
    CUTEstProblemOption.M_LINEAR_MAX.value: sys.maxsize,
    CUTEstProblemOption.M_LINEAR_INEQUALITY_MIN.value: 0,
    CUTEstProblemOption.M_LINEAR_INEQUALITY_MAX.value: sys.maxsize,
    CUTEstProblemOption.M_LINEAR_EQUALITY_MIN.value: 0,
    CUTEstProblemOption.M_LINEAR_EQUALITY_MAX.value: sys.maxsize,
    CUTEstProblemOption.M_NONLINEAR_MIN.value: 0,
    CUTEstProblemOption.M_NONLINEAR_MAX.value: sys.maxsize,
    CUTEstProblemOption.M_NONLINEAR_INEQUALITY_MIN.value: 0,
    CUTEstProblemOption.M_NONLINEAR_INEQUALITY_MAX.value: sys.maxsize,
    CUTEstProblemOption.M_NONLINEAR_EQUALITY_MIN.value: 0,
    CUTEstProblemOption.M_NONLINEAR_EQUALITY_MAX.value: sys.maxsize,
    CUTEstProblemOption.ALL_VARIABLES_CONTINUOUS.value: True,
}


class Problem:
    r"""
    Optimization problem to be used in the benchmarking.

    This class provides a uniform interface for providing general optimization
    problems. It is used to supply custom problems to the `run_benchmark`
    function, and the signature of solvers supplied to the `run_benchmark`
    function can be

        ``solver(problem) -> array_like, shape (n,)``

    where `problem` is an instance of the class `Problem`.

    Examples
    --------
    Consider the unconstrained problem of minimizing the Rosenbrock function

    .. math::

        f(x) = 100 (x_2 - x_1^2)^2 + (1 - x_1)^2.

    To create an instance of the class `Problem` for this problem, run:

    >>> from optiprofiler import Problem
    >>>
    >>> def rosen(x):
    ...     return 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2
    ...
    >>> problem = Problem(rosen, [0, 0])

    The second argument ``[0, 0]`` is the initial guess. This instance can now
    be used to evaluate the objective function at any point and access extra
    information about the problem. For example, to evaluate the objective
    function at the initial guess, run:

    >>> problem.fun(problem.x0)
    1.0

    All the attributes and methods of the class `Problem` described below are
    always well-defined. For example, since the problem ``problem`` is
    unconstrained, the components of the lower and upper bounds on the
    variables must be :math:`-\infty` and :math:`\infty`, respectively:

    >>> problem.lb
    array([-inf, -inf])
    >>> problem.ub
    array([inf, inf])

    The optional arguments of the constructor of the class `Problem` can be
    used to specify constraints. For example, to specify that the variables
    must be nonnegative, run:

    >>> problem = Problem(rosen, [0, 0], [0, 0])

    The lower bounds on the variables are now zero:

    >>> problem.lb
    array([0., 0.])
    >>> problem.ub
    array([inf, inf])

    The optional arguments ``a_ub``, ``b_ub``, ``a_eq``, and ``b_eq`` can be
    used to specify linear inequality and equality constraints, respectively.
    The optional arguments ``c_ub`` and ``c_eq`` can be used to specify
    nonlinear inequality and equality constraints, respectively. For example,
    to specify that the variables must satisfy the constraint
    :math:`x_1^2 + x_2^2 \le 1` and :math:`x_1^3 - x_2^2 \le 1`, run:

    >>> def c_ub(x):
    ...     return [x[0] ** 2 + x[1] ** 2 - 1, x[0] ** 3 - x[1] ** 2 - 1]
    ...
    >>> problem = Problem(rosen, [0, 0], c_ub=c_ub, m_nonlinear_ub=2)

    The instance ``problem`` can now be used to evaluate the nonlinear
    inequality constraints at any point. For example, to evaluate the nonlinear
    inequality constraints at the initial guess, run:

    >>> problem.c_ub(problem.x0)
    array([-1., -1.])

    If you do not provide the number of nonlinear inequality constraints in
    ``m_nonlinear_ub``, it will be inferred at the first call to ``c_ub``.
    Nonlinear equality constraints can be specified in a similar way, using the
    optional arguments ``c_eq`` and ``m_nonlinear_eq``.
    """

    def __init__(self, fun, x0, lb=None, ub=None, a_ub=None, b_ub=None, a_eq=None, b_eq=None, c_ub=None, c_eq=None, m_nonlinear_ub=None, m_nonlinear_eq=None):
        """
        Initialize an optimization problem.

        Parameters
        ----------
        fun : callable
            Objective function to be minimized.

                ``fun(x) -> float``

            where ``x`` is an array with shape (n,).
        x0 : array_like, shape (n,)
            Initial guess.
        lb : array_like, shape (n,), optional
            Lower bounds on the variables ``lb <= x``.
        ub : array_like, shape (n,), optional
            Upper bounds on the variables ``x <= ub``.
        a_ub : array_like, shape (m_linear_ub, n), optional
            Coefficient matrix of the linear constraints ``a_ub @ x <= b_ub``.
        b_ub : array_like, shape (m_linear_ub,), optional
            Right-hand side of the linear constraints ``a_ub @ x <= b_ub``.
        a_eq : array_like, shape (m_linear_eq, n), optional
            Coefficient matrix of the linear constraints ``a_eq @ x == b_eq``.
        b_eq : array_like, shape (m_linear_eq,), optional
            Right-hand side of the linear constraints ``a_eq @ x == b_eq``.
        c_ub : callable, optional
            Nonlinear inequality constraint ``c_ub(x) <= 0``.

                ``c_ub(x) -> array_like, shape (m_nonlinear_ub,)``

            where ``x`` is an array with shape (n,).
        c_eq : callable, optional
            Nonlinear equality constraint ``c_eq(x) == 0``.

                ``c_eq(x) -> array_like, shape (m_nonlinear_eq,)``

            where ``x`` is an array with shape (n,).
        m_nonlinear_ub : int, optional
            Number of nonlinear inequality constraints.
        m_nonlinear_eq : int, optional
            Number of nonlinear equality constraints.

        Raises
        ------
        TypeError
            If an argument received an invalid value.
        ValueError
            If the arguments are inconsistent.
        """
        # Preprocess the objective function.
        self._fun = fun
        if not callable(self._fun):
            raise TypeError('The argument fun must be callable.')

        # Preprocess the initial guess.
        self._x0 = _process_1d_array(x0, 'The argument x0 must be a one-dimensional array.')

        # Preprocess the bound constraints.
        self._lb = lb
        if self._lb is not None:
            self._lb = _process_1d_array(self._lb, 'The argument lb must be a one-dimensional array.')
        self._ub = ub
        if self._ub is not None:
            self._ub = _process_1d_array(self._ub, 'The argument ub must be a one-dimensional array.')

        # Preprocess the linear constraints.
        self._a_ub = a_ub
        if self._a_ub is not None:
            self._a_ub = _process_2d_array(self._a_ub, 'The argument a_ub must be a two-dimensional array.')
        self._b_ub = b_ub
        if self._b_ub is not None:
            self._b_ub = _process_1d_array(self._b_ub, 'The argument b_ub must be a one-dimensional array.')
        self._a_eq = a_eq
        if self._a_eq is not None:
            self._a_eq = _process_2d_array(self._a_eq, 'The argument a_eq must be a two-dimensional array.')
        self._b_eq = b_eq
        if self._b_eq is not None:
            self._b_eq = _process_1d_array(self._b_eq, 'The argument b_eq must be a one-dimensional array.')

        # Preprocess the nonlinear constraints.
        self._c_ub = c_ub
        if self._c_ub is not None:
            if not callable(self._c_ub):
                raise TypeError('The argument c_ub must be callable.')
        self._c_eq = c_eq
        if self._c_eq is not None:
            if not callable(self._c_eq):
                raise TypeError('The argument c_eq must be callable.')

        # Preprocess the number of nonlinear constraints.
        self._m_nonlinear_ub = m_nonlinear_ub
        if isinstance(self._m_nonlinear_ub, float) and self._m_nonlinear_ub.is_integer():
            self._m_nonlinear_ub = int(self._m_nonlinear_ub)
        if self._m_nonlinear_ub is not None and not isinstance(self._m_nonlinear_ub, int):
            raise TypeError('The argument m_nonlinear_ub must be an integer.')
        if self._m_nonlinear_ub is not None and self._m_nonlinear_ub < 0:
            raise ValueError('The argument m_nonlinear_ub must be nonnegative.')
        self._m_nonlinear_eq = m_nonlinear_eq
        if isinstance(self._m_nonlinear_eq, float) and self._m_nonlinear_eq.is_integer():
            self._m_nonlinear_eq = int(self._m_nonlinear_eq)
        if self._m_nonlinear_eq is not None and not isinstance(self._m_nonlinear_eq, int):
            raise TypeError('The argument m_nonlinear_eq must be an integer.')
        if self._m_nonlinear_eq is not None and self._m_nonlinear_eq < 0:
            raise ValueError('The argument m_nonlinear_eq must be nonnegative.')

        # Check that the arguments are consistent.
        if self.lb.size != self.n:
            raise ValueError(f'The argument lb must have size {self.n}.')
        if self.ub.size != self.n:
            raise ValueError(f'The argument ub must have size {self.n}.')
        if self.a_ub.shape != (self.m_linear_ub, self.n):
            raise ValueError(f'The argument a_ub must have shape {(self.m_linear_ub, self.n)}.')
        if self.a_eq.shape != (self.m_linear_eq, self.n):
            raise ValueError(f'The argument a_eq must have shape {(self.m_linear_eq, self.n)}.')
        if self._c_ub is None and self._m_nonlinear_ub is not None and self._m_nonlinear_ub > 0:
            raise ValueError('The argument m_nonlinear_ub must be None or zero if the argument c_ub is None.')
        if self._c_eq is None and self._m_nonlinear_eq is not None and self._m_nonlinear_eq > 0:
            raise ValueError('The argument m_nonlinear_eq must be None or zero if the argument c_eq is None.')

    @property
    def n(self):
        """
        Dimension of the problem.

        Returns
        -------
        int
            Dimension of the problem.
        """
        return self.x0.size

    @property
    def m_linear_ub(self):
        """
        Number of linear inequality constraints.

        Returns
        -------
        int
            Number of linear inequality constraints.
        """
        return self.b_ub.size

    @property
    def m_linear_eq(self):
        """
        Number of linear equality constraints.

        Returns
        -------
        int
            Number of linear equality constraints.
        """
        return self.b_eq.size

    @property
    def m_nonlinear_ub(self):
        """
        Number of nonlinear inequality constraints.

        Returns
        -------
        int
            Number of nonlinear inequality constraints.

        Raises
        ------
        ValueError
            If the number of nonlinear inequality constraints is unknown. This
            can happen if the following three conditions are met: the argument
            `m_nonlinear_ub` was not specified when the problem was
            initialized, a nonlinear inequality constraint function was
            specified when the problem was initialized, and the method `c_ub`
            has never been called.
        """
        if self._m_nonlinear_ub is None:
            if self._c_ub is None:
                self._m_nonlinear_ub = 0
            else:
                raise ValueError('The number of nonlinear inequality constraints is unknown.')
        return self._m_nonlinear_ub

    @property
    def m_nonlinear_eq(self):
        """
        Number of nonlinear equality constraints.

        Returns
        -------
        int
            Number of nonlinear equality constraints.

        Raises
        ------
        ValueError
            If the number of nonlinear equality constraints is unknown. This
            can happen if the following three conditions are met: the argument
            `m_nonlinear_eq` was not specified when the problem was
            initialized, a nonlinear equality constraint function was specified
            when the problem was initialized, and the method `c_eq` has never
            been called.
        """
        if self._m_nonlinear_eq is None:
            if self._c_eq is None:
                self._m_nonlinear_eq = 0
            else:
                raise ValueError('The number of nonlinear equality constraints is unknown.')
        return self._m_nonlinear_eq

    @property
    def type(self):
        """
        Type of the problem.

        Returns
        -------
        str
            Type of the problem.
        """
        try:
            if self.m_nonlinear_ub + self.m_nonlinear_eq > 0:
                return 'nonlinearly constrained'
            elif self.m_linear_ub + self.m_linear_eq > 0:
                return 'linearly constrained'
            elif np.any(self.lb > -np.inf) or np.any(self.ub < np.inf):
                return 'bound-constrained'
            else:
                return 'unconstrained'
        except ValueError:
            return 'nonlinearly constrained'

    @property
    def x0(self):
        """
        Initial guess.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Initial guess.
        """
        return np.copy(self._x0)

    @property
    def lb(self):
        """
        Lower bounds on the variables.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Lower bounds on the variables.
        """
        return np.copy(self._lb) if self._lb is not None else np.full(self.n, -np.inf)

    @property
    def ub(self):
        """
        Upper bounds on the variables.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Upper bounds on the variables.
        """
        return np.copy(self._ub) if self._ub is not None else np.full(self.n, np.inf)

    @property
    def a_ub(self):
        """
        Coefficient matrix of the linear constraints ``a_ub @ x <= b_ub``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_ub, n)
            Coefficient matrix of the linear inequality constraints.
        """
        return np.copy(self._a_ub) if self._a_ub is not None else np.empty((0, self.n))

    @property
    def b_ub(self):
        """
        Right-hand side of the linear constraints ``a_ub @ x <= b_ub``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_ub,)
            Right-hand side of the linear inequality constraints.
        """
        return np.copy(self._b_ub) if self._b_ub is not None else np.empty(0)

    @property
    def a_eq(self):
        """
        Coefficient matrix of the linear constraints ``a_eq @ x == b_eq``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_eq, n)
            Coefficient matrix of the linear equality constraints.
        """
        return np.copy(self._a_eq) if self._a_eq is not None else np.empty((0, self.n))

    @property
    def b_eq(self):
        """
        Right-hand side of the linear constraints ``a_eq @ x == b_eq``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_eq,)
            Right-hand side of the linear equality constraints.
        """
        return np.copy(self._b_eq) if self._b_eq is not None else np.empty(0)

    def fun(self, x):
        """
        Evaluate the objective function.

        The optimization problem is to minimize the objective function.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the objective function.

        Returns
        -------
        float
            Value of the objective function at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument x must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument x must have size {self.n}.')
        try:
            f = self._fun(x)
        except Exception as err:
            logger = get_logger(__name__)
            logger.warning(f'Failed to evaluate the objective function: {err}')
            f = np.nan
        return float(f)

    def c_ub(self, x):
        """
        Evaluate the nonlinear constraints ``c_ub(x) <= 0``.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the nonlinear inequality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (m_nonlinear_ub,)
            Values of the nonlinear inequality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `c_ub` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument x must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument x must have size {self.n}.')
        if self._c_ub is None:
            c = np.empty(0)
        else:
            try:
                c = self._c_ub(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the nonlinear inequality constraint function: {err}')
                c = np.full(self.m_nonlinear_ub, np.nan)
            c = _process_1d_array(c, 'The return value of the argument c_ub must be a one-dimensional array.')
        if self._m_nonlinear_ub is None:
            self._m_nonlinear_ub = c.size
        if c.size != self.m_nonlinear_ub:
            raise ValueError(f'The return value of the argument c_ub must have size {self.m_nonlinear_ub}.')
        return c

    def c_eq(self, x):
        """
        Evaluate the nonlinear constraints ``c_eq(x) == 0``.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the nonlinear equality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (m_nonlinear_eq,)
            Values of the nonlinear equality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `c_eq` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument x must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument x must have size {self.n}.')
        if self._c_eq is None:
            c = np.empty(0)
        else:
            try:
                c = self._c_eq(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the nonlinear equality constraint function: {err}')
                c = np.full(self.m_nonlinear_eq, np.nan)
            c = _process_1d_array(c, 'The return value of the argument c_eq must be a one-dimensional array.')
        if self._m_nonlinear_eq is None:
            self._m_nonlinear_eq = c.size
        if c.size != self.m_nonlinear_eq:
            raise ValueError(f'The return value of the argument c_eq must have size {self.m_nonlinear_eq}.')
        return c

    def maxcv(self, x):
        """
        Evaluate the maximum constraint violation.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the maximum constraint violation.

        Returns
        -------
        float
            Maximum constraint violation.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape.
        """
        return self._maxcv(x)[0]

    def project_x0(self):
        """
        Project the initial guess onto the feasible region.
        """
        if self.type == 'bound-constrained':
            self._x0 = np.clip(self._x0, self.lb, self.ub)
        elif self.type == 'linearly constrained' and self.m_linear_ub == 0 and np.all(self.lb == -np.inf) and np.all(self.ub == np.inf):
            self._x0 = lstsq(self.a_eq, self.b_eq - self.a_eq @ self.x0)[0]
        elif self.type != 'unconstrained':
            bounds = Bounds(self.lb, self.ub)
            constraints = []
            if self.m_linear_ub > 0:
                constraints.append(LinearConstraint(self.a_ub, -np.inf, self.b_ub))
            if self.m_linear_eq > 0:
                constraints.append(LinearConstraint(self.a_eq, self.b_eq, self.b_eq))
            if self.m_nonlinear_ub > 0:
                c_ub_x0 = self.c_ub(self.x0)
                constraints.append(NonlinearConstraint(self.c_ub, -np.inf, np.zeros_like(c_ub_x0)))
            if self.m_nonlinear_eq > 0:
                c_eq_x0 = self.c_eq(self.x0)
                constraints.append(NonlinearConstraint(self.c_eq, np.zeros_like(c_eq_x0), np.zeros_like(c_eq_x0)))

            def dist_x0_sq(x):
                g = x - self.x0
                return 0.5 * (g @ g), g

            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = minimize(dist_x0_sq, self.x0, jac=True, hessp=lambda p: p, bounds=bounds, constraints=constraints, tol=1e-16)
            self._x0 = res.x

    def _maxcv(self, x):
        """
        Evaluate the maximum constraint violations.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the maximum constraint violation.

        Returns
        -------
        float
            Maximum constraint violation for the bound constraints.
        float
            Maximum constraint violation for the linear constraints.
        float
            Maximum constraint violation for the nonlinear constraints.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument x must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument x must have size {self.n}.')
        cv_bounds = np.max(self.lb - x, initial=0.0)
        cv_bounds = np.max(x - self.ub, initial=cv_bounds)
        cv_linear = np.max(self.a_ub @ x - self.b_ub, initial=0.0)
        cv_linear = np.max(np.abs(self.a_eq @ x - self.b_eq), initial=cv_linear)
        cv_nonlinear = np.max(self.c_ub(x), initial=0.0)
        cv_nonlinear = np.max(np.abs(self.c_eq(x)), initial=cv_nonlinear)
        cv = max(cv_bounds, cv_linear, cv_nonlinear)
        return cv, cv_bounds, cv_linear, cv_nonlinear


class FeaturedProblem(Problem):
    """
    Optimization problem to be used in the benchmarking, with extra features.
    """

    def __init__(self, problem, feature, max_eval, seed=None):
        """
        Initialize an optimization problem.

        Parameters
        ----------
        problem : `optiprofiler.problems.Problem`
            Problem to be used in the benchmarking.
        feature : `optiprofiler.features.Feature`
            Feature to be used in the benchmarking.
        max_eval : int
            Maximum number of function evaluations.
        seed : int, optional
            Seed for the random number generator.

        Raises
        ------
        TypeError
            If an argument received an invalid value.
        """
        # Preprocess the feature.
        self._feature = feature
        if not isinstance(self._feature, Feature):
            raise TypeError('The argument feature must be an instance of the class Feature.')

        # Preprocess the maximum number of function evaluations.
        self._max_eval = max_eval
        if isinstance(self._max_eval, float) and self._max_eval.is_integer():
            self._max_eval = int(self._max_eval)
        if not isinstance(self._max_eval, int):
            raise TypeError('The argument max_eval must be an integer.')
        if self._max_eval < 1:
            raise ValueError('The argument max_eval must be positive.')

        # Preprocess the seed.
        self._seed = seed
        if self._seed is not None:
            if isinstance(self._seed, float) and self._seed.is_integer():
                self._seed = int(self._seed)
            if not isinstance(self._seed, int):
                raise TypeError('The argument seed must be an integer.')
            if self._seed < 0:
                raise ValueError('The argument seed must be nonnegative.')

        # Generate a random permutation.
        rng = self._feature.get_default_rng(self._seed)
        self._permutation = None
        if feature.name == FeatureName.PERMUTED:
            self._permutation = rng.permutation(problem.n)

        # Store the objective function values and maximum constraint violations.
        self._fun_hist = []
        self._maxcv_hist = []

    def __new__(cls, problem, feature, max_eval, seed=None):
        # Preprocess the problem.
        if not isinstance(problem, Problem):
            raise TypeError('The argument problem must be an instance of the class Problem.')

        # Create a new instance of the class `FeaturedProblem` by copying the
        # attributes of the problem passed to the __init__ method.
        instance = super().__new__(cls)
        for k, v in problem.__dict__.items():
            instance.__dict__[k] = v
        return instance

    @property
    def n_eval(self):
        """
        Number of objective function evaluations.

        Returns
        -------
        int
            Number of objective function evaluations.
        """
        return len(self._fun_hist)

    @property
    def x0(self):
        """
        Initial guess.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Initial guess.
        """
        x0 = super().x0
        if self._feature.name == FeatureName.PERMUTED:
            x0 = x0[np.argsort(self._permutation)]
        elif self._feature.name == FeatureName.PERTURBED_X0:
            rng = self._feature.get_default_rng(self._seed, *x0)
            x0 += self._feature.options[FeatureOption.DISTRIBUTION](rng, x0.size)
        return x0

    @property
    def lb(self):
        """
        Lower bounds on the variables.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Lower bounds on the variables.
        """
        lb = super().lb
        if self._feature.name == FeatureName.PERMUTED:
            lb = lb[np.argsort(self._permutation)]
        return lb

    @property
    def ub(self):
        """
        Upper bounds on the variables.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Upper bounds on the variables.
        """
        ub = super().ub
        if self._feature.name == FeatureName.PERMUTED:
            ub = ub[np.argsort(self._permutation)]
        return ub

    @property
    def a_ub(self):
        """
        Coefficient matrix of the linear constraints ``a_ub @ x <= b_ub``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_ub, n)
            Coefficient matrix of the linear inequality constraints.
        """
        a_ub = super().a_ub
        if self._feature.name == FeatureName.PERMUTED:
            a_ub = a_ub[:, np.argsort(self._permutation)]
        return a_ub

    @property
    def a_eq(self):
        """
        Coefficient matrix of the linear constraints ``a_eq @ x == b_eq``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_eq, n)
            Coefficient matrix of the linear equality constraints.
        """
        a_eq = super().a_eq
        if self._feature.name == FeatureName.PERMUTED:
            a_eq = a_eq[:, np.argsort(self._permutation)]
        return a_eq

    @property
    def fun_hist(self):
        """
        History of objective function values.

        Returns
        -------
        `numpy.ndarray`, shape (n_eval,)
            History of objective function values.
        """
        return np.array(self._fun_hist)

    @property
    def maxcv_hist(self):
        """
        History of maximum constraint violations.

        Returns
        -------
        `numpy.ndarray`, shape (n_eval,)
            History of maximum constraint violations.
        """
        return np.array(self._maxcv_hist)

    def fun(self, x):
        """
        Evaluate the objective function.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the objective function.

        Returns
        -------
        float
            Value of the objective function at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape.
        StopIteration
            If the maximum number of function evaluations has been reached.
        """
        if self.n_eval >= self._max_eval:
            raise StopIteration('The maximum number of function evaluations has been reached.')

        # Permutate the variables if necessary.
        if self._feature.name == FeatureName.PERMUTED:
            x = x[self._permutation]

        # Evaluate the objective function and store the results.
        f = super().fun(x)
        maxcv, maxcv_bounds, maxcv_linear, maxcv_nonlinear = self._maxcv(x)
        self._fun_hist.append(f)
        self._maxcv_hist.append(maxcv)

        # Modified the objective function value according to the feature and
        # return the modified value. We should not store the modified value
        # because the performance of an optimization solver should be measured
        # using the original objective function.
        return self._feature.modifier(x, f, maxcv_bounds, maxcv_linear, maxcv_nonlinear, self._seed)

    def c_ub(self, x):
        """
        Evaluate the nonlinear constraints ``c_ub(x) <= 0``.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the nonlinear inequality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (m_nonlinear_ub,)
            Values of the nonlinear inequality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `c_ub` has an invalid shape.
        """
        # Permutate the variables if necessary.
        if self._feature.name == FeatureName.PERMUTED:
            x = x[self._permutation]

        # Evaluate the nonlinear inequality constraints and return.
        return super().c_ub(x)

    def c_eq(self, x):
        """
        Evaluate the nonlinear constraints ``c_eq(x) == 0``.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the nonlinear equality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (m_nonlinear_eq,)
            Values of the nonlinear equality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `c_eq` has an invalid shape.
        """
        # Permutate the variables if necessary.
        if self._feature.name == FeatureName.PERMUTED:
            x = x[self._permutation]

        # Evaluate the nonlinear inequality constraints and return.
        return super().c_eq(x)


def get_cutest_problem_options():
    """
    Get the options used for loading CUTEst problems.

    Returns
    -------
    dict
        Options used for loading CUTEst problems.

    See Also
    --------
    set_cutest_problem_options :
        Set the options used for loading CUTEst problems.
    find_cutest_problems :
        Find the CUTEst problems that satisfy given requirements.
    """
    return dict(_cutest_problem_options)


def set_cutest_problem_options(**problem_options):
    """
    Set the options used for loading CUTEst problems.

    Other Parameters
    ----------------
    n_min : int, optional
        Minimum number of variables.
    n_max : int, optional
        Maximum number of variables.
    m_bound_min : int, optional
        Minimum number of bound constraints.
    m_bound_max : int, optional
        Maximum number of bound constraints.
    m_linear_min : int, optional
        Minimum number of linear constraints.
    m_linear_max : int, optional
        Maximum number of linear constraints.
    m_linear_inequality_min : int, optional
        Minimum number of linear inequality constraints.
    m_linear_inequality_max : int, optional
        Maximum number of linear inequality constraints.
    m_linear_equality_min : int, optional
        Minimum number of linear equality constraints.
    m_linear_equality_max : int, optional
        Maximum number of linear equality constraints.
    m_nonlinear_min : int, optional
        Minimum number of nonlinear constraints.
    m_nonlinear_max : int, optional
        Maximum number of nonlinear constraints.
    m_nonlinear_inequality_min : int, optional
        Minimum number of nonlinear inequality constraints.
    m_nonlinear_inequality_max : int, optional
        Maximum number of nonlinear inequality constraints.
    m_nonlinear_equality_min : int, optional
        Minimum number of nonlinear equality constraints.
    m_nonlinear_equality_max : int, optional
        Maximum number of nonlinear equality constraints.
    all_variables_continuous : bool, optional
        Whether all variables are continuous.

    See Also
    --------
    get_cutest_problem_options :
        Get the options used for loading CUTEst problems.
    find_cutest_problems :
        Find the CUTEst problems that satisfy given requirements.
    """
    for option_key, option_value in problem_options.items():
        if option_key in [
            CUTEstProblemOption.N_MIN,
            CUTEstProblemOption.N_MAX,
            CUTEstProblemOption.M_BOUND_MIN,
            CUTEstProblemOption.M_BOUND_MAX,
            CUTEstProblemOption.M_LINEAR_MIN,
            CUTEstProblemOption.M_LINEAR_MAX,
            CUTEstProblemOption.M_LINEAR_INEQUALITY_MIN,
            CUTEstProblemOption.M_LINEAR_INEQUALITY_MAX,
            CUTEstProblemOption.M_LINEAR_EQUALITY_MIN,
            CUTEstProblemOption.M_LINEAR_EQUALITY_MAX,
            CUTEstProblemOption.M_NONLINEAR_MIN,
            CUTEstProblemOption.M_NONLINEAR_MAX,
            CUTEstProblemOption.M_NONLINEAR_INEQUALITY_MIN,
            CUTEstProblemOption.M_NONLINEAR_INEQUALITY_MAX,
            CUTEstProblemOption.M_NONLINEAR_EQUALITY_MIN,
            CUTEstProblemOption.M_NONLINEAR_EQUALITY_MAX,
        ] and isinstance(option_value, float) and option_value.is_integer():
            option_value = int(option_value)
        if option_key in [
            CUTEstProblemOption.N_MIN,
            CUTEstProblemOption.N_MAX,
            CUTEstProblemOption.M_BOUND_MIN,
            CUTEstProblemOption.M_BOUND_MAX,
            CUTEstProblemOption.M_LINEAR_MIN,
            CUTEstProblemOption.M_LINEAR_MAX,
            CUTEstProblemOption.M_LINEAR_INEQUALITY_MIN,
            CUTEstProblemOption.M_LINEAR_INEQUALITY_MAX,
            CUTEstProblemOption.M_LINEAR_EQUALITY_MIN,
            CUTEstProblemOption.M_LINEAR_EQUALITY_MAX,
            CUTEstProblemOption.M_NONLINEAR_MIN,
            CUTEstProblemOption.M_NONLINEAR_MAX,
            CUTEstProblemOption.M_NONLINEAR_INEQUALITY_MIN,
            CUTEstProblemOption.M_NONLINEAR_INEQUALITY_MAX,
            CUTEstProblemOption.M_NONLINEAR_EQUALITY_MIN,
            CUTEstProblemOption.M_NONLINEAR_EQUALITY_MAX,
        ] and not isinstance(option_value, int):
            raise TypeError(f'The argument {option_key} must be an integer.')
        if option_key == CUTEstProblemOption.N_MIN and option_value < 1:
            raise ValueError(f'The argument {option_key} must be positive.')
        if option_key == CUTEstProblemOption.N_MAX and option_value < problem_options.get(CUTEstProblemOption.N_MIN, 1):
            raise ValueError(f'The argument {option_key} must be greater than or equal to {CUTEstProblemOption.N_MIN.value}.')
        if option_key in [
            CUTEstProblemOption.M_BOUND_MIN,
            CUTEstProblemOption.M_LINEAR_MIN,
            CUTEstProblemOption.M_LINEAR_INEQUALITY_MIN,
            CUTEstProblemOption.M_LINEAR_EQUALITY_MIN,
            CUTEstProblemOption.M_NONLINEAR_MIN,
            CUTEstProblemOption.M_NONLINEAR_INEQUALITY_MIN,
            CUTEstProblemOption.M_NONLINEAR_EQUALITY_MIN,
        ] and option_value < 0:
            raise ValueError(f'The argument {option_key} must be nonnegative.')
        if option_key == CUTEstProblemOption.M_BOUND_MAX and option_value < problem_options.get(CUTEstProblemOption.M_BOUND_MIN, 0):
            raise ValueError(f'The argument {option_key} must be greater than or equal to {CUTEstProblemOption.M_BOUND_MIN.value}.')
        if option_key == CUTEstProblemOption.M_LINEAR_MAX and option_value < problem_options.get(CUTEstProblemOption.M_LINEAR_MIN, 0):
            raise ValueError(f'The argument {option_key} must be greater than or equal to {CUTEstProblemOption.M_LINEAR_MIN.value}.')
        if option_key == CUTEstProblemOption.M_LINEAR_INEQUALITY_MAX and option_value < problem_options.get(CUTEstProblemOption.M_LINEAR_INEQUALITY_MIN, 0):
            raise ValueError(f'The argument {option_key} must be greater than or equal to {CUTEstProblemOption.M_LINEAR_INEQUALITY_MIN.value}.')
        if option_key == CUTEstProblemOption.M_LINEAR_EQUALITY_MAX and option_value < problem_options.get(CUTEstProblemOption.M_LINEAR_EQUALITY_MIN, 0):
            raise ValueError(f'The argument {option_key} must be greater than or equal to {CUTEstProblemOption.M_LINEAR_EQUALITY_MIN.value}.')
        if option_key == CUTEstProblemOption.M_NONLINEAR_MAX and option_value < problem_options.get(CUTEstProblemOption.M_NONLINEAR_MIN, 0):
            raise ValueError(f'The argument {option_key} must be greater than or equal to {CUTEstProblemOption.M_NONLINEAR_MIN.value}.')
        if option_key == CUTEstProblemOption.M_NONLINEAR_INEQUALITY_MAX and option_value < problem_options.get(CUTEstProblemOption.M_NONLINEAR_INEQUALITY_MIN, 0):
            raise ValueError(f'The argument {option_key} must be greater than or equal to {CUTEstProblemOption.M_NONLINEAR_INEQUALITY_MIN.value}.')
        if option_key == CUTEstProblemOption.M_NONLINEAR_EQUALITY_MAX and option_value < problem_options.get(CUTEstProblemOption.M_NONLINEAR_EQUALITY_MIN, 0):
            raise ValueError(f'The argument {option_key} must be greater than or equal to {CUTEstProblemOption.M_NONLINEAR_EQUALITY_MIN.value}.')
        if option_key == CUTEstProblemOption.ALL_VARIABLES_CONTINUOUS and not isinstance(option_value, bool):
            raise TypeError(f'The argument {option_key} must be a boolean.')
        if option_key in CUTEstProblemOption.__members__.values():
            _cutest_problem_options[option_key] = option_value
        else:
            raise ValueError(f'Unknown problem option: {option_key}.')


def find_cutest_problems(constraints, excluded_problem_names=(), **problem_options):
    """
    Find the CUTEst problems that satisfy given requirements.

    Note that some requirements cannot be checked without loading the problems.
    Therefore, the returned list of problem names may contain problems that do
    not satisfy the given requirements. These problems will be discarded when
    they are loaded.

    .. caution::

        To use this function, you must first install
        `PyCUTEst <https://jfowkes.github.io/pycutest/>`_. Follow the
        instructions carefully, as the CUTEst library must be installed in
        order to use `PyCUTEst <https://jfowkes.github.io/pycutest/>`_.

    Parameters
    ----------
    constraints : {str, list of str}
        Type of constraints that the CUTEst problems must have. It should
        be either: ``'unconstrained'``, ``'bound'``, ``'linear'``,
        ``'nonlinear'``, or a list containing one or more of these strings.
    excluded_problem_names: list of str, optional
        Names of the CUTEst problems to exclude from the search.

    Returns
    -------
    list of str
        Names of all the CUTEst problems that satisfy the given requirements.

    Other Parameters
    ----------------
    n_min : int, optional
        Minimum number of variables.
    n_max : int, optional
        Maximum number of variables.
    m_bound_min : int, optional
        Minimum number of bound constraints.
    m_bound_max : int, optional
        Maximum number of bound constraints.
    m_linear_min : int, optional
        Minimum number of linear constraints.
    m_linear_max : int, optional
        Maximum number of linear constraints.
    m_linear_inequality_min : int, optional
        Minimum number of linear inequality constraints.
    m_linear_inequality_max : int, optional
        Maximum number of linear inequality constraints.
    m_linear_equality_min : int, optional
        Minimum number of linear equality constraints.
    m_linear_equality_max : int, optional
        Maximum number of linear equality constraints.
    m_nonlinear_min : int, optional
        Minimum number of nonlinear constraints.
    m_nonlinear_max : int, optional
        Maximum number of nonlinear constraints.
    m_nonlinear_inequality_min : int, optional
        Minimum number of nonlinear inequality constraints.
    m_nonlinear_inequality_max : int, optional
        Maximum number of nonlinear inequality constraints.
    m_nonlinear_equality_min : int, optional
        Minimum number of nonlinear equality constraints.
    m_nonlinear_equality_max : int, optional
        Maximum number of nonlinear equality constraints.
    all_variables_continuous : bool, optional
        Whether all variables are continuous.

    Raises
    ------
    TypeError
        If an argument received an invalid value.
    ValueError
        If the arguments are inconsistent.

    See Also
    --------
    get_cutest_problem_options :
        Get the options used for loading CUTEst problems.
    set_cutest_problem_options :
        Set the options used for loading CUTEst problems.

    Notes
    -----
    The problem options supplied to this function are forwards to
    `set_cutest_problem_options`. Therefore, these options will be set
    globally and will be used when loading CUTEst problems.

    Examples
    --------
    To find all the unconstrained problems with at most 100 variables, use:

    >>> from optiprofiler import find_cutest_problems
    >>>
    >>> problem_names = find_cutest_problems('unconstrained', n_max=100)
    """
    import pycutest

    # Preprocess the constraints.
    if isinstance(constraints, str):
        constraints = [constraints]
    constraints = list(constraints)
    for constraint in constraints:
        if not isinstance(constraint, str):
            raise TypeError('The argument constraints must be either a string or a list of strings.')
        if constraint not in ['unconstrained', 'bound', 'linear', 'nonlinear']:
            raise ValueError(f'Unknown constraint: {constraint}.')

    # Preprocess the excluded problem names.
    if isinstance(excluded_problem_names, str):
        excluded_problem_names = [excluded_problem_names]
    excluded_problem_names = list(excluded_problem_names)
    for excluded_problem_name in excluded_problem_names:
        if not isinstance(excluded_problem_name, str):
            raise TypeError('The argument excluded_problem_names must be a list of strings.')

    # Set the CUTEst problem options.
    set_cutest_problem_options(**problem_options)

    # Find all the problems that satisfy the constraints.
    constraints_cutest = ''
    if 'unconstrained' in constraints:
        constraints_cutest += 'unconstrained fixed'
    if 'bound' in constraints:
        constraints_cutest += 'bound'
    if 'linear' in constraints:
        constraints_cutest += 'adjacency linear'
    if 'nonlinear' in constraints:
        constraints_cutest += 'quadratic other'
    problem_names = pycutest.find_problems(
        constraints=constraints_cutest,
        n=[
            _cutest_problem_options[CUTEstProblemOption.N_MIN],
            _cutest_problem_options[CUTEstProblemOption.N_MAX],
        ],
        userN=False,
        m=[
            _cutest_problem_options[CUTEstProblemOption.M_LINEAR_MIN] + _cutest_problem_options[CUTEstProblemOption.M_NONLINEAR_MIN],
            _cutest_problem_options[CUTEstProblemOption.M_LINEAR_MAX] + _cutest_problem_options[CUTEstProblemOption.M_NONLINEAR_MAX],
        ],
        userM=False,
    )
    return sorted(set(problem_names).difference(excluded_problem_names))


def load_cutest_problem(problem_name):
    """
    Load a CUTEst problem.

    Parameters
    ----------
    problem_name : str
        Name of the CUTEst problem to load.

    Returns
    -------
    Problem
        Loaded problem.

    Raises
    ------
    TypeError
        If an argument received an invalid value.
    ValueError
        If the arguments are inconsistent.
    ProblemError
        If the problem cannot be loaded.
    """
    import pycutest

    def _is_valid():
        """
        Check if a CUTEst problem is valid.
        """
        # Check that the dimension is within the specified range.
        is_valid = _cutest_problem_options[CUTEstProblemOption.N_MIN] <= cutest_problem.n <= _cutest_problem_options[CUTEstProblemOption.N_MAX]

        # Check that all the variable types satisfy the specified requirements.
        if _cutest_problem_options[CUTEstProblemOption.ALL_VARIABLES_CONTINUOUS]:
            is_valid = is_valid and np.all(cutest_problem.vartype == 0)

        # Ensure that the problem constraints satisfy the specified requirements.
        m_bound = np.count_nonzero(cutest_problem.bl > -1e20) + np.count_nonzero(cutest_problem.bu < 1e20)
        m_linear = np.count_nonzero(cutest_problem.is_linear_cons)
        m_linear_inequality = 0 if m_linear == 0 else np.count_nonzero(cutest_problem.is_linear_cons & ~cutest_problem.is_eq_cons)
        m_linear_equality = 0 if m_linear == 0 else np.count_nonzero(cutest_problem.is_linear_cons & cutest_problem.is_eq_cons)
        m_nonlinear = cutest_problem.m - m_linear
        m_nonlinear_inequality = 0 if m_nonlinear == 0 else np.count_nonzero(~cutest_problem.is_linear_cons & ~cutest_problem.is_eq_cons)
        m_nonlinear_equality = 0 if m_nonlinear == 0 else np.count_nonzero(~cutest_problem.is_linear_cons & cutest_problem.is_eq_cons)
        is_valid = is_valid and _cutest_problem_options[CUTEstProblemOption.M_BOUND_MIN] <= m_bound <= _cutest_problem_options[CUTEstProblemOption.M_BOUND_MAX]
        is_valid = is_valid and _cutest_problem_options[CUTEstProblemOption.M_LINEAR_MIN] <= m_linear <= _cutest_problem_options[CUTEstProblemOption.M_LINEAR_MAX]
        is_valid = is_valid and _cutest_problem_options[CUTEstProblemOption.M_LINEAR_INEQUALITY_MIN] <= m_linear_inequality <= _cutest_problem_options[CUTEstProblemOption.M_LINEAR_INEQUALITY_MAX]
        is_valid = is_valid and _cutest_problem_options[CUTEstProblemOption.M_LINEAR_EQUALITY_MIN] <= m_linear_equality <= _cutest_problem_options[CUTEstProblemOption.M_LINEAR_EQUALITY_MAX]
        is_valid = is_valid and _cutest_problem_options[CUTEstProblemOption.M_NONLINEAR_MIN] <= m_nonlinear <= _cutest_problem_options[CUTEstProblemOption.M_NONLINEAR_MAX]
        is_valid = is_valid and _cutest_problem_options[CUTEstProblemOption.M_NONLINEAR_INEQUALITY_MIN] <= m_nonlinear_inequality <= _cutest_problem_options[CUTEstProblemOption.M_NONLINEAR_INEQUALITY_MAX]
        is_valid = is_valid and _cutest_problem_options[CUTEstProblemOption.M_NONLINEAR_EQUALITY_MIN] <= m_nonlinear_equality <= _cutest_problem_options[CUTEstProblemOption.M_NONLINEAR_EQUALITY_MAX]

        return is_valid

    def _build_linear_ub():
        """
        Build the linear inequality constraints from a CUTEst problem.
        """
        idx_ub = cutest_problem.is_linear_cons & ~cutest_problem.is_eq_cons
        idx_ub_cl = cutest_problem.cl[idx_ub] > -1e20
        idx_ub_cu = cutest_problem.cu[idx_ub] < 1e20
        a_ub = []
        b_ub = []
        for i, index in enumerate(np.flatnonzero(idx_ub)):
            c_val, g_val = cutest_problem.cons(np.zeros(cutest_problem.n), index, True)
            if idx_ub_cl[i]:
                a_ub.append(-g_val)
                b_ub.append(c_val - cutest_problem.cl[index])
            if idx_ub_cu[i]:
                a_ub.append(g_val)
                b_ub.append(cutest_problem.cu[index] - c_val)
        return np.reshape(a_ub, (-1, cutest_problem.n)), np.array(b_ub)

    def _build_linear_eq():
        """
        Build the linear equality constraints from a CUTEst problem.
        """
        idx_eq = cutest_problem.is_linear_cons & cutest_problem.is_eq_cons
        a_eq = []
        b_eq = []
        for index in np.flatnonzero(idx_eq):
            c_val, g_val = cutest_problem.cons(np.zeros(cutest_problem.n), index, True)
            a_eq.append(g_val)
            b_eq.append(0.5 * (cutest_problem.cl[index] + cutest_problem.cu[index]) - c_val)
        return np.reshape(a_eq, (-1, cutest_problem.n)), np.array(b_eq)

    def _c_ub(x):
        """
        Evaluate the nonlinear inequality constraints of a CUTEst problem.
        """
        idx_ub = ~(cutest_problem.is_linear_cons | cutest_problem.is_eq_cons)
        idx_ub_cl = cutest_problem.cl[idx_ub] > -1e20
        idx_ub_cu = cutest_problem.cu[idx_ub] < 1e20
        c = []
        for i, index in enumerate(np.flatnonzero(idx_ub)):
            c_val = cutest_problem.cons(x, index)
            if idx_ub_cl[i]:
                c.append(cutest_problem.cl[index] - c_val)
            if idx_ub_cu[i]:
                c.append(c_val - cutest_problem.cu[index])
        return np.array(c)

    def _c_eq(x):
        """
        Evaluate the nonlinear equality constraints of a CUTEst problem.
        """
        idx_eq = ~cutest_problem.is_linear_cons & cutest_problem.is_eq_cons
        c = []
        for index in np.flatnonzero(idx_eq):
            c_val = cutest_problem.cons(x, index)
            c.append(0.5 * (cutest_problem.cl[index] + cutest_problem.cu[index]) - c_val)
        return np.array(c)

    # Preprocess the problem name.
    if not isinstance(problem_name, str):
        raise TypeError('The argument problem_name must be a string.')

    # Attempt to load the CUTEst problem.
    cutest_problem = None
    logger = get_logger(__name__)
    logger.info(f'Loading CUTEst problem {problem_name}.')
    try:
        with open(os.devnull, 'w') as devnull:
            with redirect_stdout(devnull), redirect_stderr(devnull):
                cutest_problem = pycutest.import_problem(problem_name)
    except Exception as err:
        logger.error(f'Failed to load CUTEst problem {problem_name}: {err}')

    # If the problem is not successfully loaded or invalid, raise an exception.
    if cutest_problem is None:
        raise ProblemError(f'Failed to load CUTEst problem {problem_name}.')
    elif cutest_problem is not None and not _is_valid():
        logger.warning(f'CUTEst problem {problem_name} successfully loaded but invalid; it is discarded.')
        raise ProblemError(f'CUTEst problem {problem_name} is invalid.')

    # The problem is successfully loaded and valid.
    lb = np.array(cutest_problem.bl)
    lb[lb <= -1e20] = -np.inf
    ub = np.array(cutest_problem.bu)
    ub[ub >= 1e20] = np.inf
    if cutest_problem.m > 0:
        constraints = {'c_ub': _c_ub, 'c_eq': _c_eq}
        constraints['a_ub'], constraints['b_ub'] = _build_linear_ub()
        constraints['a_eq'], constraints['b_eq'] = _build_linear_eq()
        idx_nonlinear_ub = ~(cutest_problem.is_linear_cons | cutest_problem.is_eq_cons)
        constraints['m_nonlinear_ub'] = np.count_nonzero(cutest_problem.cl[idx_nonlinear_ub] > -1e20) + np.count_nonzero(cutest_problem.cu[idx_nonlinear_ub] < 1e20)
        constraints['m_nonlinear_eq'] = np.count_nonzero(~cutest_problem.is_linear_cons & cutest_problem.is_eq_cons)
    else:
        constraints = {'m_nonlinear_ub': 0, 'm_nonlinear_eq': 0}
    problem = Problem(cutest_problem.obj, cutest_problem.x0, lb, ub, **constraints)
    logger.info(f'{problem.type.capitalize()} CUTEst problem {cutest_problem.name} (n = {problem.n}) successfully loaded.')
    return problem


def _process_1d_array(x, message):
    """
    Preprocess a one-dimensional array.

    Parameters
    ----------
    x : array_like
        Array to preprocess.
    message : str
        Error message to raise if the array is invalid.

    Returns
    -------
    `numpy.ndarray`
        Preprocessed array.

    Raises
    ------
    ValueError
        If the array is invalid.
    """
    x = np.atleast_1d(np.squeeze(x)).astype(float)
    if x.ndim != 1:
        raise ValueError(message)
    return x


def _process_2d_array(x, message):
    """
    Preprocess a two-dimensional array.

    Parameters
    ----------
    x : array_like
        Array to preprocess.
    message : str
        Error message to raise if the array is invalid.

    Returns
    -------
    `numpy.ndarray`
        Preprocessed array.

    Raises
    ------
    ValueError
        If the array is invalid.
    """
    x = np.atleast_2d(x).astype(float)
    if x.ndim != 2:
        raise ValueError(message)
    return x
