import os
import sys
from contextlib import redirect_stdout, redirect_stderr

import numpy as np

from .features import Feature
from .utils import FeatureName, CUTEstProblemOption, FeatureOption, ProblemError, get_logger

_cutest_problem_options = {
    CUTEstProblemOption.N_MIN.value: 1,
    CUTEstProblemOption.N_MAX.value: sys.maxsize,
    CUTEstProblemOption.M_MIN.value: 0,
    CUTEstProblemOption.M_MAX.value: sys.maxsize,
}


class Problem:
    r"""
    Optimization problem to be used in the benchmarking.

    Examples
    --------
    Consider the problem of minimizing the Rosenbrock function

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
    unconstrained, the components of the lower and upper bounds on the variables
    must be :math:`-\infty` and :math:`\infty`, respectively:

    >>> problem.lb
    array([-inf, -inf])
    >>> problem.ub
    array([inf, inf])

    The optional arguments of the constructor of the class `Problem` can be used
    to specify constraints. For example, to specify that the variables must be
    nonnegative, run:

    >>> problem = Problem(rosen, [0, 0], [0, 0])

    The lower bounds on the variables are now zero:

    >>> problem.lb
    array([0., 0.])
    >>> problem.ub
    array([inf, inf])

    The optional arguments ``a_ub``, ``b_ub``, ``a_eq``, and ``b_eq`` can be
    used to specify linear inequality and equality constraints, respectively.
    The optional arguments ``c_ub`` and ``c_eq`` can be used to specify
    nonlinear inequality and equality constraints, respectively. For example, to
    specify that the variables must satisfy the constraint
    :math:`x_1^2 + x_2^2 \le 1` and :math:`x_1^3 - x_2^2 \le 1`, run:

    >>> def c_ub(x):
    ...     return [x[0] ** 2 + x[1] ** 2 - 1, x[0] ** 3 - x[1] ** 2 - 1]
    ...
    >>> problem = Problem(rosen, [0, 0], c_ub=c_ub, num_nonlinear_ub=2)

    The instance ``problem`` can now be used to evaluate the nonlinear
    inequality constraints at any point. For example, to evaluate the nonlinear
    inequality constraints at the initial guess, run:

    >>> problem.c_ub(problem.x0)
    array([-1., -1.])

    If you do not provide the number of nonlinear inequality constraints in
    ``num_nonlinear_ub``, it will be inferred at the first call to
    ``problem.c_ub``. Nonlinear equality constraints can be specified in a
    similar way, using the optional arguments ``c_eq`` and ``num_nonlinear_eq``.
    """

    def __init__(self, fun, x0, lb=None, ub=None, a_ub=None, b_ub=None, a_eq=None, b_eq=None, c_ub=None, num_nonlinear_ub=None, c_eq=None, num_nonlinear_eq=None):
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
        a_ub : array_like, shape (num_linear_ub, n), optional
            Coefficient matrix of the linear constraints ``a_ub @ x <= b_ub``.
        b_ub : array_like, shape (num_linear_ub,), optional
            Right-hand side of the linear constraints ``a_ub @ x <= b_ub``.
        a_eq : array_like, shape (num_linear_eq, n), optional
            Coefficient matrix of the linear constraints ``a_eq @ x == b_eq``.
        b_eq : array_like, shape (num_linear_eq,), optional
            Right-hand side of the linear constraints ``a_eq @ x == b_eq``.
        c_ub : callable, optional
            Nonlinear inequality constraint ``c_ub(x) <= 0``.

                ``c_ub(x) -> array_like, shape (num_nonlinear_ub,)``

            where ``x`` is an array with shape (n,).
        c_eq : callable, optional
            Nonlinear equality constraint ``c_eq(x) == 0``.

                ``c_eq(x) -> array_like, shape (num_nonlinear_eq,)``

            where ``x`` is an array with shape (n,).
        num_nonlinear_ub : int, optional
            Number of nonlinear inequality constraints.
        num_nonlinear_eq : int, optional
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
        self._num_nonlinear_ub = num_nonlinear_ub
        if isinstance(self._num_nonlinear_ub, float) and self._num_nonlinear_ub.is_integer():
            self._num_nonlinear_ub = int(self._num_nonlinear_ub)
        if self._num_nonlinear_ub is not None and not isinstance(self._num_nonlinear_ub, int):
            raise TypeError('The argument num_nonlinear_ub must be an integer.')
        if self._num_nonlinear_ub is not None and self._num_nonlinear_ub < 0:
            raise ValueError('The argument num_nonlinear_ub must be nonnegative.')
        self._num_nonlinear_eq = num_nonlinear_eq
        if isinstance(self._num_nonlinear_eq, float) and self._num_nonlinear_eq.is_integer():
            self._num_nonlinear_eq = int(self._num_nonlinear_eq)
        if self._num_nonlinear_eq is not None and not isinstance(self._num_nonlinear_eq, int):
            raise TypeError('The argument num_nonlinear_eq must be an integer.')
        if self._num_nonlinear_eq is not None and self._num_nonlinear_eq < 0:
            raise ValueError('The argument num_nonlinear_eq must be nonnegative.')

        # Check that the arguments are consistent.
        if self.lb.size != self.dimension:
            raise ValueError(f'The argument lb must have size {self.dimension}.')
        if self.ub.size != self.dimension:
            raise ValueError(f'The argument ub must have size {self.dimension}.')
        if self.a_ub.shape != (self.num_linear_ub, self.dimension):
            raise ValueError(f'The argument a_ub must have shape {(self.num_linear_ub, self.dimension)}.')
        if self.a_eq.shape != (self.num_linear_eq, self.dimension):
            raise ValueError(f'The argument a_eq must have shape {(self.num_linear_eq, self.dimension)}.')
        if self._c_ub is None and self._num_nonlinear_ub is not None and self._num_nonlinear_ub > 0:
            raise ValueError('The argument num_nonlinear_ub must be None or zero if the argument c_ub is None.')
        if self._c_eq is None and self._num_nonlinear_eq is not None and self._num_nonlinear_eq > 0:
            raise ValueError('The argument num_nonlinear_eq must be None or zero if the argument c_eq is None.')

    @property
    def dimension(self):
        """
        Dimension of the problem.

        Returns
        -------
        int
            Dimension of the problem.
        """
        return self.x0.size

    @property
    def num_linear_ub(self):
        """
        Number of linear inequality constraints.

        Returns
        -------
        int
            Number of linear inequality constraints.
        """
        return self.b_ub.size

    @property
    def num_linear_eq(self):
        """
        Number of linear equality constraints.

        Returns
        -------
        int
            Number of linear equality constraints.
        """
        return self.b_eq.size

    @property
    def num_nonlinear_ub(self):
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
            `num_nonlinear_ub` was not specified when the problem was
            initialized, a nonlinear inequality constraint function was
            specified when the problem was initialized, and the method `c_ub`
            has never been called.
        """
        if self._num_nonlinear_ub is None:
            if self._c_ub is None:
                self._num_nonlinear_ub = 0
            else:
                raise ValueError('The number of nonlinear inequality constraints is unknown.')
        return self._num_nonlinear_ub

    @property
    def num_nonlinear_eq(self):
        """
        Number of nonlinear equality constraints.

        Returns
        -------
        int
            Number of nonlinear equality constraints.

        Raises
        ------
        ValueError
            If the number of nonlinear equality constraints is unknown. This can
            happen if the following three conditions are met: the argument
            `num_nonlinear_eq` was not specified when the problem was
            initialized, a nonlinear equality constraint function was specified
            when the problem was initialized, and the method `c_eq` has never
            been called.
        """
        if self._num_nonlinear_eq is None:
            if self._c_eq is None:
                self._num_nonlinear_eq = 0
            else:
                raise ValueError('The number of nonlinear equality constraints is unknown.')
        return self._num_nonlinear_eq

    @property
    def x0(self):
        """
        Initial guess.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Initial guess.
        """
        return self._x0

    @property
    def lb(self):
        """
        Lower bounds on the variables.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Lower bounds on the variables.
        """
        return self._lb if self._lb is not None else np.full(self.dimension, -np.inf)

    @property
    def ub(self):
        """
        Upper bounds on the variables.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Upper bounds on the variables.
        """
        return self._ub if self._ub is not None else np.full(self.dimension, np.inf)

    @property
    def a_ub(self):
        """
        Coefficient matrix of the linear constraints ``a_ub @ x <= b_ub``.

        Returns
        -------
        `numpy.ndarray`, shape (num_linear_ub, n)
            Coefficient matrix of the linear inequality constraints.
        """
        return self._a_ub if self._a_ub is not None else np.empty((0, self.dimension))

    @property
    def b_ub(self):
        """
        Right-hand side of the linear constraints ``a_ub @ x <= b_ub``.

        Returns
        -------
        `numpy.ndarray`, shape (num_linear_ub,)
            Right-hand side of the linear inequality constraints.
        """
        return self._b_ub if self._b_ub is not None else np.empty(0)

    @property
    def a_eq(self):
        """
        Coefficient matrix of the linear constraints ``a_eq @ x == b_eq``.

        Returns
        -------
        `numpy.ndarray`, shape (num_linear_eq, n)
            Coefficient matrix of the linear equality constraints.
        """
        return self._a_eq if self._a_eq is not None else np.empty((0, self.dimension))

    @property
    def b_eq(self):
        """
        Right-hand side of the linear constraints ``a_eq @ x == b_eq``.

        Returns
        -------
        `numpy.ndarray`, shape (num_linear_eq,)
            Right-hand side of the linear equality constraints.
        """
        return self._b_eq if self._b_eq is not None else np.empty(0)

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
        if x.size != self.dimension:
            raise ValueError(f'The argument x must have size {self.dimension}.')
        try:
            f = self._fun(x)
        except Exception as err:
            logger = get_logger(__name__)
            logger.warning(f'Failed to evaluate the objective function: {err}')
            f = np.nan
        f = float(f)
        return f

    def c_ub(self, x):
        """
        Evaluate the nonlinear constraints ``c_ub(x) <= 0``.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the nonlinear inequality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (num_nonlinear_ub,)
            Values of the nonlinear inequality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `c_ub` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument x must be a one-dimensional array.')
        if x.size != self.dimension:
            raise ValueError(f'The argument x must have size {self.dimension}.')
        if self._c_ub is None:
            c = np.empty(0)
        else:
            try:
                c = self._c_ub(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the nonlinear inequality constraint function: {err}')
                c = np.full(self.num_nonlinear_ub, np.nan)
            c = _process_1d_array(c, 'The return value of the argument c_ub must be a one-dimensional array.')
        if self._num_nonlinear_ub is None:
            self._num_nonlinear_ub = c.size
        if c.size != self.num_nonlinear_ub:
            raise ValueError(f'The return value of the argument c_ub must have size {self.num_nonlinear_ub}.')
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
        `numpy.ndarray`, shape (num_nonlinear_eq,)
            Values of the nonlinear equality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `c_eq` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument x must be a one-dimensional array.')
        if x.size != self.dimension:
            raise ValueError(f'The argument x must have size {self.dimension}.')
        if self._c_eq is None:
            c = np.empty(0)
        else:
            try:
                c = self._c_eq(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the nonlinear equality constraint function: {err}')
                c = np.full(self.num_nonlinear_eq, np.nan)
            c = _process_1d_array(c, 'The return value of the argument c_eq must be a one-dimensional array.')
        if self._num_nonlinear_eq is None:
            self._num_nonlinear_eq = c.size
        if c.size != self.num_nonlinear_eq:
            raise ValueError(f'The return value of the argument c_eq must have size {self.num_nonlinear_eq}.')
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
        x = _process_1d_array(x, 'The argument x must be a one-dimensional array.')
        if x.size != self.dimension:
            raise ValueError(f'The argument x must have size {self.dimension}.')
        cv = np.max(self.lb - x, initial=0.0)
        cv = np.max(x - self.ub, initial=cv)
        cv = np.max(self.a_ub @ x - self.b_ub, initial=cv)
        cv = np.max(np.abs(self.a_eq @ x - self.b_eq), initial=cv)
        cv = np.max(self.c_ub(x), initial=cv)
        cv = np.max(np.abs(self.c_eq(x)), initial=cv)
        return cv


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

        # Store the objective function values and maximum constraint violations.
        self._fun_history = []
        self._maxcv_history = []

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
        return len(self._fun_history)

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
        if self._feature.name == FeatureName.RANDOMIZE_X0:
            rng = self._feature.default_rng(self._seed, *x0)
            x0 += self._feature.options[FeatureOption.DISTRIBUTION](rng, x0.size)
        return x0

    @property
    def fun_history(self):
        """
        History of objective function values.

        Returns
        -------
        `numpy.ndarray`, shape (n_eval,)
            History of objective function values.
        """
        return np.array(self._fun_history, dtype=float)

    @property
    def maxcv_history(self):
        """
        History of maximum constraint violations.

        Returns
        -------
        `numpy.ndarray`, shape (n_eval,)
            History of maximum constraint violations.
        """
        return np.array(self._maxcv_history, dtype=float)

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
        RuntimeError
            If the maximum number of function evaluations has been reached.
        """
        if self.n_eval >= self._max_eval:
            raise RuntimeError('The maximum number of function evaluations has been reached.')

        # Evaluate the objective function and store the results.
        f = super().fun(x)
        self._fun_history.append(f)
        self._maxcv_history.append(self.maxcv(x))

        # Modified the objective function value according to the feature and
        # return the modified value. We should not store the modified value
        # because the performance of an optimization solver should be measured
        # using the original objective function.
        return self._feature.modifier(x, f, self._seed)


def get_cutest_problem_options():
    return dict(_cutest_problem_options)


def set_cutest_problem_options(**problem_options):
    for option_key, option_value in problem_options.items():
        if option_key == CUTEstProblemOption.N_MIN and isinstance(option_value, float) and option_value.is_integer():
            option_value = int(option_value)
        if option_key == CUTEstProblemOption.N_MIN and not isinstance(option_value, int):
            raise TypeError(f'The argument {CUTEstProblemOption.N_MIN.value} must be an integer.')
        if option_key == CUTEstProblemOption.N_MIN and option_value < 1:
            raise ValueError(f'The argument {CUTEstProblemOption.N_MIN.value} must be positive.')
        if option_key == CUTEstProblemOption.N_MAX and isinstance(option_value, float) and option_value.is_integer():
            option_value = int(option_value)
        if option_key == CUTEstProblemOption.N_MAX and not isinstance(option_value, int):
            raise TypeError(f'The argument {CUTEstProblemOption.N_MAX.value} must be an integer.')
        if option_key == CUTEstProblemOption.N_MAX and option_value < problem_options.get(CUTEstProblemOption.N_MIN, 1):
            raise ValueError(f'The argument {CUTEstProblemOption.N_MAX.value} must be greater than or equal to {CUTEstProblemOption.N_MIN.value}.')
        if option_key == CUTEstProblemOption.M_MIN and isinstance(option_value, float) and option_value.is_integer():
            option_value = int(option_value)
        if option_key == CUTEstProblemOption.M_MIN and not isinstance(option_value, int):
            raise TypeError(f'The argument {CUTEstProblemOption.M_MIN.value} must be an integer.')
        if option_key == CUTEstProblemOption.M_MIN and option_value < 0:
            raise ValueError(f'The argument {CUTEstProblemOption.M_MIN.value} must be nonnegative.')
        if option_key == CUTEstProblemOption.M_MAX and isinstance(option_value, float) and option_value.is_integer():
            option_value = int(option_value)
        if option_key == CUTEstProblemOption.M_MAX and not isinstance(option_value, int):
            raise TypeError(f'The argument {CUTEstProblemOption.M_MAX.value} must be an integer.')
        if option_key == CUTEstProblemOption.M_MAX and option_value < problem_options.get(CUTEstProblemOption.M_MIN, 0):
            raise ValueError(f'The argument {CUTEstProblemOption.M_MAX.value} must be greater than or equal to {CUTEstProblemOption.M_MIN.value}.')
        if option_key in CUTEstProblemOption.__members__.values():
            _cutest_problem_options[option_key] = option_value
        else:
            raise ValueError(f'Unknown problem option: {option_key}.')


def find_cutest_problems(constraints):
    """
    Find the names of all the CUTEst problems that satisfy given requirements.

    .. caution::

        To use this function, you must first install
        `PyCUTEst <https://jfowkes.github.io/pycutest/>`_. Follow the
        instructions carefully, as the CUTEst library must be installed in order
        to use `PyCUTEst <https://jfowkes.github.io/pycutest/>`_.

    Parameters
    ----------
    constraints : str
        Type of constraints that the CUTEst problems must have. It should
        contain one or more of the following substrings: 'unconstrained',
        'fixed', 'bound', 'adjacency', 'linear', 'quadratic', 'other'.

    Returns
    -------
    list of str
        Names of all the CUTEst problems that satisfy the given requirements.

    Raises
    ------
    TypeError
        If an argument received an invalid value.
    ValueError
        If the arguments are inconsistent.

    Examples
    --------
    To find all the unconstrained problems with at most 100 variables, use:


    >>> from optiprofiler import set_cutest_problem_options, find_cutest_problems
    >>>
    >>> set_cutest_problem_options(n_max=100)
    >>> problem_names = find_cutest_problems('unconstrained')
    """
    import pycutest

    # Preprocess the constraints.
    if not isinstance(constraints, str):
        raise TypeError('The argument constraints must be a string.')
    for constraint in constraints.split():
        if constraint not in ['unconstrained', 'fixed', 'bound', 'adjacency', 'linear', 'quadratic', 'other']:
            raise ValueError(f'Unknown constraint: {constraint}.')

    # Find all the problems that satisfy the constraints.
    problem_names = pycutest.find_problems(objective='constant linear quadratic sum of squares other', constraints=constraints, n=[_cutest_problem_options[CUTEstProblemOption.N_MIN], _cutest_problem_options[CUTEstProblemOption.N_MAX]], userN=False, m=[_cutest_problem_options[CUTEstProblemOption.M_MIN], _cutest_problem_options[CUTEstProblemOption.M_MAX]], userM=False)
    return sorted(set(problem_names))


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

    def _is_valid(cutest_problem):
        """
        Check if a CUTEst problem is valid.
        """
        # Check that all the variables are continuous.
        is_valid = np.all(cutest_problem.vartype == 0)

        # Check that the dimensions are within the specified range.
        is_valid = is_valid and _cutest_problem_options[CUTEstProblemOption.N_MIN] <= cutest_problem.n <= _cutest_problem_options[CUTEstProblemOption.N_MAX]
        is_valid = is_valid and _cutest_problem_options[CUTEstProblemOption.M_MIN] <= cutest_problem.m <= _cutest_problem_options[CUTEstProblemOption.M_MAX]

        return is_valid

    def _build_linear_ub(cutest_problem):
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

    def _build_linear_eq(cutest_problem):
        """
        Build the linear equality constraints from a CUTEst problem.
        """
        idx_eq = cutest_problem.is_linear_cons & cutest_problem.is_eq_cons
        a_eq = []
        b_eq = []
        for index in np.flatnonzero(idx_eq):
            c_val, g_val = cutest_problem.cons(np.zeros(cutest_problem.n), index, True)
            a_eq.append(g_val)
            b_eq.append(c_val - 0.5 * (cutest_problem.cl[index] + cutest_problem.cu[index]))
        return np.reshape(a_eq, (-1, cutest_problem.n)), np.array(b_eq)

    def _c_ub(cutest_problem, x):
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

    def _c_eq(cutest_problem, x):
        """
        Evaluate the nonlinear equality constraints of a CUTEst problem.
        """
        idx_eq = ~cutest_problem.is_linear_cons & cutest_problem.is_eq_cons
        c = []
        for index in np.flatnonzero(idx_eq):
            c_val = cutest_problem.cons(x, index)
            c.append(c_val - 0.5 * (cutest_problem.cl[index] + cutest_problem.cu[index]))
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
    elif cutest_problem is not None and not _is_valid(cutest_problem):
        logger.warning(f'CUTEst problem {problem_name} successfully loaded but invalid; it is discarded.')
        raise ProblemError(f'CUTEst problem {problem_name} is invalid.')

    # The problem is successfully loaded and valid.
    lb = np.array(cutest_problem.bl)
    lb[lb <= -1e20] = -np.inf
    ub = np.array(cutest_problem.bu)
    ub[ub >= 1e20] = np.inf
    if cutest_problem.m > 0:
        constraints = {
            'c_ub': lambda x: _c_ub(cutest_problem, x),
            'c_eq': lambda x: _c_eq(cutest_problem, x),
        }
        constraints['a_ub'], constraints['b_ub'] = _build_linear_ub(cutest_problem)
        constraints['a_eq'], constraints['b_eq'] = _build_linear_eq(cutest_problem)
        idx_ub = ~(cutest_problem.is_linear_cons | cutest_problem.is_eq_cons)
        constraints['num_nonlinear_ub'] = np.count_nonzero(cutest_problem.cl[idx_ub] > -1e20) + np.count_nonzero(cutest_problem.cu[idx_ub] < 1e20)
        constraints['num_nonlinear_eq'] = np.count_nonzero(~cutest_problem.is_linear_cons & cutest_problem.is_eq_cons)
    else:
        constraints = {
            'num_nonlinear_ub': 0,
            'num_nonlinear_eq': 0,
        }
    logger.info(f'CUTEst problem {cutest_problem.name} (n={cutest_problem.n}, m={cutest_problem.m}) successfully loaded.')
    return Problem(cutest_problem.obj, cutest_problem.x0, lb, ub, **constraints)


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
