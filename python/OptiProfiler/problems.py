import warnings
from enum import Enum

import numpy as np

from .utils import get_logger


class ProblemOptionKey(str, Enum):
    """
    Available options for loading a CUTEst problem.
    """
    N_MIN = 'n_min'
    N_MAX = 'n_max'
    M_MIN = 'm_min'
    M_MAX = 'm_max'


class Problem:
    r"""
    Optimization problem to be used in the benchmarking.

    Examples
    --------
    To illustrate the use of this class, consider the problem of minimizing the
    Rosenbrock function

    .. math::

        f(x) = 100 (x_2 - x_1^2)^2 + (1 - x_1)^2.

    To create an instance ``problem`` of the class `Problem` for this problem, run:

    >>> from OptiProfiler import Problem
    >>>
    >>> def rosen(x):
    ...     return 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2
    ...
    >>> problem = Problem(rosen, [0, 0])

    The second argument is the initial guess. The instance ``problem`` can now
    be used to evaluate the objective function at any point, and access extra
    information about the problem. For example, to evaluate the objective
    function at the initial guess, run:

    >>> problem.fun(problem.x0)
    1.0

    All the attributes and methods of the class `Problem` described below are
    always well-defined. For example, since the problem ``problem`` is
    unconstrained, the components of the lower and upper bounds on the variables
    must be :math:`-\infty` and :math:`\infty`, respectively:

    >>> problem.xl
    array([-inf, -inf])
    >>> problem.xu
    array([inf, inf])

    To optional arguments of the constructor of the class `Problem` can be used
    to specify constraints. For example, to specify that the variables must be
    nonnegative, run:

    >>> problem = Problem(rosen, [0, 0], [0, 0])

    The lower bounds on the variables are now zero:

    >>> problem.xl
    array([0., 0.])
    >>> problem.xu
    array([inf, inf])

    The optional arguments ``aub``, ``bub``, ``aeq``, and ``beq`` can be used
    to specify linear inequality and equality constraints, respectively.
    The optional arguments ``cub`` and ``ceq`` can be used to specify nonlinear
    inequality and equality constraints, respectively. For example, to specify
    that the variables must satisfy the constraint :math:`x_1^2 + x_2^2 \leq 1`
    and :math:`x_1^3 - x_2^2 \leq 1`, run:

    >>> def cub(x):
    ...     return [x[0] ** 2 + x[1] ** 2 - 1, x[0] ** 3 - x[1] ** 2 - 1]
    ...
    >>> problem = Problem(rosen, [0, 0], cub=cub, m_nonlinear_ub=2)

    The instance ``problem`` can now be used to evaluate the nonlinear
    inequality constraints at any point. For example, to evaluate the nonlinear
    inequality constraints at the initial guess, run:

    >>> problem.cub(problem.x0)
    array([-1., -1.])

    If you do not provide the number of nonlinear inequality constraints in
    ``m_nonlinear_ub``, it will be inferred at the first call to
    ``problem.cub``. Nonlinear equality constraints can be specified in a
    similar way, using the optional arguments ``ceq`` and ``m_nonlinear_eq``.
    """

    def __init__(self, fun, x0, xl=None, xu=None, aub=None, bub=None, aeq=None, beq=None, cub=None, ceq=None, m_nonlinear_ub=None, m_nonlinear_eq=None):
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
        xl : array_like, shape (n,), optional
            Lower bounds on the variables ``xl <= x``.
        xu : array_like, shape (n,), optional
            Upper bounds on the variables ``x <= xu``.
        aub : array_like, shape (m_linear_ub, n), optional
            Left-hand side matrix of the linear constraints ``aub @ x <= bub``.
        bub : array_like, shape (m_linear_ub,), optional
            Right-hand side vector of the linear constraints ``aub @ x <= bub``.
        aeq : array_like, shape (m_linear_eq, n), optional
            Left-hand side matrix of the linear constraints ``aeq @ x == beq``.
        beq : array_like, shape (m_linear_eq,), optional
            Right-hand side vector of the linear constraints ``aeq @ x == beq``.
        cub : callable, optional
            Nonlinear inequality constraint ``cub(x) <= 0``.

                ``cub(x) -> array_like, shape (m_nonlinear_ub,)``

            where ``x`` is an array with shape (n,).
        ceq : callable, optional
            Nonlinear equality constraint ``ceq(x) == 0``.

                ``ceq(x) -> array_like, shape (m_nonlinear_eq,)``

            where ``x`` is an array with shape (n,).
        m_nonlinear_ub : int, optional
            Number of nonlinear inequality constraints.
        m_nonlinear_eq : int, optional
            Number of nonlinear equality constraints.

        Raises
        ------
        TypeError
            If an argument received an invalid value.
        """
        # Preprocess the initial guess.
        self._x0 = _1d_array(x0, 'The argument x0 must be a one-dimensional array.')

        # Preprocess the objective function.
        self._fun = fun
        if not callable(self._fun):
            raise TypeError('The argument fun must be callable.')

        # Preprocess the bound constraints.
        self._xl = xl
        if self._xl is not None:
            self._xl = _1d_array(self._xl, 'The argument xl must be a one-dimensional array.')
        self._xu = xu
        if self._xu is not None:
            self._xu = _1d_array(self._xu, 'The argument xu must be a one-dimensional array.')

        # Preprocess the linear constraints.
        self._aub = aub
        if self._aub is not None:
            self._aub = _2d_array(self._aub, 'The argument aub must be a two-dimensional array.')
        self._bub = bub
        if self._bub is not None:
            self._bub = _1d_array(self._bub, 'The argument bub must be a one-dimensional array.')
        self._aeq = aeq
        if self._aeq is not None:
            self._aeq = _2d_array(self._aeq, 'The argument aeq must be a two-dimensional array.')
        self._beq = beq
        if self._beq is not None:
            self._beq = _1d_array(self._beq, 'The argument beq must be a one-dimensional array.')

        # Preprocess the nonlinear constraints.
        self._cub = cub
        if self._cub is not None:
            if not callable(self._cub):
                raise TypeError('The argument cub must be callable.')
        self._ceq = ceq
        if self._ceq is not None:
            if not callable(self._ceq):
                raise TypeError('The argument ceq must be callable.')

        # Preprocess the number of nonlinear constraints.
        self._m_nonlinear_ub = m_nonlinear_ub
        if isinstance(self._m_nonlinear_ub, float) and self._m_nonlinear_ub.is_integer():
            self._m_nonlinear_ub = int(self._m_nonlinear_ub)
        if self._m_nonlinear_ub is not None and not (isinstance(self._m_nonlinear_ub, int) and self._m_nonlinear_ub >= 0):
            raise TypeError('The argument m_nonlinear_ub must be a nonnegative integer.')
        self._m_nonlinear_eq = m_nonlinear_eq
        if isinstance(self._m_nonlinear_eq, float) and self._m_nonlinear_eq.is_integer():
            self._m_nonlinear_eq = int(self._m_nonlinear_eq)
        if self._m_nonlinear_eq is not None and not (isinstance(self._m_nonlinear_eq, int) and self._m_nonlinear_eq >= 0):
            raise TypeError('The argument m_nonlinear_eq must be a nonnegative integer.')

        # Check that the arguments are consistent.
        if self.xl.size != self.n:
            raise ValueError(f'The argument xl must have size {self.n}.')
        if self.xu.size != self.n:
            raise ValueError(f'The argument xu must have size {self.n}.')
        if self.aub.shape != (self.m_linear_ub, self.n):
            raise ValueError(f'The argument aub must have shape {(self.m_linear_ub, self.n)}.')
        if self.aeq.shape != (self.m_linear_eq, self.n):
            raise ValueError(f'The argument aeq must have shape {(self.m_linear_eq, self.n)}.')
        if self._cub is None and self._m_nonlinear_ub is not None and self._m_nonlinear_ub > 0:
            raise ValueError('The argument m_nonlinear_ub must be None or zero if the argument cub is None.')
        if self._ceq is None and self._m_nonlinear_eq is not None and self._m_nonlinear_eq > 0:
            raise ValueError('The argument m_nonlinear_eq must be None or zero if the argument ceq is None.')

    @property
    def n(self):
        """
        Number of variables.

        Returns
        -------
        int
            Number of variables.
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
        return self.bub.size

    @property
    def m_linear_eq(self):
        """
        Number of linear equality constraints.

        Returns
        -------
        int
            Number of linear equality constraints.
        """
        return self.beq.size

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
            `m_nonlinear_ub` was not specified when the problem was initialized,
            a nonlinear inequality constraint function was specified when the
            problem was initialized, and the method `cub` has never been called.
        """
        if self._m_nonlinear_ub is None:
            if self._cub is None:
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
            If the number of nonlinear equality constraints is unknown. This can
            happen if the following three conditions are met: the argument
            `m_nonlinear_eq` was not specified when the problem was initialized,
            a nonlinear equality constraint function was specified when the
            problem was initialized, and the method `ceq` has never been called.
        """
        if self._m_nonlinear_eq is None:
            if self._ceq is None:
                self._m_nonlinear_eq = 0
            else:
                raise ValueError('The number of nonlinear equality constraints is unknown.')
        return self._m_nonlinear_eq

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
    def xl(self):
        """
        Lower bounds on the variables.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Lower bounds on the variables.
        """
        return self._xl if self._xl is not None else np.full(self.n, -np.inf)

    @property
    def xu(self):
        """
        Upper bounds on the variables.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Upper bounds on the variables.
        """
        return self._xu if self._xu is not None else np.full(self.n, np.inf)

    @property
    def aub(self):
        """
        Left-hand side matrix of the linear constraints ``aub @ x <= bub``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_ub, n)
            Left-hand side matrix of the linear inequality constraints.
        """
        return self._aub if self._aub is not None else np.empty((0, self.n))

    @property
    def bub(self):
        """
        Right-hand side vector of the linear constraints ``aub @ x <= bub``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_ub,)
            Right-hand side vector of the linear inequality constraints.
        """
        return self._bub if self._bub is not None else np.empty(0)

    @property
    def aeq(self):
        """
        Left-hand side matrix of the linear constraints ``aeq @ x == beq``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_eq, n)
            Left-hand side matrix of the linear equality constraints.
        """
        return self._aeq if self._aeq is not None else np.empty((0, self.n))

    @property
    def beq(self):
        """
        Right-hand side vector of the linear constraints ``aeq @ x == beq``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_eq,)
            Right-hand side vector of the linear equality constraints.
        """
        return self._beq if self._beq is not None else np.empty(0)

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
        """
        x = _1d_array(x, 'The argument `x` must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` must have size {self.n}.')
        try:
            f = self._fun(x)
        except Exception as err:
            warnings.warn(f'Failed to evaluate the objective function: {err}', RuntimeWarning)
            f = np.nan
        f = float(f)
        return f

    def cub(self, x):
        """
        Evaluate the nonlinear constraints ``cub(x) <= 0``.

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
            the argument `cub` has an invalid shape.
        """
        x = _1d_array(x, 'The argument `x` must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` must have size {self.n}.')
        if self._cub is None:
            c = np.empty(0)
        else:
            try:
                c = self._cub(x)
            except Exception as err:
                warnings.warn(f'Failed to evaluate the nonlinear inequality constraint function: {err}', RuntimeWarning)
                c = np.full(self.m_nonlinear_ub, np.nan)
            c = _1d_array(c, 'The return value of the argument cub must be a one-dimensional array.')
        if self._m_nonlinear_ub is None:
            self._m_nonlinear_ub = c.size
        if c.size != self.m_nonlinear_ub:
            raise ValueError(f'The return value of the argument cub must have size {self.m_nonlinear_ub}.')
        return c

    def ceq(self, x):
        """
        Evaluate the nonlinear constraints ``ceq(x) == 0``.

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
            the argument `ceq` has an invalid shape.
        """
        x = _1d_array(x, 'The argument x must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument x must have size {self.n}.')
        if self._ceq is None:
            c = np.empty(0)
        else:
            try:
                c = self._ceq(x)
            except Exception as err:
                warnings.warn(f'Failed to evaluate the nonlinear equality constraint function: {err}', RuntimeWarning)
                c = np.full(self.m_nonlinear_eq, np.nan)
            c = _1d_array(c, 'The return value of the argument ceq must be a one-dimensional array.')
        if self._m_nonlinear_eq is None:
            self._m_nonlinear_eq = c.size
        if c.size != self.m_nonlinear_eq:
            raise ValueError(f'The return value of the argument ceq must have size {self.m_nonlinear_eq}.')
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
        """
        x = _1d_array(x, 'The argument x must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument x must have size {self.n}.')
        cv = np.max(self.xl - x, initial=0.0)
        cv = np.max(x - self.xu, initial=cv)
        cv = np.max(self.aub @ x - self.bub, initial=cv)
        cv = np.max(np.abs(self.aeq @ x - self.beq), initial=cv)
        cv = np.max(self.cub(x), initial=cv)
        cv = np.max(np.abs(self.ceq(x)), initial=cv)
        return cv


class FeaturedProblem(Problem):
    """
    Optimization problem to be used in the benchmarking, with extra features.
    """

    def __init__(self, problem, feature, seed=None):
        """
        Initialize an optimization problem.

        Parameters
        ----------
        problem : `OptiProfiler.problems.Problem`
            Problem to be used in the benchmarking.
        feature : `OptiProfiler.features.Feature`
            Feature to be used in the benchmarking.
        seed : int, optional
            Seed for the random number generator.
        """
        self._feature = feature
        self._seed = seed

        # Store the objective function values and maximum constraint violations.
        self._fun_values = []
        self._maxcv_values = []

    def __new__(cls, problem, feature, seed):
        """
        Initialize an optimization problem.

        Parameters
        ----------
        problem : `OptiProfiler.problems.Problem`
            Problem to be used in the benchmarking.
        feature : `OptiProfiler.features.Feature`
            Feature to be used in the benchmarking.
        seed : int, optional
            Seed for the random number generator.
        """
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
        return len(self._fun_values)

    @property
    def x0(self):
        """
        Initial guess.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Initial guess.
        """
        return super().x0

    @property
    def fun_values(self):
        """
        History of objective function values.

        Returns
        -------
        `numpy.ndarray`, shape (n_eval,)
            History of objective function values.
        """
        return np.array(self._fun_values, dtype=float)

    @property
    def maxcv_values(self):
        """
        History of maximum constraint violations.

        Returns
        -------
        `numpy.ndarray`, shape (n_eval,)
            History of maximum constraint violations.
        """
        return np.array(self._maxcv_values, dtype=float)

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
        """
        # Evaluate the objective function and store the results.
        f = super().fun(x)
        self._fun_values.append(f)
        self._maxcv_values.append(self.maxcv(x))

        # Modified the objective function value according to the feature and
        # return the modified value. We should not store the modified value
        # because the performance of an optimization solver should be measured
        # using the original objective function.
        return self._feature.modifier(x, f, self._seed)


class ProblemError(Exception):
    """
    Exception raised when a problem cannot be loaded.
    """
    pass


def load_cutest(problem_name, **problem_options):
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

    Other Parameters
    ----------------
    n_min : int, optional
        Minimum number of variables.
    n_max : int, optional
        Maximum number of variables.
    m_min : int, optional
        Minimum number of linear and nonlinear constraints.
    m_max : int, optional
        Maximum number of linear and nonlinear constraints.

    Raises
    ------
    ProblemError
        If the problem cannot be loaded.
    """
    import pycutest

    def _is_valid(cutest_problem, **problem_options):
        """
        Check if a CUTEst problem is valid.
        """
        # Check that all the variables are continuous.
        is_valid = np.all(cutest_problem.vartype == 0)

        # Check that the dimensions are within the specified range.
        if ProblemOptionKey.N_MIN in problem_options:
            is_valid = is_valid and cutest_problem.n >= problem_options[ProblemOptionKey.N_MIN]
        if ProblemOptionKey.N_MAX in problem_options:
            is_valid = is_valid and cutest_problem.n <= problem_options[ProblemOptionKey.N_MAX]
        if ProblemOptionKey.M_MIN in problem_options:
            is_valid = is_valid and cutest_problem.m >= problem_options[ProblemOptionKey.M_MIN]
        if ProblemOptionKey.M_MAX in problem_options:
            is_valid = is_valid and cutest_problem.m <= problem_options[ProblemOptionKey.M_MAX]

        return is_valid

    def _build_linear_ub(cutest_problem):
        """
        Build the linear inequality constraints from a CUTEst problem.
        """
        idx_ub = cutest_problem.is_linear_cons & ~cutest_problem.is_eq_cons
        idx_ub_cl = cutest_problem.cl[idx_ub] > -1e20
        idx_ub_cu = cutest_problem.cu[idx_ub] < 1e20
        aub = []
        bub = []
        for i, index in enumerate(np.flatnonzero(idx_ub)):
            c_val, g_val = cutest_problem.cons(np.zeros(cutest_problem.n), index, True)
            if idx_ub_cl[i]:
                aub.append(-g_val)
                bub.append(c_val - cutest_problem.cl[index])
            if idx_ub_cu[i]:
                aub.append(g_val)
                bub.append(cutest_problem.cu[index] - c_val)
        return np.reshape(aub, (-1, cutest_problem.n)), np.array(bub)

    def _build_linear_eq(cutest_problem):
        """
        Build the linear equality constraints from a CUTEst problem.
        """
        idx_eq = cutest_problem.is_linear_cons & cutest_problem.is_eq_cons
        aeq = []
        beq = []
        for index in np.flatnonzero(idx_eq):
            c_val, g_val = cutest_problem.cons(np.zeros(cutest_problem.n), index, True)
            aeq.append(g_val)
            beq.append(c_val - 0.5 * (cutest_problem.cl[index] + cutest_problem.cu[index]))
        return np.reshape(aeq, (-1, cutest_problem.n)), np.array(beq)

    def _cub(cutest_problem, x):
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

    def _ceq(cutest_problem, x):
        """
        Evaluate the nonlinear equality constraints of a CUTEst problem.
        """
        idx_eq = ~cutest_problem.is_linear_cons & cutest_problem.is_eq_cons
        c = []
        for index in np.flatnonzero(idx_eq):
            c_val = cutest_problem.cons(x, index)
            c.append(c_val - 0.5 * (cutest_problem.cl[index] + cutest_problem.cu[index]))
        return np.array(c)

    # Check the arguments.
    if not isinstance(problem_name, str):
        raise TypeError('The argument problem_name must be a string.')
    for key in problem_options:
        if key not in ProblemOptionKey.__members__.values():
            raise ValueError(f'Unknown argument: {key}.')
        if isinstance(problem_options[key], float) and problem_options[key].is_integer():
            problem_options[key] = int(problem_options[key])
        if not isinstance(problem_options[key], int) or problem_options[key] < 0:
            raise TypeError(f'The argument {key} must be a nonnegative integer.')
        if ProblemOptionKey.N_MIN in problem_options and ProblemOptionKey.N_MAX in problem_options and problem_options[ProblemOptionKey.N_MIN] > problem_options[ProblemOptionKey.N_MAX]:
            raise ValueError(f'The argument {ProblemOptionKey.N_MIN.value} must be less than or equal to the argument {ProblemOptionKey.N_MAX.value}.')
        if ProblemOptionKey.M_MIN in problem_options and ProblemOptionKey.M_MAX in problem_options and problem_options[ProblemOptionKey.M_MIN] > problem_options[ProblemOptionKey.M_MAX]:
            raise ValueError(f'The argument {ProblemOptionKey.M_MIN.value} must be less than or equal to the argument {ProblemOptionKey.M_MAX.value}.')

    # Attempt to load the CUTEst problem.
    cutest_problem = None
    logger = get_logger(__name__)
    logger.info(f'Loading CUTEst problem {problem_name}.')
    try:
        cutest_problem = pycutest.import_problem(problem_name)
    except Exception as err:
        logger.error(f'Failed to load CUTEst problem {problem_name}: {err}')

    # If the problem is not successfully loaded or invalid, raise an exception.
    if cutest_problem is None:
        raise ProblemError(f'Failed to load CUTEst problem {problem_name}.')
    elif cutest_problem is not None and not _is_valid(cutest_problem, **problem_options):
        logger.warning(f'CUTEst problem {problem_name} successfully loaded but invalid; it is discarded.')
        raise ProblemError(f'CUTEst problem {problem_name} is invalid.')

    # The problem is successfully loaded and valid. Build the bound, linear, and
    # nonlinear constraints from the CUTEst problem and return.
    xl = np.array(cutest_problem.bl)
    xl[xl <= -1e20] = -np.inf
    xu = np.array(cutest_problem.bu)
    xu[xu >= 1e20] = np.inf
    if cutest_problem.m > 0:
        constraints = {
            'cub': lambda x: _cub(cutest_problem, x),
            'ceq': lambda x: _ceq(cutest_problem, x),
        }
        constraints['aub'], constraints['bub'] = _build_linear_ub(cutest_problem)
        constraints['aeq'], constraints['beq'] = _build_linear_eq(cutest_problem)

        # Count the number of nonlinear constraints.
        idx_ub = ~(cutest_problem.is_linear_cons | cutest_problem.is_eq_cons)
        constraints['m_nonlinear_ub'] = np.count_nonzero(cutest_problem.cl[idx_ub] > -1e20) + np.count_nonzero(cutest_problem.cu[idx_ub] < 1e20)
        constraints['m_nonlinear_eq'] = np.count_nonzero(~cutest_problem.is_linear_cons & cutest_problem.is_eq_cons)
    else:
        constraints = {
            'm_nonlinear_ub': 0,
            'm_nonlinear_eq': 0,
        }
    logger.info(f'CUTEst problem {cutest_problem.name} (n={cutest_problem.n}, m={cutest_problem.m}) successfully loaded.')
    return Problem(cutest_problem.obj, cutest_problem.x0, xl, xu, **constraints)


def find_cutest_problem_names(constraints, **problem_options):
    """
    Find the names of all the CUTEst problems that satisfy given requirements.

    Parameters
    ----------
    constraints : str
        Type of constraints that the CUTEst problems must have. It should
        contain one or more of the following substrings: ``'unconstrained'``,
        ``'fixed'``, ``'bound'``, ``'adjacency'``, ``'linear'``,
        ``'quadratic'``, ``'other'``.

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
    m_min : int, optional
        Minimum number of linear and nonlinear constraints.
    m_max : int, optional
        Maximum number of linear and nonlinear constraints.

    Examples
    --------
    To find all the unconstrained problems with at most 100 variables, use:


    >>> from OptiProfiler import find_cutest_problem_names
    >>>
    >>> problem_names = find_cutest_problem_names('unconstrained', n_max=100)
    """
    import pycutest

    # Check the arguments.
    if not isinstance(constraints, str):
        raise TypeError('The argument constraints must be a string.')
    for constraint in constraints.split():
        if constraint not in ['unconstrained', 'fixed', 'bound', 'adjacency', 'linear', 'quadratic', 'other']:
            raise ValueError(f'Unknown constraint: {constraint}.')
    if ProblemOptionKey.N_MIN in problem_options and isinstance(problem_options[ProblemOptionKey.N_MIN], float) and problem_options[ProblemOptionKey.N_MIN].is_integer():
        problem_options[ProblemOptionKey.N_MIN.value] = int(problem_options[ProblemOptionKey.N_MIN])
    if ProblemOptionKey.N_MIN in problem_options and (not isinstance(problem_options[ProblemOptionKey.N_MIN], int) or problem_options[ProblemOptionKey.N_MIN] < 1):
        raise TypeError(f'The argument {ProblemOptionKey.N_MIN.value} must be a positive integer.')
    if ProblemOptionKey.N_MAX in problem_options and isinstance(problem_options[ProblemOptionKey.N_MAX], float) and problem_options[ProblemOptionKey.N_MAX].is_integer():
        problem_options[ProblemOptionKey.N_MAX.value] = int(problem_options[ProblemOptionKey.N_MAX])
    if ProblemOptionKey.N_MAX in problem_options and (not isinstance(problem_options[ProblemOptionKey.N_MAX], int) or problem_options[ProblemOptionKey.N_MAX] < problem_options.get(ProblemOptionKey.N_MIN, 1)):
        raise TypeError(f'The argument {ProblemOptionKey.N_MAX.value} must be an integer greater than or equal to {ProblemOptionKey.N_MIN.value}.')
    if ProblemOptionKey.M_MIN in problem_options and isinstance(problem_options[ProblemOptionKey.M_MIN], float) and problem_options[ProblemOptionKey.M_MIN].is_integer():
        problem_options[ProblemOptionKey.M_MIN.value] = int(problem_options[ProblemOptionKey.M_MIN])
    if ProblemOptionKey.M_MIN in problem_options and (not isinstance(problem_options[ProblemOptionKey.M_MIN], int) or problem_options[ProblemOptionKey.M_MIN] < 0):
        raise TypeError(f'The argument {ProblemOptionKey.M_MIN.value} must be a nonnegative integer.')
    if ProblemOptionKey.M_MAX in problem_options and isinstance(problem_options[ProblemOptionKey.M_MAX], float) and problem_options[ProblemOptionKey.M_MAX].is_integer():
        problem_options[ProblemOptionKey.M_MAX.value] = int(problem_options[ProblemOptionKey.M_MAX])
    if ProblemOptionKey.M_MAX in problem_options and (not isinstance(problem_options[ProblemOptionKey.M_MAX], int) or problem_options[ProblemOptionKey.M_MAX] < problem_options.get(ProblemOptionKey.M_MIN, 0)):
        raise TypeError(f'The argument {ProblemOptionKey.M_MAX.value} must be an integer greater than or equal to {ProblemOptionKey.M_MIN.value}.')

    # Find all the problems that satisfy the constraints.
    excluded_problem_names = {}
    problem_names = pycutest.find_problems(objective='constant linear quadratic sum of squares other', constraints=constraints, n=[problem_options.get(ProblemOptionKey.N_MIN, 1), problem_options.get(ProblemOptionKey.N_MAX, np.inf)], m=[problem_options.get(ProblemOptionKey.M_MIN, 0), problem_options.get(ProblemOptionKey.M_MAX, np.inf)], userM=False)
    return sorted(set(problem_names).difference(excluded_problem_names))


def _1d_array(x, message):
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


def _2d_array(x, message):
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
