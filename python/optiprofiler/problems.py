import os
import sys
import warnings
from contextlib import redirect_stdout, redirect_stderr

import numpy as np
from numpy.linalg import lstsq
from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint, minimize

from .features import Feature
from .utils import FeatureName, ProblemOption, FeatureOption, ProblemError, get_logger


class Problem:
    r"""
    Optimization problem to be used in the benchmarking.

    This class provides a uniform interface for providing general optimization
    problems. It is used to supply custom problems to the `benchmark`
    function, and the signature of solvers supplied to the `benchmark`
    function can be

        ``solver(problem) -> array_like, shape (n,)``

    where `problem` is an instance of the class `Problem`.

    Examples
    --------
    Consider the unconstrained problem of minimizing the Rosenbrock function

    .. math::

        f(x) = 100 (x_2 - x_1^2)^2 + (1 - x_1)^2.

    To create an instance of the class `Problem` for this problem, run:

    >>> from optiprofiler.problems import Problem
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

    >>> problem.xl
    array([-inf, -inf])
    >>> problem.xu
    array([inf, inf])

    The optional arguments of the constructor of the class `Problem` can be
    used to specify constraints. For example, to specify that the variables
    must be nonnegative, run:

    >>> problem = Problem(rosen, [0, 0], [0, 0])

    The lower bounds on the variables are now zero:

    >>> problem.xl
    array([0., 0.])
    >>> problem.xu
    array([inf, inf])

    The optional arguments ``aub``, ``bub``, ``aeq``, and ``beq`` can be
    used to specify linear inequality and equality constraints, respectively.
    The optional arguments ``cub`` and ``ceq`` can be used to specify
    nonlinear inequality and equality constraints, respectively. For example,
    to specify that the variables must satisfy the constraint
    :math:`x_1^2 + x_2^2 \le 1` and :math:`x_1^3 - x_2^2 \le 1`, run:

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
    ``m_nonlinear_ub``, it will be inferred at the first call to ``cub``.
    Nonlinear equality constraints can be specified in a similar way, using the
    optional arguments ``ceq`` and ``m_nonlinear_eq``.
    """

    def __init__(self, fun, x0, name=None, xl=None, xu=None, aub=None, bub=None, aeq=None, beq=None, cub=None, ceq=None, grad=None, hess=None, jcub=None, jceq=None, hcub=None, hceq=None):
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
            Coefficient matrix of the linear constraints ``aub @ x <= bub``.
        bub : array_like, shape (m_linear_ub,), optional
            Right-hand side of the linear constraints ``aub @ x <= bub``.
        aeq : array_like, shape (m_linear_eq, n), optional
            Coefficient matrix of the linear constraints ``aeq @ x == beq``.
        beq : array_like, shape (m_linear_eq,), optional
            Right-hand side of the linear constraints ``aeq @ x == beq``.
        cub : callable, optional
            Nonlinear inequality constraint ``cub(x) <= 0``.

                ``cub(x) -> array_like, shape (m_nonlinear_ub,)``

            where ``x`` is an array with shape (n,).
        ceq : callable, optional
            Nonlinear equality constraint ``ceq(x) == 0``.

                ``ceq(x) -> array_like, shape (m_nonlinear_eq,)``

            where ``x`` is an array with shape (n,).

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
            raise TypeError('The argument `fun` for problem must be callable.')

        # Preprocess the initial guess.
        self._x0 = _process_1d_array(x0, 'The argument `x0` for problem must be a one-dimensional array.')

        # Preprocess the name.
        self._name = name
        if self._name is not None and not isinstance(self._name, str):
            raise TypeError('The argument `name` for problem must be a string.')

        # Preprocess the bound constraints.
        self._xl = xl
        if self._xl is not None:
            self._xl = _process_1d_array(self._xl, 'The argument `xl` for problem must be a one-dimensional array.')
        self._xu = xu
        if self._xu is not None:
            self._xu = _process_1d_array(self._xu, 'The argument `xu` for problem must be a one-dimensional array.')

        # Preprocess the linear constraints.
        self._aub = aub
        if self._aub is not None:
            self._aub = _process_2d_array(self._aub, 'The argument `aub` for problem must be a two-dimensional array.')
        self._bub = bub
        if self._bub is not None:
            self._bub = _process_1d_array(self._bub, 'The argument `bub` for problem must be a one-dimensional array.')
        self._aeq = aeq
        if self._aeq is not None:
            self._aeq = _process_2d_array(self._aeq, 'The argument `aeq` for problem must be a two-dimensional array.')
        self._beq = beq
        if self._beq is not None:
            self._beq = _process_1d_array(self._beq, 'The argument `beq` for problem must be a one-dimensional array.')

        # Preprocess the nonlinear constraints.
        self._cub = cub
        self._m_nonlinear_ub = 0
        if self._cub is not None:
            if not callable(self._cub):
                raise TypeError('The argument `cub` for problem must be callable.')
            else:
                try:
                    self._m_nonlinear_ub = self._cub(self._x0).size
                except Exception as err:
                    logger = get_logger(__name__)
                    logger.warning(f'Failed to evaluate the nonlinear inequality constraint function: {err}')
        self._ceq = ceq
        self._m_nonlinear_eq = 0
        if self._ceq is not None:
            if not callable(self._ceq):
                raise TypeError('The argument `ceq` for problem must be callable.')
            else:
                try:
                    self._m_nonlinear_eq = self._ceq(self._x0).size
                except Exception as err:
                    logger = get_logger(__name__)
                    logger.warning(f'Failed to evaluate nonlinear equality constraints: {err}')

        # Preprocess the gradient and the Hessian of the objective function.
        self._grad = grad
        if self._grad is not None:
            if not callable(self._grad):
                raise TypeError('The argument `grad` for problem must be callable.')
        self._hess = hess
        if self._hess is not None:
            if not callable(self._hess):
                raise TypeError('The argument `hess` for problem must be callable.')

        # Preprocess the Jacobian and the Hessian of the nonlinear constraints.
        self._jcub = jcub
        if self._jcub is not None:
            if not callable(self._jcub):
                raise TypeError('The argument `jcub` for problem must be callable.')
        self._jceq = jceq
        if self._jceq is not None:
            if not callable(self._jceq):
                raise TypeError('The argument `jceq` for problem must be callable.')
        self._hcub = hcub
        if self._hcub is not None:
            if not callable(self._hcub):
                raise TypeError('The argument `hcub` for problem must be callable.')
        self._hceq = hceq
        if self._hceq is not None:
            if not callable(self._hceq):
                raise TypeError('The argument `hceq` for problem must be callable.')

        # Check that the arguments are consistent.
        if self.xl.size != self.n:
            raise ValueError(f'The argument `xl` for problem must have size {self.n}.')
        if self.xu.size != self.n:
            raise ValueError(f'The argument `xu` for problem must have size {self.n}.')
        if self.aub.shape != (self.m_linear_ub, self.n):
            raise ValueError(f'The argument `aub` for problem must have shape {(self.m_linear_ub, self.n)}.')
        if self.aeq.shape != (self.m_linear_eq, self.n):
            raise ValueError(f'The argument `aeq` for problem must have shape {(self.m_linear_eq, self.n)}.')
        if self.bub.size != self.m_linear_ub:
            raise ValueError(f'The argument `bub` for problem must have size {self.m_linear_ub}.')
        if self.beq.size != self.m_linear_eq:
            raise ValueError(f'The argument `beq` for problem must have size {self.m_linear_eq}.')

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
    def mb(self):
        """
        Number of bound constraints.

        Returns
        -------
        int
            Number of bound constraints.
        """
        return sum(self.xl > -np.inf) + sum(self.xu < np.inf)

    @property
    def m_linear_ub(self):
        """
        Number of linear inequality constraints.

        Returns
        -------
        int
            Number of linear inequality constraints.
        """
        return sum(~np.isinf(self.bub))

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
        """
        return self._m_nonlinear_ub

    @property
    def m_nonlinear_eq(self):
        """
        Number of nonlinear equality constraints.

        Returns
        -------
        int
            Number of nonlinear equality constraints.
        """
        return self._m_nonlinear_eq

    @property
    def mlcon(self):
        """
        Total number of linear constraints (inequality and equality).

        Returns
        -------
        int
            Total number of linear constraints.
        """
        return self.m_linear_ub + self.m_linear_eq

    @property
    def mnlcon(self):
        """
        Total number of nonlinear constraints (inequality and equality).

        Returns
        -------
        int
            Total number of nonlinear constraints.
        """
        return self.m_nonlinear_ub + self.m_nonlinear_eq

    @property
    def mcon(self):
        """
        Total number of constraints (linear and nonlinear).

        Returns
        -------
        int
            Total number of constraints.
        """
        return self.mlcon + self.mnlcon

    @property
    def ptype(self):
        """
        Type of the problem.

        Returns
        -------
        str
            Type of the problem.
        """
        try:
            if self.mnlcon > 0:
                return 'n'
            elif self.mlcon > 0:
                return 'l'
            elif self.mb > 0:
                return 'b'
            else:
                return 'u'
        except ValueError:
            return 'n'

    @property
    def name(self):
        """
        Name of the problem.

        Returns
        -------
        str
            Name of the problem.
        """
        return self._name if self._name is not None else 'Unnamed Problem'

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
    def xl(self):
        """
        Lower bounds on the variables.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Lower bounds on the variables.
        """
        return np.copy(self._xl) if self._xl is not None else np.full(self.n, -np.inf)

    @property
    def xu(self):
        """
        Upper bounds on the variables.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Upper bounds on the variables.
        """
        return np.copy(self._xu) if self._xu is not None else np.full(self.n, np.inf)

    @property
    def aub(self):
        """
        Coefficient matrix of the linear constraints ``aub @ x <= bub``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_ub, n)
            Coefficient matrix of the linear inequality constraints.
        """
        return np.copy(self._aub) if self._aub is not None else np.empty((0, self.n))

    @property
    def bub(self):
        """
        Right-hand side of the linear constraints ``aub @ x <= bub``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_ub,)
            Right-hand side of the linear inequality constraints.
        """
        return np.copy(self._bub) if self._bub is not None else np.empty(0)

    @property
    def aeq(self):
        """
        Coefficient matrix of the linear constraints ``aeq @ x == beq``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_eq, n)
            Coefficient matrix of the linear equality constraints.
        """
        return np.copy(self._aeq) if self._aeq is not None else np.empty((0, self.n))

    @property
    def beq(self):
        """
        Right-hand side of the linear constraints ``aeq @ x == beq``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_eq,)
            Right-hand side of the linear equality constraints.
        """
        return np.copy(self._beq) if self._beq is not None else np.empty(0)

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
        x = _process_1d_array(x, 'The argument `x` for method `fun` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `fun` in problem must have size {self.n}.')
        try:
            f = self._fun(x)
            f = float(f)
        except Exception as err:
            logger = get_logger(__name__)
            logger.warning(f'Failed to evaluate the objective function: {err}')
            f = np.nan
        return f

    def grad(self, x):
        """
        Evaluate the gradient of the objective function.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the gradient of the objective function.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Gradient of the objective function at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `grad` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `grad` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `grad` in problem must have size {self.n}.')
        if self._grad is None:
            g = np.empty(0)
        else:
            try:
                g = self._grad(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the gradient of the objective function: {err}')
                g = np.full(self.n, np.nan)
            g = _process_1d_array(g, 'The return value of the argument `grad` for problem must be a one-dimensional array.')
            if g.size != self.n:
                raise ValueError(f'The return value of the argument `grad` for problem must have size {self.n}.')
        return g

    def hess(self, x):
        """
        Evaluate the Hessian of the objective function.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the Hessian of the objective function.

        Returns
        -------
        `numpy.ndarray`, shape (n, n)
            Hessian of the objective function at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `hess` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `hess` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `hess` in problem must have size {self.n}.')
        if self._hess is None:
            h = np.empty((0, 0))
        else:
            try:
                h = self._hess(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the Hessian of the objective function: {err}')
                h = np.full((self.n, self.n), np.nan)
            h = _process_2d_array(h, 'The return value of the argument `hess` for problem must be a two-dimensional array.')
            if h.shape != (self.n, self.n):
                raise ValueError(f'The return value of the argument `hess` for problem must have shape {(self.n, self.n)}.')
        return h

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
        x = _process_1d_array(x, 'The argument `x` for method `cub` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `cub` in problem must have size {self.n}.')
        if self._cub is None:
            c = np.empty(0)
        else:
            try:
                c = self._cub(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the nonlinear inequality constraint function: {err}')
                c = np.full(self.m_nonlinear_ub, np.nan)
            c = _process_1d_array(c, 'The return value of the argument `cub` for problem must be a one-dimensional array.')
            if c.size != self.m_nonlinear_ub:
                raise ValueError(f'The return value of the argument `cub` for problem must have size {self.m_nonlinear_ub}.')
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
        x = _process_1d_array(x, 'The argument `x` for method `ceq` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `ceq` in problem must have size {self.n}.')
        if self._ceq is None:
            c = np.empty(0)
        else:
            try:
                c = self._ceq(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the nonlinear equality constraint function: {err}')
                c = np.full(self.m_nonlinear_eq, np.nan)
            c = _process_1d_array(c, 'The return value of the argument `ceq` for problem must be a one-dimensional array.')
            if c.size != self.m_nonlinear_eq:
                raise ValueError(f'The return value of the argument `ceq` for problem must have size {self.m_nonlinear_eq}.')
        return c

    def jcub(self, x):
        """
        Evaluate the Jacobian of the nonlinear inequality constraints.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the Jacobian of the nonlinear inequality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (m_nonlinear_ub, n)
            Jacobian of the nonlinear inequality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `jcub` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `jcub` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `jcub` in problem must have size {self.n}.')
        if self._jcub is None:
            j = np.empty((0, 0))
        else:
            try:
                j = self._jcub(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the Jacobian of the nonlinear inequality constraint function: {err}')
                j = np.full((self.m_nonlinear_ub, self.n), np.nan)
            j = _process_2d_array(j, 'The return value of the argument `jcub` for problem must be a two-dimensional array.')
            if j.shape != (self.m_nonlinear_ub, self.n):
                raise ValueError(f'The return value of the argument `jcub` for problem must have shape {(self.m_nonlinear_ub, self.n)}.')
        return j

    def jceq(self, x):
        """
        Evaluate the Jacobian of the nonlinear equality constraints.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the Jacobian of the nonlinear equality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (m_nonlinear_eq, n)
            Jacobian of the nonlinear equality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `jceq` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `jceq` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `jceq` in problem must have size {self.n}.')
        if self._jceq is None:
            j = np.empty((0, 0))
        else:
            try:
                j = self._jceq(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the Jacobian of the nonlinear equality constraint function: {err}')
                j = np.full((self.m_nonlinear_eq, self.n), np.nan)
            j = _process_2d_array(j, 'The return value of the argument `jceq` for problem must be a two-dimensional array.')
            if j.shape != (self.m_nonlinear_eq, self.n):
                raise ValueError(f'The return value of the argument `jceq` for problem must have shape {(self.m_nonlinear_eq, self.n)}.')
        return j

    def hcub(self, x):
        """
        Evaluate the Hessian of the nonlinear inequality constraints.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the Hessian of the nonlinear inequality constraints.

        Returns
        -------
        `list` of `numpy.ndarray`, shape (m_nonlinear_ub,)
            List of Hessians of the nonlinear inequality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `hcub` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `hcub` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `hcub` in problem must have size {self.n}.')
        if self._hcub is None:
            h = []
        else:
            try:
                h = self._hcub(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the Hessian of the nonlinear inequality constraint function: {err}')
                h = [np.full((self.n, self.n), np.nan)] * self.m_nonlinear_ub
            h = [_process_2d_array(h_i, 'Each element of the return value of the argument `hcub` for problem must be a two-dimensional array.') for h_i in h]
            if len(h) != self.m_nonlinear_ub:
                raise ValueError(f'The return value of the argument `hcub` for problem must have {self.m_nonlinear_ub} elements.')
            for h_i in h:
                if h_i.shape != (self.n, self.n):
                    raise ValueError(f'Each element of the return value of the argument `hcub` for problem must have shape {(self.n, self.n)}.')
        return h

    def hceq(self, x):
        """
        Evaluate the Hessian of the nonlinear equality constraints.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the Hessian of the nonlinear equality constraints.

        Returns
        -------
        `list` of `numpy.ndarray`, shape (m_nonlinear_eq,)
            List of Hessians of the nonlinear equality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `hceq` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `hceq` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `hceq` in problem must have size {self.n}.')
        if self._hceq is None:
            h = []
        else:
            try:
                h = self._hceq(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the Hessian of the nonlinear equality constraint function: {err}')
                h = [np.full((self.n, self.n), np.nan)] * self.m_nonlinear_eq
            h = [_process_2d_array(h_i, 'Each element of the return value of the argument `hceq` for problem must be a two-dimensional array.') for h_i in h]
            if len(h) != self.m_nonlinear_eq:
                raise ValueError(f'The return value of the argument `hceq` for problem must have {self.m_nonlinear_eq} elements.')
            for h_i in h:
                if h_i.shape != (self.n, self.n):
                    raise ValueError(f'Each element of the return value of the argument `hceq` for problem must have shape {(self.n, self.n)}.')
        return h

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
            Maximum constraint violation.
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
        x = _process_1d_array(x, 'The argument `x` for method `maxcv` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `maxcv` in problem must have size {self.n}.')
        cv_bounds = np.max(self.xl - x, initial=0.0)
        cv_bounds = np.max(x - self.xu, initial=cv_bounds)
        cv_linear = np.max(self.aub @ x - self.bub, initial=0.0)
        cv_linear = np.max(np.abs(self.aeq @ x - self.beq), initial=cv_linear)
        cv_nonlinear = np.max(self.cub(x), initial=0.0)
        cv_nonlinear = np.max(np.abs(self.ceq(x)), initial=cv_nonlinear)
        cv = max(cv_bounds, cv_linear, cv_nonlinear)
        return cv, cv_bounds, cv_linear, cv_nonlinear

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
        if self.ptype == 'b':
            self._x0 = np.clip(self._x0, self.xl, self.xu)
        elif self.ptype == 'l' and self.m_linear_ub == 0 and np.all(self.xl == -np.inf) and np.all(self.xu == np.inf):
            self._x0 = lstsq(self.aeq, self.beq - self.aeq @ self.x0)[0]
        elif self.ptype != 'u':
            bounds = Bounds(self.xl, self.xu)
            constraints = []
            if self.m_linear_ub > 0:
                constraints.append(LinearConstraint(self.aub, -np.inf, self.bub))
            if self.m_linear_eq > 0:
                constraints.append(LinearConstraint(self.aeq, self.beq, self.beq))
            if self.m_nonlinear_ub > 0:
                constraints.append(NonlinearConstraint(self.cub, -np.inf, np.zeros(self.m_nonlinear_ub)))
            if self.m_nonlinear_eq > 0:
                constraints.append(NonlinearConstraint(self.ceq, np.zeros(self.m_nonlinear_eq), np.zeros(self.m_nonlinear_eq)))

            def dist_x0_sq(x):
                g = x - self.x0
                return 0.5 * (g @ g), g

            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = minimize(dist_x0_sq, self.x0, jac=True, hessp=lambda p: p, bounds=bounds, constraints=constraints)
            self._x0 = res.x

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

        self._problem = problem
        # Preprocess the feature.
        self._feature = feature
        if not isinstance(self._feature, Feature):
            raise TypeError('The argument `feature` for featured problem must be an instance of the class Feature.')

        # Preprocess the maximum number of function evaluations.
        self._max_eval = max_eval
        if isinstance(self._max_eval, float) and self._max_eval.is_integer():
            self._max_eval = int(self._max_eval)
        if not isinstance(self._max_eval, int):
            raise TypeError('The argument `max_eval` for featured problem must be an integer.')
        if self._max_eval < 1:
            raise ValueError('The argument `max_eval` for featured problem must be positive.')

        # Preprocess the seed.
        self._seed = seed
        if self._seed is not None:
            if isinstance(self._seed, float) and self._seed.is_integer():
                self._seed = int(self._seed)
            if not isinstance(self._seed, int):
                raise TypeError('The argument seed must be an integer.')
            if self._seed < 0:
                raise ValueError('The argument seed must be nonnegative.')

        # Modify the problem according to the feature.
        self._x0 = self._feature.modifier_x0(self._seed, self._problem)
        self._xl, self._xu = self._feature.modifier_bounds(self._seed, self._problem)
        self._aub, self._bub = self._feature.modifier_linear_ub(self._seed, self._problem)
        self._aeq, self._beq = self._feature.modifier_linear_eq(self._seed, self._problem)

        # Store the histories of the objective function values, nonlinear
        # constraints, and maximum constraint violations.
        self._fun_hist = []
        self._cub_hist = []
        self._ceq_hist = []
        self._maxcv_hist = []
        self._last_fun = np.nan
        self._last_cub = np.nan
        self._last_ceq = np.nan

    def __new__(cls, problem, feature, max_eval, seed=None):
        # Preprocess the problem.
        if not isinstance(problem, Problem):
            raise TypeError('The argument `problem` for featured problem must be an instance of the class Problem.')

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
    def cub_hist(self):
        """
        History of nonlinear inequality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (n_eval, m_nonlinear_ub)
            History of nonlinear inequality constraints.
        """
        return np.array(self._cub_hist)

    @property
    def ceq_hist(self):
        """
        History of nonlinear equality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (n_eval, m_nonlinear_eq)
            History of nonlinear equality constraints.
        """
        return np.array(self._ceq_hist)

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
            Value of the objective function at `x` if maximum number of evaluations has not been reached,
            otherwise returns the last evaluated objective function value.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape.
        """
        if self.n_eval >= self._max_eval:
            # If the maximum number of evaluations has been reached, return
            # the last evaluated objective function value.
            return self._last_fun

        # Generate the affine transformation.
        A, b = self._feature.modifier_affine(self._seed, self._problem)[:2]

        # Evaluate the modified the objective function value according to the feature and return the
        # modified value.
        f = self._feature.modifier_fun(A @ x + b, self._seed, self._problem, self.n_eval)
        self._last_fun = f

        # Evaluate the objective function and store the results.
        f_true = super().fun(A @ x + b)

        # If the feature is 'quantized' and the option ``ground_truth'' is set to true, we should
        # set f_true to f.
        if self._feature.name == 'quantized' and self._feature.options[FeatureOption.GROUND_TRUTH]:
            f_true = f

        # We should not store the modified value because the performance of an optimization solver
        # should be measured using the original objective function.
        self._fun_hist.append(f_true)
        try:
            self._maxcv_hist.append(self.maxcv(x))
        except Exception:
            self._maxcv_hist.append(np.nan)

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
            Values of the nonlinear inequality constraints at `x` if maximum number of evaluations has not been reached,
            otherwise returns the last evaluated nonlinear inequality constraints.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `cub` has an invalid shape.
        """
        if len(self._cub_hist) >= self._max_eval:
            # If the maximum number of evaluations has been reached, return
            # the last evaluated nonlinear inequality constraints.
            return self._last_cub

        # Generate the affine transformation.
        A, b = self._feature.modifier_affine(self._seed, self._problem)[:2]

        # Evaluate the modified nonlinear inequality constraints and store the results.
        c = self._feature.modifier_cub(A @ x + b, self._seed, self._problem, len(self._cub_hist))
        self._last_cub = c

        # Evaluate the nonlinear inequality constraints and store the results.
        c_true = super().cub(A @ x + b)

        # If the feature is 'quantized' and the option ``ground_truth'' is set to true, we should
        # set c_true to c.
        if self._feature.name == 'quantized' and self._feature.options[FeatureOption.GROUND_TRUTH]:
            c_true = c
        self._cub_hist.append(c_true)

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
            Values of the nonlinear equality constraints at `x` if maximum number of evaluations has not been reached,
            otherwise returns the last evaluated nonlinear equality constraints.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `ceq` has an invalid shape.
        """
        if len(self._ceq_hist) >= self._max_eval:
            # If the maximum number of evaluations has been reached, return
            # the last evaluated nonlinear equality constraints.
            return self._last_ceq

        # Generate the affine transformation.
        A, b = self._feature.modifier_affine(self._seed, self._problem)[:2]

        # Evaluate the modified nonlinear equality constraints and store the results.
        c = self._feature.modifier_ceq(A @ x + b, self._seed, self._problem, len(self._ceq_hist))
        self._last_ceq = c

        # Evaluate the nonlinear equality constraints and store the results.
        c_true = super().ceq(A @ x + b)

        # If the feature is 'quantized' and the option ``ground_truth'' is set to true, we should
        # set c_true to c.
        if self._feature.name == 'quantized' and self._feature.options[FeatureOption.GROUND_TRUTH]:
            c_true = c
        self._ceq_hist.append(c_true)

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

        # Generate the affine transformation.
        A, b = self._feature.modifier_affine(self._seed, self._problem)[:2]

        # If the Feature is ``quantized'' and the option ``ground_truth'' is set to true, we should
        # use the modified constraint violation.
        if self._feature.name == 'quantized' and self._feature.options[FeatureOption.GROUND_TRUTH]:
            # `maxcv@Problem` is a method of the class `Problem`. By using this, maxcv will use self.cub_ and
            # self.ceq_ instead of self.problem.cub_ and self.problem.ceq_.
            # (We use `x` instead of `A * x + b` because this step is done inside self.ceq_ and self.cub_ methods.)
            return super().maxcv(x)
        else:
            return self._problem.maxcv(A @ x + b)

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
    if x is None or (isinstance(x, np.ndarray) and x.size == 0):
        return x
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
    if x is None or (isinstance(x, np.ndarray) and x.size == 0):
        return x
    x = np.atleast_2d(x).astype(float)
    if x.ndim != 2:
        raise ValueError(message)
    return x
