import re
import subprocess

import numpy as np


class Problem:

    def __init__(self, fun, x0, xl=None, xu=None, aub=None, bub=None, aeq=None, beq=None, cub=None, ceq=None, m_nonlinear_ub=None, m_nonlinear_eq=None):
        self._fun = fun
        self._x0 = x0
        self._xl = xl
        self._xu = xu
        self._aub = aub
        self._bub = bub
        self._aeq = aeq
        self._beq = beq
        self._cub = cub
        self._ceq = ceq
        self._m_nonlinear_ub = m_nonlinear_ub
        self._m_nonlinear_eq = m_nonlinear_eq

    @property
    def n(self):
        return self.x0.size

    @property
    def m_linear_ub(self):
        return self.bub.size

    @property
    def m_linear_eq(self):
        return self.beq.size

    @property
    def m_nonlinear_ub(self):
        if self._m_nonlinear_ub is None:
            if self._cub is None:
                self._m_nonlinear_ub = 0
            else:
                raise ValueError('The number of nonlinear inequality constraints is unknown.')
        return self._m_nonlinear_ub

    @property
    def m_nonlinear_eq(self):
        if self._m_nonlinear_eq is None:
            if self._ceq is None:
                self._m_nonlinear_eq = 0
            else:
                raise ValueError('The number of nonlinear equality constraints is unknown.')
        return self._m_nonlinear_eq

    @property
    def x0(self):
        return self._x0

    @property
    def xl(self):
        return self._xl if self._xl is not None else np.full(self.n, -np.inf)

    @property
    def xu(self):
        return self._xu if self._xu is not None else np.full(self.n, np.inf)

    @property
    def aub(self):
        return self._aub if self._aub is not None else np.empty((0, self.n))

    @property
    def bub(self):
        return self._bub if self._bub is not None else np.empty(0)

    @property
    def aeq(self):
        return self._aeq if self._aeq is not None else np.empty((0, self.n))

    @property
    def beq(self):
        return self._beq if self._beq is not None else np.empty(0)

    def fun(self, x):
        return float(self._fun(x))

    def cub(self, x):
        if self._cub is None:
            c = np.empty(0)
        else:
            c = np.array(self._cub(x), dtype=float)
        if self._m_nonlinear_ub is None:
            self._m_nonlinear_ub = c.size
        return c

    def ceq(self, x):
        if self._ceq is None:
            c = np.empty(0)
        else:
            c = np.array(self._ceq(x), dtype=float)
        if self._m_nonlinear_eq is None:
            self._m_nonlinear_eq = c.size
        return c


def load_cutest(problem_name, options=None):
    import pycutest

    # Attempt to load the CUTEst problem.
    cutest_problem = None
    try:
        if pycutest.problem_properties(problem_name)['n'] == 'variable':
            dimensions = _cutest_dimensions(problem_name, options)
            if dimensions.size > 0:
                cutest_problem = pycutest.import_problem(problem_name, sifParams={'N': dimensions[-1]})
        else:
            cutest_problem = pycutest.import_problem(problem_name)
    except Exception as err:
        pass

    # If the problem is successfully loaded and valid, return it.
    if cutest_problem is not None and _is_valid_cutest_problem(cutest_problem, options):
        x0 = np.array(cutest_problem.x0)
        xl = np.array(cutest_problem.bl)
        xl[xl <= -1e20] = -np.inf
        xu = np.array(cutest_problem.bu)
        xu[xu >= 1e20] = np.inf
        return Problem(cutest_problem.obj, x0, xl, xu)


def _cutest_dimensions(problem_name, options):
    import pycutest

    # Get all the SIF parameters.
    command = [pycutest.get_sifdecoder_path(), '-show', problem_name]
    process = subprocess.Popen(command, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    stdout = process.stdout.read()
    process.wait()

    # Extract all the dimensions that are available.
    pattern = re.compile(r'^N=(?P<dim>\d+)')
    dimensions = []
    for line in stdout.split('\n'):
        match = pattern.match(line)
        if match:
            dimensions.append(int(match.group('dim')))

    # Keep only the dimensions that are within the specified range.
    dimensions = np.sort(dimensions)
    if options is not None and 'n_min' in options:
        dimensions = dimensions[dimensions >= options['n_min']]
    if options is not None and 'n_max' in options:
        dimensions = dimensions[dimensions <= options['n_max']]

    return dimensions


def _is_valid_cutest_problem(cutest_problem, options):
    # Check that all the variables are continuous.
    is_valid = np.all(cutest_problem.vartype == 0)

    # Check that the dimensions are within the specified range.
    if options is not None and 'n_min' in options:
        is_valid = is_valid and cutest_problem.n >= options['n_min']
    if options is not None and 'n_max' in options:
        is_valid = is_valid and cutest_problem.n <= options['n_max']
    if options is not None and 'm_min' in options:
        is_valid = is_valid and cutest_problem.m >= options['m_min']
    if options is not None and 'm_max' in options:
        is_valid = is_valid and cutest_problem.m <= options['m_max']

    return is_valid
