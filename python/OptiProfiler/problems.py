import re
import subprocess

import numpy as np

from .utils import get_logger


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
            x = np.array(x, dtype=float)
            c = np.array(self._cub(x), dtype=float)
        if self._m_nonlinear_ub is None:
            self._m_nonlinear_ub = c.size
        return c

    def ceq(self, x):
        if self._ceq is None:
            c = np.empty(0)
        else:
            x = np.array(x, dtype=float)
            c = np.array(self._ceq(x), dtype=float)
        if self._m_nonlinear_eq is None:
            self._m_nonlinear_eq = c.size
        return c


class ProblemError(Exception):
    pass


def load_cutest(problem_name, options=None):
    import pycutest

    def _dimensions(problem_name, options):
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

    def _is_valid(cutest_problem, options):
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

    def _build_linear_ub(cutest_problem):
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
        idx_eq = cutest_problem.is_linear_cons & cutest_problem.is_eq_cons
        aeq = []
        beq = []
        for index in np.flatnonzero(idx_eq):
            c_val, g_val = cutest_problem.cons(np.zeros(cutest_problem.n), index, True)
            aeq.append(g_val)
            beq.append(c_val - 0.5 * (cutest_problem.cl[index] + cutest_problem.cu[index]))
        return np.reshape(aeq, (-1, cutest_problem.n)), np.array(beq)

    def _cub(cutest_problem, x):
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
        idx_eq = ~cutest_problem.is_linear_cons & cutest_problem.is_eq_cons
        c = []
        for index in np.flatnonzero(idx_eq):
            c_val = cutest_problem.cons(x, index)
            c.append(c_val - 0.5 * (cutest_problem.cl[index] + cutest_problem.cu[index]))
        return np.array(c)

    # Attempt to load the CUTEst problem.
    cutest_problem = None
    logger = get_logger(__name__)
    logger.info(f'Loading CUTEst problem {problem_name}.')
    try:
        if pycutest.problem_properties(problem_name)['n'] == 'variable':
            dimensions = _dimensions(problem_name, options)
            if dimensions.size > 0:
                logger.info(f'Loading CUTEst problem {problem_name} with N={dimensions[-1]}.')
                cutest_problem = pycutest.import_problem(problem_name, sifParams={'N': dimensions[-1]})
            else:
                logger.info(f'No valid dimensions found for CUTEst problem {problem_name}.')
        else:
            cutest_problem = pycutest.import_problem(problem_name)
    except Exception as err:
        logger.warning(f'Failed to load CUTEst problem {problem_name}: {err}')

    # If the problem is not successfully loaded or invalid, raise an exception.
    if cutest_problem is None:
        raise ProblemError(f'Failed to load CUTEst problem {problem_name}.')
    elif cutest_problem is not None and not _is_valid(cutest_problem, options):
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
