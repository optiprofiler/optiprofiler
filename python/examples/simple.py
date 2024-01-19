import numpy as np

from optiprofiler import set_cutest_problem_options, find_cutest_problems, create_profiles


def uobyqa(fun, x0):
    from pdfo import pdfo

    pdfo(fun, x0, method='uobyqa')


def newuoa(fun, x0):
    from pdfo import pdfo

    pdfo(fun, x0, method='newuoa')


def bobyqa(fun, x0, lb, ub):
    from pdfo import pdfo
    from scipy.optimize import Bounds

    bounds = Bounds(lb, ub)
    pdfo(fun, x0, method='bobyqa', bounds=bounds)


def lincoa(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq):
    from pdfo import pdfo
    from scipy.optimize import Bounds, LinearConstraint

    bounds = Bounds(lb, ub)
    constraints = []
    if b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    pdfo(fun, x0, method='lincoa', bounds=bounds, constraints=constraints)


def cobyla(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq, c_ub, c_eq):
    from pdfo import pdfo
    from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint

    bounds = Bounds(lb, ub)
    constraints = []
    if b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    c_ub_x0 = c_ub(x0)
    if c_ub_x0.size > 0:
        constraints.append(NonlinearConstraint(c_ub, -np.inf, np.zeros_like(c_ub_x0)))
    c_eq_x0 = c_eq(x0)
    if c_eq_x0.size > 0:
        constraints.append(NonlinearConstraint(c_eq, np.zeros_like(c_eq_x0), np.zeros_like(c_eq_x0)))
    pdfo(fun, x0, method='cobyla', bounds=bounds, constraints=constraints)


def cobyqa(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq, c_ub, c_eq):
    from cobyqa import minimize
    from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint

    bounds = Bounds(lb, ub)
    constraints = []
    if b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    c_ub_x0 = c_ub(x0)
    if c_ub_x0.size > 0:
        constraints.append(NonlinearConstraint(c_ub, -np.inf, np.zeros_like(c_ub_x0)))
    c_eq_x0 = c_eq(x0)
    if c_eq_x0.size > 0:
        constraints.append(NonlinearConstraint(c_eq, np.zeros_like(c_eq_x0), np.zeros_like(c_eq_x0)))
    minimize(fun, x0, bounds=bounds, constraints=constraints)


if __name__ == '__main__':
    set_cutest_problem_options(n_max=2)
    cutest_problem_names = find_cutest_problems('unconstrained')
    create_profiles([newuoa, uobyqa], ['NEWUOA', 'UOBYQA'], cutest_problem_names, feature_name='plain', benchmark_id='unconstrained')
