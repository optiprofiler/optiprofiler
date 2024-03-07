import numpy as np

from optiprofiler import set_cutest_problem_options, find_cutest_problems, run_benchmark


def uobyqa(fun, x0):
    from pdfo import pdfo

    res = pdfo(fun, x0, method='uobyqa')
    return res.x


def newuoa(fun, x0):
    from pdfo import pdfo

    res = pdfo(fun, x0, method='newuoa')
    return res.x


def bobyqa(fun, x0, lb, ub):
    from pdfo import pdfo
    from scipy.optimize import Bounds

    bounds = Bounds(lb, ub)
    res = pdfo(fun, x0, method='bobyqa', bounds=bounds)
    return res.x


def lincoa(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq):
    from pdfo import pdfo
    from scipy.optimize import Bounds, LinearConstraint

    bounds = Bounds(lb, ub)
    constraints = []
    if b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    res = pdfo(fun, x0, method='lincoa', bounds=bounds, constraints=constraints)
    return res.x


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
    res = pdfo(fun, x0, method='cobyla', bounds=bounds, constraints=constraints)
    return res.x


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
    res = minimize(fun, x0, bounds=bounds, constraints=constraints)
    return res.x


if __name__ == '__main__':
    set_cutest_problem_options(n_max=2)
    cutest_problem_names = find_cutest_problems('unconstrained')
    run_benchmark([cobyqa, newuoa], ['COBYQA', 'NEWUOA'], cutest_problem_names)
