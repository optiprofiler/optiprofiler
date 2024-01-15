import numpy as np

from OptiProfiler import find_cutest_problems, create_profiles


def uobyqa(fun, x0):
    from pdfo import pdfo

    pdfo(fun, x0, method='uobyqa')


def newuoa(fun, x0):
    from pdfo import pdfo

    pdfo(fun, x0, method='newuoa')


def bobyqa(fun, x0, xl, xu):
    from pdfo import pdfo
    from scipy.optimize import Bounds

    bounds = Bounds(xl, xu)
    pdfo(fun, x0, method='bobyqa', bounds=bounds)


def lincoa(fun, x0, xl, xu, aub, bub, aeq, beq):
    from pdfo import pdfo
    from scipy.optimize import Bounds, LinearConstraint

    bounds = Bounds(xl, xu)
    constraints = []
    if bub.size > 0:
        constraints.append(LinearConstraint(aub, -np.inf, bub))
    if beq.size > 0:
        constraints.append(LinearConstraint(aeq, beq, beq))
    pdfo(fun, x0, method='lincoa', bounds=bounds, constraints=constraints)


def cobyla(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq):
    from pdfo import pdfo
    from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint

    bounds = Bounds(xl, xu)
    constraints = []
    if bub.size > 0:
        constraints.append(LinearConstraint(aub, -np.inf, bub))
    if beq.size > 0:
        constraints.append(LinearConstraint(aeq, beq, beq))
    cub_x0 = cub(x0)
    if cub_x0.size > 0:
        constraints.append(NonlinearConstraint(cub, -np.inf, np.zeros_like(cub_x0)))
    ceq_x0 = ceq(x0)
    if ceq_x0.size > 0:
        constraints.append(NonlinearConstraint(ceq, np.zeros_like(ceq_x0), np.zeros_like(ceq_x0)))
    pdfo(fun, x0, method='cobyla', bounds=bounds, constraints=constraints)


def cobyqa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq):
    from cobyqa import minimize
    from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint

    bounds = Bounds(xl, xu)
    constraints = []
    if bub.size > 0:
        constraints.append(LinearConstraint(aub, -np.inf, bub))
    if beq.size > 0:
        constraints.append(LinearConstraint(aeq, beq, beq))
    cub_x0 = cub(x0)
    if cub_x0.size > 0:
        constraints.append(NonlinearConstraint(cub, -np.inf, np.zeros_like(cub_x0)))
    ceq_x0 = ceq(x0)
    if ceq_x0.size > 0:
        constraints.append(NonlinearConstraint(ceq, np.zeros_like(ceq_x0), np.zeros_like(ceq_x0)))

    minimize(fun, x0, bounds=bounds, constraints=constraints)


if __name__ == '__main__':
    cutest_problem_names = find_cutest_problems('unconstrained', n_max=5)
    create_profiles([newuoa, cobyqa], ['NEWUOA', 'COBYQA'], cutest_problem_names, n_max=5, subfolder='unconstrained')
