import numpy as np

from OptiProfiler import find_cutest, create_profiles


def uobyqa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    pdfo(fun, x0, method='uobyqa', options={'maxfev': max_eval})


def newuoa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    pdfo(fun, x0, method='newuoa', options={'maxfev': max_eval})


def bobyqa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    from scipy.optimize import Bounds
    bounds = Bounds(xl, xu)
    pdfo(fun, x0, method='bobyqa', bounds=bounds, options={'maxfev': max_eval})


def lincoa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    from scipy.optimize import Bounds, LinearConstraint
    bounds = Bounds(xl, xu)
    constraints = [
        LinearConstraint(aeq, beq, beq),
        LinearConstraint(aub, -np.inf, bub),
    ]
    pdfo(fun, x0, method='lincoa', bounds=bounds, constraints=constraints, options={'maxfev': max_eval})


def cobyla(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    from scipy.optimize import Bounds, LinearConstraint
    bounds = Bounds(xl, xu)
    constraints = [
        LinearConstraint(aeq, beq, beq),
        LinearConstraint(aub, -np.inf, bub),
    ]
    pdfo(fun, x0, method='cobyla', bounds=bounds, constraints=constraints, options={'maxfev': max_eval})


def cobyqa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from cobyqa import minimize
    minimize(fun, x0, options={'max_eval': max_eval})


if __name__ == '__main__':
    problem_names = find_cutest('unconstrained', n_max=5)
    create_profiles([uobyqa, newuoa, bobyqa, lincoa, cobyla], ['UOBYQA', 'NEWUOA', 'BOBYQA', 'LINCOA', 'COBYLA'], problem_names, 'plain', n_max=5)

    problem_names = find_cutest('bound', n_max=5)
    create_profiles([bobyqa, lincoa, cobyla], ['BOBYQA', 'LINCOA', 'COBYLA'], problem_names, 'plain', n_max=5)

    problem_names = find_cutest('adjacency linear', n_max=5)
    create_profiles([lincoa, cobyla], ['LINCOA', 'COBYLA'], problem_names, 'plain', n_max=5)
