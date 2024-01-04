import numpy as np

from OptiProfiler import find_cutest, create_profiles


def uobyqa(fun, x0, max_eval):
    from pdfo import pdfo
    pdfo(fun, x0, method='uobyqa', options={'maxfev': max_eval})


def newuoa(fun, x0, max_eval):
    from pdfo import pdfo
    pdfo(fun, x0, method='newuoa', options={'maxfev': max_eval})


def bobyqa(fun, x0, xl, xu, max_eval):
    from pdfo import pdfo
    from scipy.optimize import Bounds
    bounds = Bounds(xl, xu)
    pdfo(fun, x0, method='bobyqa', bounds=bounds, options={'maxfev': max_eval})


def lincoa(fun, x0, xl, xu, aub, bub, aeq, beq, max_eval):
    from pdfo import pdfo
    from scipy.optimize import Bounds, LinearConstraint
    bounds = Bounds(xl, xu)
    constraints = []
    if bub.size > 0:
        constraints.append(LinearConstraint(aub, -np.inf, bub))
    if beq.size > 0:
        constraints.append(LinearConstraint(aeq, beq, beq))
    pdfo(fun, x0, method='lincoa', bounds=bounds, constraints=constraints, options={'maxfev': max_eval})


def cobyla(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    from scipy.optimize import Bounds, LinearConstraint
    bounds = Bounds(xl, xu)
    constraints = []
    if bub.size > 0:
        constraints.append(LinearConstraint(aub, -np.inf, bub))
    if beq.size > 0:
        constraints.append(LinearConstraint(aeq, beq, beq))
    cub_x0 = cub(x0)
    if cub_x0.size > 0:
        constraints.append(LinearConstraint(cub, -np.inf, np.zeros_like(cub_x0)))
    ceq_x0 = ceq(x0)
    if ceq_x0.size > 0:
        constraints.append(LinearConstraint(ceq, np.zeros_like(ceq_x0), np.zeros_like(ceq_x0)))
    pdfo(fun, x0, method='cobyla', bounds=bounds, constraints=constraints, options={'maxfev': max_eval})


def cobyqa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from cobyqa import minimize
    minimize(fun, x0, xl=xl, xu=xu, aub=aub, bub=bub, aeq=aeq, beq=beq, cub=cub, ceq=ceq, options={'max_eval': max_eval})


if __name__ == '__main__':
    problem_names = find_cutest('unconstrained', n_max=5)
    create_profiles([newuoa, cobyqa], ['NEWUOA', 'COBYQA'], problem_names, 'randomize_x0', n_max=5)
