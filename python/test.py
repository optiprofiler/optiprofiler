from OptiProfiler import find_cutest_problem_names, create_profiles


def uobyqa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    pdfo(fun, x0, method='uobyqa', options={'maxfev': max_eval})


def newuoa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    pdfo(fun, x0, method='newuoa', options={'maxfev': max_eval})


def bobyqa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    pdfo(fun, x0, method='bobyqa', options={'maxfev': max_eval})


def lincoa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    pdfo(fun, x0, method='lincoa', options={'maxfev': max_eval})


def cobyla(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    pdfo(fun, x0, method='cobyla', options={'maxfev': max_eval})


def cobyqa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from cobyqa import minimize
    minimize(fun, x0, options={'max_eval': max_eval})


if __name__ == '__main__':
    problem_names = find_cutest_problem_names('unconstrained', n_max=5)
    create_profiles([newuoa, cobyqa], ['NEWUOA', 'COBYQA'], problem_names, 'plain', n_max=5)
