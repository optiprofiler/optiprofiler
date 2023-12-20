from OptiProfiler import find_cutest, create_profiles


def newuoa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    pdfo(fun, x0, method='newuoa', options={'maxfev': max_eval})


def cobyqa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from cobyqa import minimize
    minimize(fun, x0, options={'max_eval': max_eval})


if __name__ == '__main__':
    problem_names = find_cutest('unconstrained', n_max=5)[:30]
    create_profiles([newuoa, cobyqa], ['NEWUOA', 'COBYQA'], problem_names, 'randomize_x0', n_max=5)
