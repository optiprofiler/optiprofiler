from OptiProfiler import find_cutest_problem_names, create_profiles


def cobyqa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from cobyqa import minimize
    minimize(fun, x0, options={'max_eval': max_eval})


def lincoa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    pdfo(fun, x0, method='lincoa', options={'maxfev': max_eval})


def newuoa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    from pdfo import pdfo
    pdfo(fun, x0, method='newuoa', options={'maxfev': max_eval})


if __name__ == '__main__':
    problem_names = find_cutest_problem_names('unconstrained', n_max=5)
    create_profiles([cobyqa, newuoa], ['COBYQA', 'NEWUOA'], problem_names, 'plain', n_max=5)
