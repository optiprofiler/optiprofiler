from cobyqa import minimize
from pdfo import pdfo

from OptiProfiler import find_cutest_problem_names, create_profiles


def cobyqa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    minimize(fun, x0, options={'maxfun': max_eval})


def newuoa(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_eval):
    pdfo(fun, x0, options={'maxfev': max_eval})


if __name__ == '__main__':
    problem_names = find_cutest_problem_names('unconstrained', n_max=5)
    create_profiles([cobyqa, newuoa], ['COBYQA', 'NEWUOA'], problem_names, 'plain', n_max=5)
