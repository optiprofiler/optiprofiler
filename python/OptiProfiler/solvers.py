from contextlib import suppress

import numpy as np
from joblib import delayed

from .features import OptionKey
from .problems import ProblemError, FeaturedProblem, load_cutest
from .utils import get_logger


@delayed
def solve(problem_name, problem_options, solvers, labels, feature, max_eval):
    n_solvers = len(solvers)
    n_runs = feature.options[OptionKey.N_RUNS]
    fun_hist = np.full((n_solvers, n_runs, max_eval), np.nan)
    maxcv_hist = np.full((n_solvers, n_runs, max_eval), np.nan)
    with suppress(ProblemError):
        problem = load_cutest(problem_name, **problem_options)
        logger = get_logger(__name__)
        for i_solver in range(n_solvers):
            for i_run in range(n_runs):
                logger.info(f'Solving {problem_name} with {labels[i_solver]} ({i_run + 1}/{n_runs}).')
                featured_problem = FeaturedProblem(problem, feature)
                with suppress(Exception):
                    solvers[i_solver](lambda x: featured_problem.fun(x, i_run), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, featured_problem.cub, featured_problem.ceq, max_eval)
                n_eval = min(featured_problem.n_eval, max_eval)
                fun_hist[i_solver, i_run, :n_eval] = featured_problem.fun_hist[:n_eval]
                maxcv_hist[i_solver, i_run, :n_eval] = featured_problem.maxcv_hist[:n_eval]
                if n_eval > 0:
                    fun_hist[i_solver, i_run, n_eval:] = fun_hist[i_solver, i_run, n_eval - 1]
                    maxcv_hist[i_solver, i_run, n_eval:] = maxcv_hist[i_solver, i_run, n_eval - 1]
    return fun_hist, maxcv_hist
