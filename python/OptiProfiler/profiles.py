import numpy as np
from joblib import Parallel

from .features import Feature, OptionKey as FeatureOptionKey
from .problems import OptionKey as ProblemOptionKey
from .solvers import solve


def create_profiles(solvers, labels, problem_names, feature_name, **kwargs):
    # Get the feature and problem options from the keyword arguments.
    feature_options = {}
    problem_options = {}
    for key, value in kwargs.items():
        if key in FeatureOptionKey.__members__.values():
            feature_options[key] = value
        elif key in ProblemOptionKey.__members__.values():
            problem_options[key] = value
        else:
            raise ValueError(f'Unknown option: {key}.')

    # Build the feature.
    feature = Feature(feature_name)

    # Solve the problems.
    # TODO: Find a way to determine the maximum number of function evaluations.
    n_problems = len(problem_names)
    n_solvers = len(solvers)
    n_runs = feature.options[FeatureOptionKey.N_RUNS]
    max_eval = 500
    fun_hist, maxcv_hist = zip(*Parallel(n_jobs=-1)(solve(problem_name, problem_options, solvers, labels, feature, max_eval) for problem_name in problem_names))
    fun_hist = np.array(fun_hist)
    maxcv_hist = np.array(maxcv_hist)

    # Compute the merit values.
    is_nearly_feasible = maxcv_hist <= 1e-12
    is_very_infeasible = maxcv_hist >= 1e-6
    is_undecided = ~is_nearly_feasible & ~is_very_infeasible
    merit_values = np.empty((n_problems, n_solvers, n_runs, max_eval))
    merit_values[is_nearly_feasible] = fun_hist[is_nearly_feasible]
    merit_values[is_very_infeasible] = np.inf
    merit_values[is_undecided] = fun_hist[is_undecided] + 1e8 * maxcv_hist[is_undecided]
    print(merit_values)
