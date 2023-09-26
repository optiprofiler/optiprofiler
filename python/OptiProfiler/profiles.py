import os
import warnings
from contextlib import redirect_stderr, redirect_stdout, suppress
from enum import Enum

import numpy as np
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator

from .features import Feature, FeatureName, OptionKey as FeatureOptionKey
from .problems import FeaturedProblem, OptionKey as ProblemOptionKey, ProblemError, load_cutest
from .utils import get_logger


class OptionKey(str, Enum):
    N_JOBS = 'n_jobs'


def create_profiles(solvers, labels, problem_names, feature_name, **kwargs):
    logger = get_logger(__name__)

    # Get the different options from the keyword arguments.
    feature_options = {}
    problem_options = {}
    profile_options = {}
    for key, value in kwargs.items():
        if key in FeatureOptionKey.__members__.values():
            feature_options[key] = value
        elif key in ProblemOptionKey.__members__.values():
            problem_options[key] = value
        elif key in OptionKey.__members__.values():
            profile_options[key] = value
        else:
            raise ValueError(f'Unknown option: {key}.')

    # Set the default profile options.
    profile_options.setdefault(OptionKey.N_JOBS.value, -1)

    # Build the feature.
    feature = Feature(feature_name)
    logger.info(f'Starting the computation of the {feature.name} profiles.')

    # Solve the problems.
    max_eval_factor = 500
    fun_values, maxcv_values, problem_names, problem_dimensions = _solve_all(problem_names, problem_options, solvers, labels, feature, max_eval_factor, profile_options)

    # Compute the merit values.
    merit_values = _compute_merit_values(fun_values, maxcv_values)

    # Extract the merit values at the initial points.
    merit_init = np.nanmin(merit_values[..., 0], 1)

    # Determine the least merit value for each problem.
    merit_min = np.nanmin(merit_values, (1, 2, 3))
    if feature.name in [FeatureName.NOISY, FeatureName.TOUGH, FeatureName.TRUNCATED]:
        feature_plain = Feature('plain')
        logger.info(f'Starting the computation of the plain profiles.')
        fun_values_plain, maxcv_values_plain, _, _ = _solve_all(problem_names, problem_options, solvers, labels, feature_plain, max_eval_factor, profile_options)
        merit_values_plain = _compute_merit_values(fun_values_plain, maxcv_values_plain)
        merit_min_plain = np.nanmin(merit_values_plain, (1, 2, 3))
        merit_min = np.minimum(merit_min, merit_min_plain)

    n_problems, n_solvers, n_runs, max_eval = merit_values.shape
    tolerances = np.logspace(-1, -10, 10)
    for i_profile, tolerance in enumerate(tolerances):
        logger.info(f'Creating performance and data profiles for tolerance {tolerance}.')

        work = np.full((n_problems, n_solvers, n_runs), np.nan)
        for i_problem in range(n_problems):
            for i_solver in range(n_solvers):
                for i_run in range(n_runs):
                    if np.isfinite(merit_min[i_problem]):
                        threshold = max(tolerance * merit_init[i_problem, i_run] + (1.0 - tolerance) * merit_min[i_problem], merit_min[i_problem])
                    else:
                        threshold = -np.inf
                    if np.nanmin(merit_values[i_problem, i_solver, i_run, :]) <= threshold:
                        work[i_problem, i_solver, i_run] = np.argmax(merit_values[i_problem, i_solver, i_run, :] <= threshold) + 1

        # Calculate the x-axes of the performance profiles.
        x_perf = np.full((n_runs, n_problems, n_solvers), np.nan)
        for i_run in range(n_runs):
            for i_problem in range(n_problems):
                if not np.all(np.isnan(work[i_problem, :, i_run])):
                    x_perf[i_run, i_problem, :] = work[i_problem, :, i_run] / np.nanmin(work[i_problem, :, i_run])
        perf_ratio_max = np.nanmax(x_perf, initial=2.0 ** np.finfo(float).eps)
        x_perf[np.isnan(x_perf)] = 2.0 * perf_ratio_max
        x_perf = np.sort(x_perf, 1)
        x_perf = np.reshape(x_perf, (n_problems * n_runs, n_solvers))
        sort_perf = np.argsort(x_perf, 0, 'stable')
        x_perf = np.take_along_axis(x_perf, sort_perf, 0)

        # Calculate the y-axes of the performance profiles.
        y_perf = np.zeros((n_problems * n_runs, n_solvers))
        for i_run in range(n_runs):
            for i_solver in range(n_solvers):
                y = np.full(n_problems * n_runs, np.nan)
                y[i_run * n_problems:(i_run + 1) * n_problems] = np.linspace(1 / n_problems, 1.0, n_problems)
                y = y[sort_perf[:, i_solver]]
                for i_problem in range(n_problems * n_runs):
                    if np.isnan(y[i_problem]):
                        y[i_problem] = y[i_problem - 1] if i_problem > 0 else 0.0
                y_perf[:, i_solver] += y
        y_perf /= n_runs

        # Calculate the x-axes of the data profiles.
        x_data = np.full((n_runs, n_problems, n_solvers), np.nan)
        for i_run in range(n_runs):
            for i_problem in range(n_problems):
                if not np.all(np.isnan(work[i_problem, :, i_run])):
                    x_data[i_run, i_problem, :] = work[i_problem, :, i_run] / (problem_dimensions[i_problem] + 1)
        data_ratio_max = np.nanmax(x_data, initial=np.finfo(float).eps)
        x_data[np.isnan(x_data)] = 2.0 * data_ratio_max
        x_data = np.sort(x_data, 1)
        x_data = np.reshape(x_data, (n_problems * n_runs, n_solvers))
        sort_data = np.argsort(x_data, 0, 'stable')
        x_data = np.take_along_axis(x_data, sort_data, 0)

        # Calculate the y-axes of the performance profiles.
        y_data = np.zeros((n_problems * n_runs, n_solvers))
        for i_run in range(n_runs):
            for i_solver in range(n_solvers):
                y = np.full(n_problems * n_runs, np.nan)
                y[i_run * n_problems:(i_run + 1) * n_problems] = np.linspace(1 / n_problems, 1.0, n_problems)
                y = y[sort_data[:, i_solver]]
                for i_problem in range(n_problems * n_runs):
                    if np.isnan(y[i_problem]):
                        y[i_problem] = y[i_problem - 1] if i_problem > 0 else 0.0
                y_data[:, i_solver] += y
        y_data /= n_runs

        # Plot the performance profiles.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(MaxNLocator(5, prune='lower'))
        ax.yaxis.set_minor_locator(MaxNLocator(10))
        ax.tick_params(direction='in', which='both')
        for j in range(n_solvers):
            x = np.repeat(x_perf[:, j], 2)[1:]
            x = np.r_[0.0, x[0], x, 2.0 * perf_ratio_max]
            y = np.repeat(y_perf[:, j], 2)[:-1]
            y = np.r_[0.0, 0.0, y, y[-1]]
            plt.semilogx(x, y, base=2, label=labels[j])
        plt.xlim(1.0, 1.1 * perf_ratio_max)
        plt.ylim(0.0, 1.0)
        plt.xlabel('Performance ratio')
        plt.ylabel('Performance profiles')
        plt.legend(loc='lower right')
        plt.show()

        # Plot the data profiles.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(MaxNLocator(5, prune='lower'))
        ax.yaxis.set_minor_locator(MaxNLocator(10))
        ax.tick_params(direction='in', which='both')
        for j in range(n_solvers):
            x = np.repeat(x_data[:, j], 2)[1:]
            x = np.r_[0.0, x[0], x, 2.0 * data_ratio_max]
            y = np.repeat(y_data[:, j], 2)[:-1]
            y = np.r_[0.0, 0.0, y, y[-1]]
            plt.plot(x, y, label=labels[j])
        plt.xlim(0.0, 1.1 * data_ratio_max)
        plt.ylim(0.0, 1.0)
        plt.xlabel('Number of simplex gradients')
        plt.ylabel('Data profiles')
        plt.legend(loc='lower right')
        plt.close()


def _solve_all(problem_names, problem_options, solvers, labels, feature, max_eval_factor, profile_options):
    # Solve all problems.
    results = Parallel(n_jobs=profile_options[OptionKey.N_JOBS])(_solve_one(problem_name, problem_options, solvers, labels, feature, max_eval_factor) for problem_name in problem_names)
    _fun_values, _maxcv_values, problem_names, problem_dimensions = zip(*[result for result in results if result is not None])

    # Build
    n_runs = feature.options[FeatureOptionKey.N_RUNS]
    max_eval = max_eval_factor * max(problem_dimensions)
    fun_values = np.full((len(_fun_values), len(solvers), n_runs, max_eval), np.nan)
    maxcv_values = np.full((len(_maxcv_values), len(solvers), n_runs, max_eval), np.nan)
    for i_problem, (fun_value, maxcv_value) in enumerate(zip(_fun_values, _maxcv_values)):
        n_eval = fun_value.shape[-1]
        fun_values[i_problem, ..., :n_eval] = fun_value
        maxcv_values[i_problem, ..., :n_eval] = maxcv_value
        if n_eval > 0:
            fun_values[i_problem, ..., n_eval:] = fun_values[i_problem, ..., n_eval - 1, np.newaxis]
            maxcv_values[i_problem, ..., n_eval:] = maxcv_values[i_problem, ..., n_eval - 1, np.newaxis]
    return fun_values, maxcv_values, problem_names, problem_dimensions


@delayed
def _solve_one(problem_name, problem_options, solvers, labels, feature, max_eval_factor):
    # Load the problem and return if it cannot be loaded.
    try:
        problem = load_cutest(problem_name, **problem_options)
    except ProblemError:
        return

    # Solve the problem with each solver.
    n_solvers = len(solvers)
    n_runs = feature.options[FeatureOptionKey.N_RUNS]
    max_eval = max_eval_factor * problem.n
    fun_values = np.full((n_solvers, n_runs, max_eval), np.nan)
    maxcv_values = np.full((n_solvers, n_runs, max_eval), np.nan)
    logger = get_logger(__name__)
    for i_solver in range(n_solvers):
        for i_run in range(n_runs):
            logger.info(f'Solving {problem_name} with {labels[i_solver]} (run {i_run + 1}/{n_runs}).')
            featured_problem = FeaturedProblem(problem, feature)
            with open(os.devnull, 'w') as devnull:
                with suppress(Exception), warnings.catch_warnings(), redirect_stdout(devnull), redirect_stderr(devnull):
                    warnings.filterwarnings('ignore')
                    solvers[i_solver](lambda x: featured_problem.fun(x, i_run), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, featured_problem.cub, featured_problem.ceq, max_eval)
            n_eval = min(featured_problem.n_eval, max_eval)
            fun_values[i_solver, i_run, :n_eval] = featured_problem.fun_values[:n_eval]
            maxcv_values[i_solver, i_run, :n_eval] = featured_problem.maxcv_values[:n_eval]
            if n_eval > 0:
                fun_values[i_solver, i_run, n_eval:] = fun_values[i_solver, i_run, n_eval - 1]
                maxcv_values[i_solver, i_run, n_eval:] = maxcv_values[i_solver, i_run, n_eval - 1]
    return fun_values, maxcv_values, problem_name, problem.n


def _compute_merit_values(fun_values, maxcv_values):
    is_nearly_feasible = maxcv_values <= 1e-12
    is_very_infeasible = maxcv_values >= 1e-6
    is_undecided = ~is_nearly_feasible & ~is_very_infeasible
    merit_values = np.empty(fun_values.shape)
    merit_values[is_nearly_feasible] = fun_values[is_nearly_feasible]
    merit_values[is_very_infeasible] = np.inf
    merit_values[is_undecided] = fun_values[is_undecided] + 1e8 * maxcv_values[is_undecided]
    return merit_values
