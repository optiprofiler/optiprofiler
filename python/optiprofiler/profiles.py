from __future__ import annotations

import copy
import os
import sys
import re
import shutil
import inspect
import logging
import importlib.util
import time
import warnings
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime
from inspect import signature
from multiprocessing import Pool
from multiprocessing.reduction import ForkingPickler
from pathlib import Path
from typing import Any, Callable

import numpy as np
from cycler import cycler
from matplotlib import pyplot as plt
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D
from matplotlib.backends import backend_pdf
from matplotlib.ticker import MaxNLocator, FuncFormatter

from .opclasses import Feature, Problem, FeaturedProblem
from .utils import DEFAULT_LOG_LINE_WIDTH, FeatureName, ProfileOption, FeatureOption, ProblemOption, get_logger, print_log_message, setup_main_process_logging, setup_worker_logging, shorten_log_message, format_log_prefix
from .loader import load_results, save_results_to_h5, save_options
from .profile_utils import check_validity_problem_options, check_validity_profile_options, get_default_problem_options, get_default_profile_options, compute_merit_values, create_stamp, merge_pdfs_with_pypdf, write_report, process_results, init_readme, add_to_readme, compute_scores, _custom_problem_library_dir
from .plotting import draw_hist, set_profile_context, format_float_scientific_latex, draw_profiles, summary_legend_extra_width, latex_escape_text, format_profile_text


def _shorten_log_message(message: object, max_length: int = 180) -> str:
    """
    Shorten an exception message for log output.
    """
    return shorten_log_message(message, max_length)


_SOLVER_LOG_NAMES_KEY = '_solver_log_names'


def _standard_constrained_result_log_length(problem_name_width: int, solver_name_width: int) -> int:
    """
    Length of the standard constrained ``Output result`` log line.
    """
    problem_name = 'P' * max(1, problem_name_width)
    solver_name = 'S' * max(1, solver_name_width)
    line = (
        f'{format_log_prefix("INFO")}'
        f'Output result for {problem_name:<{problem_name_width}} with {solver_name:<{solver_name_width}} '
        f'(run {1:2d}/{1:2d}): f = {-9.5:10.4e}, maxcv = {1.1102e-16:10.4e}.'
    )
    return len(line)


def _get_solver_log_names(solver_names, problem_name_width: int, line_width: int = DEFAULT_LOG_LINE_WIDTH):
    """
    Return short solver display names for aligned terminal logs when needed.
    """
    solver_names = list(solver_names)
    solver_width = max(len(name) for name in solver_names)
    if _standard_constrained_result_log_length(problem_name_width, solver_width) <= line_width:
        return solver_names, False
    return [f'solver{i + 1}' for i in range(len(solver_names))], True


def _with_solver_log_names(profile_options, problem_name_width: int):
    """
    Copy profile options and add internal solver names used only in logs.
    """
    solver_log_names, solver_aliases_used = _get_solver_log_names(
        profile_options[ProfileOption.SOLVER_NAMES],
        problem_name_width,
    )
    profile_options_log = profile_options.copy()
    profile_options_log[_SOLVER_LOG_NAMES_KEY] = solver_log_names
    return profile_options_log, solver_aliases_used


def _log_solver_aliases(logger, solver_names, solver_log_names):
    """
    Print the mapping from short log aliases to true solver names.
    """
    if list(solver_log_names) == list(solver_names):
        return
    logger.info('Solver aliases used in the log:')
    max_alias_len = max(len(name) for name in solver_log_names)
    for alias, name in zip(solver_log_names, solver_names):
        logger.info(f'{alias:<{max_alias_len}} = {name}')


def benchmark(
    solvers: list[Callable[..., Any]] | None = None,
    /,
    **kwargs
) -> tuple[np.ndarray, np.ndarray | None, list[dict] | None]:
    """
    Benchmark optimization solvers on a set of problems with specified features.

    This function creates multiple profiles for benchmarking optimization solvers on a
    set of problems with different features. It generates performance profiles, data
    profiles, and log-ratio profiles [1]_, [2]_, [4]_, [5]_ for the given solvers on
    various test suites, returning solver scores based on the profiles.

    Parameters
    ----------
    solvers : list of callable, optional if ``load`` is provided
        Solvers to benchmark. Each solver must be a callable accepting
        corresponding arguments depending on the test suite you choose:

        - for an unconstrained problem,
          ``solver(fun, x0) -> numpy.ndarray, shape (n,)``,
          where ``fun`` is the objective function accepting a 1-D array and
          returning a float, and ``x0`` is the initial guess (1-D array);
        - for a bound-constrained problem,
          ``solver(fun, x0, xl, xu) -> numpy.ndarray, shape (n,)``,
          where ``xl`` and ``xu`` are the lower and upper bounds (1-D arrays,
          may contain ``-numpy.inf`` or ``numpy.inf``);
        - for a linearly constrained problem,
          ``solver(fun, x0, xl, xu, aub, bub, aeq, beq) -> numpy.ndarray, shape (n,)``,
          where ``aub`` and ``aeq`` are the coefficient matrices of the linear
          inequality and equality constraints, and ``bub`` and ``beq`` are the
          right-hand side vectors;
        - for a nonlinearly constrained problem,
          ``solver(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq) -> numpy.ndarray, shape (n,)``,
          where ``cub`` and ``ceq`` are the nonlinear inequality and equality
          constraint functions accepting a 1-D array and returning a 1-D array.

        All vectors and matrices mentioned above are `numpy.ndarray`.

        If the ``load`` option is provided, solvers can be None,
        in which case data from a previous experiment will be loaded to generate profiles.

    Other Parameters
    ----------------
    Feature options:

    feature_name : str, optional
        Name of the feature to apply to problems. The available features are
        'plain', 'perturbed_x0', 'noisy', 'truncated', 'permuted',
        'linearly_transformed', 'random_nan', 'unrelaxable_constraints',
        'nonquantifiable_constraints', 'quantized', and 'custom'. Default is
        'plain'.
    n_runs : int, optional
        The number of runs of the experiments with the given feature.
        Default is 5 for stochastic features and 1 for deterministic
        features.
    distribution : str or callable, optional
        The distribution of perturbation in 'perturbed_x0'
        feature or random noise in 'noisy' feature. It should be either a
        str (or char), or a callable
        ``(random_stream, dimension) -> random vector``,
        accepting a random_stream and the dimension of a problem and
        returning a random vector with the given dimension. In 'perturbed_x0'
        case, the str should be either 'spherical' or 'gaussian' (default is
        'spherical'). In 'noisy' case, the str should be either 'gaussian'
        or 'uniform' (default is 'gaussian'), and the callable should accept
        a random stream and output size.
    perturbation_level : float, optional
        The magnitude of the perturbation to the initial
        guess in the 'perturbed_x0' feature. Default is 1e-3.
    noise_level : float, optional
        The magnitude of the noise in the 'noisy' feature.
        Default is 1e-3.
    noise_type : str, optional
        The type of the noise in the 'noisy' features. It should
        be either 'absolute', 'relative', or 'mixed'. Default is 'mixed'.
    noise_mode : str, optional
        The mode of the noise in the 'noisy' feature. It should be either
        'random' or 'deterministic'. Default is 'random'. When it is
        'deterministic' and n_runs is not specified, n_runs defaults to 1.
    noise_map : str or callable, optional
        The deterministic scalar noise map in the 'noisy' feature. It should
        be either 'chebyshev' or a callable ``x -> noise`` returning a real
        scalar. It is used only when noise_mode is 'deterministic'. Default
        is 'chebyshev'. The built-in 'chebyshev' map follows the deterministic
        noise model in Moré and Wild [5]_.
    significant_digits : int, optional
        The number of significant digits in the
        'truncated' feature. Default is 6.
    perturbed_trailing_digits : bool, optional
        Whether we will randomize the trailing
        digits of the objective function value in the 'truncated' feature.
        Default is False.
    rotated : bool, optional
        Whether to use a random or given rotation matrix to rotate
        the coordinates of a problem in the 'linearly_transformed' feature.
        Default is True.
    condition_factor : float, optional
        The scaling factor of the condition number of the
        linear transformation in the 'linearly_transformed' feature. More
        specifically, the condition number of the linear transformation will
        be ``2^(condition_factor * n / 2)``, where n is the dimension of the
        problem. Default is 0.
    nan_rate : float, optional
        The probability that the evaluation of the objective
        function will return np.nan in the 'random_nan' feature. Default is
        0.05.
    unrelaxable_bounds : bool, optional
        Whether the bound constraints are unrelaxable or
        not in the 'unrelaxable_constraints' feature. Default is True.
    unrelaxable_linear_constraints : bool, optional
        Whether the linear constraints are
        unrelaxable or not in the 'unrelaxable_constraints' feature. Default
        is False.
    unrelaxable_nonlinear_constraints : bool, optional
        Whether the nonlinear constraints
        are unrelaxable or not in the 'unrelaxable_constraints' feature.
        Default is False.
    mesh_size : float, optional
        The size of the mesh in the 'quantized' feature. Default
        is 1e-3.
    mesh_type : str, optional
        The type of the mesh in the 'quantized' feature. It should
        be either 'absolute' or 'relative'. Default is 'absolute'.
    ground_truth : bool, optional
        Whether the featured problem is the ground truth or not
        in the 'quantized' feature. Default is True.
    mod_x0 : callable, optional
        The modifier function to modify the initial guess in the
        'custom' feature. It should be a callable
        ``(random_stream, problem) -> modified_x0``,
        where problem is an instance of the class Problem, and
        modified_x0 is the modified initial guess. No default.
    mod_affine : callable, optional
        The modifier function to generate the affine
        transformation applied to the variables in the 'custom' feature. It
        should be a callable
        ``(random_stream, problem) -> (A, b, inv)``,
        where problem is an instance of the class Problem, A is the
        matrix of the affine transformation, b is the vector of the affine
        transformation, and inv is the inverse of matrix A. No default.
    mod_bounds : callable, optional
        The modifier function to modify the bound constraints in
        the 'custom' feature. It should be a callable
        ``(random_stream, problem) -> (modified_xl, modified_xu)``,
        where problem is an instance of the class Problem, modified_xl is
        the modified lower bound, and modified_xu is the modified upper
        bound. No default.
    mod_linear_ub : callable, optional
        The modifier function to modify the linear inequality
        constraints in the 'custom' feature. It should be a callable
        ``(random_stream, problem) -> (modified_aub, modified_bub)``,
        where problem is an instance of the class Problem, modified_aub
        is the modified matrix of the linear inequality constraints, and
        modified_bub is the modified vector of the linear inequality
        constraints. No default.
    mod_linear_eq : callable, optional
        The modifier function to modify the linear equality
        constraints in the 'custom' feature. It should be a callable
        ``(random_stream, problem) -> (modified_aeq, modified_beq)``,
        where problem is an instance of the class Problem, modified_aeq
        is the modified matrix of the linear equality constraints, and
        modified_beq is the modified vector of the linear equality
        constraints. No default.
    mod_fun : callable, optional
        The modifier function to modify the objective function in
        the 'custom' feature. It should be a callable
        ``(x, random_stream, problem) -> modified_fun``,
        where x is the evaluation point, problem is an instance of the
        class Problem, and modified_fun is the modified objective function
        value. No default.
    mod_cub : callable, optional
        The modifier function to modify the nonlinear inequality
        constraints in the 'custom' feature. It should be a callable
        ``(x, random_stream, problem) -> modified_cub``,
        where x is the evaluation point, problem is an instance of the
        class Problem, and modified_cub is the modified vector of the
        nonlinear inequality constraints. No default.
    mod_ceq : callable, optional
        The modifier function to modify the nonlinear equality
        constraints in the 'custom' feature. It should be a callable
        ``(x, random_stream, problem) -> modified_ceq``,
        where x is the evaluation point, problem is an instance of the
        class Problem, and modified_ceq is the modified vector of the
        nonlinear equality constraints. No default.

    Profile and plot options:

    bar_colors : list or numpy.ndarray, optional
        Two different colors for the bars of two solvers in the
        log-ratio profiles. It can be a list of short names of colors
        ('r', 'g', 'b', 'c', 'm', 'y', 'k') or a 2-by-3 array with each row
        being a RGB triplet. Default is set to the first two colors in the
        'line_colors' option.
    benchmark_id : str, optional
        The identifier of the test. It is used to create the
        specific directory to store the results. Default is 'out' if the
        option 'load' is not provided, otherwise default is '.'.
    draw_hist_plots : str, optional
        Whether or how to draw the history plots of all the problems. It can
        be either 'none', 'sequential', or 'parallel'. If it is 'none', we
        will not draw the history plots. If it is 'parallel', we will draw
        the history plots at the same time when solvers are solving the
        problems. If it is 'sequential', we will draw the history plots
        after all the problems are solved. For each problem, the original
        combined history plot is saved as ``PROBLEM.pdf``; raw-only and
        cumulative-minimum-only plots are saved under the ``raw`` and
        ``cummin`` subfolders. Default is 'parallel'.
    errorbar_type : str, optional
        The type of the uncertainty interval that can be
        either 'minmax' or 'meanstd'. When 'n_runs' is greater than 1, we run
        several times of the experiments and get average curves and
        uncertainty intervals. Default is 'minmax', meaning that we takes the
        pointwise minimum and maximum of the curves.
    feature_stamp : str, optional
        The stamp of the feature with the given options. It is
        used to create the specific directory to store the results. Default
        depends on features.
    hist_aggregation : str, optional
        The aggregation method we use to reduce the number of points in the
        history plots. It can be 'min', 'mean', or 'max'. Default is 'min'.
    line_colors : list, optional
        The colors of the lines in the plots. It can be a list of any valid
        matplotlib colors (short names, hex strings, RGB tuples, etc.).
        Default line colors are from the matplotlib tab10 color cycle. Note
        that if the number of solvers is greater than the number of colors,
        we will cycle through the colors.
    line_styles : list of str, optional
        The styles of the lines in the plots. It can be a list
        of strs that are the combinations of line styles ('-', '-.',
        '--', ':') and markers ('o', '+', '*', '.', 'x', 's', 'd', '^', 'v',
        '>', '<', 'p', 'h'). Default line style order is ['-', '-.', '--',
        ':']. Note that if the number of solvers is greater than the number
        of line styles, we will cycle through the styles.
    line_widths : float or list, optional
        The widths of the lines in the plots. It should be a
        positive float or a list. Default is 1.5. Note that if the number
        of solvers is greater than the number of line widths, we will cycle
        through the widths.
    load : str, optional
        Loading the stored data from a completed experiment and draw
        profiles. It can be either 'latest' or a time stamp of an experiment
        in the format of 'yyyyMMdd_HHmmss'. No default. Note that if solvers is None,
        this key must be provided to load data from a previous experiment
        and generate profiles.
    max_eval_factor : int, optional
        The factor multiplied to each problem's dimension to
        get the maximum number of evaluations for each problem. Default is
        500.
    max_tol_order : int, optional
        The maximum order of the tolerance. In any profile
        (performance profiles, data profiles, and log-ratio profiles), we
        need to set a group of 'tolerances' to define the convergence test of
        the solvers. (Details can be found in the references.) We will set
        the tolerances as ``10**(-k)`` for ``k = 1, 2, ..., max_tol_order``.
        Default is 10.
    merit_fun : callable, optional
        The merit function to measure the quality of a point using
        the objective function value and the maximum constraint violation.
        It should be a callable
        ``(fun_value, maxcv_value, maxcv_init) -> merit_value``,
        where fun_value is the objective function value, maxcv_value is
        the maximum constraint violation, and maxcv_init is the maximum
        constraint violation at the initial guess. The default merit function
        varphi(x) is defined by the objective function f(x) and the maximum
        constraint violation v(x) as::

            varphi(x) = f(x)                        if v(x) <= v1
            varphi(x) = f(x) + 1e5 * (v(x) - v1)   if v1 < v(x) <= v2
            varphi(x) = np.inf                       if v(x) > v2

        where v1 = min(0.01, 1e-10 * max(1, v0)), v2 = max(0.1, 2 * v0),
        and v0 is the maximum constraint violation at the initial guess.
        If varphi(x_0) is inf for a problem/run, the convergence test is
        degenerate; by convention, all solvers are declared to pass that
        problem/run. These cases are listed in ``test_log/report.txt``.
    n_jobs : int, optional
        The number of parallel jobs to run the test. Default is a
        conservative number of workers, chosen as about half of the
        available workers (at least 2 when more than one worker is
        available).
    normalized_scores : bool, optional
        Whether to normalize the scores of the solvers by
        the maximum score of the solvers. Default is True.
    project_x0 : bool, optional
        Whether to project the initial point to the feasible set.
        Default is False.
    run_plain : bool, optional
        Whether to run an extra experiment with the 'plain'
        feature. Default is False.
    savepath : str, optional
        The path to store the results. Default is the current
        working directory.
    score_fun : callable, optional
        The scoring function to calculate the scores of the
        solvers. It should be a callable
        ``profile_scores -> solver_scores``,
        where profile_scores is a 4D array containing scores for all
        profiles. The first dimension of profile_scores corresponds to the
        index of the solver, the second corresponds to the index of tolerance
        starting from 1, the third represents history-based or output-based
        profiles, and the fourth represents performance profiles, data
        profiles, or log-ratio profiles. The default scoring function takes
        the average of the history-based performance profiles under all the
        tolerances.
    score_only : bool, optional
        Whether to only calculate the scores of the solvers
        without drawing the profiles and saving the data. Default is False.
    score_weight_fun : callable, optional
        The weight function to calculate the scores of the
        solvers in the performance and data profiles. It should be a callable
        representing a nonnegative function in R^+. Default is a constant
        function returning 1.
    seed : int, optional
        The seed of the random number generator. Default is 0.
    semilogx : bool, optional
        Whether to use the semilogx scale during plotting profiles
        (performance profiles and data profiles). For data profiles, the
        plotted coordinate is ``log2(1 + n_eval / (n + 1))`` but tick labels
        are displayed in the original ``n_eval / (n + 1)`` units, typically
        ``0, 1, 2, 4, ...``. Default is True.
    silent : bool, optional
        Whether to show the information of the progress. Default is
        False.
    solver_isrand : list of bool, optional
        Whether the solvers are randomized or not. Default is
        a list of bools of the same length as the number of solvers, where
        the value is True if the solver is randomized, and False otherwise.
        Note that if 'n_runs' is not specified, we will set it 5 for the
        randomized solvers.
    solver_names : list of str, optional
        The names of the solvers. Default is the names of the
        callables in solvers.
    solver_verbose : int, optional
        The level of the verbosity of the solvers. 0 means
        no verbosity, 1 means some verbosity, and 2 means full verbosity.
        Default is 1.
    solvers_to_load : list of int, optional
        The indices of the solvers to load when the 'load'
        option is provided. It can be a list of different integers selected
        from 0 to the total number of solvers minus 1 of the loading experiment. At
        least two indices should be provided. Default is all the solvers.
    summarize_data_profiles : bool, optional
        Whether to add all the data profiles to the
        summary PDF. Default is True.
    summarize_log_ratio_profiles : bool, optional
        Whether to add all the log-ratio
        profiles to the summary PDF. Default is False.
    summarize_output_based_profiles : bool, optional
        Whether to add all the output-based
        profiles of the selected profiles to the summary PDF. Default is
        True.
    summarize_performance_profiles : bool, optional
        Whether to add all the performance
        profiles to the summary PDF. Default is True.
    xlabel_data_profile : str, optional
        The label of the x-axis of the data profiles.
        Default is 'Number of simplex gradients'.
        Note: LaTeX formatting is supported. The same applies to the options
        'xlabel_log_ratio_profile', 'xlabel_performance_profile',
        'ylabel_data_profile', 'ylabel_log_ratio_profile', and
        'ylabel_performance_profile'.
    xlabel_log_ratio_profile : str, optional
        The label of the x-axis of the log-ratio
        profiles. Default is 'Problem'.
    xlabel_performance_profile : str, optional
        The label of the x-axis of the
        performance profiles. Default is 'Performance ratio'.
    ylabel_data_profile : str, optional
        The label of the y-axis of the data profiles.
        Default is 'Data profiles ($\\mathrm{tol} = %s$)', where '%s' will be
        replaced by the current tolerance in LaTeX format. You can also use
        '%s' in your custom label, and it will be replaced accordingly. The
        same applies to the options 'ylabel_log_ratio_profile' and
        'ylabel_performance_profile'.
    ylabel_log_ratio_profile : str, optional
        The label of the y-axis of the log-ratio
        profiles. Default is 'Log-ratio profiles ($\\mathrm{tol} = %s$)',
        where '%s' will be replaced by the current tolerance in LaTeX format.
    ylabel_performance_profile : str, optional
        The label of the y-axis of the
        performance profiles. Default is
        'Performance profiles ($\\mathrm{tol} = %s$)', where '%s' will be
        replaced by the current tolerance in LaTeX format.

    Problem options:

    Options in this part are used to select problems for benchmarking.
    First select which problem libraries to use based on the 'plibs'
    option. Then select problems from these libraries according to the
    given options ('problem_names', 'ptype', 'mindim', 'maxdim', 'minb',
    'maxb', 'minlcon', 'maxlcon', 'minnlcon', 'maxnlcon', 'mincon',
    'maxcon', and 'excludelist').

    plibs : list of str, optional
        The problem libraries to be used. It should be a list of strs.
        The built-in choices are ``'s2mpj'``, ``'pycutest'``,
        ``'solar'``, and
        ``'custom'``. Default setting is ``'s2mpj'``. Note that
        ``'pycutest'`` requires the separate installation of the
        ``pycutest`` package; see https://jfowkes.github.io/pycutest/
        for installation instructions. ``'solar'`` uses a slim
        SOLAR runtime and may be substantially slower than algebraic test
        problems because it calls an external simulator.
        You can also use your own problem library by specifying its name here
        together with the ``custom_problem_libs_path`` option.
    ptype : str, optional
        The type of the problems to be selected. It should be a str
        consisting of any combination of 'u' (unconstrained), 'b'
        (bound constrained), 'l' (linearly constrained), and 'n' (nonlinearly
        constrained), such as 'b', 'ul', 'ubn'. Default is 'u'.
    mindim : int, optional
        The minimum dimension of the problems to be selected. Default
        is 1.
    maxdim : int, optional
        The maximum dimension of the problems to be selected. Default
        is mindim + 1.
    minb : int, optional
        The minimum number of bound constraints of the problems to be
        selected. Default is 0.
    maxb : int, optional
        The maximum number of bound constraints of the problems to be
        selected. Default is minb + 10.
    minlcon : int, optional
        The minimum number of linear constraints of the problems to
        be selected. Default is 0.
    maxlcon : int, optional
        The maximum number of linear constraints of the problems to
        be selected. Default is minlcon + 10.
    minnlcon : int, optional
        The minimum number of nonlinear constraints of the problems
        to be selected. Default is 0.
    maxnlcon : int, optional
        The maximum number of nonlinear constraints of the problems
        to be selected. Default is minnlcon + 10.
    mincon : int, optional
        The minimum number of linear and nonlinear constraints of the
        problems to be selected. Default is min(minlcon, minnlcon).
    maxcon : int, optional
        The maximum number of linear and nonlinear constraints of the
        problems to be selected. Default is max(maxlcon, maxnlcon).
    custom_problem_libs_path : str or Path, optional
        The path to a directory containing custom problem libraries, or a direct
        path to one custom problem library. Each custom problem library must
        contain a tools module named '<library_name>_tools.py' with two functions:
        '<library_name>_load' and '<library_name>_select'. This option allows
        users to use their own problem libraries without modifying the installed
        package. If a custom library has the same name as a built-in library,
        the custom library is used. Default is None, meaning only built-in
        libraries are available.
    excludelist : list, optional
        The list of problems to be excluded. Default is not to
        exclude any problem.
    problem_names : list of str, optional
        The names of the problems to be selected. It should
        be a list of strs. Default is not to select any
        problem by name but by the options above.
    problem : Problem, optional
        A problem to be benchmarked. It should be an instance of the
        class Problem. If it is provided, we will only run the test on this
        problem with the given feature and draw the history plots. Default is
        not to set any problem.

    Returns
    -------
    solver_scores : numpy.ndarray
        Scores of the solvers based on the profiles. See 'score_fun' in
        'Other Parameters' for more details.
    profile_scores : numpy.ndarray or None
        A 4D array containing scores for all profiles. The first dimension corresponds to the
        index of the solver, the second to the index of tolerance starting from 1, the third
        represents history-based or output-based profiles, and the fourth represents
        performance profiles, data profiles, or log-ratio profiles.
    curves : list of dict or None
        A list containing the curves of all the profiles.

    Raises
    ------
    TypeError
        If an argument received an invalid value.
    ValueError
        If the arguments are inconsistent.

    See Also
    --------
    Problem : Representation of optimization problems.
    Feature : Feature applied to problems during benchmarking.
    FeaturedProblem : Problem equipped with a specific feature.

    Notes
    -----
    The current version supports benchmarking derivative-free optimization
    solvers.

    .. caution::

        The log-ratio profiles are available only when there are exactly two
        solvers. For more information on performance and data profiles, see
        [1]_, [2]_, [5]_. For that of log-ratio profiles, see [4]_, [6]_.

    .. caution::

        All callable arguments (``solvers``, ``distribution``, ``noise_map``,
        ``mod_x0``, ``mod_affine``, ``mod_bounds``, ``mod_linear_ub``,
        ``mod_linear_eq``, ``mod_fun``, ``mod_cub``, ``mod_ceq``,
        ``merit_fun``, ``score_fun``, ``score_weight_fun``) must be
        picklable for parallel execution
        (``n_jobs > 1``). In particular, **lambda functions are not
        picklable** and will cause the benchmark to fall back to sequential
        mode automatically. To take advantage of parallel execution, define
        named functions (using ``def``) instead of lambda expressions.

    1. Several problem libraries are available by default:
       `S2MPJ <https://github.com/GrattonToint/S2MPJ>`_ (see [3]_) and
       `PyCUTEst <https://jfowkes.github.io/pycutest/>`_ (Linux and macOS
       only), and SOLAR (see [7]_) through the ``solar`` adapter.
       SOLAR is distributed through its own LGPL-2.1 runtime files; see the
       adapter README and ``runtime/solar/manifest.json`` for provenance. To
       use your own problem library, see the
       ``custom_problem_libs_path`` option or the guide on our
       `website <https://www.optprof.com>`_. A custom library with the same
       name as a built-in library overrides the built-in one for that run.

    2. Each problem library has a ``config.txt`` file that controls options
       such as ``variable_size`` and ``test_feasibility_problems``. You can
       override these at runtime using `set_plib_config` or by setting the
       corresponding environment variables (e.g.,
       ``S2MPJ_VARIABLE_SIZE``). See `get_plib_config` and
       `set_plib_config` for details.

    3. When the ``load`` option is provided, the function loads data from a
       previous experiment and draws profiles using the provided options.
       Available options in this mode are:

       - *Profile and plot options*: ``benchmark_id``, ``solver_names``,
         ``feature_stamp``, ``errorbar_type``, ``savepath``,
         ``max_tol_order``, ``merit_fun``, ``run_plain``, ``score_only``,
         ``summarize_performance_profiles``,
         ``summarize_data_profiles``,
         ``summarize_log_ratio_profiles``,
         ``summarize_output_based_profiles``, ``silent``, ``semilogx``,
         ``normalized_scores``, ``score_weight_fun``, ``score_fun``,
         ``solvers_to_load``, ``line_colors``, ``line_styles``,
         ``line_widths``, ``bar_colors``.
       - *Feature options*: none.
       - *Problem options*: ``plibs``, ``ptype``, ``mindim``, ``maxdim``,
         ``minb``, ``maxb``, ``minlcon``, ``maxlcon``, ``minnlcon``,
         ``maxnlcon``, ``mincon``, ``maxcon``, ``excludelist``.

    4. More information about OptiProfiler can be found at
       https://www.optprof.com.

    References
    ----------
    .. [1] E. D. Dolan and J. J. Moré. Benchmarking optimization software with
           performance profiles. *Math. Program.*, 91(2):201–213, 2002.
           doi:10.1007/s101070100263
           <https://doi.org/10.1007/s101070100263>.
    .. [2] N. Gould and J. Scott. A note on performance profiles for
           benchmarking software. *ACM Trans. Math. Software*, 43(2):15:1–5,
           2016. doi:10.1145/2950048 <https://doi.org/10.1145/2950048>.
    .. [3] S. Gratton and Ph. L. Toint. S2MPJ and CUTEst optimization problems
           for Matlab, Python and Julia. arXiv:2407.07812, 2024.
    .. [4] J. L. Morales. A numerical study of limited memory BFGS methods.
           *Appl. Math. Lett.*, 15(4):481–487, 2002.
           doi:10.1016/S0893-9659(01)00162-8
           <https://doi.org/10.1016/S0893-9659(01)00162-8>.
    .. [5] J. J. Moré and S. M. Wild. Benchmarking derivative-free optimization
           algorithms. *SIAM J. Optim.*, 20(1):172–191, 2009.
           doi:10.1137/080724083 <https://doi.org/10.1137/080724083>.
    .. [6] H.-J. M. Shi, M. Q. Xuan, F. Oztoprak, and J. Nocedal. On the
           numerical performance of finite-difference-based methods for
           derivative-free optimization. *Optim. Methods Softw.*,
           38(2):289–311, 2023. doi:10.1080/10556788.2022.2121832
           <https://doi.org/10.1080/10556788.2022.2121832>.
    .. [7] N. Andrés-Thió, C. Audet, M. Diago, A. E. Gheribi,
           S. Le Digabel, X. Lebeuf, M. Lemyre-Garneau, and C. Tribes.
           ``solar``: A solar thermal power plant simulator for blackbox
           optimization benchmarking. *Optimization and Engineering*, 2025.
           doi:10.1007/s11081-024-09952-x
           <https://doi.org/10.1007/s11081-024-09952-x>.

    Examples
    --------
    Benchmark two toy solvers that randomly sample around the initial point.
    Save the following as a ``.py`` file and run it (the
    ``if __name__ == '__main__'`` guard is required on macOS/Windows
    because ``benchmark`` uses multiprocessing):

    .. code-block:: python

        import numpy as np
        from optiprofiler import benchmark

        def solver1(fun, x0):
            rng = np.random.default_rng(0)
            best_x, best_f = x0, fun(x0)
            for _ in range(10 * len(x0)):
                x = x0 + rng.standard_normal(len(x0))
                if fun(x) < best_f:
                    best_x, best_f = x, fun(x)
            return best_x

        def solver2(fun, x0):
            rng = np.random.default_rng(0)
            best_x, best_f = x0, fun(x0)
            for _ in range(20 * len(x0)):
                x = x0 + rng.standard_normal(len(x0))
                if fun(x) < best_f:
                    best_x, best_f = x, fun(x)
            return best_x

        if __name__ == '__main__':
            scores = benchmark([solver1, solver2])
    """
    logger = get_logger(__name__)

    # Check whether solvers or 'load' option is given.
    if solvers is None and 'load' not in kwargs:
        raise ValueError('Either solvers or the \'load\' option must be given.')

    # Preprocess the solvers if given.
    if solvers is not None:
        if not hasattr(solvers, '__len__') or not all(callable(solver) for solver in solvers):
            raise TypeError('The solvers must be a list of callables.')
        if len(solvers) < 2:
            raise ValueError('At least two solvers must be given.')
        solvers = list(solvers)

    # Save the original keyword arguments for future use.
    options_user = kwargs.copy()

    # Process the feature name.
    if 'feature_name' in kwargs:
        feature_name = kwargs.pop('feature_name')
    else:
        feature_name = FeatureName.PLAIN.value
    if feature_name not in FeatureName.__members__.values():
        raise ValueError(f'Unknown feature name: {feature_name}.')

    # Process the problem if provided.
    if 'problem' in kwargs and kwargs['problem'] is not None:
        problem = kwargs.pop('problem')

    # Get the different options from the keyword arguments.
    feature_options = {}
    profile_options = {}
    problem_options = {}
    for key, value in kwargs.items():
        if key in FeatureOption.__members__.values():
            feature_options[key] = value
        elif key in ProblemOption.__members__.values():
            problem_options[key] = value
        elif key in ProfileOption.__members__.values():
            profile_options[key] = value
        else:
            raise ValueError(f'Unknown option: {key}.')

    # Check validity of the options.
    problem_options = check_validity_problem_options(problem_options)
    profile_options = check_validity_profile_options(solvers, profile_options)

    # Whether to load the existing results.
    is_load = ProfileOption.LOAD in profile_options and profile_options[ProfileOption.LOAD] is not None

    # If `n_runs` is not specified, set it to 5 if at least one solver is randomized.
    any_solver_isrand = ProfileOption.SOLVER_ISRAND in profile_options and any(profile_options[ProfileOption.SOLVER_ISRAND])
    if FeatureOption.N_RUNS not in feature_options and any_solver_isrand and not is_load:
        if ProfileOption.SILENT not in profile_options or not profile_options[ProfileOption.SILENT]:
            print_log_message('INFO', f'We set {FeatureOption.N_RUNS} to 5 since it is not specified and at least one solver is randomized.')
        feature_options[FeatureOption.N_RUNS] = 5

    # Load the existing results if needed.
    # If 'load' is specified, we skip the solving phase and restore the results from disk.
    if is_load:
        results_plibs, profile_options = load_results(problem_options, profile_options)
    
    # Build feature.
    feature = Feature(feature_name, **feature_options)
    feature_options = feature.options
    
    # Set default values for the unspecified options.
    problem_options = get_default_problem_options(problem_options)
    profile_options = get_default_profile_options(solvers, feature, profile_options)

    # Initialize outputs.
    if 'results_plibs' in locals() and results_plibs is not None:
        n_solvers = results_plibs[0]['fun_histories'].shape[1]
    elif ProfileOption.SOLVER_NAMES in profile_options:
        n_solvers = len(profile_options[ProfileOption.SOLVER_NAMES])
    else:
        n_solvers = len(solvers)

    solver_scores = np.zeros(n_solvers)
    profile_scores = []
    curves = []
    solver_names = profile_options[ProfileOption.SOLVER_NAMES]
    many_solver_warning = None
    if n_solvers > 10 and not profile_options[ProfileOption.SILENT]:
        many_solver_warning = f'Comparing {n_solvers} solvers at once can make profiles hard to interpret; pairwise comparisons are recommended.'

    # Set the options for plotting.
    profile_context = set_profile_context(profile_options)

    # Define the directory to store the results.
    path_out = Path(profile_options[ProfileOption.SAVEPATH], profile_options[ProfileOption.BENCHMARK_ID]).resolve()

    # Record the time stamp of the current experiment. ('yyyyMMdd_HHmmss')
    time_stamp = datetime.now().astimezone().strftime('%Y%m%d_%H%M%S')

    # Set the feature stamp.
    feature_stamp = profile_options[ProfileOption.FEATURE_STAMP]

    # Create the stamp for the current experiment.
    stamp = create_stamp(solver_names, problem_options, feature_stamp, time_stamp, path_out)

    path_stamp = path_out / stamp
    path_log = path_stamp / 'test_log'
    path_report = path_log / 'report.txt'

    # Create the directory to store the results.
    if not profile_options[ProfileOption.SCORE_ONLY]:
        if not path_stamp.exists():
            path_stamp.mkdir(parents=True, exist_ok=True)
        else:
            shutil.rmtree(path_stamp)
            path_stamp.mkdir(parents=True, exist_ok=True)
    
    # Create directory to store history plots based on draw_hist_plots option.
    if profile_options[ProfileOption.DRAW_HIST_PLOTS] == 'none':
        path_hist_plots = None
    else:
        path_hist_plots = path_stamp / 'history_plots'
        if not path_hist_plots.exists():
            path_hist_plots.mkdir(parents=True, exist_ok=True)
    
    # Create the directory to store options and log files.
    if not profile_options[ProfileOption.SCORE_ONLY]:
        if not path_log.exists():
            path_log.mkdir(parents=True, exist_ok=True)
        else:
            shutil.rmtree(path_log)
            path_log.mkdir(parents=True, exist_ok=True)

    # Create a README.txt file to explain the content of the folder `path_log`.
    if not profile_options[ProfileOption.SCORE_ONLY]:
        path_readme_log = path_log / 'README.txt'
        init_readme(path_readme_log)
    else:
        path_readme_log = None

    # Create a text file named by time_stamp to record the time_stamp.
    if not profile_options[ProfileOption.SCORE_ONLY]:
        path_time_stamp = path_log / f'time_stamp_{time_stamp}.txt'
        try:
            with path_time_stamp.open('w') as f:
                f.write(f'{time_stamp}')
            add_to_readme(path_readme_log, f'time_stamp_{time_stamp}.txt', 'File, recording the time stamp of the current experiment.')
        except Exception as exc:
            if not profile_options[ProfileOption.SILENT]:
                print_log_message('WARNING', f'Failed to create the time stamp file in {path_log}.')
                print_log_message('WARNING', f'Error message: {shorten_log_message(exc)}')
        
    if not profile_options[ProfileOption.SCORE_ONLY] and 'problem' not in locals():
        # path_figs = path_log / 'profile_figs'
        # path_figs.mkdir(parents=True, exist_ok=True)
        # try:
        #     with path_readme_log.open('a') as f:
        #         f.write(f"'profile_figs': folder, containing all the FIG files of the profiles.\n")
        # except:
        #     pass
        pass

    # Save the options and record the log.
    # Initialize log_queue and listener to None for the case when score_only=True.
    log_queue = None
    listener = None
    if not profile_options[ProfileOption.SCORE_ONLY]:
        try:
            # Save the user-provided options.
            # We use pickle because options may contain function handles or other non-serializable objects.
            if 'options_user' in locals() and options_user is not None:
                save_options(options_user, path_log / 'options_user.pkl')

                add_to_readme(path_readme_log, 'options_user.pkl', 'File, storing the options provided by the user for the current experiment.')

            # Save the refined options (including defaults and internal settings).
            options_refined = profile_options.copy()
            feature_options_keys = list(feature_options.keys())
            problem_options_keys = list(problem_options.keys())
            for key in feature_options_keys:
                options_refined[key] = feature_options[key]
            for key in problem_options_keys:
                options_refined[key] = problem_options[key]
            
            save_options(options_refined, path_log / 'options_refined.pkl')
            add_to_readme(path_readme_log, 'options_refined.pkl', 'File, storing the options refined by OptiProfiler for the current experiment.')
        except Exception as exc:
            if not profile_options[ProfileOption.SILENT]:
                print_log_message('WARNING', 'Failed to save the options of the current experiment.')
                print_log_message('WARNING', f'Error message: {shorten_log_message(exc)}')

        # Set up the logger to log the information in a file.
        log_file = path_log / 'log.txt'
        log_queue, listener = setup_main_process_logging(log_file=log_file, level=logging.INFO)

        add_to_readme(path_readme_log, 'log.txt', 'File, the log file of the current experiment, recording printed information from the screen.')

    if many_solver_warning is not None:
        if listener is not None:
            logger.warning(many_solver_warning)
        else:
            print_log_message('WARNING', many_solver_warning)

    # Create a README.txt file to explain the content of the folder `path_stamp`.
    if not profile_options[ProfileOption.SCORE_ONLY]:
        path_readme_feature = path_stamp / 'README.txt'
        try:
            init_readme(path_readme_feature)
            if path_hist_plots is not None:
                add_to_readme(path_readme_feature, 'history_plots', 'Folder, containing combined history plots for each problem; raw-only and cumulative-minimum-only plots are stored in the raw and cummin subfolders of each problem-library folder.')
                add_to_readme(path_readme_feature, 'history_plots_summary.pdf', 'File, the summary PDF of combined history plots for all problems.')
            add_to_readme(path_readme_feature, 'test_log', 'Folder, containing log files and other useful experimental data.')
        except Exception as exc:
            if not profile_options[ProfileOption.SILENT]:
                logger.warning(f'Fail to create the README.txt file in {path_stamp}.')
                logger.warning(f'Error message: {shorten_log_message(exc)}')

    # We try to copy the script or function that calls the benchmark function to the log directory.
    if not profile_options[ProfileOption.SCORE_ONLY]:
        try:
            caller_path = inspect.stack()[1].filename
            if caller_path is not None and Path(caller_path).is_file():
                shutil.copy2(caller_path, path_log / os.path.basename(caller_path))
                add_to_readme(path_readme_log, os.path.basename(caller_path), 'File, the script or function that calls the benchmark function.')
                if not profile_options[ProfileOption.SILENT]:
                    logger.info(f'The script or function that calls `benchmark` function is copied to:')
                    logger.info(f'{path_log}')
        except Exception as exc:
            if not profile_options[ProfileOption.SILENT]:
                logger.warning(f'Failed to copy the script or function that calls `benchmark` function to the log directory.')
                logger.warning(f'Error message: {_shorten_log_message(exc)}')

    # Create the directories for the performance profiles, data profiles, and log-ratio profiles.
    if not profile_options[ProfileOption.SCORE_ONLY] and 'problem' not in locals():
        path_perf_hist = path_stamp / 'detailed_profiles' / 'perf_history-based'
        path_data_hist = path_stamp / 'detailed_profiles' / 'data_history-based'
        path_log_ratio_hist = path_stamp / 'detailed_profiles' / 'log-ratio_history-based'
        path_perf_out = path_stamp / 'detailed_profiles' / 'perf_output-based'
        path_data_out = path_stamp / 'detailed_profiles' / 'data_output-based'
        path_log_ratio_out = path_stamp / 'detailed_profiles' / 'log-ratio_output-based'

        if not path_perf_hist.exists():
            path_perf_hist.mkdir(parents=True, exist_ok=True)
        if not path_data_hist.exists():
            path_data_hist.mkdir(parents=True, exist_ok=True)
        if not path_perf_out.exists():
            path_perf_out.mkdir(parents=True, exist_ok=True)
        if not path_data_out.exists():
            path_data_out.mkdir(parents=True, exist_ok=True)
        if n_solvers == 2:
            if not path_log_ratio_hist.exists():
                path_log_ratio_hist.mkdir(parents=True, exist_ok=True)
            if not path_log_ratio_out.exists():
                path_log_ratio_out.mkdir(parents=True, exist_ok=True)

        try:
            add_to_readme(path_readme_feature, 'detailed_profiles', 'Folder, containing all the high-quality single profiles.')
        except:
            pass

    path_perf_hist_summary = path_stamp / 'perf_hist.pdf'
    path_perf_out_summary = path_stamp / 'perf_out.pdf'
    path_data_hist_summary = path_stamp / 'data_hist.pdf'
    path_data_out_summary = path_stamp / 'data_out.pdf'
    path_log_ratio_hist_summary = path_stamp / 'log-ratio_hist.pdf' 
    path_log_ratio_out_summary = path_stamp / 'log-ratio_out.pdf'

    # If a specific problem is provided to `problem_options`, we only solve this problem and generate the history plots for it.
    if 'problem' in locals():
        profile_options_log, _ = _with_solver_log_names(profile_options, len(problem.name))
        if not profile_options[ProfileOption.SILENT]:
            _log_solver_aliases(logger, solver_names, profile_options_log[_SOLVER_LOG_NAMES_KEY])
        result = _solve_one_problem(solvers, problem, feature, problem.name, len(problem.name), profile_options_log, True, path_hist_plots)

        if not profile_options[ProfileOption.SCORE_ONLY]:
            # We move the history plots to the feature directory.
            try:
                for file in path_hist_plots.iterdir():
                    shutil.move(file, path_stamp / file.name)
                shutil.rmtree(path_hist_plots)
                if not profile_options[ProfileOption.SILENT]:
                    logger.info(f'Detailed results stored in: {path_stamp}')
            except:
                pass
        
        # Compute merit values. ``result['fun_init']`` /
        # ``result['maxcv_init']`` are now per-run vectors of shape
        # ``(n_runs,)``; the resulting ``merit_init`` therefore also has
        # shape ``(n_runs,)`` and the per-run ``solver_scores`` are
        # averaged across runs at the end (mirroring MATLAB
        # ``benchmark.m`` lines 858-879).
        merit_fun = profile_options[ProfileOption.MERIT_FUN]
        if result is not None:
            try:
                merit_history = compute_merit_values(merit_fun, result['fun_history'], result['maxcv_history'], result['maxcv_init'])
                merit_init = compute_merit_values(merit_fun, result['fun_init'], result['maxcv_init'], result['maxcv_init'])
            except Exception as exc:
                logger.error(f'Error occurred while calculating the merit values. Please check the merit function. Error message: {shorten_log_message(exc)}')
                raise exc
            merit_init = np.atleast_1d(np.asarray(merit_init))
            local_n_solvers, local_n_runs = merit_history.shape[:2]
            # Find the least merit value for each run.
            merit_mins = np.full(local_n_runs, np.nan)
            for i_run in range(local_n_runs):
                merit_mins[i_run] = np.nanmin(merit_history[:, i_run, :])
                merit_mins[i_run] = np.nanmin([merit_mins[i_run], merit_init[i_run]])
            # Since we will not compute the profiles, we set `solver_scores`
            # to be the relative decreases in the objective function value
            # (per run, then averaged across runs).
            solver_merit_mins = np.nanmin(merit_history, axis=2)  # (n_solvers, n_runs)
            solver_scores_runs = np.zeros((local_n_solvers, local_n_runs))
            for i_run in range(local_n_runs):
                denom = max(merit_init[i_run] - merit_mins[i_run], np.finfo(float).eps)
                solver_scores_runs[:, i_run] = (merit_init[i_run] - solver_merit_mins[:, i_run]) / denom
                solver_scores_runs[:, i_run] = np.maximum(solver_scores_runs[:, i_run], 0)
            solver_scores = np.mean(solver_scores_runs, axis=1)
        else:
            solver_scores = np.zeros(n_solvers)

        if not profile_options[ProfileOption.SILENT]:
            logger.info('')
            logger.info('Scores of the solvers:')
            max_solver_name_length = max(len(name) for name in solver_names)
            for i, name in enumerate(solver_names):
                format_info_str = f'{{:<{max_solver_name_length}}}:    {{:.4f}}'
                logger.info(format_info_str.format(name, solver_scores[i]))

        # Close the listener of the logger before returning.
        if not profile_options[ProfileOption.SCORE_ONLY]:
            listener.stop()
            for h in listener.handlers:
                if hasattr(h, 'close'):
                    h.close()
            log_queue.close()
            log_queue.join_thread()
        return solver_scores, None, None

    # If 'load' option is not specified, we solve all the selected problems.
    if not is_load:
        # Print the information of the current experiment.
        if not profile_options[ProfileOption.SILENT]:
            logger.info('')
            logger.info(f'Start testing with the following options:')
            logger.info(f'- Solvers: {", ".join(solver_names)}')
            logger.info(f'- Problem libraries: {", ".join(problem_options[ProblemOption.PLIBS])}')
            logger.info(f'- Problem types: {problem_options[ProblemOption.PTYPE]}')
            logger.info(f'- Problem dimension range: [{problem_options[ProblemOption.MINDIM]}, {problem_options[ProblemOption.MAXDIM]}]')
            if any(t in 'bln' for t in problem_options[ProblemOption.PTYPE]):
                logger.info(f'- Problem mb range: [{problem_options[ProblemOption.MINB]}, {problem_options[ProblemOption.MAXB]}]')
            if any(t in 'ln' for t in problem_options[ProblemOption.PTYPE]):
                logger.info(f'- Problem mlcon range: [{problem_options[ProblemOption.MINLCON]}, {problem_options[ProblemOption.MAXLCON]}]')
            if 'n' in problem_options[ProblemOption.PTYPE]:
                logger.info(f'- Problem mnlcon range: [{problem_options[ProblemOption.MINNLCON]}, {problem_options[ProblemOption.MAXNLCON]}]')
            if any(t in 'ln' for t in problem_options[ProblemOption.PTYPE]):
                logger.info(f'- Problem mcon range: [{problem_options[ProblemOption.MINCON]}, {problem_options[ProblemOption.MAXCON]}]')
            if problem_options[ProblemOption.PROBLEM_NAMES]:
                logger.info(f'- Number of user-provided problem names: {len(problem_options[ProblemOption.PROBLEM_NAMES])}')
            if problem_options[ProblemOption.EXCLUDELIST]:
                logger.info(f'- Number of user-excluded problem names: {len(problem_options[ProblemOption.EXCLUDELIST])}')
            logger.info(f'- Feature stamp: {feature_stamp}')

        # We will solve all the problems from all the problem libraries that user specified in the `problem_options`.
        plibs = problem_options[ProblemOption.PLIBS]
        results_plibs = []
        for i, plib in enumerate(plibs):
            if not profile_options[ProfileOption.SILENT]:
                logger.info('')
                logger.info(f'Start testing problems from the problem library "{plib}" with the "{feature.name}" feature.')
                if plib == 's2mpj':
                    logger.info('')
                    logger.info('More information about the S2MPJ problem library can be found at:')
                    logger.info('https://github.com/GrattonToint/S2MPJ')
                elif plib == 'pycutest':
                    logger.info('')
                    logger.info('More information about the PyCUTEst problem library can be found at:')
                    logger.info('https://jfowkes.github.io/pycutest/_build/html/index.html')
                elif plib == 'solar':
                    logger.info('')
                    logger.info('More information about the SOLAR problem library can be found at:')
                    logger.info('https://github.com/bbopt/solar')
                    logger.info('SOLAR uses an LGPL-2.1 slim runtime and calls an external simulator; some problems can be slow.')

            # Create directory to store the history plots for each problem library.
            path_hist_plots_plib = path_hist_plots / plib if path_hist_plots is not None else None
            if path_hist_plots_plib is not None and not path_hist_plots_plib.exists():
                path_hist_plots_plib.mkdir(parents=True, exist_ok=True)

            # Determine whether to draw history plots during solving (parallel mode) or after (sequential mode).
            is_plot_parallel = profile_options[ProfileOption.DRAW_HIST_PLOTS] == 'parallel'

            # Solve all the problems from the current problem library with the specified options and get the computation results.
            results_plib = _solve_all_problems(solvers, plib, feature, problem_options, profile_options, is_plot_parallel, path_hist_plots_plib, log_queue=log_queue)

            # If there are no problems selected or solved, skip the rest of the code, and continue to the next library.
            if results_plib is None:
                continue

            # Compute the merit values.
            merit_fun = profile_options[ProfileOption.MERIT_FUN]
            try:
                merit_histories = compute_merit_values(merit_fun, results_plib['fun_histories'], results_plib['maxcv_histories'], results_plib['maxcv_inits'])
                merit_outs = compute_merit_values(merit_fun, results_plib['fun_outs'], results_plib['maxcv_outs'], results_plib['maxcv_inits'])
                merit_inits = compute_merit_values(merit_fun, results_plib['fun_inits'], results_plib['maxcv_inits'], results_plib['maxcv_inits'])
            except Exception as exc:
                logger.error(f'Error occurred while calculating the merit values. Please check the merit function. Error message: {shorten_log_message(exc)}')
                raise exc
            results_plib['merit_histories'] = merit_histories
            results_plib['merit_outs'] = merit_outs
            results_plib['merit_inits'] = merit_inits

            # Run the 'plain' feature if run_plain is true.
            if profile_options[ProfileOption.RUN_PLAIN]:
                feature_plain = Feature(FeatureName.PLAIN.value)
                if not profile_options[ProfileOption.SILENT]:
                    logger.info('')
                    logger.info(f'Start testing problems from the problem library "{plib}" with "plain" feature.')
                results_plib_plain = _solve_all_problems(solvers, plib, feature_plain, problem_options, profile_options, False, None, log_queue=log_queue)
                try:
                    results_plib_plain['merit_histories'] = compute_merit_values(merit_fun, results_plib_plain['fun_histories'], results_plib_plain['maxcv_histories'], results_plib_plain['maxcv_inits'])
                    results_plib_plain['merit_outs'] = compute_merit_values(merit_fun, results_plib_plain['fun_outs'], results_plib_plain['maxcv_outs'], results_plib_plain['maxcv_inits'])
                    results_plib_plain['merit_inits'] = compute_merit_values(merit_fun, results_plib_plain['fun_inits'], results_plib_plain['maxcv_inits'], results_plib_plain['maxcv_inits'])
                except Exception as exc:
                    logger.error(f'Error occurred while calculating the merit values for the "plain" feature. Please check the merit function. Error message: {shorten_log_message(exc)}')
                    raise exc
                
                # Store data of the 'plain' feature for later calculating merit_mins.
                results_plib['results_plib_plain'] = results_plib_plain

            # Append the results of the current problem library to the list.
            results_plibs.append(results_plib)

            # Merge the history plots for each problem library to a single pdf file.
            # Only do this in 'parallel' mode: in that case the workers have already
            # written the per-problem PDFs into `path_hist_plots_plib` while solving.
            # In 'sequential' mode, the per-problem PDFs are drawn (and merged) later
            # in the dedicated block below; in 'none' mode there are no PDFs at all.
            if (
                profile_options[ProfileOption.DRAW_HIST_PLOTS] == 'parallel'
                and not profile_options[ProfileOption.SCORE_ONLY]
                and np.any(results_plib['solvers_successes'])
            ):
                if not profile_options[ProfileOption.SILENT]:
                    logger.info(f'Merging all the history plots of problems from the "{plib}" library to a single PDF file.')
                try:
                    merge_pdfs_with_pypdf(path_hist_plots_plib, path_hist_plots / f'{plib}_history_plots_summary.pdf')
                except Exception as exc:
                    if not profile_options[ProfileOption.SILENT]:
                        logger.warning('Failed to merge the history plots to a single PDF file.')
                        logger.warning(f'Error message: {shorten_log_message(exc)}')

        # Remove the None elements from results_plibs.
        results_plibs = [results_plib for results_plib in results_plibs if results_plib is not None]
        if len(results_plibs) == 0:
            if not profile_options[ProfileOption.SILENT]:
                logger.info('')
                logger.info('No problems are selected or solved from any problem library.')
            # Close the listener of the logger before returning.
            if not profile_options[ProfileOption.SCORE_ONLY]:
                listener.stop()
                for h in listener.handlers:
                    if hasattr(h, 'close'):
                        h.close()
                log_queue.close()
                log_queue.join_thread()
            return solver_scores, None, None

    # Draw history plots sequentially if draw_hist_plots is set to 'sequential'.
    if profile_options[ProfileOption.DRAW_HIST_PLOTS] == 'sequential':
        for results_plib in results_plibs:
            plib = results_plib['plib']
            if not profile_options[ProfileOption.SILENT]:
                logger.info(f'Sequentially drawing history plots for problems from the problem library "{plib}".')
                logger.info('This may take a while if there are many problems solved from this problem library.')
                logger.info(f'The history plots for each problem will be saved in: {path_hist_plots / plib}')

            # Create directory to store the history plots for each problem library.
            path_hist_plots_plib = path_hist_plots / plib if path_hist_plots is not None else None
            if path_hist_plots_plib is not None and not path_hist_plots_plib.exists():
                path_hist_plots_plib.mkdir(parents=True, exist_ok=True)

            n_problems_solved = len(results_plib['problem_names'])
            for i_problem in range(n_problems_solved):
                problem_name = results_plib['problem_names'][i_problem]
                problem_type = results_plib['problem_types'][i_problem]
                problem_dim = results_plib['problem_dims'][i_problem]
                solvers_success = results_plib['solvers_successes'][i_problem]
                fun_history = results_plib['fun_histories'][i_problem]
                maxcv_history = results_plib['maxcv_histories'][i_problem]
                fun_init = results_plib['fun_inits'][i_problem]
                maxcv_init = results_plib['maxcv_inits'][i_problem]
                n_eval = results_plib['n_evals'][i_problem]

                # Draw the history plot for this problem.
                _draw_problem_history_plot(
                    problem_name, problem_type, problem_dim, profile_options['solver_names'],
                    solvers_success, fun_history, maxcv_history, fun_init, maxcv_init,
                    n_eval, profile_options, path_hist_plots_plib
                )

            if not profile_options[ProfileOption.SILENT]:
                logger.info(f'Finished drawing history plots for problems from the problem library "{plib}".')

            # Merge the history plots for each problem library to a single pdf file.
            if not profile_options[ProfileOption.SCORE_ONLY] and np.any(results_plib['solvers_successes']):
                if not profile_options[ProfileOption.SILENT]:
                    logger.info(f'Merging all the history plots of problems from the "{plib}" library to a single PDF file.')
                try:
                    merge_pdfs_with_pypdf(path_hist_plots_plib, path_hist_plots / f'{plib}_history_plots_summary.pdf')
                except Exception as exc:
                    if not profile_options[ProfileOption.SILENT]:
                        logger.warning('Failed to merge the history plots to a single PDF file.')
                        logger.warning(f'Error message: {shorten_log_message(exc)}')

    # Store the data for loading.
    # This HDF5 file contains all the numerical results of the experiment.
    # It can be used to reload the experiment state using the 'load' option.
    try:
        save_results_to_h5(results_plibs, path_log / 'data_for_loading.h5')
    except Exception as exc:
        if not profile_options[ProfileOption.SILENT]:
            logger.warning('Failed to save the data of the current experiment.')
            logger.warning(f'Error message: {shorten_log_message(exc)}')



    # Write the report file.
    write_report(profile_options, results_plibs, path_report, path_readme_log)

    # Process the results from all the problem libraries.
    merit_histories_merged, merit_outs_merged, merit_inits_merged, merit_mins_merged, n_evals_merged, problem_names_merged, problem_dims_merged = process_results(results_plibs, profile_options)
    n_problems, n_solvers, n_runs = merit_histories_merged.shape[:3]

    if not profile_options[ProfileOption.SILENT]:
        logger.info('')
        logger.info('Start creating profiles.')

    max_tol_order = profile_options[ProfileOption.MAX_TOL_ORDER]
    tolerances = [10**(-i) for i in range(1, max_tol_order + 1)]

    is_saving = not profile_options[ProfileOption.SCORE_ONLY]

    n_rows = 0
    is_perf = profile_options[ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES]
    is_data = profile_options[ProfileOption.SUMMARIZE_DATA_PROFILES]
    is_log_ratio = profile_options[ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES]
    is_output_based = profile_options[ProfileOption.SUMMARIZE_OUTPUT_BASED_PROFILES]
    if is_perf:
        n_rows += 1
    if is_data:
        n_rows += 1
    if is_log_ratio:
        n_rows += 1
    is_summary = n_rows > 0
    multiplier = 2 if is_output_based else 1
    default_figsize = plt.rcParams['figure.figsize']
    default_width = default_figsize[0]
    default_height = default_figsize[1]
    summary_profile_width = default_width + summary_legend_extra_width(n_solvers, default_width, solver_names)
    summary_width = len(tolerances) * summary_profile_width

    with plt.rc_context(profile_context):
        if is_summary:
            fig_summary = plt.figure(figsize=(summary_width, multiplier * n_rows * default_height), layout='constrained')
            if multiplier == 2:
                fig_summary_hist, fig_summary_out = fig_summary.subfigures(2, 1)
                subfigs_summary_hist = np.atleast_1d(fig_summary_hist.subfigures(n_rows, 1))
                subfigs_summary_out = np.atleast_1d(fig_summary_out.subfigures(n_rows, 1))
            else:
                fig_summary_hist = fig_summary.subfigures(1, 1)
                fig_summary_out = None
                subfigs_summary_hist = np.atleast_1d(fig_summary_hist.subfigures(n_rows, 1))
                subfigs_summary_out = None
            
            i_rows = 0
            if is_perf:
                ax_summary_perf_hist = np.atleast_1d(subfigs_summary_hist[i_rows].subplots(1, len(tolerances), sharey=True))
                ax_summary_perf_out = np.atleast_1d(subfigs_summary_out[i_rows].subplots(1, len(tolerances), sharey=True)) if multiplier == 2 else None
                i_rows += 1
            else:
                ax_summary_perf_hist = None
                ax_summary_perf_out = None
            if is_data:
                ax_summary_data_hist = np.atleast_1d(subfigs_summary_hist[i_rows].subplots(1, len(tolerances), sharey=True))
                ax_summary_data_out = np.atleast_1d(subfigs_summary_out[i_rows].subplots(1, len(tolerances), sharey=True)) if multiplier == 2 else None
                i_rows += 1
            else:
                ax_summary_data_hist = None
                ax_summary_data_out = None
            if is_log_ratio:
                ax_summary_log_ratio_hist = np.atleast_1d(subfigs_summary_hist[i_rows].subplots(1, len(tolerances)))
                ax_summary_log_ratio_out = np.atleast_1d(subfigs_summary_out[i_rows].subplots(1, len(tolerances))) if multiplier == 2 else None
                i_rows += 1
            else:
                ax_summary_log_ratio_hist = None
                ax_summary_log_ratio_out = None
        else:
            fig_summary = None
            fig_summary_hist = None
            fig_summary_out = None
            ax_summary_perf_hist = None
            ax_summary_perf_out = None
            ax_summary_data_hist = None
            ax_summary_data_out = None
            ax_summary_log_ratio_hist = None
            ax_summary_log_ratio_out = None

        # Find the problems that all the solvers failed to meet the convergence test for every tolerance.
        solvers_all_diverge_hist = np.zeros((n_problems, n_runs, profile_options[ProfileOption.MAX_TOL_ORDER]), dtype=bool)
        solvers_all_diverge_out = solvers_all_diverge_hist.copy()

        # Detect "merit_init = Inf" (problem, run) pairs (i.e.
        # phi(x_0) = Inf). For these the standard Moré-Wild convergence
        # threshold tau*phi(x_0) + (1-tau)*phi_min collapses to
        # Inf - Inf = NaN, so the test cannot meaningfully discriminate
        # solvers. By convention we then declare every solver to "pass"
        # the test on these (problem, run) pairs (see the threshold
        # computation below). The problem names will also be recorded in
        # the report under a dedicated section so the user can identify
        # and inspect them.
        # Note: by construction merit_min <= merit_init (see
        # ``process_results`` in ``profile_utils.py``), so the case
        # "merit_init finite, merit_min = Inf" cannot occur; only
        # "merit_init = Inf, merit_min finite" and
        # "merit_init = Inf, merit_min = Inf" need this special handling
        # and they follow the same code path.
        # Backward-compat: ``merit_inits_merged`` may be 1-D for old
        # .h5 files; broadcast to (n_problems, n_runs) here so the rest
        # of the code can treat it uniformly.
        if merit_inits_merged.ndim == 1:
            merit_inits_per_run = np.broadcast_to(
                merit_inits_merged[:, None], (n_problems, n_runs)
            )
        else:
            merit_inits_per_run = merit_inits_merged
        if merit_mins_merged.ndim == 1:
            merit_mins_per_run = np.broadcast_to(
                merit_mins_merged[:, None], (n_problems, n_runs)
            )
        else:
            merit_mins_per_run = merit_mins_merged
        merit_init_inf_mask = np.isinf(merit_inits_per_run)
        if not profile_options[ProfileOption.SILENT] and np.any(merit_init_inf_mask):
            for i_problem in range(n_problems):
                if np.any(merit_init_inf_mask[i_problem, :]):
                    logger.warning(
                        f"Problem '{problem_names_merged[i_problem]}' has merit_init = phi(x_0) = Inf at one or more runs. "
                        f"By convention, all solvers are declared to pass the convergence test for this problem at those runs."
                    )

        # Merge per-(problem, solver, run) diagnostic flags across the
        # problem libraries so they share the same problem index space as
        # ``problem_names_merged``. Default to all-False arrays for older
        # .h5 files that do not contain these fields.
        def _collect_flag(field_name):
            arrays = []
            for plib_dict in results_plibs:
                if field_name in plib_dict:
                    arrays.append(plib_dict[field_name])
                else:
                    n_p = plib_dict['fun_histories'].shape[0]
                    arrays.append(np.zeros((n_p, n_solvers, n_runs), dtype=bool))
            return np.concatenate(arrays, axis=0) if arrays else np.zeros((0, n_solvers, n_runs), dtype=bool)

        solver_abnormal_terminations_merged = _collect_flag('solver_abnormal_terminations')
        solver_output_fallbacks_merged = _collect_flag('solver_output_fallbacks')

        if is_saving:
            pdf_perf_hist_summary = backend_pdf.PdfPages(path_perf_hist_summary)
            pdf_perf_out_summary = backend_pdf.PdfPages(path_perf_out_summary)
            pdf_data_hist_summary = backend_pdf.PdfPages(path_data_hist_summary)
            pdf_data_out_summary = backend_pdf.PdfPages(path_data_out_summary)
            if n_solvers == 2:
                pdf_log_ratio_hist_summary = backend_pdf.PdfPages(path_log_ratio_hist_summary)
                pdf_log_ratio_out_summary = backend_pdf.PdfPages(path_log_ratio_out_summary)

        for i_tol, tolerance in enumerate(tolerances):
            hist = {
                'perf': [[None for _ in range(n_runs + 1)] for _ in range(n_solvers)],
                'data': [[None for _ in range(n_runs + 1)] for _ in range(n_solvers)],
            }
            if n_solvers == 2:
                hist['log_ratio'] = [None, None]
            # IMPORTANT: ``curve['hist']`` and ``curve['out']`` must be fully
            # independent containers. ``draw_profiles`` mutates the inner
            # ``perf`` / ``data`` / ``log_ratio`` lists in place when storing
            # the computed curves, and ``compute_scores`` later reads the two
            # branches separately to integrate history-based and output-based
            # profiles. A shallow copy (e.g. ``hist.copy()``) would leave the
            # nested lists shared between the two branches; the second
            # ``draw_profiles`` call would then overwrite the curves stored by
            # the first, so the history-based and output-based scores would be
            # silently identical. Use ``copy.deepcopy`` to break the sharing.
            curve = {
                'hist': hist,
                'out': copy.deepcopy(hist),
            }
            tolerance_str, tolerance_latex = format_float_scientific_latex(tolerance)
            if not profile_options[ProfileOption.SILENT]:
                logger.info(f'Creating profiles for tolerance {tolerance_str}')

            # Compute the number of function evaluations used by each
            # solver on each problem at each run to achieve convergence.
            work_hist = np.full((n_problems, n_solvers, n_runs), np.nan)
            work_out = np.full((n_problems, n_solvers, n_runs), np.nan)
            for i_problem in range(n_problems):
                for i_solver in range(n_solvers):
                    for i_run in range(n_runs):
                        if np.isinf(merit_inits_per_run[i_problem, i_run]):
                            # Degenerate case phi(x_0) = Inf. The
                            # Moré-Wild threshold
                            # tau*phi(x_0) + (1-tau)*phi_min is
                            # Inf - Inf = NaN, so the test cannot rank
                            # solvers. We declare every solver to pass
                            # by using +Inf as the threshold (any finite
                            # or infinite merit value <= +Inf). Both
                            # subcases "merit_min finite" and
                            # "merit_min = Inf" go through this branch
                            # and behave identically.
                            threshold = np.inf
                        elif np.isfinite(merit_mins_per_run[i_problem, i_run]):
                            threshold = max(tolerance * merit_inits_per_run[i_problem, i_run] + (1.0 - tolerance) * merit_mins_per_run[i_problem, i_run], merit_mins_per_run[i_problem, i_run])
                        else:
                            # Unreachable under the merit_min <=
                            # merit_init invariant established in
                            # ``process_results`` (since
                            # merit_inits_per_run is finite here). Kept
                            # defensively to preserve the previous
                            # "no convergence" behaviour.
                            threshold = -np.inf
                        if np.min(merit_histories_merged[i_problem, i_solver, i_run, :]) <= threshold:
                            work_hist[i_problem, i_solver, i_run] = np.argmax(merit_histories_merged[i_problem, i_solver, i_run, :] <= threshold) + 1
                        if merit_outs_merged[i_problem, i_solver, i_run] <= threshold:
                            work_out[i_problem, i_solver, i_run] = n_evals_merged[i_problem, i_solver, i_run]
                        
            for i_problem in range(n_problems):
                for i_run in range(n_runs):
                    solvers_all_diverge_hist[i_problem, i_run, i_tol] = np.all(np.isnan(work_hist[i_problem, :, i_run]))
                    solvers_all_diverge_out[i_problem, i_run, i_tol] = np.all(np.isnan(work_out[i_problem, :, i_run]))

            # Log if all solvers failed to meet convergence test.
            is_hist_drawable = np.any(~np.isnan(work_hist))
            is_out_drawable = np.any(~np.isnan(work_out))
            if not is_hist_drawable and not profile_options[ProfileOption.SILENT]:
                logger.info(f'All solvers failed to meet the convergence test for tolerance {tolerance_str} in history-based profiles.')
            if not is_out_drawable and not profile_options[ProfileOption.SILENT]:
                logger.info(f'All solvers failed to meet the convergence test for tolerance {tolerance_str} in output-based profiles.')

            # Draw the profiles.
            fig_perf_hist, fig_data_hist, fig_log_ratio_hist, curve['hist'] = draw_profiles(work_hist, problem_dims_merged, solver_names, tolerance_latex, i_tol, ax_summary_perf_hist, ax_summary_data_hist, ax_summary_log_ratio_hist, True, is_perf, is_data, is_log_ratio, profile_options, curve['hist'])
            fig_perf_out, fig_data_out, fig_log_ratio_out, curve['out'] = draw_profiles(work_out, problem_dims_merged, solver_names, tolerance_latex, i_tol, ax_summary_perf_out, ax_summary_data_out, ax_summary_log_ratio_out, is_output_based, is_perf, is_data, is_log_ratio, profile_options, curve['out'])
            curves.append(curve)

            # Save the profiles to files.
            if is_saving:
                if is_hist_drawable:
                    pdf_perf_hist = path_perf_hist / f'perf_hist_{i_tol + 1}.pdf'
                    fig_perf_hist.savefig(pdf_perf_hist, bbox_inches='tight')
                    pdf_perf_hist_summary.savefig(fig_perf_hist, bbox_inches='tight')

                    pdf_data_hist = path_data_hist / f'data_hist_{i_tol + 1}.pdf'
                    fig_data_hist.savefig(pdf_data_hist, bbox_inches='tight')
                    pdf_data_hist_summary.savefig(fig_data_hist, bbox_inches='tight')

                if is_out_drawable:
                    pdf_perf_out = path_perf_out / f'perf_out_{i_tol + 1}.pdf'
                    fig_perf_out.savefig(pdf_perf_out, bbox_inches='tight')
                    pdf_perf_out_summary.savefig(fig_perf_out, bbox_inches='tight')

                    pdf_data_out = path_data_out / f'data_out_{i_tol + 1}.pdf'
                    fig_data_out.savefig(pdf_data_out, bbox_inches='tight')
                    pdf_data_out_summary.savefig(fig_data_out, bbox_inches='tight')

                if n_solvers == 2:
                    if is_hist_drawable and fig_log_ratio_hist is not None:
                        pdf_log_ratio_hist = path_log_ratio_hist / f'log-ratio_hist_{i_tol + 1}.pdf'
                        fig_log_ratio_hist.savefig(pdf_log_ratio_hist, bbox_inches='tight')
                        pdf_log_ratio_hist_summary.savefig(fig_log_ratio_hist, bbox_inches='tight')
                    
                    if is_out_drawable and fig_log_ratio_out is not None:
                        pdf_log_ratio_out = path_log_ratio_out / f'log-ratio_out_{i_tol + 1}.pdf'
                        fig_log_ratio_out.savefig(pdf_log_ratio_out, bbox_inches='tight')
                        pdf_log_ratio_out_summary.savefig(fig_log_ratio_out, bbox_inches='tight')
                
            # Close the individual figures.
            plt.close(fig_perf_hist)
            plt.close(fig_perf_out)
            plt.close(fig_data_hist)
            plt.close(fig_data_out)
            if fig_log_ratio_hist is not None:
                plt.close(fig_log_ratio_hist)
            if fig_log_ratio_out is not None:
                plt.close(fig_log_ratio_out)
            
        # Close the summary pdf files.
        if is_saving:
            pdf_perf_hist_summary.close()
            pdf_perf_out_summary.close()
            pdf_data_hist_summary.close()
            pdf_data_out_summary.close()
            if n_solvers == 2:
                pdf_log_ratio_hist_summary.close()
                pdf_log_ratio_out_summary.close()
        
        if is_saving:
            try:
                if is_hist_drawable:
                    add_to_readme(path_readme_feature, 'perf_hist.pdf', 'File, the summary PDF of history-based performance profiles for all tolerances.')
                    add_to_readme(path_readme_feature, 'data_hist.pdf', 'File, the summary PDF of history-based data profiles for all tolerances.')
                    if n_solvers == 2 and fig_log_ratio_hist is not None:
                        add_to_readme(path_readme_feature, 'log-ratio_hist.pdf', 'File, the summary PDF of history-based log-ratio profiles for all tolerances.')
                if is_out_drawable:
                    add_to_readme(path_readme_feature, 'perf_out.pdf', 'File, the summary PDF of output-based performance profiles for all tolerances.')
                    add_to_readme(path_readme_feature, 'data_out.pdf', 'File, the summary PDF of output-based data profiles for all tolerances.')
                    if n_solvers == 2 and fig_log_ratio_out is not None:
                        add_to_readme(path_readme_feature, 'log-ratio_out.pdf', 'File, the summary PDF of output-based log-ratio profiles for all tolerances.')
            except Exception:
                pass
        
        # Record the names of the problems all the solvers failed to meet the convergence test for every tolerance.
        if is_saving:
            try:
                max_name_length = max(len(name) for name in problem_names_merged)
                with open(path_report, 'a') as fid:
                    fid.write('\n')
                    fid.write('## Problems among all the libraries that all the solvers failed to meet the convergence test for each tolerance and each run\n')
                    if np.any(solvers_all_diverge_hist) or np.any(solvers_all_diverge_out):
                        for i_tol in range(profile_options[ProfileOption.MAX_TOL_ORDER]):
                            tolerance = tolerances[i_tol]
                            tolerance_str = format_float_scientific_latex(tolerance, 1)[0]
                            for i_run in range(n_runs):
                                if np.any(solvers_all_diverge_hist[:, i_run, i_tol]):
                                    fid.write('\n')
                                    fid.write(f'History-based  tol = {tolerance_str:<8s} run = {i_run + 1:<3d}:\t\t')
                                    for i_problem in range(n_problems):
                                        if solvers_all_diverge_hist[i_problem, i_run, i_tol]:
                                            fid.write(f'{problem_names_merged[i_problem]:<{max_name_length}s} ')
                                if np.any(solvers_all_diverge_out[:, i_run, i_tol]):
                                    fid.write('\n')
                                    fid.write(f'Output-based   tol = {tolerance_str:<8s} run = {i_run + 1:<3d}:\t\t')
                                    for i_problem in range(n_problems):
                                        if solvers_all_diverge_out[i_problem, i_run, i_tol]:
                                            fid.write(f'{problem_names_merged[i_problem]:<{max_name_length}s} ')
                        fid.write('\n')
                    else:
                        fid.write('\n')
                        fid.write('This part is empty.\n')

                    # Record the (problem, run) pairs whose
                    # merit_init = phi(x_0) = Inf. For these the
                    # convergence test is degenerate and every solver is
                    # declared to pass; the section below lets the user
                    # identify and inspect them.
                    fid.write('\n')
                    fid.write('## Problems with merit_init = phi(x_0) = Inf (all solvers declared passing for these runs)\n')
                    if np.any(merit_init_inf_mask):
                        for i_run in range(n_runs):
                            if np.any(merit_init_inf_mask[:, i_run]):
                                fid.write('\n')
                                fid.write(f'run = {i_run + 1:<3d}:\t\t')
                                for i_problem in range(n_problems):
                                    if merit_init_inf_mask[i_problem, i_run]:
                                        fid.write(f'{problem_names_merged[i_problem]:<{max_name_length}s} ')
                        fid.write('\n')
                    else:
                        fid.write('\n')
                        fid.write('This part is empty.\n')

                    # Record (solver, problem, run) triples that
                    # terminated abnormally (solver raised an exception).
                    # Note: the evaluation history collected before the
                    # crash IS preserved and used for history-based
                    # profiles; only the output-based profile is
                    # affected (it falls back to the initial point as a
                    # penalty -- see the next section).
                    solver_names_no_esc = [s.replace('\\_', '_') for s in profile_options[ProfileOption.SOLVER_NAMES]]
                    fid.write('\n')
                    fid.write('## Solver runs that terminated abnormally (history is preserved; output uses x_0 as a penalty)\n')
                    if np.any(solver_abnormal_terminations_merged):
                        for i_solver in range(n_solvers):
                            for i_run in range(n_runs):
                                col = solver_abnormal_terminations_merged[:, i_solver, i_run]
                                if np.any(col):
                                    fid.write('\n')
                                    fid.write(f'solver = {solver_names_no_esc[i_solver]:<{len(max(solver_names_no_esc, key=len))}s}  run = {i_run + 1:<3d}:\t\t')
                                    for i_problem in range(n_problems):
                                        if col[i_problem]:
                                            fid.write(f'{problem_names_merged[i_problem]:<{max_name_length}s} ')
                        fid.write('\n')
                    else:
                        fid.write('\n')
                        fid.write('This part is empty.\n')

                    # Record (solver, problem, run) triples whose output
                    # was REPLACED by the initial point because the
                    # solver did not return a usable value (raised, or
                    # returned None / empty). This is an OUTPUT-BASED
                    # PENALTY, not a recovery of the best point from the
                    # evaluation history. See ``_solve_one_problem`` for
                    # the design rationale.
                    fid.write('\n')
                    fid.write('## Solver runs using the initial point as output fallback (output-based penalty)\n')
                    if np.any(solver_output_fallbacks_merged):
                        for i_solver in range(n_solvers):
                            for i_run in range(n_runs):
                                col = solver_output_fallbacks_merged[:, i_solver, i_run]
                                if np.any(col):
                                    fid.write('\n')
                                    fid.write(f'solver = {solver_names_no_esc[i_solver]:<{len(max(solver_names_no_esc, key=len))}s}  run = {i_run + 1:<3d}:\t\t')
                                    for i_problem in range(n_problems):
                                        if col[i_problem]:
                                            fid.write(f'{problem_names_merged[i_problem]:<{max_name_length}s} ')
                        fid.write('\n')
                    else:
                        fid.write('\n')
                        fid.write('This part is empty.\n')
            except Exception as exc:
                if not profile_options[ProfileOption.SILENT]:
                    logger.warning('Failed to record the problems that all the solvers failed to meet the convergence test.')
                    logger.warning(f'Error message: {shorten_log_message(exc)}')

    if not profile_options[ProfileOption.SILENT]:
        logger.info('')
        logger.info('All single profiles are created.')

    if is_saving and is_summary:
        if not profile_options[ProfileOption.SILENT]:
            logger.info('')
            logger.info('Start creating the summary PDF of all the profiles.')

        # Save the summary for the current feature.
        # Use rc_context to ensure LaTeX rendering is applied to the title
        with plt.rc_context(profile_context):
            fig_summary_hist.supylabel('History-based profiles', fontsize='xx-large', horizontalalignment='right')
            if fig_summary_out is not None:
                fig_summary_out.supylabel('Output-based profiles', fontsize='xx-large', horizontalalignment='right')
            fig_summary.suptitle(_format_feature_title(feature.name, profile_context), fontsize='xx-large', verticalalignment='bottom')
            path_summary = path_stamp / f'summary_{stamp}.pdf'
            fig_summary.savefig(path_summary, bbox_inches='tight')

        if not profile_options[ProfileOption.SILENT]:
            logger.info('The summary PDF of all the profiles is created.')

        plt.close(fig_summary)

    # Save curves to file.
    if is_saving:
        try:
            import pickle
            with open(path_log / 'curves.pkl', 'wb') as f:
                pickle.dump(curves, f)
            add_to_readme(path_readme_log, 'curves.pkl', 'File, storing the curves of the profiles.')
        except Exception as exc:
            if not profile_options[ProfileOption.SILENT]:
                logger.warning('Failed to save curves.')
                logger.warning(f'Error message: {shorten_log_message(exc)}')

    # Compute profile_scores from curves.
    profile_scores = compute_scores(curves, profile_options)

    # Save profile_scores to file.
    if is_saving:
        try:
            import pickle
            with open(path_log / 'profile_scores.pkl', 'wb') as f:
                pickle.dump(profile_scores, f)
            add_to_readme(path_readme_log, 'profile_scores.pkl', 'File, storing the scores of solvers on each profile.')
        except Exception as exc:
            if not profile_options[ProfileOption.SILENT]:
                logger.warning('Failed to save profile_scores.')
                logger.warning(f'Error message: {shorten_log_message(exc)}')

    # Compute solver_scores using score_fun.
    score_fun = profile_options[ProfileOption.SCORE_FUN]
    solver_scores = score_fun(profile_scores)

    # Append the solver scores to the experiment report file. Mirrors
    # MATLAB ``benchmark.m`` lines 1602-1612 so that the on-disk
    # ``report.txt`` is identical between the two implementations.
    if is_saving:
        try:
            with open(path_report, 'a') as fid:
                fid.write('\n')
                fid.write('## Scores of the solvers\n\n')
                max_solver_name_length = max(len(name) for name in solver_names)
                for i, name in enumerate(solver_names):
                    fid.write(f'{name:<{max_solver_name_length}s}:    {solver_scores[i]:.4f}\n')
        except Exception as exc:
            if not profile_options[ProfileOption.SILENT]:
                logger.warning('Failed to append the solver scores to the report file.')
                logger.warning(f'Error message: {shorten_log_message(exc)}')

    # Print solver scores.
    if not profile_options[ProfileOption.SILENT]:
        logger.info('')
        logger.info('Scores of the solvers:')
        max_solver_name_length = max(len(name) for name in solver_names)
        for i, name in enumerate(solver_names):
            format_info_str = f'{{:<{max_solver_name_length}}}:    {{:.4f}}'
            logger.info(format_info_str.format(name, solver_scores[i]))

    # Print summary of saved files.
    if not profile_options[ProfileOption.SILENT]:
        logger.info('')
        logger.info(f'Finished creating profiles with the "{feature.name}" feature.')
        if is_saving:
            logger.info('')
            logger.info('=' * 70)
            if is_summary:
                logger.info('')
                logger.info(f'Summary PDF of the profiles is saved as:')
                logger.info(f'{path_stamp / f"summary_{stamp}.pdf"}')
            logger.info('')
            logger.info(f'Single profiles are stored in:')
            logger.info(f'{path_stamp / "detailed_profiles"}')
            logger.info('')
            logger.info(f'Report of the experiment is saved as:')
            logger.info(f'{path_report}')
            logger.info('')
            logger.info(f'Check the README file for more information about the experiment and the results:')
            logger.info(f'{path_readme_feature}')
            logger.info('')
            logger.info('=' * 70)

    # Close the listener of the logger.
    if is_saving:
        listener.stop()
        for h in listener.handlers:
            if hasattr(h, 'close'):
                h.close()
        log_queue.close()
        log_queue.join_thread()

    return solver_scores, profile_scores, curves


def _solve_all_problems(solvers, plib, feature, problem_options, profile_options, is_plot, path_hist_plots, log_queue=None):
    """
    Solve all problems in plib satisfying problem_options using solvers in the solvers and stores the computing results.
    """
    logger = get_logger(__name__)
    results = []

    # Get satisfied problem names.
    option_select = problem_options.copy()
    for key in [ProblemOption.PLIBS.value, ProblemOption.PROBLEM_NAMES.value, ProblemOption.EXCLUDELIST.value, ProblemOption.CUSTOM_PROBLEM_LIBS_PATH.value]:
        option_select.pop(key, None)

    # Get the custom problem libraries path (if any).
    custom_problem_libs_path = problem_options.get(ProblemOption.CUSTOM_PROBLEM_LIBS_PATH, None)

    # Get the module file path for the problem library.
    module_file_path, module_name, _ = _get_problem_lib_module_path(plib, custom_problem_libs_path)

    # Import the problem library module.
    try:
        spec = importlib.util.spec_from_file_location(module_name, str(module_file_path))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        selector_name = plib + '_select'
        select = getattr(module, selector_name)
    except Exception as exc:
        logger.error(f'Error occurred while importing the problem library {plib}. Error message: {shorten_log_message(exc)}')
        raise exc

    # Select the problem names satisfying the options.
    try:
        if ProblemOption.PROBLEM_NAMES in problem_options:
            problem_names = problem_options[ProblemOption.PROBLEM_NAMES]
        else:
            problem_names = []
        if ProblemOption.EXCLUDELIST in problem_options:
            exclude_list = problem_options[ProblemOption.EXCLUDELIST]
        else:
            exclude_list = []
        selected_problem_names = select(option_select)
        if not problem_names:
            problem_names = selected_problem_names
        else:
            # Take an intersection of problem_names and selected_problem_names.
            problem_names = [name for name in problem_names if name in selected_problem_names]
            problem_names = list(set(problem_names))
        if exclude_list:
            problem_names = [name for name in problem_names if name not in exclude_list]
    except:
        pass

    if not problem_names:
        if not profile_options[ProfileOption.SILENT]:
            logger.info('')
            logger.info(f'No problem is selected from "{plib}".')
        return None

    n_problems = len(problem_names)
    len_problem_names = max(len(name) for name in problem_names)
    max_eval_factor = profile_options[ProfileOption.MAX_EVAL_FACTOR]
    profile_options_log, solver_aliases_used = _with_solver_log_names(profile_options, len_problem_names)
    if not profile_options[ProfileOption.SILENT]:
        logger.info('')
        logger.info(f'There are {n_problems} problems selected from "{plib}" to test.')
        if solver_aliases_used:
            _log_solver_aliases(logger, profile_options[ProfileOption.SOLVER_NAMES], profile_options_log[_SOLVER_LOG_NAMES_KEY])

    # Determine whether to use sequential mode or parallel mode.
    sequential_mode = profile_options[ProfileOption.N_JOBS] == 1 or os.cpu_count() == 1

    # Fall back to sequential mode if the arguments cannot be pickled (e.g.,
    # when the user passes lambda functions as solvers or feature modifiers).
    # Use ForkingPickler (same as multiprocessing.Pool) — plain pickle.dumps can
    # disagree and let unpicklable tasks reach starmap, which then crashes.
    if not sequential_mode:
        try:
            sample_arg = (
                solvers,
                feature,
                problem_names[0],
                len_problem_names,
                profile_options_log,
                is_plot,
                path_hist_plots,
                plib,
                custom_problem_libs_path,
            )
            ForkingPickler.dumps(sample_arg)
        except Exception:
            sequential_mode = True
            logger.info(
                'Falling back to sequential mode because worker processes cannot '
                'serialize the benchmark arguments (e.g., lambda functions in '
                'solvers, feature modifiers, or profile options). Use top-level '
                'functions (def ...) instead of lambdas to enable parallel execution.'
            )

    # Solve all problems.
    args = [(solvers, feature, problem_name, len_problem_names, profile_options_log, is_plot, path_hist_plots, plib, custom_problem_libs_path) for problem_name in problem_names]
    if sequential_mode:
        results = map(lambda arg: _solve_one_problem_wrapper(*arg), args)
    else:
        logger.info('Entering the parallel section.')
        with Pool(profile_options[ProfileOption.N_JOBS], initializer=setup_worker_logging,
        initargs=(log_queue,)) as p:
            results = p.starmap(_solve_one_problem_wrapper, args)
        logger.info('Leaving the parallel section.')

    # Delete result that is None in results
    results = [result for result in results if result is not None]
    n_problems = len(results)

    if n_problems == 0:
        logger.info('All problems from "{plib}" are not solved successfully.')
        return None

    # Process the results.
    n_solvers = len(solvers)
    n_runs = feature.options[FeatureOption.N_RUNS]
    
    problem_types = [r['problem_type'] for r in results]
    problem_dims = np.array([r['problem_dim'] for r in results])
    problem_mbs = np.array([r['problem_mb'] for r in results])
    problem_mlcons = np.array([r['problem_mlcon'] for r in results])
    problem_mnlcons = np.array([r['problem_mnlcon'] for r in results])
    problem_mcons = np.array([r['problem_mcon'] for r in results])
    fun_outs = np.array([r['fun_out'] for r in results])
    maxcv_outs = np.array([r['maxcv_out'] for r in results])
    fun_inits = np.array([r['fun_init'] for r in results])
    maxcv_inits = np.array([r['maxcv_init'] for r in results])
    n_evals = np.array([r['n_eval'] for r in results])
    problem_names = [r['problem_name'] for r in results]
    computation_times = np.array([r['computation_time'] for r in results])
    solvers_successes = np.array([r['solvers_success'] for r in results])
    # Per-(problem, solver, run) flags for the §6 report sections.
    # Default to False for backward compatibility with results that
    # predate this field (e.g. when loaded from old .h5 files).
    solver_abnormal_terminations = np.array([
        r.get('solver_abnormal_termination', np.zeros((n_solvers, n_runs), dtype=bool))
        for r in results
    ])
    solver_output_fallbacks = np.array([
        r.get('solver_output_fallback', np.zeros((n_solvers, n_runs), dtype=bool))
        for r in results
    ])

    if n_problems > 0:
        max_max_eval = int(np.ceil(max_eval_factor * np.max(problem_dims)))
    else:
        max_max_eval = 1

    fun_histories = np.full((n_problems, n_solvers, n_runs, max_max_eval), np.nan)
    maxcv_histories = np.full((n_problems, n_solvers, n_runs, max_max_eval), np.nan)

    for i_problem, r in enumerate(results):
        max_eval = int(np.ceil(max_eval_factor * problem_dims[i_problem]))
        fun_histories[i_problem, ..., :max_eval] = r['fun_history']
        maxcv_histories[i_problem, ..., :max_eval] = r['maxcv_history']
        if max_eval > 0 and max_max_eval > max_eval:
            last_fun = fun_histories[i_problem, ..., max_eval - 1, np.newaxis]
            fun_histories[i_problem, ..., max_eval:] = np.repeat(last_fun, max_max_eval - max_eval, axis=-1)
            last_maxcv = maxcv_histories[i_problem, ..., max_eval - 1, np.newaxis]
            maxcv_histories[i_problem, ..., max_eval:] = np.repeat(last_maxcv, max_max_eval - max_eval, axis=-1)

    results = {}
    results['plib'] = plib
    results['solver_names'] = profile_options[ProfileOption.SOLVER_NAMES]
    results['ptype'] = problem_options[ProblemOption.PTYPE]
    results['mindim'] = problem_options[ProblemOption.MINDIM]
    results['maxdim'] = problem_options[ProblemOption.MAXDIM]
    results['minb'] = problem_options[ProblemOption.MINB]
    results['maxb'] = problem_options[ProblemOption.MAXB]
    results['minlcon'] = problem_options[ProblemOption.MINLCON]
    results['maxlcon'] = problem_options[ProblemOption.MAXLCON]
    results['minnlcon'] = problem_options[ProblemOption.MINNLCON]
    results['maxnlcon'] = problem_options[ProblemOption.MAXNLCON]
    results['mincon'] = problem_options[ProblemOption.MINCON]
    results['maxcon'] = problem_options[ProblemOption.MAXCON]
    results['problem_names_options'] = problem_options[ProblemOption.PROBLEM_NAMES]
    results['excludelist'] = problem_options[ProblemOption.EXCLUDELIST]
    results['feature_stamp'] = profile_options[ProfileOption.FEATURE_STAMP]
    results['fun_histories'] = fun_histories
    results['maxcv_histories'] = maxcv_histories
    results['fun_outs'] = fun_outs
    results['maxcv_outs'] = maxcv_outs
    results['fun_inits'] = fun_inits
    results['maxcv_inits'] = maxcv_inits
    results['n_evals'] = n_evals
    results['problem_names'] = problem_names
    results['problem_types'] = problem_types
    results['problem_dims'] = problem_dims
    results['problem_mbs'] = problem_mbs
    results['problem_mlcons'] = problem_mlcons
    results['problem_mnlcons'] = problem_mnlcons
    results['problem_mcons'] = problem_mcons
    results['computation_times'] = computation_times
    results['solvers_successes'] = solvers_successes
    results['solver_abnormal_terminations'] = solver_abnormal_terminations
    results['solver_output_fallbacks'] = solver_output_fallbacks

    return results


def _solve_one_problem_wrapper(solvers, feature, problem_name, len_problem_names, profile_options, is_plot, path_hist_plots, plib, custom_problem_libs_path=None):
    logger = get_logger(__name__)
    
    # Get the module file path for the problem library.
    module_file_path, module_name, _ = _get_problem_lib_module_path(plib, custom_problem_libs_path)
    
    # Import the problem library module.
    try:
        spec = importlib.util.spec_from_file_location(module_name, str(module_file_path))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        loader_name = plib + '_load'
        load = getattr(module, loader_name)
    except Exception as exc:
        logger.error(f'Error occurred while importing the problem library {plib}. Error message: {shorten_log_message(exc)}')
        raise exc
    try:
        if not profile_options[ProfileOption.SILENT]:
            logger.info(f'Loading problem   {problem_name:<{len_problem_names}} from "{plib}".')
        problem = load(problem_name)
    except:
        if not profile_options[ProfileOption.SILENT]:
            logger.warning(f'Failed to load    {problem_name:<{len_problem_names}} from "{plib}".')
        return None
    result = _solve_one_problem(solvers, problem, feature, problem_name, len_problem_names, profile_options, is_plot, path_hist_plots)
    return result


def _safe_history_pdf_file_name(problem_name):
    return re.sub(r'^[-_]+', '', re.sub(r'[-_]+', '_', re.sub(r'[^a-zA-Z0-9\-_]', '', str(problem_name).replace(' ', '_')))) + '.pdf'


def _format_feature_title(feature_name, profile_context):
    feature_title = format_profile_text(feature_name, profile_context)
    return f'Profiles with "{feature_title}" feature'


def _format_problem_feature_title(problem_name, feature_stamp, profile_context):
    problem_title = format_profile_text(problem_name, profile_context)
    feature_title = format_profile_text(feature_stamp, profile_context)
    return f'Solving "{problem_title}" with "{feature_title}" feature'


def _history_title_fontsize(title_text, default_width):
    dpi = plt.rcParams.get('figure.dpi', 100)
    return min(12, 1.2 * default_width * dpi / max(len(title_text), 1))


def _save_problem_history_pdf(pdf_path, mode, problem_name, problem_type, problem_dim, solver_names,
                              fun_history, maxcv_history, merit_history, fun_init, maxcv_init, merit_init,
                              n_eval, profile_options, default_width, default_height, summary_profile_width,
                              profile_context):
    if mode not in {'combined', 'raw', 'cummin'}:
        raise ValueError(f'Unknown history plot mode: {mode}.')

    is_cum_rows = {'combined': [False, True], 'raw': [False], 'cummin': [True]}[mode]
    n_rows = len(is_cum_rows)
    n_cols = 1 if problem_type == 'u' else 3
    pdf_path = Path(pdf_path)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    with plt.rc_context(profile_context):
        fig_summary = plt.figure(figsize=(summary_profile_width * n_cols, default_height * n_rows))
        try:
            title_text = _format_problem_feature_title(problem_name, profile_options[ProfileOption.FEATURE_STAMP], profile_context)
            fig_summary.suptitle(title_text, fontsize=_history_title_fontsize(title_text, default_width))

            axs_summary = []
            for i_row in range(n_rows):
                for j_col in range(n_cols):
                    axs_summary.append(fig_summary.add_subplot(n_rows, n_cols, i_row * n_cols + j_col + 1))

            for i_row, is_cum in enumerate(is_cum_rows):
                row_start = i_row * n_cols
                if problem_type == 'u':
                    cell_axs_summary = [axs_summary[row_start]]
                else:
                    cell_axs_summary = axs_summary[row_start:row_start + n_cols]
                draw_hist(fun_history, maxcv_history, merit_history, fun_init, maxcv_init, merit_init,
                          solver_names, cell_axs_summary, is_cum, problem_type, problem_dim, n_eval,
                          profile_options, default_height)

            if mode == 'combined':
                fig_summary.tight_layout(rect=[0.05, 0.0, 1.0, 0.98])
                row1_y = 0.5 * (axs_summary[0].get_position().y0 + axs_summary[0].get_position().y1)
                row2_y = 0.5 * (axs_summary[-1].get_position().y0 + axs_summary[-1].get_position().y1)
                fig_summary.text(0.005, row1_y, "History profiles", rotation=90, va='center', fontsize=14)
                fig_summary.text(0.005, row2_y, "Cummin history profiles", rotation=90, va='center', fontsize=14)
            else:
                fig_summary.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])

            fig_summary.savefig(pdf_path, bbox_inches='tight')
        finally:
            plt.close(fig_summary)


def _export_problem_history_plots(problem_name, problem_type, problem_dim, solver_names,
                                  fun_history, maxcv_history, merit_history, fun_init, maxcv_init, merit_init,
                                  n_eval, profile_options, path_hist_plots):
    default_figsize = plt.rcParams['figure.figsize']
    default_width = default_figsize[0]
    default_height = default_figsize[1]
    summary_profile_width = default_width + summary_legend_extra_width(len(solver_names), default_width, solver_names)
    profile_context = set_profile_context(profile_options)
    processed_solver_names = [latex_escape_text(name) for name in solver_names]

    path_hist_plots = Path(path_hist_plots)
    pdf_file_name = _safe_history_pdf_file_name(problem_name)
    outputs = [
        (path_hist_plots / pdf_file_name, 'combined'),
        (path_hist_plots / 'raw' / pdf_file_name, 'raw'),
        (path_hist_plots / 'cummin' / pdf_file_name, 'cummin'),
    ]
    for pdf_path, mode in outputs:
        _save_problem_history_pdf(
            pdf_path, mode, problem_name, problem_type, problem_dim, processed_solver_names,
            fun_history, maxcv_history, merit_history, fun_init, maxcv_init, merit_init,
            n_eval, profile_options, default_width, default_height, summary_profile_width, profile_context
        )


def _solve_one_problem(solvers, problem, feature, problem_name, len_problem_names, profile_options, is_plot, path_hist_plots):
    """
    Solve a given problem.
    """
    solver_names = profile_options[ProfileOption.SOLVER_NAMES]
    solver_log_names = profile_options.get(_SOLVER_LOG_NAMES_KEY, solver_names)
    solver_isrand = profile_options[ProfileOption.SOLVER_ISRAND]
    result = None

    if type(problem) is not Problem:
        raise TypeError('The problem provided is not a valid Problem object.')

    # Project the initial point if necessary.
    if profile_options[ProfileOption.PROJECT_X0]:
        problem.project_x0()

    # Solve the problem with each solver.
    n_solvers = len(solvers)
    n_runs = feature.options[FeatureOption.N_RUNS]
    max_eval = profile_options[ProfileOption.MAX_EVAL_FACTOR] * problem.n
    max_eval = int(np.ceil(max_eval))
    n_eval = np.zeros((n_solvers, n_runs), dtype=int)
    fun_history = np.full((n_solvers, n_runs, max_eval), np.nan)
    fun_out = np.full((n_solvers, n_runs), np.nan)
    maxcv_history = np.full((n_solvers, n_runs, max_eval), np.nan)
    maxcv_out = np.full((n_solvers, n_runs), np.nan)
    # ``fun_inits`` / ``maxcv_inits`` are vectors of length ``n_runs``: the
    # initial values must be recorded per-run because the feature can
    # modify the initial guess differently at every run (e.g., the
    # ``perturbed_x0`` and ``noisy`` features), and using a single
    # problem-level value would silently misrepresent the initial point
    # actually seen by the solver. They are filled inside the run loop
    # below using ``featured_problem.fun_init`` / ``maxcv_init`` (which
    # apply the feature) instead of ``problem.fun(problem.x0)`` (which
    # would ignore the feature). This mirrors MATLAB's
    # ``solveOneProblem.m`` (lines 29-30, 63-64).
    fun_inits = np.full(n_runs, np.nan)
    maxcv_inits = np.full(n_runs, np.nan)
    computation_time = np.full((n_solvers, n_runs), np.nan)
    solvers_success = np.full((n_solvers, n_runs), False)
    # ``solver_abnormal_terminations[i_solver, i_run]`` is True when the
    # solver call raised an exception; ``solver_output_fallbacks`` is
    # True when the output ``x`` was reset to the initial guess because
    # the solver did not return a usable point (either because it raised
    # or because it returned ``None`` / an empty array). These two flags
    # power the dedicated report sections written in ``benchmark`` so the
    # user can see which (problem, solver, run) triples were affected by
    # the OUTPUT-BASED PENALTY mechanism documented in
    # ``_solve_one_problem`` above.
    solver_abnormal_terminations = np.full((n_solvers, n_runs), False)
    solver_output_fallbacks = np.full((n_solvers, n_runs), False)

    # The number of real runs for each solver, which is determined by feature and solver_isrand.
    real_n_runs = np.array([
        n_runs if feature.is_stochastic or (solver_isrand is not None and solver_isrand[i_solver]) else 1
        for i_solver in range(n_solvers)
    ], dtype=int)
    
    len_solver_log_names = max(len(name) for name in solver_log_names)

    logger = get_logger(__name__)

    # Solve the problem with each solver.
    for i_solver in range(n_solvers):
        for i_run in range(real_n_runs[i_solver]):
            if not profile_options[ProfileOption.SILENT]:
                logger.info(
                    f"Start  solving    {problem_name:<{len_problem_names}} with {solver_log_names[i_solver]:<{len_solver_log_names}} (run {i_run + 1:2d}/{real_n_runs[i_solver]:2d})."
                )

            # Construct the featured problem.
            real_seed = (23333 * profile_options[ProfileOption.SEED] + 211 * i_run) % (2**32)
            featured_problem = FeaturedProblem(problem, feature, max_eval, real_seed)

            # Save the problem information.
            problem_type = featured_problem.ptype
            problem_dim = featured_problem.n
            problem_mb = featured_problem.mb
            problem_mlcon = featured_problem.mlcon
            problem_mnlcon = featured_problem.mnlcon
            problem_mcon = featured_problem.mcon
            # Record the per-run initial values *after* the feature is
            # applied. ``real_seed`` depends only on ``i_run``, so the
            # FeaturedProblem (and hence ``fun_init`` / ``maxcv_init``)
            # is identical regardless of which solver's loop sets it
            # for a given run; assigning here keeps the loop structure
            # aligned with MATLAB's ``solveOneProblem.m``.
            fun_inits[i_run] = featured_problem.fun_init
            maxcv_inits[i_run] = featured_problem.maxcv_init

            # Define a function to call the solver.
            def call_solver():
                if problem_type == 'u':
                    return solvers[i_solver](featured_problem.fun, featured_problem.x0)
                elif problem_type == 'b':
                    return solvers[i_solver](featured_problem.fun, featured_problem.x0, featured_problem.xl, featured_problem.xu)
                elif problem_type == 'l':
                    return solvers[i_solver](featured_problem.fun, featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq)
                elif problem_type == 'n':
                    return solvers[i_solver](featured_problem.fun, featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, featured_problem.cub, featured_problem.ceq)              

            # Solve the problem with the solver.
            #
            # The control flow below mirrors the MATLAB implementation in
            # ``solveOneProblem.m`` and is built around two nested
            # ``try`` blocks:
            #
            # * The INNER ``try`` wraps only the solver call. If the solver
            #   raises, ``x`` is *not* reassigned and stays as the pre-set
            #   initial guess ``featured_problem.x0`` (assigned just below
            #   for exactly this reason). This fallback to ``x0`` is an
            #   intentional OUTPUT-BASED PENALTY for failed / crashed runs:
            #   the output-based profile (and the corresponding score) is
            #   then evaluated at the initial point, while the
            #   history-based profile is unaffected because it is built
            #   from ``featured_problem.fun_hist`` /
            #   ``featured_problem.maxcv_hist`` which already record every
            #   evaluation the solver actually performed before crashing.
            #   We deliberately do NOT try to "recover" the best evaluated
            #   point from the history for the output-based score: the
            #   user only ever sees the point a solver returns, so a
            #   crashed run must be punished there. Falling back to ``x0``
            #   is the most basic choice that is solver-independent and
            #   guarantees a finite-or-Inf value derived from the problem
            #   itself, not from the failing solver.
            #
            # * The OUTER ``try`` wraps the post-solve processing
            #   (back-transform of ``x``, evaluation of ``fun_out`` /
            #   ``maxcv_out``, logging). If anything in there raises we
            #   just log it; ``fun_out`` / ``maxcv_out`` then stay at NaN
            #   so downstream code treats them as missing data.
            #
            # Note: ``featured_problem.x0`` already returns a fresh
            # ``np.copy`` of the modified initial guess, so this assignment
            # is safe even if the solver mutates its input array.
            x = featured_problem.x0
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                time_start_solver_run = time.monotonic()
                try:
                    try:
                        if profile_options[ProfileOption.SOLVER_VERBOSE] == 2:
                            x = call_solver()
                        elif profile_options[ProfileOption.SOLVER_VERBOSE] == 1:
                            with open(os.devnull, 'w') as devnull, redirect_stdout(devnull):
                                x = call_solver()
                        else:
                            with open(os.devnull, 'w') as devnull, redirect_stdout(devnull), redirect_stderr(devnull):
                                x = call_solver()
                    except Exception as exc:
                        # Mark this (solver, run) as abnormally terminated;
                        # ``x`` keeps its pre-assigned ``featured_problem.x0``
                        # value and the fallback flag will also be set
                        # below (since ``x is featured_problem.x0`` here).
                        solver_abnormal_terminations[i_solver, i_run] = True
                        solver_output_fallbacks[i_solver, i_run] = True
                        if profile_options[ProfileOption.SOLVER_VERBOSE] >= 1:
                            logger.warning(
                                f'An error occurred while solving {problem_name} with {solver_log_names[i_solver]} '
                                f'(run {i_run + 1}/{real_n_runs[i_solver]}): {_shorten_log_message(exc)}'
                            )

                    # Mirror MATLAB's ``if isempty(x)`` guard: when the
                    # solver did not raise but returned ``None`` or an
                    # empty array, still fall back to the initial guess
                    # (same OUTPUT-BASED PENALTY as the crash case above).
                    x_arr = np.asarray(x) if x is not None else None
                    if x_arr is None or x_arr.size == 0:
                        x = featured_problem.x0
                        solver_output_fallbacks[i_solver, i_run] = True

                    # Record the END timestamp explicitly and take the
                    # straightforward float difference. Storing both
                    # endpoints (as opposed to wrapping a tic/toc-style
                    # opaque handle) means the elapsed time is just an
                    # IEEE-754 subtraction and the raw start / end
                    # values are available for forensics if either of
                    # the defensive branches below ever fires.
                    #
                    # ``time.monotonic`` is guaranteed by PEP 418 not
                    # to go backward (see
                    # https://docs.python.org/3/library/time.html#time.monotonic
                    # and https://peps.python.org/pep-0418/ ), so under
                    # normal operation neither branch below is reached.
                    # We keep them as defensive guards against
                    # pathological cases (custom solver tampering with
                    # the timer, mocked time during tests, etc.) and
                    # to keep the control flow in lockstep with the
                    # MATLAB implementation in ``solveOneProblem.m``.
                    time_end_monotonic = time.monotonic()
                    elapsed = time_end_monotonic - time_start_solver_run

                    if elapsed < 0:
                        if not profile_options[ProfileOption.SILENT]:
                            logger.warning(
                                f'Negative elapsed time for {problem_name} with {solver_log_names[i_solver]} '
                                f'(run {i_run + 1}/{real_n_runs[i_solver]}): '
                                f't_start={time_start_solver_run:.6f} t_end={time_end_monotonic:.6f} '
                                f'diff={elapsed:.6g} s (clamped to 0; time.monotonic violates PEP 418).'
                            )
                        elapsed = 0.0
                    elif not np.isfinite(elapsed) or elapsed > 1e7:
                        if not profile_options[ProfileOption.SILENT]:
                            logger.warning(
                                f'Implausibly large elapsed time {elapsed:.6g} s for {problem_name} with '
                                f'{solver_log_names[i_solver]} (run {i_run + 1}/{real_n_runs[i_solver]}): '
                                f't_start={time_start_solver_run:.6f} t_end={time_end_monotonic:.6f} '
                                '(recorded as NaN; likely a timer glitch, not a real run).'
                            )
                        elapsed = np.nan
                    computation_time[i_solver, i_run] = elapsed

                    # It is very important to transform the solution back to the one related to the original problem. (Note that the problem we solve has the objective function f(A @ x + b). Thus, if x is the output solution, then A @ x + b is the solution of the original problem.)
                    A, b = featured_problem._feature.modifier_affine(featured_problem._seed, featured_problem._problem)[:2]
                    x = A @ x + b

                    # Use problem.fun and problem.maxcv to evaluate the solution since it is possible that featured_problem.fun and featured_problem.maxcv are modified.
                    fun_out[i_solver, i_run] = problem.fun(x)
                    maxcv_out[i_solver, i_run] = problem.maxcv(x)
                    # Calculate the minimum function value and the minimum constraint violation, omitting the NaN values.
                    fun_hist = featured_problem.fun_hist
                    maxcv_hist = featured_problem.maxcv_hist
                    fun_min = np.nan if fun_hist.size == 0 else np.nanmin(fun_hist)
                    maxcv_min = np.nan if maxcv_hist.size == 0 else np.nanmin(maxcv_hist)

                    if not profile_options[ProfileOption.SILENT]:
                        logger.info(
                            f"Finish solving    {problem_name:<{len_problem_names}} with {solver_log_names[i_solver]:<{len_solver_log_names}} "
                            f"(run {i_run + 1:2d}/{real_n_runs[i_solver]:2d}) (in {computation_time[i_solver, i_run]:.2f} seconds)."
                        )

                        if problem_type == 'u':
                            logger.info(
                                f"Output result for {problem_name:<{len_problem_names}} with {solver_log_names[i_solver]:<{len_solver_log_names}} "
                                f"(run {i_run + 1:2d}/{real_n_runs[i_solver]:2d}): f = {fun_out[i_solver, i_run]:10.4e}."
                            )
                            logger.info(
                                f"Best   result for {problem_name:<{len_problem_names}} with {solver_log_names[i_solver]:<{len_solver_log_names}} "
                                f"(run {i_run + 1:2d}/{real_n_runs[i_solver]:2d}): f = {fun_min:10.4e}."
                            )
                        else:
                            logger.info(
                                f"Output result for {problem_name:<{len_problem_names}} with {solver_log_names[i_solver]:<{len_solver_log_names}} "
                                f"(run {i_run + 1:2d}/{real_n_runs[i_solver]:2d}): f = {fun_out[i_solver, i_run]:10.4e}, maxcv = {maxcv_out[i_solver, i_run]:10.4e}."
                            )
                            logger.info(
                                f"Best   result for {problem_name:<{len_problem_names}} with {solver_log_names[i_solver]:<{len_solver_log_names}} "
                                f"(run {i_run + 1:2d}/{real_n_runs[i_solver]:2d}): f = {fun_min:10.4e}, maxcv = {maxcv_min:10.4e}."
                            )
                except Exception as exc:
                    if profile_options[ProfileOption.SOLVER_VERBOSE] >= 1:
                        logger.warning(
                            f'An error occurred while processing the solution of {problem_name} from {solver_log_names[i_solver]} '
                            f'(run {i_run + 1}/{real_n_runs[i_solver]}): {_shorten_log_message(exc)}'
                        )
            n_eval[i_solver, i_run] = featured_problem.n_eval_fun
            fun_history[i_solver, i_run, :n_eval[i_solver, i_run]] = featured_problem.fun_hist[:n_eval[i_solver, i_run]]
            maxcv_history[i_solver, i_run, :n_eval[i_solver, i_run]] = featured_problem.maxcv_hist[:n_eval[i_solver, i_run]]
            if n_eval[i_solver, i_run] > 0:
                fun_history[i_solver, i_run, n_eval[i_solver, i_run]:] = fun_history[i_solver, i_run, n_eval[i_solver, i_run] - 1]
                maxcv_history[i_solver, i_run, n_eval[i_solver, i_run]:] = maxcv_history[i_solver, i_run, n_eval[i_solver, i_run] - 1]
                solvers_success[i_solver, i_run] = True
            
        
        # If real_n_runs(i_solver) == 1 ~= n_runs, then we need to copy the result to the other runs.
        if real_n_runs[i_solver] == 1 and n_runs > 1:
            for j in range(1, n_runs):
                n_eval[i_solver, j] = n_eval[i_solver, 0]
                fun_history[i_solver, j, :] = fun_history[i_solver, 0, :]
                maxcv_history[i_solver, j, :] = maxcv_history[i_solver, 0, :]
                fun_out[i_solver, j] = fun_out[i_solver, 0]
                maxcv_out[i_solver, j] = maxcv_out[i_solver, 0]
                solver_abnormal_terminations[i_solver, j] = solver_abnormal_terminations[i_solver, 0]
                solver_output_fallbacks[i_solver, j] = solver_output_fallbacks[i_solver, 0]

    # Complete fun_inits / maxcv_inits when real_n_runs[i_solver] < n_runs
    # for *every* solver (which forces real_n_runs[i_solver] == 1 for all
    # solvers, since otherwise some solver's run loop would have set every
    # run). In that degenerate case only fun_inits[0] / maxcv_inits[0] is
    # populated, so we replicate them to the remaining runs to keep the
    # output shape consistent. Mirrors MATLAB's solveOneProblem.m
    # (lines 191-198).
    if int(np.max(real_n_runs)) < n_runs:
        for i_run in range(1, n_runs):
            fun_inits[i_run] = fun_inits[0]
            maxcv_inits[i_run] = maxcv_inits[0]

    # Return the result.
    result = {
        'fun_history': fun_history,
        'maxcv_history': maxcv_history,
        'fun_out': fun_out,
        'maxcv_out': maxcv_out,
        'fun_init': fun_inits,
        'maxcv_init': maxcv_inits,
        'n_eval': n_eval,
        'problem_name': problem_name,
        'problem_type': problem_type,
        'problem_dim': problem_dim,
        'problem_mb': problem_mb,
        'problem_mlcon': problem_mlcon,
        'problem_mnlcon': problem_mnlcon,
        'problem_mcon': problem_mcon,
        'computation_time': computation_time,
        'solvers_success': solvers_success,
        'solver_abnormal_termination': solver_abnormal_terminations,
        'solver_output_fallback': solver_output_fallbacks,
    }

    # Draw the history plots if required.
    if not is_plot or path_hist_plots is None or profile_options[ProfileOption.SCORE_ONLY] or all(not success for success in solvers_success.flatten()):
        return result
    
    try:
        merit_fun = profile_options[ProfileOption.MERIT_FUN]
        try:
            # ``fun_inits`` / ``maxcv_inits`` are per-run vectors of shape
            # ``(n_runs,)`` (see the per-run initialization block above).
            merit_history = compute_merit_values(merit_fun, fun_history, maxcv_history, maxcv_inits)
            merit_init = compute_merit_values(merit_fun, fun_inits, maxcv_inits, maxcv_inits)
        except Exception as exc:
            logger.error(f'Error occurred while calculating the merit values. Please check the merit function. Error message: {shorten_log_message(exc)}')
            raise exc

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            _export_problem_history_plots(
                problem_name, problem_type, problem_dim, solver_names,
                fun_history, maxcv_history, merit_history, fun_inits, maxcv_inits, merit_init,
                n_eval, profile_options, path_hist_plots
            )

    except Exception as exc:
        if not profile_options[ProfileOption.SILENT]:
            logger.info(f'An error occurred while plotting the history plots of the problem {problem_name}.')
            logger.info(f'Error message: {shorten_log_message(exc)}')
        pass
    
    return result


def _draw_problem_history_plot(problem_name, problem_type, problem_dim, solver_names, solvers_success,
                                fun_history, maxcv_history, fun_init, maxcv_init, n_eval, 
                                profile_options, path_hist_plots):
    """
    Draw history plot for a single problem (used in sequential mode).
    
    Parameters
    ----------
    problem_name : str
        Name of the problem.
    problem_type : str
        Type of the problem ('u' for unconstrained, etc.).
    problem_dim : int
        Dimension of the problem.
    solver_names : list
        Names of the solvers.
    solvers_success : numpy.ndarray
        Boolean array indicating solver success.
    fun_history : numpy.ndarray
        History of function values.
    maxcv_history : numpy.ndarray
        History of maximum constraint violations.
    fun_init : numpy.ndarray
        Initial function values.
    maxcv_init : numpy.ndarray
        Initial maximum constraint violations.
    n_eval : numpy.ndarray
        Number of function evaluations.
    profile_options : dict
        Profile options.
    path_hist_plots : Path
        Path to save the history plots.
    """
    logger = get_logger(__name__)
    
    # Skip if history plots are disabled or no solver succeeded.
    if path_hist_plots is None or all(not success for success in solvers_success.flatten()):
        return
    
    try:
        merit_fun = profile_options[ProfileOption.MERIT_FUN]
        try:
            merit_history = compute_merit_values(merit_fun, fun_history, maxcv_history, maxcv_init)
            merit_init = compute_merit_values(merit_fun, fun_init, maxcv_init, maxcv_init)
        except Exception as exc:
            logger.error(f'Error occurred while calculating the merit values. Please check the merit function. Error message: {shorten_log_message(exc)}')
            raise exc

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            _export_problem_history_plots(
                problem_name, problem_type, problem_dim, solver_names,
                fun_history, maxcv_history, merit_history, fun_init, maxcv_init, merit_init,
                n_eval, profile_options, path_hist_plots
            )

    except Exception as exc:
        if not profile_options[ProfileOption.SILENT]:
            logger.info(f'An error occurred while plotting the history plots of the problem {problem_name}.')
            logger.info(f'Error message: {shorten_log_message(exc)}')
        pass


def _get_problem_lib_module_path(plib, custom_problem_libs_path=None):
    """
    Get the module file path and module name for a problem library.
    
    This function checks if the problem library is a built-in library (in the 
    optiprofiler/problem_libs/ directory) or a custom library (in the 
    custom_problem_libs_path directory).
    
    Parameters
    ----------
    plib : str
        The name of the problem library.
    custom_problem_libs_path : Path or None
        The path to the directory containing custom problem libraries.
        
    Returns
    -------
    module_file_path : Path
        The path to the problem library's tools module file.
    module_name : str
        The module name for importing.
    is_builtin : bool
        True if the problem library is a built-in library, False otherwise.
    """
    current_dir = Path(__file__).parent.resolve()
    builtin_module_file_path = current_dir / 'problem_libs' / plib / f'{plib}_tools.py'

    # Check explicit custom libraries before built-ins. This allows users to
    # test or override a built-in adapter by passing custom_problem_libs_path
    # with a library of the same name.
    if custom_problem_libs_path is not None:
        custom_problem_libs_path = Path(custom_problem_libs_path).expanduser().resolve()
        custom_library_dir = _custom_problem_library_dir(custom_problem_libs_path, plib)
        if custom_library_dir is not None:
            custom_module_file_path = custom_library_dir / f'{plib}_tools.py'
            if custom_module_file_path.is_file():
                module_file_path = custom_module_file_path.resolve()
                module_name = f'{plib}.{plib}_tools'
                return module_file_path, module_name, False
            raise ValueError(
                f'The custom problem library "{plib}" at {custom_library_dir} must '
                f'contain "{plib}_tools.py".'
            )

    # Check if it's a built-in library.
    if builtin_module_file_path.exists():
        module_file_path = builtin_module_file_path.resolve()
        module_name = f'optiprofiler.problem_libs.{plib}.{plib}_tools'
        return module_file_path, module_name, True
    
    # If not found in either location, return the builtin path (will fail later with a clear error).
    return builtin_module_file_path.resolve(), f'optiprofiler.problem_libs.{plib}.{plib}_tools', True





    

            
        
