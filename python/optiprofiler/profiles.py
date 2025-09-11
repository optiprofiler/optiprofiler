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
from pathlib import Path

import numpy as np
from cycler import cycler
from matplotlib import pyplot as plt
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D
from matplotlib.backends import backend_pdf
from matplotlib.ticker import MaxNLocator, FuncFormatter

from .features import Feature
from .problems import Problem, FeaturedProblem
from .utils import FeatureName, ProfileOption, FeatureOption, ProblemOption, get_logger, setup_main_process_logging, setup_worker_logging
from .loader import load_results
from .profile_utils import check_validity_problem_options, check_validity_profile_options, get_default_problem_options, get_default_profile_options, compute_merit_values, create_stamp, merge_pdfs_with_pypdf, write_report, process_results
from .plotting import draw_hist, set_profile_context, format_float_scientific_latex, draw_profiles


def benchmark(
    solvers: list[callable] | None = None,
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
    solvers : list of callable, optional if 'load' in ``**kwargs``
        Solvers to benchmark. Each solver must be a callable, as follows. For
        unconstrained problems, the signature of the callable can be

            ``solver(fun, x0) -> numpy.ndarray, shape (n,)``

        where ``fun`` is the objective function and ``x0`` is the initial
        point. The solver must return a numpy array as the solution. For
        bound constrained problems, the signature of the callable can be

            ``solver(fun, x0, lb, ub) -> numpy.ndarray, shape (n,)``

        where ``lb`` and ``ub`` are the lower and upper bounds, respectively.
        For linearly constrained problems, the signature of the callable can be

            ``solver(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq) -> numpy.ndarray, shape (n,)``

        where ``a_ub @ x <= b_ub`` and ``a_eq @ x == b_eq`` form the linear
        inequality and equality constraints, respectively. For nonlinearly
        constrained problems, the signature of the callable can be

            ``solver(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq, c_ub, c_eq) -> numpy.ndarray, shape (n,)``

        where ``c_ub(x) <= 0`` and ``c_eq(x) == 0`` form the nonlinear
        inequality and equality constraints, respectively. In all cases, the
        signature of the callable can also be

            ``solver(problem) -> numpy.ndarray, shape (n,)``

        where ``problem`` is an instance of `Problem`.

        All vectors and matrices mentioned above are `numpy.ndarray`.

        If the 'load' option is provided in ``**kwargs``, solvers can be None,
        in which case data from a previous experiment will be loaded to generate profiles.

    Other Parameters
    ----------------
    **Options for features:**

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
        feature or noise in 'noisy' feature. It should be either a str
        (or char), or a callable
            (random_stream, dimension) -> random vector,
        accepting a random_stream and the dimension of a problem and
        returning a random vector with the given dimension. In 'perturbed_x0'
        case, the str should be either 'spherical' or 'gaussian' (default is
        'spherical'). In 'noisy' case, the str should be either 'gaussian'
        or 'uniform' (default is 'gaussian').
    perturbation_level : float, optional
        The magnitude of the perturbation to the initial
        guess in the 'perturbed_x0' feature. Default is 1e-3.
    noise_level : float, optional
        The magnitude of the noise in the 'noisy' feature.
        Default is 1e-3.
    noise_type : str, optional
        The type of the noise in the 'noisy' features. It should
        be either 'absolute', 'relative', or 'mixed'. Default is 'mixed'.
    significant_digits : int, optional
        The number of significant digits in the
        'truncated' feature. Default is 6.
    perturbed_trailing_digits : bool, optional
        Whether we will randomize the trailing
        zeros of the objective function value in the 'perturbed_x0' feature.
        Default is False.
    rotated : bool, optional
        Whether to use a random or given rotation matrix to rotate
        the coordinates of a problem in the 'linearly_transformed' feature.
        Default is True.
    condition_factor : float, optional
        The scaling factor of the condition number of the
        linear transformation in the 'linearly_transformed' feature. More
        specifically, the condition number of the linear transformation will
        be 2 ** (condition_factor * n / 2), where n is the dimension of the
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
        'custom' feature. It should be a callable as follows:
            (random_stream, problem) -> modified_x0,
        where problem is an instance of the class Problem, and
        modified_x0 is the modified initial guess. No default.
    mod_affine : callable, optional
        The modifier function to generate the affine
        transformation applied to the variables in the 'custom' feature. It
        should be a callable as follows:
            (random_stream, problem) -> (A, b, inv),
        where problem is an instance of the class Problem, A is the
        matrix of the affine transformation, b is the vector of the affine
        transformation, and inv is the inverse of matrix A. No default.
    mod_bounds : callable, optional
        The modifier function to modify the bound constraints in
        the 'custom' feature. It should be a callable as follows:
            (random_stream, problem) -> (modified_xl, modified_xu),
        where problem is an instance of the class Problem, modified_xl is
        the modified lower bound, and modified_xu is the modified upper
        bound. No default.
    mod_linear_ub : callable, optional
        The modifier function to modify the linear inequality
        constraints in the 'custom' feature. It should be a callable as follows:
            (random_stream, problem) -> (modified_aub, modified_bub),
        where problem is an instance of the class Problem, modified_aub
        is the modified matrix of the linear inequality constraints, and
        modified_bub is the modified vector of the linear inequality
        constraints. No default.
    mod_linear_eq : callable, optional
        The modifier function to modify the linear equality
        constraints in the 'custom' feature. It should be a callable as follows:
            (random_stream, problem) -> (modified_aeq, modified_beq),
        where problem is an instance of the class Problem, modified_aeq
        is the modified matrix of the linear equality constraints, and
        modified_beq is the modified vector of the linear equality
        constraints. No default.
    mod_fun : callable, optional
        The modifier function to modify the objective function in
        the 'custom' feature. It should be a callable as follows:
            (x, random_stream, problem) -> modified_fun,
        where x is the evaluation point, problem is an instance of the
        class Problem, and modified_fun is the modified objective function
        value. No default.
    mod_cub : callable, optional
        The modifier function to modify the nonlinear inequality
        constraints in the 'custom' feature. It should be a callable as follows:
            (x, random_stream, problem) -> modified_cub,
        where x is the evaluation point, problem is an instance of the
        class Problem, and modified_cub is the modified vector of the
        nonlinear inequality constraints. No default.
    mod_ceq : callable, optional
        The modifier function to modify the nonlinear equality
        constraints in the 'custom' feature. It should be a callable as follows:
            (x, random_stream, problem) -> modified_ceq,
        where x is the evaluation point, problem is an instance of the
        class Problem, and modified_ceq is the modified vector of the
        nonlinear equality constraints. No default.

    **Options for profiles and plots:**

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
    line_colors : list or numpy.ndarray, optional
        The colors of the lines in the plots. It can be a list
        of short names of colors ('r', 'g', 'b', 'c', 'm', 'y', 'k') or
        an array with each row being a RGB triplet. Default line colors are
        those in the palettename named "gem" (see MATLAB documentation for
        'colororder'). Note that if the number of solvers is greater than the
        number of colors, we will cycle through the colors.
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
        the tolerances as 10**(-1:-1:-max_tol_order). Default is 10.
    merit_fun : callable, optional
        The merit function to measure the quality of a point using
        the objective function value and the maximum constraint violation.
        It should be a callable as follows:
            (fun_value, maxcv_value, maxcv_init) -> merit_value,
        where fun_value is the objective function value, maxcv_value is
        the maximum constraint violation, and maxcv_init is the maximum
        constraint violation at the initial guess. The default merit function
        varphi(x) is defined by the objective function f(x) and the maximum
        constraint violation v(x) as
          varphi(x) = f(x)                        if v(x) <= v1
          varphi(x) = f(x) + 1e5 * (v(x) - v1)   if v1 < v(x) <= v2
          varphi(x) = np.inf                      if v(x) > v2
        where v1 = max(1e-5, v0) and v2 = min(0.01, 1e-10 * max(1, v0)),
        and v0 is the maximum constraint violation at the initial guess.
    n_jobs : int, optional
        The number of parallel jobs to run the test. Default is the
        default number of workers in the default local cluster.
    normalized_scores : bool, optional
        Whether to normalize the scores of the solvers by
        the maximum score of the solvers. Default is False.
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
        solvers. It should be a callable as follows:
            profile_scores -> solver_scores,
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
        representing a nonnegative function in R^+. Default is lambda x: 1.
    seed : int, optional
        The seed of the random number generator. Default is 1.
    semilogx : bool, optional
        Whether to use the semilogx scale during plotting profiles
        (performance profiles and data profiles). Default is True.
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
        from 1 to the total number of solvers of the loading experiment. At
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

    **Options for problems:**

    Options in this part are used to select problems for benchmarking.
    First select which problem libraries to use based on the 'plibs'
    option. Then select problems from these libraries according to the
    given options ('problem_names', 'ptype', 'mindim', 'maxdim', 'minb',
    'maxb', 'minlcon', 'maxlcon', 'minnlcon', 'maxnlcon', 'mincon',
    'maxcon', and 'excludelist').

    plibs : list of str, optional
        The problem libraries to be used. It should be a list of
        strs. The available choices are subfolder names in the
        'problems' directory. There are three subfolders after installing the
        package: 's2mpj', 'matcutest', and 'custom'. Default setting is
        's2mpj'.
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

    Notes
    -----
    The current version supports benchmarking derivative-free optimization solvers.

    .. caution::

        The log-ratio profiles are available only when there are exactly two
        solvers. For more information of performance and data profiles, see
        [1]_, [2]_, [5]_. For that of log-ratio profiles, see [4]_, [6]_.

    Several points to note:

    1. The information about two problem libraries is available in the
       following links:
            S2MPJ (see [3]_) <https://github.com/GrattonToint/S2MPJ>
            MatCUTEst <https://github.com/matcutest>

    2. If you want to use your own problem library, please check the README.txt
       in the directory 'problems/' or the guidance in our website
       <https://www.optprof.com> for more details.

    3. The problem library MatCUTEst is only available when the OS is Linux.

    4. If the 'load' option is provided, we will use the provided options to
       select data from the specified experiment for plotting the profiles.
       Available options are:
       - 'benchmark_id', 'solver_names',
         'feature_stamp', 'errorbar_type', 'savepath', 'max_tol_order',
         'merit_fun', 'run_plain', 'score_only',
         'summarize_performance_profiles', 'summarize_data_profiles',
         'summarize_log_ratio_profiles', 'summarize_output_based_profiles',
         'silent', 'semilogx', 'normalized_scores', 'score_weight_fun',
         'score_fun', 'solvers_to_load', 'line_colors', 'line_styles',
         'line_widths', 'bar_colors'.
       - Options for features: none.
       - 'plibs', 'ptype', 'mindim', 'maxdim', 'minb',
         'maxb', 'minlcon', 'maxlcon', 'minnlcon', 'maxnlcon', 'mincon',
         'maxcon', 'excludelist'.

    5. More information about OptiProfiler can be found on our website:

                            https://www.optprof.com

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

    Examples
    --------
    To do.
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
            logger.info(f'We set {FeatureOption.N_RUNS} to 5 since it is not specified and at least one solver is randomized.')
        feature_options[FeatureOption.N_RUNS] = 5

    # Load the existing results if needed.
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
        n_solvers = results_plibs[0].fun_histories.shape[1]
    elif ProfileOption.SOLVER_NAMES in profile_options:
        n_solvers = len(profile_options[ProfileOption.SOLVER_NAMES])
    else:
        n_solvers = len(solvers)

    solver_scores = np.zeros(n_solvers)
    profile_scores = []
    curves = []
    solver_names = profile_options[ProfileOption.SOLVER_NAMES]

    # Set the options for plotting.
    profile_context = set_profile_context(profile_options)

    # Define the directory to store the results.
    path_out = Path(profile_options[ProfileOption.SAVEPATH], profile_options[ProfileOption.BENCHMARK_ID]).resolve()

    # Record the time stamp of the current experiment. ('yyyyMMdd_HHmmss')
    time_stamp = datetime.now().astimezone().strftime('%Y%m%d_%H%M%S')

    # Set the feature stamp.
    feature_stamp = profile_options[ProfileOption.FEATURE_STAMP]

    # Create the stamp for the current experiment.
    stamp = create_stamp(solver_names, problem_options, feature_stamp, time_stamp)

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
    
    if profile_options[ProfileOption.SCORE_ONLY] or is_load:
        # If 'load' is not empty, we will directly load the results and do not need to compute again. In this case, we do not need to create the directory to store the hist plots.
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
        try:
            with path_readme_log.open('w') as f:
                f.write("# Content of this folder\n\n")
        except Exception:
            if not profile_options[ProfileOption.SILENT]:
                logger.warning(f'Fail to create the README.txt file in {path_log}.')
    else:
        path_readme_log = None

    # Create a text file named by time_stamp to record the time_stamp.
    if not profile_options[ProfileOption.SCORE_ONLY]:
        path_time_stamp = path_log / f'time_stamp_{time_stamp}.txt'
        try:
            with path_time_stamp.open('w') as f:
                f.write(f'{time_stamp}')
            try:
                with path_readme_log.open('a') as f:
                    f.write(f"'time_stamp_{time_stamp}.txt': file, recording the time stamp of the current experiment.\n")
            except:
                pass
        except:
            if not profile_options[ProfileOption.SILENT]:
                logger.warning(f'Fail to create the time stamp file in {path_log}.')
        
    if not profile_options[ProfileOption.SCORE_ONLY] and 'problem' not in locals():
        path_figs = path_log / 'profile_figs'
        path_figs.mkdir(parents=True, exist_ok=True)
        try:
            with path_readme_log.open('a') as f:
                f.write(f"'profile_figs': folder, containing all the FIG files of the profiles.\n")
        except:
            pass

    # Save the options and record the log.
    if not profile_options[ProfileOption.SCORE_ONLY]:
        try:
            if 'options_user' in locals() and options_user is not None:
                # TODO Save the dict of user options.
                pass

                try:
                    # TODO change the file name
                    with path_readme_log.open('a') as f:
                        f.write(f"'options_user': file,storing the options provided by the user for the current experiment.\n")
                except:
                    pass
                options_refined = profile_options.copy()
                feature_options_keys = list(feature_options.keys())
                problem_options_keys = list(problem_options.keys())
                for key in feature_options_keys:
                    options_refined[key] = feature_options[key]
                for key in problem_options_keys:
                    options_refined[key] = problem_options[key]
                # TODO Save the dict of refined options.
                pass
                try:
                    # TODO change the file name
                    with path_readme_log.open('a') as f:
                        f.write(f"'options_refined': file, storing the options refined by OptiProfiler for the current experiment.\n")
                except:
                    pass
        except Exception:
            if not profile_options[ProfileOption.SILENT]:
                logger.warning(f'Failed to save the options of the current experiment.')

        # Set up the logger to log the information in a file.
        log_file = path_log / 'log.txt'
        log_queue, listener = setup_main_process_logging(log_file=log_file, level=logging.INFO)

        try:
            with path_readme_log.open('a') as f:
                f.write(f"'log.txt': file, the log file of the current experiment, recording printed information from the screen.\n")
        except:
            pass

    # Create a README.txt file to explain the content of the folder `path_stamp`.
    if not profile_options[ProfileOption.SCORE_ONLY]:
        path_readme_feature = path_stamp / 'README.txt'
        try:
            with path_readme_feature.open('w') as f:
                f.write("# Content of this folder\n\n")
                if path_hist_plots is not None:
                    f.write(f"'history_plots': folder, containing all the history plots for each problem.\n")
                    f.write(f"'history_plots_summary.pdf': file, the summary PDF of history plots for all problems.\n")
                f.write(f"'test_log': folder, containing log files and other useful experimental data.\n")
        except Exception:
            if not profile_options[ProfileOption.SILENT]:
                logger.warning(f'Fail to create the README.txt file in {path_stamp}.')

    # We try to copy the script or function that calls the benchmark function to the log directory.
    if not profile_options[ProfileOption.SCORE_ONLY]:
        try:
            caller_path = inspect.stack()[1].filename
            if caller_path is not None and Path(caller_path).is_file():
                shutil.copy2(caller_path, path_log / os.path.basename(caller_path))
                try:
                    with path_readme_log.open('a') as f:
                        f.write(f"'{os.path.basename(caller_path)}': file, the script or function that calls the benchmark function.\n")
                except:
                    pass
        except Exception:
            if not profile_options[ProfileOption.SILENT]:
                logger.warning(f'Failed to copy the script or function that calls `benchmark` function to the log directory.')

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
            with path_readme_feature.open('a') as f:
                f.write(f"'detailed_profiles': folder, containing all the high-quality single profiles.\n")
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
        result = _solve_one_problem(solvers, problem, feature, problem.name, len(problem.name), profile_options, True, path_hist_plots)

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
        
        # Compute merit values.
        merit_fun = profile_options[ProfileOption.MERIT_FUN]
        if result is not None:
            try:
                merit_history = compute_merit_values(merit_fun, result['fun_history'], result['maxcv_history'], result['maxcv_init'])
                merit_init = compute_merit_values(merit_fun, result['fun_init'], result['maxcv_init'], result['maxcv_init'])
            except Exception as exc:
                logger.error(f'Error occurred while calculating the merit values. Please check the merit function. Error message: {exc}')
                raise exc
            # Find the least merit value for each problem.
            merit_min = np.nanmin(merit_history)
            merit_min = np.nanmin([merit_min, merit_init])
            # Since we will not compute the profiles, we set `solver_scores` to be the relative decreases in the objective function value.
            solver_merit_mins = np.nanmin(np.nanmin(merit_history, axis=2), axis=1).squeeze()
            solver_scores = (merit_init - solver_merit_mins) / max(merit_init - merit_min, np.finfo(float).eps)
        else:
            solver_scores = np.zeros(n_solvers)

        if not profile_options[ProfileOption.SILENT]:
            logger.info('Scores of the solvers:')
            max_solver_name_length = max(len(name) for name in solver_names)
            for i, name in enumerate(solver_names):
                format_info_str = f'{{:<{max_solver_name_length}}}:    {{:.4f}}'
                logger.info(format_info_str.format(name, solver_scores[i]))
    
        return solver_scores, None, None

    # If 'load' option is not specified, we solve all the selected problems.
    if not is_load:
        # Print the information of the current experiment.
        if not profile_options[ProfileOption.SILENT]:
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
                logger.info(f'Start testing problems from the problem library {plib} with {feature.name} feature.')
                if plib == 's2mpj':
                    logger.info('More information about the S2MPJ problem library can be found at: https://github.com/GrattonToint/S2MPJ')
                elif plib == 'pycutest':
                    logger.info('More information about the PyCUTEst problem library can be found at: https://jfowkes.github.io/pycutest/_build/html/index.html')

            # Create directory to store the history plots for each problem library.
            path_hist_plots_plib = path_hist_plots / plib if path_hist_plots is not None else None
            if path_hist_plots_plib is not None and not path_hist_plots_plib.exists():
                path_hist_plots_plib.mkdir(parents=True, exist_ok=True)

            # Solve all the problems from the current problem library with the specified options and get the computation results.
            results_plib = _solve_all_problems(solvers, plib, feature, problem_options, profile_options, True, path_hist_plots_plib, log_queue=log_queue)

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
                logger.error(f'Error occurred while calculating the merit values. Please check the merit function. Error message: {exc}')
                raise exc
            results_plib['merit_histories'] = merit_histories
            results_plib['merit_outs'] = merit_outs
            results_plib['merit_inits'] = merit_inits

            # Run the 'plain' feature if run_plain is true.
            if profile_options[ProfileOption.RUN_PLAIN]:
                feature_plain = Feature(FeatureName.PLAIN.value)
                if not profile_options[ProfileOption.SILENT]:
                    logger.info(f'Start testing problems from the problem library {plib} with "plain" feature.')
                results_plib_plain = _solve_all_problems(solvers, plib, feature_plain, problem_options, profile_options, False, None, log_queue=log_queue)
                try:
                    results_plib_plain['merit_histories'] = compute_merit_values(merit_fun, results_plib_plain['fun_histories'], results_plib_plain['maxcv_histories'], results_plib_plain['maxcv_inits'])
                    results_plib_plain['merit_outs'] = compute_merit_values(merit_fun, results_plib_plain['fun_outs'], results_plib_plain['maxcv_outs'], results_plib_plain['maxcv_inits'])
                    results_plib_plain['merit_inits'] = compute_merit_values(merit_fun, results_plib_plain['fun_inits'], results_plib_plain['maxcv_inits'], results_plib_plain['maxcv_inits'])
                except Exception as exc:
                    logger.error(f'Error occurred while calculating the merit values for the "plain" feature. Please check the merit function. Error message: {exc}')
                    raise exc
                
                # Store data of the 'plain' feature for later calculating merit_mins.
                results_plib['results_plib_plain'] = results_plib_plain

            # Append the results of the current problem library to the list.
            results_plibs.append(results_plib)

            # Merge the history plots for each problem library to a single pdf file.
            if not profile_options[ProfileOption.SCORE_ONLY] and np.any(results_plib['solvers_successes']):
                if not profile_options[ProfileOption.SILENT]:
                    logger.info(f'Merging all the history plots of problems from the "{plib}" library to a single PDF file.')
                try:
                    merge_pdfs_with_pypdf(path_hist_plots_plib, path_hist_plots / f'{plib}_history_plots_summary.pdf')
                except Exception as exc:
                    if not profile_options[ProfileOption.SILENT]:
                        logger.warning(f'Failed to merge the history plots to a single PDF file. Error message: {exc}')

        # Remove the None elements from results_plibs.
        results_plibs = [results_plib for results_plib in results_plibs if results_plib is not None]
        if len(results_plibs) == 0:
            if not profile_options[ProfileOption.SILENT]:
                logger.warning('No problems are selected or solved from any problem library.')
            return solver_scores, None, None

    # Store the data for loading.
    # TODO
    pass



    # Write the report file.
    write_report(profile_options, results_plibs, path_report, path_readme_log)

    # Process the results from all the problem libraries.
    merit_histories_merged, merit_outs_merged, merit_inits_merged, merit_mins_merged, n_evals_merged, problem_names_merged, problem_dims_merged = process_results(results_plibs, profile_options)
    n_problems, n_solvers, n_runs = merit_histories_merged.shape[:3]

    if not profile_options[ProfileOption.SILENT]:
        logger.info('Start computing profiles.')

    max_tol_order = profile_options[ProfileOption.MAX_TOL_ORDER]
    tolerances = [10**(-i) for i in range(1, max_tol_order + 1)]

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
    multiplier = 2 if is_output_based else 1
    default_figsize = plt.rcParams['figure.figsize']
    default_width = default_figsize[0]
    default_height = default_figsize[1]




    with plt.rc_context(profile_context):
        fig_summary = plt.figure(figsize=(len(tolerances) * default_width, multiplier * n_rows * default_height), layout='constrained')
        if multiplier == 2:
            fig_summary_hist, fig_summary_out = fig_summary.subfigures(2, 1)
            subfigs_summary_hist = fig_summary_hist.subfigures(n_rows, 1)
            subfigs_summary_out = fig_summary_out.subfigures(n_rows, 1)
        else:
            fig_summary_hist = fig_summary.subfigures(1, 1)
            fig_summary_out = None
            subfigs_summary_hist = fig_summary_hist.subfigures(n_rows, 1)
            subfigs_summary_out = None
        
        i_rows = 0
        if is_perf:
            ax_summary_perf_hist = subfigs_summary_hist[i_rows].subplots(1, len(tolerances), sharey=True)
            ax_summary_perf_out = subfigs_summary_out[i_rows].subplots(1, len(tolerances), sharey=True) if multiplier == 2 else None
            i_rows += 1
        else:
            ax_summary_perf_hist = None
            ax_summary_perf_out = None
        if is_data:
            ax_summary_data_hist = subfigs_summary_hist[i_rows].subplots(1, len(tolerances), sharey=True)
            ax_summary_data_out = subfigs_summary_out[i_rows].subplots(1, len(tolerances), sharey=True) if multiplier == 2 else None
            i_rows += 1
        else:
            ax_summary_data_hist = None
            ax_summary_data_out = None
        if is_log_ratio:
            ax_summary_log_ratio_hist = subfigs_summary_hist[i_rows].subplots(1, len(tolerances))
            ax_summary_log_ratio_out = subfigs_summary_out[i_rows].subplots(1, len(tolerances)) if multiplier == 2 else None
            i_rows += 1
        else:
            ax_summary_log_ratio_hist = None
            ax_summary_log_ratio_out = None

        # Find the problems that all the solvers failed to meet the convergence test for every tolerance.
        solvers_all_diverge_hist = np.zeros((n_problems, n_runs, profile_options[ProfileOption.MAX_TOL_ORDER]), dtype=bool)
        solvers_all_diverge_out = solvers_all_diverge_hist.copy()

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
            curve = {
                'hist': hist,
                'out': hist.copy()
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
                        if np.isfinite(merit_mins_merged[i_problem]):
                            threshold = max(tolerance * merit_inits_merged[i_problem] + (1.0 - tolerance) * merit_mins_merged[i_problem], merit_mins_merged[i_problem])
                        else:
                            threshold = -np.inf
                        if np.min(merit_histories_merged[i_problem, i_solver, i_run, :]) <= threshold:
                            work_hist[i_problem, i_solver, i_run] = np.argmax(merit_histories_merged[i_problem, i_solver, i_run, :] <= threshold) + 1
                        if merit_outs_merged[i_problem, i_solver, i_run] <= threshold:
                            work_out[i_problem, i_solver, i_run] = n_evals_merged[i_problem, i_solver, i_run]
                        
            for i_problem in range(n_problems):
                for i_run in range(n_runs):
                    solvers_all_diverge_hist[i_problem, i_run, i_tol] = np.all(np.isnan(work_hist[i_problem, :, i_run]))
                    solvers_all_diverge_out[i_problem, i_run, i_tol] = np.all(np.isnan(work_out[i_problem, :, i_run]))

            # Draw the profiles.
            fig_perf_hist, fig_data_hist, fig_log_ratio_hist, curve['hist'] = draw_profiles(work_hist, problem_dims_merged, solver_names, tolerance_latex, i_tol, ax_summary_perf_hist, ax_summary_data_hist, ax_summary_log_ratio_hist, True, is_perf, is_data, is_log_ratio, profile_options, curve['hist'])
            fig_perf_out, fig_data_out, fig_log_ratio_out, curve['out'] = draw_profiles(work_out, problem_dims_merged, solver_names, tolerance_latex, i_tol, ax_summary_perf_out, ax_summary_data_out, ax_summary_log_ratio_out, is_output_based, is_perf, is_data, is_log_ratio, profile_options, curve['out'])
            curves.append(curve)

            # Save the profiles to files.
            if not profile_options[ProfileOption.SCORE_ONLY]:
                # Save the individual profiles to separate PDF files and summary files.
                pdf_perf_hist = path_perf_hist / f'perf_hist_{i_tol + 1}.pdf'
                fig_perf_hist.savefig(pdf_perf_hist, bbox_inches='tight')
                pdf_perf_hist_summary.savefig(fig_perf_hist, bbox_inches='tight')

                pdf_data_hist = path_data_hist / f'data_hist_{i_tol + 1}.pdf'
                fig_data_hist.savefig(pdf_data_hist, bbox_inches='tight')
                pdf_data_hist_summary.savefig(fig_data_hist, bbox_inches='tight')

                pdf_perf_out = path_perf_out / f'perf_out_{i_tol + 1}.pdf'
                fig_perf_out.savefig(pdf_perf_out, bbox_inches='tight')
                pdf_perf_out_summary.savefig(fig_perf_out, bbox_inches='tight')

                pdf_data_out = path_data_out / f'data_out_{i_tol + 1}.pdf'
                fig_data_out.savefig(pdf_data_out, bbox_inches='tight')
                pdf_data_out_summary.savefig(fig_data_out, bbox_inches='tight')

                if n_solvers == 2:
                    pdf_log_ratio_hist = path_log_ratio_hist / f'log-ratio_hist_{i_tol + 1}.pdf'
                    fig_log_ratio_hist.savefig(pdf_log_ratio_hist, bbox_inches='tight')
                    pdf_log_ratio_hist_summary.savefig(fig_log_ratio_hist, bbox_inches='tight')
                    
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
        pdf_perf_hist_summary.close()
        pdf_perf_out_summary.close()
        pdf_data_hist_summary.close()
        pdf_data_out_summary.close()
        if n_solvers == 2:
            pdf_log_ratio_hist_summary.close()
            pdf_log_ratio_out_summary.close()
        with open(path_readme_feature, 'a') as f:
            f.write("'data_hist.pdf': file, the summary PDF of history-based data profiles for all tolerances.\n")
            f.write("'data_out.pdf': file, the summary PDF of output-based data profiles for all tolerances.\n")
            f.write("'perf_hist.pdf': file, the summary PDF of history-based performance profiles for all tolerances.\n")
            f.write("'perf_out.pdf': file, the summary PDF of output-based performance profiles for all tolerances.\n")
            if n_solvers == 2:
                f.write("'log-ratio_hist.pdf': file, the summary PDF of history-based log-ratio profiles for all tolerances.\n")
                f.write("'log-ratio_out.pdf': file, the summary PDF of output-based log-ratio profiles for all tolerances.\n")
        
        # TODO Record the names of the problems all the solvers failed to meet the convergence test for every tolerance.
        pass

        # Save the summary for the current feature.
        fig_summary_hist.supylabel('History-based profiles', fontsize='xx-large', horizontalalignment='right')
        fig_summary_out.supylabel('Output-based profiles', fontsize='xx-large', horizontalalignment='right')
        fig_summary.suptitle(f"Profiles with the ``{feature.name}'' feature", fontsize='xx-large', verticalalignment='bottom')
        path_summary = path_stamp / f'summary_{stamp}.pdf'
        fig_summary.savefig(path_summary, bbox_inches='tight')

        plt.close(fig_summary)

    

    


    # Close the listener of the logger.
    if not profile_options[ProfileOption.SCORE_ONLY]:
        listener.stop()
    
    """

    # Run the benchmarks.
    pdf_summary = backend_pdf.PdfPages(path_summary)
    problem_names = None
    merit_init = None
    merit_histories_plain = None
    merit_out_plain = None
    for feature in features:
        # Solve the problems.
        logger.info(f'Starting the computations of the "{feature.name}" profiles.')
        max_eval_factor = 500
        if feature.name == FeatureName.PLAIN and merit_histories_plain is not None:
            merit_histories = np.copy(merit_histories_plain)
            merit_out = np.copy(merit_out_plain)
        else:
            problem_names, fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_dimensions, time_processes = _solve_all_problems(cutest_problem_names, custom_problem_loader, custom_problem_names, solvers, labels, feature, max_eval_factor, profile_options)
            merit_histories = _compute_merit_values(fun_histories, maxcv_histories, maxcv_init)
            merit_out = _compute_merit_values(fun_out, maxcv_out, maxcv_init)
            merit_init = _compute_merit_values(fun_init, maxcv_init, maxcv_init)
            if feature.name == FeatureName.PLAIN:
                merit_histories_plain = np.copy(merit_histories)
                merit_out_plain = np.copy(merit_out)

        # Determine the least merit value for each problem.
        merit_min = np.min(merit_histories, (1, 2, 3))
        if feature.is_stochastic:
            if merit_histories_plain is None:
                feature_plain = Feature('plain')
                logger.info(f'Starting the computations of the "plain" profiles.')
                problem_names, fun_histories_plain, maxcv_histories_plain, fun_out_plain, maxcv_out_plain, fun_init, maxcv_init, _, _, time_processes_plain = _solve_all_problems(cutest_problem_names, custom_problem_loader, custom_problem_names, solvers, labels, feature_plain, max_eval_factor, profile_options)
                merit_histories_plain = _compute_merit_values(fun_histories_plain, maxcv_histories_plain, maxcv_init)
                merit_out_plain = _compute_merit_values(fun_out_plain, maxcv_out_plain, maxcv_init)
                merit_init = _compute_merit_values(fun_init, maxcv_init, maxcv_init)
            merit_min_plain = np.min(merit_histories_plain, (1, 2, 3))
            merit_min = np.minimum(merit_min, merit_min_plain)

        # Paths to the individual results.
        path_feature = path_out / feature.name
        path_feature.mkdir(parents=True, exist_ok=True)
        path_problems = path_feature / 'report.txt'
        path_perf_hist = path_feature / 'perf_hist.pdf'
        path_perf_out = path_feature / 'perf_out.pdf'
        path_data_hist = path_feature / 'data_hist.pdf'
        path_data_out = path_feature / 'data_out.pdf'
        path_log_ratio_hist = path_feature / 'log-ratio_hist.pdf'
        path_log_ratio_out = path_feature / 'log-ratio_out.pdf'

        # Store the names of the problems.
        with path_problems.open('w') as f:
            f.write(os.linesep.join(problem_names))

        with plt.rc_context(profile_context):
            # Create the summary figure.
            n_rows = 0
            if profile_options[ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES]:
                n_rows += 1
            if profile_options[ProfileOption.SUMMARIZE_DATA_PROFILES]:
                n_rows += 1
            if profile_options[ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES]:
                n_rows += 1
            fig_summary = plt.figure(figsize=(len(tolerances) * 4.8, 2 * n_rows * 4.8), layout='constrained')
            fig_summary_hist, fig_summary_out = fig_summary.subfigures(2, 1)
            subfigs_summary_hist = fig_summary_hist.subfigures(n_rows, 1)
            subfigs_summary_out = fig_summary_out.subfigures(n_rows, 1)
            i_rows = 0
            if profile_options[ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES]:
                ax_summary_perf_hist = subfigs_summary_hist[i_rows].subplots(1, len(tolerances), sharey=True)
                ax_summary_perf_out = subfigs_summary_out[i_rows].subplots(1, len(tolerances), sharey=True)
                i_rows += 1
            else:
                ax_summary_perf_hist = None
                ax_summary_perf_out = None
            if profile_options[ProfileOption.SUMMARIZE_DATA_PROFILES]:
                ax_summary_data_hist = subfigs_summary_hist[i_rows].subplots(1, len(tolerances), sharey=True)
                ax_summary_data_out = subfigs_summary_out[i_rows].subplots(1, len(tolerances), sharey=True)
                i_rows += 1
            else:
                ax_summary_data_hist = None
                ax_summary_data_out = None
            if profile_options[ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES]:
                ax_summary_log_ratio_hist = subfigs_summary_hist[i_rows].subplots(1, len(tolerances))
                ax_summary_log_ratio_out = subfigs_summary_out[i_rows].subplots(1, len(tolerances))
                i_rows += 1
            else:
                ax_summary_log_ratio_hist = None
                ax_summary_log_ratio_out = None

            # Create the performance and data profiles.
            n_problems, n_solvers, n_runs, max_eval = merit_histories.shape
            pdf_perf_hist = backend_pdf.PdfPages(path_perf_hist)
            pdf_perf_out = backend_pdf.PdfPages(path_perf_out)
            pdf_data_hist = backend_pdf.PdfPages(path_data_hist)
            pdf_data_out = backend_pdf.PdfPages(path_data_out)
            pdf_log_ratio_hist = backend_pdf.PdfPages(path_log_ratio_hist, False)
            pdf_log_ratio_out = backend_pdf.PdfPages(path_log_ratio_out, False)
            for i_tolerance, tolerance in enumerate(tolerances):
                tolerance_str, tolerance_latex = _format_float_scientific_latex(tolerance)
                logger.info(f'Creating profiles for tolerance {tolerance_str}.')
                tolerance_label = f'($\\mathrm{{tol}} = {tolerance_latex}$)'

                # Compute the number of function evaluations used by each
                # solver on each problem at each run to achieve convergence.
                work_hist = np.full((n_problems, n_solvers, n_runs), np.nan)
                work_out = np.full((n_problems, n_solvers, n_runs), np.nan)
                for i_problem in range(n_problems):
                    for i_solver in range(n_solvers):
                        for i_run in range(n_runs):
                            if np.isfinite(merit_min[i_problem]):
                                threshold = max(tolerance * merit_init[i_problem] + (1.0 - tolerance) * merit_min[i_problem], merit_min[i_problem])
                            else:
                                threshold = -np.inf
                            if np.min(merit_histories[i_problem, i_solver, i_run, :]) <= threshold:
                                work_hist[i_problem, i_solver, i_run] = np.argmax(merit_histories[i_problem, i_solver, i_run, :] <= threshold) + 1
                            if merit_out[i_problem, i_solver, i_run] <= threshold:
                                work_out[i_problem, i_solver, i_run] = n_eval[i_problem, i_solver, i_run]

                # Draw and save the profiles.
                fig_perf_hist, fig_perf_out, fig_data_hist, fig_data_out, fig_log_ratio_hist, fig_log_ratio_out = _draw_profiles(work_hist, work_out, problem_dimensions, labels, tolerance_label, i_tolerance, ax_summary_perf_hist, ax_summary_data_hist, ax_summary_log_ratio_hist, ax_summary_perf_out, ax_summary_data_out, ax_summary_log_ratio_out)
                pdf_perf_hist.savefig(fig_perf_hist, bbox_inches='tight')
                pdf_perf_out.savefig(fig_perf_out, bbox_inches='tight')
                pdf_data_hist.savefig(fig_data_hist, bbox_inches='tight')
                pdf_data_out.savefig(fig_data_out, bbox_inches='tight')
                if fig_log_ratio_hist is not None:
                    pdf_log_ratio_hist.savefig(fig_log_ratio_hist, bbox_inches='tight')
                if fig_log_ratio_out is not None:
                    pdf_log_ratio_out.savefig(fig_log_ratio_out, bbox_inches='tight')

                # Close the individual figures.
                plt.close(fig_perf_hist)
                plt.close(fig_perf_out)
                plt.close(fig_data_hist)
                plt.close(fig_data_out)
                if fig_log_ratio_hist is not None:
                    plt.close(fig_log_ratio_hist)
                if fig_log_ratio_out is not None:
                    plt.close(fig_log_ratio_out)

            # Close the individual PDF files.
            pdf_perf_hist.close()
            pdf_perf_out.close()
            pdf_data_hist.close()
            pdf_data_out.close()
            pdf_log_ratio_hist.close()
            pdf_log_ratio_out.close()
            logger.info(f'Detailed results stored in {path_feature}.')

            # Save the summary for the current feature.
            fig_summary_hist.supylabel('History-based profiles', fontsize='xx-large', horizontalalignment='right')
            fig_summary_out.supylabel('Output-based profiles', fontsize='xx-large', horizontalalignment='right')
            fig_summary.suptitle(f"Profiles with the ``{feature.name}'' feature", fontsize='xx-large', verticalalignment='bottom')
            pdf_summary.savefig(fig_summary, bbox_inches='tight')

            # Close the summary figure.
            plt.close(fig_summary)

    # Close the summary PDF file.
    pdf_summary.close()
    logger.info(f'Summary stored in {path_summary}.')
    
    """


def _solve_all_problems(solvers, plib, feature, problem_options, profile_options, is_plot, path_hist_plots, log_queue=None):
    """
    Solve all problems in plib satisfying problem_options using solvers in the solvers and stores the computing results.
    """
    logger = get_logger(__name__)
    results = []

    # Get satisfied problem names.
    option_select = problem_options.copy()
    for key in [ProblemOption.PLIBS.value, ProblemOption.PROBLEM_NAMES.value, ProblemOption.EXCLUDELIST.value]:
        option_select.pop(key, None)

    current_dir = Path(__file__).parent.resolve()
    module_file_path = current_dir / '../..' / 'problems' / plib / f'{plib}_tools.py'
    module_file_path = module_file_path.resolve()
    module_name = 'problems.' + plib + '.' + plib + '_tools'

    # Import the problem library module.
    try:
        spec = importlib.util.spec_from_file_location(module_name, str(module_file_path))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        selector_name = plib + '_select'
        select = getattr(module, selector_name)
    except Exception as exc:
        logger.error(f'Error occurred while importing the problem library {plib}. Error message: {exc}')
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
            logger.info(f'No problem is selected from "{plib}".')
        return results

    n_problems = len(problem_names)
    len_problem_names = max(len(name) for name in problem_names)
    max_eval_factor = profile_options[ProfileOption.MAX_EVAL_FACTOR]
    if not profile_options[ProfileOption.SILENT]:
        logger.info(f'There are {n_problems} problems selected from "{plib}" to test.')

    # Determine whether to use sequential mode or parallel mode.
    seqential_mode = profile_options[ProfileOption.N_JOBS] == 1 or os.cpu_count() == 1

    # Solve all problems.
    args = [(solvers, feature, problem_name, len_problem_names, profile_options, is_plot, path_hist_plots, plib) for problem_name in problem_names]
    if seqential_mode:
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

    return results


def _solve_one_problem_wrapper(solvers, feature, problem_name, len_problem_names, profile_options, is_plot, path_hist_plots, plib):
    logger = get_logger(__name__)
    current_dir = Path(__file__).parent.resolve()
    module_file_path = current_dir / '../..' / 'problems' / plib / f'{plib}_tools.py'
    module_file_path = module_file_path.resolve()
    module_name = 'problems.' + plib + '.' + plib + '_tools'
    # Import the problem library module.
    try:
        spec = importlib.util.spec_from_file_location(module_name, str(module_file_path))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        loader_name = plib + '_load'
        load = getattr(module, loader_name)
    except Exception as exc:
        logger.error(f'Error occurred while importing the problem library {plib}. Error message: {exc}')
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


def _solve_one_problem(solvers, problem, feature, problem_name, len_problem_names, profile_options, is_plot, path_hist_plots):
    """
    Solve a given problem.
    """
    solver_names = profile_options[ProfileOption.SOLVER_NAMES]
    solver_isrand = profile_options[ProfileOption.SOLVER_ISRAND]
    result = None

    if type(problem) is not Problem:
        raise TypeError('The problem provided is not a valid Problem object.')

    # Project the initial point if necessary.
    if profile_options[ProfileOption.PROJECT_X0]:
        problem.project_x0()

    # Evaluate the functions at the initial point.
    fun_init = problem.fun(problem.x0)
    maxcv_init = problem.maxcv(problem.x0)

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
    computation_time = np.full((n_solvers, n_runs), np.nan)
    solvers_success = np.full((n_solvers, n_runs), False)

    # The number of real runs for each solver, which is determined by feature and solver_isrand.
    real_n_runs = np.array([
        n_runs if feature.is_stochastic or (solver_isrand is not None and solver_isrand[i_solver]) else 1
        for i_solver in range(n_solvers)
    ], dtype=int)
    
    len_solver_names = max(len(name) for name in solver_names)

    logger = get_logger(__name__)

    # Solve the problem with each solver.
    for i_solver in range(n_solvers):
        for i_run in range(real_n_runs[i_solver]):
            if not profile_options[ProfileOption.SILENT]:
                logger.info(
                    f"Start  solving    {problem_name:<{len_problem_names}} with {solver_names[i_solver]:<{len_solver_names}} (run {i_run + 1:2d}/{real_n_runs[i_solver]:2d})."
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
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                time_start_solver_run = time.monotonic()
                try:
                    if profile_options[ProfileOption.SOLVER_VERBOSE] == 2:
                        x = call_solver()
                    elif profile_options[ProfileOption.SOLVER_VERBOSE] == 1:
                        with open(os.devnull, 'w') as devnull, redirect_stdout(devnull):
                            x = call_solver()
                    else:
                        with open(os.devnull, 'w') as devnull, redirect_stdout(devnull), redirect_stderr(devnull):
                            x = call_solver()

                    computation_time[i_solver, i_run] = time.monotonic() - time_start_solver_run

                    # It is very important to transform the solution back to the one related to the original problem. (Note that the problem we solve has the objective function f(A @ x + b). Thus, if x is the output solution, then A @ x + b is the solution of the original problem.)
                    A, b = featured_problem._feature.modifier_affine(featured_problem._seed, featured_problem._problem)[:2]
                    x = A @ x + b

                    # Use problem.fun and problem.maxcv to evaluate the solution since it is possible that featured_problem.fun and featured_problem.maxcv are modified.
                    fun_out[i_solver, i_run] = problem.fun(x)
                    maxcv_out[i_solver, i_run] = problem.maxcv(x)
                    # Calculate the minimum function value and the minimum constraint violation, omitting the NaN values.
                    fun_min = np.nanmin(featured_problem.fun_hist)
                    maxcv_min = np.nanmin(featured_problem.maxcv_hist)

                    if not profile_options[ProfileOption.SILENT]:
                        logger.info(
                            f"Finish solving    {problem_name:<{len_problem_names}} with {solver_names[i_solver]:<{len_solver_names}} "
                            f"(run {i_run + 1:2d}/{real_n_runs[i_solver]:2d}) (in {computation_time[i_solver, i_run]:.2f} seconds)."
                        )

                        if problem_type == 'u':
                            logger.info(
                                f"Output result for {problem_name:<{len_problem_names}} with {solver_names[i_solver]:<{len_solver_names}} "
                                f"(run {i_run + 1:2d}/{real_n_runs[i_solver]:2d}): f = {fun_out[i_solver, i_run]:10.4e}."
                            )
                            logger.info(
                                f"Best   result for {problem_name:<{len_problem_names}} with {solver_names[i_solver]:<{len_solver_names}} "
                                f"(run {i_run + 1:2d}/{real_n_runs[i_solver]:2d}): f = {fun_min:10.4e}."
                            )
                        else:
                            logger.info(
                                f"Output result for {problem_name:<{len_problem_names}} with {solver_names[i_solver]:<{len_solver_names}} "
                                f"(run {i_run + 1:2d}/{real_n_runs[i_solver]:2d}): f = {fun_out[i_solver, i_run]:10.4e}, maxcv = {maxcv_out[i_solver, i_run]:10.4e}."
                            )
                            logger.info(
                                f"Best   result for {problem_name:<{len_problem_names}} with {solver_names[i_solver]:<{len_solver_names}} "
                                f"(run {i_run + 1:2d}/{real_n_runs[i_solver]:2d}): f = {fun_min:10.4e}, maxcv = {maxcv_min:10.4e}."
                            )
                except Exception as exc:
                    if profile_options[ProfileOption.SOLVER_VERBOSE] >= 1:
                        logger.warning(f'An error occurred while solving {problem_name} with {solver_names[i_solver]}: {exc}')
            n_eval[i_solver, i_run] = featured_problem.n_eval
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

    # Return the result.
    result = {
        'fun_history': fun_history,
        'maxcv_history': maxcv_history,
        'fun_out': fun_out,
        'maxcv_out': maxcv_out,
        'fun_init': fun_init,
        'maxcv_init': maxcv_init,
        'n_eval': n_eval,
        'problem_name': problem_name,
        'problem_type': problem_type,
        'problem_dim': problem_dim,
        'problem_mb': problem_mb,
        'problem_mlcon': problem_mlcon,
        'problem_mnlcon': problem_mnlcon,
        'problem_mcon': problem_mcon,
        'computation_time': computation_time,
        'solvers_success': solvers_success
    }

    # Draw the history plots if required.
    if not is_plot or profile_options[ProfileOption.SCORE_ONLY] or all(not success for success in solvers_success.flatten()):
        return result
    
    try:
        merit_fun = profile_options[ProfileOption.MERIT_FUN]
        try:
            merit_history = compute_merit_values(merit_fun, fun_history, maxcv_history, maxcv_init)
            merit_init = compute_merit_values(merit_fun, fun_init, maxcv_init, maxcv_init)
        except Exception as exc:
            logger.error(f'Error occurred while calculating the merit values. Please check the merit function. Error message: {exc}')
            raise exc

        warnings.filterwarnings('ignore')
        n_cols = 1 if problem_type == 'u' else 3
        default_figsize = plt.rcParams['figure.figsize']
        default_width = default_figsize[0]
        default_height = default_figsize[1]
        profile_context = set_profile_context(profile_options)
        with plt.rc_context(profile_context):
            fig_summary = plt.figure(figsize=(default_width * n_cols, default_height * 2), layout='constrained')

            F_title = profile_options['feature_stamp'].replace('_', r'\_')
            P_title = problem_name.replace('_', r'\_')
            T_title = f'Solving ``{P_title}" with ``{F_title}" feature'

            fig_summary.suptitle(T_title, verticalalignment='bottom')

            # Create axes arrays
            axs_summary = []
            
            # Create subplots - first row
            for j in range(n_cols):
                ax = fig_summary.add_subplot(2, n_cols, j+1)
                axs_summary.append(ax)
            
            # Create subplots - second row
            for j in range(n_cols):
                ax = fig_summary.add_subplot(2, n_cols, n_cols+j+1)
                axs_summary.append(ax)
            
            # Add y-axis labels to the rows
            fig_summary.text(0.04, 0.75, "History profiles", rotation=90, va='center')
            fig_summary.text(0.04, 0.25, "Cummin history profiles", rotation=90, va='center')
            
            # Create cell_axs_summary and cell_axs_summary_cum
            if problem_type == 'u':
                cell_axs_summary = [axs_summary[0]]
                cell_axs_summary_cum = [axs_summary[1]]
            else:
                cell_axs_summary = [axs_summary[0], axs_summary[1], axs_summary[2]]
                cell_axs_summary_cum = [axs_summary[3], axs_summary[4], axs_summary[5]]
            
            # Create PDF filename
            pdf_hist_file_name = re.sub(r'^[-_]+', '', re.sub(r'[-_]+', '_', re.sub(r'[^a-zA-Z0-9\-_]', '', problem_name.replace(' ', '_')))) + '.pdf'
            pdf_summary = os.path.join(path_hist_plots, pdf_hist_file_name)
            
            # Process solver names (replace underscores)
            processed_solver_names = [name.replace('_', r'\_') for name in solver_names]
            
            # Draw history plots
            draw_hist(fun_history, maxcv_history, merit_history, fun_init, maxcv_init, merit_init, processed_solver_names, cell_axs_summary, False, problem_type, problem_dim, n_eval, profile_options, default_height)
            
            # Draw cumulative minimum history plots
            draw_hist(fun_history, maxcv_history, merit_history, fun_init, maxcv_init, merit_init, processed_solver_names, cell_axs_summary_cum, True, problem_type, problem_dim, n_eval, profile_options, default_height)
            
            # Make layout tight to avoid overlapping
            fig_summary.tight_layout(rect=[0.05, 0.0, 1.0, 0.95])  # Leave space for the row labels and title
            
            # Save figure
            fig_summary.savefig(pdf_summary, bbox_inches='tight')
            plt.close(fig_summary)

    except Exception as exc:
        if not profile_options[ProfileOption.SILENT]:
            logger.info(f'An error occurred while plotting the history plots of the problem {problem_name}: {exc}')
        pass
    
    return result











    

            
        
