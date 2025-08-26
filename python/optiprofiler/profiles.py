import os
import re
import shutil
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
from matplotlib.backends import backend_pdf
from matplotlib.ticker import MaxNLocator, FuncFormatter

from .features import Feature
from .problems import FeaturedProblem
from .utils import FeatureName, ProfileOption, FeatureOption, ProblemOption, ProblemError, get_logger


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

    # Set the default profile options.
    profile_options.setdefault(ProfileOption.N_JOBS.value, None)
    profile_options.setdefault(ProfileOption.BENCHMARK_ID.value, '.')
    profile_options.setdefault(ProfileOption.PROJECT_X0.value, False)
    profile_options.setdefault(ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES.value, True)
    profile_options.setdefault(ProfileOption.SUMMARIZE_DATA_PROFILES.value, True)
    profile_options.setdefault(ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES.value, False)

    # Check whether the profile options are valid.
    if isinstance(profile_options[ProfileOption.N_JOBS], float) and profile_options[ProfileOption.N_JOBS].is_integer():
        profile_options[ProfileOption.N_JOBS] = int(profile_options[ProfileOption.N_JOBS])
    if profile_options[ProfileOption.N_JOBS] is not None and not (isinstance(profile_options[ProfileOption.N_JOBS], int) and profile_options[ProfileOption.N_JOBS] > 0):
        raise TypeError(f'Option {ProfileOption.N_JOBS} must be a positive integer or None.')
    if not isinstance(profile_options[ProfileOption.BENCHMARK_ID], str):
        raise TypeError(f'Option {ProfileOption.BENCHMARK_ID} must be a string.')
    if not isinstance(profile_options[ProfileOption.PROJECT_X0], bool):
        raise TypeError(f'Option {ProfileOption.PROJECT_X0} must be a boolean.')
    if not isinstance(profile_options[ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES], bool):
        raise TypeError(f'Option {ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES} must be a boolean.')
    if not isinstance(profile_options[ProfileOption.SUMMARIZE_DATA_PROFILES], bool):
        raise TypeError(f'Option {ProfileOption.SUMMARIZE_DATA_PROFILES} must be a boolean.')
    if not isinstance(profile_options[ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES], bool):
        raise TypeError(f'Option {ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES} must be a boolean.')

    # Build the features.
    if isinstance(feature_name, str):
        if feature_name.lower() == 'all':
            feature_name = [name.value for name in FeatureName.__members__.values()]
            feature_name.remove(FeatureName.CUSTOM.value)
        else:
            feature_name = [feature_name]
    feature_name = [name.lower() for name in feature_name]
    features = [Feature(name, **feature_options) for name in feature_name]

    # Path to the summary.
    timestamp = datetime.now().astimezone().strftime('%Y-%m-%dT%H-%M-%S%Z')
    path_out = Path('out', profile_options[ProfileOption.BENCHMARK_ID], timestamp).resolve()
    path_summary = path_out / 'summary.pdf'

    # Retrieve the tolerances of the profiles.
    tolerances = np.logspace(-1, -16, 16)

    # Customize the matplotlib style.
    prop_cycle = cycler(color=[
        '#1f77b4',
        '#ff7f0e',
        '#2ca02c',
        '#d62728',
        '#9467bd',
        '#8c564b',
        '#e377c2',
        '#7f7f7f',
        '#bcbd22',
        '#17becf',
    ])
    prop_cycle += cycler(linestyle=[
        (0, ()),
        (0, (1, 1.5)),
        (0, (3, 1.5)),
        (0, (5, 1.5, 1, 1.5)),
        (0, (5, 1.5, 1, 1.5, 1, 1.5)),
        (0, (1, 3)),
        (0, (3, 3)),
        (0, (5, 3, 1, 3)),
        (0, (5, 3, 1, 3, 1, 3)),
        (0, (1, 4.5)),
    ])
    profile_context = {
        'axes.prop_cycle': prop_cycle,
        'font.family': 'serif',
        'font.size': 16,
        'text.usetex': True if shutil.which('latex') else False,
    }

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


def _solve_all_problems(cutest_problem_names, custom_problem_loader, custom_problem_names, solvers, labels, feature, max_eval_factor, profile_options):
    """
    Solve all problems in parallel.
    """
    problem_names = cutest_problem_names + [(f'EXTRA{i_problem}', custom_problem_name) for i_problem, custom_problem_name in enumerate(custom_problem_names)]
    cutest_problem_options = get_cutest_problem_options()

    # Solve all problems.
    logger = get_logger(__name__)
    logger.info('Entering the parallel section.')
    args = [(problem_name, solvers, labels, feature, max_eval_factor, custom_problem_loader, cutest_problem_options, profile_options) for problem_name in problem_names]
    if profile_options[ProfileOption.N_JOBS] == 1 or os.cpu_count() == 1:
        results = map(lambda arg: _solve_one_problem(*arg), args)
    else:
        with Pool(profile_options[ProfileOption.N_JOBS]) as p:
            results = p.starmap(_solve_one_problem, args)
    logger.info('Leaving the parallel section.')
    if all(result is None for result in results):
        logger.critical('All problems failed to load.')
        problem_names = []
        _fun_histories, _maxcv_histories = [], []
        fun_out, maxcv_out = [], []
        fun_init, maxcv_init = [], []
        n_eval = []
        problem_dimensions = []
        time_processes = []
    else:
        problem_names, _fun_histories, _maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_dimensions, time_processes = zip(*[result for result in results if result is not None])
    fun_out = np.array(fun_out)
    maxcv_out = np.array(maxcv_out)
    fun_init = np.array(fun_init)
    maxcv_init = np.array(maxcv_init)
    n_eval = np.array(n_eval)
    problem_dimensions = np.array(problem_dimensions)
    time_processes = np.array(time_processes)

    # Build the results.
    n_problems = len(problem_names)
    n_solvers = len(solvers)
    n_runs = feature.options[FeatureOption.N_RUNS]
    max_eval = max_eval_factor * np.max(problem_dimensions) if problem_dimensions.size > 0 else 1
    fun_histories = np.full((n_problems, n_solvers, n_runs, max_eval), np.nan)
    maxcv_histories = np.full((n_problems, n_solvers, n_runs, max_eval), np.nan)
    for i_problem, (fun_hist, maxcv_hist) in enumerate(zip(_fun_histories, _maxcv_histories)):
        max_eval = max_eval_factor * problem_dimensions[i_problem]
        fun_histories[i_problem, ..., :max_eval] = fun_hist
        maxcv_histories[i_problem, ..., :max_eval] = maxcv_hist
        if max_eval > 0:
            fun_histories[i_problem, ..., max_eval:] = fun_histories[i_problem, ..., max_eval - 1, np.newaxis]
            maxcv_histories[i_problem, ..., max_eval:] = maxcv_histories[i_problem, ..., max_eval - 1, np.newaxis]
    return problem_names, fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_dimensions, time_processes


def _solve_one_problem(problem_name, solvers, labels, feature, max_eval_factor, custom_problem_loader, cutest_problem_options, profile_options):
    """
    Solve a given problem.

    If problem_name is an instance of Problem, it is used directly. In this
    case, the problem options but contain the name of the problem. Otherwise,
    the problem is loaded from CUTEst.

    Notes
    -----
    Do not call this function concurrently with the same CUTEst problem name.
    """
    # Load the problem and return if it cannot be loaded.
    set_cutest_problem_options(**cutest_problem_options)
    if isinstance(problem_name, tuple):
        problem = custom_problem_loader(problem_name[1])
        problem_name = problem_name[0]
    else:
        try:
            problem = load_cutest_problem(problem_name)
        except ProblemError:
            return

    # Project the initial point if necessary.
    if profile_options[ProfileOption.PROJECT_X0]:
        problem.project_x0()

    # Evaluate the functions at the initial point.
    fun_init = problem.fun(problem.x0)
    maxcv_init = problem.maxcv(problem.x0)

    # Solve the problem with each solver.
    time_start = time.monotonic()
    n_solvers = len(solvers)
    n_runs = feature.options[FeatureOption.N_RUNS]
    max_eval = max_eval_factor * problem.n
    n_eval = np.zeros((n_solvers, n_runs), dtype=int)
    fun_histories = np.full((n_solvers, n_runs, max_eval), np.nan)
    fun_out = np.full((n_solvers, n_runs), np.nan)
    maxcv_histories = np.full((n_solvers, n_runs, max_eval), np.nan)
    maxcv_out = np.full((n_solvers, n_runs), np.nan)
    logger = get_logger(__name__)
    for i_solver in range(n_solvers):
        for i_run in range(n_runs):
            logger.info(f'Solving {problem_name} with {labels[i_solver]} (run {i_run + 1}/{n_runs}).')
            time_start_solver_run = time.monotonic()
            featured_problem = FeaturedProblem(problem, feature, max_eval, i_run)
            sig = signature(solvers[i_solver])
            if len(sig.parameters) not in [1, 2, 4, 8, 10]:
                raise ValueError(f'Unknown signature: {sig}.')
            with open(os.devnull, 'w') as devnull:
                with warnings.catch_warnings(), redirect_stdout(devnull), redirect_stderr(devnull):
                    warnings.filterwarnings('ignore')
                    try:
                        if len(sig.parameters) == 1:
                            x = solvers[i_solver](featured_problem)
                        elif len(sig.parameters) == 2:
                            x = solvers[i_solver](featured_problem.fun, featured_problem.x0)
                        elif len(sig.parameters) == 4:
                            x = solvers[i_solver](featured_problem.fun, featured_problem.x0, featured_problem.lb, featured_problem.ub)
                        elif len(sig.parameters) == 8:
                            x = solvers[i_solver](featured_problem.fun, featured_problem.x0, featured_problem.lb, featured_problem.ub, featured_problem.a_ub, featured_problem.b_ub, featured_problem.a_eq, featured_problem.b_eq)
                        else:
                            x = solvers[i_solver](featured_problem.fun, featured_problem.x0, featured_problem.lb, featured_problem.ub, featured_problem.a_ub, featured_problem.b_ub, featured_problem.a_eq, featured_problem.b_eq, featured_problem.c_ub, featured_problem.c_eq)
                        fun_out[i_solver, i_run] = problem.fun(featured_problem.raw_variables(x))
                        maxcv_out[i_solver, i_run] = problem.maxcv(featured_problem.raw_variables(x))
                        if featured_problem.type == 'unconstrained':
                            result = f'f = {fun_out[i_solver, i_run]:.4e}'
                        else:
                            result = f'f = {fun_out[i_solver, i_run]:.4e}, maxcv = {maxcv_out[i_solver, i_run]:.4e}'
                        logger.info(f'Results for {problem_name} with {labels[i_solver]} (run {i_run + 1}/{n_runs}): {result} ({time.monotonic() - time_start_solver_run:.2f} seconds).')
                    except Exception as exc:
                        logger.warning(f'An error occurred while solving {problem_name} with {labels[i_solver]}: {exc}.')
            n_eval[i_solver, i_run] = featured_problem.n_eval
            fun_histories[i_solver, i_run, :n_eval[i_solver, i_run]] = featured_problem.fun_hist[:n_eval[i_solver, i_run]]
            maxcv_histories[i_solver, i_run, :n_eval[i_solver, i_run]] = featured_problem.maxcv_hist[:n_eval[i_solver, i_run]]
            if n_eval[i_solver, i_run] > 0:
                fun_histories[i_solver, i_run, n_eval[i_solver, i_run]:] = fun_histories[i_solver, i_run, n_eval[i_solver, i_run] - 1]
                maxcv_histories[i_solver, i_run, n_eval[i_solver, i_run]:] = maxcv_histories[i_solver, i_run, n_eval[i_solver, i_run] - 1]
    return problem_name, fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem.n, time.monotonic() - time_start


def _compute_merit_values(fun_values, maxcv_values, maxcv_init):
    """
    Compute the merit function values.
    """
    maxcv_init = maxcv_init[(...,) + (np.newaxis,) * (fun_values.ndim - 1)]
    infeasibility_thresholds = np.maximum(1e-5, maxcv_init)
    is_infeasible = maxcv_values > infeasibility_thresholds
    is_almost_feasible = (1e-10 < maxcv_values) & (maxcv_values <= infeasibility_thresholds)
    merit_values = np.copy(fun_values)
    merit_values[is_infeasible | np.isnan(merit_values)] = np.inf
    merit_values[is_almost_feasible] += 1e5 * maxcv_values[is_almost_feasible]
    return merit_values


def _format_float_scientific_latex(x):
    """
    Format a floating-point number as scientific notation in LaTeX.
    """
    raw = np.format_float_scientific(x, trim='-', exp_digits=0)
    match = re.compile(r'^(?P<coefficient>[0-9]+(\.[0-9]+)?)e(?P<exponent>(-)?[0-9]+)$').match(raw)
    if not match:
        raise ValueError(f'Cannot format {x} as scientific notation.')
    if match.group('coefficient') == '1':
        return raw, f'10^{{{match.group("exponent")}}}'
    return raw, f'{match.group("coefficient")} \\times 10^{{{match.group("exponent")}}}'


def _draw_profiles(work_hist, work_out, problem_dimensions, labels, tolerance_label, i_tolerance, ax_summary_perf_hist, ax_summary_data_hist, ax_summary_log_ratio_hist, ax_summary_perf_out, ax_summary_data_out, ax_summary_log_ratio_out):
    n_solvers = work_hist.shape[1]

    # Create the individual figures.
    fig_perf_hist, ax_perf_hist = plt.subplots()
    fig_perf_out, ax_perf_out = plt.subplots()
    fig_data_hist, ax_data_hist = plt.subplots()
    fig_data_out, ax_data_out = plt.subplots()
    if n_solvers > 2:
        fig_log_ratio_hist, ax_log_ratio_hist = None, None
        fig_log_ratio_out, ax_log_ratio_out = None, None
    else:
        fig_log_ratio_hist, ax_log_ratio_hist = plt.subplots()
        fig_log_ratio_out, ax_log_ratio_out = plt.subplots()

    # Calculate the performance and data profiles.
    x_perf_hist, y_perf_hist, ratio_max_perf_hist, x_data_hist, y_data_hist, ratio_max_data_hist = _get_extended_performances_data_profile_axes(work_hist, problem_dimensions)
    x_perf_out, y_perf_out, ratio_max_perf_out, x_data_out, y_data_out, ratio_max_data_out = _get_extended_performances_data_profile_axes(work_out, problem_dimensions)
    x_perf_hist = np.log2(x_perf_hist)
    x_perf_out = np.log2(x_perf_out)
    ratio_max_perf_hist = np.log2(ratio_max_perf_hist)
    ratio_max_perf_out = np.log2(ratio_max_perf_out)
    x_data_hist = np.log2(1.0 + x_data_hist)
    x_data_out = np.log2(1.0 + x_data_out)
    ratio_max_data_hist = np.log2(1.0 + ratio_max_data_hist)
    ratio_max_data_out = np.log2(1.0 + ratio_max_data_out)

    # Draw the performance profiles
    _draw_performance_data_profiles(ax_perf_hist, x_perf_hist, y_perf_hist, labels)
    _draw_performance_data_profiles(ax_perf_out, x_perf_out, y_perf_out, labels)
    if ax_summary_perf_hist is not None:
        if i_tolerance == 0:
            _draw_performance_data_profiles(ax_summary_perf_hist[i_tolerance], x_perf_hist, y_perf_hist, labels)
        else:
            _draw_performance_data_profiles(ax_summary_perf_hist[i_tolerance], x_perf_hist, y_perf_hist)
    if ax_summary_perf_out is not None:
        _draw_performance_data_profiles(ax_summary_perf_out[i_tolerance], x_perf_out, y_perf_out)
    perf_formatter = FuncFormatter(_perf_formatter)
    ax_perf_hist.xaxis.set_major_formatter(perf_formatter)
    ax_perf_out.xaxis.set_major_formatter(perf_formatter)
    if ax_summary_perf_hist is not None:
        ax_summary_perf_hist[i_tolerance].xaxis.set_major_formatter(perf_formatter)
    if ax_summary_perf_out is not None:
        ax_summary_perf_out[i_tolerance].xaxis.set_major_formatter(perf_formatter)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning)
        ax_perf_hist.set_xlim(0.0, 1.05 * ratio_max_perf_hist)
        ax_perf_out.set_xlim(0.0, 1.05 * ratio_max_perf_out)
        if ax_summary_perf_hist is not None:
            ax_summary_perf_hist[i_tolerance].set_xlim(0.0, 1.05 * ratio_max_perf_hist)
        if ax_summary_perf_out is not None:
            ax_summary_perf_out[i_tolerance].set_xlim(0.0, 1.05 * ratio_max_perf_out)
    ax_perf_hist.set_xlabel('Performance ratio')
    ax_perf_out.set_xlabel('Performance ratio')
    if ax_summary_perf_hist is not None:
        ax_summary_perf_hist[i_tolerance].set_xlabel('Performance ratio')
    if ax_summary_perf_out is not None:
        ax_summary_perf_out[i_tolerance].set_xlabel('Performance ratio')
    ax_perf_hist.set_ylabel(f'Performance profiles {tolerance_label}')
    ax_perf_out.set_ylabel(f'Performance profiles {tolerance_label}')
    if ax_summary_perf_hist is not None:
        ax_summary_perf_hist[i_tolerance].set_ylabel(f'Performance profiles {tolerance_label}')
    if ax_summary_perf_out is not None:
        ax_summary_perf_out[i_tolerance].set_ylabel(f'Performance profiles {tolerance_label}')

    # Draw the data profiles.
    _draw_performance_data_profiles(ax_data_hist, x_data_hist, y_data_hist, labels)
    _draw_performance_data_profiles(ax_data_out, x_data_out, y_data_out, labels)
    if ax_summary_data_hist is not None:
        _draw_performance_data_profiles(ax_summary_data_hist[i_tolerance], x_data_hist, y_data_hist)
    if ax_summary_data_out is not None:
        _draw_performance_data_profiles(ax_summary_data_out[i_tolerance], x_data_out, y_data_out)
    data_formatter = FuncFormatter(_data_formatter)
    ax_data_hist.xaxis.set_major_formatter(data_formatter)
    ax_data_out.xaxis.set_major_formatter(data_formatter)
    if ax_summary_data_hist is not None:
        ax_summary_data_hist[i_tolerance].xaxis.set_major_formatter(data_formatter)
    if ax_summary_data_out is not None:
        ax_summary_data_out[i_tolerance].xaxis.set_major_formatter(data_formatter)
    ax_data_hist.set_xlim(0.0, 1.05 * ratio_max_data_hist)
    ax_data_out.set_xlim(0.0, 1.05 * ratio_max_data_out)
    if ax_summary_data_hist is not None:
        ax_summary_data_hist[i_tolerance].set_xlim(0.0, 1.05 * ratio_max_data_hist)
    if ax_summary_data_out is not None:
        ax_summary_data_out[i_tolerance].set_xlim(0.0, 1.05 * ratio_max_data_out)
    ax_data_hist.set_xlabel('Number of simplex gradients')
    ax_data_out.set_xlabel('Number of simplex gradients')
    if ax_summary_data_hist is not None:
        ax_summary_data_hist[i_tolerance].set_xlabel('Number of simplex gradients')
    if ax_summary_data_out is not None:
        ax_summary_data_out[i_tolerance].set_xlabel('Number of simplex gradients')
    ax_data_hist.set_ylabel(f'Data profiles {tolerance_label}')
    ax_data_out.set_ylabel(f'Data profiles {tolerance_label}')
    if ax_summary_data_hist is not None:
        ax_summary_data_hist[i_tolerance].set_ylabel(f'Data profiles {tolerance_label}')
    if ax_summary_data_out is not None:
        ax_summary_data_out[i_tolerance].set_ylabel(f'Data profiles {tolerance_label}')

    # Draw the log-ratio profiles.
    if n_solvers <= 2:
        _draw_log_ratio_profiles(ax_log_ratio_hist, np.copy(work_hist), labels)
        _draw_log_ratio_profiles(ax_log_ratio_out, np.copy(work_out), labels)
        if ax_summary_log_ratio_hist is not None:
            _draw_log_ratio_profiles(ax_summary_log_ratio_hist[i_tolerance], np.copy(work_hist), labels)
        if ax_summary_log_ratio_out is not None:
            _draw_log_ratio_profiles(ax_summary_log_ratio_out[i_tolerance], np.copy(work_out), labels)
        ax_log_ratio_hist.set_xlabel('Problem')
        ax_log_ratio_out.set_xlabel('Problem')
        if ax_summary_log_ratio_hist is not None:
            ax_summary_log_ratio_hist[i_tolerance].set_xlabel('Problem')
        if ax_summary_log_ratio_out is not None:
            ax_summary_log_ratio_out[i_tolerance].set_xlabel('Problem')
        ax_log_ratio_hist.set_ylabel(f'Log-ratio profiles {tolerance_label}')
        ax_log_ratio_out.set_ylabel(f'Log-ratio profiles {tolerance_label}')
        if ax_summary_log_ratio_hist is not None:
            ax_summary_log_ratio_hist[i_tolerance].set_ylabel(f'Log-ratio profiles {tolerance_label}')
        if ax_summary_log_ratio_out is not None:
            ax_summary_log_ratio_out[i_tolerance].set_ylabel(f'Log-ratio profiles {tolerance_label}')
    return fig_perf_hist, fig_perf_out, fig_data_hist, fig_data_out, fig_log_ratio_hist, fig_log_ratio_out


def _draw_performance_data_profiles(ax, x, y, labels=None):
    n_solvers = x.shape[1]
    n_runs = y.shape[2]
    y_mean = np.mean(y, 2)
    y_min = np.min(y, 2)
    y_max = np.max(y, 2)
    for i_solver in range(n_solvers):
        x_stairs = np.repeat(x[:, i_solver], 2)[1:]
        y_mean_stairs = np.repeat(y_mean[:, i_solver], 2)[:-1]
        y_min_stairs = np.repeat(y_min[:, i_solver], 2)[:-1]
        y_max_stairs = np.repeat(y_max[:, i_solver], 2)[:-1]
        if labels is not None:
            ax.plot(x_stairs, y_mean_stairs, label=labels[i_solver])
        else:
            ax.plot(x_stairs, y_mean_stairs)
        if n_runs > 1:
            ax.fill_between(x_stairs, y_min_stairs, y_max_stairs, alpha=0.2)
    ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_major_locator(MaxNLocator(5, prune='lower'))
    ax.yaxis.set_minor_locator(MaxNLocator(10))
    ax.set_ylim(0.0, 1.0)
    ax.tick_params(which='both', direction='in')
    if labels is not None:
        ax.legend(loc='lower right')


def _draw_log_ratio_profiles(ax, work, labels):
    n_problems, n_solvers, n_runs = work.shape
    work_flat = np.reshape(np.swapaxes(work, 1, 2), (n_problems * n_runs, n_solvers))
    log_ratio = np.full(n_problems * n_runs, np.nan)
    log_ratio_finite = np.isfinite(work_flat[:, 0]) & np.isfinite(work_flat[:, 1])
    log_ratio[log_ratio_finite] = np.log2(work_flat[log_ratio_finite, 0]) - np.log2(work_flat[log_ratio_finite, 1])
    ratio_max = np.max(np.abs(log_ratio[log_ratio_finite]), initial=np.finfo(float).eps)
    log_ratio[np.isnan(work_flat[:, 0]) & np.isfinite(work_flat[:, 1])] = 2.0 * ratio_max
    log_ratio[np.isfinite(work_flat[:, 0]) & np.isnan(work_flat[:, 1])] = -2.0 * ratio_max
    log_ratio[np.isnan(work_flat[:, 0]) & np.isnan(work_flat[:, 1])] = 0.0
    log_ratio = np.sort(log_ratio)

    x = np.arange(1, n_problems * n_runs + 1)
    ax.bar(x[log_ratio < 0], log_ratio[log_ratio < 0])
    ax.bar(x[log_ratio > 0], log_ratio[log_ratio > 0])
    ax.text((n_problems * n_runs + 1) / 2, -ratio_max, labels[0], horizontalalignment='center', verticalalignment='bottom')
    ax.text((n_problems * n_runs + 1) / 2, ratio_max, labels[1], horizontalalignment='center', verticalalignment='top')
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning)
        ax.set_xlim(0.5, n_problems * n_runs + 0.5)
    ax.set_ylim(-1.1 * ratio_max, 1.1 * ratio_max)
    ax.tick_params(which='both', direction='in')


def _get_extended_performances_data_profile_axes(work, problem_dimensions):
    n_problems, n_solvers, n_runs = work.shape
    x_perf, y_perf, ratio_max_perf = _get_performance_data_profile_axes(work, lambda i_problem, i_run: np.nanmin(work[i_problem, :, i_run], initial=np.inf))
    x_perf[np.isinf(x_perf)] = ratio_max_perf ** 2.0
    x_perf = np.vstack([np.ones((1, n_solvers)), x_perf])
    y_perf = np.vstack([np.zeros((1, n_solvers, n_runs)), y_perf])
    if n_problems > 0:
        x_perf = np.vstack([x_perf, np.full((1, n_solvers), ratio_max_perf ** 2.0)])
        y_perf = np.vstack([y_perf, y_perf[-1, np.newaxis, :, :]])
    x_data, y_data, ratio_max_data = _get_performance_data_profile_axes(work, lambda i_problem, i_run: problem_dimensions[i_problem] + 1)
    x_data[np.isinf(x_data)] = ratio_max_data ** 2.0 - 1.0
    x_data = np.vstack([np.zeros((1, n_solvers)), x_data])
    y_data = np.vstack([np.zeros((1, n_solvers, n_runs)), y_data])
    if n_problems > 0:
        x_data = np.vstack([x_data, np.full((1, n_solvers), ratio_max_data ** 2.0 - 1.0)])
        y_data = np.vstack([y_data, y_data[-1, np.newaxis, :, :]])
    return x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data


def _get_performance_data_profile_axes(work, denominator):
    """
    Calculate the axes of the performance and data profiles.
    """
    n_problems, n_solvers, n_runs = work.shape

    # Calculate the x-axis values.
    x = np.full((n_runs, n_problems, n_solvers), np.nan)
    for i_run in range(n_runs):
        for i_problem in range(n_problems):
            x[i_run, i_problem, :] = work[i_problem, :, i_run] / denominator(i_problem, i_run)
    ratio_max = np.nanmax(x, initial=np.finfo(float).eps)
    x[np.isnan(x)] = np.inf
    x = np.sort(x, 1)
    x = np.reshape(x, (n_problems * n_runs, n_solvers))
    sort_x = np.argsort(x, 0, 'stable')
    x = np.take_along_axis(x, sort_x, 0)

    # Calculate the y-axis values.
    y = np.full((n_problems * n_runs, n_solvers, n_runs), np.nan)
    for i_solver in range(n_solvers):
        for i_run in range(n_runs):
            if n_problems > 0:
                y[i_run * n_problems:(i_run + 1) * n_problems, i_solver, i_run] = np.linspace(1 / n_problems, 1.0, n_problems)
                y[:, i_solver, i_run] = np.take_along_axis(y[:, i_solver, i_run], sort_x[:, i_solver], 0)
            for i_problem in range(n_problems * n_runs):
                if np.isnan(y[i_problem, i_solver, i_run]):
                    y[i_problem, i_solver, i_run] = y[i_problem - 1, i_solver, i_run] if i_problem > 0 else 0.0
    return x, y, ratio_max


def _perf_formatter(x, _):
    if x.is_integer():
        return str(int(2 ** x))
    else:
        return f'$2^{{{f"{x:.8f}".rstrip("0").rstrip(".")}}}$'


def _data_formatter(x, _):
    if x.is_integer():
        return str(int(2 ** x - 1))
    else:
        return f'$2^{{{f"{x:.8f}".rstrip("0").rstrip(".")}}}-1$'

def _check_validity_problem_options(problem_options):
    """
    Check the validity of the problem options.
    """
    if ProblemOption.PLIBS in problem_options:
        if not isinstance(problem_options[ProblemOption.PLIBS], list):
            problem_options[ProblemOption.PLIBS] = [problem_options[ProblemOption.PLIBS]]
        mydir = os.path.dirname(os.path.abspath(__file__))
        problem_dir = os.path.abspath(os.path.join(mydir, '..', '..', 'problems'))
        subfolders = [name for name in os.listdir(problem_dir) if os.path.isdir(os.path.join(problem_dir, name))]
        # Remove '__pycache__' and folders starting with '.'
        subfolder_names = [name for name in subfolders if not name.startswith('.') and name != '__pycache__']
        # Check if plibs is empty or any element is not str or not in subfolder_names
        if len(problem_options[ProblemOption.PLIBS]) == 0:
            raise ValueError(f'Option {ProblemOption.PLIBS} cannot be an empty list.')
        elif any(not isinstance(plib, str) for plib in problem_options[ProblemOption.PLIBS]):
            raise TypeError(f'Option {ProblemOption.PLIBS} must be a string or a list of strings.')
        elif any(plib not in subfolder_names for plib in problem_options[ProblemOption.PLIBS]):
            raise ValueError(f'Option {ProblemOption.PLIBS} contains invalid problem libraries. Available libraries: {subfolder_names}.')

    if ProblemOption.PTYPE in problem_options:
        if not isinstance(problem_options[ProblemOption.PTYPE], str):
            raise TypeError(f'Option {ProblemOption.PTYPE} must be a string.')
        if any(t not in 'ubln' for t in problem_options[ProblemOption.PTYPE]):
            raise ValueError(f"Option {ProblemOption.PTYPE} should be a string containing only the characters 'u', 'b', 'l', and 'n'.")

    if ProblemOption.MINDIM in problem_options:
        if isinstance(problem_options[ProblemOption.MINDIM], float) and problem_options[ProblemOption.MINDIM].is_integer():
            problem_options[ProblemOption.MINDIM] = int(problem_options[ProblemOption.MINDIM])
        if not isinstance(problem_options[ProblemOption.MINDIM], int):
            raise TypeError(f'Option {ProblemOption.MINDIM} must be an integer.')
        if problem_options[ProblemOption.MINDIM] < 1:
            raise ValueError(f'Option {ProblemOption.MINDIM} must be at least 1.')
    if ProblemOption.MAXDIM in problem_options:
        if not np.isinf(problem_options[ProblemOption.MAXDIM]):
            if isinstance(problem_options[ProblemOption.MAXDIM], float) and problem_options[ProblemOption.MAXDIM].is_integer():
                problem_options[ProblemOption.MAXDIM] = int(problem_options[ProblemOption.MAXDIM])
            if not isinstance(problem_options[ProblemOption.MAXDIM], int):
                raise TypeError(f'Option {ProblemOption.MAXDIM} must be an integer or np.inf.')
            if problem_options[ProblemOption.MAXDIM] < 1:
                raise ValueError(f'Option {ProblemOption.MAXDIM} must be at least 1 or np.inf.')
    if ProblemOption.MINDIM in problem_options and ProblemOption.MAXDIM in problem_options:
        if problem_options[ProblemOption.MINDIM] > problem_options[ProblemOption.MAXDIM]:
            raise ValueError(f'Option {ProblemOption.MINDIM} cannot be larger than option {ProblemOption.MAXDIM}.')

    if ProblemOption.MINB in problem_options:
        if isinstance(problem_options[ProblemOption.MINB], float) and problem_options[ProblemOption.MINB].is_integer():
            problem_options[ProblemOption.MINB] = int(problem_options[ProblemOption.MINB])
        if not isinstance(problem_options[ProblemOption.MINB], int):
            raise TypeError(f'Option {ProblemOption.MINB} must be an integer.')
        if problem_options[ProblemOption.MINB] < 0:
            raise ValueError(f'Option {ProblemOption.MINB} must be nonnegative.')
    if ProblemOption.MAXB in problem_options:
        if not np.isinf(problem_options[ProblemOption.MAXB]):
            if isinstance(problem_options[ProblemOption.MAXB], float) and problem_options[ProblemOption.MAXB].is_integer():
                problem_options[ProblemOption.MAXB] = int(problem_options[ProblemOption.MAXB])
            if not isinstance(problem_options[ProblemOption.MAXB], int):
                raise TypeError(f'Option {ProblemOption.MAXB} must be an integer or np.inf.')
            if problem_options[ProblemOption.MAXB] < 0:
                raise ValueError(f'Option {ProblemOption.MAXB} must be nonnegative or np.inf.')
    if ProblemOption.MINB in problem_options and ProblemOption.MAXB in problem_options:
        if problem_options[ProblemOption.MINB] > problem_options[ProblemOption.MAXB]:
            raise ValueError(f'Option {ProblemOption.MINB} cannot be larger than option {ProblemOption.MAXB}.')

    if ProblemOption.MINLCON in problem_options:
        if isinstance(problem_options[ProblemOption.MINLCON], float) and problem_options[ProblemOption.MINLCON].is_integer():
            problem_options[ProblemOption.MINLCON] = int(problem_options[ProblemOption.MINLCON])
        if not isinstance(problem_options[ProblemOption.MINLCON], int):
            raise TypeError(f'Option {ProblemOption.MINLCON} must be an integer.')
        if problem_options[ProblemOption.MINLCON] < 0:
            raise ValueError(f'Option {ProblemOption.MINLCON} must be nonnegative.')
    if ProblemOption.MAXLCON in problem_options:
        if not np.isinf(problem_options[ProblemOption.MAXLCON]):
            if isinstance(problem_options[ProblemOption.MAXLCON], float) and problem_options[ProblemOption.MAXLCON].is_integer():
                problem_options[ProblemOption.MAXLCON] = int(problem_options[ProblemOption.MAXLCON])
            if not isinstance(problem_options[ProblemOption.MAXLCON], int):
                raise TypeError(f'Option {ProblemOption.MAXLCON} must be an integer or np.inf.')
            if problem_options[ProblemOption.MAXLCON] < 0:
                raise ValueError(f'Option {ProblemOption.MAXLCON} must be nonnegative or np.inf.')
    if ProblemOption.MINLCON in problem_options and ProblemOption.MAXLCON in problem_options:
        if problem_options[ProblemOption.MINLCON] > problem_options[ProblemOption.MAXLCON]:
            raise ValueError(f'Option {ProblemOption.MINLCON} cannot be larger than option {ProblemOption.MAXLCON}.')

    if ProblemOption.MINNLCON in problem_options:
        if isinstance(problem_options[ProblemOption.MINNLCON], float) and problem_options[ProblemOption.MINNLCON].is_integer():
            problem_options[ProblemOption.MINNLCON] = int(problem_options[ProblemOption.MINNLCON])
        if not isinstance(problem_options[ProblemOption.MINNLCON], int):
            raise TypeError(f'Option {ProblemOption.MINNLCON} must be an integer.')
        if problem_options[ProblemOption.MINNLCON] < 0:
            raise ValueError(f'Option {ProblemOption.MINNLCON} must be nonnegative.')
    if ProblemOption.MAXNLCON in problem_options:
        if not np.isinf(problem_options[ProblemOption.MAXNLCON]):
            if isinstance(problem_options[ProblemOption.MAXNLCON], float) and problem_options[ProblemOption.MAXNLCON].is_integer():
                problem_options[ProblemOption.MAXNLCON] = int(problem_options[ProblemOption.MAXNLCON])
            if not isinstance(problem_options[ProblemOption.MAXNLCON], int):
                raise TypeError(f'Option {ProblemOption.MAXNLCON} must be an integer or np.inf.')
            if problem_options[ProblemOption.MAXNLCON] < 0:
                raise ValueError(f'Option {ProblemOption.MAXNLCON} must be nonnegative or np.inf.')
    if ProblemOption.MINNLCON in problem_options and ProblemOption.MAXNLCON in problem_options:
        if problem_options[ProblemOption.MINNLCON] > problem_options[ProblemOption.MAXNLCON]:
            raise ValueError(f'Option {ProblemOption.MINNLCON} cannot be larger than option {ProblemOption.MAXNLCON}.')

    if ProblemOption.MINCON in problem_options:
        if isinstance(problem_options[ProblemOption.MINCON], float) and problem_options[ProblemOption.MINCON].is_integer():
            problem_options[ProblemOption.MINCON] = int(problem_options[ProblemOption.MINCON])
        if not isinstance(problem_options[ProblemOption.MINCON], int):
            raise TypeError(f'Option {ProblemOption.MINCON} must be an integer.')
        if problem_options[ProblemOption.MINCON] < 0:
            raise ValueError(f'Option {ProblemOption.MINCON} must be nonnegative.')
    if ProblemOption.MAXCON in problem_options:
        if not np.isinf(problem_options[ProblemOption.MAXCON]):
            if isinstance(problem_options[ProblemOption.MAXCON], float) and problem_options[ProblemOption.MAXCON].is_integer():
                problem_options[ProblemOption.MAXCON] = int(problem_options[ProblemOption.MAXCON])
            if not isinstance(problem_options[ProblemOption.MAXCON], int):
                raise TypeError(f'Option {ProblemOption.MAXCON} must be an integer or np.inf.')
            if problem_options[ProblemOption.MAXCON] < 0:
                raise ValueError(f'Option {ProblemOption.MAXCON} must be nonnegative or np.inf.')
    if ProblemOption.MINCON in problem_options and ProblemOption.MAXCON in problem_options:
        if problem_options[ProblemOption.MINCON] > problem_options[ProblemOption.MAXCON]:
            raise ValueError(f'Option {ProblemOption.MINCON} cannot be larger than option {ProblemOption.MAXCON}.')

    if ProblemOption.EXCLUDELIST in problem_options:
        if not isinstance(problem_options[ProblemOption.EXCLUDELIST], list):
            problem_options[ProblemOption.EXCLUDELIST] = [problem_options[ProblemOption.EXCLUDELIST]]
        if any(not isinstance(name, str) for name in problem_options[ProblemOption.EXCLUDELIST]):
            raise TypeError(f'Option {ProblemOption.EXCLUDELIST} must be a string or a list of strings.')

    if ProblemOption.PROBLEM_NAMES in problem_options:
        if not isinstance(problem_options[ProblemOption.PROBLEM_NAMES], list):
            problem_options[ProblemOption.PROBLEM_NAMES] = [problem_options[ProblemOption.PROBLEM_NAMES]]
        if any(not isinstance(name, str) for name in problem_options[ProblemOption.PROBLEM_NAMES]):
            raise TypeError(f'Option {ProblemOption.PROBLEM_NAMES} must be a string or a list of strings.')

    return problem_options

def _check_validity_profile_options(solvers, profile_options):
    """
    Check the validity of the profile options.
    """
    if ProfileOption.N_JOBS in profile_options:
        if isinstance(profile_options[ProfileOption.N_JOBS], float) and profile_options[ProfileOption.N_JOBS].is_integer():
            profile_options[ProfileOption.N_JOBS] = int(profile_options[ProfileOption.N_JOBS])
        if not isinstance(profile_options[ProfileOption.N_JOBS], int):
            raise TypeError(f'Option {ProfileOption.N_JOBS} must be an integer.')
        if profile_options[ProfileOption.N_JOBS] < 1:
            profile_options[ProfileOption.N_JOBS] = 1
            logger = get_logger(__name__)
            logger.warning(f'Option {ProfileOption.N_JOBS} is set to 1 because it cannot be smaller than 1.')

    if ProfileOption.SEED in profile_options:
        if isinstance(profile_options[ProfileOption.SEED], float) and profile_options[ProfileOption.SEED].is_integer():
            profile_options[ProfileOption.SEED] = int(profile_options[ProfileOption.SEED])
        if not isinstance(profile_options[ProfileOption.SEED], int):
            raise TypeError(f'Option {ProfileOption.SEED} must be an integer.')
        if profile_options[ProfileOption.SEED] < 0:
            raise ValueError(f'Option {ProfileOption.SEED} must be nonnegative.')


