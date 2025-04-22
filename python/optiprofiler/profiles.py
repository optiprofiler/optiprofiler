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
from .problems import FeaturedProblem, get_cutest_problem_options, set_cutest_problem_options, load_cutest_problem
from .utils import FeatureName, ProfileOption, FeatureOption, ProblemError, get_logger


def run_benchmark(solvers, labels=(), cutest_problem_names=(), custom_problem_loader=None, custom_problem_names=(), feature_name=FeatureName.PLAIN.value, **kwargs):
    """
    Run the benchmark.

    This function generates performance/data profiles [1]_, [2]_, [4]_ and
    log-ratio profiles [3]_, [5]_ for the given solvers on the given problems.

    .. caution::

        To use CUTEst problems in your benchmark, you must first install
        `PyCUTEst <https://jfowkes.github.io/pycutest/>`_. Follow the
        instructions carefully, as the CUTEst library must be installed in
        order to use `PyCUTEst <https://jfowkes.github.io/pycutest/>`_.

    Parameters
    ----------
    solvers : list of callable
        Solvers to benchmark. Each solver must be a callable, as follows. For
        unconstrained problems, the signature of the callable can be

            ``solver(fun, x0) -> numpy.ndarray, shape (n,)``

        where ``fun`` is the objective function and ``x0`` is the initial
        point. The solver must return the best point found. The objective
        function returns a scalar and should be minimized. For
        bound-constrained problems, the signature of the callable can be

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
    labels : list of str, optional
        Labels of the solvers in the plots. By default, the labels are the
        names of the callables in `solvers`.
    cutest_problem_names : list of str, optional
        Names of the CUTEst problems to use in the benchmark. Each of these
        problems will be loaded from CUTEst with their default parameters. If a
        problem cannot be loaded, it is ignored. You can use the function
        `find_cutest_problems` to obtain a list of available CUTEst problems.
    custom_problem_loader : callable, optional
        Callable to load custom problems. The signature of the callable must be

            ``custom_problem_loader(problem_name) -> Problem``

        where ``problem_name`` is the name of the problem to load.
    custom_problem_names : list of str, optional
        Names of the custom problems to use in the benchmark. Each of these
        problems will be loaded using the custom problem loader.
    feature_name : str, optional
        Name of the feature to use in the benchmark. Available features are

            ``'plain'`` :
                The problems are not modified.
            ``'noisy'`` :
                The objective function are perturbed with noise.
            ``'permuted'`` :
                The variables are randomly permuted.
            ``'perturbed_x0'`` :
                The initial point is randomly perturbed.
            ``'random_nan'`` :
                NaN values are randomly returned by the objective function.
            ``'truncated'`` :
                The objective function values are truncated.
            ``'unrelaxable_constraints'`` :
                The objective function values are infinite if constraints are not satisfied.

    Other Parameters
    ----------------
    n_jobs : int, optional
        Number of parallel workers to use.
    benchmark_id : str, optional
        Identifier of the benchmark to use in the output paths.
    project_x0 : bool, optional
        Whether to project the initial points of all the problems included in
        the benchmark onto their feasible set.
    summarize_performance_profiles : bool, optional
        Whether to plot the performance profiles on the summary.
    summarize_data_profiles : bool, optional
        Whether to plot the data profiles on the summary.
    summarize_log_ratio_profiles : bool, optional
        Whether to plot the log-ratio profiles on the summary.
    n_runs : int, optional
        Number of runs to perform for each solver on each problem.
    distribution : callable, optional
        Only used when ``feature_name='noisy'``. Distribution of the noise to
        include in the objective function values.
    perturbed_trailing_zeros : bool, optional
        Only used when ``feature_name='truncated'``. Whether to perturb the
        trailing zeros in the objective function values.
    rate_nan : float, optional
        Only used when ``feature_name='random_nan'``. Rate of NaN values
        occurring in the objective function values.
    significant_digits : int, optional
        Only used when ``feature_name='truncated'``. Number of significant
        digits to keep in the objective function values.
    type : str, optional
        Only used when ``feature_name='noisy'``. Type of the noise to include
        in the objective function values. Available types are ``'absolute'``
        and ``'relative'``.
    unrelaxable_bounds : bool, optional
        Only used when ``feature_name='unrelaxable_constraints'``. Whether to
        make the bounds unrelaxable.
    unrelaxable_linear_constraints : bool, optional
        Only used when ``feature_name='unrelaxable_constraints'``. Whether to
        make the linear constraints unrelaxable.
    unrelaxable_nonlinear_constraints : bool, optional
        Only used when ``feature_name='unrelaxable_constraints'``. Whether to
        make the nonlinear constraints unrelaxable.

    Raises
    ------
    TypeError
        If an argument received an invalid value.
    ValueError
        If the arguments are inconsistent.

    See Also
    --------
    find_cutest_problems : Find names of CUTEst problems.
    Problem : Representation of optimization problems.

    Notes
    -----
    The current version of `optiprofiler` only supports derivative-free
    optimization solvers.

    References
    ----------
    .. [1] E. D. Dolan and J. J. Moré. Benchmarking optimization software with
           performance profiles. *Math. Program.*, 91(2):201–213, 2002.
           `doi:10.1007/s101070100263
           <https://doi.org/10.1007/s101070100263>`_.
    .. [2] N. Gould and J. Scott. A note on performance profiles for
           benchmarking software. *ACM Trans. Math. Software*, 43(2):15:1–5,
           2016. `doi:10.1145/2950048 <https://doi.org/10.1145/2950048>`_.
    .. [3] J. L. Morales. A numerical study of limited memory BFGS methods.
           *Appl. Math. Lett.*, 15(4):481–487, 2002.
           `doi:10.1016/S0893-9659(01)00162-8
           <https://doi.org/10.1016/S0893-9659(01)00162-8>`_.
    .. [4] J. J. Moré and S. M. Wild. Benchmarking derivative-free optimization
           algorithms. *SIAM J. Optim.*, 20(1):172–191, 2009.
           `doi:10.1137/080724083 <https://doi.org/10.1137/080724083>`_.
    .. [5] H.-J. M. Shi, M. Q. Xuan, F. Oztoprak, and J. Nocedal. On the
           numerical performance of finite-difference-based methods for
           derivative-free optimization. *Optim. Methods Softw.*,
           38(2):289–311, 2023. `doi:10.1080/10556788.2022.2121832
           <https://doi.org/10.1080/10556788.2022.2121832>`_.

    Examples
    --------
    To do.
    """
    logger = get_logger(__name__)

    # Preprocess the solvers.
    if not hasattr(solvers, '__len__') or not all(callable(solver) for solver in solvers):
        raise TypeError('The solvers must be a list of callables.')
    if len(solvers) < 2:
        raise ValueError('At least two solvers must be given.')
    solvers = list(solvers)

    # Preprocess the labels.
    if not hasattr(labels, '__len__') or not all(isinstance(label, str) for label in labels):
        raise TypeError('The labels must be a list of strings.')
    if len(labels) not in [0, len(solvers)]:
        raise ValueError('The number of labels must equal the number of solvers.')
    labels = list(labels)
    if len(labels) == 0:
        labels = [solver.__name__ for solver in solvers]

    # Preprocess the CUTEst problem names.
    # N.B.: Duplicate names are and MUST BE removed.
    if not hasattr(cutest_problem_names, '__len__') or not all(isinstance(problem_name, str) for problem_name in cutest_problem_names):
        raise TypeError('The CUTEst problem names must be a list of strings.')
    cutest_problem_names = list(set(problem_name.upper() for problem_name in cutest_problem_names))

    # Preprocess the custom problems.
    if custom_problem_loader is not None and not callable(custom_problem_loader):
        raise TypeError('The custom problem loader must be a callable.')
    if custom_problem_loader is not None:
        sig = signature(custom_problem_loader)
        if len(sig.parameters) != 1:
            raise ValueError('The custom problem loader must take exactly one argument.')
    if not hasattr(custom_problem_names, '__len__') or not all(isinstance(problem, str) for problem in custom_problem_names):
        raise TypeError('The custom problem names must be a list of strings.')
    custom_problem_names = list(custom_problem_names)
    if custom_problem_loader is None and len(custom_problem_names) > 0:
        raise ValueError('A custom problem loader must be given to load custom problems.')

    # Check that the number of problems is satisfactory.
    if len(cutest_problem_names) + len(custom_problem_names) < 1:
        raise ValueError('At least one problem must be given.')

    # Get the different options from the keyword arguments.
    feature_options = {}
    profile_options = {}
    for key, value in kwargs.items():
        if key in FeatureOption.__members__.values():
            feature_options[key] = value
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
