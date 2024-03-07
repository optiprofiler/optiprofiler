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
from pypdf import PdfReader, PdfWriter, Transformation

from .features import Feature
from .problems import Problem, FeaturedProblem, get_cutest_problem_options, set_cutest_problem_options, load_cutest_problem
from .utils import FeatureName, ProfileOption, FeatureOption, ProblemError, get_logger


def run_benchmark(solvers, labels=(), cutest_problem_names=(), extra_problems=(), feature_name=FeatureName.PLAIN.value, **kwargs):
    """
    Create the benchmark profiles.

    The benchmark profiles include the performance and data profiles [1]_,
    [2]_, [4]_, and the log-ratio profiles [3]_, [5]_. The log-ratio profiles
    are available only when there are exactly two solvers.

    .. caution::

        To use CUTEst problems in your benchmark, you must first install
        `PyCUTEst <https://jfowkes.github.io/pycutest/>`_. Follow the
        instructions carefully, as the CUTEst library must be installed in
        order to use `PyCUTEst <https://jfowkes.github.io/pycutest/>`_.

    Parameters
    ----------
    solvers : list of callable
        Solvers to benchmark. Each solver must be a callable, as follows. For
        unconstrained problems, the signature of the callable must be

            ``solver(fun, x0)``

        where ``fun`` is the objective function and ``x0`` is the initial point.
        The objective function returns a scalar and should be minimized. For
        bound-constrained problems, the signature of the callable must be

            ``solver(fun, x0, lb, ub)``

        where ``lb`` and ``ub`` are the lower and upper bounds, respectively.
        For linearly constrained problems, the signature of the callable must
        be

            ``solver(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq)``

        where ``a_ub @ x <= b_ub`` and ``a_eq @ x == b_eq`` form the linear
        inequality and equality constraints, respectively. For nonlinearly
        constrained problems, the signature of the callable must be

            ``solver(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq, c_ub, c_eq)``

        where ``c_ub(x) <= 0`` and ``c_eq(x) == 0`` form the nonlinear
        inequality and equality constraints, respectively.

        All vectors and matrices mentioned above are `numpy.ndarray`.
    labels : list of str, optional
        Labels of the solvers in the plots. By default, the labels are the
        names of the callables in `solvers`.
    cutest_problem_names : list of str, optional
        Names of the CUTEst problems to use in the benchmark. Each of these
        problems will be loaded from CUTEst with their default parameters. If a
        problem cannot be loaded, it is ignored. You can use the function
        `find_cutest_problems` to obtain a list of available CUTEst problems.
    extra_problems : list of `Problem`, optional
        Extra problems to use in the benchmark. If you do not want to use
        CUTEst problems, you can use this argument to provide your own
        problems.
    feature_name : str, optional
        Name of the feature to use in the benchmark.

    Other Parameters
    ----------------
    To do.

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

    # Preprocess the extra problems.
    if not hasattr(extra_problems, '__len__') or not all(isinstance(problem, Problem) for problem in extra_problems):
        raise TypeError('The extra problems must be a list of problems.')
    extra_problems = list(extra_problems)

    # Check that the number of problems is satisfactory.
    if len(cutest_problem_names) + len(extra_problems) < 1:
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
    profile_options.setdefault(ProfileOption.SUMMARIZE.value, True)
    profile_options.setdefault(ProfileOption.PROJECT_X0.value, False)

    # Check whether the profile options are valid.
    if isinstance(profile_options[ProfileOption.N_JOBS], float) and profile_options[ProfileOption.N_JOBS].is_integer():
        profile_options[ProfileOption.N_JOBS] = int(profile_options[ProfileOption.N_JOBS])
    if profile_options[ProfileOption.N_JOBS] is not None and not (isinstance(profile_options[ProfileOption.N_JOBS], int) and profile_options[ProfileOption.N_JOBS] > 0):
        raise TypeError(f'Option {ProfileOption.N_JOBS} must be a positive integer or None.')
    if not isinstance(profile_options[ProfileOption.BENCHMARK_ID], str):
        raise TypeError(f'Option {ProfileOption.BENCHMARK_ID} must be a string.')
    if not isinstance(profile_options[ProfileOption.SUMMARIZE], bool):
        raise TypeError(f'Option {ProfileOption.SUMMARIZE} must be a boolean.')
    if not isinstance(profile_options[ProfileOption.PROJECT_X0], bool):
        raise TypeError(f'Option {ProfileOption.PROJECT_X0} must be a boolean.')

    # Build the feature.
    logger = get_logger(__name__)
    feature = Feature(feature_name, **feature_options)
    logger.info(f'Starting the computations of the "{feature.name}" profiles.')

    # Solve the problems.
    max_eval_factor = 500
    problem_names, fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_dimensions, time_processes = _solve_all_problems(cutest_problem_names, extra_problems, solvers, labels, feature, max_eval_factor, profile_options)
    merit_histories = _compute_merit_values(fun_histories, maxcv_histories, maxcv_init)
    merit_out = _compute_merit_values(fun_out, maxcv_out, maxcv_init)
    merit_init = _compute_merit_values(fun_init, maxcv_init, maxcv_init)

    # Determine the least merit value for each problem.
    merit_min = np.min(merit_histories, (1, 2, 3))
    if feature.is_stochastic:
        feature_plain = Feature('plain')
        logger.info(f'Starting the computations of the "plain" profiles.')
        _, fun_histories_plain, maxcv_histories_plain, _, _, _, _, _, _, time_processes_plain = _solve_all_problems(cutest_problem_names, extra_problems, solvers, labels, feature_plain, max_eval_factor, profile_options)
        time_processes += time_processes_plain
        merit_histories_plain = _compute_merit_values(fun_histories_plain, maxcv_histories_plain, maxcv_init)
        merit_min_plain = np.min(merit_histories_plain, (1, 2, 3))
        merit_min = np.minimum(merit_min, merit_min_plain)

    # Paths to the results.
    timestamp = datetime.now().astimezone().strftime('%Y-%m-%dT%H-%M-%S%Z')
    path_out = Path('out').resolve()
    path_results = path_out / feature.name / profile_options[ProfileOption.BENCHMARK_ID] / timestamp
    path_results.mkdir(parents=True, exist_ok=True)
    path_problems = path_results / 'problems.txt'
    path_perf_hist = path_results / 'perf_hist.pdf'
    path_perf_out = path_results / 'perf_out.pdf'
    path_data_hist = path_results / 'data_hist.pdf'
    path_data_out = path_results / 'data_out.pdf'
    path_log_ratio_hist = path_results / 'log-ratio_hist.pdf'
    path_log_ratio_out = path_results / 'log-ratio_out.pdf'
    path_summary = path_out / f'summary_{timestamp}.pdf'

    # Store the names of the problems.
    with path_problems.open('w') as f:
        f.write(os.linesep.join(f'{problem_name}: {time_process:.2f} seconds' for problem_name, time_process in sorted(zip(problem_names, time_processes))))

    # Set up matplotlib for plotting the profiles.
    prop_cycle = cycler(color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])
    prop_cycle += cycler(linestyle=[(0, ()), (0, (1, 1.5)), (0, (3, 1.5)), (0, (5, 1.5, 1, 1.5)), (0, (5, 1.5, 1, 1.5, 1, 1.5)), (0, (1, 3)), (0, (3, 3)), (0, (5, 3, 1, 3)), (0, (5, 3, 1, 3, 1, 3)), (0, (1, 4.5))])
    with plt.rc_context({
        'axes.prop_cycle': prop_cycle,
        'font.family': 'serif',
        'font.size': 16,
        'text.usetex': True if shutil.which('latex') else False,
    }):
        # Create the performance and data profiles.
        n_problems, n_solvers, n_runs, max_eval = merit_histories.shape
        tolerances = np.logspace(-1, -10, 10)
        pdf_perf_hist = backend_pdf.PdfPages(path_perf_hist)
        pdf_perf_out = backend_pdf.PdfPages(path_perf_out)
        pdf_data_hist = backend_pdf.PdfPages(path_data_hist)
        pdf_data_out = backend_pdf.PdfPages(path_data_out)
        pdf_log_ratio_hist = backend_pdf.PdfPages(path_log_ratio_hist, False)
        pdf_log_ratio_out = backend_pdf.PdfPages(path_log_ratio_out, False)
        for tolerance in tolerances:
            tolerance_str, tolerance_latex = _format_float_scientific_latex(tolerance)
            logger.info(f'Creating profiles for tolerance {tolerance_str}.')
            tolerance_label = f'($\\mathrm{{tol}} = {tolerance_latex}$)'

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

            # Draw the profiles.
            fig_perf_hist, fig_perf_out, fig_data_hist, fig_data_out, fig_log_ratio_hist, fig_log_ratio_out = _draw_profiles(work_hist, work_out, problem_dimensions, labels, tolerance_label)
            pdf_perf_hist.savefig(fig_perf_hist, bbox_inches='tight')
            pdf_perf_out.savefig(fig_perf_out, bbox_inches='tight')
            pdf_data_hist.savefig(fig_data_hist, bbox_inches='tight')
            pdf_data_out.savefig(fig_data_out, bbox_inches='tight')
            if fig_log_ratio_hist is not None:
                pdf_log_ratio_hist.savefig(fig_log_ratio_hist, bbox_inches='tight')
            if fig_log_ratio_out is not None:
                pdf_log_ratio_out.savefig(fig_log_ratio_out, bbox_inches='tight')

            # Close the figures.
            plt.close(fig_perf_hist)
            plt.close(fig_perf_out)
            plt.close(fig_data_hist)
            plt.close(fig_data_out)
            if fig_log_ratio_hist is not None:
                plt.close(fig_log_ratio_hist)
            if fig_log_ratio_out is not None:
                plt.close(fig_log_ratio_out)

        # Close the PDF files.
        pdf_perf_hist.close()
        pdf_perf_out.close()
        pdf_data_hist.close()
        pdf_data_out.close()
        pdf_log_ratio_hist.close()
        pdf_log_ratio_out.close()
        logger.info(f'Results stored in {path_out}.')

    # Create the summary PDF if necessary.
    if profile_options[ProfileOption.SUMMARIZE]:
        logger.info('Creating summary PDF.')
        path_profiles = [path_perf_hist, path_data_hist]
        reader_profiles = [PdfReader(path_perf_hist), PdfReader(path_data_hist)]
        if path_log_ratio_hist.is_file():
            path_profiles.append(path_log_ratio_hist)
            reader_profiles.append(PdfReader(path_log_ratio_hist))
        path_profiles.extend([path_perf_out, path_data_out])
        reader_profiles.extend([PdfReader(path_perf_out), PdfReader(path_data_out)])
        if path_log_ratio_out.is_file():
            path_profiles.append(path_log_ratio_out)
            reader_profiles.append(PdfReader(path_log_ratio_out))
        profile_number = len(path_profiles)
        page_number = max(len(reader.pages) for reader in reader_profiles)
        page_width = max(page.mediabox.width for reader in reader_profiles for page in reader.pages)
        page_height = max(page.mediabox.height for reader in reader_profiles for page in reader.pages)
        pdf_summary = PdfWriter()
        pdf_summary.add_blank_page(page_width * page_number, page_height * profile_number)
        for i_source, pdf_source in enumerate(reader_profiles):
            for i_page, page in enumerate(pdf_source.pages):
                pdf_summary.pages[0].merge_transformed_page(
                    page,
                    Transformation().translate(
                        page_width * i_page,
                        (profile_number - i_source - 1) * page_height,
                    ),
                )
        pdf_summary.write(path_summary)
        pdf_summary.close()


def _solve_all_problems(cutest_problem_names, extra_problems, solvers, labels, feature, max_eval_factor, profile_options):
    """
    Solve all problems in parallel.
    """
    problem_names = cutest_problem_names + [(f'EXTRA{i_problem}', extra_problem) for i_problem, extra_problem in enumerate(extra_problems)]
    cutest_problem_options = get_cutest_problem_options()

    # Solve all problems.
    logger = get_logger(__name__)
    logger.info('Entering the parallel section.')
    with Pool(profile_options[ProfileOption.N_JOBS]) as p:
        results = p.starmap(_solve_one_problem, [(problem_name, solvers, labels, feature, max_eval_factor, cutest_problem_options, profile_options) for problem_name in problem_names])
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


def _solve_one_problem(problem_name, solvers, labels, feature, max_eval_factor, cutest_problem_options, profile_options):
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
        problem = problem_name[1]
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
            if len(sig.parameters) not in [2, 4, 8, 10]:
                raise ValueError(f'Unknown signature: {sig}.')
            with open(os.devnull, 'w') as devnull:
                with warnings.catch_warnings(), redirect_stdout(devnull), redirect_stderr(devnull):
                    warnings.filterwarnings('ignore')
                    try:
                        if len(sig.parameters) == 2:
                            x = solvers[i_solver](featured_problem.fun, featured_problem.x0)
                        elif len(sig.parameters) == 4:
                            x = solvers[i_solver](featured_problem.fun, featured_problem.x0, featured_problem.lb, featured_problem.ub)
                        elif len(sig.parameters) == 8:
                            x = solvers[i_solver](featured_problem.fun, featured_problem.x0, featured_problem.lb, featured_problem.ub, featured_problem.a_ub, featured_problem.b_ub, featured_problem.a_eq, featured_problem.b_eq)
                        else:
                            x = solvers[i_solver](featured_problem.fun, featured_problem.x0, featured_problem.lb, featured_problem.ub, featured_problem.a_ub, featured_problem.b_ub, featured_problem.a_eq, featured_problem.b_eq, featured_problem.c_ub, featured_problem.c_eq)
                        fun_out[i_solver, i_run] = problem.fun(x)
                        maxcv_out[i_solver, i_run] = problem.maxcv(x)
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


def _draw_profiles(work_hist, work_out, problem_dimensions, labels, tolerance_label):
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

    # Draw the performance and data profiles.
    x_perf_hist, y_perf_hist, ratio_max_perf_hist, x_data_hist, y_data_hist, ratio_max_data_hist = _get_extended_performances_data_profile_axes(work_hist, problem_dimensions)
    x_perf_out, y_perf_out, ratio_max_perf_out, x_data_out, y_data_out, ratio_max_data_out = _get_extended_performances_data_profile_axes(work_out, problem_dimensions)
    _draw_performance_data_profiles(ax_perf_hist, x_perf_hist, y_perf_hist, labels)
    _draw_performance_data_profiles(ax_perf_out, x_perf_out, y_perf_out, labels)
    _draw_performance_data_profiles(ax_data_hist, x_data_hist, y_data_hist, labels)
    _draw_performance_data_profiles(ax_data_out, x_data_out, y_data_out, labels)
    ax_perf_hist.set_xscale('log', base=2)
    ax_perf_out.set_xscale('log', base=2)
    int_formatter = FuncFormatter(lambda x, _: f'{x:.0f}')
    ax_perf_hist.xaxis.set_major_formatter(int_formatter)
    ax_perf_out.xaxis.set_major_formatter(int_formatter)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning)
        ax_perf_hist.set_xlim(1.0, ratio_max_perf_hist ** 1.1)
        ax_perf_out.set_xlim(1.0, ratio_max_perf_out ** 1.1)
    ax_data_hist.set_xlim(0.0, 1.1 * ratio_max_data_hist)
    ax_data_out.set_xlim(0.0, 1.1 * ratio_max_data_out)
    ax_perf_hist.set_xlabel('Performance ratio')
    ax_perf_out.set_xlabel('Performance ratio')
    ax_data_hist.set_xlabel('Number of simplex gradients')
    ax_data_out.set_xlabel('Number of simplex gradients')
    ax_perf_hist.set_ylabel(f'Performance profiles {tolerance_label}')
    ax_perf_out.set_ylabel(f'Performance profiles {tolerance_label}')
    ax_data_hist.set_ylabel(f'Data profiles {tolerance_label}')
    ax_data_out.set_ylabel(f'Data profiles {tolerance_label}')

    # Draw the log-ratio profiles.
    if n_solvers <= 2:
        _draw_log_ratio_profiles(ax_log_ratio_hist, np.copy(work_hist), labels)
        _draw_log_ratio_profiles(ax_log_ratio_out, np.copy(work_out), labels)
        ax_log_ratio_hist.set_xlabel('Problem')
        ax_log_ratio_out.set_xlabel('Problem')
        ax_log_ratio_hist.set_ylabel(f'Log-ratio profiles {tolerance_label}')
        ax_log_ratio_out.set_ylabel(f'Log-ratio profiles {tolerance_label}')
    return fig_perf_hist, fig_perf_out, fig_data_hist, fig_data_out, fig_log_ratio_hist, fig_log_ratio_out


def _draw_performance_data_profiles(ax, x, y, labels):
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
        ax.plot(x_stairs, y_mean_stairs, label=labels[i_solver])
        if n_runs > 1:
            ax.fill_between(x_stairs, y_min_stairs, y_max_stairs, alpha=0.2)
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_major_locator(MaxNLocator(5, prune='lower'))
    ax.yaxis.set_minor_locator(MaxNLocator(10))
    ax.set_ylim(0.0, 1.0)
    ax.tick_params(which='both', direction='in')
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
    x_data[np.isinf(x_data)] = 2.0 * ratio_max_data
    x_data = np.vstack([np.zeros((1, n_solvers)), x_data])
    y_data = np.vstack([np.zeros((1, n_solvers, n_runs)), y_data])
    if n_problems > 0:
        x_data = np.vstack([x_data, np.full((1, n_solvers), 2.0 * ratio_max_data)])
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
