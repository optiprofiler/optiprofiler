import logging
import os
import platform
import sys
from enum import Enum
from importlib.metadata import PackageNotFoundError, version


class FeatureName(str, Enum):
    """
    Feature names.
    """
    PLAIN = 'plain'
    PERTURBED_X0 = 'perturbed_x0'
    NOISY = 'noisy'
    TRUNCATED = 'truncated'
    PERMUTED = 'permuted'
    LINEARLY_TRANSFORMED = 'linearly_transformed'
    RANDOM_NAN = 'random_nan'
    UNRELAXABLE_CONSTRAINTS = 'unrelaxable_constraints'
    NONQUANTIFIABLE_CONSTRAINTS = 'nonquantifiable_constraints'
    QUANTIZED = 'quantized'
    CUSTOM = 'custom'


class ProfileOption(str, Enum):
    """
    Profile's options.
    """
    N_JOBS = 'n_jobs'
    SEED = 'seed'
    BENCHMARK_ID = 'benchmark_id'
    SOLVER_NAMES = 'solver_names'
    SOLVER_ISRAND = 'solver_isrand'
    FEATURE_STAMP = 'feature_stamp'
    ERRORBAR_TYPE = 'errorbar_type'
    SAVEPATH = 'savepath'
    MAX_TOL_ORDER = 'max_tol_order'
    MAX_EVAL_FACTOR = 'max_eval_factor'
    MERIT_FUN = 'merit_fun'
    PROJECT_X0 = 'project_x0'
    RUN_PLAIN = 'run_plain'
    SCORE_ONLY = 'score_only'
    SUMMARIZE_PERFORMANCE_PROFILES = 'summarize_performance_profiles'
    SUMMARIZE_DATA_PROFILES = 'summarize_data_profiles'
    SUMMARIZE_LOG_RATIO_PROFILES = 'summarize_log_ratio_profiles'
    SUMMARIZE_OUTPUT_BASED_PROFILES = 'summarize_output_based_profiles'
    SILENT = 'silent'
    SOLVER_VERBOSE = 'solver_verbose'
    SEMILOGX = 'semilogx'
    NORMALIZED_SCORES = 'normalized_scores'
    SCORE_WEIGHT_FUN = 'score_weight_fun'
    SCORE_FUN = 'score_fun'
    LOAD = 'load'
    SOLVERS_TO_LOAD = 'solvers_to_load'
    LINE_COLORS = 'line_colors'
    LINE_STYLES = 'line_styles'
    LINE_WIDTHS = 'line_widths'
    BAR_COLORS = 'bar_colors'
    XLABEL_PERFORMANCE_PROFILE = 'xlabel_performance_profile'
    YLABEL_PERFORMANCE_PROFILE = 'ylabel_performance_profile'
    XLABEL_DATA_PROFILE = 'xlabel_data_profile'
    YLABEL_DATA_PROFILE = 'ylabel_data_profile'
    XLABEL_LOG_RATIO_PROFILE = 'xlabel_log_ratio_profile'
    YLABEL_LOG_RATIO_PROFILE = 'ylabel_log_ratio_profile'


class ProblemOption(str, Enum):
    """
    Options for selecting problems from problem libraries.
    """
    PLIBS = 'plibs'
    PTYPE = 'ptype'
    MINDIM = 'mindim'
    MAXDIM = 'maxdim'
    MINB = 'minb'
    MAXB = 'maxb'
    MINLCON = 'minlcon'
    MAXLCON = 'maxlcon'
    MINNLCON = 'minnlcon'
    MAXNLCON = 'maxnlcon'
    MINCON = 'mincon'
    MAXCON = 'maxcon'
    EXCLUDELIST = 'excludelist'
    PROBLEM_NAMES = 'problem_names'


class FeatureOption(str, Enum):
    """
    Feature's options.
    """
    N_RUNS = 'n_runs'
    DISTRIBUTION = 'distribution'
    PERTURBATION_LEVEL = 'perturbation_level'
    NOISE_LEVEL = 'noise_level'
    NOISE_TYPE = 'noise_type'
    SIGNIFICANT_DIGITS = 'significant_digits'
    PERTURBED_TRAILING_DIGITS = 'perturbed_trailing_digits'
    ROTATED = 'rotated'
    CONDITION_FACTOR = 'condition_factor'
    NAN_RATE = 'nan_rate'
    UNRELAXABLE_BOUNDS = 'unrelaxable_bounds'
    UNRELAXABLE_LINEAR_CONSTRAINTS = 'unrelaxable_linear_constraints'
    UNRELAXABLE_NONLINEAR_CONSTRAINTS = 'unrelaxable_nonlinear_constraints'
    MESH_SIZE = 'mesh_size'
    MESH_TYPE = 'mesh_type'
    GROUND_TRUTH = 'ground_truth'
    MOD_X0 = 'mod_x0'
    MOD_AFFINE = 'mod_affine'
    MOD_BOUNDS = 'mod_bounds'
    MOD_LINEAR_UB = 'mod_linear_ub'
    MOD_LINEAR_EQ = 'mod_linear_eq'
    MOD_FUN = 'mod_fun'
    MOD_CUB = 'mod_cub'
    MOD_CEQ = 'mod_ceq'


class ProblemError(Exception):
    """
    Exception raised when a problem cannot be loaded.
    """


def show_versions():
    """
    Display useful system and dependency information.

    When reporting issues, please include this information.
    """
    # Get system and dependency information.
    sys_info = _get_sys_info()
    deps_info = _get_deps_info()

    # Display information.
    print('System settings')
    print('---------------')
    print('\n'.join(
        f'{k:>{max(map(len, sys_info.keys())) + 1}}: {v}'
        for k, v in sys_info.items()
    ))
    print()
    print('Python dependencies')
    print('-------------------')
    print('\n'.join(
        f'{k:>{max(map(len, deps_info.keys())) + 1}}: {v}'
        for k, v in deps_info.items()
    ))


def get_logger(name, level=logging.INFO):
    """
    Get a logger.

    Parameters
    ----------
    name : str
        Name of the logger. If ``None``, the root logger is returned.
    level : int, optional
        Logging level.

    Returns
    -------
    `logging.Logger`
        Logger with the given name. If a logger with the given name has already
        been created, it is returned instead of creating a new one.
    """
    logger = logging.getLogger(name)
    if len(logger.handlers) == 0:
        logger.setLevel(level)
        formatter = logging.Formatter('[%(levelname)-8s] %(message)s')

        # Attach a console handler.
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    return logger


def _get_sys_info():
    """
    Get system information.

    Returns
    -------
    dict
        System information.
    """
    return {
        'python': sys.version.replace(os.linesep, ' '),
        'executable': sys.executable,
        'machine': platform.platform(),
    }


def _get_deps_info():
    """
    Get dependency information.

    Returns
    -------
    dict
        Dependency information.
    """
    deps = [
        'optiprofiler',
        'matplotlib',
        'numpy',
        'pycutest',
        'pypdf',
        'scipy',
        'setuptools',
        'pip',
    ]
    deps_info = {}
    for module in deps:
        try:
            deps_info[module] = version(module)
        except PackageNotFoundError:
            deps_info[module] = None
    return deps_info
