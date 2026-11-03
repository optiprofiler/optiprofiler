import hashlib
import os
import re
import textwrap
from collections.abc import Mapping
from pathlib import Path
import numpy as np
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D
from pypdf import PdfWriter


from .problem_libraries import list_problem_libraries, resolve_problem_library
from .experiment import validate_n_runs
from . import scoring
# Option pickles written by 1.x and by 2.0 builds before scoring.py reference
# the default callbacks under this module.
from .scoring import _default_merit, _default_score_fun, _default_score_weight_fun  # noqa: F401
from .provenance import full_feature_stamp
from .metadata import credential_free_urls
from .utils import ProfileOption, ProblemOption, get_logger, print_log_message, shorten_log_message


def _is_positive_infinity(value):
    """Return whether *value* is a scalar positive infinity."""
    return (
        isinstance(value, (float, np.floating))
        and np.isinf(value)
        and value > 0
    )


def _get_conservative_default_n_jobs():
    """
    Get a conservative default number of parallel jobs.
    """
    n_workers = os.cpu_count() or 1
    if n_workers <= 1:
        return 1
    return max(2, n_workers // 2)


def check_validity_problem_options(problem_options, require_available=True):
    """
    Check the validity of the problem options.

    Parameters
    ----------
    problem_options : dict
        Problem options to normalize and validate.
    require_available : bool, optional
        Resolve each requested problem-library provider when true. Loading a
        saved experiment sets this to false because filtering stored results
        by library name does not require that provider to remain installed.
    """
    # Check custom_problem_libs_path first, as it affects plibs validation.
    custom_libs_path = None
    if ProblemOption.CUSTOM_PROBLEM_LIBS_PATH in problem_options:
        custom_libs_path = problem_options[ProblemOption.CUSTOM_PROBLEM_LIBS_PATH]
        if custom_libs_path is not None:
            if not isinstance(custom_libs_path, (str, Path)):
                raise TypeError(f'Option {ProblemOption.CUSTOM_PROBLEM_LIBS_PATH} must be a string or a Path.')
            # Expand ~ to user's home directory before resolving.
            custom_libs_path = Path(custom_libs_path).expanduser().resolve()
            if not custom_libs_path.exists():
                raise ValueError(f'Option {ProblemOption.CUSTOM_PROBLEM_LIBS_PATH} path does not exist: {custom_libs_path}')
            if not custom_libs_path.is_dir():
                raise ValueError(f'Option {ProblemOption.CUSTOM_PROBLEM_LIBS_PATH} must be a directory: {custom_libs_path}')
            # Store the resolved path back.
            problem_options[ProblemOption.CUSTOM_PROBLEM_LIBS_PATH] = custom_libs_path

    if ProblemOption.PLIBS in problem_options:
        if not isinstance(problem_options[ProblemOption.PLIBS], list):
            problem_options[ProblemOption.PLIBS] = [problem_options[ProblemOption.PLIBS]]
        # Check if plibs is empty or any element is not str.
        if len(problem_options[ProblemOption.PLIBS]) == 0:
            raise ValueError(f'Option {ProblemOption.PLIBS} cannot be an empty list.')
        elif any(not isinstance(plib, str) for plib in problem_options[ProblemOption.PLIBS]):
            raise TypeError(f'Option {ProblemOption.PLIBS} must be a string or a list of strings.')
        if require_available:
            all_available_libs = list_problem_libraries(custom_libs_path)
            for plib in problem_options[ProblemOption.PLIBS]:
                try:
                    resolve_problem_library(plib, custom_libs_path)
                except ValueError as exc:
                    raise ValueError(
                        f'Option {ProblemOption.PLIBS} contains invalid problem '
                        f'libraries. Available libraries: {all_available_libs}. '
                        f'{exc}'
                    ) from exc

    if ProblemOption.PLIB_OPTIONS in problem_options:
        plib_options = problem_options[ProblemOption.PLIB_OPTIONS]
        if plib_options is None:
            plib_options = {}
        if not isinstance(plib_options, Mapping):
            raise TypeError(
                f'Option {ProblemOption.PLIB_OPTIONS} must be a mapping from '
                'problem-library names to option mappings.'
            )
        plib_options = dict(plib_options)
        if any(not isinstance(plib, str) for plib in plib_options):
            raise TypeError(
                f'Option {ProblemOption.PLIB_OPTIONS} must use problem-library '
                'name strings as keys.'
            )
        for plib, options in plib_options.items():
            if not isinstance(options, Mapping):
                raise TypeError(
                    f'Options for problem library "{plib}" in '
                    f'{ProblemOption.PLIB_OPTIONS} must be a mapping.'
                )
            options = dict(options)
            if any(not isinstance(key, str) for key in options):
                raise TypeError(
                    f'Options for problem library "{plib}" in '
                    f'{ProblemOption.PLIB_OPTIONS} must use string keys.'
                )
            plib_options[plib] = options
        problem_options[ProblemOption.PLIB_OPTIONS] = plib_options

    if ProblemOption.PTYPE in problem_options:
        if not isinstance(problem_options[ProblemOption.PTYPE], str):
            raise TypeError(f'Option {ProblemOption.PTYPE} must be a string.')
        if len(problem_options[ProblemOption.PTYPE]) == 0:
            raise ValueError(f'Option {ProblemOption.PTYPE} cannot be an empty string.')
        if any(t not in 'ubln' for t in problem_options[ProblemOption.PTYPE]):
            raise ValueError(f"Option {ProblemOption.PTYPE} should be a string containing only the characters 'u', 'b', 'l', and 'n'.")

    if ProblemOption.MINDIM in problem_options:
        if isinstance(problem_options[ProblemOption.MINDIM], (float, np.floating)) and float(problem_options[ProblemOption.MINDIM]).is_integer():
            problem_options[ProblemOption.MINDIM] = int(problem_options[ProblemOption.MINDIM])
        if isinstance(problem_options[ProblemOption.MINDIM], np.integer):
            problem_options[ProblemOption.MINDIM] = int(problem_options[ProblemOption.MINDIM])
        if not isinstance(problem_options[ProblemOption.MINDIM], int):
            raise TypeError(f'Option {ProblemOption.MINDIM} must be an integer.')
        if problem_options[ProblemOption.MINDIM] < 1:
            raise ValueError(f'Option {ProblemOption.MINDIM} must be at least 1.')
    if ProblemOption.MAXDIM in problem_options:
        if not _is_positive_infinity(problem_options[ProblemOption.MAXDIM]):
            if isinstance(problem_options[ProblemOption.MAXDIM], (float, np.floating)) and float(problem_options[ProblemOption.MAXDIM]).is_integer():
                problem_options[ProblemOption.MAXDIM] = int(problem_options[ProblemOption.MAXDIM])
            if isinstance(problem_options[ProblemOption.MAXDIM], np.integer):
                problem_options[ProblemOption.MAXDIM] = int(problem_options[ProblemOption.MAXDIM])
            if not isinstance(problem_options[ProblemOption.MAXDIM], int):
                raise TypeError(f'Option {ProblemOption.MAXDIM} must be an integer or positive np.inf.')
            if problem_options[ProblemOption.MAXDIM] < 1:
                raise ValueError(f'Option {ProblemOption.MAXDIM} must be at least 1 or positive np.inf.')
    if ProblemOption.MINDIM in problem_options and ProblemOption.MAXDIM in problem_options:
        if problem_options[ProblemOption.MINDIM] > problem_options[ProblemOption.MAXDIM]:
            raise ValueError(f'Option {ProblemOption.MINDIM} cannot be larger than option {ProblemOption.MAXDIM}.')

    if ProblemOption.MINB in problem_options:
        if isinstance(problem_options[ProblemOption.MINB], (float, np.floating)) and float(problem_options[ProblemOption.MINB]).is_integer():
            problem_options[ProblemOption.MINB] = int(problem_options[ProblemOption.MINB])
        if isinstance(problem_options[ProblemOption.MINB], np.integer):
            problem_options[ProblemOption.MINB] = int(problem_options[ProblemOption.MINB])
        if not isinstance(problem_options[ProblemOption.MINB], int):
            raise TypeError(f'Option {ProblemOption.MINB} must be an integer.')
        if problem_options[ProblemOption.MINB] < 0:
            raise ValueError(f'Option {ProblemOption.MINB} must be nonnegative.')
    if ProblemOption.MAXB in problem_options:
        if not _is_positive_infinity(problem_options[ProblemOption.MAXB]):
            if isinstance(problem_options[ProblemOption.MAXB], (float, np.floating)) and float(problem_options[ProblemOption.MAXB]).is_integer():
                problem_options[ProblemOption.MAXB] = int(problem_options[ProblemOption.MAXB])
            if isinstance(problem_options[ProblemOption.MAXB], np.integer):
                problem_options[ProblemOption.MAXB] = int(problem_options[ProblemOption.MAXB])
            if not isinstance(problem_options[ProblemOption.MAXB], int):
                raise TypeError(f'Option {ProblemOption.MAXB} must be an integer or positive np.inf.')
            if problem_options[ProblemOption.MAXB] < 0:
                raise ValueError(f'Option {ProblemOption.MAXB} must be nonnegative or positive np.inf.')
    if ProblemOption.MINB in problem_options and ProblemOption.MAXB in problem_options:
        if problem_options[ProblemOption.MINB] > problem_options[ProblemOption.MAXB]:
            raise ValueError(f'Option {ProblemOption.MINB} cannot be larger than option {ProblemOption.MAXB}.')

    if ProblemOption.MINLCON in problem_options:
        if isinstance(problem_options[ProblemOption.MINLCON], (float, np.floating)) and float(problem_options[ProblemOption.MINLCON]).is_integer():
            problem_options[ProblemOption.MINLCON] = int(problem_options[ProblemOption.MINLCON])
        if isinstance(problem_options[ProblemOption.MINLCON], np.integer):
            problem_options[ProblemOption.MINLCON] = int(problem_options[ProblemOption.MINLCON])
        if not isinstance(problem_options[ProblemOption.MINLCON], int):
            raise TypeError(f'Option {ProblemOption.MINLCON} must be an integer.')
        if problem_options[ProblemOption.MINLCON] < 0:
            raise ValueError(f'Option {ProblemOption.MINLCON} must be nonnegative.')
    if ProblemOption.MAXLCON in problem_options:
        if not _is_positive_infinity(problem_options[ProblemOption.MAXLCON]):
            if isinstance(problem_options[ProblemOption.MAXLCON], (float, np.floating)) and float(problem_options[ProblemOption.MAXLCON]).is_integer():
                problem_options[ProblemOption.MAXLCON] = int(problem_options[ProblemOption.MAXLCON])
            if isinstance(problem_options[ProblemOption.MAXLCON], np.integer):
                problem_options[ProblemOption.MAXLCON] = int(problem_options[ProblemOption.MAXLCON])
            if not isinstance(problem_options[ProblemOption.MAXLCON], int):
                raise TypeError(f'Option {ProblemOption.MAXLCON} must be an integer or positive np.inf.')
            if problem_options[ProblemOption.MAXLCON] < 0:
                raise ValueError(f'Option {ProblemOption.MAXLCON} must be nonnegative or positive np.inf.')
    if ProblemOption.MINLCON in problem_options and ProblemOption.MAXLCON in problem_options:
        if problem_options[ProblemOption.MINLCON] > problem_options[ProblemOption.MAXLCON]:
            raise ValueError(f'Option {ProblemOption.MINLCON} cannot be larger than option {ProblemOption.MAXLCON}.')

    if ProblemOption.MINNLCON in problem_options:
        if isinstance(problem_options[ProblemOption.MINNLCON], (float, np.floating)) and float(problem_options[ProblemOption.MINNLCON]).is_integer():
            problem_options[ProblemOption.MINNLCON] = int(problem_options[ProblemOption.MINNLCON])
        if isinstance(problem_options[ProblemOption.MINNLCON], np.integer):
            problem_options[ProblemOption.MINNLCON] = int(problem_options[ProblemOption.MINNLCON])
        if not isinstance(problem_options[ProblemOption.MINNLCON], int):
            raise TypeError(f'Option {ProblemOption.MINNLCON} must be an integer.')
        if problem_options[ProblemOption.MINNLCON] < 0:
            raise ValueError(f'Option {ProblemOption.MINNLCON} must be nonnegative.')
    if ProblemOption.MAXNLCON in problem_options:
        if not _is_positive_infinity(problem_options[ProblemOption.MAXNLCON]):
            if isinstance(problem_options[ProblemOption.MAXNLCON], (float, np.floating)) and float(problem_options[ProblemOption.MAXNLCON]).is_integer():
                problem_options[ProblemOption.MAXNLCON] = int(problem_options[ProblemOption.MAXNLCON])
            if isinstance(problem_options[ProblemOption.MAXNLCON], np.integer):
                problem_options[ProblemOption.MAXNLCON] = int(problem_options[ProblemOption.MAXNLCON])
            if not isinstance(problem_options[ProblemOption.MAXNLCON], int):
                raise TypeError(f'Option {ProblemOption.MAXNLCON} must be an integer or positive np.inf.')
            if problem_options[ProblemOption.MAXNLCON] < 0:
                raise ValueError(f'Option {ProblemOption.MAXNLCON} must be nonnegative or positive np.inf.')
    if ProblemOption.MINNLCON in problem_options and ProblemOption.MAXNLCON in problem_options:
        if problem_options[ProblemOption.MINNLCON] > problem_options[ProblemOption.MAXNLCON]:
            raise ValueError(f'Option {ProblemOption.MINNLCON} cannot be larger than option {ProblemOption.MAXNLCON}.')

    if ProblemOption.MINCON in problem_options:
        if isinstance(problem_options[ProblemOption.MINCON], (float, np.floating)) and float(problem_options[ProblemOption.MINCON]).is_integer():
            problem_options[ProblemOption.MINCON] = int(problem_options[ProblemOption.MINCON])
        if isinstance(problem_options[ProblemOption.MINCON], np.integer):
            problem_options[ProblemOption.MINCON] = int(problem_options[ProblemOption.MINCON])
        if not isinstance(problem_options[ProblemOption.MINCON], int):
            raise TypeError(f'Option {ProblemOption.MINCON} must be an integer.')
        if problem_options[ProblemOption.MINCON] < 0:
            raise ValueError(f'Option {ProblemOption.MINCON} must be nonnegative.')
    if ProblemOption.MAXCON in problem_options:
        if not _is_positive_infinity(problem_options[ProblemOption.MAXCON]):
            if isinstance(problem_options[ProblemOption.MAXCON], (float, np.floating)) and float(problem_options[ProblemOption.MAXCON]).is_integer():
                problem_options[ProblemOption.MAXCON] = int(problem_options[ProblemOption.MAXCON])
            if isinstance(problem_options[ProblemOption.MAXCON], np.integer):
                problem_options[ProblemOption.MAXCON] = int(problem_options[ProblemOption.MAXCON])
            if not isinstance(problem_options[ProblemOption.MAXCON], int):
                raise TypeError(f'Option {ProblemOption.MAXCON} must be an integer or positive np.inf.')
            if problem_options[ProblemOption.MAXCON] < 0:
                raise ValueError(f'Option {ProblemOption.MAXCON} must be nonnegative or positive np.inf.')
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


def check_validity_profile_options(solvers, profile_options):
    """
    Check the validity of the profile options.
    """
    if ProfileOption.PLOT_OPTIONS in profile_options:
        from .plot_options import normalize_plot_options
        profile_options[ProfileOption.PLOT_OPTIONS] = normalize_plot_options(profile_options[ProfileOption.PLOT_OPTIONS])
    if ProfileOption.N_RUNS in profile_options:
        profile_options[ProfileOption.N_RUNS] = validate_n_runs(profile_options[ProfileOption.N_RUNS])
    if ProfileOption.N_JOBS in profile_options:
        if isinstance(profile_options[ProfileOption.N_JOBS], (float, np.floating)) and float(profile_options[ProfileOption.N_JOBS]).is_integer():
            profile_options[ProfileOption.N_JOBS] = int(profile_options[ProfileOption.N_JOBS])
        if isinstance(profile_options[ProfileOption.N_JOBS], np.integer):
            profile_options[ProfileOption.N_JOBS] = int(profile_options[ProfileOption.N_JOBS])
        if not isinstance(profile_options[ProfileOption.N_JOBS], int):
            raise TypeError(f'Option {ProfileOption.N_JOBS} must be an integer.')
        if profile_options[ProfileOption.N_JOBS] < 1:
            profile_options[ProfileOption.N_JOBS] = 1
            print_log_message('WARNING', f'Option {ProfileOption.N_JOBS} is set to 1 because it cannot be smaller than 1.')

    if ProfileOption.SEED in profile_options:
        if isinstance(profile_options[ProfileOption.SEED], (float, np.floating)) and float(profile_options[ProfileOption.SEED]).is_integer():
            profile_options[ProfileOption.SEED] = int(profile_options[ProfileOption.SEED])
        if isinstance(profile_options[ProfileOption.SEED], np.integer):
            profile_options[ProfileOption.SEED] = int(profile_options[ProfileOption.SEED])
        if not isinstance(profile_options[ProfileOption.SEED], int):
            raise TypeError(f'Option {ProfileOption.SEED} must be an integer.')
        if profile_options[ProfileOption.SEED] < 0:
            raise ValueError(f'Option {ProfileOption.SEED} must be nonnegative.')

    if ProfileOption.BENCHMARK_ID in profile_options:
        if not isinstance(profile_options[ProfileOption.BENCHMARK_ID], str):
            raise TypeError(f'Option {ProfileOption.BENCHMARK_ID} must be a string.')
        if len(profile_options[ProfileOption.BENCHMARK_ID]) == 0:
            raise ValueError(f'Option {ProfileOption.BENCHMARK_ID} cannot be an empty string.')
        if not re.match(r'^[a-zA-Z0-9_.-]+$', profile_options[ProfileOption.BENCHMARK_ID]):
            raise ValueError(f"Option {ProfileOption.BENCHMARK_ID} should be a string satisfying the strict file name requirements (only containing letters, numbers, underscores, hyphens, and dots).")

    if ProfileOption.SOLVER_NAMES in profile_options:
        if not isinstance(profile_options[ProfileOption.SOLVER_NAMES], list):
            raise TypeError(f'Option {ProfileOption.SOLVER_NAMES} must be a list of strings.')
        # When loading saved results, `solvers` is None and the effective solver
        # count is only known after the load; the length check is then deferred
        # to `check_post_load_profile_options`.
        if solvers is not None and len(profile_options[ProfileOption.SOLVER_NAMES]) != len(solvers):
            raise ValueError(f'Option {ProfileOption.SOLVER_NAMES} must have the same length as the number of solvers ({len(solvers)}).')
        if any(not isinstance(name, str) for name in profile_options[ProfileOption.SOLVER_NAMES]):
            raise TypeError(f'Option {ProfileOption.SOLVER_NAMES} must be a list of strings.')

    if ProfileOption.SOLVER_ISRAND in profile_options:
        if not isinstance(profile_options[ProfileOption.SOLVER_ISRAND], list):
            raise TypeError(f'Option {ProfileOption.SOLVER_ISRAND} must be a list of booleans.')
        if solvers is not None and len(profile_options[ProfileOption.SOLVER_ISRAND]) != len(solvers):
            raise ValueError(f'Option {ProfileOption.SOLVER_ISRAND} must have the same length as the number of solvers ({len(solvers)}).')
        if any(not isinstance(isrand, bool) for isrand in profile_options[ProfileOption.SOLVER_ISRAND]):
            raise TypeError(f'Option {ProfileOption.SOLVER_ISRAND} must be a list of booleans.')

    if ProfileOption.FEATURE_STAMP in profile_options:
        if not isinstance(profile_options[ProfileOption.FEATURE_STAMP], str):
            raise TypeError(f'Option {ProfileOption.FEATURE_STAMP} must be a string.')
        if len(profile_options[ProfileOption.FEATURE_STAMP]) == 0:
            raise ValueError(f'Option {ProfileOption.FEATURE_STAMP} cannot be an empty string.')
        if not re.match(r'^[a-zA-Z0-9_.-]+$', profile_options[ProfileOption.FEATURE_STAMP]):
            raise ValueError(f"Option {ProfileOption.FEATURE_STAMP} should be a string satisfying the strict file name requirements (only containing letters, numbers, underscores, hyphens, and dots).")

    if ProfileOption.ERRORBAR_TYPE in profile_options:
        if not isinstance(profile_options[ProfileOption.ERRORBAR_TYPE], str):
            raise TypeError(f'Option {ProfileOption.ERRORBAR_TYPE} must be a string.')
        if profile_options[ProfileOption.ERRORBAR_TYPE] not in ['minmax', 'meanstd']:
            raise ValueError(f"Option {ProfileOption.ERRORBAR_TYPE} must be 'minmax' or 'meanstd'.")
        
    if ProfileOption.HIST_AGGREGATION in profile_options:
        if not isinstance(profile_options[ProfileOption.HIST_AGGREGATION], str):
            raise TypeError(f'Option {ProfileOption.HIST_AGGREGATION} must be a string.')
        if profile_options[ProfileOption.HIST_AGGREGATION] not in ['min', 'mean', 'max']:
            raise ValueError(f"Option {ProfileOption.HIST_AGGREGATION} must be 'min', 'mean', or 'max'.")

    if ProfileOption.SAVEPATH in profile_options:
        if not isinstance(profile_options[ProfileOption.SAVEPATH], str):
            raise TypeError(f'Option {ProfileOption.SAVEPATH} must be a string.')
        if len(profile_options[ProfileOption.SAVEPATH]) == 0:
            raise ValueError(f'Option {ProfileOption.SAVEPATH} cannot be an empty string.')
        expanded_path = Path(profile_options[ProfileOption.SAVEPATH]).expanduser().resolve()
        profile_options[ProfileOption.SAVEPATH] = str(expanded_path)
        if not expanded_path.exists():
            try:
                expanded_path.mkdir(parents=True, exist_ok=True)
            except Exception as exc:
                raise ValueError(f'Option {ProfileOption.SAVEPATH} is not a valid path and cannot be created: {exc}.') from exc

    if ProfileOption.MAX_TOL_ORDER in profile_options:
        if isinstance(profile_options[ProfileOption.MAX_TOL_ORDER], (float, np.floating)) and float(profile_options[ProfileOption.MAX_TOL_ORDER]).is_integer():
            profile_options[ProfileOption.MAX_TOL_ORDER] = int(profile_options[ProfileOption.MAX_TOL_ORDER])
        if isinstance(profile_options[ProfileOption.MAX_TOL_ORDER], np.integer):
            profile_options[ProfileOption.MAX_TOL_ORDER] = int(profile_options[ProfileOption.MAX_TOL_ORDER])
        if not isinstance(profile_options[ProfileOption.MAX_TOL_ORDER], int):
            raise TypeError(f'Option {ProfileOption.MAX_TOL_ORDER.value} must be an integer.')
        if profile_options[ProfileOption.MAX_TOL_ORDER] <= 0 or profile_options[ProfileOption.MAX_TOL_ORDER] > 16:
            raise ValueError(f'Option {ProfileOption.MAX_TOL_ORDER.value} must be a positive integer not larger than 16.')

    if ProfileOption.MAX_EVAL_FACTOR in profile_options:
        if isinstance(profile_options[ProfileOption.MAX_EVAL_FACTOR], bool) or not isinstance(profile_options[ProfileOption.MAX_EVAL_FACTOR], (float, int)):
            raise TypeError(f'Option {ProfileOption.MAX_EVAL_FACTOR.value} must be a float or an integer.')
        if not np.isfinite(profile_options[ProfileOption.MAX_EVAL_FACTOR]) or profile_options[ProfileOption.MAX_EVAL_FACTOR] <= 0:
            raise ValueError(f'Option {ProfileOption.MAX_EVAL_FACTOR.value} must be finite and positive.')

    if ProfileOption.MERIT_FUN in profile_options:
        if not callable(profile_options[ProfileOption.MERIT_FUN]):
            raise TypeError(f'Option {ProfileOption.MERIT_FUN} must be a callable function.')

    if ProfileOption.PROJECT_X0 in profile_options:
        if not isinstance(profile_options[ProfileOption.PROJECT_X0], bool):
            raise TypeError(f'Option {ProfileOption.PROJECT_X0} must be a boolean.')

    if ProfileOption.RUN_PLAIN in profile_options:
        if not isinstance(profile_options[ProfileOption.RUN_PLAIN], bool):
            raise TypeError(f'Option {ProfileOption.RUN_PLAIN} must be a boolean.')

    if ProfileOption.SCORE_ONLY in profile_options:
        if not isinstance(profile_options[ProfileOption.SCORE_ONLY], bool):
            raise TypeError(f'Option {ProfileOption.SCORE_ONLY} must be a boolean.')

    if ProfileOption.DRAW_HIST_PLOTS in profile_options:
        if not isinstance(profile_options[ProfileOption.DRAW_HIST_PLOTS], str):
            raise TypeError(f'Option {ProfileOption.DRAW_HIST_PLOTS} must be a string.')
        valid_modes = ['none', 'parallel', 'sequential']
        if profile_options[ProfileOption.DRAW_HIST_PLOTS] not in valid_modes:
            raise ValueError(f"Option {ProfileOption.DRAW_HIST_PLOTS} must be 'none', 'parallel', or 'sequential'.")

    if ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES in profile_options:
        if not isinstance(profile_options[ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES], bool):
            raise TypeError(f'Option {ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES} must be a boolean.')

    if ProfileOption.SUMMARIZE_DATA_PROFILES in profile_options:
        if not isinstance(profile_options[ProfileOption.SUMMARIZE_DATA_PROFILES], bool):
            raise TypeError(f'Option {ProfileOption.SUMMARIZE_DATA_PROFILES} must be a boolean.')

    if ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES in profile_options:
        if not isinstance(profile_options[ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES], bool):
            raise TypeError(f'Option {ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES} must be a boolean.')
        # When loading saved results, `solvers` is None; the two-solver check is
        # deferred to `check_post_load_profile_options`, which knows the final
        # effective solver count after `solvers_to_load` is applied.
        if profile_options[ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES] and solvers is not None:
            if len(solvers) != 2:
                print_log_message('WARNING', 'The log-ratio profiles are available only when there are exactly two solvers. We will not generate the log-ratio profiles.')
                profile_options[ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES] = False

    if ProfileOption.SUMMARIZE_OUTPUT_BASED_PROFILES in profile_options:
        if not isinstance(profile_options[ProfileOption.SUMMARIZE_OUTPUT_BASED_PROFILES], bool):
            raise TypeError(f'Option {ProfileOption.SUMMARIZE_OUTPUT_BASED_PROFILES} must be a boolean.')

    if ProfileOption.SILENT in profile_options:
        if not isinstance(profile_options[ProfileOption.SILENT], bool):
            raise TypeError(f'Option {ProfileOption.SILENT} must be a boolean.')

    if ProfileOption.SOLVER_VERBOSE in profile_options:
        if isinstance(profile_options[ProfileOption.SOLVER_VERBOSE], (float, np.floating)) and float(profile_options[ProfileOption.SOLVER_VERBOSE]).is_integer():
            profile_options[ProfileOption.SOLVER_VERBOSE] = int(profile_options[ProfileOption.SOLVER_VERBOSE])
        if isinstance(profile_options[ProfileOption.SOLVER_VERBOSE], np.integer):
            profile_options[ProfileOption.SOLVER_VERBOSE] = int(profile_options[ProfileOption.SOLVER_VERBOSE])
        if not isinstance(profile_options[ProfileOption.SOLVER_VERBOSE], int):
            raise TypeError(f'Option {ProfileOption.SOLVER_VERBOSE} must be 0, 1, or 2.')
        if profile_options[ProfileOption.SOLVER_VERBOSE] < 0 or profile_options[ProfileOption.SOLVER_VERBOSE] > 2:
            raise ValueError(f'Option {ProfileOption.SOLVER_VERBOSE} must be 0, 1, or 2.')
        if (ProfileOption.SILENT in profile_options and profile_options[ProfileOption.SILENT]
                and profile_options[ProfileOption.SOLVER_VERBOSE] != 0):
            print_log_message('WARNING', f'The option `{ProfileOption.SOLVER_VERBOSE}` will be set to 0 because {ProfileOption.SILENT} is True.')
            profile_options[ProfileOption.SOLVER_VERBOSE] = 0

    if ProfileOption.SEMILOGX in profile_options:
        if not isinstance(profile_options[ProfileOption.SEMILOGX], bool):
            raise TypeError(f'Option {ProfileOption.SEMILOGX} must be a boolean.')

    if ProfileOption.NORMALIZED_SCORES in profile_options:
        if not isinstance(profile_options[ProfileOption.NORMALIZED_SCORES], bool):
            raise TypeError(f'Option {ProfileOption.NORMALIZED_SCORES} must be a boolean.')

    if ProfileOption.SCORE_WEIGHT_FUN in profile_options:
        if not callable(profile_options[ProfileOption.SCORE_WEIGHT_FUN]):
            raise TypeError(f'Option {ProfileOption.SCORE_WEIGHT_FUN} must be a callable function.')

    if ProfileOption.SCORE_FUN in profile_options:
        if not callable(profile_options[ProfileOption.SCORE_FUN]):
            raise TypeError(f'Option {ProfileOption.SCORE_FUN} must be a callable function.')

    if ProfileOption.LOAD in profile_options:
        if not isinstance(profile_options[ProfileOption.LOAD], str) and profile_options[ProfileOption.LOAD] is not None:
            raise TypeError(f'Option {ProfileOption.LOAD} must be a string or None.')
        if isinstance(profile_options[ProfileOption.LOAD], str):
            if len(profile_options[ProfileOption.LOAD]) == 0:
                profile_options[ProfileOption.LOAD] = None

    if ProfileOption.SOLVERS_TO_LOAD in profile_options:
        if ProfileOption.LOAD not in profile_options or ProfileOption.LOAD is None:
            raise ValueError(f'Option {ProfileOption.SOLVERS_TO_LOAD} cannot be used if option {ProfileOption.LOAD} is not set or is None.')
        if not isinstance(profile_options[ProfileOption.SOLVERS_TO_LOAD], list):
            raise TypeError(f'Option {ProfileOption.SOLVERS_TO_LOAD} must be a list of integers.')
        if len(profile_options[ProfileOption.SOLVERS_TO_LOAD]) < 2:
            raise ValueError(f'Option {ProfileOption.SOLVERS_TO_LOAD} must contain at least two indices since at least two solvers are needed for a comparison.')
        if ProfileOption.SOLVER_NAMES in profile_options:
            if len(profile_options[ProfileOption.SOLVERS_TO_LOAD]) != len(profile_options[ProfileOption.SOLVER_NAMES]):
                raise ValueError(f'Option {ProfileOption.SOLVERS_TO_LOAD} must have the same length as option {ProfileOption.SOLVER_NAMES} ({len(profile_options[ProfileOption.SOLVER_NAMES])}).')
        for i in range(len(profile_options[ProfileOption.SOLVERS_TO_LOAD])):
            if isinstance(profile_options[ProfileOption.SOLVERS_TO_LOAD][i], (float, np.floating)) and float(profile_options[ProfileOption.SOLVERS_TO_LOAD][i]).is_integer():
                profile_options[ProfileOption.SOLVERS_TO_LOAD][i] = int(profile_options[ProfileOption.SOLVERS_TO_LOAD][i])
            if isinstance(profile_options[ProfileOption.SOLVERS_TO_LOAD][i], np.integer):
                profile_options[ProfileOption.SOLVERS_TO_LOAD][i] = int(profile_options[ProfileOption.SOLVERS_TO_LOAD][i])
            if not isinstance(profile_options[ProfileOption.SOLVERS_TO_LOAD][i], int):
                raise TypeError(f'All elements of option {ProfileOption.SOLVERS_TO_LOAD} must be integers.')
            if profile_options[ProfileOption.SOLVERS_TO_LOAD][i] < 0:
                raise ValueError(f'All elements of option {ProfileOption.SOLVERS_TO_LOAD} must be nonnegative.')
            if profile_options[ProfileOption.SOLVERS_TO_LOAD].count(profile_options[ProfileOption.SOLVERS_TO_LOAD][i]) > 1:
                raise ValueError(f'All elements of option {ProfileOption.SOLVERS_TO_LOAD} must be different.')

    if ProfileOption.LINE_COLORS in profile_options:
        if not isinstance(profile_options[ProfileOption.LINE_COLORS], list):
            profile_options[ProfileOption.LINE_COLORS] = [profile_options[ProfileOption.LINE_COLORS]]
        if len(profile_options[ProfileOption.LINE_COLORS]) == 0:
            raise ValueError(f'Option {ProfileOption.LINE_COLORS} cannot be an empty list.')
        for i in range(len(profile_options[ProfileOption.LINE_COLORS])):
            if not is_color_like(profile_options[ProfileOption.LINE_COLORS][i]):
                    raise ValueError(f'All string elements of option {ProfileOption.LINE_COLORS} must be valid matplotlib colors.')

    if ProfileOption.LINE_STYLES in profile_options:
        if not isinstance(profile_options[ProfileOption.LINE_STYLES], list):
            profile_options[ProfileOption.LINE_STYLES] = [profile_options[ProfileOption.LINE_STYLES]]
        if len(profile_options[ProfileOption.LINE_STYLES]) == 0:
            raise ValueError(f'Option {ProfileOption.LINE_STYLES} cannot be an empty list.')
        def is_valid_linestyle(style):
            try:
                # Validate the rendered property without creating a GUI figure
                # or drawing into (and closing) the caller's current axes.
                Line2D([], [], linestyle=style)
                return True
            except Exception:
                return False
        for i in range(len(profile_options[ProfileOption.LINE_STYLES])):
            if not is_valid_linestyle(profile_options[ProfileOption.LINE_STYLES][i]):
                raise ValueError(f'All elements of option {ProfileOption.LINE_STYLES} must be valid matplotlib line styles.')

    if ProfileOption.LINE_WIDTHS in profile_options:
        if not isinstance(profile_options[ProfileOption.LINE_WIDTHS], list):
            profile_options[ProfileOption.LINE_WIDTHS] = [profile_options[ProfileOption.LINE_WIDTHS]]
        if len(profile_options[ProfileOption.LINE_WIDTHS]) == 0:
            raise ValueError(f'Option {ProfileOption.LINE_WIDTHS} cannot be an empty list.')
        for i in range(len(profile_options[ProfileOption.LINE_WIDTHS])):
            if not isinstance(profile_options[ProfileOption.LINE_WIDTHS][i], (float, int)):
                raise TypeError(f'All elements of option {ProfileOption.LINE_WIDTHS} must be floats or integers.')
            if profile_options[ProfileOption.LINE_WIDTHS][i] <= 0:
                raise ValueError(f'All elements of option {ProfileOption.LINE_WIDTHS} must be positive.')

    if ProfileOption.BAR_COLORS in profile_options:
        if not isinstance(profile_options[ProfileOption.BAR_COLORS], list):
            profile_options[ProfileOption.BAR_COLORS] = [profile_options[ProfileOption.BAR_COLORS]]
        if len(profile_options[ProfileOption.BAR_COLORS]) == 0:
            raise ValueError(f'Option {ProfileOption.BAR_COLORS} cannot be an empty list.')
        for i in range(len(profile_options[ProfileOption.BAR_COLORS])):
            if not is_color_like(profile_options[ProfileOption.BAR_COLORS][i]):
                    raise ValueError(f'All string elements of option {ProfileOption.BAR_COLORS} must be valid matplotlib colors.')

    if ProfileOption.XLABEL_PERFORMANCE_PROFILE in profile_options:
        if not isinstance(profile_options[ProfileOption.XLABEL_PERFORMANCE_PROFILE], str):
            raise TypeError(f'Option {ProfileOption.XLABEL_PERFORMANCE_PROFILE} must be a string.')

    if ProfileOption.YLABEL_PERFORMANCE_PROFILE in profile_options:
        if not isinstance(profile_options[ProfileOption.YLABEL_PERFORMANCE_PROFILE], str):
            raise TypeError(f'Option {ProfileOption.YLABEL_PERFORMANCE_PROFILE} must be a string.')

    if ProfileOption.XLABEL_DATA_PROFILE in profile_options:
        if not isinstance(profile_options[ProfileOption.XLABEL_DATA_PROFILE], str):
            raise TypeError(f'Option {ProfileOption.XLABEL_DATA_PROFILE} must be a string.')

    if ProfileOption.YLABEL_DATA_PROFILE in profile_options:
        if not isinstance(profile_options[ProfileOption.YLABEL_DATA_PROFILE], str):
            raise TypeError(f'Option {ProfileOption.YLABEL_DATA_PROFILE} must be a string.')

    if ProfileOption.XLABEL_LOG_RATIO_PROFILE in profile_options:
        if not isinstance(profile_options[ProfileOption.XLABEL_LOG_RATIO_PROFILE], str):
            raise TypeError(f'Option {ProfileOption.XLABEL_LOG_RATIO_PROFILE} must be a string.')

    if ProfileOption.YLABEL_LOG_RATIO_PROFILE in profile_options:
        if not isinstance(profile_options[ProfileOption.YLABEL_LOG_RATIO_PROFILE], str):
            raise TypeError(f'Option {ProfileOption.YLABEL_LOG_RATIO_PROFILE} must be a string.')

    return profile_options


def check_post_load_profile_options(n_solvers, profile_options):
    """
    Finish the solver-count-dependent option checks that were deferred by
    `check_validity_profile_options` because the results were loaded from disk
    (`solvers` was None there). `n_solvers` is the final effective solver count
    after the 'load' and 'solvers_to_load' options have been applied.
    """
    if ProfileOption.SOLVER_NAMES in profile_options:
        if len(profile_options[ProfileOption.SOLVER_NAMES]) != n_solvers:
            raise ValueError(f'Option {ProfileOption.SOLVER_NAMES} must have the same length as the number of loaded solvers ({n_solvers}).')

    if ProfileOption.SOLVER_ISRAND in profile_options:
        if len(profile_options[ProfileOption.SOLVER_ISRAND]) != n_solvers:
            raise ValueError(f'Option {ProfileOption.SOLVER_ISRAND} must have the same length as the number of loaded solvers ({n_solvers}).')

    if profile_options.get(ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES, False) and n_solvers != 2:
        print_log_message('WARNING', 'The log-ratio profiles are available only when there are exactly two solvers. We will not generate the log-ratio profiles.')
        profile_options[ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES] = False

    return profile_options


def get_default_problem_options(problem_options):
    """
    Get the default problem options.
    """
    problem_options.setdefault(ProblemOption.PLIBS.value, ['s2mpj'])
    problem_options.setdefault(ProblemOption.PLIB_OPTIONS.value, {})
    problem_options.setdefault(ProblemOption.PTYPE.value, 'u')
    problem_options.setdefault(ProblemOption.MINDIM.value, 1)
    problem_options.setdefault(ProblemOption.MAXDIM.value, problem_options[ProblemOption.MINDIM] + 1)
    problem_options.setdefault(ProblemOption.MINB.value, 0)
    problem_options.setdefault(ProblemOption.MAXB.value, problem_options[ProblemOption.MINB] + 10)
    problem_options.setdefault(ProblemOption.MINLCON.value, 0)
    problem_options.setdefault(ProblemOption.MAXLCON.value, problem_options[ProblemOption.MINLCON] + 10)
    problem_options.setdefault(ProblemOption.MINNLCON.value, 0)
    problem_options.setdefault(ProblemOption.MAXNLCON.value, problem_options[ProblemOption.MINNLCON] + 10)
    problem_options.setdefault(ProblemOption.MINCON.value, min(problem_options[ProblemOption.MINLCON], problem_options[ProblemOption.MINNLCON]))
    problem_options.setdefault(ProblemOption.MAXCON.value, max(problem_options[ProblemOption.MAXLCON], problem_options[ProblemOption.MAXNLCON]))
    problem_options.setdefault(ProblemOption.PROBLEM_NAMES.value, [])
    problem_options.setdefault(ProblemOption.EXCLUDELIST.value, [])

    return problem_options


def get_default_profile_options(solvers, feature, profile_options):
    """
    Get the default profile options.
    """
    # Remember whether the user explicitly chose a history-plot mode; the
    # conflict resolution below must not override an explicit choice.
    user_draw_hist_plots = profile_options.get(ProfileOption.DRAW_HIST_PLOTS.value)
    # Use a conservative default instead of all available workers.
    profile_options.setdefault(ProfileOption.N_JOBS.value, _get_conservative_default_n_jobs())
    profile_options.setdefault(ProfileOption.SEED.value, 0)
    profile_options.setdefault(ProfileOption.BENCHMARK_ID.value, 'out')
    if solvers is not None:
        profile_options.setdefault(ProfileOption.SOLVER_NAMES.value, [getattr(s, '__name__', f'Solver {i+1}') for i, s in enumerate(solvers)])
        profile_options.setdefault(ProfileOption.SOLVER_ISRAND.value, [False for _ in solvers])
    else:
        profile_options.setdefault(ProfileOption.SOLVER_NAMES.value, [])
        profile_options.setdefault(ProfileOption.SOLVER_ISRAND.value, [])
    profile_options.setdefault(ProfileOption.FEATURE_STAMP, _get_default_feature_stamp(feature))
    profile_options.setdefault(ProfileOption.ERRORBAR_TYPE.value, 'minmax')
    profile_options.setdefault(ProfileOption.HIST_AGGREGATION.value, 'min')
    profile_options.setdefault(ProfileOption.SAVEPATH.value, str(Path('.').expanduser().resolve()))
    profile_options.setdefault(ProfileOption.MAX_TOL_ORDER.value, 10)
    profile_options.setdefault(ProfileOption.MAX_EVAL_FACTOR.value, 500)
    profile_options.setdefault(ProfileOption.MERIT_FUN.value, scoring._default_merit)
    profile_options.setdefault(ProfileOption.PROJECT_X0.value, False)
    profile_options.setdefault(ProfileOption.RUN_PLAIN.value, False)
    profile_options.setdefault(ProfileOption.SCORE_ONLY.value, False)
    profile_options.setdefault(ProfileOption.DRAW_HIST_PLOTS.value, 'parallel')
    profile_options.setdefault(ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES.value, True)
    profile_options.setdefault(ProfileOption.SUMMARIZE_DATA_PROFILES.value, True)
    profile_options.setdefault(ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES.value, False)
    profile_options.setdefault(ProfileOption.SUMMARIZE_OUTPUT_BASED_PROFILES.value, True)
    profile_options.setdefault(ProfileOption.SILENT.value, False)
    profile_options.setdefault(ProfileOption.SOLVER_VERBOSE.value, 1)
    if profile_options[ProfileOption.SILENT.value]:
        profile_options[ProfileOption.SOLVER_VERBOSE.value] = 0
    profile_options.setdefault(ProfileOption.SEMILOGX.value, True)
    profile_options.setdefault(ProfileOption.NORMALIZED_SCORES.value, True)
    profile_options.setdefault(ProfileOption.SCORE_WEIGHT_FUN.value, scoring._default_score_weight_fun)
    profile_options.setdefault(ProfileOption.SCORE_FUN.value, scoring._default_score_fun)
    profile_options.setdefault(ProfileOption.LOAD.value, None)
    profile_options.setdefault(ProfileOption.LINE_COLORS.value, [
        '#1f77b4',
        '#ff7f0e',
        '#2ca02c',
        '#d62728',
        '#9467bd',
        '#8c564b',
        '#e377c2',
        '#7f7f7f',
        '#bcbd22',
        '#17becf'
    ])
    profile_options.setdefault(ProfileOption.LINE_STYLES.value, [
        '-', '-.', '--', ':'
    ])
    profile_options.setdefault(ProfileOption.LINE_WIDTHS.value, [1.5])
    profile_options.setdefault(ProfileOption.BAR_COLORS.value, profile_options[ProfileOption.LINE_COLORS])
    profile_options.setdefault(ProfileOption.PLOT_OPTIONS.value, {})
    profile_options.setdefault(ProfileOption.XLABEL_PERFORMANCE_PROFILE.value, 'Performance ratio')
    profile_options.setdefault(ProfileOption.YLABEL_PERFORMANCE_PROFILE.value, 'Performance profiles ($\\mathrm{tol} = %s$)')
    profile_options.setdefault(ProfileOption.XLABEL_DATA_PROFILE.value, 'Number of simplex gradients')
    profile_options.setdefault(ProfileOption.YLABEL_DATA_PROFILE.value, 'Data profiles ($\\mathrm{tol} = %s$)')
    profile_options.setdefault(ProfileOption.XLABEL_LOG_RATIO_PROFILE.value, 'Problem')
    profile_options.setdefault(ProfileOption.YLABEL_LOG_RATIO_PROFILE.value, 'Log-ratio profiles ($\\mathrm{tol} = %s$)')

    # Resolve potential conflicts in the options. The precedence for the
    # history-plot mode is: 'score_only' is the strongest setting and yields an
    # effective 'none'; otherwise an explicit user choice is honored exactly;
    # when the option is omitted, loading defaults to 'sequential' while fresh
    # runs keep the 'parallel' default.
    if profile_options[ProfileOption.SCORE_ONLY.value]:
        if user_draw_hist_plots is not None and user_draw_hist_plots != 'none':
            print_log_message(
                'WARNING',
                f'Since the option "score_only" is true, the requested "draw_hist_plots" mode '
                f'"{user_draw_hist_plots}" is ignored and no history plots will be drawn.',
            )
        elif not profile_options[ProfileOption.SILENT.value]:
            print_log_message(
                'INFO',
                'Since the option "score_only" is true, we will not draw any history plots of the '
                'problems (the option "draw_hist_plots" is set to "none").',
            )
        profile_options[ProfileOption.DRAW_HIST_PLOTS.value] = 'none'
    elif user_draw_hist_plots is None and profile_options[ProfileOption.LOAD.value] is not None:
        if not profile_options[ProfileOption.SILENT.value]:
            print_log_message(
                'INFO',
                'Since the option "load" is provided, we will draw history plots of the problems '
                'after loading the results (the option "draw_hist_plots" is set to "sequential").',
            )
        profile_options[ProfileOption.DRAW_HIST_PLOTS.value] = 'sequential'

    return profile_options


def create_stamp(solver_names, problem_options, feature_stamp, time_stamp, path_out):
    """
    Create a stamp for the current experiment.
    The stamp is a string that contains the solver names, problem settings, feature stamp, and time stamp.
    """

    # Set the max length of the stamp
    max_length = 100
    # Check the system os
    if os.name == 'nt':  # Windows
        max_dir_length = 250  # Windows max path length (260) minus some buffer
    else:
        max_dir_length = 4000  # Unix max path length (4096) minus some buffer
    # We want to avoid the path 'path_out/stamp/summary_stamp.pdf' exceeding the max path length
    max_length = min(max_length, (max_dir_length - len(str(path_out)) - len('summary_') - len('.pdf') - 1) // 2)


    # Get the solver stamp: replace non-alphanumeric/underscore with '_'
    solver_names_clean = [re.sub(r'[^a-zA-Z0-9_]', '_', credential_free_urls(name)) for name in solver_names]
    solver_stamp = '_'.join(solver_names_clean)

    # Get the problem stamp
    ptype = problem_options[ProblemOption.PTYPE]
    mindim = problem_options[ProblemOption.MINDIM]
    maxdim = problem_options[ProblemOption.MAXDIM]
    mincon = problem_options[ProblemOption.MINCON]
    maxcon = problem_options[ProblemOption.MAXCON]
    if ptype == 'u':
        problem_stamp = f"{ptype}_{mindim}_{maxdim}"
    else:
        problem_stamp = f"{ptype}_{mindim}_{maxdim}_{mincon}_{maxcon}"

    # Create the final stamp
    stamp = time_stamp
    if len(feature_stamp) + len(stamp) + 1 <= max_length:
        stamp = f"{feature_stamp}_{stamp}"
    if len(problem_stamp) + len(stamp) + 1 <= max_length:
        stamp = f"{problem_stamp}_{stamp}"
    if len(solver_stamp) + len(stamp) + 1 <= max_length:
        stamp = f"{solver_stamp}_{stamp}"

    return stamp


# Longest default feature stamp used as a path component. Longer composite
# stamps are cut and qualified with a digest of the full stamp, which stays
# in the archive and report metadata.
_FEATURE_STAMP_LIMIT = 64


def _bounded_feature_stamp(stamp):
    """Bound a stamp for use as a path component, keeping it collision-safe."""
    if len(stamp) <= _FEATURE_STAMP_LIMIT:
        return stamp
    digest = hashlib.sha256(stamp.encode('utf-8')).hexdigest()[:8]
    return f"{stamp[:_FEATURE_STAMP_LIMIT - 9].rstrip('_')}_{digest}"


def _get_default_feature_stamp(feature, bounded=True):
    """
    Get the default feature stamp of a specification.
    A one-stage feature keeps its established stamp; a composition joins the
    stamps of its effective stages, in order, with a double underscore and,
    with ``bounded=True`` (the default), a long composite stamp is cut and
    digest-qualified for use in paths. The identity feature is ``'plain'``.
    """
    feature_stamp = full_feature_stamp(feature)
    if len(feature.stages) <= 1:
        return feature_stamp
    return _bounded_feature_stamp(feature_stamp) if bounded else feature_stamp


def merge_pdfs_with_pypdf(input_dir, output_file):
    input_dir = Path(input_dir)
    pdfs = sorted(input_dir.glob("*.pdf"))
    if not pdfs:
        raise FileNotFoundError(f"No PDFs found in {input_dir}")
    merger = PdfWriter()
    for pdf in pdfs:
        merger.append(str(pdf))
    with open(output_file, "wb") as f:
        merger.write(f)
    merger.close()


def _write_history_display_note(stream, results_plibs):
    """Report display transforms separately from the untouched scientific data."""
    from .plotting import _HISTORY_DISPLAY_LIMIT

    stream.write(
        f"History display only: finite magnitudes above {_HISTORY_DISPLAY_LIMIT:.0e} "
        "are clipped before plotting; NaN/Inf use finite display placeholders.\n"
        "Raw data, oracle values and scores are not clipped. Counts below refer "
        "to stored history entries, including padded tails.\n"
    )
    for record in results_plibs:
        for key in ('fun_histories', 'maxcv_histories', 'merit_histories'):
            clipped = nonfinite = 0
            # Count one problem at a time rather than allocating another
            # full experiment-sized mask for a potentially multi-GB archive.
            for values in np.asarray(record.get(key, [])):
                finite = np.isfinite(values)
                clipped += np.count_nonzero(finite & (np.abs(values) > _HISTORY_DISPLAY_LIMIT))
                nonfinite += np.count_nonzero(~finite)
            if clipped or nonfinite:
                stream.write(f"  {record['plib']} / {key}: {clipped} clipped, "
                             f"{nonfinite} nonfinite entries.\n")
    stream.write('\n')


def write_report(profile_options, results_plibs, path_report, path_readme_log):
    """
    Writes the report of the current experiment to a file.
    """
    logger = get_logger(__name__)
    if profile_options[ProfileOption.SCORE_ONLY]:
        return
    try:
        with open(path_report, 'w') as fid:
            fid.write("# Report file for the current experiment\n\n")
            fid.write("This report file is generated by OptiProfiler <www.optprof.com>\n\n")
            fid.write("## Summary of the experiment\n\n")
            solver_names = results_plibs[0]['solver_names']
            fid.write(f"Solver names:            {', '.join(solver_names)}\n")
            fid.write(f"Problem types:           {results_plibs[0]['ptype']}\n")
            fid.write(f"Problem mindim:          {results_plibs[0]['mindim']}\n")
            fid.write(f"Problem maxdim:          {results_plibs[0]['maxdim']}\n")
            if any(c in results_plibs[0]['ptype'] for c in 'bln'):
                fid.write(f"Problem minb:            {results_plibs[0]['minb']}\n")
                fid.write(f"Problem maxb:            {results_plibs[0]['maxb']}\n")
            if any(c in results_plibs[0]['ptype'] for c in 'ln'):
                fid.write(f"Problem minlcon:         {results_plibs[0]['minlcon']}\n")
                fid.write(f"Problem maxlcon:         {results_plibs[0]['maxlcon']}\n")
            if 'n' in results_plibs[0]['ptype']:
                fid.write(f"Problem minnlcon:        {results_plibs[0]['minnlcon']}\n")
                fid.write(f"Problem maxnlcon:        {results_plibs[0]['maxnlcon']}\n")
            if any(c in results_plibs[0]['ptype'] for c in 'ln'):
                fid.write(f"Problem mincon:          {results_plibs[0]['mincon']}\n")
                fid.write(f"Problem maxcon:          {results_plibs[0]['maxcon']}\n")
            if results_plibs[0].get('problem_names_options'):
                fid.write(f"Problem names from user: {', '.join(results_plibs[0]['problem_names_options'])}\n")
            else:
                fid.write("Problem names from user: This part is empty.\n")
            if results_plibs[0].get('excludelist'):
                fid.write(f"Exclude list from user:  {', '.join(results_plibs[0]['excludelist'])}\n")
            else:
                fid.write("Exclude list from user:  This part is empty.\n")
            fid.write(f"Feature stamp:           {results_plibs[0]['feature_stamp']}\n\n")
            if profile_options.get(ProfileOption.DRAW_HIST_PLOTS, 'none') != 'none':
                _write_history_display_note(fid, results_plibs)

            for results_plib in results_plibs:
                plib = results_plib['plib']
                problem_names = results_plib['problem_names']
                problem_types = results_plib['problem_types']
                problem_dims = results_plib['problem_dims']
                problem_mbs = results_plib['problem_mbs']
                problem_mlcons = results_plib['problem_mlcons']
                problem_mnlcons = results_plib['problem_mnlcons']
                problem_mcons = results_plib['problem_mcons']
                solvers_successes = results_plib['solvers_successes']
                computation_times = results_plib['computation_times']

                # If plain results exist and RUN_PLAIN is set, merge computation times
                if 'results_plib_plain' in results_plib and profile_options[ProfileOption.RUN_PLAIN]:
                    results_plib_plain = results_plib['results_plib_plain']
                    problem_names_plain = results_plib_plain['problem_names']
                    computation_times_plain = results_plib_plain['computation_times']
                    if computation_times_plain.ndim != 3:
                        raise ValueError(
                            'The plain computation times must be a three-dimensional array.'
                        )
                    if computation_times_plain.ndim != computation_times.ndim:
                        raise ValueError(
                            'The plain computation times must have the same number of dimensions as '
                            'the main computation times.'
                        )
                    plain_runs = computation_times_plain.shape[2]
                    padding = np.full(computation_times.shape[:2] + (plain_runs,), np.nan)
                    computation_times = np.concatenate((computation_times, padding), axis=2)
                    for i_problem, pname in enumerate(problem_names):
                        if pname in problem_names_plain:
                            idx = problem_names_plain.index(pname)
                            computation_times[i_problem, :, -plain_runs:] = np.asarray(
                                computation_times_plain[idx, :, :]
                            )

                # Pick out unsolved problems and calculate computation times for each problem
                unsolved_problems = []
                time_processes = np.zeros(len(problem_names))
                for i_problem, pname in enumerate(problem_names):
                    if np.all(~np.array(solvers_successes[i_problem])):
                        unsolved_problems.append(pname)
                    time_process = np.array(computation_times[i_problem]).flatten()
                    time_processes[i_problem] = np.nansum(time_process)
                idx = np.argsort([name.lower() for name in problem_names])
                sorted_problem_names = [problem_names[i] for i in idx]
                sorted_problem_types = [problem_types[i] for i in idx]
                sorted_problem_dims = [problem_dims[i] for i in idx]
                sorted_problem_mbs = [problem_mbs[i] for i in idx]
                sorted_problem_mlcons = [problem_mlcons[i] for i in idx]
                sorted_problem_mnlcons = [problem_mnlcons[i] for i in idx]
                sorted_problem_mcons = [problem_mcons[i] for i in idx]
                sorted_time_processes = [time_processes[i] for i in idx]

                max_name_length = max(max(len(name) for name in sorted_problem_names), 12)
                max_type_length = max(max(len(str(t)) for t in sorted_problem_types), 4)
                max_dim_length = max(max(len(str(d)) for d in sorted_problem_dims), 9)
                max_mbs_length = max(max(len(str(mb)) for mb in sorted_problem_mbs), 2)
                max_mlcons_length = max(max(len(str(mlc)) for mlc in sorted_problem_mlcons), 5)
                max_mnlcons_length = max(max(len(str(mnlc)) for mnlc in sorted_problem_mnlcons), 6)
                max_mcons_length = max(max(len(str(mc)) for mc in sorted_problem_mcons), 4)
                max_time_length = max(max(len(f"{tp:.2f}") for tp in sorted_time_processes), 28)

                fid.write(f'## Report for the problem library "{plib}"\n\n')
                fid.write(f'Number of problems selected: {len(sorted_problem_names)}\n')
                fid.write(f'Wall-clock time spent by all the solvers: {np.nansum(time_processes):.2f} secs\n\n')

                if len(unsolved_problems) < len(sorted_problem_names):
                    header = (
                        f"{'Problem name':<{max_name_length}}    "
                        f"{'Type':<{max_type_length}}    "
                        f"{'Dimension':<{max_dim_length}}    "
                        f"{'mb':<{max_mbs_length}}    "
                        f"{'mlcon':<{max_mlcons_length}}    "
                        f"{'mnlcon':<{max_mnlcons_length}}    "
                        f"{'mcon':<{max_mcons_length}}    "
                        f"{'Time spent by solvers (secs)':<{max_time_length}}\n"
                    )
                    fid.write(header)
                    for i, name in enumerate(sorted_problem_names):
                        if name in unsolved_problems:
                            continue
                        type_ = sorted_problem_types[i]
                        dim = str(sorted_problem_dims[i])
                        mb = str(sorted_problem_mbs[i])
                        mlcon = str(sorted_problem_mlcons[i])
                        mnlcon = str(sorted_problem_mnlcons[i])
                        mcon = str(sorted_problem_mcons[i])
                        time = f"{sorted_time_processes[i]:.2f}"
                        line = (
                            f"{name:<{max_name_length}}    "
                            f"{type_:<{max_type_length}}    "
                            f"{dim:<{max_dim_length}}    "
                            f"{mb:<{max_mbs_length}}    "
                            f"{mlcon:<{max_mlcons_length}}    "
                            f"{mnlcon:<{max_mnlcons_length}}    "
                            f"{mcon:<{max_mcons_length}}    "
                            f"{time:<{max_time_length}}\n"
                        )
                        fid.write(line)
                else:
                    fid.write("This part is empty.\n")
                fid.write(f'\n## Problems from the problem library "{plib}" that all the solvers failed to evaluate a single point\n\n')
                if unsolved_problems:
                    for name in unsolved_problems:
                        fid.write(f"{name} ")
                    fid.write("\n")
                else:
                    fid.write("This part is empty.\n")
        # Optionally update the readme log
        add_to_readme(path_readme_log, 'report.txt', 'File, the report file of the current experiment, recording problem names, timing information, and diagnostics such as merit_init = phi(x_0) = inf cases, abnormal solver terminations, and output fallbacks.')
    except Exception as exc:
        if not profile_options[ProfileOption.SILENT]:
            logger.warning(f'Error occurred when writing the report to {path_report}.')
            logger.warning(f'Error message: {shorten_log_message(exc)}')


def init_readme(path_readme):
    """
    Initialize the README file with a header.
    """
    width = 100
    name_width = 33
    type_width = 8
    gap = 2
    try:
        with open(path_readme, 'w') as f:
            f.write("# Content of this folder\n\n")
            header = f"{'File/Folder Name':<{name_width}}{' ' * gap}{'Type':<{type_width}}{' ' * gap}{'Description'}\n"
            f.write(header)
            f.write('-' * width + '\n')
    except Exception:
        pass


def add_to_readme(path_readme, filename, description):
    """
    Add an entry to the README file.
    """
    width = 100
    name_width = 33
    type_width = 8
    gap = 2
    desc_width = width - name_width - type_width - 2 * gap

    # Parse type and content from description
    type_str = ""
    content_str = description
    if description.startswith("File, "):
        type_str = "File"
        content_str = description[6:]
    elif description.startswith("Folder, "):
        type_str = "Folder"
        content_str = description[8:]
    
    # Capitalize the first letter of the content
    if content_str:
        content_str = content_str[0].upper() + content_str[1:]

    try:
        with open(path_readme, 'a') as f:
            wrapped_filename = textwrap.wrap(filename, width=name_width)
            wrapped_desc = textwrap.wrap(content_str, width=desc_width)
            
            if not wrapped_filename:
                wrapped_filename = ['']
            
            n_lines = max(len(wrapped_filename), len(wrapped_desc))
            
            for i in range(n_lines):
                f_part = wrapped_filename[i] if i < len(wrapped_filename) else ''
                # Type only appears on the first line
                t_part = type_str if i == 0 else ''
                d_part = wrapped_desc[i] if i < len(wrapped_desc) else ''
                
                line = f"{f_part:<{name_width}}{' ' * gap}{t_part:<{type_width}}{' ' * gap}{d_part}\n"
                f.write(line)
    except Exception:
        pass
