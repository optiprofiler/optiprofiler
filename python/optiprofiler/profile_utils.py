import os
import re
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import is_color_like
from pypdf import PdfWriter


from .utils import FeatureName, ProfileOption, FeatureOption, ProblemOption, get_logger






def check_validity_problem_options(problem_options):
    """
    Check the validity of the problem options.
    """
    if ProblemOption.PLIBS in problem_options:
        if not isinstance(problem_options[ProblemOption.PLIBS], list):
            problem_options[ProblemOption.PLIBS] = [problem_options[ProblemOption.PLIBS]]
        mydir = Path(__file__).parent.resolve()
        problem_dir = Path(mydir, 'problem_libs').resolve()
        subfolders = [p.name for p in problem_dir.iterdir() if p.is_dir()]
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
        if len(problem_options[ProblemOption.PTYPE]) == 0:
            raise ValueError(f'Option {ProblemOption.PTYPE} cannot be an empty string.')
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


def check_validity_profile_options(solvers, profile_options):
    """
    Check the validity of the profile options.
    """
    logger = get_logger(__name__)
    if ProfileOption.N_JOBS in profile_options:
        if isinstance(profile_options[ProfileOption.N_JOBS], float) and profile_options[ProfileOption.N_JOBS].is_integer():
            profile_options[ProfileOption.N_JOBS] = int(profile_options[ProfileOption.N_JOBS])
        if not isinstance(profile_options[ProfileOption.N_JOBS], int):
            raise TypeError(f'Option {ProfileOption.N_JOBS} must be an integer.')
        if profile_options[ProfileOption.N_JOBS] < 1:
            profile_options[ProfileOption.N_JOBS] = 1
            logger.warning(f'Option {ProfileOption.N_JOBS} is set to 1 because it cannot be smaller than 1.')

    if ProfileOption.SEED in profile_options:
        if isinstance(profile_options[ProfileOption.SEED], float) and profile_options[ProfileOption.SEED].is_integer():
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
        if len(profile_options[ProfileOption.SOLVER_NAMES]) != len(solvers):
            raise ValueError(f'Option {ProfileOption.SOLVER_NAMES} must have the same length as the number of solvers ({len(solvers)}).')
        if any(not isinstance(name, str) for name in profile_options[ProfileOption.SOLVER_NAMES]):
            raise TypeError(f'Option {ProfileOption.SOLVER_NAMES} must be a list of strings.')

    if ProfileOption.SOLVER_ISRAND in profile_options:
        if not isinstance(profile_options[ProfileOption.SOLVER_ISRAND], list):
            raise TypeError(f'Option {ProfileOption.SOLVER_ISRAND} must be a list of booleans.')
        if len(profile_options[ProfileOption.SOLVER_ISRAND]) != len(solvers):
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
        if isinstance(profile_options[ProfileOption.MAX_TOL_ORDER], float) and profile_options[ProfileOption.MAX_TOL_ORDER].is_integer():
            profile_options[ProfileOption.MAX_TOL_ORDER] = int(profile_options[ProfileOption.MAX_TOL_ORDER])
        if not isinstance(profile_options[ProfileOption.MAX_TOL_ORDER], int):
            raise TypeError(f'Option {ProfileOption.MAX_TOL_ORDER} must be an integer.')
        if profile_options[ProfileOption.MAX_TOL_ORDER] < 0 or profile_options[ProfileOption.MAX_TOL_ORDER] > 16:
            raise ValueError(f'Option {ProfileOption.MAX_TOL_ORDER} must be a positive integer not larger than 16.')

    if ProfileOption.MAX_EVAL_FACTOR in profile_options:
        if not isinstance(profile_options[ProfileOption.MAX_EVAL_FACTOR], (float, int)):
            raise TypeError(f'Option {ProfileOption.MAX_EVAL_FACTOR} must be a float or an integer.')
        if profile_options[ProfileOption.MAX_EVAL_FACTOR] <= 0:
            raise ValueError(f'Option {ProfileOption.MAX_EVAL_FACTOR} must be positive.')

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

    if ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES in profile_options:
        if not isinstance(profile_options[ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES], bool):
            raise TypeError(f'Option {ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES} must be a boolean.')

    if ProfileOption.SUMMARIZE_DATA_PROFILES in profile_options:
        if not isinstance(profile_options[ProfileOption.SUMMARIZE_DATA_PROFILES], bool):
            raise TypeError(f'Option {ProfileOption.SUMMARIZE_DATA_PROFILES} must be a boolean.')

    if ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES in profile_options:
        if not isinstance(profile_options[ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES], bool):
            raise TypeError(f'Option {ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES} must be a boolean.')
        if profile_options[ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES]:
            if len(solvers) != 2:
                logger.warning(f'The log-ratio profiles are available only when there are exactly two solvers. We will not generate the log-ratio profiles.')
                profile_options[ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES] = False

    if ProfileOption.SUMMARIZE_OUTPUT_BASED_PROFILES in profile_options:
        if not isinstance(profile_options[ProfileOption.SUMMARIZE_OUTPUT_BASED_PROFILES], bool):
            raise TypeError(f'Option {ProfileOption.SUMMARIZE_OUTPUT_BASED_PROFILES} must be a boolean.')

    if ProfileOption.SILENT in profile_options:
        if not isinstance(profile_options[ProfileOption.SILENT], bool):
            raise TypeError(f'Option {ProfileOption.SILENT} must be a boolean.')

    if ProfileOption.SOLVER_VERBOSE in profile_options:
        if isinstance(profile_options[ProfileOption.SOLVER_VERBOSE], float) and profile_options[ProfileOption.SOLVER_VERBOSE].is_integer():
            profile_options[ProfileOption.SOLVER_VERBOSE] = int(profile_options[ProfileOption.SOLVER_VERBOSE])
        if not isinstance(profile_options[ProfileOption.SOLVER_VERBOSE], int):
            raise TypeError(f'Option {ProfileOption.SOLVER_VERBOSE} must be 0, 1, or 2.')
        if profile_options[ProfileOption.SOLVER_VERBOSE] < 0 or profile_options[ProfileOption.SOLVER_VERBOSE] > 2:
            raise ValueError(f'Option {ProfileOption.SOLVER_VERBOSE} must be 0, 1, or 2.')
        if ProfileOption.SILENT in profile_options and profile_options[ProfileOption.SILENT]:
            if profile_options[ProfileOption.SOLVER_VERBOSE] == 2:
                logger.warning(f'Option {ProfileOption.SOLVER_VERBOSE} is set to 1 because option {ProfileOption.SILENT} is True.')
                profile_options[ProfileOption.SOLVER_VERBOSE] = 1
            elif profile_options[ProfileOption.SOLVER_VERBOSE] == 1:
                logger.warning(f'Option {ProfileOption.SOLVER_VERBOSE} is set to 0 because option {ProfileOption.SILENT} is True.')
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
        if len(profile_options[ProfileOption.SOLVERS_TO_LOAD]) == 0:
            raise ValueError(f'Option {ProfileOption.SOLVERS_TO_LOAD} cannot be an empty list.')
        if ProfileOption.SOLVER_NAMES in profile_options:
            if len(profile_options[ProfileOption.SOLVERS_TO_LOAD]) != len(profile_options[ProfileOption.SOLVER_NAMES]):
                raise ValueError(f'Option {ProfileOption.SOLVERS_TO_LOAD} must have the same length as option {ProfileOption.SOLVER_NAMES} ({len(profile_options[ProfileOption.SOLVER_NAMES])}).')
        for i in range(len(profile_options[ProfileOption.SOLVERS_TO_LOAD])):
            if isinstance(profile_options[ProfileOption.SOLVERS_TO_LOAD][i], float) and profile_options[ProfileOption.SOLVERS_TO_LOAD][i].is_integer():
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
                plt.plot([0, 1], [0, 1], style)
                plt.close()
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


def get_default_problem_options(problem_options):
    """
    Get the default problem options.
    """
    problem_options.setdefault(ProblemOption.PLIBS.value, ['s2mpj'])
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


def _default_score_weight_fun(x):
    return 1

def _default_score_fun(x):
    return np.mean(x[:, :, 0, 0], axis=1)

def get_default_profile_options(solvers, feature, profile_options):
    """
    Get the default profile options.
    """
    profile_options.setdefault(ProfileOption.N_JOBS.value, os.cpu_count() or 1)
    profile_options.setdefault(ProfileOption.SEED.value, 0)
    profile_options.setdefault(ProfileOption.BENCHMARK_ID.value, 'out')
    profile_options.setdefault(ProfileOption.SOLVER_NAMES.value, [getattr(s, '__name__', f'Solver {i+1}') for i, s in enumerate(solvers)])
    profile_options.setdefault(ProfileOption.SOLVER_ISRAND.value, [False for _ in solvers])
    profile_options.setdefault(ProfileOption.FEATURE_STAMP, _get_default_feature_stamp(feature))
    profile_options.setdefault(ProfileOption.ERRORBAR_TYPE.value, 'minmax')
    profile_options.setdefault(ProfileOption.HIST_AGGREGATION.value, 'min')
    profile_options.setdefault(ProfileOption.SAVEPATH.value, str(Path('.').expanduser().resolve()))
    profile_options.setdefault(ProfileOption.MAX_TOL_ORDER.value, 10)
    profile_options.setdefault(ProfileOption.MAX_EVAL_FACTOR.value, 500)
    profile_options.setdefault(ProfileOption.MERIT_FUN.value, _default_merit)
    profile_options.setdefault(ProfileOption.PROJECT_X0.value, False)
    profile_options.setdefault(ProfileOption.RUN_PLAIN.value, False)
    profile_options.setdefault(ProfileOption.SCORE_ONLY.value, False)
    profile_options.setdefault(ProfileOption.SUMMARIZE_PERFORMANCE_PROFILES.value, True)
    profile_options.setdefault(ProfileOption.SUMMARIZE_DATA_PROFILES.value, True)
    profile_options.setdefault(ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES.value, False)
    profile_options.setdefault(ProfileOption.SUMMARIZE_OUTPUT_BASED_PROFILES.value, True)
    profile_options.setdefault(ProfileOption.SILENT.value, False)
    profile_options.setdefault(ProfileOption.SOLVER_VERBOSE.value, 1)
    profile_options.setdefault(ProfileOption.SEMILOGX.value, True)
    profile_options.setdefault(ProfileOption.NORMALIZED_SCORES.value, False)
    profile_options.setdefault(ProfileOption.SCORE_WEIGHT_FUN.value, _default_score_weight_fun)
    profile_options.setdefault(ProfileOption.SCORE_FUN.value, _default_score_fun)
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
    profile_options.setdefault(ProfileOption.XLABEL_PERFORMANCE_PROFILE.value, 'Performance ratio')
    profile_options.setdefault(ProfileOption.YLABEL_PERFORMANCE_PROFILE.value, 'Performance profiles ($\\mathrm{tol} = %s$)')
    profile_options.setdefault(ProfileOption.XLABEL_DATA_PROFILE.value, 'Number of simplex gradients')
    profile_options.setdefault(ProfileOption.YLABEL_DATA_PROFILE.value, 'Data profiles ($\\mathrm{tol} = %s$)')
    profile_options.setdefault(ProfileOption.XLABEL_LOG_RATIO_PROFILE.value, 'Problem')
    profile_options.setdefault(ProfileOption.YLABEL_LOG_RATIO_PROFILE.value, 'Log-ratio profiles ($\\mathrm{tol} = %s$)')

    return profile_options


def compute_merit_values(merit_fun, fun_values, maxcv_values, maxcv_init):
    """
    Compute the merit function values.

    Parameters
    ----------
    merit_fun : callable
        The merit function to be used. It must have the signature merit_fun(fun_value, maxcv_value, maxcv_init).
    fun_values : array-like
        The objective function values.
    maxcv_values : array-like
        The maximum constraint violation values.
    maxcv_init : float or array-like
        The initial maximum constraint violation value(s).

    Returns
    -------
    merit_values : np.ndarray
        The merit function values.
    """
    
    # Scale `maxcv_init` to match the dimensions of `fun_values`.
    copied_dim = np.shape(fun_values)
    if np.isscalar(maxcv_init):
        maxcv_init_expanded = np.full(copied_dim, maxcv_init)
    else:
        maxcv_init = np.asarray(maxcv_init)
        ndim = len(copied_dim)
        shaped = maxcv_init.reshape(-1, *([1] * (ndim - 1)))
        reps = (1,) + copied_dim[1:]
        maxcv_init_expanded = np.tile(shaped, reps)

    vectorized_merit_fun = np.vectorize(merit_fun)
    merit_values = vectorized_merit_fun(fun_values, maxcv_values, maxcv_init_expanded)

    return merit_values


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
    solver_names_clean = [re.sub(r'[^a-zA-Z0-9_]', '_', name) for name in solver_names]
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


def _get_default_feature_stamp(feature):
    """
    Get the default feature stamp.
    """
    name = feature.name
    if name == FeatureName.PERTURBED_X0:
        # feature_name + perturbation_level + (distribution if it is gaussian or spherical)
        feature_stamp = f"{feature.name}_{feature.options[FeatureOption.PERTURBATION_LEVEL]}"
        dist = feature.options.get(FeatureOption.DISTRIBUTION.value)
        if isinstance(dist, str) and dist in ('gaussian', 'spherical'):
            feature_stamp = f"{feature_stamp}_{dist}"
    elif name == FeatureName.NOISY:
        # feature_name + noise_level + noise_type + (distribution if it is gaussian or uniform)
        feature_stamp = f"{feature.name}_{feature.options[FeatureOption.NOISE_LEVEL]}_{feature.options[FeatureOption.NOISE_TYPE]}"
        dist = feature.options.get(FeatureOption.DISTRIBUTION.value)
        if isinstance(dist, str) and dist in ('gaussian', 'uniform'):
            feature_stamp = f"{feature_stamp}_{dist}"
    elif name == FeatureName.TRUNCATED:
        # feature_name + significant_digits + (perturbed_trailing_digits if it is true)
        feature_stamp = f"{feature.name}_{feature.options[FeatureOption.SIGNIFICANT_DIGITS]}"
        if feature.options[FeatureOption.PERTURBED_TRAILING_DIGITS]:
            feature_stamp = f"{feature_stamp}_perturbed_trailing_digits"
    elif name == FeatureName.LINEARLY_TRANSFORMED:
        # feature_name + (rotated if it is true) + (condition_factor if it is not 0)
        feature_stamp = feature.name
        if feature.options[FeatureOption.ROTATED]:
            feature_stamp = f"{feature_stamp}_rotated"
        if feature.options[FeatureOption.CONDITION_FACTOR] != 0:
            feature_stamp = f"{feature_stamp}_cond{feature.options[FeatureOption.CONDITION_FACTOR]}"
    elif name == FeatureName.RANDOM_NAN:
        # feature_name + nan_rate
        feature_stamp = f"{feature.name}_{feature.options[FeatureOption.NAN_RATE]}"
    elif name == FeatureName.UNRELAXABLE_CONSTRAINTS:
        # feature_name + (bounds if it is true) + (linear if it is true) + (nonlinear if it is true)
        feature_stamp = feature.name
        if feature.options[FeatureOption.UNRELAXABLE_BOUNDS]:
            feature_stamp = f"{feature_stamp}_bounds"
        if feature.options[FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS]:
            feature_stamp = f"{feature_stamp}_linear"
        if feature.options[FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS]:
            feature_stamp = f"{feature_stamp}_nonlinear"
    elif name == FeatureName.QUANTIZED:
        # feature_name + mesh_size + (ground_truth if is_true it is true)
        feature_stamp = f"{feature.name}_{feature.options[FeatureOption.MESH_SIZE]}"
        if feature.options[FeatureOption.GROUND_TRUTH]:
            feature_stamp = f"{feature_stamp}_ground_truth"
    else:
        feature_stamp = feature.name

    return feature_stamp


def _default_merit(fun_value, maxcv_value, maxcv_init):
    """
    The merit function varphi(x) is defined by the objective function f(x) and the maximum constraint violation v(x):
        varphi(x) = f(x)                        if v(x) <= v1
        varphi(x) = f(x) + 1e5 * (v(x) - v1)    if v1 < v(x) <= v2
        varphi(x) = np.inf                      if v(x) > v2
    where v1 = min(0.01, 1e-10 * max(1, v0)), v2 = max(0.1, 2 * v0), and v0 is the initial maximum constraint violation.

    Special case 1: _default_merit(np.nan, maxcv_value, maxcv_init) = np.inf
    Sepcial case 2: _default_merit(fun_value, np.nan, maxcv_init) = np.inf
    """
    if np.isnan(fun_value) or np.isnan(maxcv_value):
        return np.inf
    tolerance = 1e-10
    v1 = min(0.01, tolerance * max(1, maxcv_init))
    v2 = max(0.1, 2 * maxcv_init)
    if maxcv_value <= v1:
        return fun_value
    if v1 < maxcv_value <= v2:
        return fun_value + 1e5 * (maxcv_value - v1)
    else:
        return np.inf


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
                    import numpy as np
                    computation_times = np.concatenate(
                        (computation_times, np.full_like(computation_times_plain, np.nan)), axis=2
                    )
                    for i_problem, pname in enumerate(problem_names):
                        if pname in problem_names_plain:
                            idx = problem_names_plain.index(pname)
                            computation_times[i_problem, :, -1] = computation_times_plain[idx, :, :]

                # Pick out unsolved problems and calculate computation times for each problem
                unsolved_problems = []
                import numpy as np
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
        try:
            with open(path_readme_log, 'a') as fid_log:
                fid_log.write("'report.txt': file, the report file of the current experiment, recording information like problem names and time spent on solving each problem for all the problem libraries.\n")
        except Exception:
            pass
    except Exception as exc:
        if not profile_options[ProfileOption.SILENT]:
            logger.warning(f'Error occurred when writing the report to {path_report}: {exc}')


def process_results(results_plibs, profile_options):
    """
    Processes results from results_plibs.
    Returns:
        merit_histories_merged, merit_outs_merged, merit_inits_merged, merit_mins_merged, n_evals_merged, problem_names_merged, problem_dims_merged
    """
    def unify_length(results_plibs, name):
        # Unify the length of the specified field in results_plibs (merit_histories or merit_histories_plain)
        max_length = max(r[name].shape[3] for r in results_plibs)
        for r in results_plibs:
            field = r[name]
            shape = list(field.shape)
            shape[3] = max_length
            field_unified = np.full(shape, np.nan)
            field_unified[:, :, :, :field.shape[3]] = field
            if field.shape[3] < max_length:
                field_unified[:, :, :, field.shape[3]:] = np.repeat(field[:, :, :, -1:], max_length - field.shape[3], axis=3)
            r[name] = field_unified
        return results_plibs

    merit_histories_merged = []
    merit_outs_merged = []
    merit_inits_merged = []
    n_evals_merged = []
    problem_names_merged = []
    problem_dims_merged = []

    results_plibs = unify_length(results_plibs, 'merit_histories')
    for results_plib in results_plibs:
        merit_histories_merged.append(results_plib['merit_histories'])
        merit_outs_merged.append(results_plib['merit_outs'])
        merit_inits_merged.append(results_plib['merit_inits'])
        n_evals_merged.append(results_plib['n_evals'])
        problem_dims_merged.append(results_plib['problem_dims'])
        problem_names_merged.extend(results_plib['problem_names'])
        if 'results_plib_plain' in results_plib:
            results_plib['merit_histories_plain'] = results_plib['results_plib_plain']['merit_histories']

    merit_histories_merged = np.concatenate(merit_histories_merged, axis=0)
    merit_outs_merged = np.concatenate(merit_outs_merged, axis=0)
    merit_inits_merged = np.concatenate(merit_inits_merged, axis=0)
    n_evals_merged = np.concatenate(n_evals_merged, axis=0)
    problem_dims_merged = np.concatenate(problem_dims_merged, axis=0)

    # Find the least merit value for each problem
    merit_mins_merged = np.nanmin(np.nanmin(np.nanmin(merit_histories_merged, axis=3), axis=2), axis=1)
    merit_mins_merged = np.minimum(merit_mins_merged, merit_inits_merged)

    # If results_plib_plain exists and run_plain is True, redefine merit_mins_merged
    if 'results_plib_plain' in results_plibs[0] and profile_options[ProfileOption.RUN_PLAIN]:
        merit_histories_plain_merged = []
        merit_inits_plain_merged = []
        problem_names_plain_merged = []
        results_plibs = unify_length(results_plibs, 'merit_histories_plain')
        for results_plib in results_plibs:
            merit_histories_plain_merged.append(results_plib['merit_histories_plain'])
            merit_inits_plain_merged.append(results_plib['results_plib_plain']['merit_inits'])
            problem_names_plain_merged.extend(results_plib['results_plib_plain']['problem_names'])
        merit_histories_plain_merged = np.concatenate(merit_histories_plain_merged, axis=0)
        merit_inits_plain_merged = np.concatenate(merit_inits_plain_merged, axis=0)
        merit_mins_plain_merged = np.nanmin(np.nanmin(np.nanmin(merit_histories_plain_merged, axis=3), axis=2), axis=1)
        merit_mins_plain_merged = np.minimum(merit_mins_plain_merged, merit_inits_plain_merged)
        for i_problem, name in enumerate(problem_names_merged):
            if name in problem_names_plain_merged:
                idx = problem_names_plain_merged.index(name)
                merit_mins_merged[i_problem] = np.nanmin([merit_mins_merged[i_problem], merit_mins_plain_merged[idx]])

    return (merit_histories_merged, merit_outs_merged, merit_inits_merged, merit_mins_merged, n_evals_merged, problem_names_merged, problem_dims_merged)