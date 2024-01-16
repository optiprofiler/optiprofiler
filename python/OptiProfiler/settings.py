from enum import Enum


class FeatureName(str, Enum):
    """
    Feature names.
    """
    CUSTOM = 'custom'
    NOISY = 'noisy'
    PLAIN = 'plain'
    RANDOMIZE_X0 = 'randomize_x0'
    REGULARIZED = 'regularized'
    TOUGH = 'tough'
    TRUNCATED = 'truncated'


class ProfileOptionKey(str, Enum):
    """
    Profile's options.
    """
    N_JOBS = 'n_jobs'
    BENCHMARK_ID = 'benchmark_id'


class CUTEstProblemOptionKey(str, Enum):
    """
    Options for loading a CUTEst problem.
    """
    N_MIN = 'n_min'
    N_MAX = 'n_max'
    M_MIN = 'm_min'
    M_MAX = 'm_max'


class FeatureOptionKey(str, Enum):
    """
    Feature's options.
    """
    DISTRIBUTION = 'distribution'
    MODIFIER = 'modifier'
    N_RUNS = 'n_runs'
    ORDER = 'order'
    PARAMETER = 'parameter'
    RATE_ERROR = 'rate_error'
    RATE_NAN = 'rate_nan'
    SIGNIFICANT_DIGITS = 'significant_digits'
    TYPE = 'type'


class NoiseType(str, Enum):
    """
    Noise types.
    """
    ABSOLUTE = 'absolute'
    RELATIVE = 'relative'


class ProblemError(Exception):
    """
    Exception raised when a problem cannot be loaded.
    """
    pass
