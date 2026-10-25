"""Shared reference-history failure policy for feature recorders."""

import logging
import warnings
import numpy as np


def initialize_histories(owner):
    """Initialize history buffers and cached observations for one featured trial."""
    owner._fun_hist = []
    owner._cub_hist = []
    owner._ceq_hist = []
    owner._maxcv_hist = []
    owner._reference_history_warned = False
    owner._last_fun = np.nan
    owner._last_cub = np.nan
    owner._last_ceq = np.nan


def record_objective(owner, value, evaluate_maxcv, x):
    """Commit paired reference observations, including when diagnostics raise."""
    owner._fun_hist.append(value)
    # Reserve the paired slot before a warning can become an exception.
    owner._maxcv_hist.append(np.nan)
    owner._maxcv_hist[-1] = reference_maxcv_or_nan(evaluate_maxcv, x, owner, stacklevel=4)


def reference_maxcv_or_nan(evaluate, x, owner=None, *, stacklevel=3):
    """Evaluate reference violation, honoring the caller's warning policy.

    Recorders reserve a NaN history slot before calling: the diagnostic may
    raise under a warnings-as-errors filter, but histories must remain aligned.
    """
    try:
        return evaluate(x)
    except Exception as error:
        # One diagnostic per featured trial, independent of warning filters.
        if owner is not None:
            if owner._reference_history_warned:
                return np.nan
            owner._reference_history_warned = True
        message = (f'Reference constraint violation could not be recorded: '
                   f'{type(error).__name__}: {error}')
        # Solver execution suppresses Python warnings; its queue logger must
        # still retain this failure in the benchmark log.
        logging.getLogger(__name__).warning(message)
        warnings.warn(message, RuntimeWarning, stacklevel=stacklevel)
        return np.nan
