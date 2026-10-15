"""Shared reference-history failure policy for feature recorders."""

import logging
import warnings
import numpy as np


def reference_maxcv_or_nan(evaluate, x, owner=None):
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
        warnings.warn(message, RuntimeWarning, stacklevel=3)
        return np.nan
