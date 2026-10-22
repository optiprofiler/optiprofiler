"""Shared scalar normalization; each caller owns its admissible range."""

import numpy as np


def normalize_integer(value):
    """Convert integral NumPy scalars and real floats, without changing other inputs."""
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return int(value)
    if isinstance(value, np.integer):
        return int(value)
    return value
