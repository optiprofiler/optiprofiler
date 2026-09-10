"""
Private, dependency-neutral metadata encoder shared by the evaluation report
and the feature provenance.

The encoder never executes user code: callables are described by function or
class name only, text is bounded and redacted without invoking object
representations, and anything else is recorded as an explicit reason record.
This module imports nothing from the package so that both the composition
engine and the report collector can use it without a circular dependency.
"""

import math
import os
import re
import types
from enum import Enum
from pathlib import Path, PureWindowsPath

import numpy as np

_SECRET_KEY = re.compile(
    r'password|passwd|secret|credential|authorization|api[_-]?key|'
    r'(?:^|[_-])token(?:$|[_-])|access[_-]?token|refresh[_-]?token|'
    r'private[_-]?key|cookie', re.I)
_ABS_IN_TEXT = re.compile(r'(?<![\w])(?:/(?:Users|home|private|tmp|var|mnt|opt)/\S+|[A-Za-z]:[\\/][^\s]+)')


def bounded_text(value, limit=256):
    """Bound known text without invoking object representations."""
    if isinstance(value, Enum):
        value = value.value
    if not isinstance(value, str):
        return None
    value = ''.join(c if c >= ' ' else ' ' for c in value)
    if value.startswith(('http://', 'https://')):
        return '[redacted_url]'
    if os.path.isabs(value) or PureWindowsPath(value).is_absolute():
        return '[redacted_absolute_path]'
    return _ABS_IN_TEXT.sub('[redacted_absolute_path]', value)[:limit]


def describe_callback(value):
    """Read function metadata only; do not inspect descriptors or source."""
    if isinstance(value, (types.FunctionType, types.BuiltinFunctionType,
                          types.MethodType)):
        return {'kind': 'callback', 'module': bounded_text(value.__module__),
                'name': bounded_text(value.__qualname__)}
    cls = type(value)
    return {'kind': 'callback', 'module': bounded_text(cls.__module__),
            'name': bounded_text(cls.__qualname__), 'instance_state': 'not_recorded'}


def safe_metadata(value, depth=0, key=None):
    """A deliberately limited metadata encoder, not a general serializer."""
    if key and _SECRET_KEY.search(key):
        return {'value': None, 'reason': 'redacted_sensitive_option'}
    if depth > 10:
        return {'value': None, 'reason': 'metadata_depth_limit'}
    if value is None:
        return None
    if isinstance(value, Enum):
        return safe_metadata(value.value, depth, key)
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        value = float(value)
        if math.isnan(value):
            return {'value': None, 'reason': 'nan'}
        if math.isinf(value):
            return {'value': None, 'reason': 'positive_infinity' if value > 0
                    else 'negative_infinity'}
        return value
    if isinstance(value, str):
        return bounded_text(value)
    if isinstance(value, (Path, PureWindowsPath)):
        return {'name': bounded_text(value.name), 'path': None,
                'reason': 'machine_path_not_recorded'}
    if callable(value):
        return describe_callback(value)
    if type(value) is dict:
        answer = {}
        for i, (k, v) in enumerate(value.items()):
            if i >= 256:
                answer['_omission'] = {'reason': 'metadata_item_limit',
                                       'total_items': len(value)}
                break
            k = bounded_text(k)
            if k is not None:
                answer[k] = safe_metadata(v, depth + 1, k)
        return answer
    if type(value) in (list, tuple):
        if len(value) > 256:
            return {'values': [safe_metadata(v, depth + 1) for v in value[:256]],
                    'total_items': len(value), 'reason': 'metadata_item_limit'}
        return [safe_metadata(v, depth + 1) for v in value]
    if type(value) is np.ndarray:
        if value.size <= 256:
            return safe_metadata(value.tolist(), depth + 1, key)
        return {'shape': list(value.shape), 'value': None,
                'reason': 'raw_array_not_embedded'}
    return {'value': None, 'reason': 'unsupported_metadata_type',
            'type': bounded_text(type(value).__name__)}
