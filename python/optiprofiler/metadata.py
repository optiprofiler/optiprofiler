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
import zlib
from enum import Enum
from pathlib import Path, PureWindowsPath

import numpy as np

# Option names whose values never enter a report, matched without regard to
# case. `EvalReport.m` holds the same literal, doc/source/user/eval_report.rst
# describes the list, and tools/tests/test_report_redaction_contract.py keeps
# the three together. `pwd`, `auth` and `token` have to end the name or a part
# of it, which leaves `author`, `max_tokens` and `tokenizer` alone.
SENSITIVE_KEY_PATTERN = (
    r'password|passwd|passphrase|pwd(?:$|[_-])|secret|credential|'
    r'authorization|authentication|auth(?:$|[_-])|bearer|api[_-]?key|private[_-]?key|'
    r'cookie|token(?:$|[_-])')
_SECRET_KEY = re.compile(SENSITIVE_KEY_PATTERN, re.I)
_ABS_IN_TEXT = re.compile(r'(?<![\w])(?:/(?:Users|home|private|tmp|var|mnt|opt)/\S+|[A-Za-z]:[\\/][^\s]+)')
# A surrogate code point in a ``str`` (typically from ``surrogateescape``
# decoding of an undecodable file name) has no UTF-8 encoding; the report is
# UTF-8 bytes, so such code points become the replacement character.
_LONE_SURROGATE = re.compile('[\ud800-\udfff]')


TEXT_LIMIT = 256
# A clipped text must never read as the complete value (same marker in MATLAB).
TRUNCATION_MARKER = '...[truncated]'
# A URL of any scheme (`HTTPS://`, `git+https://`, `ssh://`, ...), in any case;
# the same literal as in `EvalReport.m`. Two letters at least: `C://x` is a path.
# The bound on the scheme keeps the search linear in a long text without a URL.
URL_PATTERN = r'[A-Za-z][A-Za-z0-9+.-]{1,63}://\S+'
_URL_IN_TEXT = re.compile(URL_PATTERN)


def _bounded(value, limit):
    if limit is None or len(value) <= limit:
        return value
    return value[:max(limit - len(TRUNCATION_MARKER), 0)] + TRUNCATION_MARKER


def _sanitized(value):
    value = ''.join(c if c >= ' ' else ' ' for c in value)
    return _LONE_SURROGATE.sub('\ufffd', value)


def bounded_text(value, limit=TEXT_LIMIT):
    """Bound known text without invoking object representations (``limit=None``: sanitize only)."""
    if isinstance(value, Enum):
        value = value.value
    if not isinstance(value, str):
        return None
    value = _sanitized(value)
    if _URL_IN_TEXT.match(value):
        return '[redacted_url]'
    value = _URL_IN_TEXT.sub('[redacted_url]', value)
    if os.path.isabs(value) or PureWindowsPath(value).is_absolute():
        return '[redacted_absolute_path]'
    value = _ABS_IN_TEXT.sub('[redacted_absolute_path]', value)
    return _bounded(value, limit)


def _credential_free_url(match):
    """Keep scheme, host, port and path; user information, query and fragment can carry secrets."""
    from urllib.parse import urlsplit
    text = match.group(0)
    try:
        parts = urlsplit(text)
        host = parts.hostname or ''
        port = f':{parts.port}' if parts.port is not None else ''
    except ValueError:
        return '[redacted_url]'
    # `\S+` runs on through a second URL glued to this one, whose user
    # information this split does not see. The first `://` is not in these parts.
    if not host or '://' in parts.netloc + parts.path:
        return '[redacted_url]'
    kept = f'{parts.scheme}://{host}{port}{parts.path}'
    if parts.username is not None or parts.password is not None or parts.query or parts.fragment:
        kept += ' [url parameters removed]'
    return kept


def credential_free_urls(value):
    """``value`` with the user information, query and fragment of its URLs removed."""
    return _URL_IN_TEXT.sub(_credential_free_url, value)


def credential_free_name(value):
    """``credential_free_urls(value)``, ending in a checksum of ``value`` if that changed it.

    Two names that differ only in what was removed stay two names, as the file
    names built from them must. The checksum is the CRC-32 of the UTF-8 text, as
    in `EvalReport.m`.
    """
    name = credential_free_urls(value)
    if name == value:
        return name
    return f"{name} {zlib.crc32(value.encode('utf-8', 'surrogatepass')):08x}"


def provenance_text(value, limit=TEXT_LIMIT):
    """Bounded provenance text, such as the source of a reference record.

    Unlike ``bounded_text``, a URL (for example a DOI link) is legitimate
    scientific provenance and keeps its scheme, host and path; user
    information, query and fragment are removed because they can carry
    credentials or access tokens. Absolute local paths are redacted and a
    clipped text carries ``TRUNCATION_MARKER``.
    """
    if isinstance(value, Enum):
        value = value.value
    if not isinstance(value, str):
        return None
    value = credential_free_urls(_sanitized(value))
    if os.path.isabs(value) or PureWindowsPath(value).is_absolute():
        return '[redacted_absolute_path]'
    value = _ABS_IN_TEXT.sub('[redacted_absolute_path]', value)
    return _bounded(value, limit)



def describe_callback(value):
    """Read function metadata only; do not inspect descriptors or source."""
    if isinstance(value, (types.FunctionType, types.BuiltinFunctionType,
                          types.MethodType)):
        return {'kind': 'callback', 'module': bounded_text(value.__module__),
                'name': bounded_text(value.__qualname__)}
    cls = value if isinstance(value, type) else type(value)
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
