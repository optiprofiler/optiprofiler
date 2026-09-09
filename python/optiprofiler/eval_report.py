"""Private, observation-only collector for ``optiprofiler.eval_report/1``.

This module does not call a solver, merit function, feature modifier, or other
user callback. Array values are observations supplied by the benchmark. Compact
facts go in the main JSON; a paired detail file keeps extrema-preserving history
summaries and complete numeric presentations produced by shared plot helpers.
"""

import hashlib
import json
import math
import mimetypes
import os
from pathlib import Path, PureWindowsPath
import re
import stat
import subprocess
import tempfile
import time
import types
import uuid
from datetime import datetime, timezone
from enum import Enum

import numpy as np


_SCHEMA = 'optiprofiler.eval_report/1'
_SCHEMA_NAMES = ('eval_report', 'plot_data')
_STAGES = ('numerical', 'scoring', 'persistence', 'rendering')
_STATUSES = {'not_requested', 'not_applicable', 'unknown', 'running',
             'completed', 'partial', 'failed'}
_MAX_BINS = 32
_MAX_DIAGNOSTICS = 128
_MAX_ARTIFACTS = 2048
_SECRET_KEY = re.compile(
    r'password|passwd|secret|credential|authorization|api[_-]?key|'
    r'(?:^|[_-])token(?:$|[_-])|access[_-]?token|refresh[_-]?token|'
    r'private[_-]?key|cookie', re.I)
_ABS_IN_TEXT = re.compile(r'(?<![\w])(?:/(?:Users|home|private|tmp|var|mnt|opt)/\S+|[A-Za-z]:[\\/][^\s]+)')


def _now():
    return datetime.now(timezone.utc).isoformat()


def schema_text(name):
    """Return the packaged JSON Schema source for ``'eval_report'`` or ``'plot_data'``.

    The schemas are package resources (``optiprofiler/schemas``), so an
    installed distribution, its installed test suite and the documentation
    build all read the same authoritative files. ``importlib.resources.files``
    exists from Python 3.9; on Python 3.8 the package directory is used.
    """
    if name not in _SCHEMA_NAMES:
        raise ValueError(f'Unknown EvalReport schema {name!r}; expected one of {_SCHEMA_NAMES}')
    filename = f'{name}.schema.json'
    try:
        from importlib.resources import files
    except ImportError:  # Python 3.8
        return (Path(__file__).with_name('schemas') / filename).read_text(encoding='utf-8')
    return (files('optiprofiler') / 'schemas' / filename).read_text(encoding='utf-8')


def load_schema(name):
    """Parse the packaged schema, rejecting non-standard NaN/Infinity tokens."""
    return json.loads(schema_text(name),
                      parse_constant=lambda token: (_ for _ in ()).throw(ValueError(token)))


def _hashing_capability():
    """How artifacts under the owned output tree can be hashed safely.

    ``'openat'``: POSIX descriptor-relative opens with O_NOFOLLOW (strongest;
    directory swaps during traversal are refused).
    ``'identity'``: Windows has no openat. Components are inspected with
    lstat (symlinks and junctions refused), the file is opened, and the
    opened handle must have the identity (volume, file index) recorded before
    the open, so the hashed bytes are provably the inspected file. This is the
    same guarantee level as the MATLAB collector; a reparse point inserted
    into the benchmark-owned tree between inspection and open is detected by
    the identity mismatch, not prevented.
    ``'unavailable'``: neither primitive; artifacts stay unverified with an
    explicit reason instead of a hash read through an unknown path.
    """
    if os.open in getattr(os, 'supports_dir_fd', set()):
        return 'openat'
    if os.name == 'nt':
        return 'identity'
    return 'unavailable'


def _text(value, limit=256):
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


def _callback(value):
    """Read function metadata only; do not inspect descriptors or source."""
    if isinstance(value, (types.FunctionType, types.BuiltinFunctionType,
                          types.MethodType)):
        return {'kind': 'callback', 'module': _text(value.__module__),
                'name': _text(value.__qualname__)}
    cls = type(value)
    return {'kind': 'callback', 'module': _text(cls.__module__),
            'name': _text(cls.__qualname__), 'instance_state': 'not_recorded'}


def _safe(value, depth=0, key=None):
    """A deliberately limited metadata encoder, not a general serializer."""
    if key and _SECRET_KEY.search(key):
        return {'value': None, 'reason': 'redacted_sensitive_option'}
    if depth > 10:
        return {'value': None, 'reason': 'metadata_depth_limit'}
    if value is None:
        return None
    if isinstance(value, Enum):
        return _safe(value.value, depth, key)
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
        return _text(value)
    if isinstance(value, (Path, PureWindowsPath)):
        return {'name': _text(value.name), 'path': None,
                'reason': 'machine_path_not_recorded'}
    if callable(value):
        return _callback(value)
    if type(value) is dict:
        answer = {}
        for i, (k, v) in enumerate(value.items()):
            if i >= 256:
                answer['_omission'] = {'reason': 'metadata_item_limit',
                                       'total_items': len(value)}
                break
            k = _text(k)
            if k is not None:
                answer[k] = _safe(v, depth + 1, k)
        return answer
    if type(value) in (list, tuple):
        if len(value) > 256:
            return {'values': [_safe(v, depth + 1) for v in value[:256]],
                    'total_items': len(value), 'reason': 'metadata_item_limit'}
        return [_safe(v, depth + 1) for v in value]
    if type(value) is np.ndarray:
        if value.size <= 256:
            return _safe(value.tolist(), depth + 1, key)
        return {'shape': list(value.shape), 'value': None,
                'reason': 'raw_array_not_embedded'}
    return {'value': None, 'reason': 'unsupported_metadata_type',
            'type': _text(type(value).__name__)}


def _history_bins(history, count):
    """Preserve exact per-bin extrema, not a uniform sample that misses spikes."""
    vector = _vector(history)
    if vector is None or count is None:
        return {'status': 'unavailable', 'reason': 'history_or_evaluation_count_unavailable',
                'count': 0, 'total_evaluations': count, 'representation': None, 'bins': []}
    count_retained = min(max(count, 0), len(vector))
    reason = 'history_shorter_than_evaluation_count' if count > len(vector) else ('no_evaluations' if not count else None)
    if count_retained <= 64:
        return {'status': 'available', 'reason': reason,
                'count': count_retained, 'total_evaluations': count,
                'representation': 'exact_samples',
                'values': _score_values(vector[:count_retained]), 'bins': []}
    n_bins = min(count_retained, _MAX_BINS)
    bins = []
    for i in range(n_bins):
        # Exact integer edges i*n//32. MATLAB uses floor(linspace(0, n, 33)):
        # n/32 is an exact binary fraction, so both languages partition a
        # history of the same length identically (validated by both suites).
        start, end = i * count_retained // n_bins, (i + 1) * count_retained // n_bins
        values = vector[start:end]
        finite_indices = np.flatnonzero(np.isfinite(values))
        extrema = {}
        for label, operation in (('finite_min', np.argmin), ('finite_max', np.argmax)):
            if finite_indices.size:
                index = int(finite_indices[operation(values[finite_indices])])
                extrema[label] = {'value': _safe(values[index]), 'evaluation_index': start + index + 1}
            else:
                extrema[label] = None
        bins.append({'start_index': start + 1, 'end_index': end,
                     'first': _safe(values[0]), 'last': _safe(values[-1]), **extrema,
                     'nonfinite': {'nan': int(np.count_nonzero(np.isnan(values))),
                                   'positive_infinity': int(np.count_nonzero(np.isposinf(values))),
                                   'negative_infinity': int(np.count_nonzero(np.isneginf(values)))}})
    return {'status': 'available', 'reason': reason, 'representation': 'bin_extrema',
            'count': count_retained, 'total_evaluations': count, 'bins': bins}


def _score_values(value):
    """Preserve returned score tensors, unlike bounded option metadata."""
    if type(value) is np.ndarray and value.dtype.kind in 'biuf':
        return _score_values(value.tolist())
    if type(value) in (list, tuple):
        return [_score_values(item) for item in value]
    return _safe(value)


def _numeric_payload(value):
    """Encode a known plot-data structure completely, never as option metadata."""
    if type(value) is dict:
        return {key: _numeric_payload(item) for key, item in value.items()}
    if type(value) in (list, tuple):
        return [_numeric_payload(item) for item in value]
    return _score_values(value)


def _relative_path(path, base):
    try:
        return Path(os.path.relpath(path, base)).as_posix()
    except ValueError:
        # Windows cannot express a path on another drive relatively. Never
        # fall back to exposing an absolute machine path in the report.
        return None


def _at(value, *indices):
    if value is None:
        return None
    try:
        for index in indices:
            if np.ndim(value) == 0:
                return value
            value = value[index]
        return value
    except (IndexError, KeyError, TypeError):
        return None


def _integer(value):
    if isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_)):
        return int(value)
    if isinstance(value, (float, np.floating)) and np.isfinite(value) and value == int(value):
        return int(value)
    return None


def _vector(value):
    if type(value) is np.ndarray and value.dtype.kind in 'biuf':
        return value.reshape(-1)
    if type(value) in (list, tuple) and all(
            isinstance(v, (bool, int, float, np.number)) for v in value):
        return np.asarray(value).reshape(-1)
    return None


def _metric(history, count, output, initial):
    """One objective/constraint/merit record; shared with MATLAB's metric().

    Keys that would only carry ``null`` are omitted: an absent
    ``availability_reason`` means every observation was available and an
    absent ``first_invalid_evaluation_index`` means no invalid evaluation was
    observed (``invalid_evaluations`` is then the integer 0).
    """
    vector = _vector(history)
    actual = None if vector is None or count is None else vector[:max(count, 0)]
    best = None
    best_index = None
    invalid_index = None
    invalid = None
    reason = 'history_or_evaluation_count_unavailable'
    if actual is not None:
        reason = None
        valid = actual[~np.isnan(actual)]
        if valid.size:
            best = _safe(np.min(valid))
            best_index = int(np.flatnonzero(actual == np.min(valid))[0]) + 1
        elif actual.size:
            best = _safe(float('nan'))
        invalid = {
            'nan': int(np.count_nonzero(np.isnan(actual))),
            'positive_infinity': int(np.count_nonzero(np.isposinf(actual))),
            'negative_infinity': int(np.count_nonzero(np.isneginf(actual))),
            'observed_evaluations': int(actual.size),
        }
        bad = np.flatnonzero(~np.isfinite(actual))
        invalid_index = int(bad[0]) + 1 if bad.size else None
        if count > len(vector):
            reason = 'history_shorter_than_evaluation_count'
        elif not count:
            reason = 'no_evaluations'
        if bad.size == 0 and actual.size == count:
            # Nothing is lost: observed_evaluations equals the run's
            # evaluations and every categorized count is zero.
            invalid = 0
    record = {'output': _safe(output), 'initial': _safe(initial), 'best': best,
              'best_evaluation_index': best_index, 'invalid_evaluations': invalid}
    if invalid_index is not None:
        record['first_invalid_evaluation_index'] = invalid_index
    if reason is None and (output is None or initial is None):
        reason = 'output_or_initial_unavailable'
    if reason is not None:
        record['availability_reason'] = reason
    return record


# IsReparseTagNameSurrogate: set on the reparse tags that redirect a name
# (symbolic links, junctions/mount points), clear on those that keep it
# (cloud-file placeholders, compressed or deduplicated files).
_NAME_SURROGATE = 0x20000000


def _reparse_point(status):
    """A symlink, or on Windows a junction or other name-redirecting reparse point."""
    if stat.S_ISLNK(status.st_mode):
        return True
    return bool(getattr(status, 'st_reparse_tag', 0) & _NAME_SURROGATE)


def _link_like(path):
    """True for a link (see ``_reparse_point``) or an entry that cannot be inspected."""
    try:
        return _reparse_point(os.lstat(str(path)))
    except OSError:
        return True


def _digest(path, directory=None):
    flags = os.O_RDONLY | getattr(os, 'O_NOFOLLOW', 0) | getattr(os, 'O_BINARY', 0)
    handles = []
    inspected = None
    ancestors = []
    if directory is None:
        fd = os.open(str(path), flags)
    else:
        capability = _hashing_capability()
        parts = path.relative_to(directory).parts
        if capability == 'openat':
            # Resolve every component beneath the owned directory with openat,
            # refusing symlinks even if a directory changes during traversal.
            directory_flags = flags | getattr(os, 'O_DIRECTORY', 0)
            parent = os.open(str(directory), directory_flags)
            handles.append(parent)
            try:
                for name in parts[:-1]:
                    parent = os.open(name, directory_flags, dir_fd=parent)
                    handles.append(parent)
                fd = os.open(parts[-1], flags, dir_fd=parent)
            except BaseException:
                for handle in reversed(handles):
                    os.close(handle)
                raise
        elif capability == 'identity':
            # Inspect, open, then prove the opened file is the inspected one.
            current = directory
            for name in parts[:-1]:
                current = current / name
                status = os.lstat(str(current))
                if _reparse_point(status) or not stat.S_ISDIR(status.st_mode):
                    raise OSError('artifact path component is not a plain directory')
                ancestors.append(current)
            inspected = os.lstat(str(path))
            if _reparse_point(inspected) or not stat.S_ISREG(inspected.st_mode):
                raise OSError('artifact is not a plain regular file')
            fd = os.open(str(path), flags)
        else:
            # Without openat or a verifiable identity there is no
            # race-resistant ownership check here. A missing provenance hash
            # is preferable to reading through a replaced/symlinked directory,
            # and must not abort scores.
            raise OSError('secure_artifact_hashing_unavailable')
    try:
        before = os.fstat(fd)
        if not stat.S_ISREG(before.st_mode):
            raise OSError('not a regular file')
        if inspected is not None and (before.st_dev, before.st_ino) != (inspected.st_dev, inspected.st_ino):
            raise OSError('artifact identity changed between inspection and open')
        digest = hashlib.sha256()
        with os.fdopen(fd, 'rb', closefd=False) as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(chunk)
        after = os.fstat(fd)
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise OSError('file changed during hashing')
        for ancestor in ancestors:
            # A reparse point that appeared meanwhile is reported, not trusted.
            if _reparse_point(os.lstat(str(ancestor))):
                raise OSError('artifact path component changed during hashing')
        return {'bytes': before.st_size, 'sha256': digest.hexdigest()}
    finally:
        os.close(fd)
        for handle in reversed(handles):
            os.close(handle)


def _producer():
    source = Path(__file__).resolve()
    result = {'language': 'python', 'version': None, 'revision': None,
              'dirty': None, 'exact_revision': False,
              'scope': 'current_invocation_not_original_archive_producer',
              'version_kind': 'local_package_source_declaration',
              'revision_kind': 'local_checkout_head_not_full_current_source_identity',
              'identity_note': 'Version and hashes do not establish reproducibility or comparability.'}
    try:
        with source.with_name('__init__.py').open('r', encoding='utf-8') as stream:
            match = re.search(r"^__version__\s*=\s*['\"]([^'\"]+)['\"]", stream.read(65536), re.M)
        if match:
            result['version'] = match.group(1)
    except (OSError, UnicodeError):
        pass
    root = source.parent.parent.parent
    if not (root / '.git').exists():
        result['revision_reason'] = 'source_checkout_unavailable'
        return result
    try:
        common = {'cwd': str(root), 'stderr': subprocess.DEVNULL,
                  'timeout': 2, 'check': False}
        revision = subprocess.run(['git', 'rev-parse', '--verify', 'HEAD'],
                                  stdout=subprocess.PIPE, **common)
        value = revision.stdout.decode('ascii', errors='ignore').strip()
        if revision.returncode == 0 and re.fullmatch(r'[0-9a-f]{40,64}', value):
            result['revision'] = value
        paths = ['pyproject.toml', 'python/optiprofiler']
        tracked = subprocess.run(['git', 'diff', '--quiet', 'HEAD', '--'] + paths,
                                 stdout=subprocess.DEVNULL, **common)
        current = subprocess.run(['git', 'ls-files', '--error-unmatch', '--',
                                  'python/optiprofiler/eval_report.py'],
                                 stdout=subprocess.DEVNULL, **common)
        if tracked.returncode == 1 or current.returncode == 1:
            result['dirty'] = True
            result['revision_reason'] = 'checkout_contains_uncommitted_source_changes'
        else:
            # Checking one file and tracked changes does not certify absence
            # of every untracked source; do not claim a clean exact revision.
            result['revision_reason'] = 'untracked_source_completeness_not_certified'
    except (OSError, subprocess.SubprocessError):
        result['revision_reason'] = 'local_revision_check_unavailable'
    return result


class EvalReport:
    """Collect a report owned exclusively by one benchmark invocation."""

    def __init__(self, path, request):
        self.path = Path(path).expanduser().absolute()
        self.plot_path = self.path.with_name(self.path.stem + '.plot_data.json')
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._parent_identity = self._identity(self.path.parent)
        self._started = time.monotonic()
        self._finished = False
        self._problems = {}
        self._metadata = {}
        self._options = {}
        self._output_dir = None
        self._output_identity = None
        self._initial_artifacts = set()
        self._source_path = None
        self._selection_seen = False
        self._dropped_diagnostics = 0
        self._diagnostic_keys = set()
        self._histories = {}
        self._plots = {}
        self._detail_dirty = True
        request = request if type(request) is dict else {}
        load_request = request.get('load')
        operation = 'load' if request.get('operation') == 'load' or (
            isinstance(load_request, str) and bool(load_request)) else 'benchmark'
        self.document = {
            'schema': _SCHEMA, 'evaluation_id': str(uuid.uuid4()),
            'operation': operation, 'status': 'running',
            'producer': _producer(), 'configuration': {'request': _safe(request)},
            'stages': {stage: {'status': 'unknown'} for stage in _STAGES},
            'coverage': {'selected': 0, 'loaded': 0, 'completed': 0,
                         'load_failed': 0, 'scope': 'live_selection',
                         'original_selection_known': True},
            'problems': [], 'scores': None,
            'profiles': {'plot_refs': [],
                         'convergence': [], 'work_summary': None,
                         'work_summary_reason': 'not_supplied'},
            'artifact_root': None, 'artifacts': [], 'diagnostics': [], 'source': None,
            'plot_data': None,
            # Shared vocabulary with MATLAB: the same keys are emitted by
            # EvalReport.m and pinned by eval_report.schema.json. Values are
            # prose describing THIS producer; language-specific conventions
            # (bands, seeds) are stated, not equalized.
            'semantics': {
                'index_base': 1,
                'metric_best': 'componentwise_minimum_ignoring_nan;first_tie_index;not_necessarily_a_jointly_attained_point',
                'budget': 'problems[].budget applies to every run unless runs[].budget overrides it; budget_reached is a per-run comparison and reaching the cap does not identify the termination cause.',
                'convergence': 'Never inferred from solver return values; there is no per-run convergence field. target_work in the companion observes the existing scoring predicate.',
                'run_defaults': 'evaluations, budget_reached, abnormal_termination, output_fallback, execution, oracle_seed and elapsed_seconds are per-run facts. An absent *_reason key means the observation was available; an absent first_invalid_evaluation_index means no invalid evaluation was observed and invalid_evaluations is then the integer 0.',
                'configuration': 'configuration.request lists user-supplied options; configuration.effective lists the resolved options of this invocation. In a load operation they describe reanalysis/rendering, not the archived execution.',
                'paths': {'artifacts': 'relative_to_artifact_root', 'artifact_root': 'relative_to_main_report_parent',
                          'plot_data': 'relative_to_main_report_parent', 'source': 'relative_to_main_report_parent'},
                'privacy': 'controller_private;not_an_allowlisted_agent_prompt;consumers_build_an_allowlisted_feedback_view',
            },
            'report_files': None,
            'timing': {'started_at': _now(), 'finished_at': None,
                       'elapsed_seconds': None},
        }
        self.plot_document = {'schema': 'optiprofiler.plot_data/1',
                              'evaluation_id': self.document['evaluation_id'],
                              'semantics': {'index_base': 1,
                                  'history_bins': 'exact_samples_up_to_64;otherwise_lossy_max_32_contiguous_bins_with_exact_endpoints_finite_extrema_first_tie_indices_and_nonfinite_counts;internal_order_not_retained_in_bins;excludes_padding;not_a_scoring_input',
                                  'history_plots': 'display_copy_only;across_run_band_then_shift_then_cummin_then_block_aggregation;block_aggregation_applies_when_padded_length_exceeds_1002_and_keeps_about_1000_points;array_position_is_not_an_evaluation_index_use_evaluation_indices',
                                  'profile_plots': 'full_step_vertices_and_across_run_bands;band_clamped_to_0_1;unreached_work_shown_at_failure_placeholder_not_counted_as_hit;log_ratio_zero_height_ties_retained_in_numeric_data_but_not_drawn',
                                  'error_bands': 'python_population_std_ddof_0_for_history_and_profile_bands',
                                  'target_work': 'existing_profile_work_arrays;history_first_actual_hit;output_total_evaluations_at_passing_returned_output;nan_is_nonhit_not_solver_diagnosis',
                                  'privacy': 'controller_private;not_an_allowlisted_agent_prompt'},
                              'histories': [], 'plots': [], 'target_work': [], 'diagnostics': []}
        # Preflight both names before reserving either. O_EXCL also closes the
        # race with another invocation; rollback only a still-owned empty file.
        for target in (self.path, self.plot_path):
            if target.exists() or target.is_symlink():
                raise FileExistsError(str(target))
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, 'O_NOFOLLOW', 0)
        reserved = []
        try:
            for target in (self.path, self.plot_path):
                fd = os.open(str(target), flags, 0o600)
                identity = self._stat_identity(os.fstat(fd))
                os.close(fd)
                reserved.append((target, identity))
        except BaseException:
            for target, identity in reserved:
                if not target.is_symlink() and self._identity(target) == identity:
                    target.unlink()
            raise
        self._owned_identity, self._plot_identity = (entry[1] for entry in reserved)
        # Both targets are created with mode 0600 (O_EXCL here, mkstemp for
        # every republish). This is only a best-effort owner-only policy on
        # POSIX filesystems; it is not enforced on Windows or on filesystems
        # without POSIX modes, so the directory itself must be private.
        applied = None
        if os.name == 'posix':
            try:
                applied = all(stat.S_IMODE(target.lstat().st_mode) == 0o600
                              for target, _ in reserved)
            except OSError:
                applied = None
        self.document['report_files'] = {
            'permission_policy': 'owner_read_write_only_best_effort',
            'permissions_applied': applied,
            'platform_note': 'posix_mode_0600_via_exclusive_create_and_mkstemp;not_enforced_on_windows;directory_privacy_is_the_caller_responsibility'}
        self._write()

    @staticmethod
    def _stat_identity(value):
        return value.st_dev, value.st_ino

    @classmethod
    def _identity(cls, path):
        return cls._stat_identity(path.lstat())

    def _write(self):
        self._coverage()
        self.document['problems'] = list(self._problems.values())
        self.document['timing']['elapsed_seconds'] = max(0.0, time.monotonic() - self._started)
        self.plot_document['histories'] = list(self._histories.values())
        self.plot_document['plots'] = list(self._plots.values())
        # Publish detail first. A stale main/detail pair is detectable by both
        # evaluation_id and SHA256; never publish a manifest for unwritten data.
        if self._detail_dirty:
            self._plot_identity = self._write_owned(self.plot_path, self.plot_document, self._plot_identity)
            self.document['plot_data'] = {
                'schema': self.plot_document['schema'], 'path': self.plot_path.name,
                'status': 'completed', **_digest(self.plot_path),
                'history_count': len(self._histories), 'plot_count': len(self._plots)}
            self._detail_dirty = False
        if self.plot_path.is_symlink() or self._identity(self.plot_path) != self._plot_identity:
            raise FileExistsError('EvalReport companion ownership changed')
        self._owned_identity = self._write_owned(self.path, self.document, self._owned_identity)

    def _write_owned(self, target, document, owned_identity):
        if (target.is_symlink() or self._identity(target) != owned_identity
                or self._identity(target.parent) != self._parent_identity):
            raise FileExistsError('EvalReport target ownership changed')
        fd, temp = tempfile.mkstemp(prefix='.' + target.name + '.',
                                    suffix='.tmp', dir=str(target.parent))
        try:
            with os.fdopen(fd, 'w', encoding='utf-8') as stream:
                json.dump(document, stream, allow_nan=False, ensure_ascii=False,
                          separators=(',', ':'))
                stream.write('\n')
                stream.flush()
                os.fsync(stream.fileno())
                identity = self._stat_identity(os.fstat(stream.fileno()))
            if (target.is_symlink() or self._identity(target) != owned_identity
                    or self._identity(target.parent) != self._parent_identity):
                raise FileExistsError('EvalReport target ownership changed')
            os.replace(temp, target)
            return identity
        finally:
            try:
                os.unlink(temp)
            except FileNotFoundError:
                pass

    def configure(self, problem_options, profile_options, feature, output_dir=None):
        self._options = dict(profile_options) if type(profile_options) is dict else {}
        # Feature is an internal validated object. Read stored data without
        # invoking properties, __repr__, or custom modifier callbacks.
        try:
            state = object.__getattribute__(feature, '__dict__')
        except (AttributeError, TypeError):
            state = {}
        feature_data = {'name': _safe(state.get('_name')),
                        'options': _safe(state.get('_options', {}))}
        # request = what the caller supplied; effective = the resolved options
        # of this invocation. Stated once here, never repeated per run.
        self.document['configuration']['effective'] = {
            'problem_options': _safe(problem_options),
            'profile_options': _safe(profile_options), 'feature': feature_data}
        if self.document['operation'] == 'load':
            self.document['configuration'].update({
                'scope': 'current_load_selection_reanalysis_and_rendering',
                'solver_execution_requested': False,
                'original_execution_configuration': None,
                'original_execution_configuration_reason': 'not_fully_retained_by_archive',
            })
            feature_data['scope'] = 'current_load_context_not_original_execution_feature'
        else:
            self.document['configuration'].update({
                'scope': 'current_benchmark_execution', 'solver_execution_requested': True})
            feature_data['scope'] = 'current_execution_feature'
        score_only = bool(self._options.get('score_only', False))
        for stage in ('persistence', 'rendering'):
            if self.document['stages'][stage]['status'] == 'unknown':
                self.document['stages'][stage] = {
                    'status': 'not_requested' if score_only else 'unknown',
                    'reason': 'score_only' if score_only else 'requested_not_yet_observed'}
        if output_dir is not None:
            path = Path(output_dir).absolute()
            if path != self._output_dir:
                if path.is_symlink():
                    self.add_diagnostic('artifact_directory_symlink_skipped', 'persistence')
                elif path.is_dir():
                    self._output_dir = path
                    self.document['artifact_root'] = _relative_path(path, self.path.parent)
                    self._output_identity = self._identity(path)
                    # Only files created after the benchmark establishes its
                    # owned output directory are eligible as new artifacts.
                    self._initial_artifacts = {str(p) for p in self._artifact_paths(path)}
        self._write()

    def _problem(self, plib, name, role):
        raw = tuple(value.value if isinstance(value, Enum) else value
                    if isinstance(value, str) else None for value in (plib, name, role))
        library, name, role = (_text(value) for value in raw)
        key = raw
        if key not in self._problems:
            identity = json.dumps([library, name, role], ensure_ascii=False, separators=(',', ':'))
            if raw != (library, name, role):
                private_identity = json.dumps(raw, ensure_ascii=False, separators=(',', ':'))
                identity = 'sha256:' + hashlib.sha256(private_identity.encode('utf-8')).hexdigest()
            self._problems[key] = {
                'id': identity,
                'library': library, 'name': name, 'role': role,
                'dimension': None, 'type': None,
                'selection_status': 'unknown', 'load_status': 'pending',
                'status': 'pending', 'budget': None, 'runs': [],
            }
        return self._problems[key]

    def selection(self, plib, names, provider=None, role='primary'):
        self._selection_seen = self._selection_seen or role == 'primary'
        for name in names:
            problem = self._problem(plib, name, role)
            problem['selection_status'] = 'selected'
            if provider is not None:
                if type(provider).__name__ == 'ProblemLibraryRef':
                    state = object.__getattribute__(provider, '__dict__')
                    provider = {key: state.get(key) for key in
                                ('name', 'source', 'locator', 'distribution')}
                problem['provider'] = _safe(provider)
        self._write()

    def add_diagnostic(self, code, stage, **scope):
        code, stage = _text(code, 80), _text(stage, 40)
        if code == 'problem_load_failed':
            problem = self._problem(scope.get('library'), scope.get('problem'),
                                    scope.get('role', 'primary'))
            problem.update(load_status='failed', status='failed')
        if stage == 'rendering' and code and ('failed' in code or 'failure' in code):
            self.document['stages']['rendering'] = {'status': 'failed', 'reason': code}
        if stage == 'persistence' and code in (
                'artifact_hash_unavailable', 'artifact_directory_ownership_changed'):
            if self.document['stages']['persistence']['status'] != 'failed':
                self.document['stages']['persistence'] = {'status': 'partial', 'reason': code}
        entry = {'code': code, 'stage': stage, 'scope': _safe(scope)}
        identity = json.dumps(entry, allow_nan=False, sort_keys=True, ensure_ascii=False)
        if identity in self._diagnostic_keys:
            return
        if len(self.document['diagnostics']) < _MAX_DIAGNOSTICS:
            self._diagnostic_keys.add(identity)
            self.document['diagnostics'].append(entry)
        else:
            self._dropped_diagnostics += 1

    def set_stage(self, stage, status, reason=None):
        if stage not in _STAGES or status not in _STATUSES:
            raise ValueError('Invalid EvalReport stage or status')
        if (stage in ('rendering', 'persistence') and status in ('running', 'completed')
                and self.document['stages'][stage]['status'] in ('failed', 'partial')):
            # A later successful export cannot erase an earlier missing
            # requested artifact in this same invocation.
            self._write()
            return
        self.document['stages'][stage] = {'status': status}
        if reason is not None:
            self.document['stages'][stage]['reason'] = _text(reason)
        self._write()

    def add_problem(self, result, plib, role='primary'):
        if type(result) is not dict:
            return
        self._detail_dirty = True
        name = result.get('problem_name')
        problem = self._problem(plib, name, role)
        key = problem['id']
        if type(result.get('eval_report_metadata')) is dict:
            self._metadata[key] = result['eval_report_metadata']
        metadata = self._metadata.get(key, {})
        problem.update(dimension=_safe(result.get('problem_dim')),
                       type=_safe(result.get('problem_type')),
                       load_status='loaded', status='completed')
        if problem['selection_status'] == 'unknown':
            problem['selection_status'] = ('retained' if self.document['operation'] == 'load'
                                           else 'observed_without_selection')
        n_evals = result.get('n_eval')
        if type(n_evals) is not np.ndarray:
            n_evals = np.asarray(n_evals) if type(n_evals) in (list, tuple) else None
        if n_evals is None or n_evals.ndim != 2:
            problem['status'] = 'unknown'
            problem['reason'] = 'run_evaluation_counts_unavailable'
            return
        budget = _integer(metadata.get('budget'))
        dim = _integer(result.get('problem_dim'))
        factor = self._options.get('max_eval_factor')
        if (budget is None and self.document['operation'] != 'load' and dim is not None
                and isinstance(factor, (int, float, np.number)) and np.isfinite(factor)):
            budget = int(math.ceil(float(factor) * dim))
        # The cap is a per-problem fact (same for every solver/run of this
        # invocation); each run keeps its own budget_reached comparison.
        if budget is None:
            problem['budget'] = {'evaluations': None,
                                 'reason': 'original_execution_budget_not_retained'
                                 if self.document['operation'] == 'load' else 'budget_rule_inputs_unavailable'}
        else:
            problem['budget'] = {'evaluations': budget, 'rule': 'ceil(max_eval_factor*dimension)'}
        runs = []
        for solver_index in range(n_evals.shape[0]):
            for run_index in range(n_evals.shape[1]):
                count = _integer(n_evals[solver_index, run_index])
                histories = {label: _at(result.get(field), solver_index, run_index)
                             for label, field in (('objective', 'fun_history'),
                                                  ('constraint', 'maxcv_history'),
                                                  ('merit', 'merit_history'))}
                real = _integer(_at(metadata.get('real_n_runs'), solver_index))
                execution = {'kind': 'unknown', 'source_run_index': None,
                             'reason': 'execution_metadata_not_retained'}
                if real is not None:
                    execution = {'kind': 'actual' if run_index < real else 'repeated',
                                 'source_run_index': None if run_index < real else 1}
                run = {
                    'solver_index': solver_index + 1, 'run_index': run_index + 1,
                    'evaluations': count,
                    'budget_reached': None if count is None or budget is None else count >= budget,
                    'abnormal_termination': _safe(_at(result.get('solver_abnormal_termination'), solver_index, run_index)),
                    'output_fallback': _safe(_at(result.get('solver_output_fallback'), solver_index, run_index)),
                    'execution': execution,
                    'oracle_seed': _safe(_at(metadata.get('oracle_seeds'),
                                             0 if execution['kind'] == 'repeated' else run_index)),
                    'elapsed_seconds': _safe(_at(result.get('computation_time'), solver_index, run_index)),
                }
                # Reasons are recorded only when an observation is missing;
                # see semantics.run_defaults. No convergence field exists.
                if 'solver_abnormal_termination' not in result or 'solver_output_fallback' not in result:
                    run['termination_metadata_reason'] = 'solver_termination_metadata_not_retained'
                if 'oracle_seeds' not in metadata:
                    run['oracle_seed_reason'] = 'execution_metadata_not_retained'
                elif execution['kind'] == 'repeated':
                    run['oracle_seed_reason'] = 'copied_from_source_run'
                for label, prefix in (('objective', 'fun'), ('constraint', 'maxcv'), ('merit', 'merit')):
                    run[label] = _metric(histories[label], count,
                                         _at(result.get(prefix + '_out'), solver_index, run_index),
                                         _at(result.get(prefix + '_init'), run_index))
                ordinal = next(i for i, item in enumerate(self._problems.values(), 1) if item is problem)
                history_id = f'history-{ordinal}-{solver_index + 1}-{run_index + 1}'
                self._histories[history_id] = {
                    'id': history_id, 'problem_id': key,
                    'solver_index': solver_index + 1, 'run_index': run_index + 1,
                    'channels': {label: _history_bins(value, count) for label, value in histories.items()}}
                run['history_ref'] = history_id
                runs.append(run)
        problem['runs'] = runs
        self._add_history_plots(result, problem)
        if metadata.get('render_status') in _STATUSES:
            problem['rendering'] = {'status': metadata['render_status']}
            if (metadata['render_status'] in ('failed', 'partial')
                    and self.document['stages']['rendering']['status'] != 'failed'):
                self.document['stages']['rendering'] = {'status': metadata['render_status'],
                                                       'reason': 'history_render_failed'}
            elif (metadata['render_status'] == 'completed'
                    and self.document['stages']['rendering']['status'] == 'unknown'):
                # A direct-problem run renders only history plots; observe
                # that success like MATLAB does. Profile export, when it
                # happens later, re-enters running and can still fail.
                self.document['stages']['rendering'] = {'status': 'completed'}
        for diagnostic in metadata.get('diagnostics', [])[:_MAX_DIAGNOSTICS]:
            if type(diagnostic) is dict:
                self.add_diagnostic(diagnostic.get('code'), diagnostic.get('stage'),
                                    library=plib, problem=name, role=role,
                                    exception_type=diagnostic.get('exception_type'))

    def _add_history_plots(self, result, problem):
        """Prepare display facts only: no axes, user callbacks, or raw mutation."""
        from .plotting import prepare_history_panels
        dim, counts = _integer(result.get('problem_dim')), result.get('n_eval')
        if dim is None or type(counts) is not np.ndarray:
            return
        metadata = self._metadata.get(problem['id'], {})
        panels = metadata.get('rendered_history_plots')
        scope = 'rendering_inputs' if panels is not None else 'retained_scoring_observations'
        try:
            if panels is None:
                fields = (('objective', 'fun'), ('constraint', 'maxcv'), ('merit', 'merit'))
                panels = prepare_history_panels(
                    {label: result.get(prefix + '_history') for label, prefix in fields},
                    {label: result.get(prefix + '_init') for label, prefix in fields},
                    problem['type'], dim, counts, self._options)
            self._store_history_plots(problem, panels, scope)
        except Exception as exc:
            self.add_diagnostic('history_plot_data_unavailable', 'reporting',
                                problem=problem['name'], exception_type=type(exc).__name__)

    def _store_history_plots(self, problem, panels, scope):
        ordinal = next(i for i, item in enumerate(self._problems.values(), 1) if item is problem)
        plot_refs = []
        for panel in panels:
            plot_id = f"history-plot-{ordinal}-{panel['channel']}-{panel['mode']}"
            self._plots[plot_id] = {'id': plot_id, 'problem_id': problem['id'],
                'fidelity': 'exact_rendered_data', 'status': 'numeric_prepared',
                'observation_scope': scope, **_numeric_payload(panel)}
            plot_refs.append(plot_id)
        problem['plot_refs'] = plot_refs

    def observe_rendered_history(self, plib, name, panels):
        problem = self._problem(plib, name, 'primary')
        self._metadata.setdefault(problem['id'], {})['rendered_history_plots'] = panels
        self._detail_dirty = True
        self._store_history_plots(problem, panels, 'rendering_inputs')

    def add_profile_plot(self, kind, presentation, tolerance_index, channel):
        self._detail_dirty = True
        plot_id = f'profile-{tolerance_index}-{channel}-{kind}'
        record = {'id': plot_id, 'kind': kind, 'history_or_output': channel,
                  'tolerance_index': tolerance_index,
                  'fidelity': 'exact_rendered_data', 'status': 'numeric_prepared',
                  'observation_scope': 'scoring_profile_work',
                  'target_work_ref': f'target-work-{tolerance_index}'}
        if kind != 'log_ratio':
            # Bars carry no across-run band, so band conventions are stated
            # only on performance/data records (same as MATLAB).
            record.update(errorbar_type=_safe(self._options.get('errorbar_type')), std_ddof=0)
        # Plot arrays are complete numerical data, not bounded option metadata.
        record.update(_numeric_payload(presentation))
        self._plots[plot_id] = record
        if plot_id not in self.document['profiles']['plot_refs']:
            self.document['profiles']['plot_refs'].append(plot_id)

    def add_results(self, results_plibs, operation='benchmark'):
        if operation == 'load':
            self.document['operation'] = 'load'
        mapping = {
            'problem_dim': 'problem_dims', 'problem_type': 'problem_types',
            'n_eval': 'n_evals', 'fun_history': 'fun_histories',
            'maxcv_history': 'maxcv_histories', 'merit_history': 'merit_histories',
            'fun_out': 'fun_outs', 'maxcv_out': 'maxcv_outs', 'merit_out': 'merit_outs',
            'fun_init': 'fun_inits', 'maxcv_init': 'maxcv_inits', 'merit_init': 'merit_inits',
            'computation_time': 'computation_times',
            'solver_abnormal_termination': 'solver_abnormal_terminations',
            'solver_output_fallback': 'solver_output_fallbacks',
        }

        def add_library(result, role):
            if type(result) is not dict:
                return
            if self.document['operation'] == 'load':
                retained = self.document['configuration'].setdefault('retained_result_metadata', [])
                metadata = {'library': _safe(result.get('plib')), 'role': role,
                            'feature_stamp': _safe(result.get('feature_stamp')),
                            'solver_names': _safe(result.get('solver_names')),
                            'library_options': _safe(result.get('plib_options')),
                            'scope': 'retained_result_after_load_filtering_not_complete_original_configuration'}
                if metadata not in retained:
                    retained.append(metadata)
            for index, name in enumerate(result.get('problem_names', [])):
                single = {'problem_name': name}
                for target, source in mapping.items():
                    if source in result:
                        single[target] = _at(result[source], index)
                self.add_problem(single, result.get('plib'), role)
            if role == 'primary' and type(result.get('results_plib_plain')) is dict:
                add_library(result['results_plib_plain'], 'plain_reference')

        for result in results_plibs or []:
            add_library(result, 'primary')
        self._write()

    def _coverage(self):
        primary = [problem for problem in self._problems.values()
                   if problem['role'] == 'primary']
        archived = self.document['operation'] == 'load'
        self.document['coverage'] = {
            'selected': len(primary) if archived else sum(p['selection_status'] == 'selected' for p in primary),
            'loaded': sum(p['load_status'] == 'loaded' for p in primary),
            'completed': sum(p['status'] == 'completed' for p in primary),
            'load_failed': None if archived else sum(p['load_status'] == 'failed' for p in primary),
            'scope': 'retained_archive' if archived else 'live_selection',
            'original_selection_known': False if archived else self._selection_seen,
            'reason': 'original_selection_and_load_failures_not_retained' if archived else
                      (None if self._selection_seen else 'selection_not_observed'),
        }

    def set_source(self, path):
        source = Path(path).expanduser().absolute()
        self._source_path = source
        self.document['operation'] = 'load'
        record = {'kind': 'archive', 'name': _text(source.name),
                  'path': _relative_path(source, self.path.parent),
                  'status': 'unknown', 'coverage_note': 'Only retained/replayed coverage is certified.'}
        if record['path'] is None:
            record['path_reason'] = 'cross_volume_relative_path_unavailable'
        try:
            if source.is_symlink():
                raise OSError('source is symlink')
            record.update(_digest(source), status='captured_before_filtering')
        except OSError:
            record.update(bytes=None, sha256=None, status='unavailable', reason='archive_hash_unavailable')
            self.add_diagnostic('source_hash_unavailable', 'persistence')
        self.document['source'] = record
        self._write()

    def set_profiles(self, curves, solver_scores, profile_scores, solver_names, single=False):
        n_types = profile_scores.shape[-1] if type(profile_scores) is np.ndarray and profile_scores.ndim == 4 else 0
        self.document['scores'] = {
            'solver_names': _score_values(solver_names), 'solver_scores': _score_values(solver_scores),
            'profile_scores': _score_values(profile_scores),
            'profile_axes': ['solver', 'tolerance', 'history_or_output', 'profile_type'],
            'axis_values': {'history_or_output': ['history', 'output'],
                            'profile_type': ['performance', 'data', 'log_ratio'][:n_types],
                            'tolerance': [item['tolerance'] for item in self.document['profiles']['convergence']]},
            'semantics': 'single_problem_relative_decrease_averaged_over_runs' if single else 'cohort_relative',
            'direction': 'higher_is_better_for_default_scoring; custom_callbacks_define_their_own_direction',
            'comparability_note': 'Scores depend on the retained cohort, tolerances, options, and scoring callbacks.',
        }
        self.document['stages']['scoring'] = {'status': 'completed'}
        # score_only already explains a not-requested rendering (same
        # precedence as MATLAB); only a plotting run with history plots
        # disabled needs the more specific reason.
        if (single and self._options.get('draw_hist_plots') == 'none'
                and not self._options.get('score_only', False)
                and self.document['stages']['rendering']['status'] not in ('failed', 'partial')):
            self.document['stages']['rendering'] = {'status': 'not_requested',
                                                   'reason': 'single_problem_history_plots_disabled'}
        self._write()

    def add_convergence(self, tolerance, work_hist, work_out, problem_ids):
        """Summarize work arrays actually consumed by profile construction.

        No target, merit, or convergence predicate is recalculated here. NaN
        work is a non-hit in this retained cohort, not a solver diagnosis.
        """
        self._detail_dirty = True
        records = self.document['profiles']['convergence']
        if len(records) >= 256:
            self.add_diagnostic('convergence_summary_limit_reached', 'scoring', limit=256)
            return
        entry = {'tolerance_index': len(records) + 1, 'tolerance': _safe(tolerance),
                 'semantics': 'observed_work_used_by_existing_profile_construction',
                 'denominator': 'all_retained_problem_run_pairs_per_solver',
                 'invalid_initial_separation': None,
                 'invalid_initial_reason': 'not_separately_supplied',
                 'problem_count': len(problem_ids) if type(problem_ids) in (list, tuple) else None,
                 'history_work_semantics': 'first_actual_evaluation_meeting_existing_target',
                 'output_work_semantics': 'total_evaluations_at_passing_returned_output_not_first_hit',
                 'history': [], 'output': []}
        for label, work in (('history', work_hist), ('output', work_out)):
            if type(work) is not np.ndarray or work.ndim != 3 or work.dtype.kind not in 'biuf':
                entry[label + '_reason'] = 'work_array_unavailable'
                continue
            for solver_index in range(work.shape[1]):
                observed = work[:, solver_index, :]
                finite = observed[np.isfinite(observed)]
                entry[label].append({
                    'solver_index': solver_index + 1,
                    'hits': int(finite.size), 'total': int(observed.size),
                    'not_hit': int(np.count_nonzero(np.isnan(observed))),
                    'nonfinite_other': int(np.count_nonzero(np.isinf(observed))),
                    'work_evaluations_min': _safe(np.min(finite)) if finite.size else None,
                    'work_evaluations_max': _safe(np.max(finite)) if finite.size else None,
                })
        records.append(entry)
        # The merged-primary order is passed explicitly by benchmark. It must
        # not be inferred from the collector's insertion order (which can also
        # include failed loads and optional plain reference problems).
        identities = [self._problem(item[0], item[1], 'primary')['id'] for item in problem_ids]
        self.plot_document['target_work'].append({
            'id': f'target-work-{len(records)}', 'tolerance_index': len(records),
            'tolerance': _safe(tolerance), 'axes': ['problem', 'solver', 'run'],
            'axis_values': {'problem': identities,
                            'solver': list(range(1, work_hist.shape[1] + 1)),
                            'run': list(range(1, work_hist.shape[2] + 1))},
            'history': _score_values(work_hist), 'output': _score_values(work_out)})
        self.document['profiles']['work_summary_reason'] = None
        self.document['profiles']['work_summary'] = 'see_convergence'

    def _artifact_paths(self, directory):
        count = 0
        for root, directories, files in os.walk(str(directory), followlinks=False):
            # Junctions are not symlinks to os.walk; prune every link kind.
            directories[:] = sorted(name for name in directories
                                    if not _link_like(Path(root) / name))
            for name in sorted(files):
                path = Path(root) / name
                if not _link_like(path) and path.is_file():
                    yield path
                    count += 1
                    if count >= _MAX_ARTIFACTS:
                        return

    def _harvest(self):
        if self._output_dir is None:
            return
        if (self._output_dir.is_symlink() or not self._output_dir.is_dir()
                or self._identity(self._output_dir) != self._output_identity):
            self.add_diagnostic('artifact_directory_ownership_changed', 'persistence')
            return
        artifacts = []
        for path in self._artifact_paths(self._output_dir):
            if path in (self.path, self.plot_path, self._source_path) or str(path) in self._initial_artifacts:
                continue
            try:
                media_type = mimetypes.guess_type(path.name)[0] or 'application/octet-stream'
                record = {'path': _relative_path(path, self._output_dir),
                          'media_type': media_type,
                          'kind': 'plot' if path.suffix.lower() in ('.pdf', '.png', '.svg') else
                                  ('raw_data' if path.suffix.lower() in ('.h5', '.mat', '.pkl') else 'supporting_file'),
                          'status': 'present'}
                if record['path'] is None:
                    record['path_reason'] = 'cross_volume_relative_path_unavailable'
                if _hashing_capability() == 'unavailable':
                    # The saved files still exist and may be numerically valid;
                    # only their integrity receipt is missing, and it says so.
                    record.update(bytes=None, sha256=None, status='unverified',
                                  reason='secure_artifact_hashing_unavailable_on_platform')
                    self.add_diagnostic('artifact_hash_unavailable', 'persistence',
                                        reason='secure_artifact_hashing_unavailable_on_platform')
                else:
                    record.update(_digest(path, directory=self._output_dir))
                artifacts.append(record)
            except (OSError, NotImplementedError):
                self.add_diagnostic('artifact_hash_unavailable', 'persistence', name=path.name)
        self.document['artifacts'] = artifacts
        if len(artifacts) >= _MAX_ARTIFACTS:
            self.add_diagnostic('artifact_manifest_limit_reached', 'persistence', limit=_MAX_ARTIFACTS)

    def finish(self, error=None):
        """Finalize once; secondary reporting errors never mask ``error``."""
        if self._finished:
            return
        try:
            self._coverage()
            if error is not None:
                self.add_diagnostic('benchmark_exception', 'runtime',
                                    exception_type=type(error).__name__)
                scoring_interrupted = self.document['stages']['scoring']['status'] in ('running', 'failed')
                for stage in _STAGES:
                    if self.document['stages'][stage]['status'] == 'running':
                        render_failure = self.document['stages']['rendering']['status'] in ('failed', 'partial')
                        if (stage == 'numerical' and scoring_interrupted
                                and self.document['coverage']['completed'] > 0):
                            self.document['stages'][stage] = {
                                'status': 'partial', 'reason': 'interrupted_before_all_libraries'}
                        else:
                            self.document['stages'][stage] = {
                                'status': 'unknown' if stage == 'scoring' and render_failure else 'failed',
                                'reason': 'interrupted_before_scores' if stage == 'scoring' and render_failure else 'interrupted_by_exception'}
            else:
                for stage in _STAGES:
                    if self.document['stages'][stage]['status'] == 'running':
                        self.document['stages'][stage] = {'status': 'completed'}
                coverage = self.document['coverage']
                if coverage['load_failed']:
                    self.document['stages']['numerical'] = (
                        {'status': 'partial', 'reason': 'selected_problem_load_failures'}
                        if coverage['completed'] else
                        {'status': 'failed', 'reason': 'all_selected_problem_loads_failed'})
                elif (not coverage['selected'] and not coverage['loaded']
                        and self.document['stages']['numerical']['status'] == 'completed'):
                    self.document['stages']['numerical'] = {
                        'status': 'not_applicable', 'reason': 'empty_selection'}
                if (not coverage['loaded'] and self.document['operation'] == 'benchmark'
                        and self.document['stages']['scoring']['status'] in ('unknown', 'running', 'completed')):
                    # Same vocabulary as MATLAB: nothing was loaded, so no
                    # cohort score exists even if the pipeline returned zeros.
                    self.document['stages']['scoring'] = {
                        'status': 'not_applicable', 'reason': 'no_loaded_problems'}
                if (self.document['operation'] == 'benchmark'
                        and any(problem['role'] == 'plain_reference'
                                and problem['load_status'] == 'failed'
                                for problem in self._problems.values())
                        and self.document['stages']['numerical']['status'] != 'failed'):
                    # This is a failure actually observed during a requested
                    # fresh reference execution, not an inference from absent
                    # legacy metadata. Primary coverage and completed scoring
                    # remain untouched.
                    self.document['stages']['numerical'] = {
                        'status': 'partial', 'reason': 'requested_plain_reference_incomplete'}
            self._harvest()
            if (self.document['stages']['rendering']['status'] == 'completed'
                    and not any(a['kind'] == 'plot' for a in self.document['artifacts'])):
                self.document['stages']['rendering'] = {
                    'status': 'unknown', 'reason': 'no_render_artifacts_observed'}
            # Artifact failures can change persistence/rendering outcomes;
            # compute the overall status only after harvesting their evidence.
            coverage = self.document['coverage']
            stage_failures = any(s['status'] in ('partial', 'failed')
                                 for s in self.document['stages'].values())
            if error is not None:
                self.document['status'] = 'failed'
            elif not coverage['completed'] and not coverage['selected'] and not stage_failures:
                self.document['status'] = 'empty'
            elif stage_failures:
                self.document['status'] = 'partial' if coverage['completed'] else 'failed'
            else:
                self.document['status'] = 'completed'
            self.document['timing']['finished_at'] = _now()
            if self._dropped_diagnostics:
                self.document['diagnostics_omitted'] = self._dropped_diagnostics
            self._write()
            self._finished = True
        except Exception:
            if error is None:
                raise
