#!/usr/bin/env python3
"""Validate the repository-level OptiProfiler problem-library lock.

The lock is an integration manifest, not a runtime registry.  Core package
code must continue to discover installed external libraries through the API v1
entry-point protocol rather than reading this file at runtime.
"""

import argparse
from importlib import metadata
import json
from pathlib import Path
import re
import subprocess
import sys


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LOCK = REPOSITORY_ROOT / "problem_libraries.lock"
CORE_PROTOCOL = REPOSITORY_ROOT / "python/optiprofiler/problem_libraries.py"

SHA_PATTERN = re.compile(r"^[0-9a-f]{40}$")
REPOSITORY_PATTERN = re.compile(r"^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$")
ENTRY_POINT_PATTERN = re.compile(
    r"^[A-Za-z_][A-Za-z0-9_.]*:[A-Za-z_][A-Za-z0-9_.]*$"
)
NAME_PATTERN = re.compile(r"^[a-z][a-z0-9_]*$")
ALLOWED_ROLES = {"bundled-default", "external", "external-experimental"}
TOP_LEVEL_KEYS = {"schema_version", "problem_library_api_version", "libraries"}
COMMON_LIBRARY_KEYS = {"name", "role", "repository", "commit", "api_version"}
BUNDLED_LIBRARY_KEYS = COMMON_LIBRARY_KEYS | {"gitlink"}
EXTERNAL_LIBRARY_KEYS = COMMON_LIBRARY_KEYS | {
    "distribution",
    "entry_point",
    "legacy_gitlink",
    "legacy_gitlink_commit",
}


class LockError(ValueError):
    """Raised when the integration lock is inconsistent."""


def _reject_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise LockError(f"Duplicate JSON key: {key!r}.")
        result[key] = value
    return result


def _load_json(path):
    try:
        return json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
        )
    except OSError as exc:
        raise LockError(f"Cannot read lock file {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise LockError(f"Invalid JSON in {path}: {exc}") from exc


def _require_exact_keys(value, expected, context):
    actual = set(value)
    missing = sorted(expected - actual)
    unknown = sorted(actual - expected)
    if missing or unknown:
        details = []
        if missing:
            details.append(f"missing {missing}")
        if unknown:
            details.append(f"unknown {unknown}")
        raise LockError(f"{context} has invalid keys: {', '.join(details)}.")


def _core_api_version():
    try:
        source = CORE_PROTOCOL.read_text(encoding="utf-8")
    except OSError as exc:
        raise LockError(f"Cannot read core problem-library protocol: {exc}") from exc
    match = re.search(r"^PROBLEM_LIBRARY_API_VERSION\s*=\s*(\d+)\s*$", source, re.MULTILINE)
    if match is None:
        raise LockError(f"Cannot find PROBLEM_LIBRARY_API_VERSION in {CORE_PROTOCOL}.")
    return int(match.group(1))


def _index_gitlink(path):
    if not isinstance(path, str) or not path or Path(path).is_absolute():
        raise LockError("A gitlink path must be a nonempty repository-relative string.")
    result = subprocess.run(
        ["git", "ls-files", "--stage", "--", path],
        cwd=str(REPOSITORY_ROOT),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        raise LockError(f"Cannot inspect gitlink {path}: {result.stderr.strip()}")
    line = result.stdout.strip()
    if not line:
        raise LockError(f"Expected tracked gitlink is missing: {path}.")
    metadata_part, _, indexed_path = line.partition("\t")
    fields = metadata_part.split()
    if len(fields) != 3 or indexed_path != path or fields[0] != "160000":
        raise LockError(f"Expected {path} to be a mode-160000 gitlink; got {line!r}.")
    return fields[1]


def _validate_library(entry, lock_api_version, index):
    context = f"libraries[{index}]"
    if not isinstance(entry, dict):
        raise LockError(f"{context} must be an object.")
    role = entry.get("role")
    if role not in ALLOWED_ROLES:
        raise LockError(f"{context}.role must be one of {sorted(ALLOWED_ROLES)}.")

    expected_keys = (
        BUNDLED_LIBRARY_KEYS
        if role == "bundled-default"
        else EXTERNAL_LIBRARY_KEYS
    )
    if role != "bundled-default" and "legacy_gitlink" not in entry:
        expected_keys = expected_keys - {"legacy_gitlink", "legacy_gitlink_commit"}
    _require_exact_keys(entry, expected_keys, context)

    name = entry["name"]
    if not isinstance(name, str) or NAME_PATTERN.fullmatch(name) is None:
        raise LockError(f"{context}.name must be a lowercase problem-library identifier.")
    repository = entry["repository"]
    if not isinstance(repository, str) or REPOSITORY_PATTERN.fullmatch(repository) is None:
        raise LockError(f"{context}.repository must have the form owner/repository.")
    commit = entry["commit"]
    if not isinstance(commit, str) or SHA_PATTERN.fullmatch(commit) is None:
        raise LockError(f"{context}.commit must be a full lowercase 40-character Git SHA.")
    if type(entry["api_version"]) is not int or entry["api_version"] != lock_api_version:
        raise LockError(
            f"{context}.api_version must equal problem_library_api_version {lock_api_version}."
        )

    if role == "bundled-default":
        if name != "s2mpj":
            raise LockError(
                "The bundled-default library must be s2mpj in this "
                "compatibility stage."
            )
        if _index_gitlink(entry["gitlink"]) != commit:
            raise LockError(f"{context}.commit does not match its bundled gitlink.")
        return

    distribution = entry["distribution"]
    if not isinstance(distribution, str) or not distribution.strip():
        raise LockError(f"{context}.distribution must be a nonempty string.")
    entry_point = entry["entry_point"]
    if not isinstance(entry_point, str) or ENTRY_POINT_PATTERN.fullmatch(entry_point) is None:
        raise LockError(f"{context}.entry_point must have the form module:attribute.")
    if "legacy_gitlink" in entry:
        legacy_commit = entry["legacy_gitlink_commit"]
        if not isinstance(legacy_commit, str) or SHA_PATTERN.fullmatch(legacy_commit) is None:
            raise LockError(f"{context}.legacy_gitlink_commit must be a full Git SHA.")
        if _index_gitlink(entry["legacy_gitlink"]) != legacy_commit:
            raise LockError(f"{context}.legacy_gitlink_commit does not match the Git index.")


def validate_lock(path=DEFAULT_LOCK):
    """Return a validated lock mapping or raise ``LockError``."""
    data = _load_json(Path(path))
    if not isinstance(data, dict):
        raise LockError("The lock root must be an object.")
    _require_exact_keys(data, TOP_LEVEL_KEYS, "lock")
    if data["schema_version"] != 1 or type(data["schema_version"]) is not int:
        raise LockError("schema_version must be integer 1.")
    lock_api_version = data["problem_library_api_version"]
    if type(lock_api_version) is not int or lock_api_version < 1:
        raise LockError("problem_library_api_version must be a positive integer.")
    core_api_version = _core_api_version()
    if lock_api_version != core_api_version:
        raise LockError(
            f"Lock API version {lock_api_version} does not match core API "
            f"version {core_api_version}."
        )
    libraries = data["libraries"]
    if not isinstance(libraries, list) or not libraries:
        raise LockError("libraries must be a nonempty array.")
    for index, entry in enumerate(libraries):
        _validate_library(entry, lock_api_version, index)

    names = [entry["name"] for entry in libraries]
    repositories = [entry["repository"] for entry in libraries]
    distributions = [
        entry["distribution"] for entry in libraries if entry["role"] != "bundled-default"
    ]
    for label, values in (
        ("library names", names),
        ("repositories", repositories),
        ("external distributions", distributions),
    ):
        duplicates = sorted({value for value in values if values.count(value) > 1})
        if duplicates:
            raise LockError(f"Duplicate {label}: {duplicates}.")
    bundled = [entry for entry in libraries if entry["role"] == "bundled-default"]
    if len(bundled) != 1:
        raise LockError("Exactly one bundled-default library is required.")
    return data


def external_matrix(data):
    """Return the GitHub Actions matrix for separately installed libraries."""
    include = []
    for entry in data["libraries"]:
        if entry["role"] == "bundled-default":
            continue
        include.append(
            {
                "library": entry["name"],
                "role": entry["role"],
                "repository": entry["repository"],
                "commit": entry["commit"],
                "distribution": entry["distribution"],
                "entry_point": entry["entry_point"],
            }
        )
    return {"include": include}


def _normalized_distribution_name(name):
    return re.sub(r"[-_.]+", "-", name).lower()


def verify_installed(data, library_name):
    """Verify one installed external provider against its locked metadata."""
    matches = [entry for entry in data["libraries"] if entry["name"] == library_name]
    if not matches:
        raise LockError(f"Unknown locked problem library: {library_name}.")
    entry = matches[0]
    if entry["role"] == "bundled-default":
        raise LockError(f"{library_name} is bundled and has no external entry point.")

    all_entry_points = metadata.entry_points()
    if hasattr(all_entry_points, "select"):
        candidates = list(
            all_entry_points.select(group="optiprofiler.problem_libraries", name=library_name)
        )
    else:
        candidates = [
            candidate
            for candidate in all_entry_points.get("optiprofiler.problem_libraries", ())
            if candidate.name == library_name
        ]
    if len(candidates) != 1:
        raise LockError(
            f"Expected exactly one installed entry point for {library_name}; "
            f"found {len(candidates)}."
        )
    candidate = candidates[0]
    if candidate.value != entry["entry_point"]:
        raise LockError(
            f"Installed entry point {candidate.value!r} does not match locked value "
            f"{entry['entry_point']!r}."
        )
    distribution = getattr(candidate, "dist", None)
    installed_distribution = getattr(distribution, "name", None)
    if installed_distribution is None and distribution is not None:
        installed_distribution = distribution.metadata.get("Name")
    if installed_distribution is None or _normalized_distribution_name(
        installed_distribution
    ) != _normalized_distribution_name(entry["distribution"]):
        raise LockError(
            f"Installed distribution {installed_distribution!r} does not match locked "
            f"distribution {entry['distribution']!r}."
        )

    try:
        from optiprofiler.problem_libraries import (
            load_problem_library,
            resolve_problem_library,
        )

        reference = resolve_problem_library(library_name)
        if reference.source != "entry_point":
            raise LockError(
                f"Core resolved {library_name} from {reference.source!r}, not its "
                "installed entry point."
            )
        if reference.locator != entry["entry_point"]:
            raise LockError(
                f"Core resolved {reference.locator!r}, not locked entry point "
                f"{entry['entry_point']!r}."
            )
        plugin = load_problem_library(reference)
    except LockError:
        raise
    except Exception as exc:
        raise LockError(
            f"Core could not resolve and load installed library {library_name}: {exc}"
        ) from exc
    if getattr(plugin, "name", None) != library_name:
        raise LockError(f"Factory for {library_name} returned a differently named plugin.")
    if getattr(plugin, "api_version", None) != entry["api_version"]:
        raise LockError(
            f"Plugin {library_name} API version does not match locked API version "
            f"{entry['api_version']}."
        )


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--lock",
        type=Path,
        default=DEFAULT_LOCK,
        help="Lock file to validate.",
    )
    parser.add_argument(
        "--external-matrix",
        action="store_true",
        help="Print the external-library GitHub Actions matrix as compact JSON.",
    )
    parser.add_argument(
        "--verify-installed",
        metavar="NAME",
        help="Verify an installed external provider against the lock.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)
    try:
        data = validate_lock(args.lock)
        if args.verify_installed:
            verify_installed(data, args.verify_installed)
        if args.external_matrix:
            print(json.dumps(external_matrix(data), separators=(",", ":"), sort_keys=True))
        elif not args.verify_installed:
            print(
                f"problem_libraries.lock: OK ({len(data['libraries'])} libraries, "
                f"API v{data['problem_library_api_version']})"
            )
    except LockError as exc:
        print(f"problem_libraries.lock: ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
