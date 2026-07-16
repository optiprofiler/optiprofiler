#!/usr/bin/env python3
"""Validate the MATLAB problem-library integration lock."""

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
LOCK_PATH = REPOSITORY_ROOT / "matlab_problem_libraries.lock"

SHA_PATTERN = re.compile(r"^[0-9a-f]{40}$")
NAME_PATTERN = re.compile(r"^[a-z][a-z0-9_]*$")
REPOSITORY_PATTERN = re.compile(r"^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$")
FUNCTION_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
ALLOWED_ROLES = {"bundled-default", "external"}
ALLOWED_PLATFORMS = {"linux", "macos", "windows"}
TOP_LEVEL_KEYS = {"schema_version", "problem_library_api_version", "libraries"}
COMMON_LIBRARY_KEYS = {
    "name",
    "role",
    "repository",
    "commit",
    "api_version",
    "select_function",
    "load_function",
    "collect_info_function",
    "platforms",
}
BUNDLED_LIBRARY_KEYS = COMMON_LIBRARY_KEYS | {"gitlink"}
EXTERNAL_LIBRARY_KEYS = COMMON_LIBRARY_KEYS | {"install_directory"}


class LockError(ValueError):
    """Raised when the MATLAB integration lock is inconsistent."""


def _reject_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise LockError(f"Duplicate JSON key: {key!r}.")
        result[key] = value
    return result


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


def _git_index_entry(path):
    result = subprocess.run(
        ["git", "ls-files", "--stage", "--", path],
        cwd=REPOSITORY_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        raise LockError(f"Cannot inspect {path}: {result.stderr.strip()}")
    return result.stdout.strip()


def _validate_function_name(value, context, optional=False):
    if optional and value is None:
        return
    if not isinstance(value, str) or FUNCTION_PATTERN.fullmatch(value) is None:
        suffix = " or null" if optional else ""
        raise LockError(f"{context} must be a MATLAB function name{suffix}.")


def _validate_library(entry, api_version, index):
    context = f"libraries[{index}]"
    if not isinstance(entry, dict):
        raise LockError(f"{context} must be an object.")
    role = entry.get("role")
    if role not in ALLOWED_ROLES:
        raise LockError(f"{context}.role must be one of {sorted(ALLOWED_ROLES)}.")
    expected_keys = BUNDLED_LIBRARY_KEYS if role == "bundled-default" else EXTERNAL_LIBRARY_KEYS
    _require_exact_keys(entry, expected_keys, context)

    name = entry["name"]
    if not isinstance(name, str) or NAME_PATTERN.fullmatch(name) is None:
        raise LockError(f"{context}.name must be a lowercase MATLAB library identifier.")
    repository = entry["repository"]
    if not isinstance(repository, str) or REPOSITORY_PATTERN.fullmatch(repository) is None:
        raise LockError(f"{context}.repository must have the form owner/repository.")
    if not isinstance(entry["commit"], str) or SHA_PATTERN.fullmatch(entry["commit"]) is None:
        raise LockError(f"{context}.commit must be a full lowercase 40-character Git SHA.")
    if type(entry["api_version"]) is not int or entry["api_version"] != api_version:
        raise LockError(f"{context}.api_version must equal {api_version}.")

    _validate_function_name(entry["select_function"], f"{context}.select_function")
    _validate_function_name(entry["load_function"], f"{context}.load_function")
    _validate_function_name(
        entry["collect_info_function"],
        f"{context}.collect_info_function",
        optional=True,
    )
    platforms = entry["platforms"]
    if (
        not isinstance(platforms, list)
        or not platforms
        or len(platforms) != len(set(platforms))
        or any(platform not in ALLOWED_PLATFORMS for platform in platforms)
    ):
        raise LockError(
            f"{context}.platforms must be a nonempty unique list drawn from "
            f"{sorted(ALLOWED_PLATFORMS)}."
        )

    if role == "bundled-default":
        if name != "s2mpj":
            raise LockError("The bundled-default MATLAB library must be s2mpj.")
        gitlink = entry["gitlink"]
        if not isinstance(gitlink, str) or Path(gitlink).is_absolute():
            raise LockError(f"{context}.gitlink must be a repository-relative path.")
        line = _git_index_entry(gitlink)
        metadata, separator, indexed_path = line.partition("\t")
        fields = metadata.split()
        if not separator or indexed_path != gitlink or len(fields) != 3 or fields[0] != "160000":
            raise LockError(f"Expected a mode-160000 gitlink at {gitlink}; got {line!r}.")
        if fields[1] != entry["commit"]:
            raise LockError(f"{context}.commit does not match its bundled gitlink.")
        return

    install_directory = entry["install_directory"]
    if (
        not isinstance(install_directory, str)
        or NAME_PATTERN.fullmatch(install_directory) is None
        or install_directory != name
    ):
        raise LockError(f"{context}.install_directory must equal the library name.")
    tracked_path = f"matlab/optiprofiler/problem_libs/{install_directory}"
    if _git_index_entry(tracked_path):
        raise LockError(
            f"External MATLAB provider {name!r} must not be tracked inside core at "
            f"{tracked_path}."
        )


def validate_lock(path=LOCK_PATH):
    try:
        data = json.loads(
            Path(path).read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
        )
    except OSError as exc:
        raise LockError(f"Cannot read MATLAB lock {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise LockError(f"Invalid JSON in {path}: {exc}") from exc

    if not isinstance(data, dict):
        raise LockError("The MATLAB lock root must be an object.")
    _require_exact_keys(data, TOP_LEVEL_KEYS, "lock")
    if type(data["schema_version"]) is not int or data["schema_version"] != 1:
        raise LockError("schema_version must be integer 1.")
    api_version = data["problem_library_api_version"]
    if type(api_version) is not int or api_version != 1:
        raise LockError("problem_library_api_version must be integer 1.")
    libraries = data["libraries"]
    if not isinstance(libraries, list) or not libraries:
        raise LockError("libraries must be a nonempty array.")
    for index, entry in enumerate(libraries):
        _validate_library(entry, api_version, index)

    names = [entry["name"] for entry in libraries]
    repositories = [entry["repository"] for entry in libraries]
    for label, values in (("names", names), ("repositories", repositories)):
        duplicates = sorted({value for value in values if values.count(value) > 1})
        if duplicates:
            raise LockError(f"Duplicate library {label}: {duplicates}.")
    if sum(entry["role"] == "bundled-default" for entry in libraries) != 1:
        raise LockError("Exactly one bundled-default MATLAB library is required.")
    return data


def external_matrix(data):
    """Return the GitHub Actions matrix derived from external lock entries."""

    include = []
    for entry in data["libraries"]:
        if entry["role"] != "external":
            continue
        include.append(
            {
                "library": entry["name"],
                "repository": entry["repository"],
                "commit": entry["commit"],
                "install_directory": entry["install_directory"],
                "select_function": entry["select_function"],
                "load_function": entry["load_function"],
                "collect_info_function": entry["collect_info_function"] or "",
            }
        )
    return {"include": include}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--external-matrix",
        action="store_true",
        help="print the locked external-library GitHub Actions matrix",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        data = validate_lock()
    except LockError as exc:
        print(f"MATLAB problem-library lock error: {exc}", file=sys.stderr)
        return 1
    if args.external_matrix:
        print(json.dumps(external_matrix(data), separators=(",", ":")))
        return 0
    print(f"Validated {LOCK_PATH.name}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
