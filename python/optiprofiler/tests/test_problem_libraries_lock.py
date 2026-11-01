import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
CHECKER = REPOSITORY_ROOT / "tools/check_problem_libraries_lock.py"
LOCK = REPOSITORY_ROOT / "problem_libraries.lock"
BUNDLED_DECLARATION = re.compile(r"^_BUNDLED_PROBLEM_LIBRARIES\s*=.*$", re.MULTILINE)


# The integration lock and its validator intentionally stay at repository
# level; they are not runtime package data.  Source checkouts exercise these
# tests, while wheel/sdist smoke tests skip them after installing the package.
pytestmark = pytest.mark.skipif(
    not CHECKER.is_file() or not LOCK.is_file(),
    reason="problem-library lock validation requires a repository checkout",
)


def _run_checker(*arguments):
    return subprocess.run(
        [sys.executable, str(CHECKER), *map(str, arguments)],
        cwd=str(REPOSITORY_ROOT),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_problem_libraries_lock", CHECKER)
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    return checker


def test_problem_library_lock_matches_core_protocol_and_gitlinks():
    result = _run_checker()

    assert result.returncode == 0, result.stderr
    assert "4 libraries, API v1" in result.stdout


def test_problem_library_lock_generates_external_ci_matrix():
    result = _run_checker("--external-matrix")

    assert result.returncode == 0, result.stderr
    matrix = json.loads(result.stdout)
    assert [entry["library"] for entry in matrix["include"]] == [
        "pycutest",
        "solar",
        "rs13",
    ]
    assert [entry["role"] for entry in matrix["include"]] == [
        "external",
        "external",
        "external-experimental",
    ]
    assert all(len(entry["commit"]) == 40 for entry in matrix["include"])


def test_problem_library_lock_rejects_core_api_mismatch(tmp_path):
    data = json.loads(LOCK.read_text(encoding="utf-8"))
    data["problem_library_api_version"] = 2
    invalid_lock = tmp_path / "problem_libraries.lock"
    invalid_lock.write_text(json.dumps(data), encoding="utf-8")

    result = _run_checker("--lock", invalid_lock)

    assert result.returncode == 1
    assert "does not match core API version 1" in result.stderr


def test_problem_library_lock_rejects_external_source_in_core(tmp_path):
    data = json.loads(LOCK.read_text(encoding="utf-8"))
    rs13 = next(entry for entry in data["libraries"] if entry["name"] == "rs13")
    rs13["name"] = "custom"
    invalid_lock = tmp_path / "problem_libraries.lock"
    invalid_lock.write_text(json.dumps(data), encoding="utf-8")

    result = _run_checker("--lock", invalid_lock)

    assert result.returncode == 1
    assert "must not be tracked inside the core package" in result.stderr


def test_shipped_revision_record_must_repeat_the_bundled_lock_entry(tmp_path):
    data = json.loads(LOCK.read_text(encoding="utf-8"))
    s2mpj = next(entry for entry in data["libraries"] if entry["role"] == "bundled-default")
    s2mpj["repository"] = "optiprofiler/another_s2mpj"
    invalid_lock = tmp_path / "problem_libraries.lock"
    invalid_lock.write_text(json.dumps(data), encoding="utf-8")

    result = _run_checker("--lock", invalid_lock)

    assert result.returncode == 1
    assert "bundled_revisions.json does not match the bundled-default lock entry" in result.stderr


@pytest.mark.parametrize(
    "declaration, message",
    [
        ("_BUNDLED_PROBLEM_LIBRARIES = ('custom', 's2mpj_python')", "missing from _BUNDLED_PROBLEM_LIBRARIES"),
        ("", "Cannot find _BUNDLED_PROBLEM_LIBRARIES"),
    ],
)
def test_bundled_default_must_be_declared_by_the_core(tmp_path, monkeypatch, declaration, message):
    checker = _load_checker()
    core = tmp_path / "problem_libraries.py"
    source = checker.CORE_PROTOCOL.read_text(encoding="utf-8")
    core.write_text(BUNDLED_DECLARATION.sub(declaration, source), encoding="utf-8")
    monkeypatch.setattr(checker, "CORE_PROTOCOL", core)

    with pytest.raises(checker.LockError, match=message):
        checker.validate_lock()


def test_integration_workflow_covers_the_shipped_revision_record():
    workflow = REPOSITORY_ROOT / ".github/workflows/problem-library-integration.yml"
    triggers = workflow.read_text(encoding="utf-8").split("workflow_dispatch:")[0]

    assert triggers.count("- python/optiprofiler/problem_libs/bundled_revisions.json\n") == 2
