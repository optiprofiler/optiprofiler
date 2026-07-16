import json
from pathlib import Path
import subprocess
import sys


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
CHECKER = REPOSITORY_ROOT / "tools/check_problem_libraries_lock.py"
LOCK = REPOSITORY_ROOT / "problem_libraries.lock"


def _run_checker(*arguments):
    return subprocess.run(
        [sys.executable, str(CHECKER), *map(str, arguments)],
        cwd=str(REPOSITORY_ROOT),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


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
