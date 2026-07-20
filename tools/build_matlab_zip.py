#!/usr/bin/env python3
"""Build the MATLAB-only OptiProfiler ZIP from tracked source files."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import stat
import subprocess
import sys
import zipfile


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_ROOT = PurePosixPath("optiprofiler-matlab")
S2MPJ_ROOT = Path("matlab/optiprofiler/problem_libs/s2mpj")
LOCK_PATH = Path("matlab_problem_libraries.lock")
FIXED_TIME = (1980, 1, 1, 0, 0, 0)

CORE_PATHS = (
    Path("LICENSE"),
    Path("README.rst"),
    Path("THIRD_PARTY_NOTICES.md"),
    LOCK_PATH,
    Path("setup.m"),
)
CORE_TREES = (
    Path("matlab/examples"),
    Path("matlab/optiprofiler/src"),
    Path("matlab/optiprofiler/problem_libs/README.txt"),
    Path("matlab/optiprofiler/problem_libs/custom"),
)
S2MPJ_EXCLUDED_PREFIXES = (".github/", "tests/")


class BuildError(RuntimeError):
    """Raised when a release archive cannot be built safely."""


def run_git(arguments: list[str], cwd: Path = REPOSITORY_ROOT) -> str:
    result = subprocess.run(
        ["git", *arguments],
        cwd=cwd,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise BuildError(result.stderr.strip() or "Git command failed.")
    return result.stdout


def tracked_files(paths: tuple[Path, ...], cwd: Path = REPOSITORY_ROOT) -> list[Path]:
    output = run_git(
        ["ls-files", "-z", "--", *(path.as_posix() for path in paths)], cwd=cwd
    )
    return [Path(value) for value in output.split("\0") if value]


def locked_s2mpj_commit() -> str:
    try:
        data = json.loads((REPOSITORY_ROOT / LOCK_PATH).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BuildError(f"Cannot read {LOCK_PATH}: {exc}") from exc
    matches = [entry for entry in data.get("libraries", []) if entry.get("name") == "s2mpj"]
    if len(matches) != 1 or matches[0].get("role") != "bundled-default":
        raise BuildError("The MATLAB lock must contain one bundled S2MPJ entry.")
    return str(matches[0].get("commit", ""))


def validate_s2mpj_checkout() -> str:
    expected = locked_s2mpj_commit()
    index_line = run_git(["ls-files", "--stage", "--", S2MPJ_ROOT.as_posix()]).strip()
    fields = index_line.partition("\t")[0].split()
    if len(fields) != 3 or fields[0] != "160000" or fields[1] != expected:
        raise BuildError("The bundled S2MPJ gitlink does not match the MATLAB lock.")
    checkout = REPOSITORY_ROOT / S2MPJ_ROOT
    if not checkout.is_dir():
        raise BuildError("The S2MPJ submodule is absent; initialize submodules first.")
    actual = run_git(["rev-parse", "HEAD"], cwd=checkout).strip()
    if actual != expected:
        raise BuildError(
            f"The S2MPJ checkout is {actual}, but the MATLAB lock requires {expected}."
        )
    if run_git(["status", "--porcelain"], cwd=checkout).strip():
        raise BuildError("The S2MPJ checkout contains local changes or untracked files.")
    return expected


def archive_members() -> list[tuple[Path, PurePosixPath]]:
    missing = [path for path in CORE_PATHS if not (REPOSITORY_ROOT / path).is_file()]
    if missing:
        raise BuildError(f"Required files are missing: {missing}")
    tracked_core_paths = set(tracked_files(CORE_PATHS))
    untracked = sorted(set(CORE_PATHS) - tracked_core_paths, key=lambda path: path.as_posix())
    if untracked:
        raise BuildError(f"Required distribution files are not tracked by Git: {untracked}")

    core_files = set(CORE_PATHS)
    core_files.update(tracked_files(CORE_TREES))
    if not all((REPOSITORY_ROOT / path).is_file() for path in core_files):
        raise BuildError("A selected core path is absent or is not a regular file.")

    checkout = REPOSITORY_ROOT / S2MPJ_ROOT
    submodule_files = tracked_files((Path("."),), cwd=checkout)
    submodule_files = [
        path
        for path in submodule_files
        if not path.as_posix().startswith(S2MPJ_EXCLUDED_PREFIXES)
    ]
    if not submodule_files:
        raise BuildError("No tracked S2MPJ distribution files were found.")

    members: list[tuple[Path, PurePosixPath]] = []
    for relative in sorted(core_files, key=lambda path: path.as_posix()):
        members.append((REPOSITORY_ROOT / relative, ARCHIVE_ROOT / relative.as_posix()))
    for relative in sorted(submodule_files, key=lambda path: path.as_posix()):
        source = checkout / relative
        if not source.is_file():
            raise BuildError(f"S2MPJ member is absent or not a regular file: {relative}")
        destination = ARCHIVE_ROOT / S2MPJ_ROOT.as_posix() / relative.as_posix()
        members.append((source, destination))
    return members


def build_archive(output: Path) -> tuple[str, int]:
    s2mpj_commit = validate_s2mpj_checkout()
    members = archive_members()
    output.parent.mkdir(parents=True, exist_ok=True)

    source_commit = run_git(["rev-parse", "HEAD"]).strip()
    distribution_paths = (*CORE_PATHS, *CORE_TREES, S2MPJ_ROOT)
    source_state = "dirty" if run_git(
        ["status", "--porcelain", "--", *(path.as_posix() for path in distribution_paths)]
    ).strip() else "clean"
    build_info = (
        "OptiProfiler MATLAB-only distribution\n"
        f"source_commit={source_commit}\n"
        f"source_tree={source_state}\n"
        f"s2mpj_commit={s2mpj_commit}\n"
    ).encode("utf-8")

    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        info = zipfile.ZipInfo(str(ARCHIVE_ROOT / "BUILD_INFO.txt"), FIXED_TIME)
        info.compress_type = zipfile.ZIP_DEFLATED
        info.external_attr = (stat.S_IFREG | 0o644) << 16
        archive.writestr(info, build_info)

        for source, destination in members:
            info = zipfile.ZipInfo(str(destination), FIXED_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (stat.S_IFREG | 0o644) << 16
            archive.writestr(info, source.read_bytes())

    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    checksum = output.with_suffix(output.suffix + ".sha256")
    checksum.write_text(f"{digest}  {output.name}\n", encoding="ascii")
    return digest, len(members) + 1


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("dist/optiprofiler-matlab.zip"),
        help="archive path relative to the repository root",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    output = args.output if args.output.is_absolute() else REPOSITORY_ROOT / args.output
    try:
        digest, count = build_archive(output)
    except BuildError as exc:
        print(f"MATLAB ZIP build error: {exc}", file=sys.stderr)
        return 1
    print(f"Built {output} with {count} files.")
    print(f"SHA-256: {digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
