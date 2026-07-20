#!/usr/bin/env python3
"""Validate the structure and ownership boundary of a MATLAB-only ZIP."""

from __future__ import annotations

import argparse
from pathlib import Path, PurePosixPath
import stat
import sys
import zipfile


ARCHIVE_ROOT = PurePosixPath("optiprofiler-matlab")
REQUIRED_FILES = {
    ARCHIVE_ROOT / "BUILD_INFO.txt",
    ARCHIVE_ROOT / "LICENSE",
    ARCHIVE_ROOT / "README.rst",
    ARCHIVE_ROOT / "THIRD_PARTY_NOTICES.md",
    ARCHIVE_ROOT / "matlab_problem_libraries.lock",
    ARCHIVE_ROOT / "setup.m",
    ARCHIVE_ROOT / "matlab/examples/example1.m",
    ARCHIVE_ROOT / "matlab/optiprofiler/src/Problem.m",
    ARCHIVE_ROOT / "matlab/optiprofiler/src/benchmark.m",
    ARCHIVE_ROOT / "matlab/optiprofiler/src/resolveProblemLibrary.m",
    ARCHIVE_ROOT / "matlab/optiprofiler/problem_libs/custom/custom_load.m",
    ARCHIVE_ROOT / "matlab/optiprofiler/problem_libs/s2mpj/README.md",
    ARCHIVE_ROOT / "matlab/optiprofiler/problem_libs/s2mpj/s2mpj_load.m",
    ARCHIVE_ROOT / "matlab/optiprofiler/problem_libs/s2mpj/s2mpj_select.m",
    ARCHIVE_ROOT / "matlab/optiprofiler/problem_libs/s2mpj/probinfo_matlab.mat",
    ARCHIVE_ROOT / "matlab/optiprofiler/problem_libs/s2mpj/src/s2mpjlib.m",
}
FORBIDDEN_PARTS = {
    ".git",
    ".github",
    ".idea",
    ".vscode",
    "__pycache__",
    "doc",
    "python",
    "tests",
    "venv",
}
FORBIDDEN_NAMES = {".DS_Store", "runTests.m"}


class ArchiveError(ValueError):
    """Raised when an archive violates the MATLAB distribution contract."""


def validate_archive(path: Path, release_ready: bool = False) -> int:
    try:
        archive = zipfile.ZipFile(path)
    except (OSError, zipfile.BadZipFile) as exc:
        raise ArchiveError(f"Cannot open {path}: {exc}") from exc

    with archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ArchiveError("The archive contains duplicate member names.")
        members = {PurePosixPath(name) for name in names}
        missing = sorted(REQUIRED_FILES - members, key=str)
        if missing:
            raise ArchiveError(f"Required archive members are missing: {missing}")

        for info in archive.infolist():
            member = PurePosixPath(info.filename)
            if member.is_absolute() or ".." in member.parts:
                raise ArchiveError(f"Unsafe archive path: {member}")
            if not member.parts or member.parts[0] != ARCHIVE_ROOT.name:
                raise ArchiveError(f"Member is outside {ARCHIVE_ROOT}: {member}")
            relative_parts = set(member.parts[1:])
            if relative_parts & FORBIDDEN_PARTS or member.name in FORBIDDEN_NAMES:
                raise ArchiveError(f"Development-only member is present: {member}")
            mode = info.external_attr >> 16
            if stat.S_ISLNK(mode):
                raise ArchiveError(f"Symbolic links are not allowed: {member}")

        notice = archive.read(str(ARCHIVE_ROOT / "THIRD_PARTY_NOTICES.md")).decode("utf-8")
        if "GrattonToint/S2MPJ" not in notice or "optiprofiler/s2mpj_matlab" not in notice:
            raise ArchiveError("The S2MPJ provenance notice is incomplete.")
        if release_ready:
            license_candidates = {
                ARCHIVE_ROOT / "licenses/S2MPJ-LICENSE.txt",
                ARCHIVE_ROOT / "licenses/S2MPJ-PERMISSION.txt",
            }
            if not members & license_candidates:
                raise ArchiveError(
                    "A public release requires an S2MPJ license or written "
                    "redistribution permission in the archive."
                )
        build_info = archive.read(str(ARCHIVE_ROOT / "BUILD_INFO.txt")).decode("utf-8")
        if (
            "source_commit=" not in build_info
            or "source_tree=" not in build_info
            or "s2mpj_commit=" not in build_info
        ):
            raise ArchiveError("BUILD_INFO.txt does not record both source revisions.")
        if release_ready and "source_tree=clean" not in build_info:
            raise ArchiveError("A public release archive must be built from a clean source tree.")
    return len(names)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument(
        "--release-ready",
        action="store_true",
        help="also require the upstream S2MPJ redistribution document",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        count = validate_archive(args.archive, release_ready=args.release_ready)
    except ArchiveError as exc:
        print(f"MATLAB ZIP validation error: {exc}", file=sys.stderr)
        return 1
    print(f"Validated {args.archive} with {count} files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
