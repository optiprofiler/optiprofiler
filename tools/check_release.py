#!/usr/bin/env python3
"""Check source, tag and distribution identity before an upload is permitted."""

import argparse
import ast
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import zipfile


# Revision of each bundled problem-library subset, shipped in both artifacts;
# it repeats the bundled-default entries of the lock, which the sdist ships.
BUNDLED_RECORD = 'optiprofiler/problem_libs/bundled_revisions.json'
LOCK = 'problem_libraries.lock'

REQUIRED_PACKAGE_FILES = (
    'optiprofiler/__init__.py',
    'optiprofiler/problem_libraries.py',
    BUNDLED_RECORD,
    'optiprofiler/schemas/eval_report.schema.json',
    'optiprofiler/schemas/eval_report-v2.schema.json',
    'optiprofiler/schemas/plot_data.schema.json',
    'optiprofiler/problem_libs/s2mpj/LICENCE.txt',
    'optiprofiler/problem_libs/s2mpj/THIRD_PARTY_NOTICES.md',
    'optiprofiler/problem_libs/s2mpj/probinfo_python.csv',
)

# Repository-level notices required in the sdist. The wheel carries the
# bundled S2MPJ licence and notice through REQUIRED_PACKAGE_FILES.
REQUIRED_SOURCE_FILES = (
    'THIRD_PARTY_NOTICES.md',
    'licenses/S2MPJ-LICENSE.txt',
    LOCK,
)


def source_version(root):
    # Importing the package initializes plotting and numerical dependencies.
    # Release identity must be checkable before those dependencies are installed.
    source = root / 'python' / 'optiprofiler' / '__init__.py'
    return literal_version(source.read_text(encoding='utf-8'))


def literal_version(source):
    values = [node.value.value for node in ast.parse(source).body
              if isinstance(node, ast.Assign)
              and any(isinstance(target, ast.Name) and target.id == '__version__'
                      for target in node.targets)
              and isinstance(node.value, ast.Constant)
              and isinstance(node.value.value, str)]
    if len(values) != 1:
        raise ValueError('Source must declare exactly one literal package version.')
    return values[0]


def gitlink_commit(root, path):
    result = subprocess.run(['git', 'ls-files', '--stage', '--', path], cwd=str(root),
                            capture_output=True, text=True, check=False)
    fields = result.stdout.partition('\t')[0].split()
    if result.returncode != 0 or len(fields) != 3 or fields[0] != '160000':
        raise ValueError(f'{path} is not a gitlink in the checkout {root}: '
                         f'{result.stderr.strip() or result.stdout.strip() or "not tracked"}.')
    return fields[1]


def check_bundled_revisions(root):
    """Check the shipped revision record against the lock and the gitlinks of ``root``.

    Returns the checked source bytes that both artifacts must carry unchanged,
    keyed by their source-relative names, and the recorded libraries.
    """
    lock = json.loads((root / LOCK).read_text(encoding='utf-8'))
    record = json.loads((root / 'python' / BUNDLED_RECORD).read_text(encoding='utf-8'))
    expected = [{key: entry.get(key) for key in ('name', 'repository', 'commit', 'gitlink')}
                for entry in lock.get('libraries', []) if entry.get('role') == 'bundled-default']
    if not expected or record.get('libraries') != expected:
        raise ValueError(f'{BUNDLED_RECORD} does not repeat the bundled-default entries of {LOCK}: '
                         f'{record.get("libraries")} != {expected}.')
    for entry in expected:
        commit = gitlink_commit(root, entry['gitlink'])
        if commit != entry['commit']:
            raise ValueError(f'The {entry["name"]} gitlink is {commit}, but {LOCK} and '
                             f'{BUNDLED_RECORD} record {entry["commit"]}.')
    source_files = {'python/' + BUNDLED_RECORD: (root / 'python' / BUNDLED_RECORD).read_bytes(),
                    LOCK: (root / LOCK).read_bytes()}
    return source_files, expected


def bundled_problem_libraries(source):
    """Read the shipped problem-library names declared by the core, without importing it."""
    values = [node.value for node in ast.parse(source).body
              if isinstance(node, ast.Assign)
              and any(isinstance(target, ast.Name) and target.id == '_BUNDLED_PROBLEM_LIBRARIES'
                      for target in node.targets)]
    names = ast.literal_eval(values[0]) if len(values) == 1 else None
    if not isinstance(names, tuple) or not all(isinstance(name, str) for name in names):
        raise ValueError('problem_libraries.py must declare _BUNDLED_PROBLEM_LIBRARIES once, '
                         'as a literal tuple of names.')
    return set(names)


def check_problem_library_trees(path, files, package_prefix):
    """Reject provider trees (e.g. a stale V1 checkout) outside the bundled set."""
    shipped = bundled_problem_libraries(files[package_prefix + 'optiprofiler/problem_libraries.py'].decode('utf-8'))
    libs = package_prefix + 'optiprofiler/problem_libs/'
    trees = {name[len(libs):].split('/', 1)[0] for name in files
             if name.startswith(libs) and '/' in name[len(libs):]}
    extra = sorted(trees - shipped)
    if extra:
        raise ValueError(f'{path.name} ships problem-library trees outside the bundled set '
                         f'{sorted(shipped)}: {extra}.')


def check_distribution(path, version, source_files=None):
    if path.suffix == '.whl':
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            files = {name: archive.read(name) for name in names if not name.endswith('/')}
        metadata_names = [name for name in files if name.endswith('.dist-info/METADATA')]
        license_names = [name for name in files if name.endswith('.dist-info/licenses/LICENSE')
                         or name.endswith('.dist-info/LICENSE')]
        required = REQUIRED_PACKAGE_FILES
        package_prefix = ''
        init = 'optiprofiler/__init__.py'
        if not path.name.startswith('optiprofiler-' + version + '-'):
            raise ValueError(f'{path.name} filename does not match source version {version}.')
    else:
        prefix = 'optiprofiler-' + version + '/'
        with tarfile.open(path, 'r:gz') as archive:
            members = archive.getmembers()
            names = [member.name for member in members]
            files = {member.name: archive.extractfile(member).read()
                     for member in members if member.isfile()}
        metadata_names = [prefix + 'PKG-INFO'] if prefix + 'PKG-INFO' in files else []
        license_names = [prefix + 'LICENSE'] if prefix + 'LICENSE' in files else []
        required = (tuple(prefix + 'python/' + name for name in REQUIRED_PACKAGE_FILES)
                    + tuple(prefix + name for name in REQUIRED_SOURCE_FILES))
        package_prefix = prefix + 'python/'
        init = prefix + 'python/optiprofiler/__init__.py'
        if path.name != 'optiprofiler-' + version + '.tar.gz':
            raise ValueError(f'{path.name} filename does not match source version {version}.')
    if len(names) != len(set(names)):
        raise ValueError(f'{path.name} has duplicate archive members.')
    if len(metadata_names) != 1:
        raise ValueError(f'{path.name} must contain exactly one package metadata record.')
    metadata = BytesParser().parsebytes(files[metadata_names[0]])
    if metadata.get('Name') != 'optiprofiler' or metadata.get('Version') != version:
        raise ValueError(f'{path.name} metadata does not match optiprofiler {version}.')
    for name in required:
        if not files.get(name):
            raise ValueError(f'{path.name} is missing required resource {name}.')
    check_problem_library_trees(path, files, package_prefix)
    for name, value in (source_files or {}).items():
        if path.suffix == '.whl':
            if not name.startswith('python/'):
                continue
            member = name[len('python/'):]
        else:
            member = prefix + name
        if files.get(member) != value:
            raise ValueError(f'{path.name} member {member} differs from the checked source file {name}.')
    if len(license_names) != 1 or not files[license_names[0]]:
        raise ValueError(f'{path.name} is missing the core LICENSE.')
    if literal_version(files[init].decode('utf-8')) != version:
        raise ValueError(f'{path.name} packaged version does not match source version {version}.')
    return {'file': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--tag', help='release tag, for example v2.0.0')
    parser.add_argument('--dist-dir', type=Path, help='check one wheel and one sdist')
    parser.add_argument('--manifest', type=Path, help='write the checked artifact identities')
    args = parser.parse_args()
    try:
        version = source_version(args.root)
        if args.tag is not None and args.tag != 'v' + version:
            raise ValueError(f'Tag {args.tag!r} does not match package version {version!r}.')
        source_files, bundled = check_bundled_revisions(args.root)
        artifacts = []
        if args.dist_dir is not None:
            wheels = list(args.dist_dir.glob('*.whl'))
            sources = list(args.dist_dir.glob('*.tar.gz'))
            if len(wheels) != 1 or len(sources) != 1:
                raise ValueError('Distribution directory must contain one wheel and one sdist.')
            artifacts = [check_distribution(path, version, source_files) for path in wheels + sources]
        if args.manifest is not None:
            args.manifest.write_text(json.dumps({'version': version, 'tag': args.tag,
                                                'bundled_problem_libraries': bundled,
                                                'artifacts': artifacts}, indent=2) + '\n',
                                     encoding='utf-8')
    except (OSError, SyntaxError, ValueError, tarfile.TarError, zipfile.BadZipFile) as exc:
        print(f'Release validation error: {exc}', file=sys.stderr)
        return 1
    print(f'Validated package version {version}.')
    for entry in bundled:
        print(f'Bundled {entry["name"]} revision {entry["commit"]} matches {LOCK} and its gitlink.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
