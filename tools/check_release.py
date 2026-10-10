#!/usr/bin/env python3
"""Check source, tag and distribution identity before an upload is permitted."""

import argparse
import ast
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path
import sys
import tarfile
import zipfile


REQUIRED_PACKAGE_FILES = (
    'optiprofiler/__init__.py',
    'optiprofiler/schemas/eval_report.schema.json',
    'optiprofiler/schemas/eval_report-v2.schema.json',
    'optiprofiler/schemas/plot_data.schema.json',
    'optiprofiler/problem_libs/s2mpj/LICENCE.txt',
    'optiprofiler/problem_libs/s2mpj/THIRD_PARTY_NOTICES.md',
    'optiprofiler/problem_libs/s2mpj/probinfo_python.csv',
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


def check_distribution(path, version):
    if path.suffix == '.whl':
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            files = {name: archive.read(name) for name in names if not name.endswith('/')}
        metadata_names = [name for name in files if name.endswith('.dist-info/METADATA')]
        license_names = [name for name in files if name.endswith('.dist-info/licenses/LICENSE')
                         or name.endswith('.dist-info/LICENSE')]
        required = REQUIRED_PACKAGE_FILES
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
        required = tuple(prefix + 'python/' + name for name in REQUIRED_PACKAGE_FILES)
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
        artifacts = []
        if args.dist_dir is not None:
            wheels = list(args.dist_dir.glob('*.whl'))
            sources = list(args.dist_dir.glob('*.tar.gz'))
            if len(wheels) != 1 or len(sources) != 1:
                raise ValueError('Distribution directory must contain one wheel and one sdist.')
            artifacts = [check_distribution(path, version) for path in wheels + sources]
        if args.manifest is not None:
            args.manifest.write_text(json.dumps({'version': version, 'tag': args.tag,
                                                'artifacts': artifacts}, indent=2) + '\n',
                                     encoding='utf-8')
    except (OSError, SyntaxError, ValueError, tarfile.TarError, zipfile.BadZipFile) as exc:
        print(f'Release validation error: {exc}', file=sys.stderr)
        return 1
    print(f'Validated package version {version}.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
