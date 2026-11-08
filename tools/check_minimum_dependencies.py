#!/usr/bin/env python3
"""Exercise a built wheel on the dependency floors declared in pyproject.toml.

For each requested Python this creates a clean environment with uv, installs
the lowest wheels that the declared floors allow, installs the given wheel
without its dependencies on top of them, runs ``pip check`` and then
``tools/smoke_python_distribution.py`` outside the checkout: real solvers, a
save and a solver-free load.

On Python 3.8, the oldest supported Python, every dependency has to resolve to
exactly its declared floor, so a floor is what was exercised and not merely a
number below it. A later Python has no wheels of the oldest versions; there the
floors are lower limits and the versions used are recorded in the receipt.
"""

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
FLOOR_PYTHON = '3.8'


def declared_list(text, key):
    """The strings of the first TOML array ``key = [...]`` (this needs no TOML parser on Python 3.8)."""
    match = re.search(r'^' + re.escape(key) + r'\s*=\s*\[(.*?)^\]', text, re.M | re.S)
    if match is None:
        raise SystemExit(f'pyproject.toml declares no {key}.')
    return re.findall(r'"([^"]+)"', match.group(1))


def read_floors(pyproject=ROOT / 'pyproject.toml'):
    """Map every core dependency to its declared floor. Each one must have the form ``name>=version``."""
    floors = {}
    for requirement in declared_list(pyproject.read_text(encoding='utf-8'), 'dependencies'):
        match = re.fullmatch(r'([A-Za-z0-9._-]+)>=([0-9][0-9.]*)', requirement)
        if match is None:
            raise SystemExit(f'The core dependency {requirement!r} does not declare a plain floor name>=version.')
        floors[match.group(1)] = match.group(2)
    return floors


def test_requirements(pyproject=ROOT / 'pyproject.toml'):
    """The requirements of the smoke tool that are not core dependencies."""
    return [requirement for requirement in declared_list(pyproject.read_text(encoding='utf-8'), 'tests')
            if requirement.startswith('jsonschema')]


def version_tuple(version):
    return tuple(int(part) for part in re.findall(r'\d+', version)[:3])


def floor_violations(floors, installed, exact):
    """Why ``installed`` does not exercise ``floors``: a lower version, or any other version if ``exact``."""
    problems = []
    for name, floor in floors.items():
        version = installed.get(name.lower())
        if version is None:
            problems.append(f'{name} is not installed')
        elif version_tuple(version) < version_tuple(floor) or (exact and version_tuple(version) != version_tuple(floor)):
            problems.append(f'{name} {version} instead of the declared floor {floor}')
    return problems


def run(command, **options):
    print('+', ' '.join(str(part) for part in command), flush=True)
    subprocess.run([str(part) for part in command], check=True, **options)


def check_python(uv, python, wheel, output):
    output.mkdir(parents=True)
    environment = output / 'venv'
    interpreter = environment / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python')
    floors = read_floors()
    requirements = output / 'floors.in'
    requirements.write_text(''.join(f'{name}>={floor}\n' for name, floor in floors.items()), encoding='utf-8')
    install = [uv, 'pip', 'install', '--python', interpreter]
    run([uv, 'venv', '--python', python, environment])
    run(install + ['--resolution', 'lowest-direct', '--only-binary', ':all:', '-r', requirements])
    run(install + ['--only-binary', ':all:'] + test_requirements())
    run(install + ['--no-deps', wheel])
    run([uv, 'pip', 'check', '--python', interpreter])
    listed = subprocess.run([str(uv), 'pip', 'list', '--python', str(interpreter), '--format', 'json'],
                            check=True, capture_output=True, text=True).stdout
    installed = {item['name'].lower(): item['version'] for item in json.loads(listed)}
    problems = floor_violations(floors, installed, exact=python == FLOOR_PYTHON)
    if problems:
        raise SystemExit(f'Python {python} did not exercise the declared floors: ' + '; '.join(problems) + '.')
    version = re.match(r'optiprofiler-([^-]+)-', wheel.name).group(1)
    arguments = [interpreter, ROOT / 'tools' / 'smoke_python_distribution.py', '--source-root', ROOT,
                 '--expected-version', version, '--output-dir', output / 'smoke']
    environment_variables = {key: value for key, value in os.environ.items() if key != 'PYTHONPATH'}
    environment_variables['MPLBACKEND'] = 'Agg'
    run(arguments, cwd=output, env=environment_variables)
    return {'python': python, 'exact_floors': python == FLOOR_PYTHON,
            'versions': {name: installed[name.lower()] for name in floors}, 'status': 'PASS'}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--wheel', type=Path, required=True, help='The built wheel to install.')
    parser.add_argument('--python', action='append', required=True, help='A Python version; repeat it for several.')
    parser.add_argument('--output-dir', type=Path, required=True, help='A new directory for environments and receipts.')
    parser.add_argument('--uv', default=shutil.which('uv'), help='The uv executable.')
    args = parser.parse_args()
    if args.uv is None:
        raise SystemExit('uv is required: install it or give --uv.')
    if args.output_dir.exists():
        raise SystemExit(f'The output directory {args.output_dir} already exists.')
    wheel = args.wheel.resolve()
    output = args.output_dir.resolve()
    records = []
    for python in args.python:
        try:
            records.append(check_python(args.uv, python, wheel, output / ('python-' + python)))
        except (SystemExit, subprocess.CalledProcessError) as error:
            records.append({'python': python, 'status': 'FAIL', 'reason': str(error)})
    output.mkdir(parents=True, exist_ok=True)
    (output / 'receipt.json').write_text(json.dumps(records, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(records, indent=2))
    return int(any(record['status'] != 'PASS' for record in records))


if __name__ == '__main__':
    sys.exit(main())
