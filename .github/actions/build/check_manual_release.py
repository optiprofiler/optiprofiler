"""Fail-closed checks for an explicitly requested, PyPI-only maintenance release.

Manual dispatch does not use a v* tag: the legacy tag workflow also publishes
to Anaconda. Keep the requested source/version and the built artifact hashes
bound together across jobs so a build-only dispatch cannot become a release.
"""

import ast
from email.parser import Parser
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tarfile
from urllib.error import HTTPError
from urllib.request import urlopen
import zipfile


def context():
    sha = os.environ.get('EXPECTED_SHA', '')
    version = os.environ.get('EXPECTED_VERSION', '')
    if os.environ.get('GITHUB_REF') != 'refs/heads/paper/v1.3.4-errata':
        raise ValueError('Manual maintenance publication requires the paper/v1.3.4-errata branch.')
    if not re.fullmatch(r'[0-9a-f]{40}', sha) or sha != os.environ.get('GITHUB_SHA'):
        raise ValueError('The explicitly approved SHA must equal the dispatched workflow SHA.')
    actual = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    if actual != sha:
        raise ValueError('The actual source checkout differs from the approved SHA.')
    subprocess.run(['git', 'diff', '--quiet', '--ignore-submodules=none', 'HEAD', '--'], check=True)
    if not re.fullmatch(r'1\.\d+\.\d+', version):
        raise ValueError('Manual maintenance publication requires an explicit final 1.x version.')
    tree = ast.parse(Path('python/optiprofiler/__init__.py').read_text(encoding='utf-8'))
    declared = [ast.literal_eval(node.value) for node in tree.body
                if isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == '__version__' for t in node.targets)]
    if declared != [version]:
        raise ValueError('The requested version differs from source __version__.')
    return sha, version


ROOT_NOTICES = ('LICENSE', 'THIRD_PARTY_NOTICES.md', 'licenses/S2MPJ-LICENSE.txt')
PACKAGE_NOTICES = (
    'python/optiprofiler/problem_libs/s2mpj/LICENCE.txt',
    'python/optiprofiler/problem_libs/s2mpj/THIRD_PARTY_NOTICES.md',
    'python/optiprofiler/problem_libs/solar/runtime/solar/LICENSE',
)


def artifact_receipt(sha, version):
    files = sorted(Path('dist').iterdir())
    wheels = [p for p in files if p.name.endswith('.whl')]
    sdists = [p for p in files if p.name.endswith('.tar.gz')]
    if len(wheels) != 1 or len(sdists) != 1 or len(files) != 2:
        raise ValueError('Expected exactly one wheel and one sdist, with no extra upload targets.')
    with zipfile.ZipFile(wheels[0]) as archive:
        if len(archive.namelist()) != len(set(archive.namelist())):
            raise ValueError('Duplicate wheel members are not allowed.')
        names = [name for name in archive.namelist() if name.endswith('.dist-info/METADATA')]
        if len(names) != 1:
            raise ValueError('The wheel must have one METADATA file.')
        wheel_metadata = Parser().parsestr(archive.read(names[0]).decode('utf-8'))
        # Old setuptools flattens license files directly into dist-info;
        # newer versions retain paths below dist-info/licenses. Require the
        # reviewed bytes in either supported layout, never just a filename.
        for source in ROOT_NOTICES:
            matches = [name for name in archive.namelist()
                       if '.dist-info/' in name and Path(name).name == Path(source).name
                       and archive.read(name) == Path(source).read_bytes()]
            if not matches:
                raise ValueError('Missing or changed wheel distribution notice: ' + source)
        for source in PACKAGE_NOTICES:
            if archive.read(source[len('python/'):]) != Path(source).read_bytes():
                raise ValueError('Missing or changed wheel package notice: ' + source)
    with tarfile.open(sdists[0], 'r:gz') as archive:
        if len(archive.getnames()) != len(set(archive.getnames())):
            raise ValueError('Duplicate sdist members are not allowed.')
        names = [name for name in archive.getnames() if name.count('/') == 1 and name.endswith('/PKG-INFO')]
        if len(names) != 1:
            raise ValueError('The sdist must have one root PKG-INFO file.')
        sdist_metadata = Parser().parsestr(archive.extractfile(names[0]).read().decode('utf-8'))
        prefix = names[0].split('/')[0] + '/'
        for source in ROOT_NOTICES + PACKAGE_NOTICES:
            member = archive.extractfile(prefix + source)
            if member is None or member.read() != Path(source).read_bytes():
                raise ValueError('Missing or changed sdist distribution notice: ' + source)
    for metadata in (wheel_metadata, sdist_metadata):
        if metadata['Name'].lower() != 'optiprofiler' or metadata['Version'] != version:
            raise ValueError('Artifact package name/version differs from the approved release.')
    return {'sha': sha, 'version': version, 'artifacts': {
        path.name: {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size}
        for path in files}}


def main():
    sha, version = context()
    mode = sys.argv[1]
    if mode == 'admit':
        print('Manual PyPI-only release admitted:', sha, version)
        return
    receipt = artifact_receipt(sha, version)
    proof = Path('release-proof.json')
    if mode == 'record':
        proof.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    elif mode == 'verify':
        if json.loads(proof.read_text(encoding='utf-8')) != receipt:
            raise ValueError('Downloaded artifacts differ from the approved build receipt.')
        # A timeout/service error is not proof of absence. Only a real 404
        # admits a new release; never use skip-existing to mask a partial upload.
        try:
            with urlopen('https://pypi.org/pypi/optiprofiler/' + version + '/json', timeout=30):
                raise ValueError('This version already exists on PyPI; inspect it instead of uploading again.')
        except HTTPError as error:
            if error.code != 404:
                raise
    else:
        raise ValueError('Expected admit, record or verify.')
    print(json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
