"""Release acceptance through the checker CLI, independent of runtime imports."""

from pathlib import Path
import io
import json
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
import zipfile


CHECKER = Path(__file__).resolve().parents[1] / 'check_release.py'
COMMIT = '1c8a06fd2c2516c6138bf999a54f6d3dcdf6218f'
GITLINK = 'python/optiprofiler/problem_libs/s2mpj'
RECORD = 'optiprofiler/problem_libs/bundled_revisions.json'
S2MPJ = {'name': 's2mpj', 'repository': 'optiprofiler/s2mpj_python', 'commit': COMMIT, 'gitlink': GITLINK}


def lock_bytes(commit=COMMIT):
    entry = dict(S2MPJ, commit=commit, role='bundled-default', api_version=1)
    return json.dumps({'schema_version': 1, 'problem_library_api_version': 1,
                       'libraries': [entry]}, indent=2).encode()


def record_bytes(commit=COMMIT):
    return json.dumps({'schema_version': 1, 'libraries': [dict(S2MPJ, commit=commit)]}, indent=2).encode()


class ReleaseContractTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        package = self.root / 'python' / 'optiprofiler'
        (package / 'problem_libs').mkdir(parents=True)
        # The checker must read the literal without importing plotting/runtime
        # dependencies, even when package import is deliberately impossible.
        (package / '__init__.py').write_text(
            "__version__ = '2.0.0'\nraise RuntimeError('must not import')\n")
        (self.root / 'problem_libraries.lock').write_bytes(lock_bytes())
        (self.root / 'python' / RECORD).write_bytes(record_bytes())
        self.set_gitlink(COMMIT)

    def set_gitlink(self, commit):
        subprocess.run(['git', 'init', '-q', str(self.root)], check=True)
        subprocess.run(['git', 'update-index', '--add', '--cacheinfo', '160000,' + commit + ',' + GITLINK],
                       cwd=str(self.root), check=True)

    def check(self, *args):
        return subprocess.run(
            [sys.executable, str(CHECKER), '--root', str(self.root), *args],
            text=True, capture_output=True, timeout=10)

    def test_matching_tag_accepts_literal_package_version_without_import(self):
        result = self.check('--tag', 'v2.0.0')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('2.0.0', result.stdout)
        self.assertIn('s2mpj revision ' + COMMIT, result.stdout)

    def test_record_must_repeat_the_lock(self):
        (self.root / 'python' / RECORD).write_bytes(record_bytes('0' * 40))
        result = self.check()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('does not repeat the bundled-default entries', result.stderr)

    def test_record_and_lock_must_match_the_gitlink(self):
        self.set_gitlink('f' * 40)
        result = self.check()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('gitlink is ' + 'f' * 40, result.stderr)

    def test_artifact_record_must_equal_the_checked_source(self):
        for artifact, member in (('wheel', RECORD), ('sdist', 'python/' + RECORD),
                                 ('sdist', 'problem_libraries.lock')):
            with self.subTest(artifact=artifact, member=member):
                result = self.check('--dist-dir', str(self.distributions(changed=(artifact, member))))
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('differs from the checked source file', result.stderr)
                shutil.rmtree(self.root / 'dist')

    def test_manifest_records_the_bundled_revision(self):
        manifest = self.root / 'manifest.json'
        result = self.check('--dist-dir', str(self.distributions()), '--manifest', str(manifest))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(manifest.read_text())['bundled_problem_libraries'], [S2MPJ])

    def test_wrong_tag_blocks_upload(self):
        result = self.check('--tag', 'v1.3.4')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('does not match', result.stderr)

    def distributions(self, wheel_version='2.0.0', source_version='2.0.0', omitted=None, extra=(),
                      changed=None):
        dist = self.root / 'dist'
        dist.mkdir()
        files = {
            'optiprofiler/__init__.py': b"__version__ = '2.0.0'\n",
            'optiprofiler/problem_libraries.py': b"_BUNDLED_PROBLEM_LIBRARIES = ('custom', 's2mpj')\n",
            'optiprofiler/problem_libs/__init__.py': b'',
            'optiprofiler/problem_libs/README.txt': b'guide',
            RECORD: record_bytes(),
            'optiprofiler/problem_libs/custom/custom_tools.py': b'',
            'optiprofiler/schemas/eval_report.schema.json': b'{}',
            'optiprofiler/schemas/eval_report-v2.schema.json': b'{}',
            'optiprofiler/schemas/plot_data.schema.json': b'{}',
            'optiprofiler/problem_libs/s2mpj/LICENCE.txt': b'upstream license',
            'optiprofiler/problem_libs/s2mpj/LICENSE': b'adapter license',
            'optiprofiler/problem_libs/s2mpj/THIRD_PARTY_NOTICES.md': b'provenance',
            'optiprofiler/problem_libs/s2mpj/probinfo_python.csv': b'name,n\n',
        }
        if omitted:
            files.pop(omitted, None)
        files.update({name: b'stale' for name in extra})
        wheel = dist / ('optiprofiler-' + wheel_version + '-py3-none-any.whl')
        with zipfile.ZipFile(wheel, 'w') as archive:
            for name, value in files.items():
                archive.writestr(name, record_bytes('0' * 40) if changed == ('wheel', name) else value)
            archive.writestr('optiprofiler-' + wheel_version + '.dist-info/METADATA',
                             'Name: optiprofiler\nVersion: ' + wheel_version + '\n')
            archive.writestr('optiprofiler-' + wheel_version + '.dist-info/licenses/LICENSE', 'BSD')
        source = dist / ('optiprofiler-' + source_version + '.tar.gz')
        source_files = {'python/' + name: value for name, value in files.items()}
        source_files.update({'LICENSE': b'BSD', 'PKG-INFO':
                             ('Name: optiprofiler\nVersion: ' + source_version + '\n').encode(),
                             'THIRD_PARTY_NOTICES.md': b'notices',
                             'licenses/S2MPJ-LICENSE.txt': b'upstream license',
                             'problem_libraries.lock': lock_bytes()})
        if omitted in source_files:
            source_files.pop(omitted)
        with tarfile.open(source, 'w:gz') as archive:
            for name, value in source_files.items():
                if changed == ('sdist', name):
                    value = value.replace(COMMIT.encode(), b'0' * 40)
                member = tarfile.TarInfo('optiprofiler-' + source_version + '/' + name)
                member.size = len(value)
                archive.addfile(member, io.BytesIO(value))
        return dist

    def test_matching_artifacts_are_accepted(self):
        result = self.check('--dist-dir', str(self.distributions()), '--tag', 'v2.0.0')
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_wheel_version_mismatch_is_rejected(self):
        result = self.check('--dist-dir', str(self.distributions(wheel_version='1.3.4')))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('does not match', result.stderr)

    def test_sdist_version_mismatch_is_rejected(self):
        result = self.check('--dist-dir', str(self.distributions(source_version='1.3.4')))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('does not match', result.stderr)

    def test_missing_upstream_license_blocks_distribution(self):
        result = self.check('--dist-dir', str(self.distributions(
            omitted='optiprofiler/problem_libs/s2mpj/LICENCE.txt')))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('LICENCE.txt', result.stderr)

    def test_missing_source_notice_blocks_distribution(self):
        result = self.check('--dist-dir', str(self.distributions(omitted='licenses/S2MPJ-LICENSE.txt')))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('S2MPJ-LICENSE.txt', result.stderr)

    def test_missing_schema_blocks_distribution(self):
        result = self.check('--dist-dir', str(self.distributions(
            omitted='optiprofiler/schemas/eval_report-v2.schema.json')))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('eval_report-v2.schema.json', result.stderr)

    def test_stale_provider_tree_blocks_distribution(self):
        # V1.3.4 tracked pycutest and solar submodules inside problem_libs.
        for tree in ('pycutest', 'solar'):
            with self.subTest(tree=tree):
                result = self.check('--dist-dir', str(self.distributions(
                    extra=['optiprofiler/problem_libs/' + tree + '/' + tree + '_tools.py'])))
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('outside the bundled set', result.stderr)
                self.assertIn(repr(tree), result.stderr)
                shutil.rmtree(self.root / 'dist')

    def test_empty_artifact_directory_is_rejected(self):
        dist = self.root / 'dist'
        dist.mkdir()
        result = self.check('--dist-dir', str(dist))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('one wheel and one sdist', result.stderr)


if __name__ == '__main__':
    unittest.main()
