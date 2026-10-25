"""Release acceptance through the checker CLI, independent of runtime imports."""

from pathlib import Path
import io
import subprocess
import sys
import tarfile
import tempfile
import unittest
import zipfile


CHECKER = Path(__file__).resolve().parents[1] / 'check_release.py'


class ReleaseContractTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        package = self.root / 'python' / 'optiprofiler'
        package.mkdir(parents=True)
        # The checker must read the literal without importing plotting/runtime
        # dependencies, even when package import is deliberately impossible.
        (package / '__init__.py').write_text(
            "__version__ = '2.0.0'\nraise RuntimeError('must not import')\n")

    def check(self, *args):
        return subprocess.run(
            [sys.executable, str(CHECKER), '--root', str(self.root), *args],
            text=True, capture_output=True, timeout=10)

    def test_matching_tag_accepts_literal_package_version_without_import(self):
        result = self.check('--tag', 'v2.0.0')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('2.0.0', result.stdout)

    def test_wrong_tag_blocks_upload(self):
        result = self.check('--tag', 'v1.3.4')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('does not match', result.stderr)

    def distributions(self, wheel_version='2.0.0', source_version='2.0.0', omitted=None):
        dist = self.root / 'dist'
        dist.mkdir()
        files = {
            'optiprofiler/__init__.py': b"__version__ = '2.0.0'\n",
            'optiprofiler/schemas/eval_report.schema.json': b'{}',
            'optiprofiler/schemas/eval_report-v2.schema.json': b'{}',
            'optiprofiler/schemas/plot_data.schema.json': b'{}',
            'optiprofiler/problem_libs/s2mpj/LICENCE.txt': b'upstream license',
            'optiprofiler/problem_libs/s2mpj/THIRD_PARTY_NOTICES.md': b'provenance',
            'optiprofiler/problem_libs/s2mpj/probinfo_python.csv': b'name,n\n',
        }
        if omitted:
            files.pop(omitted, None)
        wheel = dist / ('optiprofiler-' + wheel_version + '-py3-none-any.whl')
        with zipfile.ZipFile(wheel, 'w') as archive:
            for name, value in files.items():
                archive.writestr(name, value)
            archive.writestr('optiprofiler-' + wheel_version + '.dist-info/METADATA',
                             'Name: optiprofiler\nVersion: ' + wheel_version + '\n')
            archive.writestr('optiprofiler-' + wheel_version + '.dist-info/licenses/LICENSE', 'BSD')
        source = dist / ('optiprofiler-' + source_version + '.tar.gz')
        source_files = {'python/' + name: value for name, value in files.items()}
        source_files.update({'LICENSE': b'BSD', 'PKG-INFO':
                             ('Name: optiprofiler\nVersion: ' + source_version + '\n').encode(),
                             'THIRD_PARTY_NOTICES.md': b'notices',
                             'licenses/S2MPJ-LICENSE.txt': b'upstream license'})
        if omitted in source_files:
            source_files.pop(omitted)
        with tarfile.open(source, 'w:gz') as archive:
            for name, value in source_files.items():
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

    def test_empty_artifact_directory_is_rejected(self):
        dist = self.root / 'dist'
        dist.mkdir()
        result = self.check('--dist-dir', str(dist))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('one wheel and one sdist', result.stderr)


if __name__ == '__main__':
    unittest.main()
