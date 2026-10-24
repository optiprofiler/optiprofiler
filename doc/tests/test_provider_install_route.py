"""The documented provider installation must use the exact commits of the lock."""
import json
from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[2]
GUIDE = ROOT / 'doc' / 'source' / 'user' / 'problem_libraries' / 'python.rst'
INSTALL = re.compile(r'pip install "(?P<distribution>[\w-]+) @ git\+https://github\.com/(?P<repository>[\w/-]+)@(?P<commit>[0-9a-f]{40})"')


class ProviderInstallRoute(unittest.TestCase):
    def test_documented_git_installs_match_the_lock(self):
        lock = json.loads((ROOT / 'problem_libraries.lock').read_text(encoding='utf-8'))
        external = {entry['distribution']: entry for entry in lock['libraries'] if 'distribution' in entry}
        documented = {match['distribution']: match for match in INSTALL.finditer(GUIDE.read_text(encoding='utf-8'))}
        self.assertEqual(sorted(documented), sorted(external))
        for distribution, match in documented.items():
            with self.subTest(distribution=distribution):
                self.assertEqual(match['repository'], external[distribution]['repository'])
                self.assertEqual(match['commit'], external[distribution]['commit'])

    def test_unpublished_providers_are_not_presented_as_pypi_installs(self):
        text = GUIDE.read_text(encoding='utf-8')
        self.assertIsNone(re.search(r'pip install (--upgrade )?optiprofiler-(pycutest|solar|rs13)\b', text))


if __name__ == '__main__':
    unittest.main()
