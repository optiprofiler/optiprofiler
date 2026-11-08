"""PyPI and Anaconda refuse an existing file alike, so a re-run cannot make them diverge."""

from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[2]
UPLOAD = ROOT / '.github' / 'actions' / 'build' / 'upload_anaconda.sh'
WORKFLOW = ROOT / '.github' / 'workflows' / 'build.yml'


class PublishChannelTests(unittest.TestCase):
    def test_anaconda_neither_replaces_nor_skips_an_existing_file(self):
        commands = [line for line in UPLOAD.read_text(encoding='utf-8').splitlines()
                    if re.match(r'\s*anaconda\b', line)]
        self.assertEqual(len(commands), 1)
        for option in ('--force', '--skip-existing', '--skip', '--interactive', ' -i '):
            self.assertNotIn(option, commands[0])
        self.assertIn(' upload ', commands[0])

    def test_pypi_does_not_skip_an_existing_file(self):
        self.assertEqual(re.findall(r'skip-existing:\s*(\S+)', WORKFLOW.read_text(encoding='utf-8')), ['false'])


if __name__ == '__main__':
    unittest.main()
