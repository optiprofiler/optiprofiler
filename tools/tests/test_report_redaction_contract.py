"""One redaction list: the Python literal, the MATLAB literal and the documentation."""

import ast
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
METADATA = ROOT / 'python' / 'optiprofiler' / 'metadata.py'
REPORT = ROOT / 'matlab' / 'optiprofiler' / 'src' / '+optiprofiler_internal' / 'EvalReport.m'
DOCUMENT = ROOT / 'doc' / 'source' / 'user' / 'eval_report.rst'


def python_literal(name):
    for node in ast.parse(METADATA.read_text(encoding='utf-8')).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], 'id', None) == name:
            return ast.literal_eval(node.value)
    raise AssertionError(f'{name} is not defined in metadata.py')


def matlab_literal(function):
    source = REPORT.read_text(encoding='utf-8')
    body = re.search(r'function pattern = ' + function + r'\(\)(.*?)\n\s*end\b', source, re.S).group(1)
    assignment = body[body.index('pattern ='):]
    return ''.join(re.findall(r"'([^']*)'", assignment))


def documented_rules():
    text = DOCUMENT.read_text(encoding='utf-8')
    return text[text.index('Redacted metadata\n'):text.index('Hash and permission fallbacks\n')]


class RedactionContractTests(unittest.TestCase):
    def test_both_languages_hold_the_same_sensitive_key_list(self):
        self.assertEqual(python_literal('SENSITIVE_KEY_PATTERN'), matlab_literal('sensitiveKeyPattern'))

    def test_both_languages_hold_the_same_url_pattern(self):
        self.assertEqual(python_literal('URL_PATTERN'), matlab_literal('urlPattern'))

    def test_every_listed_key_is_documented_once(self):
        section = documented_rules()
        words = re.findall(r'[a-z]+(?:\[_-\]\?[a-z]+)?', python_literal('SENSITIVE_KEY_PATTERN'))
        self.assertTrue(words)
        for word in words:
            self.assertIn('``' + word.replace('[_-]?', '_') + '``', section)

    def test_matlab_omission_of_path_options_is_documented(self):
        self.assertIn("contains(key, 'path', 'IgnoreCase', true)", REPORT.read_text(encoding='utf-8'))
        self.assertIn('``savepath``, ``report_path``', ' '.join(documented_rules().split()))


if __name__ == '__main__':
    unittest.main()
