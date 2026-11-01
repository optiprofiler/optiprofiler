"""One redaction list: the Python literal, the MATLAB literal and the documentation."""

import ast
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
METADATA = ROOT / 'python' / 'optiprofiler' / 'metadata.py'
REPORT = ROOT / 'matlab' / 'optiprofiler' / 'src' / '+optiprofiler_internal' / 'EvalReport.m'
PYTHON_TESTS = ROOT / 'python' / 'optiprofiler' / 'tests' / 'test_report_text_safety.py'
MATLAB_TESTS = ROOT / 'matlab' / 'optiprofiler' / 'tests' / 'unit_tests' / 'TestReportTextSafety.m'
DOCUMENT = ROOT / 'doc' / 'source' / 'user' / 'eval_report.rst'


def python_literal(name, source=METADATA):
    for node in ast.parse(source.read_text(encoding='utf-8')).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], 'id', None) == name:
            return ast.literal_eval(node.value)
    raise AssertionError(f'{name} is not defined in {source.name}')


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

    def test_both_languages_test_the_same_glued_urls(self):
        source = MATLAB_TESTS.read_text(encoding='utf-8')
        cell = re.search(r'glued = \{(.*?)\};', source, re.S).group(1)
        self.assertEqual(python_literal('GLUED_URLS', PYTHON_TESTS), re.findall(r"'([^']*)'", cell))

    def test_both_languages_test_the_same_history_file_names(self):
        source = MATLAB_TESTS.read_text(encoding='utf-8')
        cell = re.search(r'names = \{ \.\.\.(.*?)\};\n\s*for k', source, re.S).group(1)
        texts = re.findall(r"'([^']*)'", cell)
        self.assertEqual([item for pair in python_literal('HISTORY_FILE_NAMES', PYTHON_TESTS) for item in pair], texts)

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
