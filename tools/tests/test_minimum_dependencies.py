"""The declared dependency floors are the ones the minimum-dependency check exercises."""

from pathlib import Path
import importlib.util
import json
import sys
import tempfile
import unittest
from unittest import mock

TOOL = Path(__file__).resolve().parents[1] / 'check_minimum_dependencies.py'
SPECIFICATION = importlib.util.spec_from_file_location('check_minimum_dependencies', TOOL)
tool = importlib.util.module_from_spec(SPECIFICATION)
SPECIFICATION.loader.exec_module(tool)


class MinimumDependencyTests(unittest.TestCase):
    def test_floors_cover_the_apis_in_use(self):
        floors = tool.read_floors()
        # nanmin/nanmax(initial=) and an importable numpy.testing on a stock
        # Linux Python 3.10; Figure.clear with subfigures.
        self.assertGreaterEqual(tool.version_tuple(floors['numpy']), (1, 22, 2))
        self.assertGreaterEqual(tool.version_tuple(floors['matplotlib']), (3, 5, 2))
        self.assertEqual(set(floors), {'matplotlib', 'numpy', 'scipy', 'pandas', 'h5py', 'pypdf'})

    def test_smoke_requirement_is_taken_from_the_tests_extra(self):
        self.assertEqual([item.split('>=')[0] for item in tool.test_requirements()], ['jsonschema'])

    def test_a_dependency_without_a_plain_floor_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            pyproject = Path(directory) / 'pyproject.toml'
            pyproject.write_text('dependencies = [\n    "numpy",\n]\n', encoding='utf-8')
            with self.assertRaises(SystemExit):
                tool.read_floors(pyproject)

    def test_floor_python_requires_the_exact_floors(self):
        floors = {'numpy': '1.22.2', 'scipy': '1.10.0'}
        exact = {'numpy': '1.22.2', 'scipy': '1.10.0'}
        later = {'numpy': '1.26.0', 'scipy': '1.11.2'}
        self.assertEqual(tool.floor_violations(floors, exact, exact=True), [])
        self.assertEqual(tool.floor_violations(floors, later, exact=False), [])
        self.assertEqual(len(tool.floor_violations(floors, later, exact=True)), 2)
        self.assertEqual(len(tool.floor_violations(floors, {'numpy': '1.21.6'}, exact=False)), 2)

    def test_a_failure_leaves_a_receipt_and_a_nonzero_status(self):
        def failing(uv, python, wheel, output):
            raise SystemExit('numpy 1.19.5 instead of the declared floor')

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'receipts'
            arguments = ['check', '--wheel', 'x.whl', '--python', '3.8', '--output-dir', str(output), '--uv', 'uv']
            with mock.patch.object(sys, 'argv', arguments), mock.patch.object(tool, 'check_python', failing):
                self.assertEqual(tool.main(), 1)
            records = json.loads((output / 'receipt.json').read_text(encoding='utf-8'))
        self.assertEqual([record['status'] for record in records], ['FAIL'])
        self.assertIn('declared floor', records[0]['reason'])


if __name__ == '__main__':
    unittest.main()
