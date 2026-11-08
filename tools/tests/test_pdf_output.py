"""Exercise PDF acceptance through the same CLI used by the ZIP smoke gate."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from pypdf import PdfWriter
from pypdf.generic import DecodedStreamObject, DictionaryObject, NameObject


CHECKER = Path(__file__).resolve().parents[1] / 'check_pdf_output.py'


class PdfOutputTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.receipt = self.root / 'receipt.json'

    def pdf(self, texts, name='output.pdf'):
        writer = PdfWriter()
        for text in texts:
            page = writer.add_blank_page(width=200, height=200)
            if text is not None:
                font = DictionaryObject({NameObject('/Type'): NameObject('/Font'),
                                         NameObject('/Subtype'): NameObject('/Type1'),
                                         NameObject('/BaseFont'): NameObject('/Helvetica')})
                page[NameObject('/Resources')] = DictionaryObject({NameObject('/Font'):
                    DictionaryObject({NameObject('/F1'): font})})
                content = DecodedStreamObject()
                escaped = text.replace('\\', '\\\\').replace('(', '\\(').replace(')', '\\)')
                content.set_data(('BT /F1 12 Tf 10 100 Td (' + escaped + ') Tj ET').encode('ascii'))
                page[NameObject('/Contents')] = content
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('wb') as stream:
            writer.write(stream)
        return path

    def check(self, *options):
        return subprocess.run(
            [sys.executable, str(CHECKER), str(self.root), '--report', str(self.receipt), *options],
            capture_output=True, text=True, timeout=10)

    def test_readable_page_passes(self):
        self.pdf(['solver 1'])
        result = self.check('--reject-raw-tex')
        self.assertEqual(result.returncode, 0, result.stderr)
        receipt = json.loads(self.receipt.read_text())
        self.assertTrue(receipt['passed'])
        self.assertEqual(receipt['pdfs'][0]['pages'], 1)
        self.assertGreater(receipt['pdfs'][0]['text_lengths'][0], 0)

    def test_literal_multiplication_tex_is_rejected(self):
        self.pdf([r'10 \times 2'])
        result = self.check('--reject-raw-tex')
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn('raw TeX', self.receipt.read_text())

    def test_blank_page_is_rejected_even_after_a_readable_page(self):
        self.pdf(['solver 1', None])
        result = self.check()
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn('page 2 has no readable text', self.receipt.read_text())

    def test_zero_page_pdf_is_rejected(self):
        self.pdf([])
        result = self.check()
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn('PDF has no pages', self.receipt.read_text())

    def test_literal_tex_check_is_opt_in(self):
        self.pdf([r'\rho(\alpha)'])
        result = self.check()
        self.assertEqual(result.returncode, 0, result.stderr)
        result = self.check('--reject-raw-tex')
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn('raw TeX', self.receipt.read_text())

    def test_missing_output_fails_with_a_receipt(self):
        result = self.check()
        self.assertEqual(result.returncode, 1, result.stdout)
        receipt = json.loads(self.receipt.read_text())
        self.assertFalse(receipt['passed'])
        self.assertIn('No PDF output', receipt['errors'][0])

    def test_corrupt_pdf_does_not_hide_other_results(self):
        self.pdf(['solver 1'])
        (self.root / 'broken.pdf').write_bytes(b'not a PDF')
        result = self.check()
        self.assertEqual(result.returncode, 1, result.stdout)
        receipt = json.loads(self.receipt.read_text())
        self.assertEqual(len(receipt['pdfs']), 1)
        self.assertEqual(len(receipt['errors']), 1)
        self.assertIn('broken.pdf', receipt['errors'][0])

    def test_nested_outputs_and_all_readable_pages_are_checked(self):
        self.pdf(['solver 1', 'solver 2'], 'nested outputs/profile.pdf')
        result = self.check('--reject-raw-tex')
        self.assertEqual(result.returncode, 0, result.stderr)
        receipt = json.loads(self.receipt.read_text())
        self.assertEqual(receipt['pdfs'][0]['pages'], 2)
        self.assertEqual(len(receipt['pdfs'][0]['text_lengths']), 2)


if __name__ == '__main__':
    unittest.main()
