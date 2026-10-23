"""PDF acceptance must reject outputs that pass a header-only smoke check."""

import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


TOOLS = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('check_pdf_output', TOOLS / 'check_pdf_output.py')
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


@unittest.skipUnless(importlib.util.find_spec('pypdf'), 'PDF reader is installed in the artifact gate')
class PdfOutputTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def pdf(self, text=None):
        from pypdf import PdfWriter
        from pypdf.generic import DictionaryObject, NameObject, DecodedStreamObject

        writer = PdfWriter()
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
        path = self.root / 'output.pdf'
        with path.open('wb') as stream:
            writer.write(stream)
        return path

    def test_readable_page_passes(self):
        result = checker.inspect_pdf(self.pdf('solver 1'), reject_raw_tex=True)
        self.assertEqual(result['pages'], 1)
        self.assertGreater(result['text_lengths'][0], 0)

    def test_blank_pdf_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'no readable text'):
            checker.inspect_pdf(self.pdf())

    def test_literal_tex_is_rejected_only_for_controlled_fixtures(self):
        path = self.pdf(r'\rho(\alpha)')
        self.assertEqual(checker.inspect_pdf(path)['pages'], 1)
        with self.assertRaisesRegex(ValueError, 'raw TeX'):
            checker.inspect_pdf(path, reject_raw_tex=True)

    def test_empty_output_fails_and_keeps_receipt(self):
        report = self.root / 'receipt.json'
        result = subprocess.run([sys.executable, str(TOOLS / 'check_pdf_output.py'),
                                 str(self.root), '--report', str(report)],
                                capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 1)
        self.assertIn('No PDF output', report.read_text())
