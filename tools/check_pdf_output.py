#!/usr/bin/env python3
"""Check controlled smoke-test PDFs for readable text and failed TeX rendering."""

import argparse
import json
from pathlib import Path
import re


RAW_TEX = re.compile(r'\\(?:alpha|beta|rho|tau|frac|mathrm|text|begin|end|left|right)\b')


def inspect_pdf(path, reject_raw_tex=False):
    from pypdf import PdfReader

    pages = PdfReader(str(path)).pages
    if not pages:
        raise ValueError('PDF has no pages: {}'.format(path))
    lengths = []
    for number, page in enumerate(pages, 1):
        text = page.extract_text() or ''
        if not text.strip():
            raise ValueError('PDF page {} has no readable text: {}'.format(number, path))
        raw = RAW_TEX.search(text) if reject_raw_tex else None
        if raw:
            raise ValueError('PDF page {} contains raw TeX {}: {}'.format(number, raw.group(), path))
        lengths.append(len(text.strip()))
    return {'file': str(path), 'pages': len(pages), 'text_lengths': lengths}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--reject-raw-tex', action='store_true',
                        help='Use only for fixtures whose labels should be typeset, not literal TeX.')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    files = sorted(args.root.rglob('*.pdf'))
    results, errors = [], []
    if not files:
        errors.append('No PDF output found under {}'.format(args.root))
    for path in files:
        try:
            results.append(inspect_pdf(path, args.reject_raw_tex))
        except Exception as error:
            errors.append('{}: {}'.format(path, error))
    report = {'passed': not errors, 'pdfs': results, 'errors': errors}
    encoded = json.dumps(report, indent=2) + '\n'
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(encoded, encoding='utf-8')
    print(encoded, end='')
    return 1 if errors else 0


if __name__ == '__main__':
    raise SystemExit(main())
