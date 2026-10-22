"""Report text is bounded honestly and provenance keeps no credentials."""
import json

import numpy as np
import pytest

from optiprofiler import Problem, ProblemReference, benchmark
from optiprofiler.eval_report import schema_resource
from optiprofiler.metadata import TEXT_LIMIT, TRUNCATION_MARKER, bounded_text, provenance_text
from optiprofiler.tests.eval_report_contract import assert_valid


def test_clipped_text_is_marked_and_bounded():
    text = 'x' * (TEXT_LIMIT + 50)
    for encode in (bounded_text, provenance_text):
        clipped = encode(text)
        assert len(clipped) == TEXT_LIMIT
        assert clipped.endswith(TRUNCATION_MARKER)
        assert encode('short text') == 'short text'


@pytest.mark.parametrize('source, expected', [
    ('https://doi.org/10.1080/10556788.2025.2490640', 'https://doi.org/10.1080/10556788.2025.2490640'),
    ('https://user:pw@example.org/data/set.csv?token=abc#row', 'https://example.org/data/set.csv [url parameters removed]'),
    ('Table 3 of http://example.org:8080/paper?sig=1', 'Table 3 of http://example.org:8080/paper [url parameters removed]'),
    ('best known value, see /home/someone/notes.txt', 'best known value, see [redacted_absolute_path]'),
])
def test_provenance_keeps_locations_without_secrets(source, expected):
    assert provenance_text(source) == expected
    assert 'pw' not in provenance_text(source) and 'token' not in provenance_text(source)


def test_general_metadata_still_redacts_urls():
    assert bounded_text('https://example.org/x?token=abc') == '[redacted_url]'


def finish(fun, x0):
    x = np.ones_like(x0)
    fun(x)
    return x


def stay(fun, x0):
    fun(x0)
    return x0


def test_report_shows_credential_free_and_marked_reference_sources(tmp_path):
    long_source = 'best known value reported by ' + 'y' * 300
    for index, (source, check) in enumerate((
            ('https://doi.org/10.1080/10556788.2025.2490640?utm=mail',
             lambda text: text == 'https://doi.org/10.1080/10556788.2025.2490640 [url parameters removed]'),
            (long_source, lambda text: len(text) == TEXT_LIMIT and text.endswith(TRUNCATION_MARKER)))):
        report_path = tmp_path / f'report-{index}.json'
        problem = Problem(lambda x: float(x @ x), [2.0],
                          reference=ProblemReference(0.0, 'feasible', source))
        benchmark([finish, stay], problem=problem, score_only=True, silent=True, n_jobs=1,
                  max_eval_factor=2, report_path=str(report_path))
        report = json.loads(report_path.read_text(encoding='utf-8'))
        assert_valid(report, schema_resource('eval_report', 2))
        reported = report['problems'][0]['reference_records'][0]['source']
        assert check(reported), reported
