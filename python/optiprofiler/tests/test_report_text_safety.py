"""Report text is bounded honestly and provenance keeps no credentials."""
import json
import time

import numpy as np
import pytest

from optiprofiler import Problem, ProblemReference, benchmark
from optiprofiler.eval_report import schema_resource
from optiprofiler.metadata import TEXT_LIMIT, TRUNCATION_MARKER, bounded_text, provenance_text, safe_metadata
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


CANARY = 'CANARYpw'
CREDENTIAL_URLS = [f'{scheme}://user:{CANARY}@example.org/repo.git' for scheme in (
    'https', 'HTTPS', 'Http', 'git+https', 'GIT+SSH', 'git+ssh', 'ssh', 'SSH', 'ftp', 's3', 'custom.scheme-1')]


@pytest.mark.parametrize('url', CREDENTIAL_URLS)
def test_url_is_redacted_in_any_scheme_and_letter_case(url):
    assert bounded_text(url) == '[redacted_url]'
    assert bounded_text(f'built from {url} (pinned)') == 'built from [redacted_url] (pinned)'


@pytest.mark.parametrize('url', CREDENTIAL_URLS)
def test_provenance_strips_embedded_credentials_of_any_scheme(url):
    scheme = url.split('://')[0].lower()
    assert provenance_text(url) == f'{scheme}://example.org/repo.git [url parameters removed]'
    assert provenance_text(f'Data from {url} and more') == \
        f'Data from {scheme}://example.org/repo.git [url parameters removed] and more'


# The same list as in TestReportTextSafety.m; tools/tests/test_report_redaction_contract.py compares them.
GLUED_URLS = [
    'https://user:A@example.org/x,https://u2:CANARYpw@example.org/y',
    'https://example.org/x,https://u2:CANARYpw@example.org/y',
    'https://user:A@example.org,https://u2:CANARYpw@example.org/y',
    'https://example.org,https://u2:CANARYpw@example.org/y',
    'https://example.org/x;git+ssh://u2:CANARYpw@example.org/y',
    'HTTPS://example.org/x,SSH://u2:CANARYpw@example.org/y',
]


@pytest.mark.parametrize('text', GLUED_URLS)
def test_url_glued_to_another_keeps_no_user_information(text):
    assert provenance_text(text) == '[redacted_url]'
    assert provenance_text(f'Data from {text} and more') == 'Data from [redacted_url] and more'
    assert bounded_text(text) == '[redacted_url]'


def test_url_in_a_query_leaves_with_the_query():
    assert provenance_text(f'https://example.org/x?next=https://u2:{CANARY}@example.org/y') == \
        'https://example.org/x [url parameters removed]'


@pytest.mark.parametrize('text', ['C://Users/someone/x', 'not a url: a:b', 'mailto:someone@example.org'])
def test_text_without_a_url_is_not_treated_as_one(text):
    assert bounded_text(text) != '[redacted_url]'


def test_redacted_url_keeps_text_bounded():
    text = 'https://user:pw@example.org/' + 'x' * (2 * TEXT_LIMIT)
    assert bounded_text(text) == '[redacted_url]'
    clipped = bounded_text('a' * TEXT_LIMIT + ' ' + text)
    assert 'pw' not in clipped and len(clipped) == TEXT_LIMIT


def test_long_text_without_a_url_is_handled_in_linear_time():
    text = 'a' * 200000
    start = time.perf_counter()
    assert bounded_text(text, None) == text
    assert provenance_text(text, None) == text
    assert time.perf_counter() - start < 5.0


def test_url_after_a_long_run_of_scheme_characters_is_still_found():
    assert 'pw' not in bounded_text('a' * 100 + '://user:pw@example.org/x')


SENSITIVE_KEYS = ['password', 'PassWord', 'user_passwd', 'ssh_passphrase', 'pwd', 'db-pwd', 'userPwd', 'client_secret',
                  'SECRET_KEY', 'aws_secret_access_key', 'credentials', 'Authorization', 'authentication', 'auth',
                  'http_auth', 'basicAuth', 'oauth', 'auth_value', 'bearer', 'bearer_value', 'api_key', 'apikey',
                  'API-KEY', 'x-api-key', 'private_key', 'privateKey', 'session_cookie', 'token', 'auth_token',
                  'authToken', 'access_token', 'refresh-token', 'sessiontoken', 'token_value']
BENIGN_KEYS = ['max_tokens', 'tokenizer', 'tokens', 'author', 'authors', 'authenticated', 'key', 'noise_level',
               'savepath', 'keyword']


@pytest.mark.parametrize('key', SENSITIVE_KEYS)
def test_sensitive_option_names_are_redacted_at_any_depth(key):
    secret = {'value': None, 'reason': 'redacted_sensitive_option'}
    assert safe_metadata({key: 'do-not-copy'}) == {key: secret}
    assert safe_metadata({'outer': [{'inner': {key: 'do-not-copy'}}]}) == {'outer': [{'inner': {key: secret}}]}


@pytest.mark.parametrize('key', BENIGN_KEYS)
def test_ordinary_option_names_are_kept(key):
    assert safe_metadata({key: 3}) == {key: 3}


def finish(fun, x0):
    x = np.ones_like(x0)
    fun(x)
    return x


def stay(fun, x0):
    fun(x0)
    return x0


def test_report_files_hold_no_credential_of_any_url_scheme(tmp_path):
    report_path = tmp_path / 'report.json'
    names = [CREDENTIAL_URLS[1], CREDENTIAL_URLS[3]]
    problem = Problem(lambda x: float(x @ x), [2.0], name=f'problem from {CREDENTIAL_URLS[6]}',
                      reference=ProblemReference(0.0, 'feasible', CREDENTIAL_URLS[7]))
    benchmark([finish, stay], problem=problem, solver_names=names, score_only=True, silent=True,
              n_jobs=1, max_eval_factor=2, report_path=str(report_path))
    files = [path for path in tmp_path.rglob('*') if path.is_file()]
    assert report_path in files and len(files) >= 2
    for path in files:
        assert CANARY.encode() not in path.read_bytes(), path
    assert_valid(json.loads(report_path.read_text(encoding='utf-8')), schema_resource('eval_report', 2))


def test_saved_report_and_artifact_names_hold_no_credential(tmp_path):
    report_path = tmp_path / 'report.json'
    problem = Problem(lambda x: float(x @ x), [2.0], name=f'problem from {CREDENTIAL_URLS[6]}')
    benchmark([finish, stay], problem=problem, solver_names=[f'ssh://u:{CANARY}@h', 'b'], silent=True, n_jobs=1,
              max_eval_factor=2, draw_hist_plots='sequential', savepath=str(tmp_path / 'out'),
              report_path=str(report_path))
    names = [path.name for path in (tmp_path / 'out').rglob('*')]
    assert any(name.endswith('.pdf') for name in names)
    assert not any(CANARY in name for name in names), names
    for path in (report_path, tmp_path / 'report.plot_data.json'):
        assert CANARY.encode() not in path.read_bytes(), path


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
