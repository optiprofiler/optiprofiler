"""Legacy single-feature behaviour pinned against the exact base commit."""

import json
from pathlib import Path

import pytest

from optiprofiler.tests import composition_goldens

FIXTURE = Path(__file__).parent / 'fixtures' / 'composition' / 'legacy_single_feature_goldens.json'
BASE_SHA = '4f821d8658700b8f3ec0bd480962310a74e47667'


@pytest.fixture(scope='module')
def golden_document():
    with open(FIXTURE, encoding='utf-8') as stream:
        return json.load(stream)


def test_fixture_provenance(golden_document):
    provenance = golden_document['provenance']
    assert provenance['source_sha'] == BASE_SHA
    assert provenance['fixture_version'] == composition_goldens.FIXTURE_VERSION


def test_all_single_features_match_base(golden_document):
    expected = golden_document['scenarios']
    actual = composition_goldens.run_all_scenarios()
    assert set(actual) == set(expected)
    composition_goldens.assert_matches(actual, expected)
