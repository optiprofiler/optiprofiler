"""Legacy single-feature behaviour pinned against the exact base commit."""

import json
from pathlib import Path

import pytest

from optiprofiler.tests import composition_goldens

FIXTURE = Path(__file__).parent / 'fixtures' / 'composition' / 'legacy_single_feature_goldens.json'
BASE_SHA = '4f821d8658700b8f3ec0bd480962310a74e47667'
# SHA-256 of the generator that wrote the fixture (composition_goldens.py at
# commit 269358a, run against BASE_SHA). The shipped generator has evolved
# since (it imports the experiment layer that BASE_SHA does not have), so it
# cannot regenerate the fixture at BASE_SHA; the recorded hash pins the
# provenance and detects a silently edited fixture. See fixtures/composition/README.md.
GENERATOR_SHA256 = '2ba635ad3225f9205d4854d8f9f007ffeed8f9b4d74e73e0dc1e4a373be0ef1d'


@pytest.fixture(scope='module')
def golden_document():
    with open(FIXTURE, encoding='utf-8') as stream:
        return json.load(stream)


def test_fixture_provenance(golden_document):
    provenance = golden_document['provenance']
    assert provenance['source_sha'] == BASE_SHA
    assert provenance['fixture_version'] == composition_goldens.FIXTURE_VERSION
    assert provenance['scenario_module_sha256'] == GENERATOR_SHA256


def test_all_single_features_match_base(golden_document):
    expected = golden_document['scenarios']
    actual = composition_goldens.run_all_scenarios()
    assert set(actual) == set(expected)
    composition_goldens.assert_matches(actual, expected)
