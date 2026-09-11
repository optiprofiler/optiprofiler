"""
Trusted compatibility boundary for historical serialized configurations.

Historical class paths (``optiprofiler.utils.FeatureOption`` values including
the removed ``n_runs`` member, the 1.x ``optiprofiler.opclasses.Feature``
layouts and the bridge ``optiprofiler.composition.ComposedFeature``/``Stage``)
decode to compatibility containers at the boundary; ``import_legacy_feature``
extracts the experiment count outside the canonical ``Feature``. Every fixture
is a frozen byte string with a pinned digest: the ac4 native fixtures
(``feature-single.pkl``, ``feature-composite.pkl``, ``enum-key-and-value.pkl``)
and the b93 single/composite pickles. Trusted input only: this is not a
safety promise for arbitrary pickles.
"""

import base64
import hashlib
import json
import pickle

import numpy as np
import pytest

from optiprofiler import benchmark
from optiprofiler.feature_definitions import StageRecord
from optiprofiler.legacy_compat import (LegacyComposedFeature, LegacyConfigurationError, LegacyEnumValue,
                                        LegacyFeature, LegacyStage, import_legacy_feature, load_options,
                                        loads_trusted, replay_arguments)
from optiprofiler.opclasses import Feature, FeaturedProblem, Problem
from optiprofiler.provenance import describe_feature
from optiprofiler.utils import FeatureOption

NOISY_DEFAULTS = {'noise_mode': 'random', 'distribution': 'gaussian', 'noise_map': 'chebyshev',
                  'noise_level': 0.001, 'noise_type': 'mixed'}
TRUNCATED_DEFAULTS = {'perturbed_trailing_digits': False, 'significant_digits': 6}

# Genuine ac4 native-serialization fixtures (clean ac4a67e, Python 3.10.12).
AC4_ENUM_KEY_AND_VALUE = (
    'gASVSgAAAAAAAAB9lCiMEm9wdGlwcm9maWxlci51dGlsc5SMDUZlYXR1cmVPcHRpb26Uk5SMBm5fcnVuc5SFlFKUSwOMCmVu'
    'dW1fdmFsdWWUaAZ1Lg==')
AC4_SINGLE = (
    'gASVLwEAAAAAAACMFm9wdGlwcm9maWxlci5vcGNsYXNzZXOUjAdGZWF0dXJllJOUKYGUfZQojA5fZGVjbGFyZWRfbmFtZZSM'
    'BW5vaXN5lIwFX25hbWWUaAaMDl9kZWNsYXJlZF9zcGVjlF2UfZQojARuYW1llGgGjAdvcHRpb25zlH2UjAtub2lzZV9sZXZl'
    'bJRHP+AAAAAAAABzdWGMB19zdGFnZXOUTowIX29wdGlvbnOUfZQojAZuX3J1bnOUSwKMC25vaXNlX2xldmVslEc/4AAAAAAA'
    'AIwKbm9pc2VfbW9kZZSMBnJhbmRvbZSMDGRpc3RyaWJ1dGlvbpSMCGdhdXNzaWFulIwJbm9pc2VfbWFwlIwJY2hlYnlzaGV2'
    'lIwKbm9pc2VfdHlwZZSMBW1peGVklHV1Yi4=')
AC4_COMPOSITE = (
    'gASVWQIAAAAAAACMGG9wdGlwcm9maWxlci5jb21wb3NpdGlvbpSMD0NvbXBvc2VkRmVhdHVyZZSTlCmBlH2UKIwOX2RlY2xh'
    'cmVkX3NwZWOUXZQofZQojARuYW1llIwFbm9pc3mUjAdvcHRpb25zlH2UdX2UKGgIjAl0cnVuY2F0ZWSUaAp9lHVljAZfcm91'
    'dGWUjAdmZWF0dXJllIwOX2RlY2xhcmVkX25hbWWUjA9ub2lzeSt0cnVuY2F0ZWSUjAVfbmFtZZSMD25vaXN5K3RydW5jYXRl'
    'ZJSMCF9vcHRpb25zlH2UjAZuX3J1bnOUSwNzjAdfc3RhZ2VzlGgAjAVTdGFnZZSTlCmBlE59lCiMCHBvc2l0aW9ulEsAaAho'
    'CYwEY29kZZRLAowKb2NjdXJyZW5jZZRLAGgQjBZvcHRpcHJvZmlsZXIub3BjbGFzc2VzlIwHRmVhdHVyZZSTlCmBlH2UKGgT'
    'aAloEWgJaBhOaBV9lCiMCm5vaXNlX21vZGWUjAZyYW5kb22UjAxkaXN0cmlidXRpb26UjAhnYXVzc2lhbpSMCW5vaXNlX21h'
    'cJSMCWNoZWJ5c2hldpSMC25vaXNlX2xldmVslEc/UGJN0vGp/IwKbm9pc2VfdHlwZZSMBW1peGVklHV1YnWGlGJoGimBlE59'
    'lChoHUsBaAhoDWgeSwNoH0sAaBBoIimBlH2UKGgTaA1oEWgNaBhOaBV9lCiMGXBlcnR1cmJlZF90cmFpbGluZ19kaWdpdHOU'
    'iYwSc2lnbmlmaWNhbnRfZGlnaXRzlEsGdXVidYaUYoaUdWIu')
# ``Feature('noisy+truncated', n_runs=3)`` and ``Feature('noisy', n_runs=2,
# noise_level=0.5)`` pickled (protocol 4) by the accepted candidate b93f8fb,
# whose stage dictionaries still carried their own legacy run counts (5 and 1).
B93_COMPOSITE = (
    'gASVMQIAAAAAAACMGG9wdGlwcm9maWxlci5jb21wb3NpdGlvbpSMD0NvbXBvc2VkRmVhdHVyZZSTlCmBlH2UKIwOX2RlY2xh'
    'cmVkX25hbWWUjA9ub2lzeSt0cnVuY2F0ZWSUjAVfbmFtZZSMD25vaXN5K3RydW5jYXRlZJSMCF9vcHRpb25zlH2UjAZuX3J1'
    'bnOUSwNzjAdfc3RhZ2VzlGgAjAVTdGFnZZSTlCmBlE59lCiMCHBvc2l0aW9ulEsAjARuYW1llIwFbm9pc3mUjARjb2RllEsC'
    'jApvY2N1cnJlbmNllEsAjAdmZWF0dXJllIwWb3B0aXByb2ZpbGVyLm9wY2xhc3Nlc5SMB0ZlYXR1cmWUk5QpgZR9lChoB4wF'
    'bm9pc3mUaAVoHGgMTmgJfZQojApub2lzZV9tb2RllIwGcmFuZG9tlIwMZGlzdHJpYnV0aW9ulIwIZ2F1c3NpYW6UjAlub2lz'
    'ZV9tYXCUjAljaGVieXNoZXaUaAtLBYwLbm9pc2VfbGV2ZWyURz9QYk3S8an8jApub2lzZV90eXBllIwFbWl4ZWSUdXVidYaU'
    'YmgOKYGUTn2UKGgRSwFoEowJdHJ1bmNhdGVklGgUSwNoFUsAaBZoGSmBlH2UKGgHjAl0cnVuY2F0ZWSUaAVoLWgMTmgJfZQo'
    'jBlwZXJ0dXJiZWRfdHJhaWxpbmdfZGlnaXRzlIloC0sBjBJzaWduaWZpY2FudF9kaWdpdHOUSwZ1dWJ1hpRihpR1Yi4=')
B93_SINGLE = (
    'gASV6gAAAAAAAACMFm9wdGlwcm9maWxlci5vcGNsYXNzZXOUjAdGZWF0dXJllJOUKYGUfZQojAVfbmFtZZSMBW5vaXN5lIwO'
    'X2RlY2xhcmVkX25hbWWUaAaMB19zdGFnZXOUTowIX29wdGlvbnOUfZQojAZuX3J1bnOUSwKMC25vaXNlX2xldmVslEc/4AAA'
    'AAAAAIwKbm9pc2VfbW9kZZSMBnJhbmRvbZSMDGRpc3RyaWJ1dGlvbpSMCGdhdXNzaWFulIwJbm9pc2VfbWFwlIwJY2hlYnlz'
    'aGV2lIwKbm9pc2VfdHlwZZSMBW1peGVklHV1Yi4=')
# Genuine ac4 identity declarations (``Feature(name, n_runs=2)`` for the three
# names below, pickled by the clean ac4a67e source on syu; captured by the
# controller, digests pinned).
AC4_IDENTITY_PLAIN = (
    'gASVqAAAAAAAAACMFm9wdGlwcm9maWxlci5vcGNsYXNzZXOUjAdGZWF0dXJllJOUKYGUfZQojA5fZGVjbGFyZWRfbmFtZZSM'
    'BXBsYWlulIwFX25hbWWUjAVwbGFpbpSMDl9kZWNsYXJlZF9zcGVjlF2UfZQojARuYW1llGgGjAdvcHRpb25zlH2UdWGMB19z'
    'dGFnZXOUTowIX29wdGlvbnOUfZSMBm5fcnVuc5RLAnN1Yi4=')
AC4_IDENTITY_PLAIN_PLAIN = (
    'gASVxwAAAAAAAACMFm9wdGlwcm9maWxlci5vcGNsYXNzZXOUjAdGZWF0dXJllJOUKYGUfZQojA5fZGVjbGFyZWRfbmFtZZSM'
    'C3BsYWluK3BsYWlulIwFX25hbWWUjAVwbGFpbpSMDl9kZWNsYXJlZF9zcGVjlF2UKH2UKIwEbmFtZZSMBXBsYWlulIwHb3B0'
    'aW9uc5R9lHV9lChoDIwFcGxhaW6UaA59lHVljAdfc3RhZ2VzlE6MCF9vcHRpb25zlH2UjAZuX3J1bnOUSwJzdWIu')
AC4_PLAIN_NOISY_PLAIN = (
    'gASVVAEAAAAAAACMFm9wdGlwcm9maWxlci5vcGNsYXNzZXOUjAdGZWF0dXJllJOUKYGUfZQojA5fZGVjbGFyZWRfbmFtZZSM'
    'EXBsYWluK25vaXN5K3BsYWlulIwFX25hbWWUjAVub2lzeZSMDl9kZWNsYXJlZF9zcGVjlF2UKH2UKIwEbmFtZZSMBXBsYWlu'
    'lIwHb3B0aW9uc5R9lHV9lChoDIwFbm9pc3mUaA59lHV9lChoDIwFcGxhaW6UaA59lHVljAdfc3RhZ2VzlE6MCF9vcHRpb25z'
    'lH2UKIwGbl9ydW5zlEsCjApub2lzZV9tb2RllIwGcmFuZG9tlIwMZGlzdHJpYnV0aW9ulIwIZ2F1c3NpYW6UjAlub2lzZV9t'
    'YXCUjAljaGVieXNoZXaUjAtub2lzZV9sZXZlbJRHP1BiTdLxqfyMCm5vaXNlX3R5cGWUjAVtaXhlZJR1dWIu')
DIGESTS = {
    'AC4_IDENTITY_PLAIN': '1c90257ab73d2b996a802fe43bc193cb5944c9765bf1a98e71106613e521dd69',
    'AC4_IDENTITY_PLAIN_PLAIN': 'd517c05a36b2fdfdf59528c994838627a25f5526f7edcb0b87163cf3aee3cc9a',
    'AC4_PLAIN_NOISY_PLAIN': '8985459e5b64512fda8c5466d9c29e552a946007513a08cc51aedce7f4d2ee04',
    'AC4_ENUM_KEY_AND_VALUE': 'be4aa8b7cfc0bd2823a0420c1813ed98e0db9b2348a4f6334a8b5dab5a521ffe',
    'AC4_SINGLE': 'e2e82ea875ac56308396c1d968207c19c8ec2495a975255ce39eb8d81c8d1ad0',
    'AC4_COMPOSITE': 'ccc3afc05e4490c0262efeaa46e5338c9ab8ccf915fdf4eee8849909d5da1cf2',
    'B93_COMPOSITE': '4f16a434d795f8a24e6bb49d2dd3f7791a1ffaf3ed5e406a87749503e36e8629',
    'B93_SINGLE': 'af4d93e82ab3391f693de2e8d747ea147939864fd8c1005603cc7253bb945729',
}


def fixture(name):
    raw = base64.b64decode(globals()[name])
    assert hashlib.sha256(raw).hexdigest() == DIGESTS[name]
    return raw


def sphere(x):
    return float(np.dot(x, x))


def mod_fun_plus_one(x, rng, problem):
    return problem.fun(x) + 1.0


def solver_stay(fun, x0):
    fun(x0)
    return x0


def solver_zero(fun, x0):
    fun(x0)
    return np.zeros_like(x0)


def local_options(feature):
    return [dict(stage.options) for stage in feature.stages]


class TestTrustedEnumDecoding:

    def test_plain_pickle_cannot_decode_the_removed_member(self):
        with pytest.raises(ValueError, match="'n_runs' is not a valid FeatureOption"):
            pickle.loads(fixture('AC4_ENUM_KEY_AND_VALUE'))

    def test_removed_member_decodes_to_a_string_valued_marker(self):
        decoded = loads_trusted(fixture('AC4_ENUM_KEY_AND_VALUE'))
        assert decoded['n_runs'] == 3
        marker = decoded['enum_value']
        assert isinstance(marker, LegacyEnumValue) and isinstance(marker, str)
        assert (marker.enum_name, marker.value, marker.name) == ('FeatureOption', 'n_runs', 'N_RUNS')
        assert marker == 'n_runs' and hash(marker) == hash('n_runs')
        assert [key for key in decoded if key == 'n_runs'] == ['n_runs']
        # The marker survives ordinary transport and never becomes an active member.
        clone = pickle.loads(pickle.dumps(decoded))
        assert clone['enum_value'] == 'n_runs' and isinstance(clone['enum_value'], LegacyEnumValue)
        assert 'n_runs' not in FeatureOption.__members__.values()

    def test_members_that_still_exist_decode_to_the_active_member(self):
        payload = pickle.dumps({FeatureOption.NOISE_LEVEL: 0.5})
        decoded = loads_trusted(payload)
        assert decoded == {'noise_level': 0.5}
        assert next(iter(decoded)) is FeatureOption.NOISE_LEVEL


class TestHistoricalFeatureImport:

    def test_ac4_single_layout(self):
        container = loads_trusted(fixture('AC4_SINGLE'))
        assert isinstance(container, LegacyFeature) and not isinstance(container, Feature)
        assert container.class_path == 'optiprofiler.opclasses.Feature'
        imported = import_legacy_feature(container)
        assert isinstance(imported.feature, Feature)
        assert [stage.identity for stage in imported.feature.stages] == ['noisy#0']
        assert local_options(imported.feature) == [{**NOISY_DEFAULTS, 'noise_level': 0.5}]
        assert (imported.n_runs, imported.declared_name) == (2, 'noisy')
        assert imported.declared_spec == [{'name': 'noisy', 'options': {'noise_level': 0.5}}]
        assert imported.source == 'optiprofiler.opclasses.Feature'
        # The count lives outside the Feature. The declaration the ac4 object
        # recorded is restored as the Feature's own declaration (shorthand route).
        assert imported.feature.declared.route == 'feature_name'
        assert imported.feature.declared.entries == (('noisy', {'noise_level': 0.5}),)
        assert imported.feature.declared_name == 'noisy' and imported.feature.name == 'noisy'
        assert all('n_runs' not in options for options in local_options(imported.feature))

    def test_ac4_composite_layout(self):
        container = loads_trusted(fixture('AC4_COMPOSITE'))
        assert isinstance(container, LegacyComposedFeature)
        assert all(isinstance(stage, LegacyStage) for stage in container.state['_stages'])
        imported = import_legacy_feature(container)
        assert [stage.identity for stage in imported.feature.stages] == ['noisy#0', 'truncated#0']
        assert local_options(imported.feature) == [NOISY_DEFAULTS, TRUNCATED_DEFAULTS]
        assert (imported.n_runs, imported.declared_name) == (3, 'noisy+truncated')
        assert imported.declared_spec == [{'name': 'noisy', 'options': {}}, {'name': 'truncated', 'options': {}}]
        assert imported.source == 'optiprofiler.composition.ComposedFeature'
        assert imported.feature.declared.route == 'feature' and imported.feature.declared_name == 'noisy+truncated'
        assert imported.feature.declared.entries == (('noisy', {}), ('truncated', {}))

    def test_b93_composite_root_count_wins_over_child_residue(self):
        container = loads_trusted(fixture('B93_COMPOSITE'))
        children = [stage.state['feature'].state['_options'].get('n_runs') for stage in container.state['_stages']]
        assert children == [5, 1]
        imported = import_legacy_feature(container)
        assert imported.n_runs == 3
        assert local_options(imported.feature) == [NOISY_DEFAULTS, TRUNCATED_DEFAULTS]
        assert imported.declared_spec is None and imported.declared_name == 'noisy+truncated'
        problem = Problem(sphere, np.array([1.0, 2.0]))
        x = np.array([0.3, -0.7])
        fresh = Feature('noisy+truncated')
        assert FeaturedProblem(problem, imported.feature, 5, 7).fun(x) == FeaturedProblem(problem, fresh, 5, 7).fun(x)
        assert describe_feature(imported.feature)['stages'] == describe_feature(fresh)['stages']

    def test_b93_single_layout(self):
        imported = import_legacy_feature(loads_trusted(fixture('B93_SINGLE')))
        assert imported.n_runs == 2
        assert local_options(imported.feature) == [{**NOISY_DEFAULTS, 'noise_level': 0.5}]
        assert imported.declared_spec is None and imported.declared_name == 'noisy'

    def test_public_flow_from_old_pickle_to_benchmark_and_report(self, tmp_path):
        imported = import_legacy_feature(loads_trusted(fixture('B93_COMPOSITE')))
        feature = imported.feature
        assert (feature.declared.route, feature.declared.entries, feature.declared_name) == (None, (), None)
        assert [stage.identity for stage in feature.stages] == ['noisy#0', 'truncated#0']
        assert imported.n_runs == 3
        # The current pickle form keeps both facts: effective stages and unknown declaration.
        clone = pickle.loads(pickle.dumps(feature))
        assert clone.declared.route is None and describe_feature(clone) == describe_feature(feature)
        assert local_options(clone) == [NOISY_DEFAULTS, TRUNCATED_DEFAULTS]
        problem = Problem(sphere, np.array([1.0, -2.0]), name='SPHERE')
        common = dict(problem=problem, score_only=True, silent=True, draw_hist_plots='none', savepath=str(tmp_path))
        target = tmp_path / 'imported-report.json'
        scores, _, _ = benchmark([solver_stay, solver_zero], feature=clone, n_runs=imported.n_runs,
                                 report_path=target, **common)
        fresh, _, _ = benchmark([solver_stay, solver_zero], feature_name='noisy+truncated', n_runs=3, **common)
        np.testing.assert_array_equal(scores, fresh)
        with open(target, encoding='utf-8') as stream:
            report = json.load(stream)
        block = report['configuration']['effective']['feature']
        assert (block['declaration_route'], block['declared'], block['declared_name']) == (None, [], None)
        assert block['route'] == 'feature' and block['name'] == 'noisy+truncated'
        assert [stage['identity'] for stage in block['stages']] == ['noisy#0', 'truncated#0']
        assert block['stages'][0]['options'] == NOISY_DEFAULTS and block['stages'][1]['options'] == TRUNCATED_DEFAULTS
        assert report['configuration']['effective']['experiment']['primary']['n_runs'] == 3
        assert 'n_runs' not in json.dumps(block)

    def test_public_flow_keeps_a_recorded_declaration(self, tmp_path):
        imported = import_legacy_feature(loads_trusted(fixture('AC4_SINGLE')))
        clone = pickle.loads(pickle.dumps(imported.feature))
        assert clone.declared.route == 'feature_name' and clone.declared.entries == (('noisy', {'noise_level': 0.5}),)
        assert local_options(clone) == [{**NOISY_DEFAULTS, 'noise_level': 0.5}]
        target = tmp_path / 'ac4-report.json'
        common = dict(problem=Problem(sphere, np.array([1.0, -2.0]), name='SPHERE'), score_only=True, silent=True,
                      draw_hist_plots='none', savepath=str(tmp_path))
        scores, _, _ = benchmark([solver_stay, solver_zero], feature=clone, n_runs=imported.n_runs, report_path=target,
                                 **common)
        fresh, _, _ = benchmark([solver_stay, solver_zero], feature_name='noisy', noise_level=0.5, n_runs=2, **common)
        np.testing.assert_array_equal(scores, fresh)
        with open(target, encoding='utf-8') as stream:
            block = json.load(stream)['configuration']['effective']['feature']
        assert (block['route'], block['declaration_route'], block['declared_name']) == ('feature', 'feature_name', 'noisy')
        assert block['declared'] == [{'name': 'noisy', 'options': {'noise_level': 0.5}}]
        assert block['stages'][0]['options'] == {**NOISY_DEFAULTS, 'noise_level': 0.5}

    def test_recorded_declaration_that_does_not_describe_the_stages_is_rejected(self):
        container = LegacyFeature()
        container.__setstate__({'_name': 'noisy', '_options': {'n_runs': 2, 'noise_level': 0.5},
                                '_declared_name': 'truncated', '_declared_spec': [{'name': 'truncated', 'options': {}}]})
        with pytest.raises(LegacyConfigurationError, match='does not describe'):
            import_legacy_feature(container)
        container.__setstate__({'_name': 'noisy', '_options': {'noise_level': 0.5}, '_declared_spec': 'noisy'})
        with pytest.raises(LegacyConfigurationError, match='not a list'):
            import_legacy_feature(container)

    @pytest.mark.parametrize('name, effective', [
        ('AC4_IDENTITY_PLAIN', []), ('AC4_IDENTITY_PLAIN_PLAIN', []), ('AC4_PLAIN_NOISY_PLAIN', ['noisy'])])
    def test_genuine_identity_declarations_survive_import(self, name, effective):
        imported = import_legacy_feature(loads_trusted(fixture(name)))
        assert imported.n_runs == 2
        feature = pickle.loads(pickle.dumps(imported.feature, protocol=4))
        assert isinstance(feature, Feature)
        assert [stage.name for stage in feature.stages] == effective
        assert feature.is_identity == (not effective)
        declared_name = imported.declared_name
        assert feature.declared.route == 'feature_name' and feature.declared_name == declared_name
        assert [stage for stage, _ in feature.declared.entries] == declared_name.split('+')
        assert all(not dict(options) for _, options in feature.declared.entries)
        assert describe_feature(feature)['declared'] == [{'name': token, 'options': {}} for token in declared_name.split('+')]

    def test_identity_object_options_are_validated_not_ignored(self):
        container = LegacyFeature()
        container.__setstate__({'_name': 'plain', '_options': {'n_runs': 2, 'noise_level': 0.5},
                                '_declared_name': 'plain', '_declared_spec': [{'name': 'plain', 'options': {}}]})
        with pytest.raises(LegacyConfigurationError, match="stage 'plain'"):
            import_legacy_feature(container)
        container.__setstate__({'_name': 'noisy', '_options': {'n_runs': 2, 'noise_level': 0.5},
                                '_declared_name': 'plain+noisy',
                                '_declared_spec': [{'name': 'plain', 'options': {'noise_level': 0.5}}, {'name': 'noisy', 'options': {}}]})
        with pytest.raises(LegacyConfigurationError, match="stage 'plain'"):
            import_legacy_feature(container)
        container.__setstate__({'_name': 'noisy', '_options': {'n_runs': 2},
                                '_declared_spec': [{'name': 'plain', 'options': {}}, {'name': 'noisy', 'options': {'noise_level': -1}}]})
        with pytest.raises(LegacyConfigurationError, match="stage 'noisy'"):
            import_legacy_feature(container)

    def test_pre_pipeline_layout_without_declaration(self):
        container = LegacyFeature()
        container.__setstate__({'_name': 'permuted', '_options': {'n_runs': 1}})
        imported = import_legacy_feature(container)
        assert imported.feature.name == 'permuted' and imported.n_runs == 1
        assert (imported.declared_name, imported.declared_spec) == (None, None)
        without_count = LegacyFeature()
        without_count.__setstate__({'_name': 'plain', '_options': {}})
        imported = import_legacy_feature(without_count)
        assert imported.feature.is_identity and imported.n_runs is None

    def test_invalid_historical_options_fail_explicitly(self):
        container = LegacyFeature()
        container.__setstate__({'_name': 'noisy', '_options': {'n_runs': 2, 'noise_level': -1.0}})
        with pytest.raises(LegacyConfigurationError, match='noise_level'):
            import_legacy_feature(container)
        with pytest.raises(TypeError, match='legacy'):
            import_legacy_feature(Feature('noisy'))

    def test_containers_never_reach_normal_construction(self):
        container = loads_trusted(fixture('B93_SINGLE'))
        with pytest.raises(TypeError, match='import_legacy_feature'):
            Feature(container)
        with pytest.raises(TypeError, match='import_legacy_feature'):
            benchmark([solver_stay, solver_zero], feature=container, problem=Problem(sphere, np.array([1.0])),
                      silent=True, draw_hist_plots='none', score_only=True)


class TestCurrentPickles:

    def test_current_feature_pickles_are_rebuilt_specifications(self):
        feature = Feature('custom+noisy', mod_fun=mod_fun_plus_one, noise_level=0.5)
        for loads in (pickle.loads, loads_trusted):
            clone = loads(pickle.dumps(feature))
            assert isinstance(clone, Feature)
            assert describe_feature(clone) == describe_feature(feature)
            assert clone.declared.route == 'feature_name'
            assert clone.stages[0].options['mod_fun'] is mod_fun_plus_one
            assert all(isinstance(stage, StageRecord) for stage in clone.stages)
        structured = Feature([{'name': 'plain'}, {'name': 'truncated', 'options': {'significant_digits': 3}}])
        clone = loads_trusted(pickle.dumps(structured))
        assert describe_feature(clone) == describe_feature(structured)
        assert clone.declared.entries == structured.declared.entries
        identity = pickle.loads(pickle.dumps(Feature('plain+plain')))
        assert identity.is_identity and identity.declared_name == 'plain+plain'

    def test_native_form_carries_effective_stages_and_declaration_separately(self):
        from optiprofiler.opclasses import FEATURE_NATIVE_VERSION, _rebuild_feature
        feature = Feature('noisy')
        rebuild, (version, route, declared, effective) = feature.__reduce__()
        assert rebuild is _rebuild_feature and version == FEATURE_NATIVE_VERSION == 1
        assert (route, declared) == ('feature_name', (('noisy', {}),))
        # Every validated local option is explicit in the native form, so the
        # numbers never depend on the defaults of the version that loads it.
        assert effective == (('noisy', NOISY_DEFAULTS),)
        assert 'n_runs' not in effective[0][1]
        changed = _rebuild_feature(version, route, declared, (('noisy', {**NOISY_DEFAULTS, 'noise_level': 0.5}),))
        assert changed.stages[0].options['noise_level'] == 0.5 and changed.declared.entries == (('noisy', {}),)
        with pytest.raises(ValueError, match='Unsupported native Feature form version'):
            _rebuild_feature(2, route, declared, effective)
        with pytest.raises(ValueError):
            _rebuild_feature(version, route, declared, (('noisy', {'noise_level': -1.0}),))


class TestRefinedOptionsReplay:

    def test_refined_v2_replays_the_experiment(self, tmp_path):
        kwargs = dict(problem=Problem(sphere, np.array([1.0, -2.0]), name='SPHERE'), silent=True,
                      draw_hist_plots='none', savepath=str(tmp_path), benchmark_id='refined', n_jobs=1)
        spec = [{'name': 'noisy', 'options': {'noise_level': 0.5}}, {'name': 'truncated'}]
        scores, _, _ = benchmark([solver_stay, solver_zero], feature=spec, n_runs=2, **kwargs)
        paths = list(tmp_path.rglob('options_refined.pkl'))
        assert len(paths) == 1
        refined = load_options(paths[0])
        assert refined['schema'] == 'options_refined-v2' and refined['n_runs'] == 2
        replay = replay_arguments(refined)
        assert replay['n_runs'] == 2
        assert [entry['name'] for entry in replay['feature']] == ['noisy', 'truncated']
        assert replay['feature'][0]['options']['noise_level'] == 0.5
        assert 'feature_name' not in replay and 'feature_route' not in replay
        again, _, _ = benchmark([solver_stay, solver_zero], feature=replay['feature'], n_runs=replay['n_runs'],
                                **{**kwargs, 'benchmark_id': 'replay'})
        np.testing.assert_array_equal(again, scores)

    def test_bridge_layout_with_specification_but_flat_count(self):
        options = {'feature_route': 'feature_name', 'feature_name': 'noisy', 'n_runs': 4, 'noise_level': 0.5,
                   'feature_specification': [{'name': 'noisy', 'options': {'noise_level': 0.5}}], 'ptype': 'u'}
        replay = replay_arguments(options)
        assert replay == {'feature': [{'name': 'noisy', 'options': {'noise_level': 0.5}}], 'n_runs': 4,
                          'problem_options': {'ptype': 'u'}}

    def test_flat_layout_needs_an_explicit_identity(self):
        flat = {'n_runs': 5, 'noise_level': 0.5, 'ptype': 'u', 'mindim': 2}
        with pytest.raises(LegacyConfigurationError, match='feature identity'):
            replay_arguments(flat)
        replay = replay_arguments(flat, feature_name='noisy')
        assert replay['n_runs'] == 5 and replay['problem_options'] == {'ptype': 'u', 'mindim': 2}
        assert Feature(replay['feature']).stages[0].options['noise_level'] == 0.5
        assert 'n_runs' not in replay['feature'][0]['options']

    def test_enum_keyed_options_decode(self, tmp_path):
        path = tmp_path / 'options_user.pkl'
        path.write_bytes(fixture('AC4_ENUM_KEY_AND_VALUE'))
        options = load_options(path)
        assert options['n_runs'] == 3
