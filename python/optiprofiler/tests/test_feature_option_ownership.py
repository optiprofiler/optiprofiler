"""
Ownership of feature options: the experiment-wide run count against the
options each stage owns.

``n_runs`` belongs to the experiment. A composition stores it once, on the
feature that the experiment runs; the stages of a composition own only their
local options and never carry a run count of their own. Every expectation is
a worked value; provenance is read through the public ``describe_pipeline``.
"""

import base64
import pickle

import numpy as np
import pytest

from optiprofiler.composition import describe_pipeline
from optiprofiler.opclasses import Feature, FeaturedProblem, Problem

NOISY_DEFAULTS = {'noise_mode': 'random', 'distribution': 'gaussian', 'noise_map': 'chebyshev',
                  'noise_level': 0.001, 'noise_type': 'mixed'}
TRUNCATED_DEFAULTS = {'perturbed_trailing_digits': False, 'significant_digits': 6}

# ``Feature('noisy+truncated', n_runs=3)`` and ``Feature('noisy', n_runs=2,
# noise_level=0.5)`` pickled (protocol 4) by the accepted candidate b93f8fb,
# whose stage dictionaries still carried their own legacy run counts.
B93_COMPOSITE = (
    'gASVMQIAAAAAAACMGG9wdGlwcm9maWxlci5jb21wb3NpdGlvbpSMD0NvbXBvc2VkRmVhdHVyZZSTlCmBlH2UKIwOX2RlY2xh'
    'cmVkX25hbWWUjA9ub2lzeSt0cnVuY2F0ZWSUjAVfbmFtZZSMD25vaXN5K3RydW5jYXRlZJSMCF9vcHRpb25zlH2UjAZuX3J1'
    'bnOUSwNzjAdfc3RhZ2VzlGgAjAVTdGFnZZSTlCmBlE59lCiMCHBvc2l0aW9ulEsAjARuYW1llIwFbm9pc3mUjARjb2RllEsC'
    'jApvY2N1cnJlbmNllEsAjAdmZWF0dXJllIwWb3B0aXByb2ZpbGVyLm9wY2xhc3Nlc5SMB0ZlYXR1cmWUk5QpgZR9lChoB4wF'
    'bm9pc3mUaAVoHGgMTmgJfZQojApub2lzZV9tb2RllIwGcmFuZG9tlIwMZGlzdHJpYnV0aW9ulIwIZ2F1c3NpYW6UjAlub2lz'
    'ZV9tYXCUjAljaGVieXNoZXaUaAtLBYwLbm9pc2VfbGV2ZWyURz9QYk3S8an8jApub2lzZV90eXBllIwFbWl4ZWSUdXVidYaU'
    'YmgOKYGUTn2UKGgRSwFoEowJdHJ1bmNhdGVklGgUSwNoFUsAaBZoGSmBlH2UKGgHjAl0cnVuY2F0ZWSUaAVoLWgMTmgJfZQo'
    'jBlwZXJ0dXJiZWRfdHJhaWxpbmdfZGlnaXRzlIloC0sBjBJzaWduaWZpY2FudF9kaWdpdHOUSwZ1dWJ1hpRihpR1Yi4='
)
B93_SINGLE = (
    'gASV6gAAAAAAAACMFm9wdGlwcm9maWxlci5vcGNsYXNzZXOUjAdGZWF0dXJllJOUKYGUfZQojAVfbmFtZZSMBW5vaXN5lIwO'
    'X2RlY2xhcmVkX25hbWWUaAaMB19zdGFnZXOUTowIX29wdGlvbnOUfZQojAZuX3J1bnOUSwKMC25vaXNlX2xldmVslEc/4AAA'
    'AAAAAIwKbm9pc2VfbW9kZZSMBnJhbmRvbZSMDGRpc3RyaWJ1dGlvbpSMCGdhdXNzaWFulIwJbm9pc2VfbWFwlIwJY2hlYnlz'
    'aGV2lIwKbm9pc2VfdHlwZZSMBW1peGVklHV1Yi4='
)


def mod_fun_plus_one(x, rng, problem):
    return problem.fun(x) + 1.0


def sphere(x):
    return float(np.dot(x, x))


def stage_options(feature):
    return [stage['options'] for stage in describe_pipeline(feature)['stages']]


class TestRunCountOwnership:

    def test_composite_reports_one_common_run_count_and_local_stage_options(self):
        feature = Feature('noisy+truncated', n_runs=3)
        assert feature.options == {'n_runs': 3}
        pipeline = describe_pipeline(feature)
        assert pipeline['schema'] == 'feature_pipeline-v2'
        assert pipeline['common_options'] == {'n_runs': 3}
        assert [stage['identity'] for stage in pipeline['stages']] == ['noisy#0', 'truncated#0']
        assert stage_options(feature) == [NOISY_DEFAULTS, TRUNCATED_DEFAULTS]

    def test_single_feature_reports_common_and_local_options_separately(self):
        feature = Feature('noisy', n_runs=2, noise_level=0.5)
        assert feature.options == {'n_runs': 2, **NOISY_DEFAULTS, 'noise_level': 0.5}
        pipeline = describe_pipeline(feature)
        assert pipeline['schema'] == 'feature_pipeline-v2'
        assert pipeline['seed_policy'] == 'legacy-run-seed'
        assert pipeline['common_options'] == {'n_runs': 2}
        assert stage_options(feature) == [{**NOISY_DEFAULTS, 'noise_level': 0.5}]

    def test_plain_owns_only_the_common_count(self):
        feature = Feature('plain')
        assert feature.options == {'n_runs': 1}
        pipeline = describe_pipeline(feature)
        assert pipeline['common_options'] == {'n_runs': 1}
        assert stage_options(feature) == [{}]

    @pytest.mark.parametrize('name, options, expected_runs, stochastic', [
        # custom keeps its established single run although it is stochastic.
        ('custom+truncated', {'mod_fun': mod_fun_plus_one}, 1, True),
        # An unrotated linear transformation keeps five runs although it is deterministic.
        ('linearly_transformed+quantized', {'rotated': False}, 5, False),
        ('noisy+truncated', {}, 5, True),
        ('truncated+quantized', {}, 1, False),
        ('truncated+quantized', {'perturbed_trailing_digits': True}, 5, True),
        ('noisy+truncated', {'noise_mode': 'deterministic'}, 1, False),
        ('noisy+perturbed_x0', {'noise_mode': 'deterministic'}, 5, True),
        ('plain+noisy+plain', {}, 5, True),
        ('unrelaxable_constraints+nonquantifiable_constraints', {}, 1, False),
    ])
    def test_default_count_is_the_largest_established_stage_default(self, name, options, expected_runs, stochastic):
        feature = Feature(name, **options)
        assert feature.options['n_runs'] == expected_runs
        assert feature.is_stochastic is stochastic
        assert describe_pipeline(feature)['common_options'] == {'n_runs': expected_runs}
        for options_of_stage in stage_options(feature):
            assert 'n_runs' not in options_of_stage

    def test_explicit_count_overrides_every_stage_default(self):
        feature = Feature('noisy+perturbed_x0+random_nan', n_runs=2)
        assert feature.options == {'n_runs': 2}
        assert describe_pipeline(feature)['common_options'] == {'n_runs': 2}
        assert all('n_runs' not in options_of_stage for options_of_stage in stage_options(feature))

    def test_returned_option_mappings_are_copies(self):
        feature = Feature('noisy+truncated', n_runs=3)
        feature.options['n_runs'] = 99
        describe_pipeline(feature)['common_options']['n_runs'] = 99
        assert feature.options == {'n_runs': 3}
        assert describe_pipeline(feature)['common_options'] == {'n_runs': 3}


class TestArchivedFeatureObjects:

    def test_fresh_pickle_round_trip_keeps_ownership(self):
        clone = pickle.loads(pickle.dumps(Feature('noisy+truncated', n_runs=3)))
        assert clone.options == {'n_runs': 3}
        assert describe_pipeline(clone)['common_options'] == {'n_runs': 3}
        assert stage_options(clone) == [NOISY_DEFAULTS, TRUNCATED_DEFAULTS]

    def test_b93_composite_pickle_runs_and_describes_local_options_only(self):
        old = pickle.loads(base64.b64decode(B93_COMPOSITE))
        assert old.name == 'noisy+truncated' and old.options == {'n_runs': 3}
        pipeline = describe_pipeline(old)
        assert pipeline['schema'] == 'feature_pipeline-v2'
        assert pipeline['common_options'] == {'n_runs': 3}
        # The archived stage dictionaries carried their own legacy counts (5 and
        # 1); a fresh description shows local options only. The object recorded
        # no declaration, which is reported as unknown, never fabricated.
        assert stage_options(old) == [NOISY_DEFAULTS, TRUNCATED_DEFAULTS]
        assert pipeline['declared_spec'] is None and pipeline['route'] == 'feature_name'
        problem = Problem(sphere, np.array([1.0, 2.0]))
        x = np.array([0.3, -0.7])
        fresh = Feature('noisy+truncated', n_runs=3)
        assert FeaturedProblem(problem, old, 5, 7).fun(x) == FeaturedProblem(problem, fresh, 5, 7).fun(x)

    def test_b93_single_pickle_is_unchanged(self):
        old = pickle.loads(base64.b64decode(B93_SINGLE))
        assert old.options == Feature('noisy', n_runs=2, noise_level=0.5).options
        assert describe_pipeline(old)['common_options'] == {'n_runs': 2}
        assert stage_options(old) == [{**NOISY_DEFAULTS, 'noise_level': 0.5}]
        assert describe_pipeline(old)['declared_spec'] is None
        assert describe_pipeline(Feature('noisy', n_runs=2, noise_level=0.5))['declared_spec'] == [
            {'name': 'noisy', 'options': {'noise_level': 0.5}}]
