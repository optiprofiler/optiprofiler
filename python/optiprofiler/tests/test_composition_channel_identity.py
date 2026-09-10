"""
Stream identity of genuine compositions (``seedsequence-v2``).

Every stage of a composition draws from streams seeded by
``SeedSequence(run_seed, spawn_key=(stage code, occurrence, channel tag))``
with the frozen tags fun=0, cub=1, ceq=2 and construction=3. Distinct
derivation identities remove the structural alias that omitting the channel
tag caused for channels with equal values, points and served-query indices;
32-bit seeds can still coincide by chance, and nothing here is a statistical
independence proof. Single and effective-single features keep the legacy
run-seed streams (pinned by the base goldens); this file only concerns
compositions.
"""

import json

import matplotlib
import numpy as np
import pytest

matplotlib.use('Agg')

from optiprofiler import benchmark
from optiprofiler.composition import CHANNEL_TAGS, SEED_POLICY, STAGE_CODES, Stage, stage_seed
from optiprofiler.loader import load_results_from_h5
from optiprofiler.opclasses import Feature, FeaturedProblem, Problem

X0 = np.array([0.7, -0.4])


def zero(x):
    return 0.0


def zero_cub(x):
    return np.array([0.0])


def zero_ceq(x):
    return np.array([0.0])


def constant_channels():
    """All three root channels return zero everywhere: equal values, points and indices."""
    return Problem(zero, X0, cub=zero_cub, ceq=zero_ceq)


def sequence_seed(run_seed, code, occurrence, tag):
    entropy = 0 if run_seed is None else int(run_seed)
    return int(np.random.SeedSequence(entropy, spawn_key=(code, occurrence, tag)).generate_state(1, dtype=np.uint32)[0])


class TestSeedDerivation:

    def test_frozen_tags_and_policy(self):
        assert CHANNEL_TAGS == {'fun': 0, 'cub': 1, 'ceq': 2, 'construction': 3}
        assert SEED_POLICY == 'seedsequence-v2'
        assert STAGE_CODES == {'perturbed_x0': 1, 'noisy': 2, 'truncated': 3, 'permuted': 4,
                               'linearly_transformed': 5, 'random_nan': 6, 'unrelaxable_constraints': 7,
                               'nonquantifiable_constraints': 8, 'quantized': 9, 'custom': 10}

    def test_seed_is_the_documented_conversion(self):
        for run_seed in (None, 0, 1, 2 ** 31):
            for name, code in STAGE_CODES.items():
                for occurrence in (0, 1, 3):
                    stage = Stage(occurrence, name, occurrence, Feature(name) if name != 'custom' else Feature('custom'))
                    for channel, tag in CHANNEL_TAGS.items():
                        assert stage_seed(run_seed, stage, channel) == sequence_seed(run_seed, code, occurrence, tag)

    def test_frozen_literal_seeds(self):
        noisy = Stage(0, 'noisy', 0, Feature('noisy'))
        literals = {channel: stage_seed(0, noisy, channel) for channel in CHANNEL_TAGS}
        # Frozen values, identical on NumPy 1.26.4 and 2.2.6 (SeedSequence is version-stable).
        assert literals == {'fun': 2103646603, 'cub': 3629157004, 'ceq': 2370119283, 'construction': 1055300565}
        assert stage_seed(12345, Stage(2, 'quantized', 1, Feature('quantized')), 'cub') == 1036917780

    def test_identities_are_distinct_and_position_independent(self):
        seeds = set()
        for name, code in STAGE_CODES.items():
            for occurrence in (0, 1):
                stage = Stage(0, name, occurrence, Feature(name))
                for channel in CHANNEL_TAGS:
                    seeds.add(stage_seed(0, stage, channel))
        assert len(seeds) == len(STAGE_CODES) * 2 * len(CHANNEL_TAGS)
        early = Stage(1, 'noisy', 1, Feature('noisy'))
        late = Stage(6, 'noisy', 1, Feature('noisy'))
        assert all(stage_seed(7, early, channel) == stage_seed(7, late, channel) for channel in CHANNEL_TAGS)
        assert stage_seed(0, early, 'fun') == stage_seed(None, early, 'fun')


class TestChannelsDoNotAlias:

    def test_equal_values_points_and_indices_give_different_samples(self):
        feature = Feature('noisy+truncated', noise_type='absolute', noise_level=1.0, significant_digits=12)
        featured = FeaturedProblem(constant_channels(), feature, 10, 3)
        for x in (X0, X0, X0 + 0.5):
            samples = (featured.fun(x), featured.cub(x)[0], featured.ceq(x)[0])
            assert len(set(samples)) == 3
            assert all(np.isfinite(samples))

    def test_adversarial_equal_counters_across_stages(self):
        # Two noisy stages, deterministic zero input everywhere: with the same
        # served index and the same point, the stages must still draw
        # different samples on every channel.
        feature = Feature('noisy+noisy', noise_type='absolute', noise_level=1.0)
        single = FeaturedProblem(constant_channels(), Feature('noisy+plain', noise_type='absolute', noise_level=1.0), 10, 3)
        double = FeaturedProblem(constant_channels(), feature, 10, 3)
        for x in (X0, X0):
            assert double.fun(x) != single.fun(x)
            assert double.fun(x) != 2.0 * single.fun(x)

    def test_construction_channel_drives_initialization(self):
        feature = Feature('perturbed_x0+noisy', perturbation_level=0.25, noise_level=0.0)
        featured = FeaturedProblem(constant_channels(), feature, 10, 9)
        stage = Stage(0, 'perturbed_x0', 0, Feature('perturbed_x0', perturbation_level=0.25))
        rng = Feature.get_default_rng(stage_seed(9, stage, 'construction'))
        direction = rng.standard_normal(2)
        expected = X0 + 0.25 * max(1.0, np.linalg.norm(X0)) * direction / np.linalg.norm(direction)
        np.testing.assert_allclose(featured.x0, expected, rtol=1e-12)


class TestPolicyProvenance:

    @staticmethod
    def stay(fun, x0):
        fun(x0)
        return x0

    @staticmethod
    def step(fun, x0):
        fun(x0 + 0.1)
        return x0 + 0.1

    def test_serial_and_parallel_runs_agree_and_record_v2(self, tmp_path):
        archives = []
        for n_jobs in (1, 2):
            benchmark([self.stay, self.step], feature_name='noisy+perturbed_x0', plibs=['s2mpj'], ptype='u', mindim=2,
                      maxdim=2, problem_names=['ROSENBR'], max_eval_factor=5, n_jobs=n_jobs, silent=True,
                      draw_hist_plots='none', savepath=str(tmp_path), benchmark_id=f'jobs{n_jobs}')
            archive = list((tmp_path / f'jobs{n_jobs}').rglob('data_for_loading.h5'))
            assert len(archive) == 1
            archives.append(load_results_from_h5(str(archive[0]))[0])
        for key in ('fun_histories', 'maxcv_histories', 'fun_outs', 'fun_inits'):
            np.testing.assert_array_equal(np.asarray(archives[0][key], dtype=float), np.asarray(archives[1][key], dtype=float))
        assert json.loads(archives[0]['feature_pipeline'])['seed_policy'] == 'seedsequence-v2'
