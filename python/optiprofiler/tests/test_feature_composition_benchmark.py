"""
Composed features through the public benchmark: stamps, archives, loading,
structured report metadata and notes.
"""

import json
import pickle
import re

import h5py
import matplotlib
import numpy as np
import pytest

matplotlib.use('Agg')

from optiprofiler import benchmark
from optiprofiler.loader import load_results_from_h5
from optiprofiler.opclasses import Feature, FeaturedProblem, Problem
from optiprofiler.profile_utils import get_default_profile_options

NON_PLAIN = ['perturbed_x0', 'noisy', 'truncated', 'permuted', 'linearly_transformed', 'random_nan',
             'unrelaxable_constraints', 'nonquantifiable_constraints', 'quantized']


def solver_stay(fun, x0):
    fun(x0)
    return x0


def solver_step(fun, x0):
    best_x, best_f = np.array(x0, dtype=float), fun(x0)
    for i in range(len(x0)):
        for step in (0.5, -0.5):
            trial = best_x.copy()
            trial[i] += step
            value = fun(trial)
            if value < best_f:
                best_x, best_f = trial, value
    return best_x


def common_kwargs(savepath, benchmark_id='composition'):
    return dict(plibs=['s2mpj'], ptype='u', mindim=2, maxdim=2, max_eval_factor=10, benchmark_id=benchmark_id,
                savepath=str(savepath), n_jobs=1, silent=True, draw_hist_plots='none', problem_names=['ROSENBR'])


def archive_of(directory):
    archives = list(directory.rglob('data_for_loading.h5'))
    assert len(archives) == 1
    return archives[0]


class TestStamps:

    def test_composite_stamp_joins_ordered_child_stamps(self):
        solvers = [solver_stay, solver_step]
        stamp = get_default_profile_options(solvers, Feature('noisy+truncated'), {})['feature_stamp']
        assert stamp == 'noisy_0.001_mixed_gaussian__truncated_6'
        stamp = get_default_profile_options(solvers, Feature('truncated+noisy'), {})['feature_stamp']
        assert stamp == 'truncated_6__noisy_0.001_mixed_gaussian'
        stamp = get_default_profile_options(solvers, Feature('noisy+plain'), {})['feature_stamp']
        assert stamp == 'noisy_0.001_mixed_gaussian'

    def test_long_stamps_are_bounded_and_digest_qualified(self):
        solvers = [solver_stay, solver_step]
        first = get_default_profile_options(solvers, Feature('+'.join(NON_PLAIN * 2)), {})['feature_stamp']
        second = get_default_profile_options(solvers, Feature('+'.join(reversed(NON_PLAIN * 2))), {})['feature_stamp']
        for stamp in (first, second):
            assert len(stamp) <= 64
            assert re.fullmatch(r'[a-zA-Z0-9_.-]+', stamp)
            assert re.search(r'_[0-9a-f]{8}$', stamp)
        assert first != second


class TestPickling:

    def test_composed_feature_survives_pickling(self):
        feature = Feature('noisy+quantized', noise_level=0.01, mesh_size=0.5)
        clone = pickle.loads(pickle.dumps(feature))
        assert clone.name == feature.name and clone.options == feature.options
        assert clone.is_stochastic is True
        problem = Problem(lambda x: float(x @ x), np.array([1.0, -2.0]))
        x = np.array([0.4, 0.6])
        assert FeaturedProblem(problem, clone, 5, 3).fun(x) == FeaturedProblem(problem, feature, 5, 3).fun(x)


class TestBenchmarkIntegration:

    def test_composed_feature_runs_and_records_its_pipeline(self, tmp_path):
        report_path = tmp_path / 'report.json'
        scores, profile_scores, curves = benchmark([solver_stay, solver_step], feature_name='noisy+truncated',
                                                   report_path=str(report_path), **common_kwargs(tmp_path))
        assert scores.shape == (2,)
        # The output folder name carries the feature stamp when the platform path
        # limit allows it; the archive always records the stamp itself.
        results = load_results_from_h5(str(archive_of(tmp_path / 'composition')))
        assert results[0]['feature_stamp'] == 'noisy_0.001_mixed_gaussian__truncated_6'
        pipeline = json.loads(results[0]['feature_pipeline'])
        assert pipeline['schema'] == 'feature_pipeline-v1'
        assert pipeline['declared_name'] == 'noisy+truncated' and pipeline['effective_name'] == 'noisy+truncated'
        assert pipeline['seed_policy'] == 'seedsequence-v2'
        assert [stage['identity'] for stage in pipeline['stages']] == ['noisy#0', 'truncated#0']
        assert pipeline['stages'][0]['options']['distribution'] == 'gaussian'
        assert pipeline['stages'][0]['options']['n_runs'] == 5
        assert pipeline['stages'][1]['options']['significant_digits'] == 6
        assert pipeline['full_feature_stamp'] == 'noisy_0.001_mixed_gaussian__truncated_6'
        with open(report_path, encoding='utf-8') as stream:
            report = json.load(stream)
        feature = report['configuration']['effective']['feature']
        assert feature['name'] == 'noisy+truncated'
        assert feature['declared_name'] == 'noisy+truncated'
        assert feature['seed_policy'] == 'seedsequence-v2'
        assert [stage['identity'] for stage in feature['stages']] == ['noisy#0', 'truncated#0']
        assert feature['options'] == {'n_runs': 5}

    def test_single_feature_archive_records_a_legacy_pipeline(self, tmp_path):
        benchmark([solver_stay, solver_step], feature_name='plain+noisy', **common_kwargs(tmp_path, 'single'))
        results = load_results_from_h5(str(archive_of(tmp_path / 'single')))
        assert results[0]['feature_stamp'] == 'noisy_0.001_mixed_gaussian'
        pipeline = json.loads(results[0]['feature_pipeline'])
        assert pipeline['declared_name'] == 'plain+noisy' and pipeline['effective_name'] == 'noisy'
        assert pipeline['seed_policy'] == 'legacy-run-seed'
        assert [stage['identity'] for stage in pipeline['stages']] == ['noisy#0']

    def test_unknown_token_is_rejected_before_any_output(self, tmp_path):
        with pytest.raises(ValueError, match='Unknown feature: unknown'):
            benchmark([solver_stay, solver_step], feature_name='noisy+unknown', **common_kwargs(tmp_path, 'rejected'))
        with pytest.raises(ValueError, match='empty stage token'):
            benchmark([solver_stay, solver_step], feature_name='noisy++truncated', **common_kwargs(tmp_path, 'rejected'))
        assert not (tmp_path / 'rejected').exists()

    def test_loading_restores_metadata_without_solvers_and_tolerates_old_archives(self, tmp_path, monkeypatch):
        kwargs = common_kwargs(tmp_path, 'reload')
        benchmark([solver_stay, solver_step], feature_name='quantized+noisy', **kwargs)
        monkeypatch.chdir(tmp_path)

        def no_solver(*args, **kwargs):
            pytest.fail('loading saved results must not execute solvers')

        report_path = tmp_path / 'load-report.json'
        scores, profile_scores, curves = benchmark([no_solver, no_solver], load='latest', savepath=str(tmp_path),
                                                   benchmark_id='reload', score_only=True, silent=True,
                                                   draw_hist_plots='none', report_path=str(report_path))
        assert scores.shape == (2,)
        with open(report_path, encoding='utf-8') as stream:
            report = json.load(stream)
        retained = report['configuration']['retained_result_metadata']
        assert retained[0]['feature_stamp'] == 'quantized_0.001_ground_truth__noisy_0.001_mixed_gaussian'
        assert [stage['identity'] for stage in retained[0]['feature_pipeline']['stages']] == ['quantized#0', 'noisy#0']
        # An archive written before compositions existed carries no pipeline entry.
        with h5py.File(archive_of(tmp_path / 'reload'), 'r+') as archive:
            for group in archive.values():
                del group['feature_pipeline']
        report_path = tmp_path / 'load-report-old-archive.json'
        scores, _, _ = benchmark([no_solver, no_solver], load='latest', savepath=str(tmp_path), benchmark_id='reload',
                                 score_only=True, silent=True, draw_hist_plots='none', report_path=str(report_path))
        assert scores.shape == (2,)
        with open(report_path, encoding='utf-8') as stream:
            report = json.load(stream)
        assert report['configuration']['retained_result_metadata'][0]['feature_pipeline'] is None

    def test_quantized_notes_name_each_stage(self, tmp_path):
        benchmark([solver_stay, solver_step], feature_name='quantized+noisy+quantized', ground_truth=False,
                  **common_kwargs(tmp_path, 'notes'))
        texts = [path.read_text(encoding='utf-8') for path in (tmp_path / 'notes').rglob('*.txt')]
        notes = [line for text in texts for line in text.splitlines() if line.startswith('Quantized truth')]
        assert any("(stage 1 'quantized#0')" in line and 'original' in line for line in notes)
        assert any("(stage 3 'quantized#1')" in line and 'ground_truth=false' in line for line in notes)
