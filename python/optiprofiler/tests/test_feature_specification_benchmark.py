"""
The structured ``feature`` entry through the public ``benchmark``: one route
per call, stage options inside the entries, the run count experiment-wide,
provenance and refined configuration for both routes, and loading.
"""

import json
import pickle

import h5py
import matplotlib
import numpy as np
import pytest

matplotlib.use('Agg')

from optiprofiler import benchmark, profiles
from optiprofiler.loader import load_results_from_h5

NOISY_EFFECTIVE = {'noise_mode': 'random', 'distribution': 'gaussian', 'noise_map': 'chebyshev',
                   'noise_level': 0.01, 'noise_type': 'mixed'}
TRUNCATED_EFFECTIVE = {'perturbed_trailing_digits': False, 'significant_digits': 6}
SPEC = [{'name': 'noisy', 'options': {'noise_level': 1e-2}}, 'truncated']
ARRAYS = ('fun_histories', 'maxcv_histories', 'fun_outs', 'maxcv_outs', 'fun_inits', 'maxcv_inits', 'n_evals')


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


def no_solver(*args, **kwargs):
    pytest.fail('loading saved results must not execute solvers')


def common_kwargs(savepath, benchmark_id='composition'):
    return dict(plibs=['s2mpj'], ptype='u', mindim=2, maxdim=2, max_eval_factor=10, benchmark_id=benchmark_id,
                savepath=str(savepath), n_jobs=1, silent=True, draw_hist_plots='none', problem_names=['ROSENBR'])


def archive_of(directory):
    archives = list(directory.rglob('data_for_loading.h5'))
    assert len(archives) == 1
    return archives[0]


def refined_options_of(directory):
    files = list(directory.rglob('options_refined.pkl'))
    assert len(files) == 1
    with open(files[0], 'rb') as stream:
        return pickle.load(stream)


def assert_same_numbers(first, second):
    for key in ARRAYS:
        np.testing.assert_array_equal(first[key], second[key], err_msg=key)


class TestStructuredRoute:

    def test_structured_feature_runs_and_records_route_and_specification(self, tmp_path):
        report_path = tmp_path / 'report.json'
        scores, _, _ = benchmark([solver_stay, solver_step], feature=SPEC, report_path=str(report_path),
                                 **common_kwargs(tmp_path))
        assert scores.shape == (2,)
        results = load_results_from_h5(str(archive_of(tmp_path / 'composition')))
        assert results[0]['feature_stamp'] == 'noisy_0.01_mixed_gaussian__truncated_6'
        payload = json.loads(results[0]['feature_pipeline'])
        assert payload['schema'] == 'feature_pipeline-v3'
        pipeline = payload['feature']
        assert pipeline['route'] == 'feature'
        assert pipeline['declared_name'] == 'noisy+truncated'
        assert pipeline['declared'] == [{'name': 'noisy', 'options': {'noise_level': 0.01}},
                                        {'name': 'truncated', 'options': {}}]
        assert payload['experiment'] == {'role': 'primary', 'n_runs': 5, 'origin': 'stage_hints',
                                         'run_policy': 'legacy-hints-v1', 'execution_strategy': 'composed-views'}
        assert pipeline['stages'][0]['options'] == NOISY_EFFECTIVE
        assert all('n_runs' not in stage['options'] for stage in pipeline['stages'])
        refined = refined_options_of(tmp_path / 'composition')
        assert refined['schema'] == 'options_refined-v2'
        assert refined['feature_route'] == 'feature'
        assert refined['feature_name'] == 'noisy+truncated'
        assert refined['n_runs'] == 5
        assert refined['feature_specification'] == [{'name': 'noisy', 'options': NOISY_EFFECTIVE},
                                                    {'name': 'truncated', 'options': TRUNCATED_EFFECTIVE}]
        # The structured route has no flat stage keys to record.
        assert 'noise_level' not in refined and 'feature' not in refined
        with open(report_path, encoding='utf-8') as stream:
            report = json.load(stream)
        feature = report['configuration']['effective']['feature']
        assert feature['name'] == 'noisy+truncated' and feature['route'] == 'feature'
        assert 'common_options' not in feature and 'options' not in feature
        assert [entry['name'] for entry in feature['declared']] == ['noisy', 'truncated']
        assert all('n_runs' not in stage['options'] for stage in feature['stages'])
        assert report['configuration']['request']['feature'] == [{'name': 'noisy', 'options': {'noise_level': 0.01}},
                                                                 'truncated']

    def test_shorthand_run_records_the_same_refined_specification(self, tmp_path):
        benchmark([solver_stay, solver_step], feature_name='noisy+truncated', noise_level=1e-2,
                  **common_kwargs(tmp_path, 'shorthand'))
        refined = refined_options_of(tmp_path / 'shorthand')
        assert refined['feature_route'] == 'feature_name'
        assert refined['feature_name'] == 'noisy+truncated'
        assert refined['n_runs'] == 5
        # No flat stage-option projection in the refined configuration.
        assert 'noise_level' not in refined
        assert refined['feature_specification'] == [{'name': 'noisy', 'options': NOISY_EFFECTIVE},
                                                    {'name': 'truncated', 'options': TRUNCATED_EFFECTIVE}]
        payload = json.loads(load_results_from_h5(str(archive_of(tmp_path / 'shorthand')))[0]['feature_pipeline'])
        assert payload['feature']['route'] == 'feature_name'
        # The declaration keeps the supplied broadcast options per owning token.
        assert payload['feature']['declared'] == [{'name': 'noisy', 'options': {'noise_level': 0.01}},
                                                  {'name': 'truncated', 'options': {}}]
        assert payload['experiment']['n_runs'] == 5 and payload['experiment']['origin'] == 'stage_hints'

    def test_single_stage_specification_records_the_legacy_pipeline(self, tmp_path):
        benchmark([solver_stay, solver_step], feature={'name': 'noisy', 'options': {'noise_level': 1e-2}}, n_runs=2,
                  **common_kwargs(tmp_path, 'single'))
        results = load_results_from_h5(str(archive_of(tmp_path / 'single')))
        assert results[0]['feature_stamp'] == 'noisy_0.01_mixed_gaussian'
        # Histories are (problems, solvers, runs, evaluations).
        assert results[0]['fun_histories'].shape[2] == 2
        payload = json.loads(results[0]['feature_pipeline'])
        assert payload['feature']['route'] == 'feature' and payload['feature']['seed_policy'] == 'legacy-run-seed'
        assert payload['experiment']['n_runs'] == 2 and payload['experiment']['origin'] == 'explicit'
        assert payload['feature']['stages'][0]['options'] == NOISY_EFFECTIVE
        refined = refined_options_of(tmp_path / 'single')
        assert refined['feature_specification'] == [{'name': 'noisy', 'options': NOISY_EFFECTIVE}]
        assert refined['n_runs'] == 2 and refined['feature_name'] == 'noisy'

    def test_equal_settings_give_identical_numbers_stamps_and_scores(self, tmp_path):
        structured = benchmark([solver_stay, solver_step], feature=SPEC, **common_kwargs(tmp_path, 'a'))
        shorthand = benchmark([solver_stay, solver_step], feature_name='noisy+truncated', noise_level=1e-2,
                              **common_kwargs(tmp_path, 'b'))
        np.testing.assert_array_equal(structured[0], shorthand[0])
        first = load_results_from_h5(str(archive_of(tmp_path / 'a')))[0]
        second = load_results_from_h5(str(archive_of(tmp_path / 'b')))[0]
        assert_same_numbers(first, second)
        assert first['feature_stamp'] == second['feature_stamp']
        pipelines = [json.loads(result['feature_pipeline']) for result in (first, second)]
        assert [pipeline['feature'].pop('route') for pipeline in pipelines] == ['feature', 'feature_name']
        # Everything else, the declared specification and the plan included, is identical.
        assert pipelines[0] == pipelines[1]

    def test_refined_specification_replays_the_experiment(self, tmp_path):
        original = benchmark([solver_stay, solver_step], feature=SPEC, n_runs=2, **common_kwargs(tmp_path, 'original'))
        refined = refined_options_of(tmp_path / 'original')
        replay = benchmark([solver_stay, solver_step], feature=refined['feature_specification'],
                           n_runs=refined['n_runs'], **common_kwargs(tmp_path, 'replay'))
        np.testing.assert_array_equal(original[0], replay[0])
        assert_same_numbers(load_results_from_h5(str(archive_of(tmp_path / 'original')))[0],
                            load_results_from_h5(str(archive_of(tmp_path / 'replay')))[0])

    def test_randomized_solver_rule_applies_to_the_structured_route(self, tmp_path):
        benchmark([solver_stay, solver_step], feature=['truncated', 'quantized'], solver_isrand=[True, False],
                  **common_kwargs(tmp_path, 'isrand'))
        results = load_results_from_h5(str(archive_of(tmp_path / 'isrand')))
        assert results[0]['fun_histories'].shape[2] == 5
        assert json.loads(results[0]['feature_pipeline'])['experiment']['origin'] == 'randomized_solvers'


class TestErrorsBeforeOutput:

    @pytest.mark.parametrize('kwargs, message', [
        (dict(feature_name='noisy', feature=['noisy']), 'cannot both be given'),
        (dict(feature_name='plain', feature=['noisy']), 'cannot both be given'),
        (dict(feature=None), 'cannot be None'),
        (dict(feature=[]), 'at least one stage entry'),
        (dict(feature=['noisy'], noise_level=1e-2), r"Unexpected keyword\(s\): \['noise_level'\]"),
        (dict(feature=[{'name': 'noisy', 'options': {'n_runs': 3}}]), 'experiment-wide'),
        (dict(feature='noisy+truncated'), 'belongs to `feature_name`'),
        (dict(feature=[{'name': 'perturbed_x0', 'options': {'noise_level': 1e-2}}]),
         r"entry 1 of the feature specification \(stage 'perturbed_x0'\)"),
        (dict(feature=[{'name': 'plain', 'options': {'noise_level': 1.0}}, 'noisy']),
         r"entry 1 of the feature specification \(stage 'plain'\)"),
        (dict(feature=[{'name': 'noisy+truncated'}]), 'stage names are atomic'),
        (dict(feature=[{'name': 'noisy', 'level': 1e-2}]), 'unknown key'),
    ])
    def test_invalid_input_is_rejected_before_any_output(self, tmp_path, kwargs, message):
        with pytest.raises(ValueError, match=message):
            benchmark([solver_stay, solver_step], **kwargs, **common_kwargs(tmp_path, 'rejected'))
        assert not (tmp_path / 'rejected').exists()

    @pytest.mark.parametrize('route', [dict(feature=['truncated', 'quantized']), dict(feature_name='truncated+quantized')])
    def test_explicit_none_run_count_is_rejected_before_output_even_with_randomized_solvers(self, tmp_path, route):
        # A present ``n_runs=None`` bypasses the randomized-solver default and
        # must then fail validation, never run a silently defaulted experiment.
        with pytest.raises(TypeError, match='must be an integer'):
            benchmark([solver_stay, solver_step], n_runs=None, solver_isrand=[True, False], **route,
                      **common_kwargs(tmp_path, 'rejected'))
        assert not (tmp_path / 'rejected').exists()

    def test_requested_report_records_the_failure(self, tmp_path):
        report_path = tmp_path / 'failure.json'
        with pytest.raises(ValueError, match='cannot both be given'):
            benchmark([solver_stay, solver_step], feature_name='noisy', feature=['noisy'], report_path=str(report_path),
                      **common_kwargs(tmp_path, 'rejected'))
        assert not (tmp_path / 'rejected').exists()
        with open(report_path, encoding='utf-8') as stream:
            report = json.load(stream)
        assert report['status'] == 'failed'
        diagnostics = json.dumps(report['diagnostics'])
        assert 'benchmark_exception' in diagnostics and 'ValueError' in diagnostics
        assert report['coverage']['completed'] == 0


class TestLoad:

    def test_structured_archive_loads_without_a_new_feature_argument(self, tmp_path, monkeypatch):
        scores, _, _ = benchmark([solver_stay, solver_step], feature=SPEC, **common_kwargs(tmp_path, 'reload'))
        monkeypatch.chdir(tmp_path)
        report_path = tmp_path / 'load-report.json'
        loaded, _, _ = benchmark([no_solver, no_solver], load='latest', savepath=str(tmp_path), benchmark_id='reload',
                                 score_only=True, silent=True, draw_hist_plots='none', report_path=str(report_path))
        np.testing.assert_array_equal(loaded, scores)
        with open(report_path, encoding='utf-8') as stream:
            report = json.load(stream)
        retained = report['configuration']['retained_result_metadata'][0]['feature_pipeline']
        assert retained['schema'] == 'feature_pipeline-v3' and retained['feature']['route'] == 'feature'
        assert retained['experiment']['n_runs'] == 5 and retained['experiment']['role'] == 'primary'
        assert [stage['identity'] for stage in retained['feature']['stages']] == ['noisy#0', 'truncated#0']
        assert all('n_runs' not in stage['options'] for stage in retained['feature']['stages'])
        effective = report['configuration']['effective']['feature']
        assert effective['scope'] == 'current_load_context_not_original_execution_feature'
        assert effective['name'] == 'plain' and effective['route'] == 'feature_name'

    def test_feature_name_still_labels_a_load(self, tmp_path, monkeypatch):
        scores, _, _ = benchmark([solver_stay, solver_step], feature=SPEC, **common_kwargs(tmp_path, 'reload'))
        monkeypatch.chdir(tmp_path)
        report_path = tmp_path / 'load-report.json'
        loaded, _, _ = benchmark([no_solver, no_solver], load='latest', feature_name='noisy', savepath=str(tmp_path),
                                 benchmark_id='reload', score_only=True, silent=True, draw_hist_plots='none',
                                 report_path=str(report_path))
        np.testing.assert_array_equal(loaded, scores)
        with open(report_path, encoding='utf-8') as stream:
            report = json.load(stream)
        assert report['configuration']['effective']['feature']['name'] == 'noisy'
        assert report['configuration']['retained_result_metadata'][0]['feature_pipeline']['feature']['route'] == 'feature'

    def test_feature_is_rejected_in_load_mode_before_loading(self, tmp_path, monkeypatch):
        benchmark([solver_stay, solver_step], feature=SPEC, **common_kwargs(tmp_path, 'reload'))
        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(profiles, 'load_results',
                            lambda *args, **kwargs: pytest.fail('a rejected feature must not start loading'))
        before = sorted(path.name for path in tmp_path.iterdir())
        with pytest.raises(ValueError, match='cannot be used with `load`'):
            benchmark([no_solver, no_solver], load='latest', feature=['noisy'], savepath=str(tmp_path),
                      benchmark_id='reload', silent=True, draw_hist_plots='none')
        assert sorted(path.name for path in tmp_path.iterdir()) == before

    def test_v1_archive_pipeline_is_retained_verbatim(self, tmp_path, monkeypatch):
        scores, _, _ = benchmark([solver_stay, solver_step], feature_name='quantized+noisy',
                                 **common_kwargs(tmp_path, 'reload'))
        v1 = {'schema': 'feature_pipeline-v1', 'declared_name': 'quantized+noisy', 'effective_name': 'quantized+noisy',
              'seed_policy': 'seedsequence-v2', 'feature_stamp': 'quantized_0.001_ground_truth__noisy_0.001_mixed_gaussian',
              'full_feature_stamp': 'quantized_0.001_ground_truth__noisy_0.001_mixed_gaussian',
              'stages': [{'position': 0, 'name': 'quantized', 'code': 9, 'occurrence': 0, 'identity': 'quantized#0',
                          'options': {'n_runs': 1, 'mesh_size': 0.001, 'mesh_type': 'absolute', 'ground_truth': True}},
                         {'position': 1, 'name': 'noisy', 'code': 2, 'occurrence': 0, 'identity': 'noisy#0',
                          'options': {'noise_mode': 'random', 'distribution': 'gaussian', 'noise_map': 'chebyshev',
                                      'n_runs': 5, 'noise_level': 0.001, 'noise_type': 'mixed'}}]}
        with h5py.File(archive_of(tmp_path / 'reload'), 'r+') as archive:
            for group in archive.values():
                del group['feature_pipeline']
                group.create_dataset('feature_pipeline',
                                     data=np.array(json.dumps(v1, sort_keys=True), dtype=h5py.special_dtype(vlen=str)))
        monkeypatch.chdir(tmp_path)
        report_path = tmp_path / 'load-report.json'
        loaded, _, _ = benchmark([no_solver, no_solver], load='latest', savepath=str(tmp_path), benchmark_id='reload',
                                 score_only=True, silent=True, draw_hist_plots='none', report_path=str(report_path))
        np.testing.assert_array_equal(loaded, scores)
        with open(report_path, encoding='utf-8') as stream:
            report = json.load(stream)
        assert report['configuration']['retained_result_metadata'][0]['feature_pipeline'] == v1
