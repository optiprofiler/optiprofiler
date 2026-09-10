"""
Provenance of features must never execute user code.

A callable option is described by its function or class name only. Objects
with hostile or stateful representation and attribute hooks are legal
options: their ``__call__`` is what the feature uses, and the metadata
encoder must not touch anything else.
"""

import json

import matplotlib
import numpy as np
import pytest

matplotlib.use('Agg')

from optiprofiler import benchmark
from optiprofiler.composition import describe_pipeline
from optiprofiler.loader import load_results_from_h5
from optiprofiler.opclasses import Feature, Problem


class HostileNoiseMap:
    """Valid noise map whose representation and attribute hooks raise."""

    def __call__(self, x):
        return 0.5

    def __repr__(self):
        raise RuntimeError('metadata_called_repr')

    def __str__(self):
        raise RuntimeError('metadata_called_str')

    def __getattr__(self, name):
        raise RuntimeError(f'metadata_called_getattr:{name}')


class StatefulNoiseMap:
    """Valid noise map whose representation records that it was invoked."""

    def __init__(self):
        self.representations = 0

    def __call__(self, x):
        return 0.5

    def __repr__(self):
        self.representations += 1
        return 'StatefulNoiseMap'


class HostileModFun:

    def __call__(self, x, rng, problem):
        return problem.fun(x) + 1.0

    def __repr__(self):
        raise RuntimeError('metadata_called_repr')

    def __getattr__(self, name):
        raise RuntimeError(f'metadata_called_getattr:{name}')


def stay(fun, x0):
    fun(x0)
    return x0


def step(fun, x0):
    fun(x0 + 0.1)
    return x0 + 0.1


def direct_kwargs(tmp_path, **extra):
    options = dict(problem=Problem(lambda x: float(x @ x), np.array([1.0, -1.0]), name='QUAD'), score_only=True,
                   draw_hist_plots='none', silent=True, savepath=str(tmp_path), n_jobs=1)
    options.update(extra)
    return options


def library_kwargs(tmp_path, **extra):
    options = dict(plibs=['s2mpj'], ptype='u', mindim=2, maxdim=2, problem_names=['ROSENBR'], max_eval_factor=5,
                   n_jobs=1, silent=True, draw_hist_plots='none', savepath=str(tmp_path), benchmark_id='meta')
    options.update(extra)
    return options


def callback_record(record):
    assert record['kind'] == 'callback' and record['instance_state'] == 'not_recorded'
    return record['name']


class TestBenchmarkNeverExecutesOptionRepresentations:

    def test_direct_problem_single_feature_with_report(self, tmp_path):
        report_path = tmp_path / 'report.json'
        scores, _, _ = benchmark([stay, step], feature_name='noisy', noise_mode='deterministic',
                                 noise_map=HostileNoiseMap(), report_path=str(report_path), **direct_kwargs(tmp_path))
        assert scores.shape == (2,)
        with open(report_path, encoding='utf-8') as stream:
            feature = json.load(stream)['configuration']['effective']['feature']
        assert callback_record(feature['options']['noise_map']) == 'HostileNoiseMap'
        assert callback_record(feature['stages'][0]['options']['noise_map']) == 'HostileNoiseMap'

    def test_library_single_feature_score_only_without_report(self, tmp_path):
        scores, _, _ = benchmark([stay, step], feature_name='noisy', noise_mode='deterministic',
                                 noise_map=HostileNoiseMap(), score_only=True, **library_kwargs(tmp_path))
        assert scores.shape == (2,)

    def test_composite_with_report_and_archive(self, tmp_path):
        report_path = tmp_path / 'report.json'
        scores, _, _ = benchmark([stay, step], feature_name='noisy+truncated', noise_mode='deterministic',
                                 noise_map=HostileNoiseMap(), report_path=str(report_path), **library_kwargs(tmp_path))
        assert scores.shape == (2,)
        archives = list((tmp_path / 'meta').rglob('data_for_loading.h5'))
        assert len(archives) == 1
        pipeline = json.loads(load_results_from_h5(str(archives[0]))[0]['feature_pipeline'])
        assert callback_record(pipeline['stages'][0]['options']['noise_map']) == 'HostileNoiseMap'
        assert pipeline['stages'][1]['options']['significant_digits'] == 6
        with open(report_path, encoding='utf-8') as stream:
            feature = json.load(stream)['configuration']['effective']['feature']
        assert callback_record(feature['stages'][0]['options']['noise_map']) == 'HostileNoiseMap'

    def test_custom_stage_with_callable_instance(self, tmp_path):
        report_path = tmp_path / 'report.json'
        scores, _, _ = benchmark([stay, step], feature_name='custom+noisy', mod_fun=HostileModFun(),
                                 noise_level=0.0, report_path=str(report_path), **direct_kwargs(tmp_path))
        assert scores.shape == (2,)
        with open(report_path, encoding='utf-8') as stream:
            feature = json.load(stream)['configuration']['effective']['feature']
        assert callback_record(feature['stages'][0]['options']['mod_fun']) == 'HostileModFun'

    def test_stateful_representation_is_never_invoked(self, tmp_path):
        noise_map = StatefulNoiseMap()
        benchmark([stay, step], feature_name='plain+noisy', noise_mode='deterministic', noise_map=noise_map,
                  report_path=str(tmp_path / 'report.json'), **direct_kwargs(tmp_path))
        benchmark([stay, step], feature_name='noisy+truncated', noise_mode='deterministic', noise_map=noise_map,
                  report_path=str(tmp_path / 'report2.json'), **direct_kwargs(tmp_path))
        assert noise_map.representations == 0


class TestDescribePipeline:

    def test_callables_are_described_by_name_only(self):
        def plain_function(x):
            return 1.0

        feature = Feature('noisy+custom', noise_mode='deterministic', noise_map=HostileNoiseMap(),
                          mod_fun=HostileModFun(), mod_x0=plain_function)
        pipeline = describe_pipeline(feature)
        options = pipeline['stages'][0]['options']
        assert callback_record(options['noise_map']) == 'HostileNoiseMap'
        custom = pipeline['stages'][1]['options']
        assert callback_record(custom['mod_fun']) == 'HostileModFun'
        record = custom['mod_x0']
        assert record['kind'] == 'callback' and record['name'].endswith('plain_function')
        assert 'instance_state' not in record

    def test_non_finite_values_become_reason_records(self):
        pipeline = describe_pipeline(Feature('noisy+truncated', noise_level=float('nan')))
        assert pipeline['stages'][0]['options']['noise_level'] == {'value': None, 'reason': 'nan'}

    def test_report_and_pipeline_share_one_encoder(self):
        import optiprofiler.eval_report as eval_report
        import optiprofiler.metadata as metadata
        assert eval_report._safe is metadata.safe_metadata
        assert eval_report._callback is metadata.describe_callback
