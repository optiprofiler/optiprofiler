"""
Report contract ``optiprofiler.eval_report/2``: explicit schema selection,
the immutable version 1 resource, the canonical feature block and the
experiment plans, and the separation between the unchanged archive bytes and
the sanitized metadata retained by a load report.
"""

import hashlib
import json

import h5py
import matplotlib
import numpy as np
import pytest

matplotlib.use('Agg')

from optiprofiler import benchmark
from optiprofiler import eval_report as module
from optiprofiler.opclasses import Problem
from optiprofiler.tests.eval_report_contract import assert_valid

# eval_report.schema.json is version 1 of the main report: the contract of the
# reports already written and of the MATLAB producer. It never changes.
V1_SHA256 = '9b8f76fd9b2d3eff323304b758ff61c5c456d35b1d81fcf05d3ab9d64453a8a3'
PLOT_DATA_SHA256 = '76f83e1b7922a1479d8770888bd417c86fe5f837a450789891618960582d323c'


def stay(fun, x0):
    fun(x0)
    return np.asarray(x0, dtype=float)


def zero(fun, x0):
    fun(x0)
    return np.zeros_like(np.asarray(x0, dtype=float))


def forbidden(fun, x0):
    raise SystemExit('loading saved results must not execute solvers')


def read(path):
    return json.loads(path.read_text(encoding='utf-8'),
                      parse_constant=lambda token: (_ for _ in ()).throw(ValueError(token)))


def quad():
    return Problem(lambda x: float(x @ x), [1.0, 2.0], name='QUAD')


class TestSchemaSelection:

    def test_version_one_resources_are_immutable(self):
        text = module.schema_text('eval_report', 1)
        assert hashlib.sha256(text.encode('utf-8')).hexdigest() == V1_SHA256
        assert module.schema_text('eval_report.schema.json') == text
        assert module.schema_text('optiprofiler.eval_report/1') == text
        assert module.schema_text('eval_report-v1') == text
        v1 = module.load_schema('eval_report', 1)
        assert v1['$id'] == 'urn:optiprofiler:eval_report:1'
        assert v1['properties']['schema']['const'] == 'optiprofiler.eval_report/1'
        assert hashlib.sha256(module.schema_text('plot_data').encode('utf-8')).hexdigest() == PLOT_DATA_SHA256

    def test_current_main_schema_is_version_two(self):
        v2 = module.load_schema('eval_report')
        assert v2['$id'] == 'urn:optiprofiler:eval_report:2'
        assert v2['properties']['schema']['const'] == 'optiprofiler.eval_report/2'
        assert module.load_schema('eval_report-v2') == v2 == module.load_schema('eval_report-v2.schema.json')
        assert module.load_schema('optiprofiler.eval_report/2') == v2 == module.load_schema('eval_report', 2)
        assert module.schema_resource('eval_report') == 'eval_report-v2.schema.json'
        assert module.schema_resource('eval_report', 1) == 'eval_report.schema.json'
        assert module.schema_identifier('eval_report') == 'optiprofiler.eval_report/2'
        assert module.schema_identifier('plot_data') == 'optiprofiler.plot_data/1'
        # The numeric companion keeps its contract; v2 references the same identifier.
        assert module.load_schema('plot_data')['$id'] == 'urn:optiprofiler:plot_data:1'
        assert v2['$defs']['plotDataReference']['properties']['schema']['const'] == 'optiprofiler.plot_data/1'
        effective = v2['$defs']['configuration']['properties']['effective']
        assert effective['required'] == ['problem_options', 'profile_options', 'feature', 'experiment']
        assert v2['$defs']['experimentPlan']['required'] == ['role', 'n_runs', 'origin', 'run_policy', 'execution_strategy',
                                                           'runtime_policy']
        assert v2['$defs']['experimentPlans']['additionalProperties'] is False

    def test_dispatch_by_document_identity_rejects_unknown_versions(self):
        assert module.schema_for_document({'schema': 'optiprofiler.eval_report/1'}) == ('eval_report', 1)
        assert module.schema_for_document({'schema': 'optiprofiler.eval_report/2'}) == ('eval_report', 2)
        assert module.schema_for_document({'schema': 'optiprofiler.plot_data/1'}) == ('plot_data', 1)
        assert module.load_schema_for({'schema': 'optiprofiler.eval_report/1'})['$id'] == 'urn:optiprofiler:eval_report:1'
        for identifier in ('optiprofiler.eval_report/3', 'optiprofiler.eval_report', 'eval_report/2',
                           'optiprofiler.agent_report/1', 'optiprofiler.plot_data/2', '', None, 2):
            with pytest.raises(ValueError):
                module.schema_for_document({'schema': identifier})
        with pytest.raises(ValueError):
            module.schema_for_document({})
        with pytest.raises(ValueError, match='Unknown EvalReport schema'):
            module.schema_text('eval_report', 3)
        with pytest.raises(ValueError, match='Unknown EvalReport schema'):
            module.schema_text('eval_report.schema.json', 2)
        with pytest.raises(TypeError):
            module.schema_text(None)


class TestProducedReports:

    def test_benchmark_report_states_feature_and_plans_once(self, tmp_path):
        target = tmp_path / 'report.json'
        benchmark([stay, zero], report_path=target, plibs=['s2mpj'], ptype='u', mindim=2, maxdim=2,
                  problem_names=['ROSENBR'], max_eval_factor=10, benchmark_id='plans', n_jobs=1, score_only=True,
                  draw_hist_plots='none', silent=True, savepath=str(tmp_path), n_runs=2, run_plain=True,
                  feature=[{'name': 'noisy', 'options': {'noise_level': 0.5}}, 'truncated'])
        report = read(target)
        assert report['schema'] == 'optiprofiler.eval_report/2' and report['status'] == 'completed'
        assert_valid(report, 'eval_report-v2.schema.json')
        effective = report['configuration']['effective']
        feature = effective['feature']
        assert (feature['scope'], feature['route'], feature['name']) == ('current_execution_feature', 'feature', 'noisy+truncated')
        assert feature['declaration_route'] == 'feature' and feature['feature_stamp_origin'] == 'generated'
        assert feature['feature_stamp'] == feature['full_feature_stamp'] == 'noisy_0.5_mixed_gaussian__truncated_6'
        assert [stage['identity'] for stage in feature['stages']] == ['noisy#0', 'truncated#0']
        assert feature['declared'] == [{'name': 'noisy', 'options': {'noise_level': 0.5}}, {'name': 'truncated', 'options': {}}]
        assert feature['stages'][0]['options']['noise_level'] == 0.5 and feature['stages'][1]['options']['significant_digits'] == 6
        assert 'options' not in feature and 'n_runs' not in json.dumps(feature)
        assert effective['experiment'] == {
            'primary': {'role': 'primary', 'n_runs': 2, 'origin': 'explicit', 'run_policy': 'legacy-hints-v1',
                        'execution_strategy': 'composed-views', 'runtime_policy': 'python-featured-problem-v1'},
            'plain_reference': {'role': 'plain_reference', 'n_runs': 1, 'origin': 'reference_policy',
                                'run_policy': 'legacy-hints-v1', 'execution_strategy': 'identity',
                                'runtime_policy': 'python-featured-problem-v1'}}
        roles = {(problem['role'], len(problem['runs'])) for problem in report['problems']}
        assert roles == {('primary', 4), ('plain_reference', 2)}
        # Runtime facts are stated once per role; Python emits no per-run receipt.
        for problem in report['problems']:
            for run in problem['runs']:
                assert run['execution']['kind'] == 'actual' and 'runtime' not in run

    def test_default_count_and_identity_feature_in_the_experiment_block(self, tmp_path):
        target = tmp_path / 'report.json'
        benchmark([stay, zero], report_path=target, problem=quad(), score_only=True, draw_hist_plots='none',
                  silent=True, savepath=str(tmp_path))
        report = read(target)
        feature = report['configuration']['effective']['feature']
        assert (feature['name'], feature['effective_name'], feature['stages']) == ('plain', 'plain', [])
        assert feature['declared'] == [{'name': 'plain', 'options': {}}] and feature['declaration_route'] == 'feature_name'
        assert feature['route'] is None  # neither keyword was given: the default plain feature
        assert report['configuration']['effective']['experiment'] == {
            'primary': {'role': 'primary', 'n_runs': 1, 'origin': 'stage_hints', 'run_policy': 'legacy-hints-v1',
                        'execution_strategy': 'identity', 'runtime_policy': 'python-featured-problem-v1'}}

    def test_reused_feature_object_keeps_its_declaration_route_and_records_the_input(self, tmp_path):
        from optiprofiler import Feature
        declared_once = Feature('noisy+truncated', noise_level=0.5)
        target = tmp_path / 'report.json'
        benchmark([stay, zero], report_path=target, problem=quad(), score_only=True, draw_hist_plots='none',
                  silent=True, savepath=str(tmp_path), feature=declared_once)
        report = read(target)
        assert_valid(report, 'eval_report-v2.schema.json')
        feature = report['configuration']['effective']['feature']
        assert (feature['declaration_route'], feature['route'], feature['declared_name']) == ('feature_name', 'feature', 'noisy+truncated')
        assert feature['declared'] == [{'name': 'noisy', 'options': {'noise_level': 0.5}}, {'name': 'truncated', 'options': {}}]
        target = tmp_path / 'shorthand.json'
        benchmark([stay, zero], report_path=target, problem=quad(), score_only=True, draw_hist_plots='none',
                  silent=True, savepath=str(tmp_path), feature_name='noisy+truncated', noise_level=0.5)
        shorthand = read(target)['configuration']['effective']['feature']
        assert (shorthand['declaration_route'], shorthand['route']) == ('feature_name', 'feature_name')
        for key in ('declared', 'declared_name', 'effective_name', 'stages', 'seed_policy'):
            assert shorthand[key] == feature[key]

    def test_explicit_stamp_and_copied_runs(self, tmp_path):
        target = tmp_path / 'report.json'
        benchmark([stay, zero], report_path=target, problem=quad(), score_only=True, draw_hist_plots='none',
                  silent=True, savepath=str(tmp_path), feature_name='truncated', n_runs=3, feature_stamp='mine')
        report = read(target)
        assert_valid(report, 'eval_report-v2.schema.json')
        feature = report['configuration']['effective']['feature']
        assert (feature['feature_stamp'], feature['feature_stamp_origin']) == ('mine', 'explicit')
        assert report['configuration']['effective']['experiment']['primary']['execution_strategy'] == 'legacy-single'
        [problem] = report['problems']
        kinds = [(run['run_index'], run['execution']['kind'], 'runtime' in run) for run in problem['runs']]
        # Deterministic feature and solvers: one actual run per solver, the rest copied; no per-run receipts.
        assert kinds == [(1, 'actual', False), (2, 'repeated', False), (3, 'repeated', False)] * 2

    def test_long_chain_text_is_omitted_with_length_and_reason_not_clipped(self, tmp_path):
        target = tmp_path / 'report.json'
        chain = '+'.join(['noisy'] * 50)
        benchmark([stay, zero], report_path=target, problem=quad(), score_only=True, draw_hist_plots='none',
                  silent=True, savepath=str(tmp_path), feature_name=chain, n_runs=1)
        report = read(target)
        assert_valid(report, 'eval_report-v2.schema.json')
        feature = report['configuration']['effective']['feature']
        assert len(feature['stages']) == 50 and feature['stages'][-1]['identity'] == 'noisy#49'
        full = '__'.join(['noisy_0.001_mixed_gaussian'] * 50)
        assert len(full) > 256 and len(chain) > 256
        for field, native in (('full_feature_stamp', full), ('effective_name', chain), ('declared_name', chain), ('name', chain)):
            assert feature[field] is None
            assert feature[field + '_bytes'] == len(native.encode('utf-8'))
            assert feature[field + '_reason'] == 'omitted_from_bounded_metadata_projection'
        # The bounded stamp is short and present; nothing is clipped silently.
        assert len(feature['feature_stamp']) <= 64 and 'feature_stamp_reason' not in feature
        assert not any(isinstance(value, str) and len(value) == 256 for value in feature.values())

    def test_failed_report_before_configuration_selects_the_current_schema(self, tmp_path):
        target = tmp_path / 'failed.json'
        with pytest.raises(ValueError, match='at least one stage entry'):
            benchmark([stay, zero], report_path=target, problem=quad(), feature=[], score_only=True,
                      draw_hist_plots='none', silent=True, savepath=str(tmp_path))
        report = read(target)
        assert report['schema'] == 'optiprofiler.eval_report/2' and report['status'] == 'failed'
        assert 'effective' not in report['configuration']
        assert_valid(report, 'eval_report-v2.schema.json')

    def test_load_report_separates_archive_bytes_from_sanitized_metadata(self, tmp_path, monkeypatch):
        kwargs = dict(plibs=['s2mpj'], ptype='u', mindim=2, maxdim=2, problem_names=['ROSENBR'], max_eval_factor=10,
                      benchmark_id='v2load', savepath=str(tmp_path), n_jobs=1, silent=True, draw_hist_plots='none')
        benchmark([stay, zero], feature_name='noisy+truncated', n_runs=2, **kwargs)
        archives = list(tmp_path.rglob('data_for_loading.h5'))
        assert len(archives) == 1
        archive_bytes = archives[0].read_bytes()
        monkeypatch.chdir(tmp_path)
        target = tmp_path / 'load-report.json'
        benchmark([forbidden, forbidden], load='latest', benchmark_id='v2load', savepath=str(tmp_path), score_only=True,
                  silent=True, draw_hist_plots='none', report_path=target)
        report = read(target)
        assert report['schema'] == 'optiprofiler.eval_report/2' and report['operation'] == 'load'
        assert_valid(report, 'eval_report-v2.schema.json')
        configuration = report['configuration']
        assert configuration['effective']['experiment'] == {}
        assert all('runtime' not in run for problem in report['problems'] for run in problem['runs'])
        # The raw observations are recorded before the display options are
        # resolved and again afterwards: the prepared plots exist and no stale
        # 'unavailable' diagnostic from the first pass survives.
        assert all(problem['plot_refs'] for problem in report['problems'])
        assert [d for d in report['diagnostics'] if d['code'] == 'history_plot_data_unavailable'] == []
        assert report['status'] == 'completed'
        assert configuration['effective']['feature']['scope'] == 'current_load_context_not_original_execution_feature'
        assert configuration['retained_result_metadata_encoding'].startswith('sanitized_metadata_copy;not_byte_identical')
        # The archive itself is untouched and identified by its digest; the
        # retained metadata is a sanitized copy of what the archive recorded.
        assert archives[0].read_bytes() == archive_bytes
        assert report['source']['sha256'] == hashlib.sha256(archive_bytes).hexdigest()
        assert report['source']['bytes'] == len(archive_bytes)
        retained = configuration['retained_result_metadata']
        assert len(retained) == 1 and retained[0]['scope'].startswith('retained_result_after_load_filtering')
        pipeline = retained[0]['feature_pipeline']
        assert pipeline['schema'] == 'feature_pipeline-v3'
        assert pipeline['experiment'] == {'role': 'primary', 'n_runs': 2, 'origin': 'explicit',
                                          'run_policy': 'legacy-hints-v1', 'execution_strategy': 'composed-views',
                                          'runtime_policy': 'python-featured-problem-v1'}
        with h5py.File(archives[0], 'r') as archive:
            stored = archive['plib_0']['feature_pipeline'][()]
        stored = json.loads(stored.decode('utf-8') if isinstance(stored, bytes) else stored)
        assert pipeline == stored
