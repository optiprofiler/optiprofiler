"""Tests for problem-library discovery and protocol validation."""

import json
import pickle
from pathlib import Path

import pytest

import optiprofiler
from optiprofiler.profile_utils import check_validity_problem_options
from optiprofiler.problem_libraries import (
    PROBLEM_LIBRARY_API_VERSION,
    PROBLEM_LIBRARY_ENTRY_POINT_GROUP,
    ProblemLibraryPlugin,
    ProblemLibraryRef,
    _resolve_problem_library_options,
    list_problem_libraries,
    load_problem_library,
    resolve_problem_library,
)
from optiprofiler import problem_libraries


def _select(options, library_options):
    return ['TOY']


def _load(name, library_options):
    return name


class _FakeDistribution:
    def __init__(self, name):
        self.name = name


class _FakeEntryPoint:
    def __init__(self, name, value, loaded, distribution='test-distribution'):
        self.name = name
        self.value = value
        self.dist = _FakeDistribution(distribution)
        self._loaded = loaded
        self.load_calls = 0

    def load(self):
        self.load_calls += 1
        return self._loaded


def _plugin(name='externaltoy', api_version=PROBLEM_LIBRARY_API_VERSION):
    return ProblemLibraryPlugin(
        name=name,
        api_version=api_version,
        select=_select,
        load=_load,
    )


def _write_legacy_library(root, name='localtoy'):
    library_dir = root / name
    library_dir.mkdir()
    tools_file = library_dir / f'{name}_tools.py'
    tools_file.write_text(
        '\n'.join([
            f'def {name}_select(options):',
            "    return ['TOY']",
            '',
            f'def {name}_load(problem_name):',
            '    return problem_name',
            '',
            f'def {name}_collect_info():',
            "    return 'collected'",
            '',
            f'def {name}_check_available():',
            "    return 'available'",
        ]),
        encoding='utf-8',
    )
    return library_dir, tools_file


def _use_entry_point_for_loading(monkeypatch, entry_point):
    monkeypatch.setattr(
        problem_libraries,
        '_entry_point_from_reference',
        lambda reference: entry_point,
    )


def test_protocol_api_is_public():
    assert optiprofiler.ProblemLibraryPlugin is ProblemLibraryPlugin
    assert optiprofiler.list_problem_libraries is list_problem_libraries


def test_list_problem_libraries_does_not_import_entry_points(monkeypatch):
    entry_point = _FakeEntryPoint(
        'externaltoy',
        'optiprofiler_externaltoy:get_problem_library',
        lambda: _plugin(),
    )
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: [entry_point],
    )

    names = list_problem_libraries()

    assert 's2mpj' in names
    assert 'externaltoy' in names
    assert entry_point.load_calls == 0


def test_python_38_entry_point_mapping_is_supported(monkeypatch):
    entry_point = _FakeEntryPoint(
        'externaltoy',
        'optiprofiler_externaltoy:get_problem_library',
        lambda: _plugin(),
    )
    monkeypatch.setattr(
        problem_libraries.metadata,
        'entry_points',
        lambda: {PROBLEM_LIBRARY_ENTRY_POINT_GROUP: [entry_point]},
    )

    assert 'externaltoy' in list_problem_libraries()
    assert resolve_problem_library('externaltoy').source == 'entry_point'
    assert entry_point.load_calls == 0


def test_invalid_entry_point_name_is_not_listed(monkeypatch):
    entry_point = _FakeEntryPoint(
        '../externaltoy',
        'optiprofiler_externaltoy:get_problem_library',
        lambda: _plugin(),
    )
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: [entry_point],
    )

    assert '../externaltoy' not in list_problem_libraries()


def test_invalid_custom_library_name_is_not_listed(tmp_path):
    library_dir = tmp_path / 'local-toy'
    library_dir.mkdir()
    (library_dir / 'local-toy_tools.py').write_text('', encoding='utf-8')

    assert 'local-toy' not in list_problem_libraries(tmp_path)


def test_problem_option_validation_does_not_import_entry_points(monkeypatch):
    entry_point = _FakeEntryPoint(
        'externaltoy',
        'optiprofiler_externaltoy:get_problem_library',
        lambda: _plugin(),
    )
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: [entry_point],
    )

    options = check_validity_problem_options({'plibs': ['externaltoy']})

    assert options['plibs'] == ['externaltoy']
    assert entry_point.load_calls == 0


def test_resolve_and_load_entry_point(monkeypatch):
    entry_point = _FakeEntryPoint(
        'externaltoy',
        'optiprofiler_externaltoy:get_problem_library',
        lambda: _plugin(),
    )
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: [entry_point],
    )
    _use_entry_point_for_loading(monkeypatch, entry_point)

    reference = resolve_problem_library('externaltoy')
    roundtripped = pickle.loads(pickle.dumps(reference))
    plugin = load_problem_library(roundtripped)

    assert reference.source == 'entry_point'
    assert reference.distribution == 'test-distribution'
    assert plugin.name == 'externaltoy'
    assert plugin.select({}, {}) == ['TOY']
    assert plugin.load('TOY', {}) == 'TOY'
    assert entry_point.load_calls == 1


def test_resolved_entry_point_load_does_not_rescan_packages(monkeypatch):
    entry_point = _FakeEntryPoint(
        'externaltoy',
        'optiprofiler_externaltoy:get_problem_library',
        lambda: _plugin(),
    )
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: [entry_point],
    )
    reference = resolve_problem_library('externaltoy')
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: pytest.fail('loading a resolved entry point must not rescan packages'),
    )
    _use_entry_point_for_loading(monkeypatch, entry_point)

    assert load_problem_library(reference).name == 'externaltoy'


def test_duplicate_entry_points_are_rejected(monkeypatch):
    entry_points = [
        _FakeEntryPoint('externaltoy', 'package_a:get_plugin', lambda: _plugin(), 'package-a'),
        _FakeEntryPoint('externaltoy', 'package_b:get_plugin', lambda: _plugin(), 'package-b'),
    ]
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: entry_points,
    )

    with pytest.raises(ValueError, match='multiple installed distributions'):
        resolve_problem_library('externaltoy')


def test_bundled_and_installed_name_collision_is_rejected(monkeypatch):
    entry_point = _FakeEntryPoint(
        's2mpj',
        'optiprofiler_s2mpj:get_problem_library',
        lambda: _plugin(name='s2mpj'),
    )
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: [entry_point],
    )

    with pytest.raises(ValueError, match='both bundled and installed'):
        resolve_problem_library('s2mpj')


def test_explicit_custom_path_overrides_other_providers(tmp_path, monkeypatch):
    _, tools_file = _write_legacy_library(tmp_path, 's2mpj')
    entry_point = _FakeEntryPoint(
        's2mpj',
        'optiprofiler_s2mpj:get_problem_library',
        lambda: _plugin(name='s2mpj'),
    )
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: [entry_point],
    )

    reference = resolve_problem_library('s2mpj', tmp_path)

    assert reference.source == 'custom'
    assert Path(reference.locator) == tools_file.resolve()
    assert reference.shadows == ('bundled', 'test-distribution')


def test_distribution_provenance_reads_version_and_vcs_commit(monkeypatch):
    commit = 'a' * 40

    class _Distribution:
        version = '1.2.3'

        def read_text(self, name):
            if name == 'direct_url.json':
                return json.dumps({'url': 'https://example.invalid/toy.git',
                                   'vcs_info': {'vcs': 'git', 'commit_id': commit}})
            return None

    def distribution(name):
        if name == 'missing-dist':
            raise problem_libraries.metadata.PackageNotFoundError(name)
        return _Distribution()

    monkeypatch.setattr(problem_libraries.metadata, 'distribution', distribution)
    provenance = problem_libraries.distribution_provenance

    assert provenance(ProblemLibraryRef('toy', 'entry_point', 'toy:get', distribution='toy-dist')) == ('1.2.3', commit)
    # A bundled library is identified by the imported core's version only;
    # the core's own VCS commit is not the bundled library's revision.
    assert provenance(ProblemLibraryRef('s2mpj', 'builtin', 'x', distribution='optiprofiler')) == (optiprofiler.__version__, None)
    assert provenance(ProblemLibraryRef('local', 'custom', 'x')) == (None, None)
    assert provenance(ProblemLibraryRef('toy', 'entry_point', 'toy:get', distribution='missing-dist')) == (None, None)


def test_unique_custom_library_shadows_nothing(tmp_path):
    _write_legacy_library(tmp_path)

    assert resolve_problem_library('localtoy', tmp_path).shadows == ()
    assert resolve_problem_library('s2mpj').shadows == ()


def _write_quadratic_library(root, name):
    library_dir = root / name
    library_dir.mkdir()
    (library_dir / f'{name}_tools.py').write_text(
        '\n'.join([
            'import numpy as np',
            'from optiprofiler import Problem',
            f'def {name}_select(options):',
            "    return ['LOCALQ']",
            f'def {name}_load(problem_name):',
            '    return Problem(lambda x: float(np.sum((x - 1.0) ** 2)), [0.0, 0.0], name=problem_name)',
        ]),
        encoding='utf-8',
    )


def _one_step_solver(fun, x0):
    fun(x0)
    return x0


@pytest.mark.parametrize('name, warned', [('s2mpj', True), ('custom', False)])
def test_benchmark_reports_a_local_library_that_replaces_a_provider(tmp_path, caplog, name, warned):
    _write_quadratic_library(tmp_path, name)
    report_path = tmp_path / 'report.json'

    optiprofiler.benchmark([_one_step_solver, _one_step_solver], plibs=[name], custom_problem_libs_path=str(tmp_path),
                           n_jobs=1, n_runs=1, max_eval_factor=2, max_tol_order=1, run_plain=True,
                           score_only=True, silent=True, report_path=str(report_path))

    report = json.loads(report_path.read_text(encoding='utf-8'))
    codes = [entry['code'] for entry in report['diagnostics']]
    assert {problem['provider']['source'] for problem in report['problems']} == {'custom'}
    assert {tuple(problem['provider']['shadows']) for problem in report['problems']} == {('bundled',)}
    assert {(problem['provider']['version'], problem['provider']['revision']) for problem in report['problems']} == {(None, None)}
    warnings = [record.getMessage() for record in caplog.records
                if record.levelname == 'WARNING' and 'replaces the provider' in record.getMessage()]
    assert len(warnings) == (1 if warned else 0)
    assert all('custom_problem_libs_path' in message and f'"{name}"' in message for message in warnings)
    assert codes.count('problem_library_shadowed') == (1 if warned else 0)


def test_load_legacy_custom_library_from_parent_and_direct_path(tmp_path):
    library_dir, tools_file = _write_legacy_library(tmp_path)

    parent_ref = resolve_problem_library('localtoy', tmp_path)
    direct_ref = resolve_problem_library('localtoy', library_dir)
    plugin = load_problem_library(parent_ref)

    assert Path(parent_ref.locator) == tools_file.resolve()
    assert direct_ref == parent_ref
    assert plugin.select({}, {}) == ['TOY']
    assert plugin.load('TOY', {}) == 'TOY'
    assert plugin.collect_info() == 'collected'
    assert plugin.check_available() == 'available'


def test_legacy_get_info_alias_is_adapted(tmp_path):
    library_dir = tmp_path / 'localtoy'
    library_dir.mkdir()
    (library_dir / 'localtoy_tools.py').write_text(
        '\n'.join([
            'def localtoy_select(options):',
            "    return ['TOY']",
            '',
            'def localtoy_load(problem_name):',
            '    return problem_name',
            '',
            'def localtoy_get_info():',
            "    return 'legacy info'",
        ]),
        encoding='utf-8',
    )

    plugin = load_problem_library(resolve_problem_library('localtoy', tmp_path))

    assert plugin.collect_info() == 'legacy info'


def test_ambiguous_custom_path_is_rejected(tmp_path):
    direct_dir, _ = _write_legacy_library(tmp_path, 'localtoy')
    _write_legacy_library(direct_dir, 'localtoy')

    with pytest.raises(ValueError, match='ambiguous'):
        resolve_problem_library('localtoy', direct_dir)


def test_custom_library_filename_is_strict(tmp_path):
    library_dir = tmp_path / 'localtoy'
    library_dir.mkdir()
    (library_dir / 'other_tools.py').write_text('', encoding='utf-8')

    with pytest.raises(ValueError, match='must contain "localtoy_tools.py"'):
        resolve_problem_library('localtoy', library_dir)


def test_public_discovery_validates_custom_path(tmp_path):
    missing_path = tmp_path / 'missing'

    with pytest.raises(ValueError, match='does not exist'):
        list_problem_libraries(missing_path)
    with pytest.raises(ValueError, match='does not exist'):
        resolve_problem_library('localtoy', missing_path)


@pytest.mark.parametrize(
    'name',
    [
        None,
        42,
        '',
        '2localtoy',
        '../localtoy',
        'local/toy',
        r'local\toy',
        '.localtoy',
        'local-toy',
        'local.toy',
        'local toy',
    ],
)
def test_problem_library_name_is_path_safe(name):
    with pytest.raises(ValueError, match='ASCII Python identifier'):
        resolve_problem_library(name)


@pytest.mark.parametrize(
    ('factory', 'error_type', 'message'),
    [
        (lambda: object(), TypeError, 'must return a ProblemLibraryPlugin'),
        (lambda: _plugin(name='wrong'), ValueError, 'returned a plugin named'),
        (lambda: _plugin(api_version=999), ValueError, 'uses API version'),
        (lambda: _plugin(api_version=True), ValueError, 'uses API version'),
        (lambda: _plugin(api_version='1'), ValueError, 'uses API version'),
    ],
)
def test_entry_point_plugin_contract_is_validated(
    monkeypatch, factory, error_type, message
):
    entry_point = _FakeEntryPoint(
        'externaltoy',
        'optiprofiler_externaltoy:get_problem_library',
        factory,
    )
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: [entry_point],
    )
    _use_entry_point_for_loading(monkeypatch, entry_point)

    reference = resolve_problem_library('externaltoy')
    with pytest.raises(error_type, match=message):
        load_problem_library(reference)


def test_entry_point_must_resolve_to_factory(monkeypatch):
    entry_point = _FakeEntryPoint(
        'externaltoy',
        'optiprofiler_externaltoy:plugin',
        _plugin(),
    )
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: [entry_point],
    )
    _use_entry_point_for_loading(monkeypatch, entry_point)

    reference = resolve_problem_library('externaltoy')
    with pytest.raises(TypeError, match='zero-argument factory'):
        load_problem_library(reference)


def test_plugin_functions_must_be_callable(monkeypatch):
    invalid_plugin = ProblemLibraryPlugin(
        name='externaltoy',
        api_version=PROBLEM_LIBRARY_API_VERSION,
        select=None,
        load=_load,
    )
    entry_point = _FakeEntryPoint(
        'externaltoy',
        'optiprofiler_externaltoy:get_problem_library',
        lambda: invalid_plugin,
    )
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: [entry_point],
    )
    _use_entry_point_for_loading(monkeypatch, entry_point)

    reference = resolve_problem_library('externaltoy')
    with pytest.raises(TypeError, match='callable select and load'):
        load_problem_library(reference)


@pytest.mark.parametrize(
    'field',
    [
        'collect_info',
        'check_available',
        'get_default_options',
        'validate_options',
    ],
)
def test_optional_plugin_functions_must_be_callable(monkeypatch, field):
    kwargs = {field: 'not callable'}
    invalid_plugin = ProblemLibraryPlugin(
        name='externaltoy',
        api_version=PROBLEM_LIBRARY_API_VERSION,
        select=_select,
        load=_load,
        **kwargs,
    )
    entry_point = _FakeEntryPoint(
        'externaltoy',
        'optiprofiler_externaltoy:get_problem_library',
        lambda: invalid_plugin,
    )
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: [entry_point],
    )
    _use_entry_point_for_loading(monkeypatch, entry_point)

    reference = resolve_problem_library('externaltoy')
    with pytest.raises(TypeError, match=field):
        load_problem_library(reference)


@pytest.mark.parametrize(
    'config_callbacks',
    [
        {'get_default_options': lambda: {}},
        {'validate_options': lambda options: options},
    ],
)
def test_plugin_must_declare_both_config_callbacks(monkeypatch, config_callbacks):
    invalid_plugin = ProblemLibraryPlugin(
        name='externaltoy',
        api_version=PROBLEM_LIBRARY_API_VERSION,
        select=_select,
        load=_load,
        **config_callbacks,
    )
    entry_point = _FakeEntryPoint(
        'externaltoy',
        'optiprofiler_externaltoy:get_problem_library',
        lambda: invalid_plugin,
    )
    _use_entry_point_for_loading(monkeypatch, entry_point)

    reference = ProblemLibraryRef(
        'externaltoy',
        'entry_point',
        entry_point.value,
    )
    with pytest.raises(TypeError, match='must provide both'):
        load_problem_library(reference)


def test_problem_library_options_are_merged_and_normalized():
    plugin = ProblemLibraryPlugin(
        name='externaltoy',
        api_version=PROBLEM_LIBRARY_API_VERSION,
        select=_select,
        load=_load,
        get_default_options=lambda: {'variant': 1, 'enabled': True},
        validate_options=lambda options: {
            'variant': int(options['variant']),
            'enabled': bool(options['enabled']),
        },
    )

    options = _resolve_problem_library_options(plugin, {'variant': '3'})

    assert options == {'variant': 3, 'enabled': True}


def test_plugin_validation_cannot_mutate_caller_options():
    def validate(options):
        options['trace'].append('validated')
        return options

    plugin = ProblemLibraryPlugin(
        name='externaltoy',
        api_version=PROBLEM_LIBRARY_API_VERSION,
        select=_select,
        load=_load,
        get_default_options=lambda: {'trace': []},
        validate_options=validate,
    )
    overrides = {'trace': []}

    resolved = _resolve_problem_library_options(plugin, overrides)

    assert overrides == {'trace': []}
    assert resolved == {'trace': ['validated']}


def test_nonconfigurable_plugin_rejects_nonempty_options():
    with pytest.raises(ValueError, match='does not declare configurable options'):
        _resolve_problem_library_options(_plugin(), {'variant': 1})


@pytest.mark.parametrize(
    ('get_defaults', 'validator', 'error_type', 'message'),
    [
        (lambda: [], None, TypeError, 'must be a mapping'),
        (lambda: {1: 'bad'}, None, TypeError, 'string keys'),
        (lambda: {}, lambda options: [], TypeError, 'must be a mapping'),
        (lambda: {}, lambda options: {1: 'bad'}, TypeError, 'string keys'),
        (lambda: {}, lambda options: {'bad': lambda: None}, TypeError, 'pickleable'),
    ],
)
def test_problem_library_option_contract_is_validated(
    get_defaults,
    validator,
    error_type,
    message,
):
    plugin = ProblemLibraryPlugin(
        name='externaltoy',
        api_version=PROBLEM_LIBRARY_API_VERSION,
        select=_select,
        load=_load,
        get_default_options=get_defaults,
        validate_options=validator,
    )

    with pytest.raises(error_type, match=message):
        _resolve_problem_library_options(plugin)


def test_load_rejects_invalid_reference():
    with pytest.raises(TypeError, match='ProblemLibraryRef'):
        load_problem_library('externaltoy')

    reference = ProblemLibraryRef('externaltoy', 'unknown', 'nowhere')
    with pytest.raises(ValueError, match='Unknown problem-library source'):
        load_problem_library(reference)
