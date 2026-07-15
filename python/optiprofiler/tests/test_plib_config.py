"""Tests for problem-library configuration resolution."""

import pytest

from optiprofiler import get_plib_config, set_plib_config
from optiprofiler import plib_config, problem_libraries
from optiprofiler.plib_config import _resolve_plib_options
from optiprofiler.problem_libraries import ProblemLibraryPlugin, ProblemLibraryRef


@pytest.fixture(autouse=True)
def _clean_process_overrides(monkeypatch):
    # These are unit tests for the bundled and explicit custom providers. Keep
    # them independent of optional problem-library distributions that may be
    # installed in the developer's environment; provider-collision behavior is
    # covered separately in test_problem_libraries.py.
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: [],
    )
    plib_config._PLIB_CONFIG_OVERRIDES.clear()
    monkeypatch.delenv('S2MPJ_VARIABLE_SIZE', raising=False)
    monkeypatch.delenv('S2MPJ_TEST_FEASIBILITY_PROBLEMS', raising=False)
    monkeypatch.delenv('PYCUTEST_VARIABLE_SIZE', raising=False)
    monkeypatch.delenv('PYCUTEST_TEST_FEASIBILITY_PROBLEMS', raising=False)
    yield
    plib_config._PLIB_CONFIG_OVERRIDES.clear()


def _write_configurable_library(root, default_variant):
    library_dir = root / 'localtoy'
    library_dir.mkdir(parents=True)
    (library_dir / 'localtoy_tools.py').write_text(
        '\n'.join([
            'def localtoy_get_default_options():',
            f"    return {{'variant': {default_variant}}}",
            '',
            'def localtoy_validate_options(options):',
            "    if set(options) != {'variant'}:",
            "        raise ValueError('variant is required')",
            "    return {'variant': int(options['variant'])}",
            '',
            'def localtoy_select(problem_options, library_options):',
            "    return ['TOY']",
            '',
            'def localtoy_load(problem_name, *, library_options=None):',
            '    return problem_name',
        ]),
        encoding='utf-8',
    )
    return library_dir


def test_get_plib_config_reads_effective_s2mpj_defaults():
    assert get_plib_config('s2mpj') == {
        'variable_size': 'default',
        'test_feasibility_problems': 0,
    }


def test_get_pycutest_config_does_not_require_the_optional_runtime(monkeypatch):
    monkeypatch.setattr(
        plib_config,
        'load_problem_library',
        lambda reference: pytest.fail(
            'reading built-in config must not import the adapter'
        ),
    )

    assert get_plib_config('pycutest') == {
        'variable_size': 'default',
        'test_feasibility_problems': 0,
    }
    set_plib_config('pycutest', variable_size='all')
    assert get_plib_config('pycutest')['variable_size'] == 'all'


def test_set_plib_config_remains_a_process_level_default():
    set_plib_config(
        's2mpj',
        variable_size='all',
        test_feasibility_problems=2,
    )

    assert get_plib_config('s2mpj') == {
        'variable_size': 'all',
        'test_feasibility_problems': 2,
    }
    assert plib_config.os.environ['S2MPJ_VARIABLE_SIZE'] == 'all'
    assert plib_config.os.environ['S2MPJ_TEST_FEASIBILITY_PROBLEMS'] == '2'


def test_per_run_options_override_session_defaults_without_mutating_them(
    tmp_path,
):
    library = _write_configurable_library(tmp_path, 1)
    set_plib_config(
        'localtoy',
        custom_problem_libs_path=library,
        variant=2,
    )

    resolved = _resolve_plib_options(
        'localtoy',
        {'variant': 3},
        library,
    )

    assert resolved == {'variant': 3}
    assert get_plib_config(
        'localtoy',
        custom_problem_libs_path=library,
    ) == {'variant': 2}


def test_per_run_options_mask_an_invalid_library_default(tmp_path):
    library = _write_configurable_library(tmp_path, repr('invalid'))

    resolved = _resolve_plib_options(
        'localtoy',
        {'variant': 4},
        library,
    )

    assert resolved == {'variant': 4}


@pytest.mark.parametrize(
    'kwargs',
    [
        {'unknown': 1},
        {'variant': 'not-an-integer'},
    ],
)
def test_set_plib_config_uses_library_validation(tmp_path, kwargs):
    library = _write_configurable_library(tmp_path, 1)
    with pytest.raises(ValueError):
        set_plib_config(
            'localtoy',
            custom_problem_libs_path=library,
            **kwargs,
        )


def test_process_overrides_are_scoped_to_the_resolved_provider(tmp_path):
    first_root = tmp_path / 'first'
    second_root = tmp_path / 'second'
    first_library = _write_configurable_library(first_root, 1)
    second_library = _write_configurable_library(second_root, 2)

    set_plib_config(
        'localtoy',
        custom_problem_libs_path=first_library,
        variant=10,
    )

    assert get_plib_config(
        'localtoy',
        custom_problem_libs_path=first_library,
    ) == {'variant': 10}
    assert get_plib_config(
        'localtoy',
        custom_problem_libs_path=second_library,
    ) == {'variant': 2}


def test_custom_provider_override_is_not_mirrored_to_name_only_environment(
    tmp_path,
):
    library = _write_configurable_library(tmp_path, 1)
    config_path = library / 'config.txt'
    config_path.write_text('variant=1\n', encoding='utf-8')

    set_plib_config(
        'localtoy',
        custom_problem_libs_path=library,
        variant=10,
    )

    assert 'LOCALTOY_VARIANT' not in plib_config.os.environ


def test_legacy_file_config_remains_available(tmp_path):
    library_dir = tmp_path / 'legacytoy'
    library_dir.mkdir()
    (library_dir / 'legacytoy_tools.py').write_text(
        '\n'.join([
            'def legacytoy_select(options):',
            "    return ['TOY']",
            '',
            'def legacytoy_load(problem_name):',
            '    return problem_name',
        ]),
        encoding='utf-8',
    )
    (library_dir / 'config.txt').write_text('variant=1\n', encoding='utf-8')

    assert get_plib_config(
        'legacytoy',
        custom_problem_libs_path=library_dir,
    ) == {'variant': 1}

    with pytest.raises(ValueError, match='configuration callbacks'):
        set_plib_config(
            'legacytoy',
            custom_problem_libs_path=library_dir,
            variant=3,
        )
    with pytest.raises(ValueError, match='does not support explicit'):
        _resolve_plib_options(
            'legacytoy',
            {'variant': 4},
            library_dir,
        )


def test_nonconfigurable_installed_plugin_has_an_empty_config(monkeypatch):
    reference = ProblemLibraryRef(
        'externaltoy',
        'entry_point',
        'optiprofiler_externaltoy:get_problem_library',
        distribution='optiprofiler-externaltoy',
    )
    plugin = ProblemLibraryPlugin(
        name='externaltoy',
        api_version=1,
        select=lambda problem_options, library_options: ['TOY'],
        load=lambda problem_name, library_options: problem_name,
    )
    monkeypatch.setattr(
        plib_config,
        'resolve_problem_library',
        lambda name, custom_path=None: reference,
    )
    monkeypatch.setattr(
        plib_config,
        'load_problem_library',
        lambda resolved_reference: plugin,
    )

    assert get_plib_config('externaltoy') == {}
    assert set_plib_config('externaltoy') == {}
    with pytest.raises(ValueError, match='does not declare configurable options'):
        set_plib_config('externaltoy', unknown=1)


def test_process_overrides_are_isolated_from_nested_caller_mutation(monkeypatch):
    reference = ProblemLibraryRef(
        'externaltoy',
        'entry_point',
        'optiprofiler_externaltoy:get_problem_library',
        distribution='optiprofiler-externaltoy',
    )
    plugin = ProblemLibraryPlugin(
        name='externaltoy',
        api_version=1,
        select=lambda problem_options, library_options: ['TOY'],
        load=lambda problem_name, library_options: problem_name,
        get_default_options=lambda: {'trace': []},
        validate_options=lambda options: {'trace': list(options['trace'])},
    )
    monkeypatch.setattr(
        plib_config,
        'resolve_problem_library',
        lambda name, custom_path=None: reference,
    )
    monkeypatch.setattr(
        plib_config,
        'load_problem_library',
        lambda resolved_reference: plugin,
    )
    trace = []

    set_plib_config('externaltoy', trace=trace)
    trace.append('caller mutation')
    first_read = get_plib_config('externaltoy')
    first_read['trace'].append('reader mutation')

    assert get_plib_config('externaltoy') == {'trace': []}
