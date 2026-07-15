"""Tests for the benchmark function and solve helpers in profiles module."""
import os
import logging
import pickle
import shutil
import tempfile
from pathlib import Path
from unittest.mock import patch

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pytest

import optiprofiler.profiles as profiles_module
import optiprofiler.profile_utils as profile_utils_module
from optiprofiler import benchmark
from optiprofiler.opclasses import Feature, FeaturedProblem, Problem
from optiprofiler.problem_libraries import ProblemLibraryPlugin
from optiprofiler import plib_config, problem_libraries
from optiprofiler.loader import load_results_from_h5
from optiprofiler.profile_utils import get_default_problem_options, get_default_profile_options
from optiprofiler.profiles import _format_feature_title, _format_problem_feature_title, _resolve_benchmark_plib_options, _solve_all_problems
from optiprofiler.utils import ProblemOption


def simple_solver_1(fun, x0):
    """A trivial solver that returns x0."""
    return x0


def simple_solver_2(fun, x0):
    """A trivial solver that does one step of coordinate search."""
    n = len(x0)
    best_x = x0.copy()
    best_f = fun(x0)
    for i in range(n):
        for step in [0.5, -0.5, 0.1, -0.1]:
            x_trial = best_x.copy()
            x_trial[i] += step
            f_trial = fun(x_trial)
            if f_trial < best_f:
                best_f = f_trial
                best_x = x_trial.copy()
    return best_x


def simple_solver_zero_no_eval(fun, x0):
    """A trivial solver that returns the zero vector without evaluating the objective."""
    return np.zeros_like(x0)


def history_solver_1(fun, x0, *args):
    """A solver that records at least one history entry for any problem type."""
    fun(x0)
    return x0


def history_solver_2(fun, x0, *args):
    """A solver that records two history entries for any problem type."""
    fun(x0)
    x = x0.copy()
    x[0] += 0.1
    fun(x)
    return x


def externaltoy_select(options, library_options):
    return ['EXTERNAL_TOY']


def externaltoy_load(problem_name, library_options):
    return Problem(
        lambda x: float(np.dot(x, x)),
        np.array([1.0, -1.0]),
        name=problem_name,
    )


def externaltoy_get_default_options():
    return {'variant': 'default', 'include_special': 0}


def externaltoy_validate_options(options):
    known = {'variant', 'include_special'}
    unknown = sorted(set(options) - known)
    if unknown:
        raise ValueError(f'Unknown externaltoy library options: {unknown}')
    variant = options.get('variant', 'default')
    if variant not in {'default', 'min', 'max', 'all'}:
        raise ValueError('Invalid externaltoy variant')
    include_special = options.get('include_special', 0)
    if include_special not in {0, 1}:
        raise ValueError('Invalid externaltoy include_special value')
    return {
        'variant': variant,
        'include_special': include_special,
    }


def externaltoy_plugin_factory():
    return ProblemLibraryPlugin(
        name='externaltoy',
        api_version=1,
        select=externaltoy_select,
        load=externaltoy_load,
        get_default_options=externaltoy_get_default_options,
        validate_options=externaltoy_validate_options,
    )


class ExternalToyEntryPoint:
    name = 'externaltoy'
    value = 'optiprofiler_externaltoy:get_problem_library'
    dist = None

    @staticmethod
    def load():
        return externaltoy_plugin_factory


def use_externaltoy_entry_point(monkeypatch, entry_point=None):
    if entry_point is None:
        entry_point = ExternalToyEntryPoint()
    monkeypatch.setattr(
        problem_libraries,
        '_problem_library_entry_points',
        lambda: [entry_point],
    )
    monkeypatch.setattr(
        problem_libraries,
        '_entry_point_from_reference',
        lambda reference: entry_point,
    )
    return entry_point


def _common_kwargs(tmpdir, benchmark_id='test'):
    return dict(
        plibs=['s2mpj'],
        ptype='u',
        mindim=2,
        maxdim=2,
        max_eval_factor=10,
        benchmark_id=benchmark_id,
        savepath=tmpdir,
        n_jobs=1,
        silent=True,
        draw_hist_plots='none',
        problem_names=['ROSENBR'],
    )


class TestBenchmarkBasic:

    def _assert_history_plot_files(self, tmpdir, problem_name='ROSENBR'):
        root = Path(tmpdir)
        assert len(list(root.glob(f'**/history_plots/s2mpj/{problem_name}.pdf'))) == 1
        assert len(list(root.glob(f'**/history_plots/s2mpj/raw/{problem_name}.pdf'))) == 1
        assert len(list(root.glob(f'**/history_plots/s2mpj/cummin/{problem_name}.pdf'))) == 1
        assert len(list(root.glob('**/history_plots/s2mpj_history_plots_summary.pdf'))) == 1

    def test_two_solvers_plain(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, profile_scores, curves = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                **_common_kwargs(tmpdir),
            )
            assert isinstance(scores, np.ndarray)
            assert scores.shape[0] == 2

    def test_plib_options_are_refined_and_stored(self, monkeypatch):
        use_externaltoy_entry_point(monkeypatch)
        plib_config._PLIB_CONFIG_OVERRIDES.clear()
        with tempfile.TemporaryDirectory() as tmpdir:
            kwargs = _common_kwargs(tmpdir, benchmark_id='plib-options')
            kwargs.update({
                'plibs': ['externaltoy'],
                'problem_names': ['EXTERNAL_TOY'],
            })
            benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                plib_options={
                    'externaltoy': {
                        'variant': 'max',
                    },
                },
                **kwargs,
            )

            log_dir = next(Path(tmpdir).glob('plib-options/*/test_log'))
            with (log_dir / 'options_user.pkl').open('rb') as f:
                options_user = pickle.load(f)
            with (log_dir / 'options_refined.pkl').open('rb') as f:
                options_refined = pickle.load(f)
            results = load_results_from_h5(log_dir / 'data_for_loading.h5')

            assert options_user['plib_options'] == {
                'externaltoy': {'variant': 'max'},
            }
            expected = {
                'variant': 'max',
                'include_special': 0,
            }
            assert options_refined['plib_options'] == {'externaltoy': expected}
            assert results[0]['plib_options'] == expected

    def test_plib_options_for_unselected_library_are_rejected(self, tmp_path):
        kwargs = _common_kwargs(str(tmp_path), benchmark_id='unselected-options')
        with pytest.raises(ValueError, match='unselected problem libraries'):
            benchmark(
                [simple_solver_1, simple_solver_2],
                plib_options={'pycutest': {'variable_size': 'all'}},
                **kwargs,
            )

    def test_unknown_plib_option_is_rejected_by_the_library(
        self,
        tmp_path,
        monkeypatch,
    ):
        use_externaltoy_entry_point(monkeypatch)
        kwargs = _common_kwargs(str(tmp_path), benchmark_id='unknown-options')
        kwargs.update({
            'plibs': ['externaltoy'],
            'problem_names': ['EXTERNAL_TOY'],
        })
        with pytest.raises(ValueError, match='Unknown externaltoy library options'):
            benchmark(
                [simple_solver_1, simple_solver_2],
                plib_options={'externaltoy': {'varient': 'all'}},
                **kwargs,
            )
        assert not (tmp_path / 'unknown-options').exists()

    def test_public_benchmark_closes_logging_on_plugin_error(
        self,
        tmp_path,
        monkeypatch,
    ):
        def failing_select(options, library_options):
            raise RuntimeError('selection failed')

        plugin = ProblemLibraryPlugin(
            name='externaltoy',
            api_version=1,
            select=failing_select,
            load=externaltoy_load,
        )
        entry_point = ExternalToyEntryPoint()
        monkeypatch.setattr(entry_point, 'load', lambda: lambda: plugin)
        use_externaltoy_entry_point(monkeypatch, entry_point)

        class FakeListener:
            handlers = []

            def __init__(self):
                self.stopped = 0

            def stop(self):
                self.stopped += 1

        class FakeQueue:
            def __init__(self):
                self.closed = 0
                self.joined = 0

            def close(self):
                self.closed += 1

            def join_thread(self):
                self.joined += 1

        listener = FakeListener()
        log_queue = FakeQueue()

        class FakeQueueHandler(logging.Handler):
            def __init__(self, queue):
                super().__init__()
                self.queue = queue

            def emit(self, record):
                pass

        queue_handler = FakeQueueHandler(log_queue)
        logging.getLogger().addHandler(queue_handler)
        monkeypatch.setattr(
            profiles_module,
            'setup_main_process_logging',
            lambda **kwargs: (log_queue, listener),
        )
        kwargs = _common_kwargs(str(tmp_path), benchmark_id='plugin-error')
        kwargs.update({
            'plibs': ['externaltoy'],
            'problem_names': ['EXTERNAL_TOY'],
        })

        with pytest.raises(RuntimeError, match='selection failed'):
            benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                **kwargs,
            )

        assert listener.stopped == 1
        assert log_queue.closed == 1
        assert log_queue.joined == 1
        assert queue_handler not in logging.getLogger().handlers

    def test_loading_saved_results_does_not_require_plugin(
        self,
        tmp_path,
        monkeypatch,
    ):
        use_externaltoy_entry_point(monkeypatch)
        kwargs = _common_kwargs(str(tmp_path), benchmark_id='retired-plugin')
        kwargs.update({
            'plibs': ['externaltoy'],
            'problem_names': ['EXTERNAL_TOY'],
        })
        benchmark(
            [simple_solver_1, simple_solver_2],
            feature_name='plain',
            **kwargs,
        )
        monkeypatch.setattr(
            profile_utils_module,
            'resolve_problem_library',
            lambda *args, **kwargs: pytest.fail(
                'loading saved results must not resolve a provider'
            ),
        )
        monkeypatch.chdir(tmp_path)

        scores, profile_scores, curves = benchmark(
            None,
            load='latest',
            plibs=['externaltoy'],
            savepath=str(tmp_path),
            benchmark_id='retired-plugin',
            score_only=True,
            silent=True,
            draw_hist_plots='none',
        )

        assert scores.shape == (2,)
        assert profile_scores is not None
        assert curves is not None

    def test_plib_options_are_rejected_for_direct_problem(self, tmp_path):
        problem = Problem(lambda x: float(np.dot(x, x)), np.ones(2), name='TOY')
        with pytest.raises(ValueError, match='cannot be used with the `problem`'):
            benchmark(
                [simple_solver_1, simple_solver_2],
                problem=problem,
                plib_options={'s2mpj': {'variable_size': 'all'}},
                savepath=str(tmp_path),
                benchmark_id='direct-problem-options',
                silent=True,
            )
        assert not (tmp_path / 'direct-problem-options').exists()

    def test_plib_options_are_rejected_when_loading_saved_results(self, tmp_path):
        with pytest.raises(ValueError, match='cannot be used with the `load`'):
            benchmark(
                None,
                load='latest',
                plib_options={'s2mpj': {'variable_size': 'all'}},
                savepath=str(tmp_path),
                benchmark_id='load-options',
                silent=True,
            )
        assert not (tmp_path / 'load-options').exists()

    def test_returns_curves(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, profile_scores, curves = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                **_common_kwargs(tmpdir),
            )
            assert curves is not None
            assert isinstance(curves, list)
            assert len(curves) > 0

    def test_silent_false(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, profile_scores, curves = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                plibs=['s2mpj'],
                ptype='u',
                mindim=2,
                maxdim=2,
                max_eval_factor=10,
                benchmark_id='test',
                savepath=tmpdir,
                n_jobs=1,
                silent=False,
                draw_hist_plots='none',
                problem_names=['ROSENBR'],
            )
            assert isinstance(scores, np.ndarray)
            assert scores.shape[0] == 2

    def test_with_problem_option(self):
        def rosen(x):
            return np.sum(1e2 * (x[1:] - x[:-1] ** 2) ** 2 + (1.0 - x[:-1]) ** 2)
        problem = Problem(rosen, np.zeros(2), name='rosen')
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, profile_scores, curves = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                problem=problem,
                savepath=tmpdir,
                benchmark_id='test',
                n_jobs=1,
                silent=True,
                draw_hist_plots='none',
            )
            assert isinstance(scores, np.ndarray)

    @pytest.mark.parametrize('n_jobs', [1, 2])
    def test_custom_problem_library_with_direct_path(self, tmp_path, n_jobs):
        lib_dir = tmp_path / 'solar'
        lib_dir.mkdir()
        (lib_dir / 'solar_tools.py').write_text(
            "\n".join([
                "import numpy as np",
                "from optiprofiler.opclasses import Problem",
                "",
                "def solar_select(options):",
                "    return ['SOLAR_TOY']",
                "",
                "def solar_load(problem_name):",
                "    return Problem(lambda x: float(np.dot(x, x)), np.array([1.0, -1.0]), name=problem_name)",
            ]),
            encoding='utf-8',
        )

        solvers = [simple_solver_1, simple_solver_zero_no_eval]
        feature = Feature('plain')
        problem_options = get_default_problem_options(
            {
                'plibs': ['solar'],
                'custom_problem_libs_path': lib_dir,
                'problem_names': ['SOLAR_TOY'],
            }
        )
        profile_options = get_default_profile_options(
            solvers,
            feature,
            {
                'max_eval_factor': 2,
                'n_jobs': n_jobs,
                'silent': True,
                'solver_names': ['solver1', 'solver2'],
            },
        )
        results = _solve_all_problems(
            solvers,
            'solar',
            feature,
            problem_options,
            profile_options,
            False,
            None,
        )

        assert results['plib'] == 'solar'
        assert results['problem_names'] == ['SOLAR_TOY']
        assert results['fun_outs'].shape == (1, 2, 1)

    @pytest.mark.parametrize('n_jobs', [1, 2])
    def test_custom_library_options_reach_select_and_load(
        self,
        tmp_path,
        n_jobs,
    ):
        lib_dir = tmp_path / 'configtoy'
        lib_dir.mkdir()
        (lib_dir / 'configtoy_tools.py').write_text(
            "\n".join([
                'import numpy as np',
                'from optiprofiler.opclasses import Problem',
                '',
                'def configtoy_get_default_options():',
                "    return {'x0_value': 1.0}",
                '',
                'def configtoy_validate_options(options):',
                "    if set(options) != {'x0_value'}:",
                "        raise ValueError('x0_value is the only option')",
                "    return {'x0_value': float(options['x0_value'])}",
                '',
                'def configtoy_select(problem_options, library_options):',
                "    if library_options['x0_value'] != 3.0:",
                "        raise ValueError('selector did not receive options')",
                "    return ['CONFIG_TOY']",
                '',
                'def configtoy_load(problem_name, *, library_options=None):',
                "    x0 = np.array([library_options['x0_value'], 0.0])",
                '    return Problem(lambda x: float(np.dot(x, x)), x0, name=problem_name)',
            ]),
            encoding='utf-8',
        )

        solvers = [simple_solver_1, simple_solver_zero_no_eval]
        feature = Feature('plain')
        problem_options = get_default_problem_options({
            'plibs': ['configtoy'],
            'custom_problem_libs_path': lib_dir,
            'problem_names': ['CONFIG_TOY'],
            'plib_options': {'configtoy': {'x0_value': 3}},
        })
        problem_options[ProblemOption.PLIB_OPTIONS] = (
            _resolve_benchmark_plib_options(problem_options)
        )
        profile_options = get_default_profile_options(
            solvers,
            feature,
            {
                'max_eval_factor': 2,
                'n_jobs': n_jobs,
                'silent': True,
                'solver_names': ['solver1', 'solver2'],
            },
        )

        results = _solve_all_problems(
            solvers,
            'configtoy',
            feature,
            problem_options,
            profile_options,
            False,
            None,
        )

        assert results['plib_options'] == {'x0_value': 3.0}
        assert np.all(results['fun_inits'] == 9.0)

    def test_installed_problem_library_entry_point(self, monkeypatch):
        use_externaltoy_entry_point(monkeypatch)
        solvers = [simple_solver_1, simple_solver_zero_no_eval]
        feature = Feature('plain')
        problem_options = get_default_problem_options({
            'plibs': ['externaltoy'],
            'problem_names': ['EXTERNAL_TOY'],
        })
        profile_options = get_default_profile_options(
            solvers,
            feature,
            {
                'max_eval_factor': 2,
                'n_jobs': 1,
                'silent': True,
                'solver_names': ['solver1', 'solver2'],
            },
        )

        results = _solve_all_problems(
            solvers,
            'externaltoy',
            feature,
            problem_options,
            profile_options,
            False,
            None,
        )

        assert results['plib'] == 'externaltoy'
        assert results['problem_names'] == ['EXTERNAL_TOY']
        assert results['fun_outs'].shape == (1, 2, 1)

    @pytest.mark.parametrize('n_jobs', [1, 2])
    def test_real_entry_point_config_reaches_worker(
        self,
        tmp_path,
        monkeypatch,
        n_jobs,
    ):
        library_name = 'entryconfigtoy'
        module_name = 'optiprofiler_entryconfigtoy'
        (tmp_path / f'{module_name}.py').write_text(
            "\n".join([
                'import numpy as np',
                'from optiprofiler import Problem, ProblemLibraryPlugin',
                '',
                'def get_default_options():',
                "    return {'x0_value': 1.0, 'trace': []}",
                '',
                'def validate_options(options):',
                "    if set(options) != {'x0_value', 'trace'}:",
                "        raise ValueError('wrong option names')",
                "    return {'x0_value': float(options['x0_value']), 'trace': list(options['trace'])}",
                '',
                'def select(problem_options, library_options):',
                "    if library_options != {'x0_value': 4.0, 'trace': []}:",
                "        raise ValueError('selector received the wrong options')",
                "    library_options['trace'].append('selected')",
                "    return ['ENTRY_CONFIG_TOY']",
                '',
                'def load(problem_name, library_options):',
                "    if library_options['trace']:",
                "        raise ValueError('selector mutation leaked into loader')",
                "    x0 = np.array([library_options['x0_value'], 0.0])",
                '    return Problem(lambda x: float(np.dot(x, x)), x0, name=problem_name)',
                '',
                'def get_problem_library():',
                '    return ProblemLibraryPlugin(',
                "        name='entryconfigtoy',",
                '        api_version=1,',
                '        select=select,',
                '        load=load,',
                '        get_default_options=get_default_options,',
                '        validate_options=validate_options,',
                '    )',
            ]),
            encoding='utf-8',
        )
        dist_info = tmp_path / 'optiprofiler_entryconfigtoy-1.0.dist-info'
        dist_info.mkdir()
        (dist_info / 'METADATA').write_text(
            '\n'.join([
                'Metadata-Version: 2.1',
                'Name: optiprofiler-entryconfigtoy',
                'Version: 1.0',
                '',
            ]),
            encoding='utf-8',
        )
        (dist_info / 'entry_points.txt').write_text(
            '\n'.join([
                '[optiprofiler.problem_libraries]',
                f'{library_name} = {module_name}:get_problem_library',
                '',
            ]),
            encoding='utf-8',
        )
        monkeypatch.syspath_prepend(str(tmp_path))

        assert library_name in problem_libraries.list_problem_libraries()
        problem_options = get_default_problem_options({
            'plibs': [library_name],
            'problem_names': ['ENTRY_CONFIG_TOY'],
            'plib_options': {library_name: {'x0_value': 4}},
        })
        problem_options[ProblemOption.PLIB_OPTIONS] = (
            _resolve_benchmark_plib_options(problem_options)
        )
        profile_options = get_default_profile_options(
            [simple_solver_1, simple_solver_zero_no_eval],
            Feature('plain'),
            {
                'max_eval_factor': 2,
                'n_jobs': n_jobs,
                'silent': True,
                'solver_names': ['solver1', 'solver2'],
            },
        )

        results = _solve_all_problems(
            [simple_solver_1, simple_solver_zero_no_eval],
            library_name,
            Feature('plain'),
            problem_options,
            profile_options,
            False,
            None,
        )

        assert results['plib_options'] == {'x0_value': 4.0, 'trace': []}
        assert np.all(results['fun_inits'] == 16.0)

    def test_plugin_availability_is_checked_once(self, monkeypatch):
        checks = []
        plugin = ProblemLibraryPlugin(
            name='externaltoy',
            api_version=1,
            select=externaltoy_select,
            load=externaltoy_load,
            check_available=lambda: checks.append('checked'),
        )
        entry_point = ExternalToyEntryPoint()
        monkeypatch.setattr(entry_point, 'load', lambda: lambda: plugin)
        use_externaltoy_entry_point(monkeypatch, entry_point)
        solvers = [simple_solver_1, simple_solver_zero_no_eval]
        feature = Feature('plain')
        problem_options = get_default_problem_options({
            'plibs': ['externaltoy'],
            'problem_names': ['EXTERNAL_TOY'],
        })
        profile_options = get_default_profile_options(
            solvers,
            feature,
            {
                'max_eval_factor': 2,
                'n_jobs': 1,
                'silent': True,
                'solver_names': ['solver1', 'solver2'],
            },
        )

        results = _solve_all_problems(
            solvers,
            'externaltoy',
            feature,
            problem_options,
            profile_options,
            False,
            None,
        )

        assert results['problem_names'] == ['EXTERNAL_TOY']
        assert checks == ['checked']

    def test_plugin_selection_error_is_not_hidden(self, monkeypatch):
        def failing_select(options, library_options):
            raise RuntimeError('selection failed')

        plugin = ProblemLibraryPlugin(
            name='externaltoy',
            api_version=1,
            select=failing_select,
            load=externaltoy_load,
        )
        entry_point = ExternalToyEntryPoint()
        monkeypatch.setattr(entry_point, 'load', lambda: lambda: plugin)
        use_externaltoy_entry_point(monkeypatch, entry_point)
        solvers = [simple_solver_1, simple_solver_zero_no_eval]
        feature = Feature('plain')
        problem_options = get_default_problem_options({
            'plibs': ['externaltoy'],
        })
        profile_options = get_default_profile_options(
            solvers,
            feature,
            {
                'max_eval_factor': 2,
                'n_jobs': 1,
                'silent': True,
                'solver_names': ['solver1', 'solver2'],
            },
        )

        with pytest.raises(RuntimeError, match='selection failed'):
            _solve_all_problems(
                solvers,
                'externaltoy',
                feature,
                problem_options,
                profile_options,
                False,
                None,
            )

    @pytest.mark.parametrize(
        ('selected_problem_names', 'message'),
        [
            ('EXTERNAL_TOY', 'ordered sequence'),
            (['EXTERNAL_TOY', 1], 'non-string problem name'),
        ],
    )
    def test_plugin_selection_result_is_validated(
        self, monkeypatch, selected_problem_names, message
    ):
        plugin = ProblemLibraryPlugin(
            name='externaltoy',
            api_version=1,
            select=lambda options, library_options: selected_problem_names,
            load=externaltoy_load,
        )
        entry_point = ExternalToyEntryPoint()
        monkeypatch.setattr(entry_point, 'load', lambda: lambda: plugin)
        use_externaltoy_entry_point(monkeypatch, entry_point)
        solvers = [simple_solver_1, simple_solver_zero_no_eval]
        feature = Feature('plain')
        problem_options = get_default_problem_options({
            'plibs': ['externaltoy'],
        })
        profile_options = get_default_profile_options(
            solvers,
            feature,
            {
                'max_eval_factor': 2,
                'n_jobs': 1,
                'silent': True,
                'solver_names': ['solver1', 'solver2'],
            },
        )

        with pytest.raises(TypeError, match=message):
            _solve_all_problems(
                solvers,
                'externaltoy',
                feature,
                problem_options,
                profile_options,
                False,
                None,
            )

    def test_score_only(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, profile_scores, curves = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                plibs=['s2mpj'],
                ptype='u',
                mindim=2,
                maxdim=2,
                max_eval_factor=10,
                benchmark_id='test',
                savepath=tmpdir,
                n_jobs=1,
                silent=True,
                score_only=True,
                problem_names=['ROSENBR'],
            )
            assert isinstance(scores, np.ndarray)

    def test_run_plain(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='noisy',
                n_runs=1,
                run_plain=True,
                **_common_kwargs(tmpdir),
            )
            assert isinstance(scores, np.ndarray)

    def test_sequential_hist_plots(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                plibs=['s2mpj'],
                ptype='u',
                mindim=2,
                maxdim=2,
                max_eval_factor=10,
                benchmark_id='test',
                savepath=tmpdir,
                n_jobs=1,
                silent=True,
                draw_hist_plots='sequential',
                problem_names=['ROSENBR'],
            )
            assert isinstance(scores, np.ndarray)
            self._assert_history_plot_files(tmpdir)

    def test_constrained_hist_plots(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [history_solver_1, history_solver_2],
                feature_name='plain',
                plibs=['s2mpj'],
                ptype='n',
                mindim=2,
                maxdim=2,
                max_eval_factor=10,
                maxnlcon=2,
                maxcon=2,
                benchmark_id='test',
                savepath=tmpdir,
                n_jobs=1,
                silent=True,
                draw_hist_plots='sequential',
                problem_names=['BT1'],
                solver_names=['solver_1', 'solver_2'],
            )
            assert isinstance(scores, np.ndarray)
            self._assert_history_plot_files(tmpdir, 'BT1')

    def test_parallel_hist_plots(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                plibs=['s2mpj'],
                ptype='u',
                mindim=2,
                maxdim=2,
                max_eval_factor=10,
                benchmark_id='test',
                savepath=tmpdir,
                n_jobs=1,
                silent=True,
                draw_hist_plots='parallel',
                problem_names=['ROSENBR'],
            )
            assert isinstance(scores, np.ndarray)
            self._assert_history_plot_files(tmpdir)

    def test_profile_titles_escape_only_for_latex_backend(self):
        assert _format_feature_title('plain_feature', {'text.usetex': False}) == 'Profiles with "plain_feature" feature'
        assert _format_feature_title('plain_feature', {'text.usetex': True}) == r"Profiles with ``plain\_feature'' feature"
        assert _format_problem_feature_title('PROB_1', 'feat_1', {'text.usetex': False}) == 'Solving "PROB_1" with "feat_1" feature'
        assert _format_problem_feature_title('PROB_1', 'feat_1', {'text.usetex': True}) == r"Solving ``PROB\_1'' with ``feat\_1'' feature"


class TestBenchmarkFeatures:

    @pytest.fixture(autouse=True)
    def setup_tmpdir(self):
        self.tmpdir = tempfile.mkdtemp()
        yield
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _run_benchmark(self, feature_name, **extra_kwargs):
        kwargs = _common_kwargs(self.tmpdir)
        kwargs.update(extra_kwargs)
        return benchmark(
            [simple_solver_1, simple_solver_2],
            feature_name=feature_name,
            **kwargs,
        )

    def test_noisy(self):
        scores, _, _ = self._run_benchmark('noisy', n_runs=2)
        assert isinstance(scores, np.ndarray)
        assert scores.shape[0] == 2

    def test_noisy_deterministic(self):
        scores, _, curves = self._run_benchmark('noisy', noise_mode='deterministic')
        assert isinstance(scores, np.ndarray)
        assert scores.shape[0] == 2
        assert len(curves[0]['hist']['perf'][0]) == 2

        scores, _, curves = self._run_benchmark(
            'noisy',
            noise_mode='deterministic',
            solver_isrand=[True, False],
        )
        assert isinstance(scores, np.ndarray)
        assert scores.shape[0] == 2
        assert len(curves[0]['hist']['perf'][0]) == 6

    def test_truncated(self):
        scores, _, _ = self._run_benchmark('truncated')
        assert isinstance(scores, np.ndarray)

    def test_perturbed_x0(self):
        scores, _, _ = self._run_benchmark('perturbed_x0', n_runs=2)
        assert isinstance(scores, np.ndarray)

    def test_quantized(self):
        scores, _, _ = self._run_benchmark('quantized')
        assert isinstance(scores, np.ndarray)

    def test_permuted(self):
        scores, _, _ = self._run_benchmark('permuted', n_runs=2)
        assert isinstance(scores, np.ndarray)

    def test_linearly_transformed(self):
        scores, _, _ = self._run_benchmark('linearly_transformed', n_runs=2)
        assert isinstance(scores, np.ndarray)


class TestBenchmarkErrors:

    def test_no_solvers_no_load(self):
        with pytest.raises((TypeError, ValueError)):
            benchmark(None)

    def test_invalid_feature_name(self):
        with pytest.raises((TypeError, ValueError)):
            benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='invalid_feature',
            )

    def test_single_solver_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(ValueError, match='At least two'):
                benchmark(
                    [simple_solver_1],
                    feature_name='plain',
                    **_common_kwargs(tmpdir),
                )

    def test_solvers_not_callable(self):
        with pytest.raises(TypeError):
            benchmark(['not_callable', 'also_not_callable'])

    def test_copy_failure_logs_reason(self, capsys):
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch('optiprofiler.profiles.shutil.copy2', side_effect=OSError('copy failed for testing')):
                benchmark(
                    [simple_solver_1, simple_solver_2],
                    feature_name='plain',
                    plibs=['s2mpj'],
                    ptype='u',
                    mindim=2,
                    maxdim=2,
                    max_eval_factor=10,
                    benchmark_id='test',
                    savepath=tmpdir,
                    n_jobs=1,
                    silent=False,
                    draw_hist_plots='none',
                    problem_names=['ROSENBR'],
                )
        captured = capsys.readouterr()
        normalized_output = ' '.join(captured.out.split())
        assert 'Failed to copy the script or function that calls `benchmark` function to the log directory.' in normalized_output
        assert 'Error message: copy failed for testing' in normalized_output


class TestBenchmarkOptions:

    def test_custom_solver_names(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                solver_names=['MySolver1', 'MySolver2'],
                **_common_kwargs(tmpdir),
            )
            assert isinstance(scores, np.ndarray)

    def test_custom_max_tol_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                max_tol_order=3,
                **_common_kwargs(tmpdir),
            )
            assert isinstance(scores, np.ndarray)

    def test_custom_seed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                seed=42,
                **_common_kwargs(tmpdir),
            )
            assert isinstance(scores, np.ndarray)

    def test_project_x0(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                project_x0=True,
                **_common_kwargs(tmpdir),
            )
            assert isinstance(scores, np.ndarray)

    def test_custom_merit_fun(self):
        def my_merit(f, v, v0):
            return f
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                merit_fun=my_merit,
                **_common_kwargs(tmpdir),
            )
            assert isinstance(scores, np.ndarray)

    def test_normalized_scores_false(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                normalized_scores=False,
                **_common_kwargs(tmpdir),
            )
            assert isinstance(scores, np.ndarray)

    def test_semilogx_false(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scores, _, _ = benchmark(
                [simple_solver_1, simple_solver_2],
                feature_name='plain',
                semilogx=False,
                **_common_kwargs(tmpdir),
            )
            assert isinstance(scores, np.ndarray)
