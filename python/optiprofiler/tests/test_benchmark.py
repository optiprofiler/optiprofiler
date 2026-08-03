"""Tests for the benchmark function and solve helpers in profiles module."""
import os
import shutil
import tempfile
from pathlib import Path
from unittest.mock import patch

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pytest

from optiprofiler import benchmark
from optiprofiler.opclasses import Feature, FeaturedProblem, Problem
from optiprofiler.profile_utils import get_default_problem_options, get_default_profile_options
from optiprofiler.profiles import _format_feature_title, _format_problem_feature_title, _solve_all_problems


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

    def test_custom_problem_library_with_direct_path(self, tmp_path):
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
                'n_jobs': 1,
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


def _load_solver_1(fun, x0):
    """Records at least one history entry so loaded history plots have data."""
    fun(x0)
    return x0


def _load_solver_2(fun, x0):
    """Records two history entries so loaded history plots have data."""
    fun(x0)
    x = x0.copy()
    x[0] += 0.1
    fun(x)
    return x


def _load_solver_3(fun, x0):
    """A third distinct solver for the multi-solver load fixtures."""
    fun(x0)
    x = x0.copy()
    x[0] -= 0.2
    fun(x)
    return x


@pytest.fixture(scope='module')
def saved_two_solver_experiment(tmp_path_factory):
    """A small saved experiment with two solvers, reused by the load tests."""
    tmpdir = str(tmp_path_factory.mktemp('load_two'))
    benchmark(
        [_load_solver_1, _load_solver_2],
        feature_name='plain',
        **_common_kwargs(tmpdir, benchmark_id='loadtest'),
    )
    return tmpdir


@pytest.fixture(scope='module')
def saved_three_solver_experiment(tmp_path_factory):
    """A small saved experiment with three solvers, reused by the load tests."""
    tmpdir = str(tmp_path_factory.mktemp('load_three'))
    benchmark(
        [_load_solver_1, _load_solver_2, _load_solver_3],
        feature_name='plain',
        **_common_kwargs(tmpdir, benchmark_id='loadtest3'),
    )
    return tmpdir


class TestBenchmarkLoad:
    """End-to-end checks of the load path: deferred solver-dependent option
    validation and the post-load two-solver decision for log-ratio profiles."""

    @staticmethod
    def _load_and_diff(tmpdir, benchmark_id, **load_kwargs):
        """Run a load call and return (scores, stamp dirs created by this call)."""
        bench_root = Path(tmpdir) / benchmark_id
        before = {p for p in bench_root.iterdir() if p.is_dir()}
        scores, _, _ = benchmark(
            None,
            load='latest',
            benchmark_id=benchmark_id,
            savepath=tmpdir,
            silent=True,
            n_jobs=1,
            **load_kwargs,
        )
        new_dirs = [p for p in bench_root.iterdir() if p.is_dir()] 
        new_dirs = [p for p in new_dirs if p not in before]
        return scores, new_dirs

    def test_load_two_solvers_generates_log_ratio(self, saved_two_solver_experiment, monkeypatch):
        tmpdir = saved_two_solver_experiment
        monkeypatch.chdir(tmpdir)
        scores, new_dirs = self._load_and_diff(
            tmpdir, 'loadtest',
            summarize_log_ratio_profiles=True, max_tol_order=2,
        )
        assert scores.shape[0] == 2
        assert new_dirs
        assert any(list(d.glob('**/log-ratio_hist_*.pdf')) for d in new_dirs), \
            'log-ratio profile PDFs missing after a two-solver load'

    def test_load_three_solvers_warns_and_skips_log_ratio(self, saved_three_solver_experiment, monkeypatch, capsys):
        tmpdir = saved_three_solver_experiment
        monkeypatch.chdir(tmpdir)
        scores, new_dirs = self._load_and_diff(
            tmpdir, 'loadtest3',
            summarize_log_ratio_profiles=True, max_tol_order=2,
        )
        assert scores.shape[0] == 3
        assert 'exactly two solvers' in capsys.readouterr().out
        for d in new_dirs:
            for pdf in d.glob('**/log-ratio_hist_*.pdf'):
                pytest.fail(f'log-ratio profile {pdf} generated for a three-solver load')

    def test_load_select_two_of_three_generates_log_ratio(self, saved_three_solver_experiment, monkeypatch):
        tmpdir = saved_three_solver_experiment
        monkeypatch.chdir(tmpdir)
        scores, new_dirs = self._load_and_diff(
            tmpdir, 'loadtest3',
            solvers_to_load=[0, 2], summarize_log_ratio_profiles=True, max_tol_order=2,
        )
        assert scores.shape[0] == 2
        assert any(list(d.glob('**/log-ratio_hist_*.pdf')) for d in new_dirs), \
            'log-ratio profile PDFs missing after selecting two of three solvers'

    def test_load_with_solver_names_renames(self, saved_two_solver_experiment, monkeypatch):
        tmpdir = saved_two_solver_experiment
        monkeypatch.chdir(tmpdir)
        scores, _, _ = benchmark(
            None,
            load='latest',
            benchmark_id='loadtest',
            savepath=tmpdir,
            solver_names=['first', 'second'],
            max_tol_order=1,
            silent=True,
            n_jobs=1,
        )
        assert scores.shape[0] == 2

    def test_load_with_solver_isrand_accepted(self, saved_two_solver_experiment, monkeypatch):
        tmpdir = saved_two_solver_experiment
        monkeypatch.chdir(tmpdir)
        scores, _, _ = benchmark(
            None,
            load='latest',
            benchmark_id='loadtest',
            savepath=tmpdir,
            solver_isrand=[True, False],
            max_tol_order=1,
            silent=True,
            n_jobs=1,
        )
        assert scores.shape[0] == 2

    def test_load_with_wrong_solver_names_length_raises(self, saved_two_solver_experiment, monkeypatch):
        tmpdir = saved_two_solver_experiment
        monkeypatch.chdir(tmpdir)
        with pytest.raises(ValueError, match='loaded solvers'):
            benchmark(
                None,
                load='latest',
                benchmark_id='loadtest',
                savepath=tmpdir,
                solver_names=['only_one'],
                silent=True,
            )

    def test_load_with_single_solvers_to_load_raises(self, saved_two_solver_experiment, monkeypatch):
        tmpdir = saved_two_solver_experiment
        monkeypatch.chdir(tmpdir)
        with pytest.raises(ValueError, match='at least two indices'):
            benchmark(
                None,
                load='latest',
                benchmark_id='loadtest',
                savepath=tmpdir,
                solvers_to_load=[0],
                silent=True,
            )


class TestBenchmarkLoadStamp:
    """The stamp of a reloaded experiment must describe the loaded subset."""

    def test_load_without_options_keeps_saved_selection_in_stamp(self, saved_two_solver_experiment, monkeypatch):
        tmpdir = saved_two_solver_experiment
        monkeypatch.chdir(tmpdir)
        bench_root = Path(tmpdir) / 'loadtest'
        before = {p for p in bench_root.iterdir() if p.is_dir()}
        benchmark(
            None,
            load='latest',
            benchmark_id='loadtest',
            savepath=tmpdir,
            max_tol_order=1,
            silent=True,
            n_jobs=1,
        )
        new_dirs = [p for p in bench_root.iterdir() if p.is_dir() and p not in before]
        assert new_dirs
        # The fixture experiment was saved with ptype 'u', mindim 2, maxdim 2.
        assert any('_u_2_2_' in d.name for d in new_dirs), \
            f'stamp does not describe the loaded subset: {[d.name for d in new_dirs]}'
        for d in new_dirs:
            assert '_u_1_2_' not in d.name, f'default problem selection leaked into stamp {d.name}'


class TestBenchmarkLoadDrawHistPrecedence:
    """Precedence of draw_hist_plots on the load path: score_only is strongest,
    an explicit user choice is honored, and only an omitted option defaults to
    sequential drawing after loading."""

    @staticmethod
    def _history_pdfs(dirs):
        return [pdf for d in dirs for pdf in d.glob('history_plots/**/*.pdf')]

    def test_explicit_none_is_honored_on_load(self, saved_two_solver_experiment, monkeypatch, capsys):
        tmpdir = saved_two_solver_experiment
        monkeypatch.chdir(tmpdir)
        scores, new_dirs = TestBenchmarkLoad._load_and_diff(
            tmpdir, 'loadtest',
            draw_hist_plots='none', max_tol_order=1,
        )
        assert scores.shape[0] == 2
        assert not self._history_pdfs(new_dirs), 'history plots were drawn despite an explicit draw_hist_plots="none"'
        assert 'is set to "sequential"' not in ' '.join(capsys.readouterr().out.split())

    def test_omitted_defaults_to_sequential_on_load(self, saved_two_solver_experiment, monkeypatch, capsys):
        tmpdir = saved_two_solver_experiment
        monkeypatch.chdir(tmpdir)
        scores, new_dirs = TestBenchmarkLoad._load_and_diff(
            tmpdir, 'loadtest',
            max_tol_order=1,
        )
        assert scores.shape[0] == 2
        assert self._history_pdfs(new_dirs), 'history plots missing although draw_hist_plots was omitted on load'
        assert 'is set to "sequential"' in ' '.join(capsys.readouterr().out.split())

    def test_explicit_parallel_still_draws_on_load(self, saved_two_solver_experiment, monkeypatch):
        tmpdir = saved_two_solver_experiment
        monkeypatch.chdir(tmpdir)
        scores, new_dirs = TestBenchmarkLoad._load_and_diff(
            tmpdir, 'loadtest',
            draw_hist_plots='parallel', max_tol_order=1,
        )
        assert scores.shape[0] == 2
        assert self._history_pdfs(new_dirs), 'history plots missing for an explicit draw_hist_plots="parallel" on load'

    def test_score_only_with_explicit_drawing_warns_and_suppresses(self, saved_two_solver_experiment, monkeypatch, capsys):
        tmpdir = saved_two_solver_experiment
        monkeypatch.chdir(tmpdir)
        bench_root = Path(tmpdir) / 'loadtest'
        before = {p for p in bench_root.iterdir() if p.is_dir()}
        scores, _, _ = benchmark(
            None,
            load='latest',
            benchmark_id='loadtest',
            savepath=tmpdir,
            score_only=True,
            draw_hist_plots='sequential',
            silent=True,
            n_jobs=1,
        )
        assert scores.shape[0] == 2
        out = ' '.join(capsys.readouterr().out.split())
        assert 'is ignored and no history plots will be drawn' in out
        new_dirs = [p for p in bench_root.iterdir() if p.is_dir() and p not in before]
        assert not self._history_pdfs(new_dirs), 'history plots were drawn despite score_only=True'


class TestBenchmarkRendererFallback:
    """Without a latex executable, plot text must fall back to plain rendering:
    identity text such as solver names must not carry LaTeX escapes."""

    def test_no_latex_pdf_has_unescaped_solver_names(self, capsys):
        from pypdf import PdfReader
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch('optiprofiler.plotting.shutil.which', return_value=None):
                benchmark(
                    [simple_solver_1, simple_solver_2],
                    feature_name='plain',
                    plibs=['s2mpj'],
                    ptype='u',
                    mindim=2,
                    maxdim=2,
                    max_eval_factor=10,
                    max_tol_order=1,
                    benchmark_id='render',
                    savepath=tmpdir,
                    n_jobs=1,
                    silent=False,
                    draw_hist_plots='none',
                    problem_names=['ROSENBR'],
                )
            assert 'Rendering plot text without LaTeX' in ' '.join(capsys.readouterr().out.split())
            summary_pdfs = list(Path(tmpdir).glob('render/**/summary_*.pdf'))
            assert summary_pdfs
            text = ''.join(page.extract_text() or '' for page in PdfReader(str(summary_pdfs[0])).pages)
            assert 'simple_solver_1' in text.replace('\n', '')
            assert 'simple\\_solver' not in text, 'LaTeX escapes leaked into a non-LaTeX PDF'
