"""Tests for the utils module."""
import logging
import multiprocessing as mp

import numpy as np
import pytest

from optiprofiler.utils import (
    FeatureName,
    FeatureOption,
    ProblemError,
    ProblemOption,
    ProfileOption,
    get_logger,
    setup_main_process_logging,
    setup_worker_logging,
    show_versions,
    _get_deps_info,
    _get_sys_info,
)


class TestFeatureName:

    def test_enumeration_values(self):
        assert FeatureName.PLAIN == 'plain'
        assert FeatureName.PERTURBED_X0 == 'perturbed_x0'
        assert FeatureName.NOISY == 'noisy'
        assert FeatureName.TRUNCATED == 'truncated'
        assert FeatureName.PERMUTED == 'permuted'
        assert FeatureName.LINEARLY_TRANSFORMED == 'linearly_transformed'
        assert FeatureName.RANDOM_NAN == 'random_nan'
        assert FeatureName.UNRELAXABLE_CONSTRAINTS == 'unrelaxable_constraints'
        assert FeatureName.NONQUANTIFIABLE_CONSTRAINTS == 'nonquantifiable_constraints'
        assert FeatureName.QUANTIZED == 'quantized'
        assert FeatureName.CUSTOM == 'custom'

    def test_constructor(self):
        assert FeatureName('plain') == FeatureName.PLAIN
        assert FeatureName('noisy') == FeatureName.NOISY
        with pytest.raises(ValueError):
            FeatureName('invalid')

    def test_membership(self):
        assert 'plain' in FeatureName.__members__.values()
        assert 'noisy' in FeatureName.__members__.values()
        assert 'invalid' not in FeatureName.__members__.values()


class TestProfileOption:

    def test_enumeration_values(self):
        assert ProfileOption.N_JOBS == 'n_jobs'
        assert ProfileOption.SEED == 'seed'
        assert ProfileOption.BENCHMARK_ID == 'benchmark_id'
        assert ProfileOption.SOLVER_NAMES == 'solver_names'
        assert ProfileOption.SOLVER_ISRAND == 'solver_isrand'
        assert ProfileOption.FEATURE_STAMP == 'feature_stamp'
        assert ProfileOption.MAX_TOL_ORDER == 'max_tol_order'
        assert ProfileOption.MAX_EVAL_FACTOR == 'max_eval_factor'
        assert ProfileOption.MERIT_FUN == 'merit_fun'
        assert ProfileOption.PROJECT_X0 == 'project_x0'
        assert ProfileOption.RUN_PLAIN == 'run_plain'
        assert ProfileOption.SCORE_ONLY == 'score_only'
        assert ProfileOption.DRAW_HIST_PLOTS == 'draw_hist_plots'
        assert ProfileOption.SILENT == 'silent'
        assert ProfileOption.SEMILOGX == 'semilogx'
        assert ProfileOption.LOAD == 'load'

    def test_constructor(self):
        assert ProfileOption('n_jobs') == ProfileOption.N_JOBS
        with pytest.raises(ValueError):
            ProfileOption('invalid')


class TestProblemOption:

    def test_enumeration_values(self):
        assert ProblemOption.PLIBS == 'plibs'
        assert ProblemOption.PTYPE == 'ptype'
        assert ProblemOption.MINDIM == 'mindim'
        assert ProblemOption.MAXDIM == 'maxdim'
        assert ProblemOption.MINB == 'minb'
        assert ProblemOption.MAXB == 'maxb'
        assert ProblemOption.MINLCON == 'minlcon'
        assert ProblemOption.MAXLCON == 'maxlcon'
        assert ProblemOption.MINNLCON == 'minnlcon'
        assert ProblemOption.MAXNLCON == 'maxnlcon'
        assert ProblemOption.MINCON == 'mincon'
        assert ProblemOption.MAXCON == 'maxcon'
        assert ProblemOption.EXCLUDELIST == 'excludelist'
        assert ProblemOption.PROBLEM_NAMES == 'problem_names'
        assert ProblemOption.CUSTOM_PROBLEM_LIBS_PATH == 'custom_problem_libs_path'

    def test_constructor(self):
        assert ProblemOption('plibs') == ProblemOption.PLIBS
        with pytest.raises(ValueError):
            ProblemOption('invalid')


class TestFeatureOption:

    def test_enumeration_values(self):
        assert FeatureOption.N_RUNS == 'n_runs'
        assert FeatureOption.DISTRIBUTION == 'distribution'
        assert FeatureOption.PERTURBATION_LEVEL == 'perturbation_level'
        assert FeatureOption.NOISE_LEVEL == 'noise_level'
        assert FeatureOption.NOISE_TYPE == 'noise_type'
        assert FeatureOption.SIGNIFICANT_DIGITS == 'significant_digits'
        assert FeatureOption.PERTURBED_TRAILING_DIGITS == 'perturbed_trailing_digits'
        assert FeatureOption.ROTATED == 'rotated'
        assert FeatureOption.CONDITION_FACTOR == 'condition_factor'
        assert FeatureOption.NAN_RATE == 'nan_rate'
        assert FeatureOption.MESH_SIZE == 'mesh_size'
        assert FeatureOption.MESH_TYPE == 'mesh_type'
        assert FeatureOption.GROUND_TRUTH == 'ground_truth'
        assert FeatureOption.MOD_FUN == 'mod_fun'
        assert FeatureOption.MOD_X0 == 'mod_x0'

    def test_constructor(self):
        assert FeatureOption('n_runs') == FeatureOption.N_RUNS
        with pytest.raises(ValueError):
            FeatureOption('invalid')


class TestProblemError:

    def test_basic(self):
        err = ProblemError('test error')
        assert str(err) == 'test error'
        assert err.message == 'test error'

    def test_raise(self):
        with pytest.raises(ProblemError, match='something went wrong'):
            raise ProblemError('something went wrong')

    def test_inheritance(self):
        assert issubclass(ProblemError, Exception)


class TestShowVersions:

    def test_runs_without_error(self, capsys):
        show_versions()
        captured = capsys.readouterr()
        assert 'System settings' in captured.out
        assert 'Python dependencies' in captured.out
        assert 'python' in captured.out
        assert 'numpy' in captured.out

    def test_get_sys_info(self):
        info = _get_sys_info()
        assert 'python' in info
        assert 'executable' in info
        assert 'machine' in info

    def test_get_deps_info(self):
        info = _get_deps_info()
        assert 'numpy' in info
        assert 'matplotlib' in info
        assert 'scipy' in info
        assert info['numpy'] is not None


class TestLogging:

    def test_get_logger(self):
        logger = get_logger('test_logger')
        assert isinstance(logger, logging.Logger)
        assert logger.name == 'test_logger'
        assert logger.level == logging.INFO
        assert logger.propagate is True

    def test_get_logger_custom_level(self):
        logger = get_logger('test_debug', level=logging.DEBUG)
        assert logger.level == logging.DEBUG

    def test_setup_main_process_logging(self):
        log_queue, listener = setup_main_process_logging()
        assert isinstance(log_queue, mp.queues.Queue)
        listener.stop()

    def test_setup_worker_logging_none(self):
        setup_worker_logging(None)

    def test_setup_worker_logging(self):
        log_queue = mp.Queue(-1)
        setup_worker_logging(log_queue, level=logging.DEBUG)
        root = logging.getLogger()
        assert root.level == logging.DEBUG
