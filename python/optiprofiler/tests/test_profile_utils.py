"""Tests for the profile_utils module."""
import os
import tempfile

import numpy as np
import pytest

from optiprofiler.opclasses import Feature
from optiprofiler.profile_utils import (
    _default_merit,
    _default_score_fun,
    _default_score_weight_fun,
    _get_default_feature_stamp,
    add_to_readme,
    check_validity_problem_options,
    check_validity_profile_options,
    compute_merit_values,
    create_stamp,
    get_default_problem_options,
    get_default_profile_options,
    init_readme,
    integrate,
    merge_pdfs_with_pypdf,
)
from optiprofiler.utils import FeatureOption, ProblemOption, ProfileOption


class TestCheckValidityProblemOptions:

    def test_empty_options(self):
        opts = check_validity_problem_options({})
        assert isinstance(opts, dict)

    def test_ptype_valid(self):
        opts = check_validity_problem_options({ProblemOption.PTYPE: 'u'})
        assert opts[ProblemOption.PTYPE] == 'u'
        opts = check_validity_problem_options({ProblemOption.PTYPE: 'ubln'})
        assert opts[ProblemOption.PTYPE] == 'ubln'

    def test_ptype_invalid_type(self):
        with pytest.raises(TypeError):
            check_validity_problem_options({ProblemOption.PTYPE: 123})

    def test_ptype_empty(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.PTYPE: ''})

    def test_ptype_invalid_char(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.PTYPE: 'x'})

    def test_mindim_valid(self):
        opts = check_validity_problem_options({ProblemOption.MINDIM: 5})
        assert opts[ProblemOption.MINDIM] == 5

    def test_mindim_float_integer(self):
        opts = check_validity_problem_options({ProblemOption.MINDIM: 5.0})
        assert opts[ProblemOption.MINDIM] == 5
        assert isinstance(opts[ProblemOption.MINDIM], int)

    def test_mindim_numpy_integer(self):
        opts = check_validity_problem_options({ProblemOption.MINDIM: np.int64(5)})
        assert opts[ProblemOption.MINDIM] == 5
        assert isinstance(opts[ProblemOption.MINDIM], int)

    def test_mindim_invalid_type(self):
        with pytest.raises(TypeError):
            check_validity_problem_options({ProblemOption.MINDIM: 'five'})

    def test_mindim_too_small(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.MINDIM: 0})

    def test_maxdim_valid(self):
        opts = check_validity_problem_options({ProblemOption.MAXDIM: 10})
        assert opts[ProblemOption.MAXDIM] == 10

    def test_maxdim_inf(self):
        opts = check_validity_problem_options({ProblemOption.MAXDIM: np.inf})
        assert np.isinf(opts[ProblemOption.MAXDIM])

    def test_maxdim_invalid_type(self):
        with pytest.raises(TypeError):
            check_validity_problem_options({ProblemOption.MAXDIM: 'ten'})

    def test_maxdim_too_small(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.MAXDIM: 0})

    def test_mindim_greater_than_maxdim(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.MINDIM: 10, ProblemOption.MAXDIM: 5})

    def test_minb_valid(self):
        opts = check_validity_problem_options({ProblemOption.MINB: 0})
        assert opts[ProblemOption.MINB] == 0

    def test_minb_negative(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.MINB: -1})

    def test_maxb_valid(self):
        opts = check_validity_problem_options({ProblemOption.MAXB: np.inf})
        assert np.isinf(opts[ProblemOption.MAXB])

    def test_minb_greater_than_maxb(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.MINB: 10, ProblemOption.MAXB: 5})

    def test_minlcon_valid(self):
        opts = check_validity_problem_options({ProblemOption.MINLCON: 0})
        assert opts[ProblemOption.MINLCON] == 0

    def test_minlcon_negative(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.MINLCON: -1})

    def test_minlcon_greater_than_maxlcon(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.MINLCON: 10, ProblemOption.MAXLCON: 5})

    def test_minnlcon_valid(self):
        opts = check_validity_problem_options({ProblemOption.MINNLCON: 0})
        assert opts[ProblemOption.MINNLCON] == 0

    def test_minnlcon_negative(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.MINNLCON: -1})

    def test_minnlcon_greater_than_maxnlcon(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.MINNLCON: 10, ProblemOption.MAXNLCON: 5})

    def test_mincon_valid(self):
        opts = check_validity_problem_options({ProblemOption.MINCON: 0})
        assert opts[ProblemOption.MINCON] == 0

    def test_mincon_negative(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.MINCON: -1})

    def test_mincon_greater_than_maxcon(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.MINCON: 10, ProblemOption.MAXCON: 5})

    def test_excludelist_string(self):
        opts = check_validity_problem_options({ProblemOption.EXCLUDELIST: 'PROB1'})
        assert opts[ProblemOption.EXCLUDELIST] == ['PROB1']

    def test_excludelist_list(self):
        opts = check_validity_problem_options({ProblemOption.EXCLUDELIST: ['PROB1', 'PROB2']})
        assert opts[ProblemOption.EXCLUDELIST] == ['PROB1', 'PROB2']

    def test_excludelist_invalid_type(self):
        with pytest.raises(TypeError):
            check_validity_problem_options({ProblemOption.EXCLUDELIST: [1, 2]})

    def test_problem_names_string(self):
        opts = check_validity_problem_options({ProblemOption.PROBLEM_NAMES: 'PROB1'})
        assert opts[ProblemOption.PROBLEM_NAMES] == ['PROB1']

    def test_problem_names_invalid_type(self):
        with pytest.raises(TypeError):
            check_validity_problem_options({ProblemOption.PROBLEM_NAMES: [1]})

    def test_plibs_empty(self):
        with pytest.raises(ValueError):
            check_validity_problem_options({ProblemOption.PLIBS: []})

    def test_plibs_invalid_type(self):
        with pytest.raises(TypeError):
            check_validity_problem_options({ProblemOption.PLIBS: [123]})


class TestCheckValidityProfileOptions:

    @staticmethod
    def _dummy_solver(problem):
        return problem.x0

    def test_empty_options(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {})
        assert isinstance(opts, dict)

    def test_n_jobs_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.N_JOBS: 4})
        assert opts[ProfileOption.N_JOBS] == 4

    def test_n_jobs_float_integer(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.N_JOBS: 4.0})
        assert opts[ProfileOption.N_JOBS] == 4
        assert isinstance(opts[ProfileOption.N_JOBS], int)

    def test_n_jobs_invalid_type(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.N_JOBS: 'four'})

    def test_n_jobs_too_small(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.N_JOBS: 0})
        assert opts[ProfileOption.N_JOBS] == 1

    def test_seed_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.SEED: 42})
        assert opts[ProfileOption.SEED] == 42

    def test_seed_invalid_type(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.SEED: 'abc'})

    def test_seed_negative(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.SEED: -1})

    def test_benchmark_id_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.BENCHMARK_ID: 'test_run'})
        assert opts[ProfileOption.BENCHMARK_ID] == 'test_run'

    def test_benchmark_id_invalid_type(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.BENCHMARK_ID: 123})

    def test_benchmark_id_empty(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.BENCHMARK_ID: ''})

    def test_benchmark_id_invalid_chars(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.BENCHMARK_ID: 'test run!'})

    def test_solver_names_valid(self):
        solvers = [self._dummy_solver, self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.SOLVER_NAMES: ['s1', 's2']})
        assert opts[ProfileOption.SOLVER_NAMES] == ['s1', 's2']

    def test_solver_names_wrong_length(self):
        solvers = [self._dummy_solver, self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.SOLVER_NAMES: ['s1']})

    def test_solver_names_invalid_type(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.SOLVER_NAMES: 'not_a_list'})

    def test_solver_isrand_valid(self):
        solvers = [self._dummy_solver, self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.SOLVER_ISRAND: [True, False]})
        assert opts[ProfileOption.SOLVER_ISRAND] == [True, False]

    def test_solver_isrand_wrong_length(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.SOLVER_ISRAND: [True, False]})

    def test_solver_isrand_invalid_type(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.SOLVER_ISRAND: [1]})

    def test_feature_stamp_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.FEATURE_STAMP: 'plain'})
        assert opts[ProfileOption.FEATURE_STAMP] == 'plain'

    def test_feature_stamp_empty(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.FEATURE_STAMP: ''})

    def test_feature_stamp_invalid_chars(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.FEATURE_STAMP: 'a b c'})

    def test_errorbar_type_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.ERRORBAR_TYPE: 'minmax'})
        assert opts[ProfileOption.ERRORBAR_TYPE] == 'minmax'
        opts = check_validity_profile_options(solvers, {ProfileOption.ERRORBAR_TYPE: 'meanstd'})
        assert opts[ProfileOption.ERRORBAR_TYPE] == 'meanstd'

    def test_errorbar_type_invalid(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.ERRORBAR_TYPE: 'invalid'})

    def test_hist_aggregation_valid(self):
        solvers = [self._dummy_solver]
        for mode in ['min', 'mean', 'max']:
            opts = check_validity_profile_options(solvers, {ProfileOption.HIST_AGGREGATION: mode})
            assert opts[ProfileOption.HIST_AGGREGATION] == mode

    def test_hist_aggregation_invalid(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.HIST_AGGREGATION: 'invalid'})

    def test_max_tol_order_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.MAX_TOL_ORDER: 10})
        assert opts[ProfileOption.MAX_TOL_ORDER] == 10

    def test_max_tol_order_out_of_range(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.MAX_TOL_ORDER: 17})

    def test_max_eval_factor_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.MAX_EVAL_FACTOR: 500.0})
        assert opts[ProfileOption.MAX_EVAL_FACTOR] == 500.0

    def test_max_eval_factor_invalid_type(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.MAX_EVAL_FACTOR: 'big'})

    def test_max_eval_factor_nonpositive(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.MAX_EVAL_FACTOR: 0})

    def test_merit_fun_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.MERIT_FUN: lambda f, v, v0: f})
        assert callable(opts[ProfileOption.MERIT_FUN])

    def test_merit_fun_invalid(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.MERIT_FUN: 'not_callable'})

    def test_project_x0_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.PROJECT_X0: True})
        assert opts[ProfileOption.PROJECT_X0] is True

    def test_project_x0_invalid(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.PROJECT_X0: 1})

    def test_draw_hist_plots_valid(self):
        solvers = [self._dummy_solver]
        for mode in ['none', 'parallel', 'sequential']:
            opts = check_validity_profile_options(solvers, {ProfileOption.DRAW_HIST_PLOTS: mode})
            assert opts[ProfileOption.DRAW_HIST_PLOTS] == mode

    def test_draw_hist_plots_invalid(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.DRAW_HIST_PLOTS: 'invalid'})

    def test_silent_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.SILENT: True})
        assert opts[ProfileOption.SILENT] is True

    def test_silent_invalid(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.SILENT: 1})

    def test_semilogx_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.SEMILOGX: False})
        assert opts[ProfileOption.SEMILOGX] is False

    def test_semilogx_invalid(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.SEMILOGX: 1})

    def test_score_weight_fun_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.SCORE_WEIGHT_FUN: lambda x: 1})
        assert callable(opts[ProfileOption.SCORE_WEIGHT_FUN])

    def test_score_weight_fun_invalid(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.SCORE_WEIGHT_FUN: 'not_callable'})

    def test_load_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.LOAD: 'latest'})
        assert opts[ProfileOption.LOAD] == 'latest'

    def test_load_none(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.LOAD: None})
        assert opts[ProfileOption.LOAD] is None

    def test_load_empty_string(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.LOAD: ''})
        assert opts[ProfileOption.LOAD] is None

    def test_load_invalid_type(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.LOAD: 123})

    def test_line_colors_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.LINE_COLORS: ['red', 'blue']})
        assert opts[ProfileOption.LINE_COLORS] == ['red', 'blue']

    def test_line_colors_single(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.LINE_COLORS: 'red'})
        assert opts[ProfileOption.LINE_COLORS] == ['red']

    def test_line_colors_invalid(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.LINE_COLORS: ['not_a_color_xyz']})

    def test_line_widths_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.LINE_WIDTHS: [1.5, 2.0]})
        assert opts[ProfileOption.LINE_WIDTHS] == [1.5, 2.0]

    def test_line_widths_nonpositive(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.LINE_WIDTHS: [0]})

    def test_line_widths_invalid_type(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.LINE_WIDTHS: ['a']})

    def test_bar_colors_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.BAR_COLORS: ['red']})
        assert opts[ProfileOption.BAR_COLORS] == ['red']

    def test_xlabel_ylabel_valid(self):
        solvers = [self._dummy_solver]
        opts = check_validity_profile_options(solvers, {
            ProfileOption.XLABEL_PERFORMANCE_PROFILE: 'X',
            ProfileOption.YLABEL_PERFORMANCE_PROFILE: 'Y',
            ProfileOption.XLABEL_DATA_PROFILE: 'X',
            ProfileOption.YLABEL_DATA_PROFILE: 'Y',
            ProfileOption.XLABEL_LOG_RATIO_PROFILE: 'X',
            ProfileOption.YLABEL_LOG_RATIO_PROFILE: 'Y',
        })
        assert opts[ProfileOption.XLABEL_PERFORMANCE_PROFILE] == 'X'

    def test_xlabel_invalid_type(self):
        solvers = [self._dummy_solver]
        with pytest.raises(TypeError):
            check_validity_profile_options(solvers, {ProfileOption.XLABEL_PERFORMANCE_PROFILE: 123})

    def test_solver_verbose_valid(self):
        solvers = [self._dummy_solver]
        for v in [0, 1, 2]:
            opts = check_validity_profile_options(solvers, {ProfileOption.SOLVER_VERBOSE: v})
            assert opts[ProfileOption.SOLVER_VERBOSE] == v

    def test_solver_verbose_out_of_range(self):
        solvers = [self._dummy_solver]
        with pytest.raises(ValueError):
            check_validity_profile_options(solvers, {ProfileOption.SOLVER_VERBOSE: 3})

    def test_summarize_log_ratio_not_two_solvers(self):
        solvers = [self._dummy_solver, self._dummy_solver, self._dummy_solver]
        opts = check_validity_profile_options(solvers, {ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES: True})
        assert opts[ProfileOption.SUMMARIZE_LOG_RATIO_PROFILES] is False


class TestGetDefaultProblemOptions:

    def test_defaults(self):
        opts = get_default_problem_options({})
        assert opts[ProblemOption.PLIBS.value] == ['s2mpj']
        assert opts[ProblemOption.PTYPE.value] == 'u'
        assert opts[ProblemOption.MINDIM.value] == 1
        assert opts[ProblemOption.MAXDIM.value] == 2
        assert opts[ProblemOption.MINB.value] == 0
        assert opts[ProblemOption.MINLCON.value] == 0
        assert opts[ProblemOption.MINNLCON.value] == 0

    def test_preserves_existing(self):
        opts = get_default_problem_options({ProblemOption.MINDIM.value: 5})
        assert opts[ProblemOption.MINDIM.value] == 5
        assert opts[ProblemOption.MAXDIM.value] == 6


class TestGetDefaultProfileOptions:

    @staticmethod
    def _dummy_solver(problem):
        return problem.x0

    def test_defaults(self):
        solvers = [self._dummy_solver]
        feature = Feature('plain')
        opts = get_default_profile_options(solvers, feature, {})
        assert opts[ProfileOption.SEED.value] == 0
        assert opts[ProfileOption.BENCHMARK_ID.value] == 'out'
        assert opts[ProfileOption.MAX_TOL_ORDER.value] == 10
        assert opts[ProfileOption.MAX_EVAL_FACTOR.value] == 500
        assert opts[ProfileOption.PROJECT_X0.value] is False
        assert opts[ProfileOption.RUN_PLAIN.value] is False
        assert opts[ProfileOption.SCORE_ONLY.value] is False
        assert opts[ProfileOption.SEMILOGX.value] is True
        assert callable(opts[ProfileOption.MERIT_FUN.value])

    def test_solver_names_auto(self):
        def my_solver(problem):
            return problem.x0
        solvers = [my_solver]
        feature = Feature('plain')
        opts = get_default_profile_options(solvers, feature, {})
        assert opts[ProfileOption.SOLVER_NAMES.value] == ['my_solver']


class TestDefaultMerit:

    def test_feasible(self):
        assert _default_merit(1.0, 0.0, 0.0) == 1.0

    def test_nan_fun(self):
        assert np.isinf(_default_merit(np.nan, 0.0, 0.0))

    def test_nan_maxcv(self):
        assert np.isinf(_default_merit(1.0, np.nan, 0.0))

    def test_small_violation(self):
        result = _default_merit(1.0, 0.005, 0.0)
        assert result > 1.0

    def test_large_violation(self):
        result = _default_merit(1.0, 1.0, 0.0)
        assert np.isinf(result)

    def test_with_nonzero_maxcv_init(self):
        result = _default_merit(1.0, 0.0, 1.0)
        assert result == 1.0


class TestDefaultScoreFunctions:

    def test_default_score_weight_fun(self):
        assert _default_score_weight_fun(0) == 1
        assert _default_score_weight_fun(100) == 1

    def test_default_score_fun(self):
        x = np.random.rand(3, 5, 2, 2)
        result = _default_score_fun(x)
        expected = np.mean(x[:, :, 0, 0], axis=1)
        np.testing.assert_array_equal(result, expected)


class TestComputeMeritValues:

    def test_scalar_maxcv_init(self):
        fun_values = np.array([1.0, 2.0, 3.0])
        maxcv_values = np.array([0.0, 0.0, 0.0])
        result = compute_merit_values(_default_merit, fun_values, maxcv_values, 0.0)
        np.testing.assert_array_equal(result, fun_values)

    def test_array_maxcv_init(self):
        fun_values = np.array([[1.0, 2.0], [3.0, 4.0]])
        maxcv_values = np.array([[0.0, 0.0], [0.0, 0.0]])
        maxcv_init = np.array([0.0, 0.0])
        result = compute_merit_values(_default_merit, fun_values, maxcv_values, maxcv_init)
        np.testing.assert_array_equal(result, fun_values)

    def test_nan_returns_inf(self):
        fun_values = np.array([np.nan])
        maxcv_values = np.array([0.0])
        result = compute_merit_values(_default_merit, fun_values, maxcv_values, 0.0)
        assert np.isinf(result[0])


class TestCreateStamp:

    def test_basic(self):
        solver_names = ['solver1', 'solver2']
        problem_options = {
            ProblemOption.PTYPE: 'u',
            ProblemOption.MINDIM: 1,
            ProblemOption.MAXDIM: 10,
            ProblemOption.MINCON: 0,
            ProblemOption.MAXCON: 5,
        }
        stamp = create_stamp(solver_names, problem_options, 'plain', '20240101_120000', '/tmp')
        assert isinstance(stamp, str)
        assert '20240101_120000' in stamp

    def test_unconstrained_stamp(self):
        solver_names = ['s1']
        problem_options = {
            ProblemOption.PTYPE: 'u',
            ProblemOption.MINDIM: 1,
            ProblemOption.MAXDIM: 10,
            ProblemOption.MINCON: 0,
            ProblemOption.MAXCON: 0,
        }
        stamp = create_stamp(solver_names, problem_options, 'plain', '20240101_120000', '/tmp')
        assert 'u_1_10' in stamp

    def test_constrained_stamp(self):
        solver_names = ['s1']
        problem_options = {
            ProblemOption.PTYPE: 'b',
            ProblemOption.MINDIM: 1,
            ProblemOption.MAXDIM: 10,
            ProblemOption.MINCON: 0,
            ProblemOption.MAXCON: 5,
        }
        stamp = create_stamp(solver_names, problem_options, 'plain', '20240101_120000', '/tmp')
        assert 'b_1_10_0_5' in stamp


class TestGetDefaultFeatureStamp:

    def test_plain(self):
        feature = Feature('plain')
        assert _get_default_feature_stamp(feature) == 'plain'

    def test_noisy(self):
        feature = Feature('noisy')
        stamp = _get_default_feature_stamp(feature)
        assert 'noisy' in stamp
        assert '0.001' in stamp

    def test_truncated(self):
        feature = Feature('truncated')
        stamp = _get_default_feature_stamp(feature)
        assert 'truncated' in stamp
        assert '6' in stamp

    def test_perturbed_x0(self):
        feature = Feature('perturbed_x0')
        stamp = _get_default_feature_stamp(feature)
        assert 'perturbed_x0' in stamp

    def test_random_nan(self):
        feature = Feature('random_nan')
        stamp = _get_default_feature_stamp(feature)
        assert 'random_nan' in stamp
        assert '0.05' in stamp

    def test_linearly_transformed_rotated(self):
        feature = Feature('linearly_transformed', rotated=True)
        stamp = _get_default_feature_stamp(feature)
        assert 'linearly_transformed' in stamp
        assert 'rotated' in stamp

    def test_linearly_transformed_conditioned(self):
        feature = Feature('linearly_transformed', rotated=False, condition_factor=2.0)
        stamp = _get_default_feature_stamp(feature)
        assert 'cond' in stamp

    def test_quantized(self):
        feature = Feature('quantized')
        stamp = _get_default_feature_stamp(feature)
        assert 'quantized' in stamp
        assert 'ground_truth' in stamp

    def test_unrelaxable_constraints(self):
        feature = Feature('unrelaxable_constraints', unrelaxable_bounds=True)
        stamp = _get_default_feature_stamp(feature)
        assert 'unrelaxable_constraints' in stamp
        assert 'bounds' in stamp

    def test_permuted(self):
        feature = Feature('permuted')
        stamp = _get_default_feature_stamp(feature)
        assert stamp == 'permuted'

    def test_custom(self):
        feature = Feature('custom', mod_fun=lambda x, rng, p: p.fun(x))
        stamp = _get_default_feature_stamp(feature)
        assert stamp == 'custom'


class TestIntegrate:

    def test_none_curve(self):
        assert integrate(None, 'perf', {}) == 0.0

    def test_empty_tuple(self):
        assert integrate((), 'perf', {}) == 0.0

    def test_empty_arrays(self):
        assert integrate((np.array([]), np.array([])), 'perf', {}) == 0.0

    def test_perf_step_function(self):
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([0.5, 0.8, 1.0])
        result = integrate((x, y), 'perf', {})
        expected = 0.5 * 1.0 + 0.8 * 1.0
        np.testing.assert_allclose(result, expected)

    def test_data_step_function(self):
        x = np.array([0.0, 1.0, 2.0])
        y = np.array([0.0, 0.5, 1.0])
        result = integrate((x, y), 'data', {})
        expected = 0.0 * 1.0 + 0.5 * 1.0
        np.testing.assert_allclose(result, expected)

    def test_log_ratio(self):
        x = np.array([1, 2, 3])
        y = np.array([0.5, -0.3, 0.2])
        result = integrate((x, y), 'log_ratio', {})
        expected = np.sum(np.abs(y))
        np.testing.assert_allclose(result, expected)

    def test_ndarray_input(self):
        curve = np.array([[1.0, 2.0, 3.0], [0.5, 0.8, 1.0]])
        result = integrate(curve, 'perf', {})
        expected = 0.5 * 1.0 + 0.8 * 1.0
        np.testing.assert_allclose(result, expected)

    def test_single_point(self):
        x = np.array([1.0])
        y = np.array([0.5])
        result = integrate((x, y), 'perf', {})
        assert result == 0.0

    def test_custom_weight_fun(self):
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([0.5, 0.8, 1.0])
        opts = {ProfileOption.SCORE_WEIGHT_FUN: lambda t: 2.0}
        result = integrate((x, y), 'perf', opts)
        expected = 2.0 * (0.5 * 1.0 + 0.8 * 1.0)
        np.testing.assert_allclose(result, expected)

    def test_unknown_type_returns_zero(self):
        x = np.array([1.0, 2.0])
        y = np.array([0.5, 0.8])
        result = integrate((x, y), 'unknown', {})
        assert result == 0.0


class TestInitReadme:

    def test_creates_file(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.md', delete=False) as f:
            path = f.name
        try:
            init_readme(path)
            with open(path) as f:
                content = f.read()
            assert '# Content of this folder' in content
            assert 'File/Folder Name' in content
        finally:
            os.unlink(path)


class TestAddToReadme:

    def test_adds_entry(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.md', delete=False) as f:
            path = f.name
        try:
            init_readme(path)
            add_to_readme(path, 'test.pdf', 'File, a test file.')
            with open(path) as f:
                content = f.read()
            assert 'test.pdf' in content
            assert 'File' in content
        finally:
            os.unlink(path)

    def test_folder_entry(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.md', delete=False) as f:
            path = f.name
        try:
            init_readme(path)
            add_to_readme(path, 'data/', 'Folder, contains data files.')
            with open(path) as f:
                content = f.read()
            assert 'data/' in content
            assert 'Folder' in content
        finally:
            os.unlink(path)


class TestMergePdfs:

    def test_no_pdfs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(FileNotFoundError):
                merge_pdfs_with_pypdf(tmpdir, os.path.join(tmpdir, 'out.pdf'))
