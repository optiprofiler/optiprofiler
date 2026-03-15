"""Tests for the loader module."""
import os
import pickle
import tempfile

import h5py
import numpy as np
import pytest

from optiprofiler.loader import (
    load_results_from_h5,
    merit_fun_compute,
    recompute_merits,
    save_options,
    save_results_to_h5,
    search_in_dir,
    truncate_problems,
    truncate_solvers,
)
from optiprofiler.profile_utils import _default_merit, compute_merit_values
from optiprofiler.utils import ProblemOption


def _make_dummy_results_plib(n_problems=5, n_solvers=2, n_runs=3, n_evals=100):
    """Create a dummy results_plib dictionary for testing."""
    return {
        'plib': 's2mpj',
        'solver_names': [f'solver_{i}' for i in range(n_solvers)],
        'fun_histories': np.random.rand(n_problems, n_solvers, n_runs, n_evals),
        'maxcv_histories': np.zeros((n_problems, n_solvers, n_runs, n_evals)),
        'fun_outs': np.random.rand(n_problems, n_solvers, n_runs),
        'maxcv_outs': np.zeros((n_problems, n_solvers, n_runs)),
        'fun_inits': np.random.rand(n_problems),
        'maxcv_inits': np.zeros(n_problems),
        'n_evals': np.full((n_problems, n_solvers, n_runs), n_evals),
        'computation_times': np.random.rand(n_problems, n_solvers, n_runs),
        'solvers_successes': np.ones((n_problems, n_solvers, n_runs), dtype=bool),
        'merit_histories': np.random.rand(n_problems, n_solvers, n_runs, n_evals),
        'merit_outs': np.random.rand(n_problems, n_solvers, n_runs),
        'merit_inits': np.random.rand(n_problems),
        'problem_names': [f'PROB_{i}' for i in range(n_problems)],
        'problem_types': ['u'] * n_problems,
        'problem_dims': np.array([2, 5, 10, 20, 50][:n_problems]),
        'problem_mbs': np.zeros(n_problems, dtype=int),
        'problem_mlcons': np.zeros(n_problems, dtype=int),
        'problem_mnlcons': np.zeros(n_problems, dtype=int),
        'problem_mcons': np.zeros(n_problems, dtype=int),
        'feature_stamp': 'plain',
        'ptype': 'u',
        'mindim': 1,
        'maxdim': 100,
        'minb': 0,
        'maxb': 10,
        'minlcon': 0,
        'maxlcon': 10,
        'minnlcon': 0,
        'maxnlcon': 10,
        'mincon': 0,
        'maxcon': 10,
    }


class TestSearchInDir:

    def test_finds_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for i in range(3):
                with open(os.path.join(tmpdir, f'test_{i}.txt'), 'w') as f:
                    f.write('hello')
            files = search_in_dir(tmpdir, '*.txt', 1, 0, [])
            assert len(files) == 3
            for fi in files:
                assert fi['name'].endswith('.txt')
                assert fi['isdir'] is False

    def test_respects_max_depth(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            subdir = os.path.join(tmpdir, 'sub')
            os.makedirs(subdir)
            with open(os.path.join(subdir, 'deep.txt'), 'w') as f:
                f.write('hello')
            files = search_in_dir(tmpdir, '*.txt', 0, 0, [])
            assert len(files) == 0

    def test_recursive_search(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            subdir = os.path.join(tmpdir, 'sub')
            os.makedirs(subdir)
            with open(os.path.join(subdir, 'deep.txt'), 'w') as f:
                f.write('hello')
            files = search_in_dir(tmpdir, '*.txt', 2, 0, [])
            assert len(files) == 1
            assert files[0]['name'] == 'deep.txt'

    def test_no_match(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with open(os.path.join(tmpdir, 'test.txt'), 'w') as f:
                f.write('hello')
            files = search_in_dir(tmpdir, '*.csv', 1, 0, [])
            assert len(files) == 0


class TestTruncateSolvers:

    def test_select_first_solver(self):
        results = _make_dummy_results_plib(n_problems=3, n_solvers=3)
        truncated = truncate_solvers(results, [0])
        assert len(truncated['solver_names']) == 1
        assert truncated['fun_histories'].shape[1] == 1
        assert truncated['merit_histories'].shape[1] == 1

    def test_select_multiple_solvers(self):
        results = _make_dummy_results_plib(n_problems=3, n_solvers=4)
        truncated = truncate_solvers(results, [1, 3])
        assert len(truncated['solver_names']) == 2
        assert truncated['solver_names'] == ['solver_1', 'solver_3']
        assert truncated['fun_histories'].shape[1] == 2

    def test_does_not_modify_original(self):
        results = _make_dummy_results_plib(n_problems=3, n_solvers=3)
        original_names = list(results['solver_names'])
        truncate_solvers(results, [0])
        assert results['solver_names'] == original_names


class TestTruncateProblems:

    def test_filter_by_ptype(self):
        results = _make_dummy_results_plib(n_problems=5)
        results['problem_types'] = ['u', 'b', 'u', 'l', 'n']
        truncated = truncate_problems(results, {ProblemOption.PTYPE.value: 'u'})
        assert len(truncated['problem_names']) == 2

    def test_filter_by_mindim(self):
        results = _make_dummy_results_plib(n_problems=5)
        truncated = truncate_problems(results, {ProblemOption.MINDIM.value: 10})
        for dim in truncated['problem_dims']:
            assert dim >= 10

    def test_filter_by_maxdim(self):
        results = _make_dummy_results_plib(n_problems=5)
        truncated = truncate_problems(results, {ProblemOption.MAXDIM.value: 10})
        for dim in truncated['problem_dims']:
            assert dim <= 10

    def test_filter_by_problem_names(self):
        results = _make_dummy_results_plib(n_problems=5)
        truncated = truncate_problems(results, {ProblemOption.PROBLEM_NAMES.value: ['PROB_0', 'PROB_2']})
        assert len(truncated['problem_names']) == 2
        assert 'PROB_0' in truncated['problem_names']
        assert 'PROB_2' in truncated['problem_names']

    def test_filter_by_excludelist(self):
        results = _make_dummy_results_plib(n_problems=5)
        truncated = truncate_problems(results, {ProblemOption.EXCLUDELIST.value: ['PROB_0', 'PROB_1']})
        assert len(truncated['problem_names']) == 3
        assert 'PROB_0' not in truncated['problem_names']
        assert 'PROB_1' not in truncated['problem_names']

    def test_filter_by_minb(self):
        results = _make_dummy_results_plib(n_problems=5)
        results['problem_mbs'] = np.array([0, 2, 4, 6, 8])
        truncated = truncate_problems(results, {ProblemOption.MINB.value: 4})
        for mb in truncated['problem_mbs']:
            assert mb >= 4

    def test_filter_by_maxb(self):
        results = _make_dummy_results_plib(n_problems=5)
        results['problem_mbs'] = np.array([0, 2, 4, 6, 8])
        truncated = truncate_problems(results, {ProblemOption.MAXB.value: 4})
        for mb in truncated['problem_mbs']:
            assert mb <= 4

    def test_empty_result(self):
        results = _make_dummy_results_plib(n_problems=5)
        truncated = truncate_problems(results, {ProblemOption.MINDIM.value: 1000})
        assert len(truncated['problem_names']) == 0

    def test_does_not_modify_original(self):
        results = _make_dummy_results_plib(n_problems=5)
        original_names = list(results['problem_names'])
        truncate_problems(results, {ProblemOption.MINDIM.value: 10})
        assert results['problem_names'] == original_names


class TestMeritFunCompute:

    def test_basic(self):
        def my_merit(f, v, v0):
            return f + v
        fun_values = np.array([1.0, 2.0])
        maxcv_values = np.array([0.1, 0.2])
        maxcv_inits = np.array([0.0, 0.0])
        result = merit_fun_compute(my_merit, fun_values, maxcv_values, maxcv_inits)
        np.testing.assert_array_almost_equal(result, [1.1, 2.2])


class TestRecomputeMerits:

    def test_basic(self):
        results = _make_dummy_results_plib(n_problems=3, n_solvers=2, n_runs=1, n_evals=10)
        results['fun_histories'] = np.ones((3, 2, 1, 10))
        results['maxcv_histories'] = np.zeros((3, 2, 1, 10))
        results['fun_outs'] = np.ones((3, 2, 1))
        results['maxcv_outs'] = np.zeros((3, 2, 1))
        results['fun_inits'] = np.ones(3)
        results['maxcv_inits'] = np.zeros(3)

        def my_merit(fun_values, maxcv_values, maxcv_inits):
            return fun_values + maxcv_values

        recomputed = recompute_merits(results, my_merit)
        assert 'merit_histories' in recomputed
        assert 'merit_outs' in recomputed
        assert 'merit_inits' in recomputed

    def test_invalid_merit_fun(self):
        results = _make_dummy_results_plib(n_problems=3)
        def bad_merit(fun_values, maxcv_values, maxcv_inits):
            raise RuntimeError("bad")
        with pytest.raises(ValueError, match='not valid'):
            recompute_merits(results, bad_merit)


class TestSaveLoadH5:

    def test_round_trip(self):
        results_plibs = [_make_dummy_results_plib(n_problems=3, n_solvers=2)]
        with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as f:
            path = f.name
        try:
            save_results_to_h5(results_plibs, path)
            loaded = load_results_from_h5(path)
            assert len(loaded) == 1
            assert loaded[0]['plib'] == 's2mpj'
            assert len(loaded[0]['solver_names']) == 2
            assert len(loaded[0]['problem_names']) == 3
            np.testing.assert_array_almost_equal(
                loaded[0]['fun_histories'],
                results_plibs[0]['fun_histories'],
            )
        finally:
            os.unlink(path)

    def test_multiple_plibs(self):
        plib1 = _make_dummy_results_plib(n_problems=2)
        plib1['plib'] = 'lib1'
        plib2 = _make_dummy_results_plib(n_problems=3)
        plib2['plib'] = 'lib2'
        results_plibs = [plib1, plib2]
        with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as f:
            path = f.name
        try:
            save_results_to_h5(results_plibs, path)
            loaded = load_results_from_h5(path)
            assert len(loaded) == 2
        finally:
            os.unlink(path)

    def test_nested_dict(self):
        results_plibs = [_make_dummy_results_plib(n_problems=2)]
        results_plibs[0]['nested'] = {
            'key1': np.array([1.0, 2.0]),
            'key2': ['a', 'b'],
        }
        with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as f:
            path = f.name
        try:
            save_results_to_h5(results_plibs, path)
            loaded = load_results_from_h5(path)
            assert 'nested' in loaded[0]
            np.testing.assert_array_equal(loaded[0]['nested']['key1'], [1.0, 2.0])
            assert loaded[0]['nested']['key2'] == ['a', 'b']
        finally:
            os.unlink(path)

    def test_scalar_string(self):
        results_plibs = [_make_dummy_results_plib(n_problems=2)]
        results_plibs[0]['feature_stamp'] = 'plain'
        with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as f:
            path = f.name
        try:
            save_results_to_h5(results_plibs, path)
            loaded = load_results_from_h5(path)
            assert loaded[0]['feature_stamp'] == 'plain'
        finally:
            os.unlink(path)


class TestSaveOptions:

    def test_basic(self):
        opts = {'key1': 'value1', 'key2': 42}
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as f:
            path = f.name
        try:
            save_options(opts, path)
            with open(path, 'rb') as f:
                loaded = pickle.load(f)
            assert loaded == opts
        finally:
            os.unlink(path)

    def test_with_non_picklable(self):
        opts = {'merit_fun': lambda f, v, v0: f}
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as f:
            path = f.name
        try:
            save_options(opts, path)
            assert os.path.exists(path)
        finally:
            os.unlink(path)
