"""Native replay entry point, independent of historical decoding internals."""

import pickle

import numpy as np
import pytest

from optiprofiler import Problem, benchmark, load_benchmark_options, legacy_compat


def objective(x):
    return float(x @ x)


def stay(fun, x0):
    fun(x0)
    return x0


@pytest.mark.parametrize('schema', ['options_user-v2', 'options_refined-v2'])
@pytest.mark.parametrize('path_type', [str, lambda path: path])
def test_file_and_mapping_replay_directly_into_benchmark(tmp_path, schema, path_type):
    problem = Problem(objective, [1., 2.], name='replay')
    source = dict(schema=schema, problem=problem, problem_route='direct_problem',
                  n_runs=1, silent=False, seed=19, savepath='old-output')
    if schema == 'options_user-v2':
        source['feature_name'] = 'plain'
    else:
        source['feature_specification'] = [{'name': 'plain', 'options': {}}]
    path = tmp_path / 'options.pkl'
    original = pickle.dumps(source)
    path.write_bytes(original)

    mapping_options = load_benchmark_options(source)
    file_options = load_benchmark_options(path_type(path))
    assert set(mapping_options) == set(file_options) == {'feature', 'problem', 'n_runs'}
    assert mapping_options['problem'] is problem
    assert mapping_options['feature'] == file_options['feature']
    settings = dict(silent=True, score_only=True, n_jobs=1, run_plain=False,
                    max_eval_factor=1, max_tol_order=1)
    expected = benchmark([stay, stay], **mapping_options, **settings)[0]
    actual = benchmark([stay, stay], **file_options, **settings)[0]
    np.testing.assert_array_equal(actual, expected)
    assert pickle.dumps(source) == original
    assert path.read_bytes() == original


def test_old_replay_entry_points_are_not_retained():
    assert not hasattr(legacy_compat, 'load_options')
    assert not hasattr(legacy_compat, 'replay_arguments')
