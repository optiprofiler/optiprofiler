"""Public new-run and saved-run regressions for undefined evaluations."""

import hashlib

import h5py
import numpy as np
import pytest

from optiprofiler import benchmark
from optiprofiler.loader import load_results_from_h5


def visit_one(fun, x0, *constraints):
    x = np.ones_like(x0)
    fun(x)
    return x


def visit_one_again(fun, x0, *constraints):
    return visit_one(fun, x0, *constraints)


def invalid_then_valid(fun, x0, *constraints):
    fun(np.ones_like(x0))
    fun(x0)
    return x0


def ignore_constraint_merit(fun, cv, cv_init):
    return 0.0


@pytest.fixture
def nan_library(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    root = tmp_path / 'nanvalidity'
    root.mkdir()
    (root / 'nanvalidity_tools.py').write_text('''
import numpy as np
from optiprofiler import Problem

def nanvalidity_select(options):
    return ['FINITE', 'NAN_INIT_CV', 'NAN_INIT_FUN', 'NAN_LATER_CV', 'NAN_LATER_FUN', 'VALID_INF', 'INF_INIT_NAN_LATER']

def nanvalidity_load(name):
    def fun(x):
        if name == 'NAN_INIT_FUN' and x[0] == 0:
            return np.nan
        if name == 'NAN_LATER_FUN' and x[0] != 0:
            return np.nan
        if name == 'VALID_INF' or (name == 'INF_INIT_NAN_LATER' and x[0] == 0):
            return np.inf
        return 7.0
    def cub(x):
        if (name == 'NAN_INIT_CV' and x[0] == 0) or (name in ('NAN_LATER_CV', 'INF_INIT_NAN_LATER') and x[0] != 0):
            return np.array([np.nan])
        return np.array([-1.0])
    return Problem(fun, [0.0], cub=cub, name=name)
''', encoding='utf-8')
    return dict(plibs=['nanvalidity'], custom_problem_libs_path=str(tmp_path),
                ptype='n', mindim=1, maxdim=1, feature_name='plain', n_runs=2,
                n_jobs=1, max_eval_factor=3, max_tol_order=1,
                draw_hist_plots='none', solver_names=['one', 'other'],
                benchmark_id='nan-scoring', savepath=str(tmp_path), silent=True)


def passing_fraction(curves, channel, solver=0):
    curve = curves[0][channel]['perf'][solver][-1]
    return 0.0 if curve is None else float(curve[1][-1])


@pytest.mark.parametrize('name', ['NAN_INIT_CV', 'NAN_INIT_FUN', 'NAN_LATER_CV', 'NAN_LATER_FUN', 'INF_INIT_NAN_LATER'])
@pytest.mark.parametrize('custom_merit', [False, True])
def test_undefined_evaluations_cannot_pass(nan_library, name, custom_merit):
    options = dict(nan_library, problem_names=[name], score_only=True)
    if custom_merit:
        options['merit_fun'] = ignore_constraint_merit
    _, _, curves = benchmark([visit_one, visit_one_again], **options)
    assert passing_fraction(curves, 'hist') == 0.0
    assert passing_fraction(curves, 'out') == 0.0


@pytest.mark.parametrize('name', ['FINITE', 'VALID_INF'])
def test_valid_initial_conventions_are_unchanged(nan_library, name):
    _, _, curves = benchmark([visit_one, visit_one_again],
                             **dict(nan_library, problem_names=[name], score_only=True))
    assert passing_fraction(curves, 'hist') == 1.0
    assert passing_fraction(curves, 'out') == 1.0


def test_valid_points_after_nan_still_count(nan_library):
    _, _, curves = benchmark([invalid_then_valid, visit_one],
                             **dict(nan_library, problem_names=['INF_INIT_NAN_LATER'], score_only=True))
    assert passing_fraction(curves, 'hist', 0) == 1.0
    assert passing_fraction(curves, 'out', 0) == 1.0
    assert passing_fraction(curves, 'hist', 1) == 0.0
    assert passing_fraction(curves, 'out', 1) == 0.0
    x, y = curves[0]['hist']['data'][0][-1]
    assert x[np.flatnonzero(y > 0)[0]] == 1.0  # second evaluation / (dimension + 1)


@pytest.mark.parametrize('legacy', [False, True])
def test_saved_raw_nan_validity_survives_load(nan_library, tmp_path, legacy):
    options = dict(nan_library, problem_names=['FINITE', 'NAN_INIT_CV', 'INF_INIT_NAN_LATER'],
                   score_only=False)
    _, _, original_curves = benchmark([visit_one, visit_one_again], **options)
    source = next(tmp_path.glob('nan-scoring/**/data_for_loading.h5'))
    if legacy:
        # Emulate a pre-fix saved file: one initial value per problem, and
        # a custom merit that had assigned zero to undefined evaluations.
        with h5py.File(source, 'r+') as saved:
            group = next(iter(saved.values()))
            for field in ['fun_inits', 'maxcv_inits', 'merit_inits']:
                values = group[field][...][:, 0]
                del group[field]
                group.create_dataset(field, data=values)
            for suffix in ['histories', 'outs', 'inits']:
                invalid = np.isnan(group['fun_' + suffix][...]) | np.isnan(group['maxcv_' + suffix][...])
                values = group['merit_' + suffix][...]
                values[invalid] = 0.0
                group['merit_' + suffix][...] = values
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    data = load_results_from_h5(source)
    row = data[0]['problem_names'].index('NAN_INIT_CV')
    assert np.isnan(data[0]['maxcv_inits'][row]).all()
    _, _, loaded_curves = benchmark(None, **dict(options, load='latest', score_only=True))
    for channel in ['hist', 'out']:
        assert passing_fraction(original_curves, channel) == pytest.approx(1 / 3)
        assert passing_fraction(loaded_curves, channel) == pytest.approx(1 / 3)
    assert hashlib.sha256(source.read_bytes()).hexdigest() == digest
    report = next(source.parent.glob('*report*.txt')).read_text()
    assert 'undefined initial objective or constraint' in report
