"""Independent-review follow-ups: actual worker transport and failure boundaries."""
import json
import multiprocessing
import os
from pathlib import Path

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem, benchmark


def objective(x):
    return float(x[0])


def stay(fun, x0):
    fun(x0)
    return x0


def worker_grid(feature):
    import optiprofiler.opclasses as classes
    fp = FeaturedProblem(Problem(objective, [2.6]), feature, 5, 17)
    return os.getpid(), classes.__file__, fp.fun(np.array([2.6]))


@pytest.mark.parametrize('composed', [False, True])
def test_fresh_relative_feature_uses_current_grid_in_spawned_worker(composed):
    import optiprofiler.opclasses as classes
    entries = [{'name': 'quantized', 'options': {
        'mesh_size': 0.5, 'mesh_type': 'RELATIVE', 'ground_truth': False}}]
    if composed:
        entries.append({'name': 'noisy', 'options': {'noise_level': 0.0}})
    with multiprocessing.get_context('spawn').Pool(1) as pool:
        pid, origin, value = pool.apply_async(worker_grid, (Feature(entries),)).get(timeout=30)
    assert pid != os.getpid()
    assert Path(origin).resolve() == Path(classes.__file__).resolve()
    assert value == pytest.approx(2.6)


@pytest.mark.parametrize('level', [[0.1, 0.2], (0.1, 0.2), np.array([0.1, 0.2])])
def test_vector_amplitude_composition_has_real_benchmark_stamp_and_report(level, tmp_path):
    spec = Feature([{'name': 'perturbed_x0', 'options': {'perturbation_level': level}},
                    {'name': 'noisy', 'options': {'noise_level': 0.0}}])
    base = Problem(objective, [1.0, 2.0])
    fp = FeaturedProblem(base, spec, 5, 17)
    control = Feature([{'name': 'perturbed_x0', 'options': {
        'perturbation_level': np.array([0.1, 0.2])}},
        {'name': 'noisy', 'options': {'noise_level': 0.0}}])
    np.testing.assert_array_equal(fp.x0, FeaturedProblem(base, control, 5, 17).x0)
    report_path = tmp_path / 'report.json'
    # A provider benchmark exercises persisted H5/stamp output; the direct
    # problem route does not promise a data_for_loading.h5 archive.
    benchmark([stay, stay], plibs=['s2mpj'], problem_names=['ROSENBR'],
              ptype='u', mindim=2, maxdim=2, max_eval_factor=2,
              feature=spec, n_runs=1, n_jobs=1,
              silent=True, draw_hist_plots='none', benchmark_id='vector-report',
              savepath=str(tmp_path), report_path=report_path)
    report = json.loads(report_path.read_text())
    assert report['status'] == 'completed'
    assert report['coverage']['completed'] == 1
    assert len(list(tmp_path.rglob('data_for_loading.h5'))) == 1
    options = report['configuration']['effective']['feature']['stages'][0]['options']
    assert options['perturbation_level'] == [0.1, 0.2]


@pytest.mark.parametrize('n_jobs', [1, 2])
def test_mixed_dimension_library_rejects_incompatible_vector_without_silent_skip(n_jobs, tmp_path):
    name = f'dimensionprobe{n_jobs}'
    library = tmp_path / 'libraries' / name
    library.mkdir(parents=True)
    (library / f'{name}_tools.py').write_text(
        'from optiprofiler import Problem\n'
        f'def {name}_select(options):\n    return ["DIM2", "DIM3"]\n'
        f'def {name}_load(name):\n'
        '    n = int(name[-1])\n'
        '    return Problem(lambda x: float(x @ x), [1.0] * n, name=name)\n')
    report_path = tmp_path / 'failed.json'
    with pytest.raises(ValueError, match='perturbation_level'):
        benchmark([stay, stay], plibs=[name], custom_problem_libs_path=library.parent,
                  feature=Feature('perturbed_x0', perturbation_level=[0.1, 0.2]),
                  n_runs=1, n_jobs=n_jobs, mindim=2, maxdim=3, ptype='u',
                  silent=True, score_only=True, draw_hist_plots='none',
                  savepath=str(tmp_path), report_path=report_path)
    assert json.loads(report_path.read_text())['status'] != 'completed'
