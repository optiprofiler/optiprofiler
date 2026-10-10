#!/usr/bin/env python3
"""Exercise an installed distribution, outside its source checkout.

This is a headless release gate, not a native GUI-backend test. It exercises
actual SciPy solvers, bundled selection, composed features, spawned workers,
PDF output and solver-free reload. Full unit suites remain separate gates.
"""

import argparse
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import sys

os.environ['MPLBACKEND'] = 'Agg'

import numpy as np
from jsonschema import Draft202012Validator
from pypdf import PdfReader
from scipy.optimize import minimize

import optiprofiler
from optiprofiler import Feature, FeaturedProblem, Problem, ProblemReference, benchmark
from optiprofiler.eval_report import load_schema_for


def powell(fun, x0):
    return minimize(fun, x0, method='Powell', options={'maxfev': 20}).x


def nelder_mead(fun, x0):
    return minimize(fun, x0, method='Nelder-Mead', options={'maxfev': 20}).x


def sphere(x):
    return float(x @ x)


def halfway(fun, x0, *constraints):
    x = x0 / 2
    fun(x)
    return x


def unchanged(fun, x0, *constraints):
    fun(x0)
    return x0


def linear_merit(fun, maxcv, initial_maxcv):
    return fun + 2 * maxcv


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_report(path):
    report = json.loads(path.read_text(encoding='utf-8'))
    Draft202012Validator(load_schema_for(report)).validate(report)
    assert report['status'] == 'completed', report
    return report


def smoke(output, jobs):
    root = output / ('jobs-' + str(jobs))
    root.mkdir()
    report_path = root / 'saved.json'
    options = dict(plibs=['s2mpj'], problem_names=['ROSENBR'], ptype='u',
                   mindim=2, maxdim=2, n_runs=2, n_jobs=jobs, max_eval_factor=20,
                   max_tol_order=2, benchmark_id='installed', savepath=str(root),
                   silent=True, draw_hist_plots='parallel' if jobs == 2 else 'sequential',
                   run_plain=True, report_path=report_path)
    original = benchmark([powell, nelder_mead],
                         feature=Feature([{'name': 'noisy', 'options': {'noise_level': 0.001}},
                                          'permuted']), **options)
    report = read_report(report_path)
    primary = [problem for problem in report['problems'] if problem['role'] == 'primary']
    plain = [problem for problem in report['problems'] if problem['role'] == 'plain_reference']
    assert len(primary) == len(plain) == 1, report['problems']
    assert primary[0]['name'] == 'ROSENBR' and len(primary[0]['runs']) == 4
    assert len(plain[0]['runs']) == 2
    for problem in report['problems']:
        for run in problem['runs']:
            assert run['evaluations'] > 1, run
            assert not run['abnormal_termination'] and not run['output_fallback'], run
    assert np.all(np.isfinite(original[0])), original[0]
    archives = list((root / 'installed').rglob('data_for_loading.h5'))
    assert len(archives) == 1, archives
    preserved = {path: sha256(path) for path in (root / 'installed').rglob('*') if path.is_file()}
    pdfs = list((root / 'installed').rglob('*.pdf'))
    assert any('summary' in path.name for path in pdfs), pdfs
    assert any('hist' in str(path).lower() for path in pdfs), pdfs
    assert all(len(PdfReader(str(path)).pages) > 0 for path in pdfs)
    # No callables are supplied: a successful load cannot accidentally rerun
    # the solvers. Verify scores and immutable original artifacts, not a banner.
    loaded_path = root / 'loaded.json'
    # load searches cwd/benchmark_id; savepath controls new output, not the
    # input archive. Exercise the documented route rather than assuming that
    # an output option also redirects archive discovery.
    previous_cwd = Path.cwd()
    try:
        os.chdir(root)
        loaded = benchmark(None, load='latest', benchmark_id='installed', savepath=str(root),
                           score_only=True, silent=True, draw_hist_plots='none',
                           max_tol_order=2, max_eval_factor=20, run_plain=True,
                           report_path=loaded_path)
    finally:
        os.chdir(previous_cwd)
    np.testing.assert_array_equal(loaded[0], original[0])
    np.testing.assert_array_equal(loaded[1], original[1])
    loaded_report = read_report(loaded_path)
    assert loaded_report['operation'] == 'load'
    assert loaded_report['source']['sha256'] == preserved[archives[0]]
    assert all(path.is_file() and sha256(path) == value for path, value in preserved.items())
    return {'jobs': jobs, 'archive_sha256': preserved[archives[0]],
            'primary_runs': 4, 'plain_runs': 2, 'readable_pdfs': len(pdfs),
            'solver_free_reload': 'PASS', 'source_artifacts_unchanged': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--expected-version', required=True)
    args = parser.parse_args()
    installed = Path(optiprofiler.__file__).resolve()
    source = args.source_root.resolve()
    assert source not in installed.parents, 'The smoke imported the source checkout, not the installed artifact.'
    assert optiprofiler.__version__ == metadata.version('optiprofiler') == args.expected_version
    output = args.output_dir.resolve()
    output.mkdir(parents=True)
    os.chdir(output)
    # Check scoring, not just metadata. The same measured step has relative
    # decrease 1 under the empirical target and 3/4 under a supplied zero.
    reference = ProblemReference(0.0, 'optimum', 'analytic sphere', 'feasible_objective/1')
    problem = Problem(sphere, [1.0, 2.0], reference=reference)
    trial = FeaturedProblem(problem, Feature('noisy+permuted'), 10, 7)
    assert trial.reference == reference
    assert np.isfinite(trial.fun(trial.x0))
    reference_options = dict(score_only=True, silent=True, n_jobs=1, max_eval_factor=2)
    empirical = benchmark([halfway, unchanged], problem=Problem(sphere, [2.0]), **reference_options)
    explicit = benchmark([halfway, unchanged], problem=Problem(sphere, [2.0], reference=reference),
                         **reference_options)
    np.testing.assert_array_equal(empirical[0], [1.0, 0.0])
    np.testing.assert_array_equal(explicit[0], [0.75, 0.0])
    violated = ProblemReference(fun=0.0, maxcv=0.5, kind='target', source='test target')
    converted = benchmark([halfway, unchanged], problem=Problem(sphere, [2.0], xl=[0.5], reference=violated),
                          merit_fun=linear_merit, **reference_options)
    np.testing.assert_array_equal(converted[0], [1.0, 0.0])
    records = [smoke(output, jobs) for jobs in (1, 2)]
    receipt = {'status': 'PASS', 'version': args.expected_version, 'python': sys.version,
               'installed_package': str(installed), 'backend': 'Agg',
               'reference_contract': 'raw objective/violation propagation and current-merit target scoring',
               'cases': records}
    (output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
