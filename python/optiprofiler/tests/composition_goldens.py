"""
Legacy single-feature golden scenarios.

This module is shared by the golden generator (run against the exact base
commit) and by ``test_composition_goldens.py`` (run against the current
code). Every scenario constructs a small problem whose root callbacks log
each invocation, applies one feature through the public ``Feature`` and
``FeaturedProblem`` classes, performs a fixed query sequence, and records
every observable outcome: structure, oracle outputs, reference histories,
budget behaviour and the exact sequence of root callback invocations.

The recorded fixture pins supported single-feature behaviour, including the
callback multiplicity of each feature, so that the composition work cannot
change it silently.

Run ``python composition_goldens.py --write <path>`` to (re)generate the
fixture with whatever ``optiprofiler`` is importable on ``PYTHONPATH``.
"""

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np

from optiprofiler.opclasses import Feature, FeaturedProblem, Problem

FIXTURE_VERSION = 1
MAX_EVAL = 4
SEEDS = (0, 12345)
KINDS = ('u', 'b', 'l', 'n')


class CallLog:
    """Record every root callback invocation per channel, in call order."""

    def __init__(self):
        self.calls = {'fun': [], 'cub': [], 'ceq': []}

    def wrap(self, channel, function):
        def wrapped(x):
            self.calls[channel].append(np.array(x, dtype=float).tolist())
            return function(x)
        return wrapped

    def clear(self):
        for channel in self.calls:
            self.calls[channel] = []


def root_fun(x):
    x = np.asarray(x, dtype=float)
    return float(np.sum(np.array([1.0, 2.0, 0.5]) * (x - np.array([0.5, -1.0, 2.0])) ** 2) + 0.5)


def root_cub(x):
    return np.array([x[0] ** 2 + x[1] - 1.5, x[2] - 3.0])


def root_ceq(x):
    return np.array([x[0] * x[1] - 0.2])


def custom_mod_fun(x, rng, problem):
    # Two genuine predecessor reads plus a stream draw: pins the legacy
    # custom contract, including the base-value dependent stream.
    return problem.fun(x) + 0.5 * problem.fun(x) + 1e-3 * rng.standard_normal()


def custom_mod_x0(rng, problem):
    return problem.x0 + 0.1 * rng.standard_normal(problem.n)


def custom_mod_affine(rng, problem):
    d = np.array([2.0, 0.5, 1.0])
    return np.diag(d), np.array([0.1, -0.2, 0.3]), np.diag(1.0 / d)


def custom_mod_cub(x, rng, problem):
    return problem.cub(x) + 0.25


FEATURE_CASES = [
    ('plain', {}),
    ('perturbed_x0', {}),
    ('perturbed_x0', {'distribution': 'gaussian', 'perturbation_level': 0.5}),
    ('noisy', {}),
    ('noisy', {'distribution': 'uniform', 'noise_type': 'relative', 'noise_level': 0.1}),
    ('noisy', {'noise_type': 'absolute', 'noise_level': 0.5}),
    ('noisy', {'noise_mode': 'deterministic'}),
    ('truncated', {}),
    ('truncated', {'significant_digits': 3, 'perturbed_trailing_digits': True}),
    ('permuted', {}),
    ('linearly_transformed', {}),
    ('linearly_transformed', {'rotated': False, 'condition_factor': 2.0}),
    ('linearly_transformed', {'condition_factor': 1.0}),
    ('random_nan', {}),
    ('random_nan', {'nan_rate': 0.5}),
    ('unrelaxable_constraints', {}),
    ('unrelaxable_constraints', {'unrelaxable_bounds': True, 'unrelaxable_linear_constraints': True,
                                 'unrelaxable_nonlinear_constraints': True}),
    ('nonquantifiable_constraints', {}),
    ('quantized', {}),
    ('quantized', {'mesh_size': 0.3, 'ground_truth': False}),
    ('quantized', {'mesh_size': 0.3, 'mesh_type': 'relative', 'ground_truth': True}),
    ('custom', {'mod_fun': custom_mod_fun, 'mod_x0': custom_mod_x0, 'mod_affine': custom_mod_affine,
                'mod_cub': custom_mod_cub}),
]


def case_id(name, options):
    if not options:
        return name
    return name + '[' + ','.join(f'{k}={getattr(v, "__name__", v)}' for k, v in sorted(options.items())) + ']'


def make_problem(kind, log):
    x0 = np.array([0.3, -1.2, 2.5])
    fun = log.wrap('fun', root_fun)
    if kind == 'u':
        return Problem(fun, x0, name='golden_u')
    xl = np.array([-2.0, -3.0, -1.0])
    xu = np.array([2.0, 3.0, 4.0])
    if kind == 'b':
        return Problem(fun, x0, name='golden_b', xl=xl, xu=xu)
    aub = np.array([[1.0, 1.0, 0.0]])
    bub = np.array([2.5])
    aeq = np.array([[0.0, 1.0, -1.0]])
    beq = np.array([-0.5])
    if kind == 'l':
        return Problem(fun, x0, name='golden_l', xl=xl, xu=xu, aub=aub, bub=bub, aeq=aeq, beq=beq)
    return Problem(fun, x0, name='golden_n', xl=xl, xu=xu, aub=aub, bub=bub, aeq=aeq, beq=beq,
                   cub=log.wrap('cub', root_cub), ceq=log.wrap('ceq', root_ceq))


def _listify(value):
    return np.asarray(value, dtype=float).tolist()


def _call(method, x):
    try:
        return _listify(method(x))
    except StopIteration:
        return 'StopIteration'


def run_scenario(kind, name, options, seed):
    log = CallLog()
    problem = make_problem(kind, log)
    log.clear()
    feature = Feature(name, **options)
    featured = FeaturedProblem(problem, feature, MAX_EVAL, seed)
    record = {
        'is_stochastic': bool(feature.is_stochastic),
        'options': {k: getattr(v, '__name__', v) for k, v in feature.options.items()},
        'n': int(featured.n),
        'ptype': featured.ptype,
        'm_nonlinear_ub': int(featured.m_nonlinear_ub),
        'm_nonlinear_eq': int(featured.m_nonlinear_eq),
        'x0': _listify(featured.x0),
        'xl': _listify(featured.xl),
        'xu': _listify(featured.xu),
        'aub': _listify(featured.aub),
        'bub': _listify(featured.bub),
        'aeq': _listify(featured.aeq),
        'beq': _listify(featured.beq),
        'fun_init': float(featured.fun_init),
        'maxcv_init': float(featured.maxcv_init),
        'construction_calls': {k: list(v) for k, v in log.calls.items()},
    }
    log.clear()
    delta = np.array([0.1, -0.1, 0.1])
    points = [featured.x0, featured.x0, featured.x0 + delta, np.array([5.0, -7.0, 9.0])]
    points.extend([featured.x0] * (2 * MAX_EVAL + 1 - len(points)))
    record['fun_outputs'] = [_call(featured.fun, x) for x in points]
    record['maxcv_outputs'] = [_call(featured.maxcv, x) for x in points[:4]]
    record['cub_outputs'] = [_call(featured.cub, x) for x in points[:2]]
    record['ceq_outputs'] = [_call(featured.ceq, x) for x in points[:2]]
    record['fun_hist'] = _listify(featured.fun_hist)
    record['maxcv_hist'] = _listify(featured.maxcv_hist)
    record['cub_hist'] = _listify(featured.cub_hist)
    record['ceq_hist'] = _listify(featured.ceq_hist)
    record['n_eval_fun'] = int(featured.n_eval_fun)
    record['n_eval_cub'] = int(featured.n_eval_cub)
    record['n_eval_ceq'] = int(featured.n_eval_ceq)
    record['query_calls'] = {k: list(v) for k, v in log.calls.items()}
    return record


def run_all_scenarios():
    scenarios = {}
    for kind in KINDS:
        for name, options in FEATURE_CASES:
            for seed in SEEDS:
                scenarios[f'{kind}/{case_id(name, options)}/seed={seed}'] = run_scenario(kind, name, options, seed)
    return scenarios


def provenance():
    import scipy
    source = Path(sys.modules['optiprofiler'].__file__).resolve().parents[2]
    try:
        sha = subprocess.run(['git', '-C', str(source), 'rev-parse', 'HEAD'], check=True,
                             capture_output=True, text=True).stdout.strip()
    except Exception:
        sha = None
    return {
        'fixture_version': FIXTURE_VERSION,
        'source_sha': sha,
        'python': platform.python_version(),
        'numpy': np.__version__,
        'scipy': scipy.__version__,
        'hostname': platform.node(),
        'scenario_module_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def _assert_close(actual, expected, path):
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and set(actual) == set(expected), f'{path}: keys differ'
        for key in expected:
            _assert_close(actual[key], expected[key], f'{path}/{key}')
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), f'{path}: length differs'
        if expected and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in expected) \
                and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in actual):
            np.testing.assert_allclose(np.array(actual, dtype=float), np.array(expected, dtype=float),
                                       rtol=1e-12, atol=0.0, equal_nan=True, err_msg=path)
        else:
            for i, (a, e) in enumerate(zip(actual, expected)):
                _assert_close(a, e, f'{path}[{i}]')
    elif isinstance(expected, float):
        np.testing.assert_allclose(float(actual), expected, rtol=1e-12, atol=0.0, equal_nan=True, err_msg=path)
    else:
        assert actual == expected, f'{path}: {actual!r} != {expected!r}'


def assert_matches(actual, expected):
    _assert_close(actual, expected, 'scenarios')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', required=True, help='Destination JSON file.')
    args = parser.parse_args(argv)
    document = {'provenance': provenance(), 'scenarios': run_all_scenarios()}
    Path(args.write).write_text(json.dumps(document, indent=1, sort_keys=True) + '\n', encoding='utf-8')
    print(f"wrote {args.write}: {len(document['scenarios'])} scenarios, source {document['provenance']['source_sha']}")


if __name__ == '__main__':
    main()
