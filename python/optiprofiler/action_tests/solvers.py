"""Solvers for action tests."""
import json
from pathlib import Path
import tempfile

import numpy as np
from scipy import __version__ as scipy_version
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint, show_options

from optiprofiler import benchmark


def scipy_cobyla(fun, x0, xl=None, xu=None, aub=None, bub=None, aeq=None, beq=None, cub=None, ceq=None):
    """Solver using scipy.optimize.minimize with COBYLA method."""
    # Build constraints list for COBYLA
    constraints = []
    
    # Bound constraints as inequality constraints for COBYLA
    if xl is not None:
        for i in range(len(xl)):
            if np.isfinite(xl[i]):
                constraints.append({'type': 'ineq', 'fun': lambda x, i=i, lb=xl[i]: x[i] - lb})
    if xu is not None:
        for i in range(len(xu)):
            if np.isfinite(xu[i]):
                constraints.append({'type': 'ineq', 'fun': lambda x, i=i, ub=xu[i]: ub - x[i]})
    
    # Linear inequality constraints: aub @ x <= bub
    if aub is not None and bub is not None and len(bub) > 0:
        for i in range(len(bub)):
            constraints.append({'type': 'ineq', 'fun': lambda x, i=i: bub[i] - aub[i] @ x})
    
    # Linear equality constraints: aeq @ x == beq
    if aeq is not None and beq is not None and len(beq) > 0:
        for i in range(len(beq)):
            constraints.append({'type': 'eq', 'fun': lambda x, i=i: aeq[i] @ x - beq[i]})
    
    # Nonlinear inequality constraints: cub(x) <= 0
    if cub is not None:
        def cub_constraint(x):
            val = cub(x)
            return -np.atleast_1d(val)  # COBYLA wants ineq >= 0
        constraints.append({'type': 'ineq', 'fun': cub_constraint})
    
    # Nonlinear equality constraints: ceq(x) == 0
    if ceq is not None:
        constraints.append({'type': 'eq', 'fun': ceq})
    
    result = minimize(fun, x0, method='COBYLA', constraints=constraints if constraints else ())
    return result.x


def scipy_cobyqa(fun, x0, xl=None, xu=None, aub=None, bub=None, aeq=None, beq=None, cub=None, ceq=None):
    """Solver using scipy.optimize.minimize with COBYQA method."""
    # Set up bounds
    bounds = None
    if xl is not None or xu is not None:
        lb = xl if xl is not None else np.full(len(x0), -np.inf)
        ub = xu if xu is not None else np.full(len(x0), np.inf)
        bounds = Bounds(lb, ub)
    
    # Build constraints
    constraints = []
    
    # Linear inequality constraints
    if aub is not None and bub is not None and len(bub) > 0:
        constraints.append(LinearConstraint(aub, -np.inf, bub))
    
    # Linear equality constraints
    if aeq is not None and beq is not None and len(beq) > 0:
        constraints.append(LinearConstraint(aeq, beq, beq))
    
    # Nonlinear inequality constraints
    if cub is not None:
        constraints.append(NonlinearConstraint(cub, -np.inf, 0))
    
    # Nonlinear equality constraints
    if ceq is not None:
        constraints.append(NonlinearConstraint(ceq, 0, 0))
    
    result = minimize(fun, x0, method='COBYQA', bounds=bounds, 
                     constraints=constraints if constraints else ())
    return result.x


def scipy_nelder_mead(fun, x0, xl=None, xu=None, aub=None, bub=None, aeq=None, beq=None, cub=None, ceq=None):
    """Solver using scipy.optimize.minimize with Nelder-Mead method.
    
    Note: Nelder-Mead only supports unconstrained problems. Bounds and constraints are ignored.
    """
    # Nelder-Mead doesn't support bounds or constraints
    result = minimize(fun, x0, method='Nelder-Mead')
    return result.x


def get_solver_configuration(unconstrained=False):
    """Select methods provided by the installed SciPy, without solving at import."""
    try:
        show_options(solver='minimize', method='cobyqa', disp=False)
    except ValueError:
        # SciPy releases supporting older Python versions may lack COBYQA.
        # Keep two real solvers instead of benchmarking an unknown method.
        if unconstrained:
            return [scipy_nelder_mead, scipy_cobyla], ['Nelder-Mead', 'COBYLA'], ['COBYQA']
        return [scipy_cobyla, scipy_nelder_mead], ['COBYLA', 'Nelder-Mead'], ['COBYQA']
    if unconstrained:
        return [scipy_nelder_mead, scipy_cobyqa], ['Nelder-Mead', 'COBYQA'], []
    return [scipy_cobyla, scipy_cobyqa, scipy_nelder_mead], ['COBYLA', 'COBYQA', 'Nelder-Mead'], []


def verify_solver_support(solvers, solver_names):
    """Fail before benchmark can catch a missing method as a solver failure."""
    for solver, name in zip(solvers, solver_names):
        evaluations = 0

        def objective(x):
            nonlocal evaluations
            evaluations += 1
            return float(np.dot(x, x))

        try:
            point = np.asarray(solver(objective, np.array([1., -1.])))
            if not evaluations:
                raise ValueError('no objective evaluations')
            if point.shape != (2,) or not np.all(np.isfinite(point)):
                raise ValueError('invalid returned point')
        except Exception as exc:
            raise RuntimeError(f'Action-test solver preflight failed for {name}: {exc}') from exc
        print(f'Solver preflight {name}: {evaluations} objective evaluations.')


def verify_benchmark_coverage(report, solver_names):
    """Require measured primary work from every solver, not just a return code."""
    runs = [run for problem in report['problems']
            if problem['role'] == 'primary' and problem['status'] == 'completed'
            for run in problem['runs']
            if run.get('execution', {}).get('kind') != 'repeated']
    for index, name in enumerate(solver_names, start=1):
        measured = [run for run in runs if run['solver_index'] == index
                    and (run['evaluations'] or 0) > 0]
        if not measured:
            raise RuntimeError(f'Action-test solver {name} has no recorded objective evaluations.')
        abnormal = sum(run['abnormal_termination'] is True for run in measured)
        # Budget exhaustion and deliberately difficult features can cause
        # abnormal termination. Retain those facts rather than rejecting them
        # or presenting them as normal solver returns.
        print(f'Benchmark coverage {name}: {len(measured)} evaluated runs; '
              f'{abnormal} abnormal terminations.')


def checked_benchmark(solvers, **options):
    """Run a benchmark with explicit solver capabilities and history coverage."""
    names = options['solver_names']
    if len(solvers) != len(names):
        raise ValueError('Each action-test solver must have one solver name.')
    for name in SKIPPED_SOLVERS:
        print(f'SKIP {name}: SciPy {scipy_version} does not provide this minimize method.')
    verify_solver_support(solvers, names)
    if not options.get('report_path'):
        report_root = Path('action_reports')
        report_root.mkdir(exist_ok=True)
        options['report_path'] = Path(tempfile.mkdtemp(dir=report_root, prefix='benchmark_')) / 'report.json'
    result = benchmark(solvers, **options)
    report = json.loads(Path(options['report_path']).read_text(encoding='utf-8'))
    verify_benchmark_coverage(report, names)
    return result


# Exported selections remain available to all three action scripts.
SOLVERS, SOLVER_NAMES, SKIPPED_SOLVERS = get_solver_configuration()

# Unconstrained solvers
UNCONSTRAINED_SOLVERS, UNCONSTRAINED_SOLVER_NAMES, _ = get_solver_configuration(True)
