"""Solvers for action tests."""
import numpy as np
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint


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


# Exported solvers for action tests
SOLVERS = [scipy_cobyla, scipy_cobyqa, scipy_nelder_mead]
SOLVER_NAMES = ['COBYLA', 'COBYQA', 'Nelder-Mead']

# Unconstrained solvers
UNCONSTRAINED_SOLVERS = [scipy_nelder_mead, scipy_cobyqa]
UNCONSTRAINED_SOLVER_NAMES = ['Nelder-Mead', 'COBYQA']
