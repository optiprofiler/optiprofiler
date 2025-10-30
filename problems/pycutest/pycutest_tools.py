import sys, os, io, re, importlib
from contextlib import redirect_stdout
import numpy as np
import pandas as pd

# Set the destination directory for pycutest cache.
current_dir = os.path.dirname(os.path.abspath(__file__))
cache_dir = os.path.join(current_dir, 'pycutest_cache')
os.makedirs(cache_dir, exist_ok=True)
if cache_dir not in sys.path:
    sys.path.append(cache_dir)
os.environ['PYCUTEST_CACHE'] = cache_dir
import pycutest

# Import Problem class from optiprofiler
from ..utils import add_optiprofiler
add_optiprofiler()
from optiprofiler.problems import Problem

def pycutest_load(problem_name):
    """
    Load a problem from pycutest and return a Problem instance.
    """

    
    name = problem_name

    # Initialize CUTEst problem
    p = pycutest.import_problem(name, destination=problem_name)

    fun = lambda x: p.obj(x)
    grad = lambda x: p.grad(x)
    hess = lambda x: p.ihess(x)

    # We replace 1.0e+20 (which represents infinity in pycutest) bounds with np.inf.
    p.bl = np.where(p.bl <= -1.0e+20, -np.inf, p.bl)
    p.bu = np.where(p.bu >= 1.0e+20, np.inf, p.bu)
    p.cl = np.where(p.cl <= -1.0e+20, -np.inf, p.cl) if p.m > 0 else p.cl
    p.cu = np.where(p.cu >= 1.0e+20, np.inf, p.cu) if p.m > 0 else p.cu

    x0 = p.x0
    xl = p.bl
    xu = p.bu

    mask_linear = getattr(p, "is_linear_cons", None)
    if mask_linear is None:
        mask_linear = np.zeros(p.m, dtype=bool)
    mask_eq = getattr(p, "is_eq_cons", None)
    if mask_eq is None:
        mask_eq = np.zeros(p.m, dtype=bool)
    mask_linear_eq = mask_linear & mask_eq
    mask_linear_ineq = mask_linear & ~mask_eq
    mask_nonlinear_eq = ~mask_linear & mask_eq
    mask_nonlinear_ineq = ~mask_linear & ~mask_eq

    buf = io.StringIO()
    with redirect_stdout(buf):
        try:
            cx, jx = p.cons(x0, gradient=True)
            bx = jx @ x0 - cx
        except Exception:
            jx = None
            bx = None

    # Create the linear equality constraints if any.
    aeq = jx[mask_linear_eq, :] if np.any(mask_linear_eq) else np.zeros((0, p.n))
    beq = bx[mask_linear_eq] if np.any(mask_linear_eq) else np.zeros(0)

    # Create the linear inequality constraints if any.
    idx_upper = np.where(mask_linear_ineq & (p.cu < np.inf))[0] if p.m > 0 else np.array([], dtype=int)
    idx_lower = np.where(mask_linear_ineq & (p.cl > -np.inf))[0] if p.m > 0 else np.array([], dtype=int)

    aub_upper = jx[idx_upper, :] if np.any(idx_upper) else np.zeros((0, p.n))
    bub_upper = p.cu[idx_upper] - bx[idx_upper] if np.any(idx_upper) else np.zeros(0)

    aub_lower = -jx[idx_lower, :] if np.any(idx_lower) else np.zeros((0, p.n))
    bub_lower = -(p.cl[idx_lower] - bx[idx_lower]) if np.any(idx_lower) else np.zeros(0)

    if idx_upper.size > 0 or idx_lower.size > 0:
        aub = np.vstack([aub_upper, aub_lower])
        bub = np.concatenate([bub_upper, bub_lower])
    else:
        aub = np.zeros((0, p.n))
        bub = np.zeros(0)

    # Handle nonlinear constraints.
    def _process_nonlinear_ineq(x, mode="value"):
        """Helper for nonlinear inequality constraints
        mode in {"value", "jacobian", "hessian"}
        """
        if not np.any(mask_nonlinear_ineq):
            if mode == "value":
                return np.zeros(0)
            elif mode == "jacobian":
                return np.zeros((0, p.n))
            elif mode == "hessian":
                return []

        ineq_idx = np.where(mask_nonlinear_ineq)[0]
        cu = np.asarray(p.cu)
        cl = np.asarray(p.cl)
        # Indices of upper and lower bounds for nonlinear inequalities.
        upper_idx = ineq_idx[cu[ineq_idx] < np.inf]
        lower_idx = ineq_idx[cl[ineq_idx] > -np.inf]

        if mode == "value":
            c_all = p.cons(x)
            upper_vals = c_all[upper_idx] - cu[upper_idx]
            lower_vals = -(c_all[lower_idx] - cl[lower_idx])
            return np.concatenate([upper_vals, lower_vals])
        elif mode == "jacobian":
            _, j_all = p.cons(x, gradient=True)
            j_upper = j_all[upper_idx, :]
            j_lower = -j_all[lower_idx, :]
            return np.vstack([j_upper, j_lower]) if (upper_idx.size or lower_idx.size) else np.zeros((0, p.n))
        elif mode == "hessian":
            hlist = []
            for i in upper_idx:
                H = p.ihess(x, cons_index=i)
                hlist.append(H)
            for i in lower_idx:
                H = -p.ihess(x, cons_index=i)
                hlist.append(H)
            return hlist

    # Create the nonlinear equality constraints.
    def ceq(x):
        if not np.any(mask_nonlinear_eq):
            return np.zeros(0)
        c_all = p.cons(x)
        return c_all[mask_nonlinear_eq]
    
    # Create the nonlinear inequality constraints.
    def cub(x):
        return _process_nonlinear_ineq(x, "value")
    
    # Create Jacobian functions for nonlinear constraints.
    def jceq(x):
        if not np.any(mask_nonlinear_eq):
            return np.zeros((0, p.n))
        _, j_all = p.cons(x, gradient=True)
        return j_all[mask_nonlinear_eq, :]
    
    def jcub(x):
        return _process_nonlinear_ineq(x, "jacobian")
    
    # Create Hessian functions for nonlinear constraints.
    def hceq(x):
        if not np.any(mask_nonlinear_eq):
            return []
        idx_eq = np.where(mask_nonlinear_eq)[0]
        hlist = []
        for i in idx_eq:
            H = p.ihess(x, cons_index=i)
            hlist.append(H)
        return hlist

    def hcub(x):
        return _process_nonlinear_ineq(x, "hessian")

    problem = Problem(fun, x0, name=name, xl=xl, xu=xu, aub=aub, bub=bub, aeq=aeq, beq=beq, cub=cub, ceq=ceq, grad=grad, hess=hess, jcub=jcub, jceq=jceq, hcub=hcub, hceq=hceq)

    return problem

