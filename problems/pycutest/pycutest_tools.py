import sys, os, io, re, importlib
from contextlib import redirect_stdout
import numpy as np
import pandas as pd

import builtins
builtins.np = np  # Make numpy available globally as 'np'

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

def pycutest_load(problem_name, **kwargs):
    """
    Load a problem from pycutest and return a Problem instance.

    Parameters
    ----------
    problem_name : str
        The name of the problem in pycutest to load.
    **kwargs : dict
        Additional keyword arguments (only for problems with available SIF parameters).
    
    Returns
    -------
    `Problem`
        An instance of the Problem class.
    """

    # Check if 'problem_name' has the pattern '_n_m' or 'n'. If it has, find the position of the pattern and return the dimension 'n' and the number of constraints 'm'.
    # Note that when 'm' is 0, the pattern is '_n' instead of '_n_0'.
    pattern = r'_(\d+)_(\d+)$|_(\d+)$'
    match = re.search(pattern, problem_name)
    name = problem_name
    if match:
        name = problem_name[:match.start()]
        if match.group(1) and match.group(2):
            dim = int(match.group(1))
            mcon = int(match.group(2))
        elif match.group(3):
            dim = int(match.group(3))
            mcon = 0
        else:
            raise ValueError(f"Invalid problem name format: {problem_name}")
    else:
        dim = None
        mcon = None

    # Initialize CUTEst problem
    p = pycutest.import_problem(name, destination=problem_name, sifParams=kwargs)

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

def pycutest_select(**kwargs):
    problem_names = pycutest.find_problems(**kwargs)
    return problem_names

# Define a function to collect SIF parameters
def pycutest_get_sif_params(problem_name):
    """
    Parse the printed output of pycutest.print_available_sif_params(problem_name)
    and extract available SIF parameters.

    Parameters
    ----------
    problem_name : str
        The name of the problem in pycutest.

    Returns
    -------
    para_names : list of str
        list of parameter names (e.g., ['A', 'N'])
    para_values : list of list
        list of possible values, each is a list (e.g., [[1, 2, 3], [10, 20, 40]])
    para_defaults : list
        list of default values (e.g., [2, 20])
    If the problem has no SIF parameters, all three lists are empty.
    """

    # Capture the printed output of pycutest.print_available_sif_params
    buf = io.StringIO()
    sys_stdout = sys.stdout
    sys.stdout = buf
    try:
        pycutest.print_available_sif_params(problem_name)
    finally:
        sys.stdout = sys_stdout

    lines = buf.getvalue().splitlines()
    buf.close()

    # Updated pattern — now we capture also the optional type (inside parentheses)
    # Example line: "N = 40 (int) [default]"
    # Exclude lines where parentheses contain 'for' (e.g., "(for SIF parameter N=1)")
    pattern = re.compile(
        r"^\s*([^\s=]+)\s*=\s*([^\s\(]+)"   # param name and value
        r"(?:\s*\((?!.*\bfor\b)([^)]*)\))?" # optional '(int)' or '(float)' group, negative lookahead to avoid 'for'
        r"\s*(\[default\])?\s*$"            # optional [default]
    )

    params = {}    # dict: param_name -> list of values
    defaults = {}  # dict: param_name -> default value

    for line in lines:
        m = pattern.search(line)
        if not m:
            continue

        key = m.group(1)                # parameter name
        val_raw = m.group(2)            # numeric/string value
        type_hint = m.group(3) or ""    # 'int', 'float', maybe empty
        default_flag = bool(m.group(4)) # whether [default] exists

        # Determine type conversion priority:
        # 1. explicit type hint in parentheses (best)
        # 2. fallback heuristic based on the string itself
        val = val_raw
        try:
            if "float" in type_hint.lower():
                val = float(val_raw)
            elif "int" in type_hint.lower():
                val = int(val_raw)
            else:
                # fallback heuristic if no type specified
                if "." in val_raw or "e" in val_raw.lower():
                    val = float(val_raw)
                else:
                    val = int(val_raw)
        except ValueError:
            # fallback to original string if conversion fails
            val = val_raw

        params.setdefault(key, []).append(val)
        if default_flag:
            defaults[key] = val

    # Return empty lists if no parameters found
    if not params:
        return [], [], []
    
    # Remove duplicate values while preserving order
    for key in params:
        seen = set()
        unique_values = []
        for val in params[key]:
            if val not in seen:
                seen.add(val)
                unique_values.append(val)
        params[key] = unique_values

    # Keep order of appearance
    para_names = list(params.keys())
    para_values = [params[k] for k in para_names]
    para_defaults = [defaults.get(k) for k in para_names]

    # Filter out parameters with only one value
    filtered_names = []
    filtered_values = []
    filtered_defaults = []
    
    for name, values, default in zip(para_names, para_values, para_defaults):
        if len(values) > 1:
            filtered_names.append(name)
            filtered_values.append(values)
            filtered_defaults.append(default)

    return filtered_names, filtered_values, filtered_defaults

def pycutest_clear_cache(problem_name, **kwargs):
    pycutest.clear_cache(problem_name, sifParams=kwargs)