import sys, os, io, importlib
from contextlib import redirect_stdout
import numpy as np

from ..utils import add_optiprofiler

add_optiprofiler()

from optiprofiler.problems import Problem, FeaturedProblem
from optiprofiler.features import Feature

def s2mpj_load(problem_name, **kwargs):
    """
    Load the S2MPJ problem.

    Parameters
    ----------
    problem_name : str
        Name of the problem in S2MPJ.
    **kwargs : dict
        Additional keyword arguments for the problem.

    Returns
    -------
    `Problem`
        An instance of the Problem class.
    """

    # Add the path of the problem to the system path.
    try:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        src_dir = os.path.join(current_dir, 'src')
        if src_dir not in sys.path:
            # Add the subfolder 'src' under the current directory to the system path (containing 's2mpjlib.py').
            sys.path.insert(0, src_dir)
        python_problems_dir = os.path.join(src_dir, 'python_problems')
        if python_problems_dir not in sys.path:
            # Add the subfolder 'python_problems' under the 'src' directory to the system path.
            sys.path.insert(0, python_problems_dir)
    except Exception as err:
        print(f"Failed to add the path of Python problems of S2MPJ to the system path: {err}")
        raise

    # Import the problem module dynamically.
    problem_module = importlib.import_module(f'python_problems.{problem_name}')

    # Get the class from the module.
    problem_class = getattr(problem_module, problem_name)

    # Create an instance of the class from the module.
    p = problem_class()

    # List of feasibility problems.
    feasibility_list = [
        'AIRCRFTA', 'ARGAUSS', 'ARGLALE', 'ARGLBLE'
        ]
    is_feasibility = problem_name in feasibility_list

    # Collect the problem parameters for transforming into a Problem instance.
    name = problem_name

    fun = lambda x: getfun(p, is_feasibility, x)
    grad = lambda x: getgrad(p, is_feasibility, x)
    hess = lambda x: gethess(p, is_feasibility, x)

    x0 = p.x0
    xl = p.xlower
    xu = p.xupper

    if not hasattr(p, 'lincons'):
        p.lincons = None
    if not hasattr(p, 'nle'):
        p.nle = 0
    if not hasattr(p, 'neq'):
        p.neq = 0
    if not hasattr(p, 'nge'):
        p.nge = 0

    buf = io.StringIO()
    with redirect_stdout(buf):
        try:
            cx, jx = p.cJHx(x0)[2]
            if hasattr(cx, 'toarray'):
                cx = cx.toarray()
            if hasattr(jx, 'toarray'):
                jx = jx.toarray()
            bx = jx @ x0 - cx
        except Exception:
            jx = None
            bx = None

    cx, jx = p.cJHx(x0)[2]
    if hasattr(cx, 'toarray'):
        cx = cx.toarray()
    if hasattr(jx, 'toarray'):
        jx = jx.toarray()
    bx = jx @ x0 - cx

    nonlincons = np.setdiff1d(np.arange(p.m), p.lincons)
    idx_le = np.arange(p.nle)
    idx_eq = np.arange(p.nle, p.nle + p.neq)
    idx_ge = np.arange(p.nle + p.neq, p.nle + p.neq + p.nge)
    idx_aeq = np.intersect1d(idx_eq, p.lincons)
    idx_aub_le = np.intersect1d(idx_le, p.lincons)
    idx_aub_ge = np.intersect1d(idx_ge, p.lincons)
    idx_cle = np.intersect1d(idx_le, nonlincons)
    idx_ceq = np.intersect1d(idx_eq, nonlincons)
    idx_cge = np.intersect1d(idx_ge, nonlincons)
    aeq = jx[idx_aeq, :]
    aub =np.vstack([jx[idx_aub_le, :], -jx[idx_aub_ge, :]])
    beq = bx[idx_aeq]
    bub = np.concatenate([bx[idx_aub_le], -bx[idx_aub_ge]])

    getidx = lambda y, idx: y[idx]
    ceq = lambda x: getidx(getcx(p, x), idx_ceq)
    cub = lambda x: np.concatenate([getidx(getcx(p, x), idx_cle),
                                    -getidx(getcx(p, x), idx_cge)])
    hceq = lambda x: getidx(getHx(p, x), idx_ceq)
    hcub = lambda x: np.hstack([getidx(getHx(p, x), idx_cle),
                                getidx(getHx(p, x), idx_cge)])

    getidx_mat = lambda y, idx: y[idx, :]
    jceq = lambda x: getidx_mat(getJx(p, x), idx_ceq)
    jcub = lambda x: np.vstack([getidx_mat(getJx(p, x), idx_cle),
                                -getidx_mat(getJx(p, x), idx_cge)])

    # Construct the Problem instance.
    problem = Problem(fun, x0, name=name, xl=xl, xu=xu, aub=aub, bub=bub, aeq=aeq, beq=beq, cub=cub, ceq=ceq, grad=grad, hess=hess, jcub=jcub, jceq=jceq, hcub=hcub, hceq=hceq)

    return problem

def s2mpj_select():
    pass

def getfun(p, is_feasibility, x):
    if is_feasibility:
        f = 0
    else:
        buf = io.StringIO()
        with redirect_stdout(buf):
            try:
                f = p.fx(x)
            except Exception:
                f = None
    return f

def getgrad(p, is_feasibility, x):
    if is_feasibility:
        g = np.zeros_like(x)
    else:
        buf = io.StringIO()
        with redirect_stdout(buf):
            try:
                _, g = p.fgx(x)
                g = g.toarray() if hasattr(g, 'toarray') else g
            except Exception:
                g = None
    return g

def gethess(p, is_feasibility, x):
    if is_feasibility:
        h = np.zeros((len(x), len(x)))
    else:
        buf = io.StringIO()
        with redirect_stdout(buf):
            try:
                _, _, h = p.fgHx(x)
                h = h.toarray() if hasattr(h, 'toarray') else h
            except Exception:
                h = None
    return h

def getcx(p, x):
    buf = io.StringIO()
    with redirect_stdout(buf):
        try:
            c = p.cx(x)
            c = c.toarray() if hasattr(c, 'toarray') else c
        except Exception:
            c = None
    return c

def getJx(p, x):
    buf = io.StringIO()
    with redirect_stdout(buf):
        try:
            _, j = p.cJHx(x)[2]
            j = j.toarray() if hasattr(j, 'toarray') else j
        except Exception:
            j = None
    return j

def getHx(p, x):
    # buf = io.StringIO()
    # with redirect_stdout(buf):
    try:
        _, _, h = p.cJHx(x)
        for i in range(len(h)):
            if hasattr(h[i], 'toarray'):
                h[i] = h[i].toarray()
    except Exception:
            h = None
    return h