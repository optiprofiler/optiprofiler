"""Example: wrap SciPy COBYQA for OptiProfiler.

OptiProfiler passes nonlinear inequality and equality callbacks separately:
``cub(x) <= 0`` and ``ceq(x) = 0``. SciPy's ``minimize`` interface expects
constraints to be supplied as ``Bounds``, ``LinearConstraint``, and
``NonlinearConstraint`` objects. This example shows the signature conversion
for the COBYQA method in SciPy: the wrapper translates OptiProfiler's
separate callbacks into SciPy constraint objects before calling ``minimize``.
"""

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint, minimize

from optiprofiler import benchmark


def scipy_cobyqa_wrapper(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, maxfev=200):
    """Solve an OptiProfiler problem with SciPy's COBYQA method."""
    constraints = []

    if bub.size > 0:
        # OptiProfiler gives aub @ x <= bub; SciPy represents it with
        # lower bound -inf and upper bound bub.
        constraints.append(LinearConstraint(aub, -np.inf, bub))
    if beq.size > 0:
        # Equality constraints use identical lower and upper bounds.
        constraints.append(LinearConstraint(aeq, beq, beq))

    c_ub_x0 = np.atleast_1d(cub(x0))
    if c_ub_x0.size > 0:
        # Convert cub(x) <= 0 to SciPy's object-based constraint interface.
        constraints.append(NonlinearConstraint(cub, -np.inf, np.zeros_like(c_ub_x0)))

    c_eq_x0 = np.atleast_1d(ceq(x0))
    if c_eq_x0.size > 0:
        # Convert ceq(x) = 0 by setting both bounds to zero.
        constraints.append(NonlinearConstraint(ceq, np.zeros_like(c_eq_x0), np.zeros_like(c_eq_x0)))

    result = minimize(
        fun,
        x0,
        method="COBYQA",
        bounds=Bounds(xl, xu),
        constraints=constraints,
        options={"maxfev": maxfev},
    )
    return result.x


def scipy_cobyqa_short(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq):
    """A short-budget COBYQA wrapper for comparison in this example."""
    return scipy_cobyqa_wrapper(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, maxfev=100)


def scipy_cobyqa_long(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq):
    """A longer-budget COBYQA wrapper for comparison in this example."""
    return scipy_cobyqa_wrapper(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, maxfev=200)


if __name__ == "__main__":
    scores = benchmark(
        [scipy_cobyqa_short, scipy_cobyqa_long],
        solver_names=["SciPy COBYQA short", "SciPy COBYQA long"],
        ptype="n",
        problem_names=["HS10", "HS11", "HS12"],
        mindim=2,
        maxdim=5,
        max_eval_factor=500,
        plibs=["s2mpj"],
        draw_hist_plots="none",
        n_jobs=1,
    )
    print(f"Solver scores: {scores[0]}")
