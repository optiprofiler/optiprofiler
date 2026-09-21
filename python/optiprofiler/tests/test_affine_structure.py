"""
Affine structure safeguard: one decision for bounds and linear constraints.

A feature that changes variables by ``x = A @ y + b`` (``linearly_transformed``
and ``custom`` with ``mod_affine``) must hand the solver the *same* problem in
the new coordinates. There are two representations of the bounds:

- the diagonal shortcut, valid when ``A`` is diagonal: the bounds stay bounds,
  scaled by ``diag(inv)``;
- the generic representation, valid for every invertible ``A``: the solver's
  bounds become infinite and every finite bound becomes a linear row of ``A``.

The defect these tests pin: the bounds were classified by looking at ``inv``
and the linear rows by looking at ``A``, both with an exact test. A diagonal
``A`` whose supplied inverse carries any off-diagonal entry, even at roundoff
level, made the bounds infinite (``inv`` is "not diagonal") while no bound row
was added (``A`` "is diagonal"). The bounds vanished from the posed problem,
silently, while the truth channel went on scoring the original bounds.

What is asserted here is the structure the solver receives (``xl``, ``xu``,
``aub``, ``bub``, ``aeq``, ``beq``), not merely that construction returns, and
that a point satisfies that structure exactly when the truth channel calls it
feasible.

The follow-up these tests also pin. The feasible set ``xl <= A @ y + b <= xu``
is a box exactly when ``A`` is diagonal, so the decision reads ``A`` exactly: an
entry of ``4e-16`` next to bounds of ``1e16`` moves the set by 4, and a
tolerance on the entries called it roundoff. The inverse cannot change the set
and decides nothing, except that the shortcut scales by ``diag(inv)`` and
therefore requires it to be the reciprocal of ``diag(A)`` to roundoff. One
transformation is produced per problem and seed, so that user code with a state
cannot hand the bounds one map and the linear constraints another. A supplied
inverse has to invert ``A`` from both sides, measured so that the units of
neither set of variables matter. And every finite quantity that is transported
(shifted bounds, right-hand sides, composed rows, the initial point) raises if
it leaves the floating-point range, instead of being posed as "no constraint".

The second follow-up. The range has two ends: a nonzero bound that falls below
the smallest normal number, or a transported entry whose terms all do, raises
as well (``1e-200 * [1e-200, 2e-200]`` was posed as the point ``[0, 0]``). The
initial point is verified where it is used, because no tolerance on the
matrices bounds an error at a point: ``inv @ (x0 - b)`` has to be mapped back
to ``x0`` by ``A`` to the rounding of that evaluation, or the equation is solved
from ``A``, or construction raises. The derivative methods of a single-feature
problem follow the change of variables by the chain rule (they returned the
derivatives of the original callbacks at the solver's point), and a composition
still provides none. The consistency of a supplied inverse allows the rounding
of its products and nothing more. And three decisions are pinned as they are:
what the deprecated conveniences keep (nothing), what ``mod_bounds`` replaces
(the original bounds in every representation), and how an integer beyond ``2**53`` is
converted (rounded, as every number is).
"""

import copy
import copyreg
import io
import pickle

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem
from optiprofiler.composition import AffineView, ComposedFeaturedProblem
from optiprofiler.opclasses import (_StageRuntime, _affine_is_diagonal, _checked_affine, _pulled_back,
                                    _restore_featured_problem)


EPS = np.finfo(float).eps
MAPPING = 'feasible_objective/1'
REFERENCE = {'merit': 0.0, 'kind': 'lower_bound', 'source': 'sum of squares', 'mapping': MAPPING}

# A diagonal change of variables with a negative entry (bounds must swap) and
# three different scales, and its exact inverse (all entries are powers of two).
D = np.array([2.0, -4.0, 0.5])
B = np.array([0.5, -0.5, 1.0])


def sphere(x):
    return float(np.sum(np.asarray(x, dtype=float) ** 2))


def bounded_problem():
    """Finite and infinite bounds, no linear constraints."""
    return Problem(sphere, [0.5, 0.5, 0.5], xl=[-1.0, -2.0, -np.inf], xu=[3.0, np.inf, 4.0], reference=REFERENCE)


def linear_problem():
    """The same bounds with one linear inequality and one linear equality."""
    return Problem(sphere, [0.5, 0.5, 0.5], xl=[-1.0, -2.0, -np.inf], xu=[3.0, np.inf, 4.0],
                   aub=[[1.0, 1.0, 0.0]], bub=[5.0], aeq=[[1.0, 0.0, -1.0]], beq=[0.25], reference=REFERENCE)


def fixed_variable_problem():
    """The first variable is fixed (``xl == xu``): the generic path must turn it into an equality row."""
    return Problem(sphere, [1.0, 0.5, 0.5], xl=[1.0, -2.0, -np.inf], xu=[1.0, np.inf, 4.0],
                   aub=[[1.0, 1.0, 0.0]], bub=[5.0], reference=REFERENCE)


PROBLEMS = {'bounded': bounded_problem, 'linear': linear_problem, 'fixed': fixed_variable_problem}


# ---------------------------------------------------------------- transforms (module level: they are pickled)

def exact_diagonal(rng, problem):
    return np.diag(D), B.copy(), np.diag(1.0 / D)


def roundoff_inverse(rng, problem):
    """Diagonal ``A``; the inverse carries off-diagonal entries of roundoff size (a fraction of one ulp of its diagonal)."""
    inv = np.diag(1.0 / D)
    inv[0, 1] = 1e-17
    inv[2, 0] = -0.5 * EPS * 0.5
    return np.diag(D), B.copy(), inv


def roundoff_matrix(rng, problem):
    """The mirror case: the inverse is exactly diagonal and ``A`` carries the roundoff."""
    A = np.diag(D)
    A[1, 0] = 0.5 * EPS * 2.0
    A[0, 2] = -1e-17
    return A, B.copy(), np.diag(1.0 / D)


def sloppy_inverse(rng, problem):
    """
    Diagonal ``A``; the inverse is accurate to about 1e-12, far above roundoff
    and far below the 1e-8 at which the pair would be called inconsistent.
    """
    inv = np.diag(1.0 / D)
    inv[0, 1] = 1e-12
    return np.diag(D), B.copy(), inv


def sloppy_matrix(rng, problem):
    """The mirror case: the inverse is exactly diagonal and ``A`` carries a coupling of 1e-12, which is posed."""
    A = np.diag(D)
    A[1, 0] = 1e-12
    return A, B.copy(), np.diag(1.0 / D)


def coupling_hidden_by_scale(rng, problem):
    """
    Diagonal ``A``; the inverse (diagonal 0.5, -0.25, 2) couples the rows of
    scales 0.25 and 2 by 2 eps. That is above ``n * eps`` times the smaller of
    the two scales (0.75 eps), although below ``n * eps`` times the largest
    entry of the matrix (6 eps): a tolerance scaled by a norm would overlook it.
    """
    inv = np.diag(1.0 / D)
    inv[1, 2] = 2.0 * EPS
    return np.diag(D), B.copy(), inv


DENSE = np.array([[2.0, 0.0, 1.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]])


def dense(rng, problem):
    return DENSE.copy(), B.copy(), np.linalg.inv(DENSE)


def singular(rng, problem):
    """A singular ``A`` with its pseudo-inverse, the most plausible stand-in for an inverse that does not exist."""
    A = np.array([[1.0, 2.0, 0.0], [2.0, 4.0, 0.0], [0.0, 0.0, 1.0]])
    return A, np.zeros(3), np.array([[0.04, 0.08, 0.0], [0.08, 0.16, 0.0], [0.0, 0.0, 1.0]])


def nan_in_matrix(rng, problem):
    A = np.diag(D)
    A[0, 1] = np.nan
    return A, B.copy(), np.diag(1.0 / D)


def inf_in_inverse(rng, problem):
    inv = np.diag(1.0 / D)
    inv[1, 2] = np.inf
    return np.diag(D), B.copy(), inv


def nan_in_shift(rng, problem):
    return np.diag(D), np.array([0.5, np.nan, 1.0]), np.diag(1.0 / D)


def wrong_matrix_shape(rng, problem):
    return np.eye(2), np.zeros(3), np.eye(3)


def wrong_inverse_shape(rng, problem):
    return np.eye(3), np.zeros(3), np.eye(2)


def wrong_shift_size(rng, problem):
    return np.eye(3), np.zeros(2), np.eye(3)


def complex_matrix(rng, problem):
    return np.diag(D).astype(complex), B.copy(), np.diag(1.0 / D)


def identity_as_inverse(rng, problem):
    return DENSE.copy(), B.copy(), np.eye(3)


def materially_wrong_inverse(rng, problem):
    inv = np.diag(1.0 / D)
    inv[0, 1] = 1e-3
    return np.diag(D), B.copy(), inv


ROTATION = np.array([[2.0, -2.0, 1.0], [1.0, 2.0, 2.0], [2.0, 1.0, -2.0]]) / 3.0


def scaled_rotation(ratio):
    """A rotation after scaling one NEW variable by ``ratio``, with its inverse (consistent to roundoff)."""
    scale = np.array([1.0, ratio, 1.0])
    return ROTATION * scale, np.zeros(3), (ROTATION / scale).T


def numerically_singular(rng, problem):
    """
    Consistent to roundoff, so the residual test cannot see it, and singular
    to working precision: no scaling of the original variables repairs 1e17.
    """
    return scaled_rotation(1e17)


def badly_scaled_rotation(rng, problem):
    """The same with 1e13: badly scaled, invertible, and to be accepted."""
    return scaled_rotation(1e13)


def extreme_diagonal(rng, problem):
    """An exact diagonal scaling over 300 orders of magnitude: its usual condition number is 1e300, its posed problem exact."""
    d = np.array([1e-150, -4.0, 1e150])
    return np.diag(d), B.copy(), np.diag(1.0 / d)


def overflowing_scale(rng, problem):
    """
    Exactly consistent and perfectly conditioned, yet the scaled bound
    ``2**1023 * 3`` overflows: a finite bound would become infinite.
    """
    return np.diag([2.0 ** -1023] * 3), np.zeros(3), np.diag([2.0 ** 1023] * 3)


def rounding_level_coupling(rng, problem):
    """
    ``A`` differs from the identity by one entry of 4e-16, which is below
    ``n * eps``, and the inverse is exact. With bounds of 1e16 that entry moves
    the feasible set by 4: it is a coupling, however small it looks.
    """
    return np.array([[1.0, 4e-16], [0.0, 1.0]]), np.zeros(2), np.array([[1.0, -4e-16], [0.0, 1.0]])


def inaccurate_diagonal_inverse(rng, problem):
    """
    Exactly diagonal, and consistent to 1e-9, which the contract accepts (1e-8).
    Bounds scaled by this ``diag(inv)`` would be off by 1e-9 of their size.
    """
    return np.diag(D), B.copy(), np.diag((1.0 + 1e-9) / D)


def one_sided_inverse(rng, problem):
    """``A @ inv`` is the identity to 1e-16, and ``inv @ A`` misses it by 1: ``inv @ (A @ y)`` is not ``y``."""
    inv = np.diag([1e8, 1e-8, 1.0])
    inv[0, 1] = 1e-8
    return np.diag([1e-8, 1e8, 1.0]), np.zeros(3), inv


def one_sided_inverse_mirror(rng, problem):
    """The same with the large and the small scale exchanged, so that the stray entry sits below the diagonal."""
    inv = np.diag([1e-8, 1e8, 1.0])
    inv[1, 0] = 1e-8
    return np.diag([1e8, 1e-8, 1.0]), np.zeros(3), inv


def generic_rotation():
    """
    A rotation that is orthogonal to roundoff only. The entries of ``ROTATION``
    are ``t`` and ``2 * t``, so that its products cancel exactly in plain
    floating-point arithmetic and not in fused arithmetic: whether the base
    accepted a scaled ``ROTATION`` depended on the kernel that multiplied it.
    """
    c1, s1, c2, s2 = np.cos(0.3), np.sin(0.3), np.cos(0.5), np.sin(0.5)
    return np.array([[c1, -s1, 0.0], [s1, c1, 0.0], [0.0, 0.0, 1.0]]) @ np.array([[1.0, 0.0, 0.0], [0.0, c2, -s2], [0.0, s2, c2]])


def row_scaled_rotation(rng, problem):
    """
    The mirror of ``badly_scaled_rotation``: one ORIGINAL variable is in units
    of 1e13 (a row of ``A``). Consistent to roundoff and well conditioned in
    the sense that is checked, so it is to be accepted as well.
    """
    rotation, scale = generic_rotation(), np.array([1.0, 1e13, 1.0])
    return rotation * scale[:, None], np.zeros(3), rotation.T / scale


def column_scaled_rotation(rng, problem):
    """The same rotation with one NEW variable in units of 1e13 (a column of ``A``)."""
    rotation, scale = generic_rotation(), np.array([1.0, 1e13, 1.0])
    return rotation * scale, np.zeros(3), (rotation / scale).T


def shifted_out_of_range(rng, problem):
    return np.eye(3), np.full(3, -1e308), np.eye(3)


def dense_shifted_out_of_range(rng, problem):
    return DENSE.copy(), np.full(3, -1e308), np.linalg.inv(DENSE)


def shifted_the_other_way(rng, problem):
    return np.eye(3), np.full(3, 1e308), np.eye(3)


def huge_scaling(rng, problem):
    """Exactly consistent and perfectly conditioned; a coefficient of 1e200 times 1e200 overflows."""
    return np.diag([1e200, 1.0, 1.0]), np.zeros(3), np.diag([1e-200, 1.0, 1.0])


class Alternating:
    """User code with a state: a valid diagonal transform and a valid dense one in turn."""

    def __init__(self, calls=0):
        self.calls = calls

    def __call__(self, rng, problem):
        self.calls += 1
        return exact_diagonal(rng, problem) if self.calls % 2 == 1 else dense(rng, problem)


class Drifting:
    """
    Another valid rotation on every call, as when user code draws from a
    global generator instead of the stream it is handed.
    """

    def __init__(self):
        self.calls = 0

    def __call__(self, rng, problem):
        self.calls += 1
        return self.rotation(self.calls)

    @staticmethod
    def rotation(call):
        c, s = np.cos(0.3 * call), np.sin(0.3 * call)
        Q = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        return Q, B.copy(), Q.T


class Counting:
    """A pure transform that counts how often it is asked."""

    def __init__(self):
        self.calls = 0

    def __call__(self, rng, problem):
        self.calls += 1
        return dense(rng, problem)


def custom_stage(transform):
    return {'name': 'custom', 'options': {'mod_affine': transform}}


# ---------------------------------------------------------------- what the solver receives

def affine_of(featured):
    """The map ``x = A @ y + b`` from solver coordinates to original coordinates, composed over the stages."""
    n = featured.n
    if isinstance(featured, ComposedFeaturedProblem):
        A, b = np.eye(n), np.zeros(n)
        for view in featured._views:
            if isinstance(view, AffineView):
                A, b = A @ view._A, A @ view._b + b
        return A, b
    A, b, _ = featured._runtime.modifier_affine(featured._seed, featured._problem)
    return np.asarray(A, dtype=float), np.asarray(b, dtype=float)


def posed_violation(featured, y):
    """Largest violation of the structure handed to the solver: bounds, linear inequalities, linear equalities."""
    violation = 0.0
    lower, upper = np.isfinite(featured.xl), np.isfinite(featured.xu)
    if lower.any():
        violation = max(violation, np.max(featured.xl[lower] - y[lower]))
    if upper.any():
        violation = max(violation, np.max(y[upper] - featured.xu[upper]))
    if featured.aub.size:
        violation = max(violation, np.max(featured.aub @ y - featured.bub))
    if featured.aeq.size:
        violation = max(violation, np.max(np.abs(featured.aeq @ y - featured.beq)))
    return violation


def assert_posed_problem_is_the_scored_problem(featured, problem, seed=0):
    """
    A point satisfies the solver-facing structure exactly when the truth
    channel (the original problem at ``A @ y + b``) calls it feasible. A
    dropped bound breaks this at once: far outside the original box the posed
    structure reports no violation and the truth reports a large one.
    """
    A, b = affine_of(featured)
    rng = np.random.default_rng(seed)
    agree = {'feasible': 0, 'infeasible': 0}
    for x in rng.uniform(-6.0, 8.0, size=(400, problem.n)):
        if problem.m_linear_eq:  # sample on the equality so that feasible points exist
            x = x - problem.aeq.T @ np.linalg.solve(problem.aeq @ problem.aeq.T, problem.aeq @ x - problem.beq)
        if np.any(problem.xl == problem.xu):  # and on the fixed variables
            fixed = problem.xl == problem.xu
            x[fixed] = problem.xl[fixed]
        y = np.linalg.solve(A, x - b)
        truth, posed = featured.maxcv(y), posed_violation(featured, y)
        if 1e-9 < truth < 1e-3 or 1e-9 < posed < 1e-3:
            continue  # too close to a boundary to classify in floating point
        assert (truth <= 1e-9) == (posed <= 1e-9), (x, truth, posed)
        agree['feasible' if truth <= 1e-9 else 'infeasible'] += 1
    assert agree['feasible'] > 20 and agree['infeasible'] > 20, agree  # the sample covers both sides


def expected_diagonal_structure(problem, A, b, inv):
    """The diagonal shortcut: scaled bounds (swapped where the scale is negative) and linear maps composed with A."""
    scale = np.diagonal(inv)
    lower, upper = scale * (problem.xl - b), scale * (problem.xu - b)
    return {'xl': np.minimum(lower, upper), 'xu': np.maximum(lower, upper),
            'aub': problem.aub @ A if problem.aub.size else problem.aub, 'bub': problem.bub - problem.aub @ b,
            'aeq': problem.aeq @ A if problem.aeq.size else problem.aeq, 'beq': problem.beq - problem.aeq @ b}


def expected_generic_structure(problem, A, b):
    """The generic representation: no solver bounds; every finite bound is a row of A (an equality row if fixed)."""
    fixed = problem.xl == problem.xu
    upper, lower = np.isfinite(problem.xu) & ~fixed, np.isfinite(problem.xl) & ~fixed
    aub = [A[upper], -A[lower]] + ([problem.aub @ A] if problem.aub.size else [])
    bub = [problem.xu[upper] - b[upper], -(problem.xl[lower] - b[lower])] + \
        ([problem.bub - problem.aub @ b] if problem.aub.size else [])
    aeq = [A[fixed]] + ([problem.aeq @ A] if problem.aeq.size else [])
    beq = [problem.xu[fixed] - b[fixed]] + ([problem.beq - problem.aeq @ b] if problem.aeq.size else [])
    return {'xl': np.full(problem.n, -np.inf), 'xu': np.full(problem.n, np.inf),
            'aub': np.vstack(aub), 'bub': np.concatenate(bub), 'aeq': np.vstack(aeq), 'beq': np.concatenate(beq)}


def assert_structure(featured, expected):
    for key, value in expected.items():
        actual = getattr(featured, key)
        if value.size == 0:
            assert actual.size == 0, key
        else:
            np.testing.assert_array_equal(actual, value, err_msg=key)


def assert_no_bound_is_lost_or_doubled(featured, problem):
    """Every finite original bound is posed exactly once: as a solver bound or as a linear row, never neither or both."""
    fixed = problem.xl == problem.xu
    finite = int(np.isfinite(problem.xl).sum() + np.isfinite(problem.xu).sum())
    as_bounds = int(np.isfinite(featured.xl).sum() + np.isfinite(featured.xu).sum())
    as_rows = (featured.aub.shape[0] if featured.aub.size else 0) - problem.m_linear_ub \
        + 2 * ((featured.aeq.shape[0] if featured.aeq.size else 0) - problem.m_linear_eq)
    assert as_bounds + as_rows == finite, (as_bounds, as_rows, finite)
    assert as_bounds in (0, finite), 'bounds must be all solver bounds or all linear rows'
    assert int(fixed.sum()) * 2 <= finite


# ---------------------------------------------------------------- custom with mod_affine

class TestCustomAffine:

    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    def test_exact_diagonal_keeps_bounds_as_bounds(self, name):
        problem = PROBLEMS[name]()
        featured = FeaturedProblem(problem, Feature('custom', mod_affine=exact_diagonal), 10, 3)
        assert_structure(featured, expected_diagonal_structure(problem, np.diag(D), B, np.diag(1.0 / D)))
        # The negative scale swapped the second pair of bounds: (-2, inf) became (-inf, 0.375).
        np.testing.assert_array_equal(featured.xl, [-0.75, -np.inf, -np.inf] if name != 'fixed' else [0.25, -np.inf, -np.inf])
        np.testing.assert_array_equal(featured.xu[1:], [0.375, 6.0])
        assert_no_bound_is_lost_or_doubled(featured, problem)
        assert_posed_problem_is_the_scored_problem(featured, problem)
        np.testing.assert_allclose(np.diag(D) @ featured.x0 + B, problem.x0, atol=1e-15)

    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    @pytest.mark.parametrize('transform', [roundoff_inverse, sloppy_inverse, coupling_hidden_by_scale],
                             ids=lambda f: f.__name__)
    def test_stray_entries_of_the_inverse_keep_bounds_as_bounds(self, transform, name):
        # The original defect: with a stray entry in the inverse the bounds
        # became infinite and no bound row was added, so the bounds were gone.
        # ``A`` is exactly diagonal in all three, so the feasible set is a box
        # whatever the inverse carries off its diagonal, at roundoff level or
        # above it: the inverse cannot change the set and decides nothing.
        problem = PROBLEMS[name]()
        featured = FeaturedProblem(problem, Feature('custom', mod_affine=transform), 10, 3)
        A, b, inv = transform(None, problem)
        assert_structure(featured, expected_diagonal_structure(problem, A, b, inv))
        # Every finite bound is still a finite solver bound, and no bound row was added.
        assert np.isfinite(featured.xl).sum() + np.isfinite(featured.xu).sum() == \
            np.isfinite(problem.xl).sum() + np.isfinite(problem.xu).sum()
        assert_no_bound_is_lost_or_doubled(featured, problem)
        assert_posed_problem_is_the_scored_problem(featured, problem)
        # The reference is retained, and rightly so: the posed problem is the scored one.
        assert featured.reference == problem.reference

    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    @pytest.mark.parametrize('transform', [roundoff_matrix, sloppy_matrix], ids=lambda f: f.__name__)
    def test_any_coupling_in_the_matrix_takes_the_generic_path(self, transform, name):
        # An off-diagonal entry of ``A`` couples two variables, so the feasible
        # set is no box, and how far it is from one depends on the size of the
        # other variable, which no tolerance on the entry knows. It is posed,
        # exactly, by the representation that needs only ``A``. On the base a
        # coupling in the matrix posed every bound twice.
        problem = PROBLEMS[name]()
        featured = FeaturedProblem(problem, Feature('custom', mod_affine=transform), 10, 3)
        A, b, _ = transform(None, problem)
        assert_structure(featured, expected_generic_structure(problem, A, b))
        assert_no_bound_is_lost_or_doubled(featured, problem)
        assert_posed_problem_is_the_scored_problem(featured, problem)
        assert featured.reference == problem.reference

    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    def test_badly_scaled_transformations_keep_working(self, name):
        # Scaling is not singularity. The condition number that is checked
        # does not change with the units of the original variables, so an
        # exact diagonal scaling passes whatever its entries, and a rotation
        # of badly scaled variables passes while it is invertible in double
        # precision. Structures are compared exactly; sampling points at
        # these scales would measure roundoff, not structure.
        problem = PROBLEMS[name]()
        featured = FeaturedProblem(problem, Feature('custom', mod_affine=extreme_diagonal), 10, 3)
        assert_structure(featured, expected_diagonal_structure(problem, *extreme_diagonal(None, problem)))
        assert_no_bound_is_lost_or_doubled(featured, problem)
        featured = FeaturedProblem(problem, Feature('custom', mod_affine=badly_scaled_rotation), 10, 3)
        A, b, _ = badly_scaled_rotation(None, problem)
        assert_structure(featured, expected_generic_structure(problem, A, b))
        assert_no_bound_is_lost_or_doubled(featured, problem)
        # The units of the ORIGINAL variables matter as little as those of the
        # new ones. ``A @ inv`` misses the identity by 5e-4 here, roundoff
        # times the ratio of the units, and the base called that inconsistent
        # although the condition number it checks is 2.8.
        featured = FeaturedProblem(problem, Feature('custom', mod_affine=row_scaled_rotation), 10, 3)
        A, b, inv = row_scaled_rotation(None, problem)
        assert np.linalg.norm(A @ inv - np.eye(3)) > 1e-8 * 3
        assert_structure(featured, expected_generic_structure(problem, A, b))
        assert_no_bound_is_lost_or_doubled(featured, problem)

    @pytest.mark.parametrize('composed', [False, True], ids=['single', 'composed'])
    def test_a_problem_without_variables_is_still_built(self, composed):
        # Older NumPy cannot take the 1-norm or the infinity norm of an empty
        # matrix; the validation must not be what breaks the empty problem.
        def empty(rng, problem):
            return np.zeros((0, 0)), np.zeros(0), np.zeros((0, 0))

        stages = [{'name': 'custom', 'options': {'mod_affine': empty}}] + (['noisy'] if composed else [])
        featured = FeaturedProblem(Problem(sphere, np.zeros(0)), Feature(stages), 5, 3)
        assert featured.n == 0 and featured.xl.size == 0 and featured.aub.size == 0 and featured.aeq.size == 0

    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    def test_dense_transform_keeps_every_bound_as_a_linear_row(self, name):
        problem = PROBLEMS[name]()
        featured = FeaturedProblem(problem, Feature('custom', mod_affine=dense), 10, 3)
        assert_structure(featured, expected_generic_structure(problem, DENSE, B))
        assert not np.isfinite(featured.xl).any() and not np.isfinite(featured.xu).any()
        assert_no_bound_is_lost_or_doubled(featured, problem)
        assert_posed_problem_is_the_scored_problem(featured, problem)

    @pytest.mark.parametrize('transform, message', [
        (singular, 'not an identity matrix'),
        (identity_as_inverse, 'not an identity matrix'),
        (materially_wrong_inverse, 'not an identity matrix'),
        (nan_in_matrix, 'must be finite'),
        (inf_in_inverse, 'must be finite'),
        (nan_in_shift, 'must be finite'),
        (wrong_matrix_shape, r'shape \(3, 3\)'),
        (wrong_inverse_shape, r'shape \(3, 3\)'),
        (wrong_shift_size, 'size 3'),
        (complex_matrix, 'real'),
        (numerically_singular, 'numerically singular'),
        (one_sided_inverse, 'not an identity matrix'),
        (one_sided_inverse_mirror, 'not an identity matrix'),
    ], ids=lambda value: value.__name__ if callable(value) else None)
    @pytest.mark.parametrize('composed', [False, True], ids=['single', 'composed'])
    def test_invalid_transforms_fail_closed(self, transform, message, composed):
        # No featured problem exists afterwards, so nothing can be posed,
        # scored or claimed. NaN used to pass: `norm(...) > tolerance` is
        # false for NaN.
        stages = [custom_stage(transform)] + (['noisy'] if composed else [])
        for make in PROBLEMS.values():
            with pytest.raises(ValueError, match=message):
                FeaturedProblem(make(), Feature(stages), 10, 3)

    def test_a_finite_bound_never_becomes_infinite_by_overflow(self):
        # Every check on the matrices passes (the pair is exactly consistent
        # and perfectly conditioned); only the scaled bound overflows. Letting
        # it through would drop the bound without a word.
        with pytest.raises(ValueError, match='finite bound'):
            FeaturedProblem(bounded_problem(), Feature('custom', mod_affine=overflowing_scale), 10, 3)

    @pytest.mark.parametrize('composed', [False, True], ids=['single', 'composed'])
    def test_supplied_linear_rows_cannot_silently_replace_the_bound_rows(self, composed):
        # A supplied linear modifier replaces the linear constraints verbatim
        # (an established rule). Under a dense map the framework poses the
        # bounds as exactly those constraints, so the combination had no place
        # left for the bounds and dropped them. It now raises, and only when a
        # bound is really at stake.
        def rows(rng, problem):
            return np.array([[1.0, 1.0, 0.0]]), np.array([5.0])

        def equalities(rng, problem):
            return np.array([[1.0, 0.0, -1.0]]), np.array([0.25])

        def box(rng, problem):
            return np.full(3, -9.0), np.full(3, 9.0)

        def build(problem, **options):
            return FeaturedProblem(problem, Feature([{'name': 'custom', 'options': options}] + (['noisy'] if composed else [])), 10, 3)

        with pytest.raises(ValueError, match='would be dropped silently'):
            build(bounded_problem(), mod_affine=dense, mod_linear_ub=rows)
        with pytest.raises(ValueError, match='would be dropped silently'):
            build(fixed_variable_problem(), mod_affine=dense, mod_linear_eq=equalities)
        # The same decision is read here: a diagonal map whose inverse is not
        # good enough for the shortcut poses the bounds as rows as well.
        with pytest.raises(ValueError, match='would be dropped silently'):
            build(bounded_problem(), mod_affine=inaccurate_diagonal_inverse, mod_linear_ub=rows)

        # Nothing fixed, so replacing the equalities loses nothing: the bounds are inequality rows.
        featured = build(bounded_problem(), mod_affine=dense, mod_linear_eq=equalities)
        np.testing.assert_array_equal(featured.aub, expected_generic_structure(bounded_problem(), DENSE, B)['aub'])
        np.testing.assert_array_equal(featured.aeq, [[1.0, 0.0, -1.0]])
        # No finite bound, so nothing is at stake.
        unbounded = Problem(sphere, [0.5, 0.5, 0.5], aub=[[1.0, 1.0, 0.0]], bub=[5.0])
        np.testing.assert_array_equal(build(unbounded, mod_affine=dense, mod_linear_ub=rows).aub, [[1.0, 1.0, 0.0]])
        # The user takes over the bounds as well: both are theirs, verbatim.
        featured = build(bounded_problem(), mod_affine=dense, mod_linear_ub=rows, mod_bounds=box)
        np.testing.assert_array_equal(featured.xl, [-9.0, -9.0, -9.0])
        np.testing.assert_array_equal(featured.aub, [[1.0, 1.0, 0.0]])
        # A diagonal map keeps the bounds as bounds, so supplied rows replace nothing of theirs.
        featured = build(bounded_problem(), mod_affine=exact_diagonal, mod_linear_ub=rows)
        assert np.isfinite(featured.xl).sum() + np.isfinite(featured.xu).sum() == 4
        np.testing.assert_array_equal(featured.aub, [[1.0, 1.0, 0.0]])

    def test_structure_does_not_depend_on_which_modifier_is_asked_first(self):
        # One decision, read by all three modifiers: asking in any order, any
        # number of times, gives the same representation.
        problem = linear_problem()
        for transform in (exact_diagonal, roundoff_inverse, roundoff_matrix, sloppy_inverse, dense):
            runtime = _StageRuntime('custom', {'mod_affine': transform})
            first = [runtime.modifier_linear_eq(3, problem), runtime.modifier_linear_ub(3, problem),
                     runtime.modifier_bounds(3, problem)]
            second = [runtime.modifier_bounds(3, problem), runtime.modifier_linear_ub(3, problem),
                      runtime.modifier_linear_eq(3, problem)]
            for left, right in zip(first, reversed(second)):
                for a, b in zip(left, right):
                    np.testing.assert_array_equal(a, b)
            bounds_kept = np.isfinite(second[0][0]).any() or np.isfinite(second[0][1]).any()
            rows_added = second[1][0].shape[0] > problem.m_linear_ub
            assert bounds_kept != rows_added, transform.__name__


# ---------------------------------------------------------------- linearly_transformed

def perturb_linearly_transformed(monkeypatch, perturb):
    """
    The framework builds both matrices of ``linearly_transformed`` itself, so
    roundoff in its inverse cannot be requested through an option. The pair is
    perturbed where it is produced; the classification code under test is
    untouched.
    """
    original = _StageRuntime.modifier_affine

    def modifier_affine(self, seed, problem):
        A, b, inv = original(self, seed, problem)
        if self._name == 'linearly_transformed':
            A, inv = perturb(np.array(A), np.array(inv))
        return A, b, inv

    monkeypatch.setattr(_StageRuntime, 'modifier_affine', modifier_affine)


class TestLinearlyTransformed:

    UNROTATED = dict(rotated=False, condition_factor=4)

    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    def test_unrotated_is_the_exact_diagonal_case(self, name):
        problem = PROBLEMS[name]()
        featured = FeaturedProblem(problem, Feature('linearly_transformed', **self.UNROTATED), 10, 3)
        A, b, inv = featured._runtime.modifier_affine(3, problem)
        assert np.count_nonzero(A - np.diag(np.diagonal(A))) == 0 and np.linalg.cond(A) > 5.0
        assert_structure(featured, expected_diagonal_structure(problem, A, b, inv))
        assert_no_bound_is_lost_or_doubled(featured, problem)
        assert_posed_problem_is_the_scored_problem(featured, problem)

    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    @pytest.mark.parametrize('size', [0.25 * EPS, 1e-12], ids=['roundoff', 'sloppy'])
    def test_a_stray_entry_of_the_inverse_keeps_bounds_as_bounds(self, monkeypatch, size, name):
        # The matrix stays exactly diagonal, so the set is a box whatever the
        # inverse carries off its diagonal.
        def perturb(A, inv):
            inv[0, 1] = size * min(abs(inv[0, 0]), abs(inv[1, 1]))
            return A, inv

        perturb_linearly_transformed(monkeypatch, perturb)
        problem = PROBLEMS[name]()
        featured = FeaturedProblem(problem, Feature('linearly_transformed', **self.UNROTATED), 10, 3)
        A, b, inv = featured._runtime.modifier_affine(3, problem)
        assert inv[0, 1] != 0.0  # the perturbation is in place
        assert_structure(featured, expected_diagonal_structure(problem, A, b, inv))
        assert_no_bound_is_lost_or_doubled(featured, problem)
        assert_posed_problem_is_the_scored_problem(featured, problem)
        assert featured.reference == problem.reference

    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    def test_roundoff_in_the_matrix_takes_the_generic_path(self, monkeypatch, name):
        # A quarter of a unit in the last place of the diagonal is still a
        # coupling of two variables: it is posed, not ignored.
        def perturb(A, inv):
            A[0, 1] = 0.25 * EPS * min(abs(A[0, 0]), abs(A[1, 1]))
            return A, inv

        perturb_linearly_transformed(monkeypatch, perturb)
        problem = PROBLEMS[name]()
        featured = FeaturedProblem(problem, Feature('linearly_transformed', **self.UNROTATED), 10, 3)
        A, b, _ = featured._runtime.modifier_affine(3, problem)
        assert A[0, 1] != 0.0  # the perturbation is in place
        assert_structure(featured, expected_generic_structure(problem, A, b))
        assert_no_bound_is_lost_or_doubled(featured, problem)
        assert_posed_problem_is_the_scored_problem(featured, problem)
        assert featured.reference == problem.reference

    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    def test_rotated_is_the_dense_case(self, name):
        problem = PROBLEMS[name]()
        featured = FeaturedProblem(problem, Feature('linearly_transformed', rotated=True, condition_factor=4), 10, 3)
        A, b, _ = featured._runtime.modifier_affine(3, problem)
        assert np.count_nonzero(A - np.diag(np.diagonal(A))) > 0
        assert_structure(featured, expected_generic_structure(problem, A, b))
        assert_no_bound_is_lost_or_doubled(featured, problem)
        assert_posed_problem_is_the_scored_problem(featured, problem)

    @pytest.mark.parametrize('rotated', [False, True])
    def test_unusable_condition_factors_fail_closed(self, rotated):
        # 2**(+-1937) overflows. NumPy warns about the overflow it is asked to
        # produce here; the point of the test is that the framework then
        # refuses the result.
        with np.errstate(over='ignore', invalid='ignore'), pytest.raises(ValueError, match='must be finite'):
            FeaturedProblem(linear_problem(), Feature('linearly_transformed', rotated=rotated, condition_factor=1e7), 10, 3)

    def test_extreme_but_exact_scaling_keeps_working(self):
        # cond(A) = 2**95, unrotated. Every entry is exact and so is the posed
        # problem: a large condition number alone is no reason to refuse.
        problem = linear_problem()
        feature = Feature('linearly_transformed', rotated=False, condition_factor=6000.0)
        featured = FeaturedProblem(problem, feature, 10, 3)
        A, b, inv = featured._runtime.modifier_affine(featured._seed, featured._problem)
        assert np.linalg.cond(A) > 1e28
        assert_structure(featured, expected_diagonal_structure(problem, A, b, inv))
        assert_no_bound_is_lost_or_doubled(featured, problem)

    @pytest.mark.parametrize('condition_factor', [0, 2, 5, 40])
    @pytest.mark.parametrize('rotated', [False, True])
    def test_ordinary_condition_factors_keep_working(self, rotated, condition_factor):
        # The framework's own inverse is not held to the residual test of a
        # supplied one: a rotation with a large condition number is accurate
        # only to roundoff times that number, which is not an inconsistency.
        problem = linear_problem()
        featured = FeaturedProblem(problem, Feature('linearly_transformed', rotated=rotated,
                                                    condition_factor=condition_factor), 10, 3)
        assert_no_bound_is_lost_or_doubled(featured, problem)
        assert_posed_problem_is_the_scored_problem(featured, problem)

    def test_the_frameworks_own_inverse_is_not_held_to_the_residual_test(self):
        # cond(A) = 2**40 for n = 3. Both factors are built directly (nothing
        # is inverted), so each is accurate to roundoff, but their product
        # misses the identity by roundoff times the condition number: more
        # than a supplied inverse is allowed, and no inconsistency. The posed
        # rows need A only.
        problem = linear_problem()
        feature = Feature('linearly_transformed', rotated=True, condition_factor=2.0 * 40.0 ** 2 / 3.0)
        featured = FeaturedProblem(problem, feature, 10, 3)
        A, b, inv = featured._runtime.modifier_affine(featured._seed, featured._problem)
        assert np.linalg.norm(A @ inv - np.eye(3)) > 1e-8 * 3
        assert_structure(featured, expected_generic_structure(problem, A, b))
        assert_no_bound_is_lost_or_doubled(featured, problem)


# ---------------------------------------------------------------- composition, load and save

class TestCompositionAndPersistence:

    @pytest.mark.parametrize('stages, bounds_kept', [
        ([custom_stage(roundoff_inverse), 'noisy'], True),
        (['perturbed_x0', custom_stage(roundoff_matrix), 'truncated'], False),  # a coupling in A is posed
        ([{'name': 'linearly_transformed', 'options': {'rotated': False, 'condition_factor': 4}},
          custom_stage(roundoff_inverse)], True),
        (['permuted', custom_stage(exact_diagonal)], True),
        ([custom_stage(sloppy_inverse), 'noisy'], True),  # A is exactly diagonal: the set is a box
        ([custom_stage(dense), {'name': 'linearly_transformed', 'options': {'rotated': False}}], False),
    ], ids=['roundoff-inverse+noisy', 'x0+roundoff-matrix+truncated', 'scaling+roundoff-inverse',
            'permuted+exact', 'sloppy+noisy', 'dense+scaling'])
    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    def test_compositions_pose_the_scored_problem(self, stages, bounds_kept, name):
        problem = PROBLEMS[name]()
        featured = FeaturedProblem(problem, Feature(stages), 10, 3)
        assert isinstance(featured, ComposedFeaturedProblem)
        assert bool(np.isfinite(featured.xl).any() or np.isfinite(featured.xu).any()) is bounds_kept
        assert_posed_problem_is_the_scored_problem(featured, problem)
        assert featured.reference == problem.reference  # every stage here is a safe one

    def test_composed_problem_survives_pickling_with_its_structure(self):
        problem = linear_problem()
        for transform in (exact_diagonal, roundoff_inverse, sloppy_inverse, dense):
            featured = FeaturedProblem(problem, Feature([custom_stage(transform), 'noisy']), 10, 3)
            restored = pickle.loads(pickle.dumps(featured))
            for key in ('x0', 'xl', 'xu', 'aub', 'bub', 'aeq', 'beq'):
                np.testing.assert_array_equal(getattr(restored, key), getattr(featured, key), err_msg=key)
            assert_posed_problem_is_the_scored_problem(restored, problem)

    def test_single_feature_survives_pickling_without_resampling_affine(self):
        """A saved single-stage trial must keep its validated transform and history."""
        problem = linear_problem()
        callback = Drifting()
        featured = FeaturedProblem(problem, Feature('custom', mod_affine=callback), 10, 3)
        assert callback.calls == 1
        featured.fun(np.array([0.2, 0.3, 0.4]))
        restored = pickle.loads(pickle.dumps(featured))
        assert isinstance(restored, FeaturedProblem)
        assert callback.calls == 1  # unpickling must not call the user's callback
        restored_callback = restored._runtime._options['mod_affine']
        assert restored_callback.calls == 1
        for key in ('x0', 'xl', 'xu', 'aub', 'bub', 'aeq', 'beq', 'fun_hist', 'maxcv_hist'):
            np.testing.assert_array_equal(getattr(restored, key), getattr(featured, key), err_msg=key)
        restored_affine = restored._runtime.modifier_affine(restored._seed, restored._problem)
        original_affine = featured._runtime.modifier_affine(featured._seed, featured._problem)
        for actual, expected in zip(restored_affine, original_affine):
            np.testing.assert_array_equal(actual, expected)
            assert not actual.flags.writeable
            with pytest.raises(ValueError):
                actual.flat[0] = actual.flat[0]
        assert restored._runtime._kept_affine[1] == restored._seed
        assert restored_callback.calls == 1

    @pytest.mark.parametrize('transform', [nan_in_matrix, wrong_inverse_shape, materially_wrong_inverse,
                                           numerically_singular], ids=lambda f: f.__name__)
    def test_saved_invalid_transform_still_fails_closed_after_loading(self, transform):
        # A specification is data: it can be saved whatever its callbacks
        # return. The safeguard runs where the problem is built, so loading
        # cannot smuggle an invalid transform past it.
        for stages in ([custom_stage(transform)], [custom_stage(transform), 'noisy']):
            feature = pickle.loads(pickle.dumps(Feature(stages)))
            with pytest.raises(ValueError):
                FeaturedProblem(linear_problem(), feature, 10, 3)

    def test_saved_valid_transform_poses_the_same_problem_after_loading(self):
        problem = linear_problem()
        for transform in (exact_diagonal, roundoff_inverse, sloppy_inverse, dense):
            feature = Feature([custom_stage(transform)])
            before = FeaturedProblem(problem, feature, 10, 3)
            after = FeaturedProblem(problem, pickle.loads(pickle.dumps(feature)), 10, 3)
            for key in ('x0', 'xl', 'xu', 'aub', 'bub', 'aeq', 'beq'):
                np.testing.assert_array_equal(getattr(after, key), getattr(before, key), err_msg=key)


# ---------------------------------------------------------------- follow-up: the decision is exact

def wide_box():
    """Two variables between 0 and 1e16: next to bounds that large, an entry of 4e-16 is worth 4."""
    return Problem(sphere, [1.0, 1.0], xl=[0.0, 0.0], xu=[1e16, 1e16], reference=REFERENCE)


def build(problem, transform, composed=False, **options):
    stages = [{'name': 'custom', 'options': dict(options, mod_affine=transform)}] + (['noisy'] if composed else [])
    return FeaturedProblem(problem, Feature(stages), 20, 3)


COMPOSED = pytest.mark.parametrize('composed', [False, True], ids=['single', 'composed'])


def runtime_and_callback_of(featured):
    """The runtime of the custom stage built by `build`, and the ``mod_affine`` of the user in it."""
    runtime = featured._views[0]._runtime if isinstance(featured, ComposedFeaturedProblem) else featured._runtime
    callback = runtime._options['mod_affine']
    return runtime, getattr(callback, 'user', callback)  # a composition wraps the callbacks of a custom stage


class EarlierLayout(pickle.Pickler):
    """Writes a featured problem the way an earlier version of the package did."""

    def __init__(self, stream, layout):
        super().__init__(stream, protocol=4)
        self.layout = layout

    def reducer_override(self, obj):
        if not isinstance(obj, FeaturedProblem):
            return NotImplemented
        if self.layout == 'state among the arguments':  # the first reducer of the class
            return _restore_featured_problem, (type(obj), obj.__dict__)
        # A composition before the class had a reducer: ``__getnewargs__`` and the instance dictionary.
        return copyreg.__newobj__, (type(obj), obj._problem, obj._runtime, obj._max_eval, obj._seed), obj.__dict__


class TestPicklingRestoresWhatWasBuilt:

    @COMPOSED
    @pytest.mark.parametrize('how', ['protocol 4', 'protocol 5', 'deepcopy'])
    def test_a_callback_that_keeps_its_featured_problem_survives_pickling_and_copying(self, how, composed):
        # A reference back to the object from within its own state. pickle and
        # copy register an object after its reconstructor returned, so a state
        # handed over among the arguments of the reconstructor is restored
        # before the object exists: ``pickle.loads`` returned an object without
        # a single attribute here, silently, and ``copy.deepcopy`` a copy whose
        # callback kept another, half-built object.
        callback = Drifting()
        featured = build(linear_problem(), callback, composed)
        featured.fun(np.array([0.2, 0.3, 0.4]))
        callback.owner = featured
        restored = copy.deepcopy(featured) if how == 'deepcopy' else pickle.loads(pickle.dumps(featured, protocol=int(how[-1])))
        assert type(restored) is type(featured) and callback.calls == 1
        for key in ('x0', 'xl', 'xu', 'aub', 'bub', 'aeq', 'beq', 'fun_hist', 'maxcv_hist'):
            np.testing.assert_array_equal(getattr(restored, key), getattr(featured, key), err_msg=key)
        runtime, kept_callback = runtime_and_callback_of(restored)
        assert kept_callback.owner is restored and kept_callback.calls == 1
        assert all(not array.flags.writeable for array in runtime._kept_affine[2])

    @pytest.mark.parametrize('layout, composed', [('state among the arguments', False), ('state among the arguments', True),
                                                  ('constructor arguments and a dictionary', True)],
                             ids=['first reducer, single', 'first reducer, composed', 'before any reducer, composed'])
    def test_pickles_of_the_earlier_layouts_are_still_read(self, layout, composed):
        # Nothing of what is restored is produced again: the callback is not
        # called, the kept transformation is read-only, and the restored trial
        # continues as the live one does.
        callback = Drifting()
        featured = build(linear_problem(), callback, composed)
        featured.fun(np.array([0.2, 0.3, 0.4]))
        stream = io.BytesIO()
        EarlierLayout(stream, layout).dump(featured)
        restored = pickle.loads(stream.getvalue())
        assert type(restored) is type(featured) and callback.calls == 1
        for key in ('x0', 'xl', 'xu', 'aub', 'bub', 'aeq', 'beq', 'fun_hist', 'maxcv_hist'):
            np.testing.assert_array_equal(getattr(restored, key), getattr(featured, key), err_msg=key)
        runtime, kept_callback = runtime_and_callback_of(restored)
        assert kept_callback.calls == 1
        assert all(not array.flags.writeable for array in runtime._kept_affine[2])
        point = np.array([0.1, -0.2, 0.3])
        assert restored.fun(point) == featured.fun(point) and restored.maxcv(point) == featured.maxcv(point)
        assert kept_callback.calls == 1 and callback.calls == 1


class TestExactStructuralDecision:

    @COMPOSED
    def test_a_coupling_below_every_tolerance_is_posed_not_ignored(self, composed):
        # The base called 4e-16 negligible next to a diagonal of 1 (n * eps is
        # 4.4e-16) and posed the box [0, 1e16]**2 in the new variables.
        problem = wide_box()
        featured = build(problem, rounding_level_coupling, composed)
        A, b, _ = rounding_level_coupling(None, problem)
        assert_structure(featured, expected_generic_structure(problem, A, b))
        assert_no_bound_is_lost_or_doubled(featured, problem)
        # y = (-4, 1e16) is mapped to x = (0, 1e16), a vertex of the original
        # box, and that box in the new variables rejected it by 4; it accepted
        # (1e16, 1e16), which is mapped 4 outside.
        for y, violation in (([-4.0, 1e16], 0.0), ([1e16, 1e16], 4.0), ([0.0, 0.0], 0.0), ([3.0, 5e15], 0.0)):
            y = np.array(y)
            assert featured.maxcv(y) == pytest.approx(violation, abs=1e-9), y
            assert posed_violation(featured, y) == pytest.approx(violation, abs=1e-9), y
        assert featured.reference == problem.reference

    def test_the_decision_reads_the_matrix_exactly_and_the_inverse_for_its_scale_only(self):
        A, _, inv = rounding_level_coupling(None, None)
        assert _affine_is_diagonal(A, inv) is False
        assert _affine_is_diagonal(A.T, inv.T) is False
        for transform, diagonal in ((exact_diagonal, True), (roundoff_inverse, True), (sloppy_inverse, True),
                                    (coupling_hidden_by_scale, True), (extreme_diagonal, True),
                                    (overflowing_scale, True), (roundoff_matrix, False), (sloppy_matrix, False),
                                    (dense, False), (badly_scaled_rotation, False),
                                    (inaccurate_diagonal_inverse, False)):
            A, _, inv = transform(None, None)
            assert _affine_is_diagonal(A, inv) is diagonal, transform.__name__
        assert _affine_is_diagonal(np.zeros((0, 0)), np.zeros((0, 0))) is True
        # The reciprocal to roundoff, whichever way it was rounded, and no further.
        d = np.array([3.0, -7.0, 0.1, 1e-150, 2.0 ** 0.37, 1e150])
        assert _affine_is_diagonal(np.diag(d), np.diag(1.0 / d)) is True
        assert _affine_is_diagonal(np.diag(d), np.diag(np.nextafter(1.0 / d, np.inf))) is True
        assert _affine_is_diagonal(np.diag(d), np.diag(np.nextafter(1.0 / d, -np.inf))) is True
        assert _affine_is_diagonal(np.diag(d), np.diag((1.0 + 1e-13) / d)) is False

    @COMPOSED
    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    def test_a_diagonal_inverse_that_is_not_the_reciprocal_is_not_used_for_the_bounds(self, name, composed):
        # The pair is accepted: it is consistent to 1e-9 and the contract asks
        # for 1e-8. The shortcut would scale the bounds by that inverse; the
        # rows of A are exact whatever the inverse is.
        problem = PROBLEMS[name]()
        featured = build(problem, inaccurate_diagonal_inverse, composed)
        A, b, _ = inaccurate_diagonal_inverse(None, problem)
        assert_structure(featured, expected_generic_structure(problem, A, b))
        assert_no_bound_is_lost_or_doubled(featured, problem)
        assert_posed_problem_is_the_scored_problem(featured, problem)

    def test_bounds_of_any_size_are_posed_to_roundoff(self):
        problem = Problem(sphere, [1.0, 1.0, 1.0], xl=[0.0] * 3, xu=[1e16] * 3)
        featured = build(problem, inaccurate_diagonal_inverse)
        _, b, inv = inaccurate_diagonal_inverse(None, problem)
        # The upper bounds that diag(inv) gives: on the base the corner of the
        # posed box, and 1e7 outside the original one.
        scale = np.diagonal(inv)
        y = np.maximum(scale * (problem.xl - b), scale * (problem.xu - b))
        assert featured.maxcv(y) > 1e6
        assert posed_violation(featured, y) == pytest.approx(featured.maxcv(y), rel=1e-5)


# ---------------------------------------------------------------- follow-up: one transform per problem

class TestOneTransformPerProblem:

    @pytest.mark.parametrize('calls', [0, 1], ids=['diagonal-first', 'dense-first'])
    @COMPOSED
    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    def test_a_callback_with_a_state_cannot_mix_two_transforms(self, name, composed, calls):
        # On the base every modifier asked the callback again. With the dense
        # answer for the bounds and the diagonal one for the rows, the bounds
        # became infinite and no bound row was added: every finite bound was
        # lost. The other way round every bound was posed twice.
        problem = PROBLEMS[name]()
        callback = Alternating(calls)
        featured = build(problem, callback, composed)
        assert callback.calls == calls + 1  # asked once
        if calls == 0:
            assert_structure(featured, expected_diagonal_structure(problem, np.diag(D), B, np.diag(1.0 / D)))
        else:
            assert_structure(featured, expected_generic_structure(problem, DENSE, B))
        assert_no_bound_is_lost_or_doubled(featured, problem)
        assert_posed_problem_is_the_scored_problem(featured, problem)
        assert callback.calls == calls + 1  # and not again, 400 truth reads later
        assert featured.reference == problem.reference

    @COMPOSED
    def test_evaluations_use_the_transform_the_structure_was_built_with(self, composed):
        problem = linear_problem()
        callback = Drifting()
        featured = build(problem, callback, composed)
        A, b, _ = Drifting.rotation(1)  # the first answer, and the only one that may be used
        assert_structure(featured, expected_generic_structure(problem, A, b))
        np.testing.assert_allclose(A @ featured.x0 + b, problem.x0, atol=1e-14)
        assert featured.fun_init == pytest.approx(sphere(problem.x0), rel=1e-13)
        for y in np.random.default_rng(1).uniform(-2.0, 2.0, size=(6, 3)):
            featured.fun(y)
            assert featured.fun_hist[-1] == pytest.approx(sphere(A @ y + b), rel=1e-13)
            assert featured.maxcv_hist[-1] == pytest.approx(problem.maxcv(A @ y + b), rel=1e-13, abs=1e-13)
            assert featured.maxcv(y) == pytest.approx(problem.maxcv(A @ y + b), rel=1e-13, abs=1e-13)
        assert callback.calls == 1

    def test_the_transform_is_kept_per_seed_and_per_problem(self):
        first, second = linear_problem(), bounded_problem()
        callback = Counting()
        runtime = _StageRuntime('custom', {'mod_affine': callback})
        for _ in range(3):
            runtime.modifier_x0(3, first)
            runtime.modifier_bounds(3, first)
            runtime.modifier_linear_ub(3, first)
            runtime.modifier_linear_eq(3, first)
            runtime.modifier_affine(3, first)
        assert callback.calls == 1
        runtime.modifier_affine(4, first)
        assert callback.calls == 2  # another seed is another stream for the callback
        runtime.modifier_bounds(4, second)
        assert callback.calls == 3  # and another problem another question
        runtime.modifier_linear_ub(4, second)
        assert callback.calls == 3
        # A specification keeps no such state: every featured problem asks afresh.
        feature = Feature('custom', mod_affine=callback)
        FeaturedProblem(first, feature, 10, 3)
        FeaturedProblem(first, feature, 10, 3)
        assert callback.calls == 5

    def test_the_kept_transform_cannot_be_changed_from_outside(self):
        problem = linear_problem()
        for runtime in (_StageRuntime('custom', {'mod_affine': dense}),
                        _StageRuntime('linearly_transformed', {'rotated': True, 'condition_factor': 4}),
                        _StageRuntime('permuted', {})):
            for array in runtime.modifier_affine(3, problem):
                with pytest.raises(ValueError, match='read-only'):
                    array[...] = 0.0
        # What the callback returned is still the user's to change, without changing the kept transform.
        own = (DENSE.copy(), B.copy(), np.linalg.inv(DENSE))
        runtime = _StageRuntime('custom', {'mod_affine': lambda rng, problem: own})
        runtime.modifier_affine(3, problem)
        own[0][0, 0] = 5.0
        assert runtime.modifier_affine(3, problem)[0][0, 0] == 2.0

    def test_a_pickled_composition_keeps_the_transform_it_was_built_with(self):
        problem = linear_problem()
        featured = build(problem, Drifting(), composed=True)
        restored = pickle.loads(pickle.dumps(featured))
        A, b, _ = Drifting.rotation(1)
        assert_structure(restored, expected_generic_structure(problem, A, b))
        for y in np.random.default_rng(2).uniform(-2.0, 2.0, size=(4, 3)):
            featured.fun(y)
            restored.fun(y)
            assert restored.fun_hist[-1] == pytest.approx(sphere(A @ y + b), rel=1e-13)
        np.testing.assert_array_equal(restored.fun_hist, featured.fun_hist)
        np.testing.assert_array_equal(restored.maxcv_hist, featured.maxcv_hist)


# ---------------------------------------------------------------- follow-up: the inverse, from both sides

class TestTwoSidedInverse:

    @pytest.mark.parametrize('transform', [one_sided_inverse, one_sided_inverse_mirror], ids=lambda f: f.__name__)
    def test_the_inverse_has_to_invert_from_both_sides(self, transform):
        # Badly scaled, so that one product is the identity to 1e-16 while the
        # other misses it by 1: the base looked at the first only.
        A, b, inv = transform(None, None)
        assert np.linalg.norm(A @ inv - np.eye(3)) <= 1e-8 * 3 < 0.9 < np.linalg.norm(inv @ A - np.eye(3))
        y = np.array([0.0, 1.0, 0.0]) if transform is one_sided_inverse else np.array([1.0, 0.0, 0.0])
        assert np.max(np.abs(inv @ (A @ y) - y)) == pytest.approx(1.0)
        with pytest.raises(ValueError, match='not an identity matrix'):
            _checked_affine(A, b, inv, 3, supplied=True)
        # The transposed pairs fail on the other side, as they did before.
        with pytest.raises(ValueError, match='not an identity matrix'):
            _checked_affine(A.T, b, inv.T, 3, supplied=True)

    def test_neither_set_of_units_makes_a_consistent_pair_inconsistent(self):
        # A rotation with one variable in units of 1e13: new variable
        # (columns of A) or original one (rows). Roundoff times the ratio of
        # the units, 1e-3, is in one product or in the other; measured against
        # the terms the entries are summed from, both are roundoff.
        for transform, large in ((column_scaled_rotation, 1), (row_scaled_rotation, 0)):
            A, b, inv = transform(None, None)
            residuals = np.linalg.norm(A @ inv - np.eye(3)), np.linalg.norm(inv @ A - np.eye(3))
            assert residuals[large] > 1e-8 * 3 > 1e-12 > residuals[1 - large], transform.__name__  # the plain rule refuses one side
            checked = _checked_affine(A, b, inv, 3, supplied=True)
            for kept, given in zip(checked, (A, b, inv)):
                np.testing.assert_array_equal(kept, given)

    def test_what_was_inconsistent_stays_inconsistent(self):
        # The scale of the terms never excuses an error of the size of the terms.
        for transform in (singular, identity_as_inverse, materially_wrong_inverse):
            A, b, inv = transform(None, None)
            for pair in ((A, inv), (A.T, inv.T)):
                with pytest.raises(ValueError, match='not an identity matrix'):
                    _checked_affine(pair[0], b, pair[1], 3, supplied=True)
        A, b, wrong = row_scaled_rotation(None, None)
        wrong[0, 1] *= 1.0 + 1e-6  # one entry of the inverse off by 1e-6 of its size
        with pytest.raises(ValueError, match='not an identity matrix'):
            _checked_affine(A, b, wrong, 3, supplied=True)


# ---------------------------------------------------------------- follow-up: nothing finite leaves the range

def far_box():
    return Problem(sphere, [0.0, 0.0, 0.0], xl=[0.0] * 3, xu=[1e308] * 3)


class TestTransportOverflow:

    @COMPOSED
    @pytest.mark.parametrize('transform', [shifted_out_of_range, dense_shifted_out_of_range], ids=lambda f: f.__name__)
    def test_a_finite_bound_never_becomes_infinite_in_the_shift(self, transform, composed):
        # xu - b is 1e308 + 1e308. The scale guard of the base looked at the
        # product only, and an infinite factor is "no bound" to it: the upper
        # bounds were posed as infinite, or as rows with an infinite right-hand
        # side, which count as no constraint.
        with pytest.raises(ValueError, match='finite bound'):
            build(far_box(), transform, composed)

    @COMPOSED
    def test_a_fixed_variable_never_loses_its_equality_in_the_shift(self, composed):
        problem = Problem(sphere, [0.0, 0.0, 0.0], xl=[1e308, 0.0, 0.0], xu=[1e308, 1.0, 1.0])
        with pytest.raises(ValueError, match='finite bound'):
            build(problem, dense_shifted_out_of_range, composed)

    @COMPOSED
    def test_a_fixed_variable_is_guarded_where_its_equality_is_posed(self, composed):
        # With supplied inequality rows the inequality modifier of the
        # framework does not run (and nothing of its is replaced here: the only
        # finite bounds are those of the fixed variable). The equality of that
        # variable is still posed by the framework, and has its own guard.
        def rows(rng, problem):
            return np.array([[1.0, 1.0, 0.0]]), np.array([5.0])

        problem = Problem(sphere, [0.0, 0.0, 0.0], xl=[1e308, -np.inf, -np.inf], xu=[1e308, np.inf, np.inf])
        with pytest.raises(ValueError, match='finite bound'):
            build(problem, dense_shifted_out_of_range, composed, mod_linear_ub=rows)

    @COMPOSED
    @pytest.mark.parametrize('kind', ['ub', 'eq'])
    def test_a_right_hand_side_never_leaves_the_range(self, kind, composed):
        # bub - aub @ b with aub @ b = 2e318: on the base a right-hand side of
        # -inf (no point satisfies the row) or NaN.
        rows = {'aub': [[1e10, 1e10, 0.0]], 'bub': [1.0]} if kind == 'ub' else {'aeq': [[1e10, -3e10, 0.0]], 'beq': [2.0]}
        with pytest.raises(ValueError, match='linear constraint'):
            build(Problem(sphere, [0.0, 0.0, 0.0], **rows), shifted_the_other_way, composed)

    @COMPOSED
    @pytest.mark.parametrize('kind', ['ub', 'eq'])
    def test_a_composed_row_never_leaves_the_range(self, kind, composed):
        # aub @ A with 1e200 * 1e200: on the base a coefficient of inf.
        rows = {'aub': [[1e200, 1.0, 0.0]], 'bub': [1.0]} if kind == 'ub' else {'aeq': [[1e200, 0.0, 1.0]], 'beq': [0.0]}
        with pytest.raises(ValueError, match='linear constraint'):
            build(Problem(sphere, [0.0, 0.0, 0.0], **rows), huge_scaling, composed)

    @COMPOSED
    def test_the_initial_point_never_leaves_the_range(self, composed):
        with pytest.raises(ValueError, match='initial point'):
            build(Problem(sphere, [1e308, 1e308, 1e308]), shifted_out_of_range, composed)

    def test_the_scaling_feature_is_guarded_as_well(self):
        # diag(A) reaches 2**47 for this condition factor.
        feature = Feature('linearly_transformed', rotated=False, condition_factor=6000.0)
        with pytest.raises(ValueError, match='linear constraint'):
            FeaturedProblem(Problem(sphere, [0.0, 0.0, 0.0], aub=[[1e300, 0.0, 1e300]], bub=[1.0]), feature, 10, 3)
        with pytest.raises(ValueError, match='initial point'):
            FeaturedProblem(Problem(sphere, [1e300, 0.0, 0.0]), feature, 10, 3)

    @COMPOSED
    @pytest.mark.parametrize('transform', [exact_diagonal, dense], ids=lambda f: f.__name__)
    def test_a_bound_that_is_infinite_already_stays_as_it_is(self, transform, composed):
        # An infinite bound is no bound: nothing finite is lost, so nothing is
        # refused, however large the finite data next to it. (A `Problem` does
        # not accept an infinite right-hand side, so an infinite one in a
        # featured problem can only be an overflow.)
        problem = Problem(sphere, [0.5, 0.5, 0.5], xl=[-np.inf, -1e300, -np.inf], xu=[np.inf, np.inf, 1e300],
                          aub=[[1.0, 1.0, 0.0]], bub=[1e300])
        featured = build(problem, transform, composed)
        A, b, inv = transform(None, problem)
        expected = expected_diagonal_structure(problem, A, b, inv) if transform is exact_diagonal \
            else expected_generic_structure(problem, A, b)
        assert_structure(featured, expected)
        assert np.all(np.isfinite(featured.bub)) and np.all(np.isfinite(featured.aub))
        assert np.isfinite(featured.xl).sum() + np.isfinite(featured.xu).sum() == (2 if transform is exact_diagonal else 0)


# ---------------------------------------------------------------- follow-up: numbers from a callback

class TestCallbackNumbers:

    @COMPOSED
    @pytest.mark.parametrize('where', ['matrix', 'shift', 'inverse'])
    def test_an_integer_beyond_the_float_range_is_refused_as_every_other_invalid_output(self, where, composed):
        # In a composition the base let the OverflowError of the conversion
        # through; the contract of a custom stage is a ValueError that names it.
        huge, ones = 10 ** 400, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        triple = {'matrix': ([[huge, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0], ones),
                  'shift': (ones, [huge, 0, 0], ones),
                  'inverse': (ones, [0, 0, 0], [[huge, 0, 0], [0, 1, 0], [0, 0, 1]])}[where]
        with pytest.raises(ValueError, match='mod_affine' if composed else 'must be a real array'):
            build(bounded_problem(), lambda rng, problem: triple, composed)

    @COMPOSED
    def test_integers_are_numbers(self, composed):
        problem = linear_problem()
        integers = build(problem, lambda rng, problem: (np.diag([2, -4, 1]), np.array([1, -1, 2]), np.diag([0.5, -0.25, 1.0])), composed)
        floats = build(problem, lambda rng, problem: (np.diag([2.0, -4.0, 1.0]), np.array([1.0, -1.0, 2.0]), np.diag([0.5, -0.25, 1.0])), composed)
        for key in ('x0', 'xl', 'xu', 'aub', 'bub', 'aeq', 'beq'):
            assert getattr(integers, key).dtype == float, key
            np.testing.assert_array_equal(getattr(integers, key), getattr(floats, key), err_msg=key)


# ---------------------------------------------------------------- second follow-up: nothing finite is lost to underflow

TINY = np.finfo(float).tiny  # the smallest normal number, 2**-1022: the documented boundary


def power_scaling(exponent):
    """``x = 2**exponent * y`` in the first variable, with its exact inverse (powers of two: nothing is rounded)."""
    def transform(rng, problem):
        return np.diag([2.0 ** exponent, 1.0, 1.0]), np.zeros(3), np.diag([2.0 ** -exponent, 1.0, 1.0])
    return transform


class TestTransportUnderflow:

    @COMPOSED
    def test_a_nonempty_interval_never_collapses_to_a_point(self, composed):
        # The finding: xl = 1e-200, xu = 2e-200, A = 1e200, inv = 1e-200. Both
        # scaled bounds are 1e-400 and 2e-400, which round to 0: the base posed
        # the single point y = 0, which is mapped to x = 0, outside the interval.
        def huge(rng, problem):
            return np.array([[1e200]]), np.zeros(1), np.array([[1e-200]])

        stages = [custom_stage(huge)] + (['noisy'] if composed else [])
        with pytest.raises(ValueError, match='finite bound.*underflow'):
            # (x0 = 0 is pulled back exactly, so that the bounds are what is refused.)
            FeaturedProblem(Problem(sphere, [0.0], xl=[1e-200], xu=[2e-200]), Feature(stages), 10, 3)
        # From inside the interval the initial point is lost in the same way,
        # and it is asked for first: refused as well, for that reason.
        with pytest.raises(ValueError, match='initial point.*underflow'):
            FeaturedProblem(Problem(sphere, [1.5e-200], xl=[1e-200], xu=[2e-200]), Feature(stages), 10, 3)

    @COMPOSED
    @pytest.mark.parametrize('exponent, representable', [(1022 - 500, True), (1023 - 500, False), (1050 - 500, False), (1100 - 500, False)],
                             ids=['realmin', 'half-realmin', 'subnormal', 'zero'])
    def test_the_boundary_is_the_smallest_normal_number(self, exponent, representable, composed):
        # The bound 2**-500 scaled by 2**-exponent: exactly realmin is posed,
        # anything below it (a subnormal number, whose spacing is absolute, or
        # zero) is refused. Powers of two, so that nothing else is rounded.
        problem = Problem(sphere, [0.0, 0.0, 0.0], xl=[2.0 ** -500, -1.0, -1.0], xu=[1.0, 1.0, 1.0])
        if representable:
            featured = build(problem, power_scaling(exponent), composed)
            assert featured.xl[0] == TINY == 2.0 ** -500 * 2.0 ** -exponent
        else:
            with pytest.raises(ValueError, match='finite bound.*underflow'):
                build(problem, power_scaling(exponent), composed)

    @COMPOSED
    def test_a_zero_bound_is_zero_whatever_the_scale(self, composed):
        # 0 times anything is 0 exactly: nothing is lost, so nothing is refused.
        problem = Problem(sphere, [0.0, 0.0, 0.0], xl=[0.0, -1.0, -1.0], xu=[np.inf, 1.0, 1.0])
        featured = build(problem, power_scaling(900), composed)
        np.testing.assert_array_equal(featured.xl, [0.0, -1.0, -1.0])

    @COMPOSED
    @pytest.mark.parametrize('kind', ['ub', 'eq'])
    def test_a_coefficient_row_never_vanishes(self, kind, composed):
        # aub @ A with 2**-600 * 2**-600: on the base the row [0, 0, 0], which
        # is no constraint for bub >= 0 and satisfied nowhere for bub < 0.
        rows = {'aub': [[2.0 ** -600, 0.0, 0.0]], 'bub': [1.0]} if kind == 'ub' else {'aeq': [[2.0 ** -600, 0.0, 0.0]], 'beq': [0.0]}
        with pytest.raises(ValueError, match='linear constraint.*underflow'):
            build(Problem(sphere, [0.0, 0.0, 0.0], **rows), power_scaling(-600), composed)

    @COMPOSED
    @pytest.mark.parametrize('kind', ['ub', 'eq'])
    def test_a_right_hand_side_never_vanishes(self, kind, composed):
        # bub - aub @ b with bub = 0 and aub @ b = 2**-1200: the right-hand
        # side is -2**-1200 and the base posed 0.
        def shifted(rng, problem):
            return np.eye(3), np.array([2.0 ** -600, 0.0, 0.0]), np.eye(3)

        rows = {'aub': [[2.0 ** -600, 0.0, 0.0]], 'bub': [0.0]} if kind == 'ub' else {'aeq': [[2.0 ** -600, 0.0, 0.0]], 'beq': [0.0]}
        with pytest.raises(ValueError, match='linear constraint.*underflow'):
            build(Problem(sphere, [0.0, 0.0, 0.0], **rows), shifted, composed)
        # Next to a right-hand side that is a normal number the same shift is
        # below its rounding, and the posed right-hand side is the exact one rounded.
        rows['bub' if kind == 'ub' else 'beq'] = [5.0]
        featured = build(Problem(sphere, [0.0, 0.0, 0.0], **rows), shifted, composed)
        assert (featured.bub if kind == 'ub' else featured.beq)[-1] == 5.0

    @COMPOSED
    def test_the_initial_point_never_vanishes(self, composed):
        # inv @ (x0 - b) with 2**-600 * 2**-600: the base started at 0, where
        # this objective is -inf, for a problem whose fun(x0) is finite.
        def log_objective(x):
            return float(np.log(x[0]) + x[1] ** 2 + x[2] ** 2)

        problem = Problem(log_objective, [2.0 ** -600, 1.0, 1.0])
        assert np.isfinite(problem.fun(problem.x0))
        with pytest.raises(ValueError, match='initial point.*underflow'):
            build(problem, power_scaling(600), composed)

    def test_the_scaling_feature_is_guarded_as_well(self):
        # diag(inv) reaches 2**-47 for this condition factor: 1e-300 becomes a subnormal number.
        feature = Feature('linearly_transformed', rotated=False, condition_factor=6000.0)
        with pytest.raises(ValueError, match='finite bound.*underflow'):
            FeaturedProblem(Problem(sphere, [0.0, 0.0, 0.0], xl=[-1.0, -1.0, 1e-300], xu=[1.0, 1.0, 1.0]), feature, 10, 3)
        with pytest.raises(ValueError, match='linear constraint.*underflow'):
            FeaturedProblem(Problem(sphere, [0.0, 0.0, 0.0], aub=[[1e-300, 0.0, 0.0]], bub=[1.0]), feature, 10, 3)

    @COMPOSED
    def test_cancellation_is_not_underflow(self, composed):
        # An entry that is zero because its terms cancel has lost nothing: the
        # terms are normal numbers. Only an entry whose terms are all lost is refused.
        def mixing(rng, problem):
            A = np.array([[2.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            return A, np.array([1.0, 1.0, 0.0]), np.linalg.inv(A)

        problem = Problem(sphere, [1.0, 1.0, 0.0], aub=[[1.0, -1.0, 0.0]], bub=[0.0], aeq=[[1.0, -2.0, 0.0]], beq=[-1.0])
        featured = build(problem, mixing, composed)
        np.testing.assert_array_equal(featured.aub, [[1.0, 0.0, 0.0]])  # 1 * 1 - 1 * 1 in the second entry
        np.testing.assert_array_equal(featured.bub, [0.0])  # 0 - (1 * 1 - 1 * 1)
        np.testing.assert_array_equal(featured.aeq, [[0.0, -1.0, 0.0]])
        np.testing.assert_array_equal(featured.beq, [0.0])  # -1 - (1 - 2)
        np.testing.assert_array_equal(featured.x0, [0.0, 0.0, 0.0])  # inv @ (x0 - b) with x0 = b


# ---------------------------------------------------------------- second follow-up: derivatives in the solver's variables

def smooth(x):
    return float(np.sum(x ** 4) + x[0] * x[1] + np.sin(x[2]))


def smooth_grad(x):
    return 4.0 * x ** 3 + np.array([x[1], x[0], np.cos(x[2])])


def smooth_hess(x):
    H = np.diag(12.0 * x ** 2)
    H[0, 1] = H[1, 0] = 1.0
    H[2, 2] -= np.sin(x[2])
    return H


def smooth_cub(x):
    return np.array([x[0] ** 2 + x[1] * x[2] - 1.0, np.exp(x[0]) - x[2]])


def smooth_jcub(x):
    return np.array([[2.0 * x[0], x[2], x[1]], [np.exp(x[0]), 0.0, -1.0]])


def smooth_hcub(x):
    first = np.zeros((3, 3))
    first[0, 0], first[1, 2], first[2, 1] = 2.0, 1.0, 1.0
    second = np.zeros((3, 3))
    second[0, 0] = np.exp(x[0])
    return [first, second]


def smooth_ceq(x):
    return np.array([x[0] * x[1] * x[2] - 0.5])


def smooth_jceq(x):
    return np.array([[x[1] * x[2], x[0] * x[2], x[0] * x[1]]])


def smooth_hceq(x):
    return [np.array([[0.0, x[2], x[1]], [x[2], 0.0, x[0]], [x[1], x[0], 0.0]])]


def smooth_problem(**without):
    callbacks = dict(grad=smooth_grad, hess=smooth_hess, cub=smooth_cub, jcub=smooth_jcub, hcub=smooth_hcub,
                     ceq=smooth_ceq, jceq=smooth_jceq, hceq=smooth_hceq)
    for name in without:
        del callbacks[name]
    return Problem(smooth, [0.3, -0.4, 0.5], xl=[-2.0, -2.0, -2.0], xu=[2.0, 2.0, 2.0], **callbacks)


def central_differences(function, y, h=1e-5):
    columns = [(np.atleast_1d(function(y + h * e)) - np.atleast_1d(function(y - h * e))) / (2.0 * h) for e in np.eye(y.size)]
    return np.array(columns).T  # one row per output, one column per variable


def dense_with_shift(rng, problem):
    return DENSE.copy(), B.copy(), np.linalg.inv(DENSE)


VARIABLE_CHANGES = {
    'custom-dense': lambda: Feature('custom', mod_affine=dense_with_shift),
    'custom-diagonal': lambda: Feature('custom', mod_affine=exact_diagonal),
    'linearly_transformed': lambda: Feature('linearly_transformed', rotated=True, condition_factor=4),
    'permuted': lambda: Feature('permuted'),
}
Y = np.array([0.2, -0.3, 0.4])


class TestDerivativesInSolverVariables:

    def test_the_gradient_of_the_finding(self):
        # f(x) = ||x||**2, A = diag(2, 1), y = (1, 1): the function the solver
        # sees is 4 * y1**2 + y2**2 with gradient (8, 2). The base returned
        # grad(y) = (2, 2), the gradient of another function at another point.
        def scaling(rng, problem):
            return np.diag([2.0, 1.0]), np.zeros(2), np.diag([0.5, 1.0])

        problem = Problem(sphere, [1.0, 1.0], grad=lambda x: 2.0 * np.asarray(x, dtype=float), hess=lambda x: 2.0 * np.eye(2))
        featured = FeaturedProblem(problem, Feature('custom', mod_affine=scaling), 100, 3)
        y = np.array([1.0, 1.0])
        np.testing.assert_array_equal(featured.grad(y), [8.0, 2.0])
        np.testing.assert_array_equal(featured.hess(y), [[8.0, 0.0], [0.0, 2.0]])
        np.testing.assert_allclose(featured.grad(y), central_differences(featured.fun, y)[0], rtol=1e-6)

    @pytest.mark.parametrize('name', sorted(VARIABLE_CHANGES))
    def test_every_derivative_follows_the_chain_rule(self, name):
        problem = smooth_problem()
        featured = FeaturedProblem(problem, VARIABLE_CHANGES[name](), 10000, 5)
        A, b, _ = featured._runtime.modifier_affine(featured._seed, problem)
        assert np.count_nonzero(A - np.eye(3)) > 0, 'the change of variables of this seed is the identity'
        x = A @ Y + b
        # Exactly the chain rule of x = A @ y + b ...
        np.testing.assert_allclose(featured.grad(Y), A.T @ smooth_grad(x), rtol=1e-13, atol=1e-13)
        np.testing.assert_allclose(featured.hess(Y), A.T @ smooth_hess(x) @ A, rtol=1e-13, atol=1e-13)
        np.testing.assert_allclose(featured.jcub(Y), smooth_jcub(x) @ A, rtol=1e-13, atol=1e-13)
        np.testing.assert_allclose(featured.jceq(Y), smooth_jceq(x) @ A, rtol=1e-13, atol=1e-13)
        for returned, original in ((featured.hcub(Y), smooth_hcub(x)), (featured.hceq(Y), smooth_hceq(x))):
            assert len(returned) == len(original)
            for kept, H in zip(returned, original):
                np.testing.assert_allclose(kept, A.T @ H @ A, rtol=1e-13, atol=1e-13)
        # ... which is what the functions the solver evaluates vary by.
        np.testing.assert_allclose(featured.grad(Y), central_differences(featured.fun, Y)[0], rtol=1e-6, atol=1e-7)
        np.testing.assert_allclose(featured.hess(Y), central_differences(featured.grad, Y), rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(featured.jcub(Y), central_differences(featured.cub, Y), rtol=1e-6, atol=1e-7)
        np.testing.assert_allclose(featured.jceq(Y), central_differences(featured.ceq, Y), rtol=1e-6, atol=1e-7)
        for i, H in enumerate(featured.hcub(Y)):
            np.testing.assert_allclose(H, central_differences(lambda y: featured.jcub(y)[i], Y), rtol=1e-6, atol=1e-6)
        for i, H in enumerate(featured.hceq(Y)):
            np.testing.assert_allclose(H, central_differences(lambda y: featured.jceq(y)[i], Y), rtol=1e-6, atol=1e-6)

    def test_derivatives_cost_nothing_and_read_the_kept_transformation(self):
        callback = Counting()
        featured = FeaturedProblem(smooth_problem(), Feature('custom', mod_affine=callback), 10, 3)
        for method in ('grad', 'hess', 'jcub', 'jceq', 'hcub', 'hceq'):
            getattr(featured, method)(Y)
        assert callback.calls == 1
        assert (featured.n_eval_fun, featured.n_eval_cub, featured.n_eval_ceq) == (0, 0, 0)
        assert featured.fun_hist.size == 0 and featured.maxcv_hist.size == 0

    @pytest.mark.parametrize('name', sorted(VARIABLE_CHANGES))
    def test_an_absent_derivative_stays_absent(self, name):
        problem = Problem(smooth, [0.3, -0.4, 0.5], cub=smooth_cub, ceq=smooth_ceq)
        featured = FeaturedProblem(problem, VARIABLE_CHANGES[name](), 10, 5)
        for method in ('grad', 'hess', 'jcub', 'jceq'):
            assert getattr(featured, method)(Y).size == 0 == getattr(problem, method)(Y).size, method
        assert featured.hcub(Y) == [] == problem.hcub(Y) and featured.hceq(Y) == []

    @pytest.mark.parametrize('name', sorted(VARIABLE_CHANGES))
    def test_the_point_is_checked_as_before(self, name):
        featured = FeaturedProblem(smooth_problem(), VARIABLE_CHANGES[name](), 10, 5)
        for method in ('grad', 'hess', 'jcub', 'jceq', 'hcub', 'hceq'):
            with pytest.raises(ValueError, match=f'`x` for method `{method}` in problem must have size 3'):
                getattr(featured, method)([0.1, 0.2])

    @pytest.mark.parametrize('feature', ['plain', 'noisy', 'truncated', 'perturbed_x0', 'random_nan', 'quantized',
                                         'unrelaxable_constraints', 'nonquantifiable_constraints'])
    def test_without_a_change_of_variables_the_derivatives_are_those_of_the_problem(self, feature):
        # The established passthrough: these features change values or the
        # initial point, never the variables, and the derivative methods
        # describe the callbacks of the problem, not the observed values.
        problem = smooth_problem()
        featured = FeaturedProblem(problem, Feature(feature), 10, 5)
        for method in ('grad', 'hess', 'jcub', 'jceq'):
            np.testing.assert_array_equal(getattr(featured, method)(Y), getattr(problem, method)(Y), err_msg=method)
        for method in ('hcub', 'hceq'):
            for kept, original in zip(getattr(featured, method)(Y), getattr(problem, method)(Y)):
                np.testing.assert_array_equal(kept, original, err_msg=method)

    def test_a_composition_still_provides_no_derivative(self):
        featured = FeaturedProblem(smooth_problem(), Feature([custom_stage(dense_with_shift), 'noisy']), 10, 3)
        for method in ('grad', 'hess', 'jcub', 'jceq', 'hcub', 'hceq'):
            with pytest.raises(NotImplementedError, match='composed feature'):
                getattr(featured, method)(Y)


# ---------------------------------------------------------------- second follow-up: decisions made explicit

class TestExplicitDecisions:

    def test_consistency_allows_the_rounding_of_the_products_and_nothing_more(self):
        # 330805f measured each residual against max(1, terms), which allowed
        # 1e-8 OF THE TERMS. For this matrix the terms are 2e15, and an inverse
        # with one entry off by 1e-9 of its size was accepted: the initial point
        # was pulled back to a point that is mapped 1e6 away from x0. 3a43a19
        # refused it. What the terms justify is the rounding of the products.
        A = np.array([[1.0, 1e15], [0.0, 1.0]])
        exact, sloppy = np.array([[1.0, -1e15], [0.0, 1.0]]), np.array([[1.0, -1e15 + 1e6], [0.0, 1.0]])
        problem = Problem(sphere, [1.0, 1.0])
        featured = FeaturedProblem(problem, Feature('custom', mod_affine=lambda rng, problem: (A, np.zeros(2), exact)), 10, 3)
        np.testing.assert_array_equal(A @ featured.x0, problem.x0)
        for pair in ((A, sloppy), (A.T, sloppy.T)):
            with pytest.raises(ValueError, match='not an identity matrix'):
                _checked_affine(pair[0], np.zeros(2), pair[1], 2, supplied=True)
        with pytest.raises(ValueError, match='not an identity matrix'):
            FeaturedProblem(problem, Feature('custom', mod_affine=lambda rng, problem: (A, np.zeros(2), sloppy)), 10, 3)

    def test_a_pair_that_is_consistent_to_roundoff_is_accepted_at_any_condition_number(self):
        # Built as linearly_transformed builds its own, D @ Q.T and Q @ D**-1,
        # with a condition number of 1e12: each product misses the identity by
        # 1e-5, which is a few units of rounding of its terms.
        rng = np.random.default_rng(7)
        Q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
        d = np.array([1e-6, 1.0, 1e6])
        A, inv = np.diag(d) @ Q.T, Q @ np.diag(1.0 / d)
        assert max(np.linalg.norm(A @ inv - np.eye(3)), np.linalg.norm(inv @ A - np.eye(3))) > 1e-8 * 3
        kept = _checked_affine(A, np.zeros(3), inv, 3, supplied=True)
        np.testing.assert_array_equal(kept[0], A)

    @COMPOSED
    def test_mod_bounds_replaces_the_bounds_and_nothing_else(self, composed):
        # Replacement owns the logical original bounds, even when an affine
        # map would otherwise turn them into linear rows. The original linear
        # constraints must still be transported, including fixed-bound cases.
        def box(rng, problem):
            return np.full(3, -9.0), np.full(3, 9.0)

        for transform in (exact_diagonal, dense, roundoff_matrix):
            problem = linear_problem()
            problem._xu[0] = problem._xl[0]  # a fixed bound would become an equality row
            featured = build(problem, transform, composed, mod_bounds=box)
            A, b, _ = transform(None, problem)
            np.testing.assert_array_equal(featured.xl, np.full(3, -9.0))
            np.testing.assert_array_equal(featured.xu, np.full(3, 9.0))
            np.testing.assert_array_equal(featured.aub, problem.aub @ A)
            np.testing.assert_array_equal(featured.bub, problem.bub - problem.aub @ b)
            np.testing.assert_array_equal(featured.aeq, problem.aeq @ A)
            np.testing.assert_array_equal(featured.beq, problem.beq - problem.aeq @ b)
            assert featured.reference is None

    @COMPOSED
    def test_an_integer_beyond_2_53_is_rounded_as_every_number_is(self, composed):
        # Conversion to double precision is rounding to nearest, for an integer
        # as for a decimal literal, the same in both paths and both languages.
        # The rounded triple is the one that is validated, kept and used, by
        # the structure and by the truth alike.
        big = 2 ** 53 + 1
        integers = (np.diag(np.array([big, 1, 1], dtype=np.int64)), np.zeros(3, dtype=np.int64), np.diag([1.0 / big, 1.0, 1.0]))
        problem = bounded_problem()
        featured = build(problem, lambda rng, problem: integers, composed)
        rounded = build(problem, lambda rng, problem: (np.diag([2.0 ** 53, 1.0, 1.0]), np.zeros(3), np.diag([1.0 / big, 1.0, 1.0])), composed)
        for key in ('x0', 'xl', 'xu', 'aub', 'bub', 'aeq', 'beq'):
            np.testing.assert_array_equal(getattr(featured, key), getattr(rounded, key), err_msg=key)
        y = np.array([5e-16, 0.5, 0.5])  # mapped to x1 = 4.5, which violates xu = 3 by 1.5
        assert featured.maxcv(y) == rounded.maxcv(y) == problem.maxcv(np.array([2.0 ** 53 * 5e-16, 0.5, 0.5])) > 1.0

    def test_the_deprecated_conveniences_keep_nothing(self):
        # Each call of a deprecated Feature.modifier_* builds a runtime of its
        # own, so a callback with a state is asked again by each of them: one
        # transformation per problem is a property of FeaturedProblem.
        callback = Counting()
        feature = Feature('custom', mod_affine=callback)
        with pytest.warns(DeprecationWarning):
            feature.modifier_bounds(3, linear_problem())
            feature.modifier_linear_ub(3, linear_problem())
        assert callback.calls == 2


# ---------------------------------------------------------------- second follow-up: the initial point

def stray_inverse(rng, problem):
    """``A = I``; the inverse is off by 1e-12 in one entry, far inside every matrix tolerance."""
    return np.eye(2), np.zeros(2), np.array([[1.0, 1e-12], [0.0, 1.0]])


def ill_conditioned(entry):
    def transform(rng, problem):
        return np.array([[1.0, 1e15], [0.0, 1.0]]), np.zeros(2), np.array([[1.0, entry], [0.0, 1.0]])
    return transform


class TestInitialPoint:

    @COMPOSED
    def test_the_point_is_found_from_the_matrix_not_from_the_supplied_inverse(self, composed):
        # The base built x0 = (100, 1e14) for the feasible (0, 1e14): maxcv_init 99.
        problem = Problem(sphere, [0.0, 1e14], xl=[-1.0, 0.0], xu=[1.0, 1e14])
        featured = build(problem, stray_inverse, composed)
        np.testing.assert_array_equal(featured.x0, [0.0, 1e14])
        assert featured.maxcv_init == 0.0
        assert featured.fun_init == sphere(problem.x0)

    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    @pytest.mark.parametrize('transform', [exact_diagonal, dense, sloppy_inverse, badly_scaled_rotation,
                                           row_scaled_rotation, column_scaled_rotation], ids=lambda f: f.__name__)
    def test_the_point_is_mapped_back_to_the_original_one(self, transform, name):
        problem = PROBLEMS[name]()
        featured = build(problem, transform)
        A, b, _ = transform(None, problem)
        allowed = 64 * 3 * EPS * (np.abs(A) @ np.abs(featured.x0) + np.abs(b) + np.abs(problem.x0))
        assert np.all(np.abs(A @ featured.x0 + b - problem.x0) <= allowed)

    def test_an_inverse_that_is_good_at_the_point_gives_the_point_it_always_gave(self):
        # Bitwise: the framework's own inverse (measured at most 3.1 units of
        # the allowance of 64 over n up to 200 and condition factors up to
        # 6000), and a supplied one that is an inverse to roundoff.
        problem = linear_problem()
        for options in (dict(rotated=False, condition_factor=4), dict(rotated=True, condition_factor=40),
                        dict(rotated=True, condition_factor=6000.0)):
            featured = FeaturedProblem(problem, Feature('linearly_transformed', **options), 10, 3)
            _, _, inv = featured._runtime.modifier_affine(3, problem)
            np.testing.assert_array_equal(featured.x0, inv @ problem.x0)
        for transform in (exact_diagonal, dense, badly_scaled_rotation, row_scaled_rotation, column_scaled_rotation):
            _, b, inv = transform(None, problem)
            np.testing.assert_array_equal(build(problem, transform).x0, inv @ (problem.x0 - b), err_msg=transform.__name__)

    def test_the_allowance_is_the_documented_one_and_is_taken_by_component(self):
        from optiprofiler.opclasses import _ROUNDING  # (imported here: the module has to load on the sources before it)
        assert _ROUNDING == 64
        # By component: next to a component of 1e14 an error of 1e-3 in the
        # other one is 1e-17 of the norm, and it is not excused.
        A, x0 = np.eye(2), np.array([1.0, 1e14])
        inv = np.array([[1.0, 1e-17], [0.0, 1.0]])  # moves the first component by 1e-3
        assert (inv @ x0)[0] != 1.0
        np.testing.assert_array_equal(_pulled_back(A, inv, x0), x0)  # solved from A instead
        # Within the allowance the supplied inverse is believed: one unit in the last place.
        nudged = np.array([[np.nextafter(1.0, 2.0), 0.0], [0.0, 1.0]])
        np.testing.assert_array_equal(_pulled_back(A, nudged, x0), nudged @ x0)

    def test_the_allowance_is_the_rounding_of_the_evaluation_at_that_point(self):
        # x = (y1 + 1e15 * y2, y2): at x0 = (1/3, 1/7) the first component is
        # evaluated as 1.4e14 - 1.4e14 + 1/3, which double precision resolves to
        # 0.03, whatever y is. The exact inverse is therefore accepted, its
        # point is mapped back within the allowance, and that allowance is NOT
        # the rounding of x0: it is what the truth, which goes through the same
        # map at every point, can resolve there. Measuring against x0 alone
        # would refuse the rotations of linearly_transformed at ordinary
        # condition factors. What is refused is a point that is lost: overflow
        # (TestTransportOverflow) and underflow (TestTransportUnderflow).
        A = np.array([[1.0, 1e15], [0.0, 1.0]])
        inv = np.array([[1.0, -1e15], [0.0, 1.0]])
        x0 = np.array([1.0 / 3.0, 1.0 / 7.0])
        point = _pulled_back(A, inv, x0)
        np.testing.assert_array_equal(point, inv @ x0)
        error = np.abs(A @ point - x0)
        assert np.all(error <= 64 * 2 * EPS * (np.abs(A) @ np.abs(point) + np.abs(x0)))
        assert error[0] <= 0.0625 and 64 * 2 * EPS * abs(x0[0]) < 1e-14  # resolved to 2 ulp of 1.4e14, not to the rounding of 1/3

    def test_a_sloppy_inverse_of_an_ill_conditioned_matrix_is_refused_again(self):
        # 3a43a19 refused it; a tolerance of 1e-8 relative to the terms (2e15
        # here) accepted it. The allowance is now the rounding of the products.
        problem = Problem(sphere, [1.0, 1.0])
        np.testing.assert_allclose(build(problem, ill_conditioned(-1e15)).x0, [1.0 - 1e15, 1.0])
        with pytest.raises(ValueError, match='not an identity matrix'):
            build(problem, ill_conditioned(-1e15 + 1e6))
