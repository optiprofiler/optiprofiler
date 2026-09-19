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
"""

import pickle

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem
from optiprofiler.composition import AffineView, ComposedFeaturedProblem
from optiprofiler.opclasses import _StageRuntime


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
    @pytest.mark.parametrize('transform', [roundoff_inverse, roundoff_matrix], ids=lambda f: f.__name__)
    def test_roundoff_in_either_matrix_keeps_bounds_as_bounds(self, transform, name):
        # The defect. With roundoff in the inverse the bounds became infinite
        # and no bound row was added: the bounds were gone. With roundoff in
        # the matrix the bounds stayed AND were added again as rows.
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
    @pytest.mark.parametrize('transform', [sloppy_inverse, sloppy_matrix, coupling_hidden_by_scale],
                             ids=lambda f: f.__name__)
    def test_pair_that_is_not_diagonal_to_roundoff_takes_the_generic_path(self, transform, name):
        # Not roundoff, not inconsistent: the two matrices do not tell the same
        # structural story, so the representation that needs only A is used.
        # On the base a coupling in the inverse lost every bound, and one in
        # the matrix posed every bound twice.
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
    @pytest.mark.parametrize('where', ['inverse', 'matrix'])
    def test_roundoff_in_either_matrix_keeps_bounds_as_bounds(self, monkeypatch, where, name):
        def perturb(A, inv):
            target = inv if where == 'inverse' else A
            target[0, 1] = 0.25 * EPS * min(abs(target[0, 0]), abs(target[1, 1]))
            return A, inv

        perturb_linearly_transformed(monkeypatch, perturb)
        problem = PROBLEMS[name]()
        featured = FeaturedProblem(problem, Feature('linearly_transformed', **self.UNROTATED), 10, 3)
        A, b, inv = featured._runtime.modifier_affine(3, problem)
        assert (inv if where == 'inverse' else A)[0, 1] != 0.0  # the perturbation is in place
        assert_structure(featured, expected_diagonal_structure(problem, A, b, inv))
        assert_no_bound_is_lost_or_doubled(featured, problem)
        assert_posed_problem_is_the_scored_problem(featured, problem)
        assert featured.reference == problem.reference

    @pytest.mark.parametrize('name', sorted(PROBLEMS))
    def test_inverse_that_is_not_diagonal_to_roundoff_takes_the_generic_path(self, monkeypatch, name):
        def perturb(A, inv):
            inv[0, 1] = 1e-12 * abs(inv[0, 0])
            return A, inv

        perturb_linearly_transformed(monkeypatch, perturb)
        problem = PROBLEMS[name]()
        featured = FeaturedProblem(problem, Feature('linearly_transformed', **self.UNROTATED), 10, 3)
        A, b, _ = featured._runtime.modifier_affine(3, problem)
        assert_structure(featured, expected_generic_structure(problem, A, b))
        assert_no_bound_is_lost_or_doubled(featured, problem)
        assert_posed_problem_is_the_scored_problem(featured, problem)

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
        (['perturbed_x0', custom_stage(roundoff_matrix), 'truncated'], True),
        ([{'name': 'linearly_transformed', 'options': {'rotated': False, 'condition_factor': 4}},
          custom_stage(roundoff_inverse)], True),
        (['permuted', custom_stage(exact_diagonal)], True),
        ([custom_stage(sloppy_inverse), 'noisy'], False),
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
