"""An approximate-inverse residual alone must not certify a change of variables."""

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem
from optiprofiler.opclasses import _checked_affine, _equilibrated_affine_condition


@pytest.mark.parametrize('composed', [False, True])
def test_an_annihilating_claimed_inverse_cannot_hide_a_singular_map(composed):
    # A maps every solver point to (s, s). Its alleged inverse makes both
    # products zero, but the terms in those products give an allowance > 1.
    # Accepting it would change this problem's optimum from 0 to 2 while
    # retaining its optimum reference as 0.
    A = np.ones((2, 2))
    inverse = 1e14 * np.array([[1.0, -1.0], [-1.0, 1.0]])
    problem = Problem(lambda x: (x[0] - 1.0)**2 + (x[1] + 1.0)**2, [0.0, 0.0],
                      reference={'merit': 0.0, 'kind': 'optimum', 'source': 'analytic',
                                 'mapping': 'feasible_objective/1'})
    stage = {'name': 'custom', 'options': {'mod_affine': lambda rng, p: (A, np.zeros(2), inverse)}}
    feature = Feature([stage, 'noisy'] if composed else [stage])
    np.testing.assert_array_equal(A @ inverse, np.zeros((2, 2)))
    np.testing.assert_array_equal(inverse @ A, np.zeros((2, 2)))
    with pytest.raises(ValueError, match='numerically singular after equilibration'):
        FeaturedProblem(problem, feature, 10, 0)


@pytest.mark.parametrize('n', [2, 10, 50, 200])
def test_exact_rank_deficiency_is_detected_in_dense_integer_matrices(n):
    rng = np.random.default_rng(26 + n)
    A = rng.integers(-9, 10, size=(n, n)).astype(float)
    # The exact dependence is representable; no SVD-generated near-singular
    # input is used as evidence of a rank-deficient matrix.
    A[-1, :] = np.sum(A[:-1, :], axis=0)
    assert not _equilibrated_affine_condition(A) * np.finfo(float).eps < 1.0


@pytest.mark.parametrize('n', [2, 10, 50, 200])
def test_valid_dense_condition_1e14_is_not_rejected_only_for_its_scale(n):
    rng = np.random.default_rng(26)
    Q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    R, _ = np.linalg.qr(rng.standard_normal((n, n)))
    diagonal = np.geomspace(1.0, 1e-14, n)
    A = (Q * diagonal) @ R.T
    inverse = (R / diagonal) @ Q.T
    kept, _, _ = _checked_affine(A, np.zeros(n), inverse, n, supplied=True)
    np.testing.assert_array_equal(kept, A)


@pytest.mark.parametrize('scale_rows', [False, True])
def test_independent_check_preserves_changes_of_row_and_column_units(scale_rows):
    Q, _ = np.linalg.qr(np.random.default_rng(7).standard_normal((3, 3)))
    diagonal = np.array([1.0, 1e13, 1.0])
    A = diagonal[:, None] * Q if scale_rows else Q * diagonal
    inverse = Q.T / diagonal if scale_rows else Q.T / diagonal[:, None]
    kept, _, _ = _checked_affine(A, np.zeros(3), inverse, 3, supplied=True)
    np.testing.assert_array_equal(kept, A)


def test_extreme_diagonal_shear_and_empty_transform_remain_supported():
    diagonal = np.array([1e-150, -4.0, 1e150])
    pairs = [(np.diag(diagonal), np.diag(1.0 / diagonal)),
             (np.array([[1.0, 1e15], [0.0, 1.0]]), np.array([[1.0, -1e15], [0.0, 1.0]])),
             (np.empty((0, 0)), np.empty((0, 0)))]
    for A, inverse in pairs:
        n = A.shape[0]
        kept, _, _ = _checked_affine(A, np.zeros(n), inverse, n, supplied=True)
        np.testing.assert_array_equal(kept, A)


def test_equilibration_must_not_silently_drop_a_nonzero_entry():
    # A rank estimate for a rounded-to-zero coefficient would concern a
    # different matrix. Such an extreme scaling is refused explicitly.
    A = np.array([[1e300, 1e-300], [1e300, 2e-300]])
    with pytest.raises(ValueError, match='losing data during equilibration'):
        _equilibrated_affine_condition(A)
