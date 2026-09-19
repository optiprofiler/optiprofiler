import numpy as np
from scipy.linalg import qr
import re
import sys
import warnings
from collections.abc import Mapping
from decimal import Context, Decimal, InvalidOperation, ROUND_HALF_UP
from numpy.linalg import lstsq
from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint, minimize
from scipy import __version__ as _SCIPY_VERSION

import warnings
from .utils import FeatureName, FeatureOption, get_logger, shorten_log_message
from .feature_definitions import (_SPEC_TYPE_MESSAGE, Declaration, StageRecord, fold_option_names, normalize_entries,
                                  normalize_shorthand, reject_experiment_options, reject_flat_stage_options)
from .feature_definitions import is_stochastic as _stage_is_stochastic
from .feature_definitions import retains_reference as _stage_retains_reference
from .experiment import STRATEGY_COMPOSED, select_execution_strategy
from .legacy_compat import LegacyObject, historical_effective_options

def _round_truncated(value, digits):
    """Round decimal ties using MATLAB's default away-from-zero direction."""
    # Use the decimal representation rather than scaling a binary float:
    # this preserves decimal ties (e.g. 1.005) and avoids scale overflow.
    decimal_value = Decimal(str(float(value)))
    if decimal_value.as_tuple().exponent >= -int(digits):
        return float(value)
    # A float has at most 17 significant decimal digits; allow a carry and
    # isolate rounding from the caller's global decimal context.
    return float(decimal_value.quantize(
        Decimal((0, (1,), -int(digits))), rounding=ROUND_HALF_UP,
        context=Context(prec=18, rounding=ROUND_HALF_UP, Emin=-999999, Emax=999999,
                        capitals=1, clamp=0, flags=[], traps=[InvalidOperation]),
    ))


def _scipy_version_less_than(major, minor):
    """
    Check whether the installed SciPy version is older than ``major.minor``.
    """
    match = re.match(r'(\d+)\.(\d+)', _SCIPY_VERSION)
    if match is None:
        return False
    scipy_major, scipy_minor = (int(match.group(1)), int(match.group(2)))
    return (scipy_major, scipy_minor) < (major, minor)


class _StageRuntime:
    """
    Runtime implementation of one stage kind (private).

    Holds the validated stage-local options of one stage and implements the
    numerical modifiers the recorder and the composition views call. It is
    created fresh for every trial from a stage record of a ``Feature``
    specification; it never owns experiment state and is never a public type.
    """

    def __init__(self, name, options):
        self._name = name
        self._options = dict(options)

    @property
    def name(self):
        """
        Name of the feature.

        Returns
        -------
        str
            Name of the feature.
        """
        return self._name

    @property
    def options(self):
        """
        Options of the feature.

        Returns
        -------
        dict
            Options of the feature.
        """
        return dict(self._options)

    @property
    def is_stochastic(self):
        return _stage_is_stochastic(self._name, self._options)

    def modifier_x0(self, seed, problem):
        """
        Modify the initial point.

        Parameters
        ----------
        seed : int
            Seed used to generate random numbers.
        problem : `Problem`
            Problem for which the initial point is modified.

        Returns
        -------
        `numpy.ndarray`
            Modified initial point.
        """

        if self._name == FeatureName.CUSTOM:
            # If the user specifies a custom modifier for the initial point, use it.
            if FeatureOption.MOD_X0 in self._options:
                rng_custom = self.get_default_rng(seed)
                return self._options[FeatureOption.MOD_X0](rng_custom, problem)
            # If the user does not specify a custom modifier for the initial point but specifies
            # a custom affine transformation, we need to apply the inverse of the affine
            # transformation to the initial point.
            if FeatureOption.MOD_AFFINE in self._options:
                _, b, inv = self.modifier_affine(seed, problem)
                return inv @ (problem.x0 - b)
            else:
                return problem.x0
        elif self._name == FeatureName.PERTURBED_X0:
            # Use max(1, norm(x0)) to avoid no perturbation when x0 is zero.
            rng_perturbed_x0 = self.get_default_rng(seed)
            level = self._options[FeatureOption.PERTURBATION_LEVEL]
            # NumPy vectors historically give coordinatewise amplitudes. Keep
            # that behavior, validate their dimension here (unknown in Feature),
            # and make finite list/tuple vectors work like their array form.
            if isinstance(level, (list, tuple, np.ndarray)):
                level = np.asarray(level)
                if level.ndim == 1 and level.size not in (1, problem.n):
                    raise ValueError('Option `perturbation_level` must have length 1 or problem.n.')
            perturbation_level = level * max(1, np.linalg.norm(problem.x0))
            if self._options[FeatureOption.DISTRIBUTION] == 'gaussian':
                return problem.x0 + perturbation_level * rng_perturbed_x0.standard_normal(problem.n)
            elif self._options[FeatureOption.DISTRIBUTION] == 'spherical':
                perturbation = rng_perturbed_x0.standard_normal(problem.n)
                return problem.x0 + perturbation_level * perturbation / np.linalg.norm(perturbation)
            else:
                return problem.x0 + perturbation_level * self._options[FeatureOption.DISTRIBUTION](rng_perturbed_x0, problem.n)
        elif self._name == FeatureName.PERMUTED:
            # Note that we need to apply the reverse permutation to the initial point so that
            # the new problem is mathematically equivalent to the original one.
            rng_permuted = self.get_default_rng(seed)
            permutation = rng_permuted.permutation(problem.n)
            reverse_permutation = np.argsort(permutation)
            return problem.x0[reverse_permutation]
        elif self._name == FeatureName.LINEARLY_TRANSFORMED:
            # Apply the inverse of the affine transformation to the initial point.
            _, __, inv = self.modifier_affine(seed, problem)
            return inv @ problem.x0
        else:
            return problem.x0

    def modifier_affine(self, seed, problem):
        """
        Generate an invertible matrix A and a vector b for the affine transformation applied to the variables.

        Parameters
        ----------
        seed : int
            Seed used to generate random numbers.
        problem : `Problem`
            Problem for which the affine transformation is generated.

        Returns
        -------
        `numpy.ndarray`, `numpy.ndarray`, `numpy.ndarray`
            Matrix A, vector b, and inverse of the matrix A.
        """

        # Default values
        A = np.eye(problem.n)
        b = np.zeros(problem.n)
        inv = np.eye(problem.n)

        if self._name == FeatureName.CUSTOM:
            if FeatureOption.MOD_AFFINE in self._options:
                rng_custom = self.get_default_rng(seed)
                A, b, inv = self._options[FeatureOption.MOD_AFFINE](rng_custom, problem)
            # Check whether A * inv is an identity matrix.
            if np.linalg.norm(A @ inv - np.eye(problem.n)) > 1e-8 * problem.n:
                raise ValueError('The multiplication of the affine transformation matrix and its inverse is not an identity matrix.')
        elif self._name == FeatureName.PERMUTED:
            # Generate a random permutation matrix.
            rng_permuted = self.get_default_rng(seed)
            permutation = rng_permuted.permutation(problem.n)
            A = np.eye(problem.n)
            A = A[permutation, :]
            inv = A.T
            b = np.zeros(problem.n)
        elif self._name == FeatureName.LINEARLY_TRANSFORMED:
            # Generate A in the form D * Q' with D being a diagonal matrix and Q being an
            # orthogonal matrix.

            # Generate a random rotation matrix Q if the 'rotated' option is set to true.
            if self._options[FeatureOption.ROTATED]:
                """
                We generate a random orthogonal matrix Q following the uniform distribution
                on O(n). The method refers to a note written by the late Professor Nicholas
                Higham, see:
                https://nhigham.com/2020/04/22/what-is-a-random-orthogonal-matrix/
                and an answer from MathStackExchange, see:
                https://math.stackexchange.com/a/4891933/1088047
                """
                rng_linearly_transformed = self.get_default_rng(seed)
                rand_matrix = rng_linearly_transformed.standard_normal((problem.n, problem.n))
                q, r = qr(rand_matrix)
                q[:, np.diag(r) < 0] *= -1
            else:
                q = np.eye(problem.n)
            # The extreme exponents differ by sqrt(condition_factor * n / 2),
            # so cond(A) is 2**sqrt(condition_factor * n / 2) for n >= 2.
            # Orthogonal rotation does not change the singular values; n=1 has cond(A)=1.
            log_condition_number = np.sqrt(self._options[FeatureOption.CONDITION_FACTOR] * problem.n / 2)
            power = np.linspace(-log_condition_number/2, log_condition_number/2, problem.n)
            A = np.diag(2**power) @ q.T
            inv = q @ np.diag(2**(-power))
        return A, b, inv

    def modifier_bounds(self, seed, problem):
        """
        Modify the bounds.

        Parameters
        ----------
        seed : int
            Seed used to generate random numbers.
        problem : `Problem`
            Problem for which the bounds are modified.

        Returns
        -------
        `numpy.ndarray`, `numpy.ndarray`
            Lower bounds and upper bounds of the modified problem.
        """

        if self._name == FeatureName.CUSTOM:
            # If the user specifies a custom modifier for the bounds, use it.
            if FeatureOption.MOD_BOUNDS in self._options:
                rng_custom = self.get_default_rng(seed)
                return self._options[FeatureOption.MOD_BOUNDS](rng_custom, problem)
            if FeatureOption.MOD_AFFINE not in self._options:
                return problem.xl, problem.xu
            # If the user does not specify a custom modifier for the bounds but specifies a
            # custom affine transformation, we need to specially handle the bounds.
            _, b, inv = self.modifier_affine(seed, problem)
            if np.count_nonzero(inv - np.diag(np.diagonal(inv))) != 0:  # Check if inv is a diagonal matrix.
                return np.full(problem.n, -np.inf), np.full(problem.n, np.inf)
            # If the inverse of the affine transformation is diagonal, we can apply it to get
            # the modified bounds.
            xl_tmp = np.diag(inv) * (problem.xl - b)
            xu_tmp = np.diag(inv) * (problem.xu - b)
            return np.minimum(xl_tmp, xu_tmp), np.maximum(xl_tmp, xu_tmp)
        elif self._name == FeatureName.PERMUTED:
            # Note that we need to apply the reverse permutation to the bounds so that the new
            # problem is mathematically equivalent to the original one.
            rng_permuted = self.get_default_rng(seed)
            permutation = rng_permuted.permutation(problem.n)
            reverse_permutation = np.argsort(permutation)
            return problem.xl[reverse_permutation], problem.xu[reverse_permutation]
        elif self._name == FeatureName.LINEARLY_TRANSFORMED:
            # Apply the inverse of the affine transformation to the bounds.
            _, __, inv = self.modifier_affine(seed, problem)
            if np.count_nonzero(inv - np.diag(np.diagonal(inv))) != 0:  # Check if inv is a diagonal matrix.
                return np.full(problem.n, -np.inf), np.full(problem.n, np.inf)
            xl_tmp = np.diag(inv) * (problem.xl - np.zeros(problem.n))
            xu_tmp = np.diag(inv) * (problem.xu - np.zeros(problem.n))
            return np.minimum(xl_tmp, xu_tmp), np.maximum(xl_tmp, xu_tmp)
        else:
            return problem.xl, problem.xu

    def modifier_linear_ub(self, seed, problem):
        """
        Modify the linear inequality constraints.

        Parameters
        ----------
        seed : int
            Seed used to generate random numbers.
        problem : `Problem`
            Problem for which the linear inequality constraints are modified.

        Returns
        -------
        `numpy.ndarray`, `numpy.ndarray`
            Coefficients of the linear inequality constraints and the right-hand side values of the modified problem.
        """

        if self._name == FeatureName.CUSTOM:
            # If the user specifies a custom modifier for the linear inequality constraints, use it.
            if FeatureOption.MOD_LINEAR_UB in self._options:
                rng_custom = self.get_default_rng(seed)
                return self._options[FeatureOption.MOD_LINEAR_UB](rng_custom, problem)
            if FeatureOption.MOD_AFFINE not in self._options:
                return problem.aub, problem.bub
            # If the user does not specify a custom modifier for the linear inequality
            # constraints but specifies a custom affine transformation, we need to specially
            # handle the linear inequality constraints.
            A, b = self.modifier_affine(seed, problem)[:2]
            if np.count_nonzero(A - np.diag(np.diagonal(A))) == 0:
                return problem.aub @ A, problem.bub - problem.aub @ b
            """
            We need to specially handle bound constraints and linear inequality constraints.

            Bound constraints
            xl <= A @ x + b <= xu
            should be modified to
            A @ x <= xu - b
            and
            -A @ x <= -xl + b
            when A is not diagonal.

            Linear inequality constraints
            turn out to be
            (aub @ A) @ x <= bub - aub @ b
            """

            # Pick out the indices of lower bounds who are not -Inf and upper bounds who are
            # not Inf since later we will not transform them into linear inequality constraints.
            idx_lb = ~np.isinf(problem.xl)
            idx_ub = ~np.isinf(problem.xu)
            # Remove the indices, of which the lower and upper bound are equal, since later
            # we will put them into the linear equality constraints.
            idx_eq = np.where(problem.xl == problem.xu)[0]
            idx_lb[idx_eq] = False
            idx_ub[idx_eq] = False
            if problem.aub.size == 0:
                return np.vstack([A[idx_ub, :], -A[idx_lb, :]]), \
                       np.concatenate([problem.xu[idx_ub] - b[idx_ub], -(problem.xl[idx_lb] - b[idx_lb])])
            else:
                return np.vstack([A[idx_ub, :], -A[idx_lb, :], problem.aub @ A]), \
                       np.concatenate([problem.xu[idx_ub] - b[idx_ub], -(problem.xl[idx_lb] - b[idx_lb]), problem.bub - problem.aub @ b])
        elif self._name == FeatureName.PERMUTED:
            rng_permuted = self.get_default_rng(seed)
            permutation = rng_permuted.permutation(problem.n)
            reverse_permutation = np.argsort(permutation)
            return problem.aub[:, reverse_permutation], problem.bub
        elif self._name == FeatureName.LINEARLY_TRANSFORMED:
            # Similar to the case in the custom feature where a custom affine transformation
            # is specified.
            A = self.modifier_affine(seed, problem)[0]
            if np.count_nonzero(A - np.diag(np.diagonal(A))) == 0:
                return problem.aub @ A, problem.bub
            idx_lb = ~np.isinf(problem.xl)
            idx_ub = ~np.isinf(problem.xu)
            idx_eq = np.where(problem.xl == problem.xu)[0]
            idx_lb[idx_eq] = False
            idx_ub[idx_eq] = False
            if problem.aub.size == 0:
                return np.vstack([A[idx_ub, :], -A[idx_lb, :]]), \
                       np.concatenate([problem.xu[idx_ub], -problem.xl[idx_lb]])
            else:
                return np.vstack([A[idx_ub, :], -A[idx_lb, :], problem.aub @ A]), \
                       np.concatenate([problem.xu[idx_ub], -problem.xl[idx_lb], problem.bub])
        else:
            return problem.aub, problem.bub

    def modifier_linear_eq(self, seed, problem):
        """
        Modify the linear equality constraints.

        Parameters
        ----------
        seed : int
            Seed used to generate random numbers.
        problem : `Problem`
            Problem for which the linear equality constraints are modified.

        Returns
        -------
        `numpy.ndarray`, `numpy.ndarray`
            Coefficients of the linear equality constraints and the right-hand side values of the modified problem.
        """

        if self._name == FeatureName.CUSTOM:
            # If the user specifies a custom modifier for the linear equality constraints, use it.
            if FeatureOption.MOD_LINEAR_EQ in self._options:
                rng_custom = self.get_default_rng(seed)
                return self._options[FeatureOption.MOD_LINEAR_EQ](rng_custom, problem)
            if FeatureOption.MOD_AFFINE not in self._options:
                return problem.aeq, problem.beq
            # If the user does not specify a custom modifier for the linear equality
            # constraints but specifies a custom affine transformation, we need to specially
            # handle the linear equality constraints.
            A, b = self.modifier_affine(seed, problem)[:2]
            if np.count_nonzero(A - np.diag(np.diagonal(A))) == 0:
                return problem.aeq @ A, problem.beq - problem.aeq @ b
            """
            We need to specially handle bound constraints and linear equality constraints.

            Bound constraints
            xl <= A @ x + b <= xu
            with xl = xu should be modified to
            A @ x = xu - b
            when A is not diagonal.

            Linear equality constraints
            turn out to be
            (aeq @ A) @ x = beq - aeq @ b
            """

            # Pick out the indices, of which the lower and upper bound are equal.
            idx_eq = np.where(problem.xl == problem.xu)[0]
            if problem.aeq.size == 0:
                return A[idx_eq, :], problem.xu[idx_eq] - b[idx_eq]
            else:
                return np.vstack([A[idx_eq, :], problem.aeq @ A]), \
                       np.concatenate([problem.xu[idx_eq] - b[idx_eq], problem.beq - problem.aeq @ b])
        elif self._name == FeatureName.PERMUTED:
            rng_permuted = self.get_default_rng(seed)
            permutation = rng_permuted.permutation(problem.n)
            reverse_permutation = np.argsort(permutation)
            return problem.aeq[:, reverse_permutation], problem.beq
        elif self._name == FeatureName.LINEARLY_TRANSFORMED:
            # Similar to the case in the custom feature where a custom affine transformation
            # is specified.
            A = self.modifier_affine(seed, problem)[0]
            if np.count_nonzero(A - np.diag(np.diagonal(A))) == 0:
                return problem.aeq @ A, problem.beq
            idx_eq = np.where(problem.xl == problem.xu)[0]
            if problem.aeq.size == 0:
                return A[idx_eq, :], problem.xu[idx_eq]
            else:
                return np.vstack([A[idx_eq, :], problem.aeq @ A]), \
                       np.concatenate([problem.xu[idx_eq], problem.beq])
        else:
            return problem.aeq, problem.beq

    def modifier_fun(self, x, seed, problem, n_eval):
        """
        Modify the objective function value.

        Parameters
        ----------
        x : `numpy.ndarray`, shape (n,)
            Point at which the objective function is evaluated.
        seed : int
            Seed used to generate random numbers.
        problem : `Problem`
            Problem for which the objective function is modified.
        n_eval: int
            Number of evaluations of the objective function.
            (We will use it to generate random streams so that evaluating
            the same point multiple times will not lead to the same
            random numbers.)

        Returns
        -------
        float
            Modified objective function value.
        """

        f = problem.fun(x)
        if self._name == FeatureName.CUSTOM:
            if FeatureOption.MOD_FUN in self._options:
                rng_custom = self.get_default_rng(seed, f, *x, n_eval)
                return self._options[FeatureOption.MOD_FUN](x, rng_custom, problem)
            else:
                return f
        elif self._name == FeatureName.NOISY:
            noise = self._compute_noise(x, seed, n_eval, f)
            return self._apply_noise(f, noise)
        elif self._name == FeatureName.RANDOM_NAN:
            return self._random_nan_scalar(f, x, seed, n_eval)
        elif self._name == FeatureName.TRUNCATED:
            return self._truncate_scalar(f, x, seed, n_eval)
        elif self._name == FeatureName.UNRELAXABLE_CONSTRAINTS:
            _, maxcv_bounds, maxcv_linear, maxcv_nonlinear = problem._maxcv(x)
            if self._options[FeatureOption.UNRELAXABLE_BOUNDS] and maxcv_bounds > 0.0:
                return np.inf
            elif self._options[FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS] and maxcv_linear > 0.0:
                return np.inf
            elif self._options[FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS] and maxcv_nonlinear > 0.0:
                return np.inf
            else:
                return f
        elif self._name == FeatureName.QUANTIZED:
            return problem.fun(self._quantize_point(x))
        else:
            return f

    def modifier_cub(self, x, seed, problem, n_eval_cub):
        """
        Modify the values of the nonlinear inequality constraints.

        Parameters
        ----------
        x : `numpy.ndarray`, shape (n,)
            Point at which the nonlinear inequality constraints are evaluated.
        seed : int
            Seed used to generate random numbers.
        problem : `Problem`
            Problem for which the nonlinear inequality constraints are modified.
        n_eval_cub : int
            Number of evaluations of the nonlinear inequality constraints.
            (We will use it to generate random streams so that evaluating
            the same point multiple times will not lead to the same
            random numbers.)

        Returns
        -------
        `numpy.ndarray`
            Modified values of the nonlinear inequality constraints.
        """

        cub = problem.cub(x)
        if cub.size == 0:
            return cub
        elif self._name == FeatureName.CUSTOM:
            if FeatureOption.MOD_CUB in self._options:
                rng_custom = self.get_default_rng(seed, *cub, *x, n_eval_cub)
                return self._options[FeatureOption.MOD_CUB](x, rng_custom, problem)
            else:
                return cub
        elif self._name == FeatureName.NOISY:
            # Similar to the case in the modifier_fun method.
            noise = self._compute_noise(x, seed, n_eval_cub, cub, cub.size)
            return self._apply_noise(cub, noise)
        elif self._name == FeatureName.RANDOM_NAN:
            # Similar to the case in the modifier_fun method.
            return self._random_nan_vector(cub, x, seed, n_eval_cub)
        elif self._name == FeatureName.TRUNCATED:
            # Similar to the case in the modifier_fun method.
            return self._truncate_vector(cub, x, seed, n_eval_cub)
        elif self._name == FeatureName.NONQUANTIFIABLE_CONSTRAINTS:
            return self._nonquantifiable_cub(cub)
        elif self._name == FeatureName.QUANTIZED:
            # Similar to the case in the modifier_fun method.
            return problem.cub(self._quantize_point(x))
        else:
            return cub

    def modifier_ceq(self, x, seed, problem, n_eval_ceq):
        """
        Modify the values of the nonlinear equality constraints.

        Parameters
        ----------
        x : `numpy.ndarray`, shape (n,)
            Point at which the nonlinear equality constraints are evaluated.
        seed : int
            Seed used to generate random numbers.
        problem : `Problem`
            Problem for which the nonlinear equality constraints are modified.
        n_eval_ceq : int
            Number of evaluations of the nonlinear equality constraints.
            (We will use it to generate random streams so that evaluating
            the same point multiple times will not lead to the same
            random numbers.)

        Returns
        -------
        `numpy.ndarray`
            Modified values of the nonlinear equality constraints.
        """

        ceq = problem.ceq(x)
        if ceq.size == 0:
            return ceq
        elif self._name == FeatureName.CUSTOM:
            if FeatureOption.MOD_CEQ in self._options:
                rng_custom = self.get_default_rng(seed, *ceq, *x, n_eval_ceq)
                return self._options[FeatureOption.MOD_CEQ](x, rng_custom, problem)
            else:
                return ceq
        elif self._name == FeatureName.NOISY:
            # Similar to the case in the modifier_fun method.
            noise = self._compute_noise(x, seed, n_eval_ceq, ceq, ceq.size)
            return self._apply_noise(ceq, noise)
        elif self._name == FeatureName.RANDOM_NAN:
            # Similar to the case in the modifier_fun method.
            return self._random_nan_vector(ceq, x, seed, n_eval_ceq)
        elif self._name == FeatureName.TRUNCATED:
            # Similar to the case in the modifier_fun method.
            return self._truncate_vector(ceq, x, seed, n_eval_ceq)
        elif self._name == FeatureName.NONQUANTIFIABLE_CONSTRAINTS:
            return self._nonquantifiable_ceq(ceq)
        elif self._name == FeatureName.QUANTIZED:
            # Similar to the case in the modifier_fun method.
            return problem.ceq(self._quantize_point(x))
        else:
            return ceq

    # The value helpers below hold the arithmetic of the value features. The
    # modifiers above apply them to a value they evaluated themselves; the
    # views of ``optiprofiler.composition`` apply them to the value served by
    # a predecessor. Vector helpers modify their argument in place.

    def _random_nan_scalar(self, f, x, seed, n_eval):
        rng_random_nan = self.get_default_rng(seed, f, *x, n_eval)
        if rng_random_nan.random() < self._options[FeatureOption.NAN_RATE]:
            return np.nan
        else:
            return f

    def _random_nan_vector(self, values, x, seed, n_eval):
        rng_random_nan = self.get_default_rng(seed, *values, *x, n_eval)
        values[rng_random_nan.random(values.size) < self._options[FeatureOption.NAN_RATE]] = np.nan
        return values

    def _truncate_scalar(self, f, x, seed, n_eval):
        if np.isnan(f) or np.isinf(f):
            # If f is NaN or Inf, we do not need to truncate it.
            # Note that if f is NaN or Inf, digits will be set to NaN or Inf respectively, which will lead
            # to an error when calling 'round(f, digits)'.
            return f
        rng_truncated = self.get_default_rng(seed, f, *x, n_eval)
        if f == 0.0:
            digits = self._options[FeatureOption.SIGNIFICANT_DIGITS] - 1
        else:
            # Floor matters below one: int would truncate a negative logarithm toward zero.
            digits = self._options[FeatureOption.SIGNIFICANT_DIGITS] - int(np.floor(np.log10(np.abs(f)))) - 1
        f = _round_truncated(f, digits)
        # Round f to the desired number of significant digits.
        if self._options[FeatureOption.PERTURBED_TRAILING_DIGITS]:
            if f >= 0.0:
                f += rng_truncated.uniform(0.0, 10.0 ** (-digits))
            else:
                f -= rng_truncated.uniform(0.0, 10.0 ** (-digits))
        return f

    def _truncate_vector(self, values, x, seed, n_eval):
        rng_truncated = self.get_default_rng(seed, *values, *x, n_eval)
        digits = np.zeros(values.size, dtype=int)
        finite = np.isfinite(values)
        nonzero = finite & (values != 0.0)
        digits[values == 0.0] = self._options[FeatureOption.SIGNIFICANT_DIGITS] - 1
        # Do not cast NaN/Inf exponents to integers or perturb them.
        digits[nonzero] = self._options[FeatureOption.SIGNIFICANT_DIGITS] - np.floor(np.log10(np.abs(values[nonzero]))).astype(int) - 1
        for i in range(values.size):
            if not np.isnan(values[i]) and not np.isinf(values[i]):
                values[i] = _round_truncated(values[i], digits[i])
        if self._options[FeatureOption.PERTURBED_TRAILING_DIGITS]:
            positive = finite & (values >= 0.0)
            negative = finite & (values < 0.0)
            values[positive] += rng_truncated.uniform(0.0, 10.0 ** (-digits[positive]))
            values[negative] -= rng_truncated.uniform(0.0, 10.0 ** (-digits[negative]))
        return values

    @staticmethod
    def _nonquantifiable_cub(values):
        # Set the elements whose value are less than or equal to 0 to 0.
        values[values <= 0.0] = 0.0
        # Set the rest to 1.
        values[values > 0.0] = 1.0
        return values

    @staticmethod
    def _nonquantifiable_ceq(values):
        # Set the elements whose absolute value are less than or equal to 10^(-6) to 0.
        values[np.abs(values) <= 1e-6] = 0.0
        # Set the rest to 1.
        values[np.abs(values) > 1e-6] = 1.0
        return values

    def _quantize_point(self, x):
        mesh_size = self._options[FeatureOption.MESH_SIZE]
        if self._options[FeatureOption.MESH_TYPE] == 'relative':
            mesh_size *= np.maximum(1, np.abs(x))
        return mesh_size * np.round(x / mesh_size)

    def _compute_noise(self, x, seed, n_eval, base_values, noise_size=None):
        if self._options[FeatureOption.NOISE_MODE] == 'deterministic':
            noise = self._evaluate_noise_map(x)
            if noise_size is None:
                return noise
            return np.full(noise_size, noise)

        base_values = np.asarray(base_values, dtype=float).ravel()
        rng_noisy = self.get_default_rng(seed, *base_values, *x, n_eval)
        if self._options[FeatureOption.DISTRIBUTION] == 'gaussian':
            if noise_size is None:
                return rng_noisy.standard_normal()
            return rng_noisy.standard_normal(noise_size)
        elif self._options[FeatureOption.DISTRIBUTION] == 'uniform':
            if noise_size is None:
                return rng_noisy.uniform(-1, 1)
            return rng_noisy.uniform(-1, 1, noise_size)
        elif noise_size is None:
            return self._sample_scalar_distribution(rng_noisy)
        else:
            return self._options[FeatureOption.DISTRIBUTION](rng_noisy, noise_size)

    def _apply_noise(self, value, noise):
        if self._options[FeatureOption.NOISE_TYPE] == 'absolute':
            return value + self._options[FeatureOption.NOISE_LEVEL] * noise
        elif self._options[FeatureOption.NOISE_TYPE] == 'relative':
            return value * (1.0 + self._options[FeatureOption.NOISE_LEVEL] * noise)
        else:
            return value + np.maximum(1, np.abs(value)) * self._options[FeatureOption.NOISE_LEVEL] * noise

    def _sample_scalar_distribution(self, rng):
        distribution = self._options[FeatureOption.DISTRIBUTION]
        try:
            noise = distribution(rng, 1)
        except TypeError:
            noise = distribution(rng)

        noise = np.asarray(noise)
        if noise.shape == ():
            return float(noise)
        if noise.size == 1:
            return float(noise.ravel()[0])
        raise ValueError('The output of `distribution` must be scalar for objective noise.')

    def _evaluate_noise_map(self, x):
        noise_map = self._options[FeatureOption.NOISE_MAP]
        if callable(noise_map):
            noise = noise_map(x)
        elif noise_map == 'chebyshev':
            noise = self.chebyshev_noise_map(x)
        else:
            raise ValueError('Option `noise_map` must be "chebyshev" or it must be callable.')

        noise = np.asarray(noise)
        if noise.shape != () or not np.isrealobj(noise):
            raise ValueError('The output of `noise_map` must be a real scalar.')
        return float(noise)

    @staticmethod
    def chebyshev_noise_map(x):
        # Deterministic noise map from Moré and Wild, "Benchmarking
        # derivative-free optimization algorithms" (2009).
        alpha = 0.9 * np.sin(100 * np.linalg.norm(x, 1)) * np.cos(100 * np.linalg.norm(x, np.inf)) + 0.1 * np.cos(np.linalg.norm(x, 2))
        return alpha * (4 * alpha ** 2 - 3)

    @staticmethod
    def get_default_rng(seed, *args):
        """
        Generate a random number generator.

        Parameters
        ----------
        seed : int
            Seed used to generate an initial random number generator.
        *args : tuple of int or float
            Arguments used to generate the returned random number generator.

        Returns
        -------
        `numpy.random.Generator`
            Random number generator.
        """
        # Preprocess the seed.
        if seed is not None:
            if np.isnan(seed) or np.isinf(seed):
                seed = 0
            if isinstance(seed, (float, np.floating)) and float(seed).is_integer():
                seed = int(seed)
            if isinstance(seed, np.integer):
                seed = int(seed)
            if not isinstance(seed, int):
                raise TypeError('The argument seed must be an integer.')
            if seed < 0:
                raise ValueError('The argument seed must be nonnegative.')

        # Convert the arguments to numbers.
        args = [float(arg) for arg in args if isinstance(arg, (int, float)) and not np.isnan(arg) and not np.isinf(arg)]

        # Create a random number generator based on the specified seed.
        rng = np.random.default_rng(seed if seed is not None else 0)

        # Generate a new seed based on the initial rand_stream and the additional arguments
        try:
            rand_val = rng.standard_normal()
            if not np.isfinite(rand_val):
                rand_val = 0.0
            sin_rand = np.sin(1e5 * rand_val)
            sin_args = sum(np.sin(1e5 * arg) for arg in args if np.isfinite(arg))

            new_seed = 1e9 * np.abs(sin_rand + sin_args)
            if not np.isfinite(new_seed):
                new_seed = 42
            new_seed = int(new_seed) % 2**32
        except:
            new_seed = 42

        return np.random.default_rng(new_seed)


FEATURE_NATIVE_VERSION = 1


def _rebuild_feature(version, route, declared_entries, effective_entries):
    """
    Rebuild a pickled ``Feature`` from its native form: the canonical effective
    stages (every validated local option explicit, so the defaults of the
    loading version cannot change the numbers) and, separately, the declaration
    as provenance. The effective entries are validated again; an identity
    pipeline has none.
    """
    if version != FEATURE_NATIVE_VERSION:
        raise ValueError(f'Unsupported native Feature form version {version!r} (this version reads {FEATURE_NATIVE_VERSION}).')
    if effective_entries:
        _, stages = normalize_entries([
            {'name': name, 'options': historical_effective_options(name, options)}
            for name, options in effective_entries])
    else:
        stages = ()
    feature = Feature.__new__(Feature)
    object.__setattr__(feature, '_declared', Declaration(route, declared_entries))
    object.__setattr__(feature, '_stages', stages)
    return feature


class Feature:
    """
    Specification of the feature applied to the benchmarked problems: an
    ordered pipeline of stages, each with its own validated stage-local options.

    A single feature is a one-stage pipeline; ``'plain'`` is the identity with
    no effective stage. The specification is immutable and carries no experiment
    state: the number of runs belongs to ``benchmark(..., n_runs=N)``, and every
    trial builds fresh runtime state from the records.

    Parameters
    ----------
    name : str, dict, list, tuple or Feature
        Either a feature name, several names joined with ``+`` to apply them in
        order (``'noisy+truncated'`` adds noise first and truncates the noisy
        value afterwards; keyword options are broadcast to every stage that
        accepts them), or a structured specification: one stage entry
        ``{'name': ..., 'options': {...}}`` or an ordered list/tuple of such
        entries and bare names, each stage owning the options inside its entry
        (no keyword option is accepted then). An existing ``Feature`` is
        accepted as it is, without reparsing. The available stage names are
        ``'plain'``, ``'perturbed_x0'``, ``'noisy'``, ``'truncated'``,
        ``'permuted'``, ``'linearly_transformed'``, ``'random_nan'``,
        ``'unrelaxable_constraints'``, ``'nonquantifiable_constraints'``,
        ``'quantized'`` and ``'custom'``.

    Other Parameters
    ----------------
    distribution : str or callable, optional
        Distribution used by the 'noisy' feature ('gaussian' or 'uniform')
        and the 'perturbed_x0' feature ('spherical' or 'gaussian'), or a
        callable ``distribution(rng, size)``.
    noise_level, noise_type, noise_mode, noise_map : optional
        Options of the 'noisy' feature; ``noise_map`` is 'chebyshev' or a
        callable deterministic scalar map. ``noise_level`` is a finite,
        nonnegative real scalar.
    perturbation_level : float or array_like, optional
        Finite, nonnegative amplitude for 'perturbed_x0': a scalar or a
        one-dimensional vector of length 1 or problem.n (coordinatewise).
        Zero leaves the initial point unchanged; NaN and infinity are invalid.
    significant_digits, perturbed_trailing_digits : optional
        Options of the 'truncated' feature.
    nan_rate : int or float, optional
        Finite real scalar in [0, 1] for the 'random_nan' feature.
    rotated, condition_factor : optional
        Options of 'linearly_transformed'; ``condition_factor`` is a finite,
        nonnegative real scalar.
    unrelaxable_bounds, unrelaxable_linear_constraints, unrelaxable_nonlinear_constraints : bool, optional
        Options of the 'unrelaxable_constraints' feature.
    mesh_size, mesh_type, ground_truth : optional
        Options of 'quantized'. ``mesh_size`` is a finite, positive real scalar;
        ``mesh_type`` is 'absolute' or 'relative' (case-insensitive, stored in
        lowercase). Booleans are not numeric magnitude/count options. NumPy
        real scalar magnitudes are accepted. Only ``perturbation_level`` also
        accepts vectors; other magnitude options must be scalars. Option names
        are case-insensitive; two spellings of one name in the same call or
        stage entry are rejected as a duplicate.
    mod_x0, mod_affine, mod_bounds, mod_linear_ub, mod_linear_eq, mod_fun, mod_cub, mod_ceq : callable, optional
        Callbacks of the 'custom' feature.

    Attributes
    ----------
    stages : tuple
        The effective stages in order. Each record has ``name``, ``occurrence``,
        ``identity`` (``'noisy#0'``), ``code`` and a read-only ``options`` view
        of its validated, defaulted stage-local options.
    name : str
        Effective name (``'noisy+truncated'``; ``'plain'`` for the identity).
    declared_name, declared : str or None, Declaration
        The declaration as given (declaration route and entries before
        defaults). A specification imported from a historical object that
        recorded no declaration has route ``None``, no entries and
        ``declared_name`` ``None``; its effective stages are complete.
    is_stochastic, is_identity : bool

    Raises
    ------
    TypeError
        If an argument received an invalid value.
    ValueError
        If the arguments are inconsistent, or if an experiment option such as
        ``n_runs`` is given here instead of to ``benchmark``.

    Notes
    -----
    ``options`` and the ``modifier_*`` methods remain as deprecated one-stage
    conveniences only; the benchmark, the recorder, the archives and the
    reports use ``stages``.

    Ownership of option values: exact built-in data containers (``list``,
    ``tuple``, ``dict``) and exact ``numpy.ndarray`` values are copied at
    construction (arrays as read-only copies) and the inspection views hand
    out isolated read-only array copies, fresh container copies and read-only
    mappings, so changing such a container or array after construction, or
    changing what a view returned, never changes the specification. Opaque
    values are retained by identity and are not copied: callables, container
    subclasses, any other object and the elements of object-dtype arrays; the
    specification never mutates them, but their own state stays the caller's.
    No copy, iteration or reduction hook of a user object is invoked.

    Examples
    --------
    .. code-block:: python

        from optiprofiler import Feature

        feature = Feature('noisy', noise_level=1e-2)
        feature.stages[0].options['noise_level']   # 0.01

        pipeline = Feature([
            {'name': 'noisy', 'options': {'distribution': 'uniform'}},
            {'name': 'perturbed_x0', 'options': {'distribution': 'gaussian'}},
        ])
        [stage.identity for stage in pipeline.stages]   # ['noisy#0', 'perturbed_x0#0']
    """

    def __init__(self, name, **options):
        lowered = fold_option_names(options, 'Feature options')
        if isinstance(name, Feature):
            reject_experiment_options(lowered)
            reject_flat_stage_options(lowered)
            declared, stages = name._declared, name._stages
        elif isinstance(name, LegacyObject):
            raise TypeError(f'{type(name).__name__} is a historical feature object decoded by the trusted loader; '
                            'convert it with optiprofiler.legacy_compat.import_legacy_feature(...) and pass the '
                            'returned feature and n_runs to benchmark.')
        elif isinstance(name, str):
            declared, stages = normalize_shorthand(name, lowered)
        elif name is None:
            raise TypeError(_SPEC_TYPE_MESSAGE)
        else:
            reject_experiment_options(lowered)
            reject_flat_stage_options(lowered)
            declared, stages = normalize_entries(name)
        object.__setattr__(self, '_declared', declared)
        object.__setattr__(self, '_stages', stages)

    def __setattr__(self, key, value):
        raise AttributeError('Feature is immutable; build a new Feature instead.')

    def __delattr__(self, key):
        raise AttributeError('Feature is immutable; build a new Feature instead.')

    def __reduce__(self):
        # The native form is versioned: the canonical effective stages with all
        # validated local options explicit, plus the declaration as provenance.
        # No run count and no runtime state; unpickling validates again.
        effective = tuple((stage.name, stage.native_options()) for stage in self._stages)
        return _rebuild_feature, (FEATURE_NATIVE_VERSION, self._declared.route, self._declared.native_entries(), effective)

    def __setstate__(self, state):
        # Current pickles rebuild through ``__reduce__`` and never carry a
        # state dictionary; one only arrives when plain pickle restores an
        # OptiProfiler 1.x layout. Refuse it here instead of returning a
        # half-built object whose first attribute access fails.
        raise TypeError('This pickle holds an OptiProfiler 1.x Feature layout, which plain pickle cannot restore; read '
                        'it with optiprofiler.legacy_compat.load_trusted (or loads_trusted) and convert it with '
                        'optiprofiler.legacy_compat.import_legacy_feature.')

    def __repr__(self):
        return f'Feature({self.name!r})'

    @property
    def stages(self):
        """The effective stages, in order (empty for the identity)."""
        return self._stages

    @property
    def name(self):
        """The effective name; ``'plain'`` for the identity."""
        return '+'.join(stage.name for stage in self._stages) or FeatureName.PLAIN.value

    @property
    def declared_name(self):
        """The declared name, ``plain`` tokens included."""
        return self._declared.name

    @property
    def declared(self):
        """The declaration as given: input route and entries before defaults."""
        return self._declared

    @property
    def is_identity(self):
        return not self._stages

    @property
    def is_stochastic(self):
        """Whether any stage draws random numbers."""
        return any(stage.is_stochastic for stage in self._stages)

    def _single_stage_runtime(self):
        """Fresh runtime of the only stage (identity: a plain runtime); compositions have none."""
        if len(self._stages) > 1:
            raise ValueError('This feature has several stages; use `stages` (a composition has no single modifier).')
        if self._stages:
            return _StageRuntime(self._stages[0].name, self._stages[0].options)
        return _StageRuntime(FeatureName.PLAIN.value, {})

    def _deprecated_runtime(self, what):
        if len(self._stages) > 1:
            raise ValueError(f'`Feature.{what}` is available for one-stage features only; use `stages` for a composition.')
        warnings.warn(f'`Feature.{what}` is a deprecated single-stage convenience; the specification is `Feature.stages` '
                      f'(OptiProfiler 2.0).', DeprecationWarning, stacklevel=3)
        return self._single_stage_runtime()

    @property
    def options(self):
        """
        Deprecated: the stage-local options of a one-stage feature (``{}`` for
        the identity). Use ``stages``; a composition raises ``ValueError``.
        """
        runtime = self._deprecated_runtime('options')
        return dict(runtime._options)

    def modifier_x0(self, seed, problem):
        """Deprecated one-stage convenience; see ``stages``."""
        return self._deprecated_runtime('modifier_x0').modifier_x0(seed, problem)

    def modifier_affine(self, seed, problem):
        """Deprecated one-stage convenience; see ``stages``."""
        return self._deprecated_runtime('modifier_affine').modifier_affine(seed, problem)

    def modifier_bounds(self, seed, problem):
        """Deprecated one-stage convenience; see ``stages``."""
        return self._deprecated_runtime('modifier_bounds').modifier_bounds(seed, problem)

    def modifier_linear_ub(self, seed, problem):
        """Deprecated one-stage convenience; see ``stages``."""
        return self._deprecated_runtime('modifier_linear_ub').modifier_linear_ub(seed, problem)

    def modifier_linear_eq(self, seed, problem):
        """Deprecated one-stage convenience; see ``stages``."""
        return self._deprecated_runtime('modifier_linear_eq').modifier_linear_eq(seed, problem)

    def modifier_fun(self, x, seed, problem, n_eval):
        """Deprecated one-stage convenience; see ``stages``."""
        return self._deprecated_runtime('modifier_fun').modifier_fun(x, seed, problem, n_eval)

    def modifier_cub(self, x, seed, problem, n_eval_cub):
        """Deprecated one-stage convenience; see ``stages``."""
        return self._deprecated_runtime('modifier_cub').modifier_cub(x, seed, problem, n_eval_cub)

    def modifier_ceq(self, x, seed, problem, n_eval_ceq):
        """Deprecated one-stage convenience; see ``stages``."""
        return self._deprecated_runtime('modifier_ceq').modifier_ceq(x, seed, problem, n_eval_ceq)

    # Shared numerical utilities of the stage kinds (stateless).
    chebyshev_noise_map = staticmethod(_StageRuntime.chebyshev_noise_map)
    get_default_rng = staticmethod(_StageRuntime.get_default_rng)


def _reference_merit(value):
    """
    The ``merit`` of a reference record as a finite ``float``.

    Booleans, strings, arrays with a dimension, complex and object values are
    rejected with ``TypeError``; NaN, infinities and numbers that a float
    cannot hold exactly are rejected with ``ValueError``. The stored scalar is
    therefore exactly the scalar that was given: nothing is rounded, clipped
    or replaced.
    """
    message = 'The field `merit` of a problem reference must be a real scalar.'
    if isinstance(value, (bool, np.bool_, str, bytes)):
        raise TypeError(message)
    try:
        array = np.asarray(value)
    except Exception as exc:
        raise TypeError(message) from exc
    if array.ndim != 0 or array.dtype.kind not in 'iuf':
        raise TypeError(message)
    merit = float(array)
    if not np.isfinite(merit):
        raise ValueError('The field `merit` of a problem reference must be finite.')
    # An integer beyond 2**53 or an extended-precision value would be stored
    # as a different number; that is a silent reinterpretation, so reject it.
    exact = int(array) == int(merit) if array.dtype.kind in 'iu' else bool(array.dtype.type(merit) == array)
    if not exact:
        raise ValueError('The field `merit` of a problem reference cannot be represented exactly as a float.')
    return merit


def _restore_problem_reference(record):
    """
    Rebuild a serialized `ProblemReference` (the callable named by its pickle).

    The record is validated again under the contract of the running version.
    A record this version does not accept (an unknown mapping token written by
    another version, a non-finite merit, a foreign field) is read as *unknown*:
    ``None`` is returned with a ``RuntimeWarning``. It is never reinterpreted,
    for example by guessing a mapping or by reading another field as the merit.
    """
    try:
        if not isinstance(record, Mapping):
            raise TypeError('the serialized record is not a mapping')
        return ProblemReference.from_record(record)
    except (TypeError, ValueError) as exc:
        warnings.warn('A serialized problem reference is not valid under the reference contract of this version; it is '
                      f'read as unknown and is not reinterpreted ({exc}).', RuntimeWarning, stacklevel=2)
        return None


class ProblemReference:
    """
    Feasible reference fact of a `Problem`: one scalar stated by the author
    or provider of the problem, with its kind, its provenance and the mapping
    that fixes how the scalar is read.

    The record has exactly four fields and holds no point, no callable and no
    constraint violation.

    Parameters
    ----------
    merit : float
        The reference value. Must be a finite real scalar.
    kind : {'lower_bound', 'optimum', 'best_known', 'target'}
        The claim carried by ``merit``. Every kind is a claim over the
        *feasible* points of the problem:

        - ``'lower_bound'``: ``f(x) >= merit`` for every feasible ``x``.
        - ``'optimum'``: ``merit`` is the exact optimal value of the objective
          over the feasible set (an attained lower bound).
        - ``'best_known'``: ``merit`` is the objective value of a known
          feasible point, hence an upper bound on the optimum.
        - ``'target'``: ``merit`` is a level of the objective that the author
          chose as a target for feasible points; no mathematical claim.
    source : str
        Non-empty provenance of the fact (for example ``'author'`` or a
        citation). A reference without provenance is rejected.
    mapping : str
        Token of the closed registry `ProblemReference.MAPPINGS` that fixes
        how ``merit`` is read. The only token is ``'feasible_objective/1'``:
        ``merit`` is an objective value over feasible points, so it equals the
        merit of a feasible point under every merit function with the
        *feasible identity* ``merit_fun(f, 0, maxcv_init) == f`` for every
        ``maxcv_init``. The default merit function has this identity. A
        mapping is never a callable and users cannot register one; an unknown
        token is rejected, not guessed.

    Notes
    -----
    What the record is not:

    - It is not a run-history minimum and not the dynamic cohort minimum of a
      benchmark. The profile baseline is the least merit observed over the
      selected solvers, runs and evaluations; it changes with the cohort, is
      recomputed on every load and is never stored in a problem. Nothing
      derives a reference from solver output.
    - It is not a floor for run merits. On a constrained problem the merit of
      a run may be *below* the reference, because a merit function tolerates
      or penalizes small violations: an infeasible point can have a lower
      objective value than every feasible point. Such values are legitimate
      and are never clamped to the reference.
    - It is stated for feasible points only. Before a consumer compares run
      merits with the record under a custom ``merit_fun``, that function must
      be known to preserve the feasible identity above; otherwise the scalar
      and the run merits are not in the same space.

    The record is immutable and validated structurally. Validation never
    evaluates the objective or the constraints, so a reference adds no
    callback call to building or loading a problem, and a malformed record is
    rejected before any callback of the problem is touched. An omitted
    reference means unknown, not a bound at the initial point.

    See Also
    --------
    Problem : Optimization problem carrying an optional reference.
    FeaturedProblem : Which features retain the reference.
    """

    #: The four kinds; each is a claim over feasible points.
    KINDS = ('lower_bound', 'optimum', 'best_known', 'target')
    #: Closed registry of mapping tokens. A token names one fixed, versioned
    #: reading of ``merit``; a changed reading gets a new token, never a new
    #: meaning for an old one.
    MAPPINGS = ('feasible_objective/1',)
    #: The fields of the record, in order.
    FIELDS = ('merit', 'kind', 'source', 'mapping')
    # Fields of the superseded record layout (objective value, violation and
    # point). They are named so that the rejection can say why.
    _LEGACY_FIELDS = ('fun', 'maxcv', 'point')

    __slots__ = ('_merit', '_kind', '_source', '_mapping')

    def __init__(self, merit, kind, source, mapping):
        merit = _reference_merit(merit)
        if not isinstance(kind, str):
            raise TypeError('The field `kind` of a problem reference must be a string.')
        if kind not in self.KINDS:
            raise ValueError(f'The field `kind` of a problem reference must be one of {self.KINDS}, not {str(kind)!r}.')
        if not isinstance(source, str):
            raise TypeError('The field `source` of a problem reference must be a string.')
        if not source.strip():
            raise ValueError('The field `source` of a problem reference must be a non-empty provenance string.')
        if callable(mapping):
            raise TypeError('The field `mapping` of a problem reference must be a token of the closed registry '
                            f'{self.MAPPINGS}; callables and user-defined mappings are not accepted.')
        if not isinstance(mapping, str):
            raise TypeError('The field `mapping` of a problem reference must be a string token of the closed registry '
                            f'{self.MAPPINGS}.')
        if mapping not in self.MAPPINGS:
            raise ValueError(f'Unknown problem reference mapping {str(mapping)!r}; the closed registry is {self.MAPPINGS}. '
                             'An unknown mapping is rejected, never guessed.')
        object.__setattr__(self, '_merit', merit)
        object.__setattr__(self, '_kind', str(kind))
        object.__setattr__(self, '_source', str(source))
        object.__setattr__(self, '_mapping', str(mapping))

    @classmethod
    def from_record(cls, value):
        """
        Normalize a constructor input: a `ProblemReference` is returned as is
        and a mapping with exactly the four fields is converted. Anything
        else is rejected: a naked scalar (it has no kind, provenance or
        mapping), a record with a missing or unknown field, and a record of
        the superseded layout (``fun``, ``maxcv``, ``point``), which is never
        reinterpreted as a merit.
        """
        if isinstance(value, cls):
            return value
        if not isinstance(value, Mapping):
            raise TypeError('The argument `reference` for problem must be a ProblemReference or a mapping with the '
                            f'fields {cls.FIELDS}; a naked scalar has no kind, provenance or mapping and is rejected.')
        legacy = [key for key in cls._LEGACY_FIELDS if key in value]
        if legacy:
            raise ValueError(f'The problem reference field(s) {legacy} belong to a superseded record layout. Such a '
                             f'record is rejected and never reinterpreted; state the fields {cls.FIELDS}.')
        unknown = [key for key in value if key not in cls.FIELDS]
        if unknown:
            raise ValueError(f'Unknown problem reference field(s) {unknown!r}; the fields are {cls.FIELDS}.')
        missing = [key for key in cls.FIELDS if key not in value]
        if missing:
            raise ValueError(f'The problem reference is missing the required field(s) {missing}; all of {cls.FIELDS} '
                             'are required.')
        return cls(value['merit'], value['kind'], value['source'], value['mapping'])

    def __setattr__(self, key, value):
        raise AttributeError('ProblemReference is immutable.')

    def __delattr__(self, key):
        raise AttributeError('ProblemReference is immutable.')

    def __reduce__(self):
        # Serialized as the four native fields. Loading validates them again
        # and yields unknown (``None``) for a record this version rejects.
        return _restore_problem_reference, (self.as_dict(),)

    def __setstate__(self, state):
        # Never reached by a stream this class writes (see ``__reduce__``). A
        # stream that restores instance state directly was written for another
        # record layout, so its content is not read at all.
        raise TypeError('This serialized ProblemReference uses a superseded record layout; it is rejected and never '
                        f'reinterpreted. Build the record again with the fields {self.FIELDS}.')

    def __eq__(self, other):
        if not isinstance(other, ProblemReference):
            return NotImplemented
        return self._key() == other._key()

    def __ne__(self, other):
        result = self.__eq__(other)
        return result if result is NotImplemented else not result

    def __hash__(self):
        return hash(self._key())

    def __repr__(self):
        return (f'ProblemReference(merit={self._merit!r}, kind={self._kind!r}, source={self._source!r}, '
                f'mapping={self._mapping!r})')

    def _key(self):
        return self._merit, self._kind, self._source, self._mapping

    @property
    def merit(self):
        """The reference value, read as stated by ``mapping``."""
        return self._merit

    @property
    def kind(self):
        """The claim carried by ``merit``; always a claim over feasible points."""
        return self._kind

    @property
    def source(self):
        """Provenance of the record."""
        return self._source

    @property
    def mapping(self):
        """Registry token that fixes how ``merit`` is read."""
        return self._mapping

    def as_dict(self):
        """Native copy of the record: the four fields, in order."""
        return {'merit': self._merit, 'kind': self._kind, 'source': self._source, 'mapping': self._mapping}


def _propagated_reference(reference, stages):
    """
    The reference fact behind ``stages``, given as ``(name, options)`` pairs in
    order: ``reference`` itself if every stage retains it, otherwise ``None``.
    One stage that may change values, constraints, bounds or the truth makes
    the fact unknown for the whole pipeline; see
    `optiprofiler.feature_definitions.retains_reference`.
    """
    if reference is None or not all(_stage_retains_reference(name, options) for name, options in stages):
        return None
    return reference


class Problem:
    r"""
    Optimization problem to be used in the benchmarking.

    ``Problem`` describes an optimization problem with the following
    structure:

    .. math::

        \min \quad & \mathrm{fun}(x) \\
        \text{s.t.} \quad & x_l \le x \le x_u, \\
        & A_{\mathrm{ub}} x \le b_{\mathrm{ub}}, \\
        & A_{\mathrm{eq}} x = b_{\mathrm{eq}}, \\
        & c_{\mathrm{ub}}(x) \le 0, \\
        & c_{\mathrm{eq}}(x) = 0, \\
        & \text{with initial point } x_0,

    where ``fun`` is the objective function, ``x`` is the variable to
    optimize, ``xl`` and ``xu`` are the lower and upper bounds, ``aub`` and
    ``bub`` are the coefficient matrix and right-hand side vector of the
    linear inequality constraints, ``aeq`` and ``beq`` are the coefficient
    matrix and right-hand side vector of the linear equality constraints,
    ``cub`` is the function of nonlinear inequality constraints, and ``ceq``
    is the function of nonlinear equality constraints.

    Parameters
    ----------
    fun : callable
        Objective function to be minimized: ``fun(x) -> float``, where
        ``x`` is an array with shape ``(n,)``.
    x0 : array_like, shape (n,)
        Initial guess.
    name : str, optional
        Name of the problem. Default is ``'Unnamed Problem'``.
    xl : array_like, shape (n,), optional
        Lower bounds on the variables ``xl <= x``. Default is ``-numpy.inf``
        for each component.
    xu : array_like, shape (n,), optional
        Upper bounds on the variables ``x <= xu``. Default is ``numpy.inf``
        for each component.
    aub : array_like, shape (m_linear_ub, n), optional
        Coefficient matrix of the linear inequality constraints
        ``aub @ x <= bub``. Default is an empty matrix.
    bub : array_like, shape (m_linear_ub,), optional
        Right-hand side of the linear inequality constraints
        ``aub @ x <= bub``. Default is an empty vector.
    aeq : array_like, shape (m_linear_eq, n), optional
        Coefficient matrix of the linear equality constraints
        ``aeq @ x == beq``. Default is an empty matrix.
    beq : array_like, shape (m_linear_eq,), optional
        Right-hand side of the linear equality constraints
        ``aeq @ x == beq``. Default is an empty vector.
    cub : callable, optional
        Nonlinear inequality constraints ``cub(x) <= 0``:
        ``cub(x) -> array_like, shape (m_nonlinear_ub,)``. Default returns
        an empty array.
    ceq : callable, optional
        Nonlinear equality constraints ``ceq(x) == 0``:
        ``ceq(x) -> array_like, shape (m_nonlinear_eq,)``. Default returns
        an empty array.
    grad : callable, optional
        Gradient of the objective function: ``grad(x) -> array, shape (n,)``.
        Default returns an empty array.
    hess : callable, optional
        Hessian of the objective function:
        ``hess(x) -> array, shape (n, n)``. Default returns an empty matrix.
    jcub : callable, optional
        Jacobian of the nonlinear inequality constraints:
        ``jcub(x) -> array, shape (m_nonlinear_ub, n)``. The number of
        columns must equal ``n`` and the number of rows must equal
        ``m_nonlinear_ub``. Default returns an empty matrix.
    jceq : callable, optional
        Jacobian of the nonlinear equality constraints:
        ``jceq(x) -> array, shape (m_nonlinear_eq, n)``. Default returns
        an empty matrix.
    hcub : callable, optional
        Hessians of the nonlinear inequality constraints:
        ``hcub(x) -> list of arrays, each shape (n, n)``. The *i*-th
        element is the Hessian of the *i*-th constraint in ``cub``. Default
        returns an empty list.
    hceq : callable, optional
        Hessians of the nonlinear equality constraints:
        ``hceq(x) -> list of arrays, each shape (n, n)``. Default returns
        an empty list.
    reference : ProblemReference or dict, optional
        Optional feasible reference fact of the problem, stated by its author
        or provider: a `ProblemReference` or a dict with exactly the fields
        ``merit``, ``kind``, ``source`` and ``mapping``. A naked scalar, a
        record with other fields and an unknown mapping are rejected. Omitted
        means unknown. The record is validated structurally, first of all the
        arguments and without evaluating anything, so a reference adds no
        callback call to building or loading a problem.

    Attributes
    ----------
    name : str
        Name of the problem.
    x0 : numpy.ndarray, shape (n,)
        Initial guess.
    xl : numpy.ndarray, shape (n,)
        Lower bounds on the variables.
    xu : numpy.ndarray, shape (n,)
        Upper bounds on the variables.
    aub : numpy.ndarray, shape (m_linear_ub, n)
        Coefficient matrix of the linear inequality constraints.
    bub : numpy.ndarray, shape (m_linear_ub,)
        Right-hand side of the linear inequality constraints.
    aeq : numpy.ndarray, shape (m_linear_eq, n)
        Coefficient matrix of the linear equality constraints.
    beq : numpy.ndarray, shape (m_linear_eq,)
        Right-hand side of the linear equality constraints.
    ptype : str
        Type of the problem: ``'u'`` (unconstrained), ``'b'``
        (bound-constrained), ``'l'`` (linearly constrained), or ``'n'``
        (nonlinearly constrained).
    n : int
        Dimension of the problem (length of ``x``).
    mb : int
        Number of finite bound constraints.
    m_linear_ub : int
        Number of linear inequality constraints.
    m_linear_eq : int
        Number of linear equality constraints.
    m_nonlinear_ub : int
        Number of nonlinear inequality constraints.
    m_nonlinear_eq : int
        Number of nonlinear equality constraints.
    mlcon : int
        Total number of linear constraints (``m_linear_ub + m_linear_eq``).
    mnlcon : int
        Total number of nonlinear constraints
        (``m_nonlinear_ub + m_nonlinear_eq``).
    mcon : int
        Total number of constraints (``mlcon + mnlcon``).
    reference : ProblemReference or None
        Feasible reference fact of the problem, or ``None`` when unknown. It
        is author/provider metadata about the problem itself. It is not the
        profile baseline, which is the least merit observed over the selected
        solver histories, changes with the solver cohort and is never stored
        here; and it is not a floor for run merits, which may fall below it
        on a constrained problem (see `ProblemReference`).

    Methods
    -------
    fun(x)
        Evaluate the objective function.
    grad(x)
        Evaluate the gradient of the objective function.
    hess(x)
        Evaluate the Hessian of the objective function.
    cub(x)
        Evaluate the nonlinear inequality constraints.
    ceq(x)
        Evaluate the nonlinear equality constraints.
    jcub(x)
        Evaluate the Jacobian of the nonlinear inequality constraints.
    jceq(x)
        Evaluate the Jacobian of the nonlinear equality constraints.
    hcub(x)
        Evaluate the Hessians of the nonlinear inequality constraints.
    hceq(x)
        Evaluate the Hessians of the nonlinear equality constraints.
    maxcv(x)
        Compute the maximum constraint violation at ``x``, defined as the
        maximum of the infinity norms of ``max(xl - x, 0)``,
        ``max(x - xu, 0)``, ``max(aub @ x - bub, 0)``,
        ``aeq @ x - beq``, ``max(cub(x), 0)``, and ``ceq(x)``.
    project_x0()
        Attempt to project the initial guess ``x0`` onto the feasible
        region (may fail).

    See Also
    --------
    Feature : Feature applied to problems during benchmarking.
    FeaturedProblem : Problem equipped with a specific feature.
    benchmark : Main benchmarking function.

    Examples
    --------
    Consider the unconstrained problem of minimizing the Rosenbrock function

    .. math::

        f(x) = 100 (x_2 - x_1^2)^2 + (1 - x_1)^2.

    To create an instance of the class ``Problem`` for this problem:

    .. code-block:: python

        from optiprofiler import Problem

        def rosen(x):
            return 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2

        problem = Problem(rosen, [0, 0])

    The second argument ``[0, 0]`` is the initial guess. This instance can
    now be used to evaluate the objective function at any point and access
    extra information about the problem:

    .. code-block:: python

        problem.fun(problem.x0)  # returns 1.0
        problem.xl               # array([-inf, -inf])
        problem.xu               # array([inf, inf])

    The optional arguments of the constructor can be used to specify
    constraints. For example, to specify that the variables must be
    nonnegative:

    .. code-block:: python

        problem = Problem(rosen, [0, 0], xl=[0, 0])
        problem.xl  # array([0., 0.])

    Nonlinear inequality constraints can be specified using ``cub``. For
    example, to require :math:`x_1^2 + x_2^2 \le 1` and
    :math:`x_1^3 - x_2^2 \le 1`:

    .. code-block:: python

        def cub(x):
            return [x[0] ** 2 + x[1] ** 2 - 1, x[0] ** 3 - x[1] ** 2 - 1]

        problem = Problem(rosen, [0, 0], cub=cub)
        problem.cub(problem.x0)  # array([-1., -1.])

    The number of nonlinear inequality constraints is inferred
    automatically from the return value of ``cub`` at ``x0``. Nonlinear
    equality constraints can be specified in a similar way using ``ceq``.
    """

    def __init__(self, fun, x0, name=None, xl=None, xu=None, aub=None, bub=None, aeq=None, beq=None, cub=None, ceq=None, grad=None, hess=None, jcub=None, jceq=None, hcub=None, hceq=None, reference=None):
        """
        Initialize an optimization problem.

        Parameters
        ----------
        fun : callable
            Objective function to be minimized.

                ``fun(x) -> float``

            where ``x`` is an array with shape (n,).
        x0 : array_like, shape (n,)
            Initial guess.
        xl : array_like, shape (n,), optional
            Lower bounds on the variables ``xl <= x``.
        xu : array_like, shape (n,), optional
            Upper bounds on the variables ``x <= xu``.
        aub : array_like, shape (m_linear_ub, n), optional
            Coefficient matrix of the linear constraints ``aub @ x <= bub``.
        bub : array_like, shape (m_linear_ub,), optional
            Right-hand side of the linear constraints ``aub @ x <= bub``.
        aeq : array_like, shape (m_linear_eq, n), optional
            Coefficient matrix of the linear constraints ``aeq @ x == beq``.
        beq : array_like, shape (m_linear_eq,), optional
            Right-hand side of the linear constraints ``aeq @ x == beq``.
        cub : callable, optional
            Nonlinear inequality constraint ``cub(x) <= 0``.

                ``cub(x) -> array_like, shape (m_nonlinear_ub,)``

            where ``x`` is an array with shape (n,).
        ceq : callable, optional
            Nonlinear equality constraint ``ceq(x) == 0``.

                ``ceq(x) -> array_like, shape (m_nonlinear_eq,)``

            where ``x`` is an array with shape (n,).
        reference : `optiprofiler.opclasses.ProblemReference` or dict, optional
            Feasible reference fact of the problem: exactly the fields
            ``merit``, ``kind``, ``source`` and ``mapping``. Omitted means
            unknown. See the class documentation.

        Raises
        ------
        TypeError
            If an argument received an invalid value.
        ValueError
            If the arguments are inconsistent.
        """

        # Preprocess the optional feasible reference fact, before anything
        # else: the check is purely structural (four fields, a finite scalar,
        # a registry token) and evaluates nothing, so a malformed record is
        # rejected before any callback of the problem is touched. The
        # constraint callbacks are probed at `x0` further below, for their
        # dimensions; that happens with or without a reference, which
        # therefore adds no callback call.
        self._reference = None if reference is None else ProblemReference.from_record(reference)

        # Preprocess the objective function.
        self._fun = fun
        if not callable(self._fun):
            raise TypeError('The argument `fun` for problem must be callable.')

        # Preprocess the initial guess.
        self._x0 = _process_1d_array(x0, 'The argument `x0` for problem must be a one-dimensional array.')

        # Preprocess the name.
        self._name = name
        if self._name is not None and not isinstance(self._name, str):
            raise TypeError('The argument `name` for problem must be a string.')

        # Preprocess the bound constraints.
        self._xl = xl
        if self._xl is not None:
            self._xl = _process_1d_array(self._xl, 'The argument `xl` for problem must be a one-dimensional array.')
        self._xu = xu
        if self._xu is not None:
            self._xu = _process_1d_array(self._xu, 'The argument `xu` for problem must be a one-dimensional array.')

        # Preprocess the linear constraints.
        self._aub = aub
        if self._aub is not None:
            self._aub = _process_2d_array(self._aub, 'The argument `aub` for problem must be a two-dimensional array.')
        self._bub = bub
        if self._bub is not None:
            self._bub = _process_1d_array(self._bub, 'The argument `bub` for problem must be a one-dimensional array.')
        self._aeq = aeq
        if self._aeq is not None:
            self._aeq = _process_2d_array(self._aeq, 'The argument `aeq` for problem must be a two-dimensional array.')
        self._beq = beq
        if self._beq is not None:
            self._beq = _process_1d_array(self._beq, 'The argument `beq` for problem must be a one-dimensional array.')

        # Preprocess the nonlinear constraints.
        self._m_nonlinear_ub = 0
        self._cub = cub
        if self._cub is not None:
            if not callable(self._cub):
                raise TypeError('The argument `cub` for problem must be callable.')
            else:
                try:
                    c = _process_1d_array(self._cub(self._x0), 'The return value of the argument `cub` for problem must be a one-dimensional array.')
                    self._m_nonlinear_ub = c.size
                except Exception as err:
                    raise ValueError('Failed to determine the number of nonlinear inequality constraints at the initial guess.') from err
        self._m_nonlinear_eq = 0
        self._ceq = ceq
        if self._ceq is not None:
            if not callable(self._ceq):
                raise TypeError('The argument `ceq` for problem must be callable.')
            else:
                try:
                    c = _process_1d_array(self._ceq(self._x0), 'The return value of the argument `ceq` for problem must be a one-dimensional array.')
                    self._m_nonlinear_eq = c.size
                except Exception as err:
                    raise ValueError('Failed to determine the number of nonlinear equality constraints at the initial guess.') from err

        # Preprocess the gradient and the Hessian of the objective function.
        self._grad = grad
        if self._grad is not None:
            if not callable(self._grad):
                raise TypeError('The argument `grad` for problem must be callable.')
        self._hess = hess
        if self._hess is not None:
            if not callable(self._hess):
                raise TypeError('The argument `hess` for problem must be callable.')

        # Preprocess the Jacobian and the Hessian of the nonlinear constraints.
        self._jcub = jcub
        if self._jcub is not None:
            if not callable(self._jcub):
                raise TypeError('The argument `jcub` for problem must be callable.')
        self._jceq = jceq
        if self._jceq is not None:
            if not callable(self._jceq):
                raise TypeError('The argument `jceq` for problem must be callable.')
        self._hcub = hcub
        if self._hcub is not None:
            if not callable(self._hcub):
                raise TypeError('The argument `hcub` for problem must be callable.')
        self._hceq = hceq
        if self._hceq is not None:
            if not callable(self._hceq):
                raise TypeError('The argument `hceq` for problem must be callable.')

        # Check that the arguments are consistent.
        if self.xl.size != self.n:
            raise ValueError(f'The argument `xl` for problem must have size {self.n}.')
        if self.xu.size != self.n:
            raise ValueError(f'The argument `xu` for problem must have size {self.n}.')
        if self.aub.shape != (self.m_linear_ub, self.n):
            raise ValueError(f'The argument `aub` for problem must have shape {(self.m_linear_ub, self.n)}.')
        if self.aeq.shape != (self.m_linear_eq, self.n):
            raise ValueError(f'The argument `aeq` for problem must have shape {(self.m_linear_eq, self.n)}.')
        if self.bub.size != self.m_linear_ub:
            raise ValueError(f'The argument `bub` for problem must have size {self.m_linear_ub}.')
        if self.beq.size != self.m_linear_eq:
            raise ValueError(f'The argument `beq` for problem must have size {self.m_linear_eq}.')

    @property
    def n(self):
        """
        Dimension of the problem.

        Returns
        -------
        int
            Dimension of the problem.
        """
        return self.x0.size

    @property
    def mb(self):
        """
        Number of bound constraints.

        Returns
        -------
        int
            Number of bound constraints.
        """
        return sum(self.xl > -np.inf) + sum(self.xu < np.inf)

    @property
    def m_linear_ub(self):
        """
        Number of linear inequality constraints.

        Returns
        -------
        int
            Number of linear inequality constraints.
        """
        return sum(~np.isinf(self.bub))

    @property
    def m_linear_eq(self):
        """
        Number of linear equality constraints.

        Returns
        -------
        int
            Number of linear equality constraints.
        """
        return self.beq.size

    @property
    def m_nonlinear_ub(self):
        """
        Number of nonlinear inequality constraints.

        Returns
        -------
        int
            Number of nonlinear inequality constraints.
        """
        return self._m_nonlinear_ub

    @property
    def m_nonlinear_eq(self):
        """
        Number of nonlinear equality constraints.

        Returns
        -------
        int
            Number of nonlinear equality constraints.
        """
        return self._m_nonlinear_eq

    @property
    def mlcon(self):
        """
        Total number of linear constraints (inequality and equality).

        Returns
        -------
        int
            Total number of linear constraints.
        """
        return self.m_linear_ub + self.m_linear_eq

    @property
    def mnlcon(self):
        """
        Total number of nonlinear constraints (inequality and equality).

        Returns
        -------
        int
            Total number of nonlinear constraints.
        """
        return self.m_nonlinear_ub + self.m_nonlinear_eq

    @property
    def mcon(self):
        """
        Total number of constraints (linear and nonlinear).

        Returns
        -------
        int
            Total number of constraints.
        """
        return self.mlcon + self.mnlcon

    @property
    def ptype(self):
        """
        Type of the problem.

        Returns
        -------
        str
            Type of the problem.
        """
        try:
            if self.mnlcon > 0:
                return 'n'
            elif self.mlcon > 0:
                return 'l'
            elif self.mb > 0:
                return 'b'
            else:
                return 'u'
        except ValueError:
            return 'n'

    @property
    def name(self):
        """
        Name of the problem.

        Returns
        -------
        str
            Name of the problem.
        """
        return self._name if self._name is not None else 'Unnamed Problem'

    @property
    def x0(self):
        """
        Initial guess.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Initial guess.
        """
        return np.copy(self._x0)

    @property
    def xl(self):
        """
        Lower bounds on the variables.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Lower bounds on the variables.
        """
        return np.copy(self._xl) if self._xl is not None else np.full(self.n, -np.inf)

    @property
    def xu(self):
        """
        Upper bounds on the variables.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Upper bounds on the variables.
        """
        return np.copy(self._xu) if self._xu is not None else np.full(self.n, np.inf)

    @property
    def aub(self):
        """
        Coefficient matrix of the linear constraints ``aub @ x <= bub``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_ub, n)
            Coefficient matrix of the linear inequality constraints.
        """
        return np.copy(self._aub) if self._aub is not None else np.empty((0, self.n))

    @property
    def bub(self):
        """
        Right-hand side of the linear constraints ``aub @ x <= bub``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_ub,)
            Right-hand side of the linear inequality constraints.
        """
        return np.copy(self._bub) if self._bub is not None else np.empty(0)

    @property
    def aeq(self):
        """
        Coefficient matrix of the linear constraints ``aeq @ x == beq``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_eq, n)
            Coefficient matrix of the linear equality constraints.
        """
        return np.copy(self._aeq) if self._aeq is not None else np.empty((0, self.n))

    @property
    def beq(self):
        """
        Right-hand side of the linear constraints ``aeq @ x == beq``.

        Returns
        -------
        `numpy.ndarray`, shape (m_linear_eq,)
            Right-hand side of the linear equality constraints.
        """
        return np.copy(self._beq) if self._beq is not None else np.empty(0)

    @property
    def reference(self):
        """
        Feasible reference fact of the problem.

        Returns
        -------
        `ProblemReference` or None
            The author/provider record, or ``None`` when unknown. An object
            serialized before the record existed, or holding anything that is
            not a validated record, reads as unknown.
        """
        # Only a validated, immutable record is ever handed out: whatever a
        # foreign or older serialization left in the instance reads as
        # unknown and is never reinterpreted.
        reference = getattr(self, '_reference', None)
        return reference if isinstance(reference, ProblemReference) else None

    def fun(self, x):
        """
        Evaluate the objective function.

        The optimization problem is to minimize the objective function.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the objective function.

        Returns
        -------
        float
            Value of the objective function at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `fun` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `fun` in problem must have size {self.n}.')
        try:
            f = self._fun(x)
            f = float(f)
        except Exception as err:
            logger = get_logger(__name__)
            logger.warning(f'Failed to evaluate the objective function: {shorten_log_message(err)}')
            f = np.nan
        return f

    def grad(self, x):
        """
        Evaluate the gradient of the objective function.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the gradient of the objective function.

        Returns
        -------
        `numpy.ndarray`, shape (n,)
            Gradient of the objective function at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `grad` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `grad` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `grad` in problem must have size {self.n}.')
        if self._grad is None:
            g = np.empty(0)
        else:
            try:
                g = self._grad(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the gradient of the objective function: {shorten_log_message(err)}')
                g = np.full(self.n, np.nan)
            g = _process_1d_array(g, 'The return value of the argument `grad` for problem must be a one-dimensional array.')
            if g.size != 0 and g.size != self.n:
                raise ValueError(f'The return value of the argument `grad` for problem must have size {self.n}.')
        return g

    def hess(self, x):
        """
        Evaluate the Hessian of the objective function.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the Hessian of the objective function.

        Returns
        -------
        `numpy.ndarray`, shape (n, n)
            Hessian of the objective function at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `hess` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `hess` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `hess` in problem must have size {self.n}.')
        if self._hess is None:
            h = np.empty((0, 0))
        else:
            try:
                h = self._hess(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the Hessian of the objective function: {shorten_log_message(err)}')
                h = np.full((self.n, self.n), np.nan)
            h = _process_2d_array(h, 'The return value of the argument `hess` for problem must be a two-dimensional array.')
            if h.size != 0 and h.shape != (self.n, self.n):
                raise ValueError(f'The return value of the argument `hess` for problem must have shape {(self.n, self.n)}.')
        return h

    def cub(self, x):
        """
        Evaluate the nonlinear constraints ``cub(x) <= 0``.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the nonlinear inequality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (m_nonlinear_ub,)
            Values of the nonlinear inequality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `cub` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `cub` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `cub` in problem must have size {self.n}.')
        if self._cub is None:
            c = np.empty(0)
        else:
            try:
                c = self._cub(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the nonlinear inequality constraint function: {shorten_log_message(err)}')
                c = np.full(self.m_nonlinear_ub, np.nan)
            c = _process_1d_array(c, 'The return value of the argument `cub` for problem must be a one-dimensional array.')
            if c.size != 0 and c.size != self.m_nonlinear_ub:
                raise ValueError(f'The return value of the argument `cub` for problem must have size {self.m_nonlinear_ub}.')
        return c

    def ceq(self, x):
        """
        Evaluate the nonlinear constraints ``ceq(x) == 0``.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the nonlinear equality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (m_nonlinear_eq,)
            Values of the nonlinear equality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `ceq` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `ceq` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `ceq` in problem must have size {self.n}.')
        if self._ceq is None:
            c = np.empty(0)
        else:
            try:
                c = self._ceq(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the nonlinear equality constraint function: {shorten_log_message(err)}')
                c = np.full(self.m_nonlinear_eq, np.nan)
            c = _process_1d_array(c, 'The return value of the argument `ceq` for problem must be a one-dimensional array.')
            if c.size != 0 and c.size != self.m_nonlinear_eq:
                raise ValueError(f'The return value of the argument `ceq` for problem must have size {self.m_nonlinear_eq}.')
        return c

    def jcub(self, x):
        """
        Evaluate the Jacobian of the nonlinear inequality constraints.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the Jacobian of the nonlinear inequality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (m_nonlinear_ub, n)
            Jacobian of the nonlinear inequality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `jcub` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `jcub` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `jcub` in problem must have size {self.n}.')
        if self._jcub is None:
            j = np.empty((0, 0))
        else:
            try:
                j = self._jcub(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the Jacobian of the nonlinear inequality constraint function: {shorten_log_message(err)}')
                j = np.full((self.m_nonlinear_ub, self.n), np.nan)
            j = _process_2d_array(j, 'The return value of the argument `jcub` for problem must be a two-dimensional array.')
            if j.size != 0 and j.shape != (self.m_nonlinear_ub, self.n):
                raise ValueError(f'The return value of the argument `jcub` for problem must have shape {(self.m_nonlinear_ub, self.n)}.')
        return j

    def jceq(self, x):
        """
        Evaluate the Jacobian of the nonlinear equality constraints.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the Jacobian of the nonlinear equality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (m_nonlinear_eq, n)
            Jacobian of the nonlinear equality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `jceq` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `jceq` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `jceq` in problem must have size {self.n}.')
        if self._jceq is None:
            j = np.empty((0, 0))
        else:
            try:
                j = self._jceq(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the Jacobian of the nonlinear equality constraint function: {shorten_log_message(err)}')
                j = np.full((self.m_nonlinear_eq, self.n), np.nan)
            j = _process_2d_array(j, 'The return value of the argument `jceq` for problem must be a two-dimensional array.')
            if j.size != 0 and j.shape != (self.m_nonlinear_eq, self.n):
                raise ValueError(f'The return value of the argument `jceq` for problem must have shape {(self.m_nonlinear_eq, self.n)}.')
        return j

    def hcub(self, x):
        """
        Evaluate the Hessian of the nonlinear inequality constraints.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the Hessian of the nonlinear inequality constraints.

        Returns
        -------
        `list` of `numpy.ndarray`, shape (m_nonlinear_ub,)
            List of Hessians of the nonlinear inequality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `hcub` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `hcub` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `hcub` in problem must have size {self.n}.')
        if self._hcub is None:
            h = []
        else:
            try:
                h = self._hcub(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the Hessian of the nonlinear inequality constraint function: {shorten_log_message(err)}')
                h = [np.full((self.n, self.n), np.nan)] * self.m_nonlinear_ub
            h = [_process_2d_array(h_i, 'Each element of the return value of the argument `hcub` for problem must be a two-dimensional array.') for h_i in h]
            if len(h) != self.m_nonlinear_ub:
                raise ValueError(f'The return value of the argument `hcub` for problem must have {self.m_nonlinear_ub} elements.')
            for h_i in h:
                if h_i.size != 0 and h_i.shape != (self.n, self.n):
                    raise ValueError(f'Each element of the return value of the argument `hcub` for problem must have shape {(self.n, self.n)}.')
        return h

    def hceq(self, x):
        """
        Evaluate the Hessian of the nonlinear equality constraints.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the Hessian of the nonlinear equality constraints.

        Returns
        -------
        `list` of `numpy.ndarray`, shape (m_nonlinear_eq,)
            List of Hessians of the nonlinear equality constraints at `x`.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `hceq` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `hceq` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `hceq` in problem must have size {self.n}.')
        if self._hceq is None:
            h = []
        else:
            try:
                h = self._hceq(x)
            except Exception as err:
                logger = get_logger(__name__)
                logger.warning(f'Failed to evaluate the Hessian of the nonlinear equality constraint function: {shorten_log_message(err)}')
                h = [np.full((self.n, self.n), np.nan)] * self.m_nonlinear_eq
            h = [_process_2d_array(h_i, 'Each element of the return value of the argument `hceq` for problem must be a two-dimensional array.') for h_i in h]
            if len(h) != self.m_nonlinear_eq:
                raise ValueError(f'The return value of the argument `hceq` for problem must have {self.m_nonlinear_eq} elements.')
            for h_i in h:
                if h_i.size != 0 and h_i.shape != (self.n, self.n):
                    raise ValueError(f'Each element of the return value of the argument `hceq` for problem must have shape {(self.n, self.n)}.')
        return h

    def _maxcv(self, x):
        """
        Evaluate the maximum constraint violations.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the maximum constraint violation.

        Returns
        -------
        float
            Maximum constraint violation.
        float
            Maximum constraint violation for the bound constraints.
        float
            Maximum constraint violation for the linear constraints.
        float
            Maximum constraint violation for the nonlinear constraints.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape.
        """
        x = _process_1d_array(x, 'The argument `x` for method `maxcv` in problem must be a one-dimensional array.')
        if x.size != self.n:
            raise ValueError(f'The argument `x` for method `maxcv` in problem must have size {self.n}.')
        
        cv = 0.0
        cv_bounds = 0.0
        cv_linear = 0.0
        cv_nonlinear = 0.0

        if self.ptype == 'u':
            return cv, cv_bounds, cv_linear, cv_nonlinear

        if np.any(np.isfinite(self.xl)):
            cv_bounds = np.max(self.xl - x, initial=0.0)
        if np.any(np.isfinite(self.xu)):
            cv_bounds = np.max(x - self.xu, initial=cv_bounds)
        if self.ptype == 'b':
            cv = cv_bounds
            return cv, cv_bounds, cv_linear, cv_nonlinear
        
        if self.aub.size > 0:
            cv_linear = np.max(self.aub @ x - self.bub, initial=0.0)
        if self.aeq.size > 0:
            cv_linear = np.max(np.abs(self.aeq @ x - self.beq), initial=cv_linear)
        if self.ptype == 'l':
            cv = np.maximum(cv_bounds, cv_linear)
            return cv, cv_bounds, cv_linear, cv_nonlinear
        
        if self.m_nonlinear_ub > 0:
            cv_nonlinear = np.max(self.cub(x), initial=0.0)
        if self.m_nonlinear_eq > 0:
            cv_nonlinear = np.max(np.abs(self.ceq(x)), initial=cv_nonlinear)
        cv = np.max([cv_bounds, cv_linear, cv_nonlinear])
        return cv, cv_bounds, cv_linear, cv_nonlinear

    def maxcv(self, x):
        """
        Evaluate the maximum constraint violation.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the maximum constraint violation.

        Returns
        -------
        float
            Maximum constraint violation.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape.
        """
        return self._maxcv(x)[0]

    def project_x0(self):
        """
        Project the initial guess onto the feasible region.
        """
        if self.ptype == 'b':
            self._x0 = np.clip(self._x0, self.xl, self.xu)
        elif self.ptype == 'l' and self.m_linear_ub == 0 and np.all(self.xl == -np.inf) and np.all(self.xu == np.inf):
            self._x0 = self.x0 + lstsq(self.aeq, self.beq - self.aeq @ self.x0)[0]
        elif self.ptype != 'u':
            bounds = Bounds(self.xl, self.xu)
            constraints = []
            if self.m_linear_ub > 0:
                constraints.append(LinearConstraint(self.aub, -np.inf, self.bub))
            if self.m_linear_eq > 0:
                constraints.append(LinearConstraint(self.aeq, self.beq, self.beq))
            if self.m_nonlinear_ub > 0:
                constraints.append(NonlinearConstraint(self.cub, -np.inf, np.zeros(self.m_nonlinear_ub)))
            if self.m_nonlinear_eq > 0:
                constraints.append(NonlinearConstraint(self.ceq, np.zeros(self.m_nonlinear_eq), np.zeros(self.m_nonlinear_eq)))

            def dist_x0_sq(x):
                g = x - self.x0
                return 0.5 * (g @ g), g
            def hessp(x, p):
                return p

            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                minimize_options = {}
                # In SciPy < 1.16, constrained ``minimize`` chooses SLSQP
                # when ``method`` is not specified.  On Python 3.10, SciPy
                # cannot be upgraded past 1.15.x, and the old f2py/Fortran
                # SLSQP wrapper can corrupt native memory for overdetermined
                # equality systems (``m_eq > n``).  The solver still reports
                # status 2 ("More equality constraints than independent
                # variables"), but the process may later abort or segfault
                # during garbage collection.  This is exactly the failure
                # chain observed for PyCUTEst problem JENSMPNE
                # (n = 2, m_eq = 10) when ``project_x0=True``.  Keep the
                # workaround narrow: only this Python/SciPy/version/shape
                # combination avoids the unsafe default SLSQP path.  In the
                # same overdetermined equality case, ``trust-constr`` also
                # needs the SVD projection factorization; otherwise SciPy can
                # fail before returning a projected point with
                # "expected square matrix".
                if (
                    sys.version_info[:2] == (3, 10)
                    and _scipy_version_less_than(1, 16)
                    and self.m_linear_eq + self.m_nonlinear_eq > self.n
                ):
                    minimize_options['method'] = 'trust-constr'
                    minimize_options['options'] = {'factorization_method': 'SVDFactorization'}
                res = minimize(
                    dist_x0_sq,
                    self.x0,
                    jac=True,
                    hessp=hessp,
                    bounds=bounds,
                    constraints=constraints,
                    tol=1e-16,
                    **minimize_options,
                )
            self._x0 = res.x

class FeaturedProblem(Problem):
    r"""
    Subclass of `Problem` that equips an optimization problem with a feature.

    ``Problem`` and its subclass ``FeaturedProblem`` describe the following
    optimization problem:

    .. math::

        \min \quad & \mathrm{fun}(x) \\
        \text{s.t.} \quad & x_l \le x \le x_u, \\
        & A_{\mathrm{ub}} x \le b_{\mathrm{ub}}, \\
        & A_{\mathrm{eq}} x = b_{\mathrm{eq}}, \\
        & c_{\mathrm{ub}}(x) \le 0, \\
        & c_{\mathrm{eq}}(x) = 0, \\
        & \text{with initial point } x_0.

    Parameters
    ----------
    problem : Problem
        The original optimization problem.
    feature : Feature
        The feature to apply to the problem.
    max_eval : int
        Maximum number of function evaluations.
    seed : int, optional
        Nonnegative integer seed for the random number generator.

    Attributes
    ----------
    problem : Problem
        The original (unmodified) optimization problem.
    feature : Feature
        The feature applied to the optimization problem.
    max_eval : int
        Maximum number of function evaluations.
    seed : int
        Seed for the random number generator.
    fun_hist : list of float
        History of evaluated objective function values.
    cub_hist : list of numpy.ndarray
        History of evaluated nonlinear inequality constraint values.
    ceq_hist : list of numpy.ndarray
        History of evaluated nonlinear equality constraint values.
    maxcv_hist : list of float
        History of maximum constraint violations.
    n_eval_fun : int
        Minimum of the number of objective function evaluations and
        ``max_eval``.
    n_eval_cub : int
        Minimum of the number of nonlinear inequality constraint
        evaluations and ``max_eval``.
    n_eval_ceq : int
        Minimum of the number of nonlinear equality constraint evaluations
        and ``max_eval``.
    fun_init : float
        Objective function value at the initial point.
    maxcv_init : float
        Maximum constraint violation at the initial point.
    reference : ProblemReference or None
        The original problem's feasible reference fact if the feature retains
        it, otherwise ``None`` (unknown). The record is retained unchanged or
        dropped; nothing is transported or derived, and the rule is the same
        for every kind. Retaining stages: those that change only observations
        (``noisy``, ``truncated``, ``random_nan``,
        ``nonquantifiable_constraints``, ``unrelaxable_constraints``, and
        ``quantized`` with ``ground_truth=False``), ``perturbed_x0``,
        ``permuted`` and ``linearly_transformed``. ``custom`` retains it only
        if its options are a subset of ``mod_x0`` and ``mod_affine``; any
        other custom option (``mod_fun``, ``mod_cub``, ``mod_ceq``,
        ``mod_bounds``, ``mod_linear_ub``, ``mod_linear_eq``) may change
        values, constraints or bounds and makes it unknown. ``quantized`` with
        ``ground_truth=True`` (the default) makes it unknown, because the
        truth is then the mesh problem. A composition retains the record only
        if every stage does.

    Notes
    -----
    ``FeaturedProblem`` inherits all methods of ``Problem``, but the
    methods ``fun``, ``cub``, ``ceq``, and ``maxcv`` are modified by the
    input ``Feature``.

    1. When the number of function evaluations reaches ``max_eval``, the
       methods ``fun``, ``cub``, and ``ceq`` will return the values at the
       point where the maximum number of function evaluations was reached.
    2. When the number of function evaluations reaches
       ``termination_eval`` (used internally by ``benchmark``), the
       methods ``fun``, ``cub``, and ``ceq`` will raise an error to
       terminate the optimization process.

    .. note::

        For consistency with the rest of OptiProfiler, we recommend
        defining callables (such as ``fun``) with ``def`` rather than
        ``lambda``. Lambda expressions are not picklable, which prevents
        parallel execution when such callables are eventually passed to
        :func:`~optiprofiler.benchmark` with ``n_jobs > 1``. See
        :ref:`py_callable_picklability` for details.

    See Also
    --------
    Problem : Optimization problem.
    Feature : Feature applied to problems during benchmarking.
    benchmark : Main benchmarking function.

    Examples
    --------
    Create a featured problem with a noisy feature:

    .. code-block:: python

        import numpy as np
        from optiprofiler import Problem, Feature, FeaturedProblem

        def fun(x):
            return np.sum(x ** 2)

        problem = Problem(fun, [1.0, 2.0])
        feature = Feature('noisy')
        fp = FeaturedProblem(problem, feature, max_eval=100, seed=0)
    """

    def __init__(self, problem, feature, max_eval, seed=None):
        """
        Initialize an optimization problem.

        Parameters
        ----------
        problem : `optiprofiler.opclasses.Problem`
            Problem to be used in the benchmarking.
        feature : `optiprofiler.opclasses.Feature`
            Feature to be used in the benchmarking.
        max_eval : int
            Maximum number of function evaluations.
        seed : int, optional
            Seed for the random number generator.

        Raises
        ------
        TypeError
            If an argument received an invalid value.
        """

        self._problem = problem
        # Preprocess the feature.
        self._feature = feature
        if not isinstance(self._feature, Feature):
            raise TypeError('The argument `feature` for featured problem must be an instance of the class Feature.')
        # Explicit execution strategy (recorded in provenance) and a fresh
        # runtime of the only stage for this trial; the specification itself
        # never carries runtime state.
        self.execution_strategy = select_execution_strategy(feature)
        self._runtime = feature._single_stage_runtime()

        # Preprocess the maximum number of function evaluations.
        self._max_eval = _validate_max_eval(max_eval)

        # Preprocess the seed.
        self._seed = _validate_seed(seed)

        # Record the real evaluation numbers.
        self._real_n_eval_fun = 0
        self._real_n_eval_cub = 0
        self._real_n_eval_ceq = 0

        # Modify the problem according to the feature.
        self._x0 = self._runtime.modifier_x0(self._seed, self._problem)
        self._xl, self._xu = self._runtime.modifier_bounds(self._seed, self._problem)
        self._aub, self._bub = self._runtime.modifier_linear_ub(self._seed, self._problem)
        self._aeq, self._beq = self._runtime.modifier_linear_eq(self._seed, self._problem)

        # The reference fact is never copied blindly (``__new__`` copied the
        # original problem's attributes): the stage either retains the record
        # unchanged or makes it unknown. Only the stage name and its option
        # names are read, so no callback is called for this.
        self._reference = _propagated_reference(self._problem.reference, [(self._runtime.name, self._runtime._options)])

        # Store the histories of the objective function values, nonlinear
        # constraints, and maximum constraint violations.
        self._fun_hist = []
        self._cub_hist = []
        self._ceq_hist = []
        self._maxcv_hist = []
        self._last_fun = np.nan
        self._last_cub = np.nan
        self._last_ceq = np.nan

        self._fun_init, self._maxcv_init = self._evaluate_truth(self._x0)

    def _evaluate_truth(self, x):
        """Evaluate scoring truth at solver coordinates, without recording an oracle call."""
        A, b = self._runtime.modifier_affine(self._seed, self._problem)[:2]
        # True means the quantized problem itself is the ground truth. False
        # (and all other features) keeps base truth. Never snap the returned
        # solver point or call self.fun here: bookkeeping must not consume its
        # budget, append a synthetic history entry, or update last-value caches.
        if self._runtime.name == 'quantized' and self._runtime.options[FeatureOption.GROUND_TRUTH]:
            f = self._runtime.modifier_fun(A @ x + b, self._seed, self._problem, self.n_eval_fun)
            cv = self.maxcv(x)
        else:
            f = self._problem.fun(A @ x + b)
            cv = self._problem.maxcv(A @ x + b)
        return f, cv


    def __new__(cls, problem, feature, max_eval, seed=None):
        # Preprocess the problem.
        if not isinstance(problem, Problem):
            raise TypeError('The argument `problem` for featured problem must be an instance of the class Problem.')

        # Single dispatch boundary: a feature with at least two effective
        # stages is applied through lazy problem views and the recorder of
        # ``optiprofiler.composition``; everything else keeps this class.
        if cls is FeaturedProblem and isinstance(feature, Feature) and select_execution_strategy(feature) == STRATEGY_COMPOSED:
            from .composition import ComposedFeaturedProblem
            return object.__new__(ComposedFeaturedProblem)

        # Create a new instance of the class `FeaturedProblem` by copying the
        # attributes of the problem passed to the __init__ method.
        instance = super().__new__(cls)
        for k, v in problem.__dict__.items():
            instance.__dict__[k] = v
        return instance

    @property
    def fun_init(self):
        """
        Objective function value at the initial point.

        Returns
        -------
        float
            Objective function value at the initial point.
        """
        return self._fun_init
    
    @property
    def maxcv_init(self):
        """
        Maximum constraint violation at the initial point.

        Returns
        -------
        float
            Maximum constraint violation at the initial point.
        """
        return self._maxcv_init

    @property
    def n_eval_fun(self):
        """
        Number of objective function evaluations.

        Returns
        -------
        int
            Number of objective function evaluations.
        """
        return len(self._fun_hist)
    
    @property
    def n_eval_cub(self):
        """
        Number of nonlinear inequality constraint evaluations.

        Returns
        -------
        int
            Number of nonlinear inequality constraint evaluations.
        """
        return len(self._cub_hist)
    
    @property
    def n_eval_ceq(self):
        """
        Number of nonlinear equality constraint evaluations.

        Returns
        -------
        int
            Number of nonlinear equality constraint evaluations.
        """
        return len(self._ceq_hist)

    @property
    def fun_hist(self):
        """
        History of objective function values.

        Returns
        -------
        `numpy.ndarray`, shape (n_eval_fun,)
            History of objective function values.
        """
        return np.array(self._fun_hist)

    @property
    def cub_hist(self):
        """
        History of nonlinear inequality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (n_eval_cub, m_nonlinear_ub)
            History of nonlinear inequality constraints.
        """
        return np.array(self._cub_hist)

    @property
    def ceq_hist(self):
        """
        History of nonlinear equality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (n_eval_ceq, m_nonlinear_eq)
            History of nonlinear equality constraints.
        """
        return np.array(self._ceq_hist)

    @property
    def maxcv_hist(self):
        """
        History of maximum constraint violations.

        Returns
        -------
        `numpy.ndarray`, shape (n_eval_fun,)
            History of maximum constraint violations.
        """
        return np.array(self._maxcv_hist)

    def fun(self, x):
        """
        Evaluate the objective function.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the objective function.

        Returns
        -------
        float
            Value of the objective function at `x` if maximum number of evaluations has not been reached,
            otherwise returns the last evaluated objective function value.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape.
        StopIteration
            If the number of the objective function evaluations has reached two times the maximum function evaluations.
        """
        if self._real_n_eval_fun >= 2 * self._max_eval:
            raise StopIteration(f'The number of the objective function evaluations has reached {2 * self._max_eval} (two times the maximum function evaluations).')
        self._real_n_eval_fun += 1

        if self.n_eval_fun >= self._max_eval:
            # If the maximum number of evaluations has been reached, return
            # the last evaluated objective function value.
            return self._last_fun

        # Generate the affine transformation.
        A, b = self._runtime.modifier_affine(self._seed, self._problem)[:2]

        # Evaluate the modified the objective function value according to the feature and return the
        # modified value.
        f = self._runtime.modifier_fun(A @ x + b, self._seed, self._problem, self.n_eval_fun)
        self._last_fun = f

        # Evaluate the objective function and store the results.
        f_true = super().fun(A @ x + b)

        # If the feature is 'quantized' and the option ``ground_truth'' is set to true, we should
        # set f_true to f.
        if self._runtime.name == 'quantized' and self._runtime.options[FeatureOption.GROUND_TRUTH]:
            f_true = f

        # We should not store the modified value because the performance of an optimization solver
        # should be measured using the original objective function.
        self._fun_hist.append(f_true)
        try:
            self._maxcv_hist.append(self.maxcv(x))
        except Exception:
            self._maxcv_hist.append(np.nan)

        return f

    def cub(self, x, record_hist=True):
        """
        Evaluate the nonlinear constraints ``cub(x) <= 0``.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the nonlinear inequality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (m_nonlinear_ub,)
            Values of the nonlinear inequality constraints at `x` if maximum number of evaluations has not been reached,
            otherwise returns the last evaluated nonlinear inequality constraints.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `cub` has an invalid shape.
        StopIteration
            If the number of the nonlinear inequality constraint evaluations has reached two times the maximum function evaluations.
        """
        if self._real_n_eval_cub >= 2 * self._max_eval:
            raise StopIteration(f'The number of the nonlinear inequality constraint evaluations has reached {2 * self._max_eval} (two times the maximum function evaluations).')
        self._real_n_eval_cub += 1

        if self.n_eval_cub >= self._max_eval:
            # If the maximum number of evaluations has been reached, return
            # the last evaluated nonlinear inequality constraints.
            return self._last_cub

        # Generate the affine transformation.
        A, b = self._runtime.modifier_affine(self._seed, self._problem)[:2]

        # Evaluate the modified nonlinear inequality constraints and store the results.
        c = self._runtime.modifier_cub(A @ x + b, self._seed, self._problem, len(self._cub_hist))
        self._last_cub = c

        # Evaluate the nonlinear inequality constraints and store the results.
        c_true = super().cub(A @ x + b)

        # If the feature is 'quantized' and the option ``ground_truth'' is set to true, we should
        # set c_true to c.
        if self._runtime.name == 'quantized' and self._runtime.options[FeatureOption.GROUND_TRUTH]:
            c_true = c

        # Record the history of the nonlinear inequality constraints only when `record_hist` is true.
        if record_hist:
            self._cub_hist.append(c_true)

        return c

    def ceq(self, x, record_hist=True):
        """
        Evaluate the nonlinear constraints ``ceq(x) == 0``.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the nonlinear equality constraints.

        Returns
        -------
        `numpy.ndarray`, shape (m_nonlinear_eq,)
            Values of the nonlinear equality constraints at `x` if maximum number of evaluations has not been reached,
            otherwise returns the last evaluated nonlinear equality constraints.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape or if the return value of
            the argument `ceq` has an invalid shape.
        StopIteration
            If the number of the nonlinear equality constraint evaluations has reached two times the maximum function evaluations.
        """
        if self._real_n_eval_ceq >= 2 * self._max_eval:
            raise StopIteration(f'The number of the nonlinear equality constraint evaluations has reached {2 * self._max_eval} (two times the maximum function evaluations).')
        self._real_n_eval_ceq += 1

        if self.n_eval_ceq >= self._max_eval:
            # If the maximum number of evaluations has been reached, return
            # the last evaluated nonlinear equality constraints.
            return self._last_ceq

        # Generate the affine transformation.
        A, b = self._runtime.modifier_affine(self._seed, self._problem)[:2]

        # Evaluate the modified nonlinear equality constraints and store the results.
        c = self._runtime.modifier_ceq(A @ x + b, self._seed, self._problem, len(self._ceq_hist))
        self._last_ceq = c

        # Evaluate the nonlinear equality constraints and store the results.
        c_true = super().ceq(A @ x + b)

        # If the feature is 'quantized' and the option ``ground_truth'' is set to true, we should
        # set c_true to c.
        if self._runtime.name == 'quantized' and self._runtime.options[FeatureOption.GROUND_TRUTH]:
            c_true = c
        
        # Record the history of the nonlinear equality constraints only when `record_hist` is true.
        if record_hist:
            self._ceq_hist.append(c_true)

        return c

    def maxcv(self, x):
        """
        Evaluate the maximum constraint violation.

        Parameters
        ----------
        x : array_like, shape (n,)
            Point at which to evaluate the maximum constraint violation.

        Returns
        -------
        float
            Maximum constraint violation.

        Raises
        ------
        ValueError
            If the argument `x` has an invalid shape.
        """

        # If the Feature is ``quantized'' and the option ``ground_truth'' is set to true, we should
        # use the modified constraint violation.
        if self._runtime.name == 'quantized' and self._runtime.options[FeatureOption.GROUND_TRUTH]:
            if self.ptype == 'u':
                cv = 0.0
                return cv
            
            if np.any(np.isfinite(self.xl)):
                cv_bounds = np.max(self.xl - x, initial=0.0)
            else:
                cv_bounds = 0.0
            if np.any(np.isfinite(self.xu)):
                cv_bounds = np.max(x - self.xu, initial=cv_bounds)
            if self.ptype == 'b':
                cv = cv_bounds
                return cv
            
            if self.aub.size > 0:
                cv_linear = np.max(self.aub @ x - self.bub, initial=0.0)
            else:
                cv_linear = 0.0
            if self.aeq.size > 0:
                cv_linear = np.max(np.abs(self.aeq @ x - self.beq), initial=cv_linear)
            if self.ptype == 'l':
                cv = np.maximum(cv_bounds, cv_linear)
                return cv
            
            # Quantized modifiers are deterministic. Calling the public cub/ceq
            # oracle even with record_hist=False would still consume its real
            # evaluation budget and may reuse a cached last value after exhaustion.
            if self.m_nonlinear_ub > 0:
                cub = self._runtime.modifier_cub(x, self._seed, self._problem, self.n_eval_cub)
                cv_nonlinear = np.max(cub, initial=0.0)
            else:
                cv_nonlinear = 0.0
            if self.m_nonlinear_eq > 0:
                ceq = self._runtime.modifier_ceq(x, self._seed, self._problem, self.n_eval_ceq)
                cv_nonlinear = np.max(np.abs(ceq), initial=cv_nonlinear)
            cv = np.max([cv_bounds, cv_linear, cv_nonlinear])
            return cv
        else:
            # Generate the affine transformation.
            A, b = self._runtime.modifier_affine(self._seed, self._problem)[:2]
            return self._problem.maxcv(A @ x + b)
        
    # Note: We need to add methods `grad`, `hess`, `jcub`, and `jceq` to the FeaturedProblem class in the future.

def _validate_max_eval(max_eval):
    """Validate the evaluation budget of a featured problem and return it as an ``int``."""
    if isinstance(max_eval, (float, np.floating)) and float(max_eval).is_integer():
        max_eval = int(max_eval)
    if isinstance(max_eval, np.integer):
        max_eval = int(max_eval)
    if not isinstance(max_eval, int):
        raise TypeError('The argument `max_eval` for featured problem must be an integer.')
    if max_eval < 1:
        raise ValueError('The argument `max_eval` for featured problem must be positive.')
    return max_eval


def _validate_seed(seed):
    """Validate the seed of a featured problem (``None`` is allowed) and return it as an ``int``."""
    if seed is not None:
        if isinstance(seed, (float, np.floating)) and float(seed).is_integer():
            seed = int(seed)
        if isinstance(seed, np.integer):
            seed = int(seed)
        if not isinstance(seed, int):
            raise TypeError('The argument seed must be an integer.')
        if seed < 0:
            raise ValueError('The argument seed must be nonnegative.')
    return seed


def _process_1d_array(x, message):
    """
    Preprocess a one-dimensional array.

    Parameters
    ----------
    x : array_like
        Array to preprocess.
    message : str
        Error message to raise if the array is invalid.

    Returns
    -------
    `numpy.ndarray`
        Preprocessed array.

    Raises
    ------
    ValueError
        If the array is invalid.
    """
    if x is None:
        return np.empty(0)
    elif isinstance(x, np.ndarray) and x.size == 0:
        return x
    x = np.atleast_1d(np.squeeze(x)).astype(float)
    if x.ndim != 1:
        raise ValueError(message)
    return x

def _process_2d_array(x, message):
    """
    Preprocess a two-dimensional array.

    Parameters
    ----------
    x : array_like
        Array to preprocess.
    message : str
        Error message to raise if the array is invalid.

    Returns
    -------
    `numpy.ndarray`
        Preprocessed array.

    Raises
    ------
    ValueError
        If the array is invalid.
    """
    if x is None:
        return np.empty(0)
    elif isinstance(x, np.ndarray) and x.size == 0:
        return x
    x = np.atleast_2d(x).astype(float)
    if x.ndim != 2:
        raise ValueError(message)
    return x
