import numpy as np
from scipy.linalg import qr
import warnings
from numpy.linalg import lstsq
from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint, minimize

from .utils import FeatureName, FeatureOption, get_logger


class Feature:
    """
    Feature is a class that defines a mapping from an optimization problem to a new one with specified features.
    We are interested to test solvers on problems with different features. For example, we want to test the performance of solvers under the case where the objective function is noisy. For this purpose, we define `Feature` class.
    """

    def __init__(self, name, **feature_options):
        """
        Initialize a feature.

        Parameters
        ----------
        name : str
            Name of the feature.

        Other Parameters
        ----------------
        distribution : callable, optional
            Distribution used by the 'noisy' and 'randomize_x0' feature.
        modifier : callable, optional
            Custom modifier used by the 'custom' feature.
        n_runs : int, optional
            Number of runs for all features.
        rate_nan : int or float, optional
            Rate of NaNs used by the 'tough' feature.
        significant_digits : int, optional
            Number of significant digits used by the 'truncated' feature.
        type : str, optional
            Type of the noise used by the 'noisy' feature.

        Raises
        ------
        TypeError
            If an argument received an invalid value.
        ValueError
            If the arguments are inconsistent.
        """
        # Preprocess the feature name.
        self._name = name
        if not isinstance(self._name, str):
            raise TypeError('The first input argument for `Feature` must be a string.')
        self._name = self._name.lower()
        if self._name not in FeatureName.__members__.values():
            raise ValueError(f'Unknown feature: {self._name}.')

        # Preprocess the feature options.
        self._options = {k.lower(): v for k, v in feature_options.items()}
        for key in self._options:
            # Check whether the option is known.
            if key not in FeatureOption.__members__.values():
                raise ValueError(f'Unknown option for feature: {key}.')

            # Check whether the options are valid for the feature.
            known_options = [FeatureOption.N_RUNS]
            if self._name == FeatureName.CUSTOM:
                known_options.extend([FeatureOption.MOD_X0, FeatureOption.MOD_BOUNDS, FeatureOption.MOD_LINEAR_UB, FeatureOption.MOD_LINEAR_EQ, FeatureOption.MOD_AFFINE, FeatureOption.MOD_FUN, FeatureOption.MOD_CUB, FeatureOption.MOD_CEQ])
            elif self._name == FeatureName.NOISY:
                known_options.extend([FeatureOption.DISTRIBUTION, FeatureOption.NOISE_LEVEL, FeatureOption.NOISE_TYPE])
            elif self._name == FeatureName.PERTURBED_X0:
                known_options.extend([FeatureOption.DISTRIBUTION, FeatureOption.PERTURBATION_LEVEL])
            elif self._name == FeatureName.RANDOM_NAN:
                known_options.extend([FeatureOption.NAN_RATE])
            elif self._name == FeatureName.TRUNCATED:
                known_options.extend([FeatureOption.PERTURBED_TRAILING_DIGITS, FeatureOption.SIGNIFICANT_DIGITS])
            elif self._name == FeatureName.UNRELAXABLE_CONSTRAINTS:
                known_options.extend([FeatureOption.UNRELAXABLE_BOUNDS, FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS, FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS])
            elif self._name == FeatureName.LINEARLY_TRANSFORMED:
                known_options.extend([FeatureOption.ROTATED, FeatureOption.CONDITION_FACTOR])
            elif self._name == FeatureName.QUANTIZED:
                known_options.extend([FeatureOption.MESH_SIZE, FeatureOption.MESH_TYPE, FeatureOption.GROUND_TRUTH])
            elif self._name not in [FeatureName.PERMUTED, FeatureName.NONQUANTIFIABLE_CONSTRAINTS, FeatureName.PLAIN]:
                raise NotImplementedError(f'Unknown feature: {self._name}.')
            if key not in known_options:
                raise ValueError(f"Option `{key}` is not valid for feature '{self._name}'.")

            # Check whether the options are valid.
            if key == FeatureOption.N_RUNS:
                if isinstance(self._options[key], float) and self._options[key].is_integer():
                    self._options[key] = int(self._options[key])
                if not isinstance(self._options[key], int):
                    raise TypeError(f'Option `{key}` must be an integer.')
                if self._options[key] <= 0:
                    raise ValueError(f'Option `{key}` must be positive.')
            elif key == FeatureOption.DISTRIBUTION:
                if isinstance(self._options[key], str):
                    if self._name == FeatureName.NOISY and self._options[key] not in ['gaussian', 'uniform']:
                        raise ValueError(f'Option `{key}` for feature `{self._name}` must be either "gaussian" or "uniform" when specified as a string.')
                    elif self._name == FeatureName.PERTURBED_X0 and self._options[key] not in ['gaussian', 'spherical']:
                        raise ValueError(f'Option `{key}` for feature `{self._name}` must be either "gaussian" or "spherical" when specified as a string.')
                elif not callable(self._options[key]):
                    raise TypeError(f'Option `{key}` must be a string or it must be callable.')
            elif key == FeatureOption.NAN_RATE:
                if not isinstance(self._options[key], (int, float)):
                    raise TypeError(f'Option `{key}` must be a number.')
                if not (0.0 <= self._options[key] <= 1.0):
                    raise ValueError(f'Option `{key}` must be between 0 and 1.')
            elif key == FeatureOption.SIGNIFICANT_DIGITS:
                if isinstance(self._options[key], float) and self._options[key].is_integer():
                    self._options[key] = int(self._options[key])
                if not isinstance(self._options[key], int):
                    raise TypeError(f'Option `{key}` must be an integer.')
                if self._options[key] <= 0:
                    raise ValueError(f'Option `{key}` must be positive.')
            elif key in [FeatureOption.NOISE_LEVEL, FeatureOption.CONDITION_FACTOR]:
                if not isinstance(self._options[key], (int, float)):
                    raise TypeError(f'Option `{key}` must be a number.')
                if self._options[key] < 0.0:
                    raise ValueError(f'Option `{key}` must be nonnegative.')
            elif key == FeatureOption.NOISE_TYPE:
                if not isinstance(self._options[key], str):
                    raise TypeError(f'Option {key} must be a string.')
                if self._options[key].lower() not in ['absolute', 'relative', 'mixed']:
                    raise ValueError(f"Option `{key}` must be one of 'absolute', 'relative', or 'mixed'.")
            elif key in [FeatureOption.PERTURBED_TRAILING_DIGITS, FeatureOption.ROTATED, FeatureOption.UNRELAXABLE_BOUNDS, FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS, FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS, FeatureOption.GROUND_TRUTH]:
                if not isinstance(self._options[key], bool):
                    raise TypeError(f'Option `{key}` must be a boolean.')
            elif key == FeatureOption.MESH_SIZE:
                if not isinstance(self._options[key], (int, float)):
                    raise TypeError(f'Option `{key}` must be a number.')
                if self._options[key] <= 0.0:
                    raise ValueError(f'Option `{key}` must be positive.')
            elif key == FeatureOption.MESH_TYPE:
                if not isinstance(self._options[key], str):
                    raise TypeError(f'Option `{key}` must be a string.')
                if self._options[key].lower() not in ['absolute', 'relative']:
                    raise ValueError(f"Option `{key}` must be 'absolute' or 'relative'.")
            elif key in [FeatureOption.MOD_X0, FeatureOption.MOD_BOUNDS, FeatureOption.MOD_LINEAR_UB, FeatureOption.MOD_LINEAR_EQ, FeatureOption.MOD_AFFINE, FeatureOption.MOD_FUN, FeatureOption.MOD_CUB, FeatureOption.MOD_CEQ]:
                if not callable(self._options[key]):
                    raise TypeError(f'Option `{key}` must be callable.')

        # Set default options for the unspecified options.
        self._set_default_options()

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
        """
        Whether the feature is stochastic.

        Returns
        -------
        bool
            Whether the feature is stochastic.
        """
        if self._name in [FeatureName.PERTURBED_X0, FeatureName.NOISY, FeatureName.PERMUTED, FeatureName.RANDOM_NAN, FeatureName.CUSTOM]:
            return True
        elif self._name == FeatureName.TRUNCATED:
            return self._options[FeatureOption.PERTURBED_TRAILING_DIGITS]
        elif self._name == FeatureName.LINEARLY_TRANSFORMED:
            return self._options[FeatureOption.ROTATED]
        else:
            return False

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
            perturbation_level = self._options[FeatureOption.PERTURBATION_LEVEL] * max(1, np.linalg.norm(problem.x0))
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
            # Generate a positive definite diagonal matrix D with condition number equal to
            # 2^(condition_factor * n / 2), where n is the dimension of the problem. In this
            # way, the condition number of Q * D^2 * Q^T is 2^(condition_factor * n).
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
            if np.count_nonzero(A - np.diag(np.diagonal(A))) != 0:
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
            rng_noisy = self.get_default_rng(seed, f, *x, n_eval)
            if self._options[FeatureOption.DISTRIBUTION] == 'gaussian':
                noise = rng_noisy.standard_normal()
            elif self._options[FeatureOption.DISTRIBUTION] == 'uniform':
                noise = rng_noisy.uniform(-1, 1)
            else:
                noise = self._options[FeatureOption.DISTRIBUTION](rng_noisy)
            if self._options[FeatureOption.NOISE_TYPE] == 'absolute':
                return f + self._options[FeatureOption.NOISE_LEVEL] * noise
            elif self._options[FeatureOption.NOISE_TYPE] == 'relative':
                return f * (1.0 + self._options[FeatureOption.NOISE_LEVEL] * noise)
            else:
                return f + max(1, np.abs(f)) * self._options[FeatureOption.NOISE_LEVEL] * noise
        elif self._name == FeatureName.RANDOM_NAN:
            rng_random_nan = self.get_default_rng(seed, f, *x, n_eval)
            if rng_random_nan.random() < self._options[FeatureOption.NAN_RATE]:
                return np.nan
            else:
                return f
        elif self._name == FeatureName.TRUNCATED:
            if np.isnan(f) or np.isinf(f):
                # If f is NaN or Inf, we do not need to truncate it.
                # Note that if f is NaN or Inf, digits will be set to NaN or Inf respectively, which will lead
                # to an error when calling 'round(f, digits)'.
                return f
            rng_truncated = self.get_default_rng(seed, f, *x, n_eval)
            if f == 0.0:
                digits = self._options[FeatureOption.SIGNIFICANT_DIGITS] - 1
            else:
                digits = self._options[FeatureOption.SIGNIFICANT_DIGITS] - int(np.log10(np.abs(f))) - 1
            f = round(f, digits)
            # Round f to the desired number of significant digits.
            if self._options[FeatureOption.PERTURBED_TRAILING_DIGITS]:
                if f >= 0.0:
                    f += rng_truncated.uniform(0.0, 10.0 ** (-digits))
                else:
                    f -= rng_truncated.uniform(0.0, 10.0 ** (-digits))
            return f
        elif self._name == FeatureName.UNRELAXABLE_CONSTRAINTS:
            _, maxcv_bounds, maxcv_linear, maxcv_nonlinear = problem.maxcv(x)
            if self._options[FeatureOption.UNRELAXABLE_BOUNDS] and maxcv_bounds > 0.0:
                return np.inf
            elif self._options[FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS] and maxcv_linear > 0.0:
                return np.inf
            elif self._options[FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS] and maxcv_nonlinear > 0.0:
                return np.inf
            else:
                return f
        elif self._name == FeatureName.QUANTIZED:
            mesh_size = self._options[FeatureOption.MESH_SIZE]
            if self._options[FeatureOption.MESH_TYPE] == 'relative':
                mesh_size *= np.maximum(1, np.abs(x))
            x = mesh_size * np.round(x / mesh_size)
            return problem.fun(x)
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
            rng_noisy = self.get_default_rng(seed, *cub, *x, n_eval_cub)
            if self._options[FeatureOption.DISTRIBUTION] == 'gaussian':
                noise = rng_noisy.standard_normal(cub.size)
            elif self._options[FeatureOption.DISTRIBUTION] == 'uniform':
                noise = rng_noisy.uniform(-1, 1, cub.size)
            else:
                noise = self._options[FeatureOption.DISTRIBUTION](rng_noisy, cub.size)
            if self._options[FeatureOption.NOISE_TYPE] == 'absolute':
                return cub + self._options[FeatureOption.NOISE_LEVEL] * noise
            elif self._options[FeatureOption.NOISE_TYPE] == 'relative':
                return cub * (1.0 + self._options[FeatureOption.NOISE_LEVEL] * noise)
            else:
                return cub + np.maximum(1, np.abs(cub)) * self._options[FeatureOption.NOISE_LEVEL] * noise
        elif self._name == FeatureName.RANDOM_NAN:
            # Similar to the case in the modifier_fun method.
            rng_random_nan = self.get_default_rng(seed, *cub, *x, n_eval_cub)
            cub[rng_random_nan.random(cub.size) < self._options[FeatureOption.NAN_RATE]] = np.nan
            return cub
        elif self._name == FeatureName.TRUNCATED:
            # Similar to the case in the modifier_fun method.
            rng_truncated = self.get_default_rng(seed, *cub, *x, n_eval_cub)
            digits = np.zeros(cub.size)
            digits[cub == 0.0] = self._options[FeatureOption.SIGNIFICANT_DIGITS] - 1
            digits[cub != 0.0] = self._options[FeatureOption.SIGNIFICANT_DIGITS] - np.int_(np.log10(np.abs(cub[cub != 0.0]))) - 1
            for i in range(cub.size):
                if not np.isnan(cub[i]) and not np.isinf(cub[i]):
                    cub[i] = round(cub[i], digits[i])
            if self._options[FeatureOption.PERTURBED_TRAILING_DIGITS]:
                cub[cub >= 0.0] += rng_truncated.uniform(0.0, 10.0 ** (-digits[cub >= 0.0]))
                cub[cub < 0.0] -= rng_truncated.uniform(0.0, 10.0 ** (-digits[cub < 0.0]))
            return cub
        elif self._name == FeatureName.NONQUANTIFIABLE_CONSTRAINTS:
            # Set the elements whose value are less than or equal to 0 to 0.
            cub[cub <= 0.0] = 0.0
            # Set the rest to 1.
            cub[cub > 0.0] = 1.0
            return cub
        elif self._name == FeatureName.QUANTIZED:
            # Similar to the case in the modifier_fun method.
            mesh_size = self._options[FeatureOption.MESH_SIZE]
            if self._options[FeatureOption.MESH_TYPE] == 'relative':
                mesh_size *= np.maximum(1, np.abs(x))
            x = mesh_size * np.round(x / mesh_size)
            return problem.cub(x)
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
            rng_noisy = self.get_default_rng(seed, *ceq, *x, n_eval_ceq)
            if self._options[FeatureOption.DISTRIBUTION] == 'gaussian':
                noise = rng_noisy.standard_normal(ceq.size)
            elif self._options[FeatureOption.DISTRIBUTION] == 'uniform':
                noise = rng_noisy.uniform(-1, 1, ceq.size)
            else:
                noise = self._options[FeatureOption.DISTRIBUTION](rng_noisy, ceq.size)
            if self._options[FeatureOption.NOISE_TYPE] == 'absolute':
                return ceq + self._options[FeatureOption.NOISE_LEVEL] * noise
            elif self._options[FeatureOption.NOISE_TYPE] == 'relative':
                return ceq * (1.0 + self._options[FeatureOption.NOISE_LEVEL] * noise)
            else:
                return ceq + np.maximum(1, np.abs(ceq)) * self._options[FeatureOption.NOISE_LEVEL] * noise
        elif self._name == FeatureName.RANDOM_NAN:
            # Similar to the case in the modifier_fun method.
            rng_random_nan = self.get_default_rng(seed, *ceq, *x, n_eval_ceq)
            ceq[rng_random_nan.random(ceq.size) < self._options[FeatureOption.NAN_RATE]] = np.nan
            return ceq
        elif self._name == FeatureName.TRUNCATED:
            # Similar to the case in the modifier_fun method.
            rng_truncated = self.get_default_rng(seed, *ceq, *x, n_eval_ceq)
            digits = np.zeros(ceq.size)
            digits[ceq == 0.0] = self._options[FeatureOption.SIGNIFICANT_DIGITS] - 1
            digits[ceq != 0.0] = self._options[FeatureOption.SIGNIFICANT_DIGITS] - np.int_(np.log10(np.abs(ceq[ceq != 0.0]))) - 1
            for i in range(ceq.size):
                if not np.isnan(ceq[i]) and not np.isinf(ceq[i]):
                    ceq[i] = round(ceq[i], digits[i])
            if self._options[FeatureOption.PERTURBED_TRAILING_DIGITS]:
                ceq[ceq >= 0.0] += rng_truncated.uniform(0.0, 10.0 ** (-digits[ceq >= 0.0]))
                ceq[ceq < 0.0] -= rng_truncated.uniform(0.0, 10.0 ** (-digits[ceq < 0.0]))
            return ceq
        elif self._name == FeatureName.NONQUANTIFIABLE_CONSTRAINTS:
            # Set the elements whose absolute value are less than or equal to 10^(-6) to 0.
            ceq[np.abs(ceq) <= 1e-6] = 0.0
            # Set the rest to 1.
            ceq[np.abs(ceq) > 1e-6] = 1.0
            return ceq
        elif self._name == FeatureName.QUANTIZED:
            # Similar to the case in the modifier_fun method.
            mesh_size = self._options[FeatureOption.MESH_SIZE]
            if self._options[FeatureOption.MESH_TYPE] == 'relative':
                mesh_size *= np.maximum(1, np.abs(x))
            x = mesh_size * np.round(x / mesh_size)
            return problem.ceq(x)
        else:
            return ceq

    def _set_default_options(self):
        """
        Set default options.

        Notes
        -----
        The default distribution are defined as static methods of the class and
        not using lambda functions because the latter are not picklable.
        """

        if self._name in [FeatureName.PLAIN, FeatureName.CUSTOM, FeatureName.NONQUANTIFIABLE_CONSTRAINTS]:
            self._options.setdefault(FeatureOption.N_RUNS.value, 1)
        elif self._name == FeatureName.NOISY:
            self._options.setdefault(FeatureOption.DISTRIBUTION.value, 'gaussian')
            self._options.setdefault(FeatureOption.N_RUNS.value, 5)
            self._options.setdefault(FeatureOption.NOISE_LEVEL.value, 1e-3)
            self._options.setdefault(FeatureOption.NOISE_TYPE.value, 'mixed')
        elif self._name == FeatureName.PERMUTED:
            self._options.setdefault(FeatureOption.N_RUNS.value, 5)
        elif self._name == FeatureName.LINEARLY_TRANSFORMED:
            self._options.setdefault(FeatureOption.ROTATED.value, True)
            self._options.setdefault(FeatureOption.CONDITION_FACTOR.value, 0)
            self._options.setdefault(FeatureOption.N_RUNS.value, 5)
        elif self._name == FeatureName.PERTURBED_X0:
            self._options.setdefault(FeatureOption.DISTRIBUTION.value, 'spherical')
            self._options.setdefault(FeatureOption.PERTURBATION_LEVEL.value, 1e-3)
            self._options.setdefault(FeatureOption.N_RUNS.value, 5)
        elif self._name == FeatureName.RANDOM_NAN:
            self._options.setdefault(FeatureOption.N_RUNS.value, 5)
            self._options.setdefault(FeatureOption.NAN_RATE.value, 0.05)
        elif self._name == FeatureName.TRUNCATED:
            self._options.setdefault(FeatureOption.PERTURBED_TRAILING_DIGITS.value, False)
            if self._options[FeatureOption.PERTURBED_TRAILING_DIGITS]:
                self._options.setdefault(FeatureOption.N_RUNS.value, 5)
            else:
                self._options.setdefault(FeatureOption.N_RUNS.value, 1)
            self._options.setdefault(FeatureOption.SIGNIFICANT_DIGITS.value, 6)
        elif self._name == FeatureName.UNRELAXABLE_CONSTRAINTS:
            self._options.setdefault(FeatureOption.N_RUNS.value, 1)
            self._options.setdefault(FeatureOption.UNRELAXABLE_BOUNDS.value, True)
            self._options.setdefault(FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS.value, False)
            self._options.setdefault(FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value, False)
        elif self._name == FeatureName.QUANTIZED:
            self._options.setdefault(FeatureOption.N_RUNS.value, 1)
            self._options.setdefault(FeatureOption.MESH_SIZE.value, 1e-3)
            self._options.setdefault(FeatureOption.MESH_TYPE.value, 'absolute')
            self._options.setdefault(FeatureOption.GROUND_TRUTH.value, True)
        else:
            raise NotImplementedError(f'Unknown feature: {self._name}.')

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
            if isinstance(seed, float) and seed.is_integer():
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

class Problem:
    r"""
    Optimization problem to be used in the benchmarking.

    This class provides a uniform interface for providing general optimization
    problems. It is used to supply custom problems to the `benchmark`
    function, and the signature of solvers supplied to the `benchmark`
    function can be

        ``solver(problem) -> array_like, shape (n,)``

    where `problem` is an instance of the class `Problem`.

    Examples
    --------
    Consider the unconstrained problem of minimizing the Rosenbrock function

    .. math::

        f(x) = 100 (x_2 - x_1^2)^2 + (1 - x_1)^2.

    To create an instance of the class `Problem` for this problem, run:

    >>> from optiprofiler.modules import Problem
    >>>
    >>> def rosen(x):
    ...     return 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2
    ...
    >>> problem = Problem(rosen, [0, 0])

    The second argument ``[0, 0]`` is the initial guess. This instance can now
    be used to evaluate the objective function at any point and access extra
    information about the problem. For example, to evaluate the objective
    function at the initial guess, run:

    >>> problem.fun(problem.x0)
    1.0

    All the attributes and methods of the class `Problem` described below are
    always well-defined. For example, since the problem ``problem`` is
    unconstrained, the components of the lower and upper bounds on the
    variables must be :math:`-\infty` and :math:`\infty`, respectively:

    >>> problem.xl
    array([-inf, -inf])
    >>> problem.xu
    array([inf, inf])

    The optional arguments of the constructor of the class `Problem` can be
    used to specify constraints. For example, to specify that the variables
    must be nonnegative, run:

    >>> problem = Problem(rosen, [0, 0], [0, 0])

    The lower bounds on the variables are now zero:

    >>> problem.xl
    array([0., 0.])
    >>> problem.xu
    array([inf, inf])

    The optional arguments ``aub``, ``bub``, ``aeq``, and ``beq`` can be
    used to specify linear inequality and equality constraints, respectively.
    The optional arguments ``cub`` and ``ceq`` can be used to specify
    nonlinear inequality and equality constraints, respectively. For example,
    to specify that the variables must satisfy the constraint
    :math:`x_1^2 + x_2^2 \le 1` and :math:`x_1^3 - x_2^2 \le 1`, run:

    >>> def cub(x):
    ...     return [x[0] ** 2 + x[1] ** 2 - 1, x[0] ** 3 - x[1] ** 2 - 1]
    ...
    >>> problem = Problem(rosen, [0, 0], cub=cub, m_nonlinear_ub=2)

    The instance ``problem`` can now be used to evaluate the nonlinear
    inequality constraints at any point. For example, to evaluate the nonlinear
    inequality constraints at the initial guess, run:

    >>> problem.cub(problem.x0)
    array([-1., -1.])

    If you do not provide the number of nonlinear inequality constraints in
    ``m_nonlinear_ub``, it will be inferred at the first call to ``cub``.
    Nonlinear equality constraints can be specified in a similar way, using the
    optional arguments ``ceq`` and ``m_nonlinear_eq``.
    """

    def __init__(self, fun, x0, name=None, xl=None, xu=None, aub=None, bub=None, aeq=None, beq=None, cub=None, ceq=None, grad=None, hess=None, jcub=None, jceq=None, hcub=None, hceq=None):
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

        Raises
        ------
        TypeError
            If an argument received an invalid value.
        ValueError
            If the arguments are inconsistent.
        """

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
                    if self._cub(self._x0) is not None:
                        self._m_nonlinear_ub = self._cub(self._x0).size
                    else:
                        self._m_nonlinear_ub = 0
                except Exception as err:
                    logger = get_logger(__name__)
                    logger.warning(f'Failed to evaluate the nonlinear inequality constraint function: {err}')
        self._m_nonlinear_eq = 0
        self._ceq = ceq
        if self._ceq is not None:
            if not callable(self._ceq):
                raise TypeError('The argument `ceq` for problem must be callable.')
            else:
                try:
                    if self._ceq(self._x0) is not None:
                        self._m_nonlinear_eq = self._ceq(self._x0).size
                    else:
                        self._m_nonlinear_eq = 0
                except Exception as err:
                    logger = get_logger(__name__)
                    logger.warning(f'Failed to evaluate nonlinear equality constraints: {err}')

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
            logger.warning(f'Failed to evaluate the objective function: {err}')
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
                logger.warning(f'Failed to evaluate the gradient of the objective function: {err}')
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
                logger.warning(f'Failed to evaluate the Hessian of the objective function: {err}')
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
                logger.warning(f'Failed to evaluate the nonlinear inequality constraint function: {err}')
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
                logger.warning(f'Failed to evaluate the nonlinear equality constraint function: {err}')
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
                logger.warning(f'Failed to evaluate the Jacobian of the nonlinear inequality constraint function: {err}')
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
                logger.warning(f'Failed to evaluate the Jacobian of the nonlinear equality constraint function: {err}')
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
                logger.warning(f'Failed to evaluate the Hessian of the nonlinear inequality constraint function: {err}')
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
                logger.warning(f'Failed to evaluate the Hessian of the nonlinear equality constraint function: {err}')
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
            cv = max(cv_bounds)
            return cv, cv_bounds, cv_linear, cv_nonlinear
        
        if self.aub.size > 0:
            cv_linear = np.max(self.aub @ x - self.bub, initial=0.0)
        if self.aeq.size > 0:
            cv_linear = np.max(np.abs(self.aeq @ x - self.beq), initial=cv_linear)
        if self.ptype == 'l':
            cv = max(cv_bounds, cv_linear)
            return cv, cv_bounds, cv_linear, cv_nonlinear
        
        if self.m_nonlinear_ub > 0:
            cv_nonlinear = np.max(self.cub(x), initial=0.0)
        if self.m_nonlinear_eq > 0:
            cv_nonlinear = np.max(np.abs(self.ceq(x)), initial=cv_nonlinear)
        cv = max(cv_bounds, cv_linear, cv_nonlinear)
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
            self._x0 = lstsq(self.aeq, self.beq - self.aeq @ self.x0)[0]
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
            def hessp(p):
                return p

            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = minimize(dist_x0_sq, self.x0, jac=True, hessp=hessp, bounds=bounds, constraints=constraints, tol=1e-16)
            self._x0 = res.x

class FeaturedProblem(Problem):
    """
    Optimization problem to be used in the benchmarking, with extra features.
    """

    def __init__(self, problem, feature, max_eval, seed=None):
        """
        Initialize an optimization problem.

        Parameters
        ----------
        problem : `optiprofiler.modules.Problem`
            Problem to be used in the benchmarking.
        feature : `optiprofiler.features.Feature`
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

        # Preprocess the maximum number of function evaluations.
        self._max_eval = max_eval
        if isinstance(self._max_eval, float) and self._max_eval.is_integer():
            self._max_eval = int(self._max_eval)
        if not isinstance(self._max_eval, int):
            raise TypeError('The argument `max_eval` for featured problem must be an integer.')
        if self._max_eval < 1:
            raise ValueError('The argument `max_eval` for featured problem must be positive.')

        # Preprocess the seed.
        self._seed = seed
        if self._seed is not None:
            if isinstance(self._seed, float) and self._seed.is_integer():
                self._seed = int(self._seed)
            if not isinstance(self._seed, int):
                raise TypeError('The argument seed must be an integer.')
            if self._seed < 0:
                raise ValueError('The argument seed must be nonnegative.')
            
        # Record the real evaluation numbers.
        self._real_n_eval_fun = 0
        self._real_n_eval_cub = 0
        self._real_n_eval_ceq = 0

        # Modify the problem according to the feature.
        self._x0 = self._feature.modifier_x0(self._seed, self._problem)
        self._xl, self._xu = self._feature.modifier_bounds(self._seed, self._problem)
        self._aub, self._bub = self._feature.modifier_linear_ub(self._seed, self._problem)
        self._aeq, self._beq = self._feature.modifier_linear_eq(self._seed, self._problem)

        # Store the histories of the objective function values, nonlinear
        # constraints, and maximum constraint violations.
        self._fun_hist = []
        self._cub_hist = []
        self._ceq_hist = []
        self._maxcv_hist = []
        self._last_fun = np.nan
        self._last_cub = np.nan
        self._last_ceq = np.nan

        # Evaluate the objective function and the maximum constraint violation at the initial point.
        # Pay attention to the case when the feature is 'quantized' and the option ``ground_truth'' is set to true.
        A, b = self._feature.modifier_affine(self._seed, self._problem)[:2]
        if self._feature.name == 'quantized' and self._feature.options[FeatureOption.GROUND_TRUTH]:
            self._fun_init = self._feature.modifier_fun(A @ self._x0 + b, self._seed, self._problem, self.n_eval_fun)
        else:
            self._fun_init = self._problem.fun(A @ self._x0 + b)
        self._maxcv_init = self._problem.maxcv(A @ self._x0 + b)


    def __new__(cls, problem, feature, max_eval, seed=None):
        # Preprocess the problem.
        if not isinstance(problem, Problem):
            raise TypeError('The argument `problem` for featured problem must be an instance of the class Problem.')

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
        A, b = self._feature.modifier_affine(self._seed, self._problem)[:2]

        # Evaluate the modified the objective function value according to the feature and return the
        # modified value.
        f = self._feature.modifier_fun(A @ x + b, self._seed, self._problem, self.n_eval_fun)
        self._last_fun = f

        # Evaluate the objective function and store the results.
        f_true = super().fun(A @ x + b)

        # If the feature is 'quantized' and the option ``ground_truth'' is set to true, we should
        # set f_true to f.
        if self._feature.name == 'quantized' and self._feature.options[FeatureOption.GROUND_TRUTH]:
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
        A, b = self._feature.modifier_affine(self._seed, self._problem)[:2]

        # Evaluate the modified nonlinear inequality constraints and store the results.
        c = self._feature.modifier_cub(A @ x + b, self._seed, self._problem, len(self._cub_hist))
        self._last_cub = c

        # Evaluate the nonlinear inequality constraints and store the results.
        c_true = super().cub(A @ x + b)

        # If the feature is 'quantized' and the option ``ground_truth'' is set to true, we should
        # set c_true to c.
        if self._feature.name == 'quantized' and self._feature.options[FeatureOption.GROUND_TRUTH]:
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
        A, b = self._feature.modifier_affine(self._seed, self._problem)[:2]

        # Evaluate the modified nonlinear equality constraints and store the results.
        c = self._feature.modifier_ceq(A @ x + b, self._seed, self._problem, len(self._ceq_hist))
        self._last_ceq = c

        # Evaluate the nonlinear equality constraints and store the results.
        c_true = super().ceq(A @ x + b)

        # If the feature is 'quantized' and the option ``ground_truth'' is set to true, we should
        # set c_true to c.
        if self._feature.name == 'quantized' and self._feature.options[FeatureOption.GROUND_TRUTH]:
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
        if self._feature.name == 'quantized' and self._feature.options[FeatureOption.GROUND_TRUTH]:
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
                cv = max(cv_bounds)
                return cv
            
            if self.aub.size > 0:
                cv_linear = np.max(self.aub @ x - self.bub, initial=0.0)
            else:
                cv_linear = 0.0
            if self.aeq.size > 0:
                cv_linear = np.max(np.abs(self.aeq @ x - self.beq), initial=cv_linear)
            if self.ptype == 'l':
                cv = max(cv_bounds, cv_linear)
                return cv
            
            if self.m_nonlinear_ub > 0:
                cv_nonlinear = np.max(self.cub(x, record_hist=False), initial=0.0)
            else:
                cv_nonlinear = 0.0
            if self.m_nonlinear_eq > 0:
                cv_nonlinear = np.max(np.abs(self.ceq(x, record_hist=False)), initial=cv_nonlinear)
            cv = max(cv_bounds, cv_linear, cv_nonlinear)
            return cv
        else:
            # Generate the affine transformation.
            A, b = self._feature.modifier_affine(self._seed, self._problem)[:2]
            return self._problem.maxcv(A @ x + b)
        
    # Note: We need to add methods `grad`, `hess`, `jcub`, and `jceq` to the FeaturedProblem class in the future.

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