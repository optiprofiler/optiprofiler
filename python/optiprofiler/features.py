import numpy as np
from scipy.linalg import qr

from .utils import FeatureName, FeatureOption


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