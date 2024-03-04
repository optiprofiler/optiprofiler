import numpy as np

from .utils import FeatureName, FeatureOption, NoiseType


class Feature:
    """
    Feature used to modify the objective function.
    """

    def __init__(self, feature_name, **feature_options):
        """
        Initialize a feature.

        Parameters
        ----------
        feature_name : str
            Name of the feature.

        Other Parameters
        ----------------
        distribution : callable, optional
            Distribution used in the 'noisy' and 'randomize_x0' feature.
        modifier : callable, optional
            Modifier used in the 'custom' feature.
        n_runs : int, optional
            Number of runs for all features.
        order : int or float, optional
            Order of the 'regularized' feature.
        parameter : int or float, optional
            Regularization parameter of the 'regularized' feature.
        rate_nan : int or float, optional
            Rate of NaNs of the 'tough' feature.
        significant_digits : int, optional
            Number of significant digits of the 'truncated' feature.
        type : str, optional
            Type of the 'noisy' feature.

        Raises
        ------
        TypeError
            If an argument received an invalid value.
        ValueError
            If the arguments are inconsistent.
        """
        # Preprocess the feature name.
        self._name = feature_name
        if not isinstance(self._name, str):
            raise TypeError('The feature name must be a string.')
        self._name = self._name.lower()
        if self._name not in FeatureName.__members__.values():
            raise ValueError(f'Unknown feature: {self._name}.')

        # Preprocess the feature options.
        self._options = {k.lower(): v for k, v in feature_options.items()}
        for key in self._options:
            # Check whether the option is known.
            if key not in FeatureOption.__members__.values():
                raise ValueError(f'Unknown option: {key}.')

            # Check whether the options are valid for the feature.
            known_options = [FeatureOption.N_RUNS]
            if self._name == FeatureName.CUSTOM:
                known_options.extend([FeatureOption.MODIFIER])
            elif self._name == FeatureName.NOISY:
                known_options.extend([FeatureOption.DISTRIBUTION, FeatureOption.TYPE])
            elif self._name == FeatureName.RANDOMIZE_X0:
                known_options.extend([FeatureOption.DISTRIBUTION])
            elif self._name == FeatureName.REGULARIZE:
                known_options.extend([FeatureOption.ORDER, FeatureOption.PARAMETER])
            elif self._name == FeatureName.TOUGH:
                known_options.extend([FeatureOption.RATE_NAN])
            elif self._name == FeatureName.TRUNCATE:
                known_options.extend([FeatureOption.SIGNIFICANT_DIGITS])
            elif self._name == FeatureName.UNRELAXABLE_CONSTRAINTS:
                known_options.extend([FeatureOption.UNRELAXABLE_BOUNDS, FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS, FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS])
            elif self._name != FeatureName.PLAIN:
                raise NotImplementedError(f'Unknown feature: {self._name}.')
            if key not in known_options:
                raise ValueError(f'Option {key} is not valid for feature {self._name}.')

            # Check whether the options are valid.
            if key == FeatureOption.N_RUNS:
                if isinstance(self._options[key], float) and self._options[key].is_integer():
                    self._options[key] = int(self._options[key])
                if not isinstance(self._options[key], int):
                    raise TypeError(f'Option {key} must be an integer.')
                if self._options[key] <= 0:
                    raise ValueError(f'Option {key} must be integer.')
            elif key == FeatureOption.MODIFIER:
                if not callable(self._options[key]):
                    raise TypeError(f'Option {key} must be callable.')
            elif key == FeatureOption.DISTRIBUTION:
                if not callable(self._options[key]):
                    raise TypeError(f'Option {key} must be callable.')
            elif key == FeatureOption.ORDER:
                if not isinstance(self._options[key], (int, float)):
                    raise TypeError(f'Option {key} must be a number.')
            elif key == FeatureOption.PARAMETER:
                if not isinstance(self._options[key], (int, float)):
                    raise TypeError(f'Option {key} must be a number.')
                if self._options[key] < 0.0:
                    raise ValueError(f'Option {key} must be nonnegative.')
            elif key == FeatureOption.RATE_NAN:
                if not isinstance(self._options[key], (int, float)):
                    raise TypeError(f'Option {key} must be a number.')
                if not (0.0 <= self._options[key] <= 1.0):
                    raise ValueError(f'Option {key} must be between 0 and 1.')
            elif key == FeatureOption.SIGNIFICANT_DIGITS:
                if isinstance(self._options[key], float) and self._options[key].is_integer():
                    self._options[key] = int(self._options[key])
                if not isinstance(self._options[key], int):
                    raise TypeError(f'Option {key} must be an integer.')
                if self._options[key] <= 0:
                    raise ValueError(f'Option {key} must be positive.')
            elif key == FeatureOption.TYPE:
                if not isinstance(self._options[key], str):
                    raise TypeError(f'Option {key} must be a string.')
                if self._options[key].lower() not in NoiseType.__members__.values():
                    raise ValueError(f'Option {key} must be either "{NoiseType.ABSOLUTE.value}" or "{NoiseType.RELATIVE.value}".')
            elif key in [FeatureOption.UNRELAXABLE_BOUNDS, FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS, FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS]:
                if not isinstance(self._options[key], bool):
                    raise TypeError(f'Option {key} must be a boolean.')

        # Set default options.
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
        return self._options

    def modifier(self, x, f, maxcv_bounds=0.0, maxcv_linear=0.0, maxcv_nonlinear=0.0, seed=None):
        """
        Modify the objective function value according to the feature.

        Parameters
        ----------
        x : `numpy.ndarray`, shape (n,)
            Point at which the objective function is evaluated.
        f : float
            Objective function value at `x`.
        seed : int, optional
            Seed used to generate random numbers.

        Returns
        -------
        float
            Modified objective function value.
        """
        # Preprocess the point at which the objective function is evaluated.
        if not isinstance(x, np.ndarray):
            raise TypeError('The argument x must be a numpy array.')

        # Preprocess the objective function value.
        if not isinstance(f, float):
            raise TypeError('The argument f must be a float.')

        # Preprocess the seed.
        if seed is not None:
            if isinstance(seed, float) and seed.is_integer():
                seed = int(seed)
            if not isinstance(seed, int):
                raise TypeError('The argument seed must be an integer.')
            if seed < 0:
                raise ValueError('The argument seed must be nonnegative.')

        # Modify the objective function value.
        if self._name == FeatureName.CUSTOM:
            f = self._options[FeatureOption.MODIFIER](x, f, seed)
        elif self._name == FeatureName.NOISY:
            rng = self.get_default_rng(seed, f, sum(ord(letter) for letter in self._options[FeatureOption.TYPE]), *x)
            if self._options[FeatureOption.TYPE] == NoiseType.ABSOLUTE:
                f += self._options[FeatureOption.DISTRIBUTION](rng)
            else:
                f *= 1.0 + self._options[FeatureOption.DISTRIBUTION](rng)
        elif self._name == FeatureName.REGULARIZE:
            f += self._options[FeatureOption.PARAMETER] * np.linalg.norm(x, self._options[FeatureOption.ORDER])
        elif self._name == FeatureName.TOUGH:
            rng = self.get_default_rng(seed, f, self._options[FeatureOption.RATE_NAN], *x)
            if rng.uniform() < self._options[FeatureOption.RATE_NAN]:
                f = np.nan
        elif self._name == FeatureName.TRUNCATE:
            rng = self.get_default_rng(seed, f, self._options[FeatureOption.SIGNIFICANT_DIGITS], *x)
            if f == 0.0:
                digits = self._options[FeatureOption.SIGNIFICANT_DIGITS] - 1
            else:
                digits = self._options[FeatureOption.SIGNIFICANT_DIGITS] - int(np.floor(np.log10(np.abs(f)))) - 1
            if f >= 0.0:
                f = round(f, digits) + rng.uniform(0.0, 10.0 ** (-digits))
            else:
                f = round(f, digits) - rng.uniform(0.0, 10.0 ** (-digits))
        elif self._name == FeatureName.UNRELAXABLE_CONSTRAINTS:
            if self._options[FeatureOption.UNRELAXABLE_BOUNDS] and maxcv_bounds > 0.0:
                f = np.inf
            elif self._options[FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS] and maxcv_linear > 0.0:
                f = np.inf
            elif self._options[FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS] and maxcv_nonlinear > 0.0:
                f = np.inf
        elif self._name not in [FeatureName.PLAIN, FeatureName.RANDOMIZE_X0]:
            raise NotImplementedError(f'Unknown feature: {self._name}.')
        return f

    def _set_default_options(self):
        """
        Set default options.

        Notes
        -----
        The default distribution are defined as static methods of the class and
        not using lambda functions because the latter are not picklable.
        """

        if self._name == FeatureName.PLAIN:
            self._options.setdefault(FeatureOption.N_RUNS.value, 1)
        elif self._name == FeatureName.CUSTOM:
            if FeatureOption.MODIFIER not in self._options:
                raise ValueError(f'When using a custom feature, you must specify the {FeatureOption.MODIFIER} option.')
            self._options.setdefault(FeatureOption.N_RUNS.value, 1)
        elif self._name == FeatureName.NOISY:
            self._options.setdefault(FeatureOption.DISTRIBUTION.value, self._default_distribution)
            self._options.setdefault(FeatureOption.N_RUNS.value, 10)
            self._options.setdefault(FeatureOption.TYPE.value, NoiseType.RELATIVE.value)
        elif self._name == FeatureName.RANDOMIZE_X0:
            self._options.setdefault(FeatureOption.DISTRIBUTION.value, self._default_distribution)
            self._options.setdefault(FeatureOption.N_RUNS.value, 10)
        elif self._name == FeatureName.REGULARIZE:
            self._options.setdefault(FeatureOption.N_RUNS.value, 1)
            self._options.setdefault(FeatureOption.ORDER.value, 2)
            self._options.setdefault(FeatureOption.PARAMETER.value, 1.0)
        elif self._name == FeatureName.TOUGH:
            self._options.setdefault(FeatureOption.N_RUNS.value, 10)
            self._options.setdefault(FeatureOption.RATE_NAN.value, 0.05)
        elif self._name == FeatureName.TRUNCATE:
            self._options.setdefault(FeatureOption.N_RUNS.value, 10)
            self._options.setdefault(FeatureOption.SIGNIFICANT_DIGITS.value, 6)
        elif self._name == FeatureName.UNRELAXABLE_CONSTRAINTS:
            self._options.setdefault(FeatureOption.N_RUNS.value, 1)
            self._options.setdefault(FeatureOption.UNRELAXABLE_BOUNDS.value, True)
            self._options.setdefault(FeatureOption.UNRELAXABLE_LINEAR_CONSTRAINTS.value, False)
            self._options.setdefault(FeatureOption.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value, False)
        else:
            raise NotImplementedError(f'Unknown feature: {self._name}.')

    @staticmethod
    def _default_distribution(rng, n=None):
        return 1e-3 * rng.standard_normal(n)

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
            if isinstance(seed, float) and seed.is_integer():
                seed = int(seed)
            if not isinstance(seed, int):
                raise TypeError('The argument seed must be an integer.')
            if seed < 0:
                raise ValueError('The argument seed must be nonnegative.')

        # Preprocess the other arguments.
        for arg in args:
            if not isinstance(arg, (int, float)):
                raise TypeError('The arguments must be numbers.')

        # Generate the random number generator.
        rng = np.random.default_rng(seed)
        return np.random.default_rng(int(1e9 * abs(np.sin(1e5 * rng.standard_normal()) + sum(np.sin(1e5 * arg) for arg in args))))
