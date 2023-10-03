from enum import Enum

import numpy as np


class FeatureName(str, Enum):
    """
    Enumeration of the available features.
    """
    CUSTOM = 'custom'
    NOISY = 'noisy'
    PLAIN = 'plain'
    REGULARIZED = 'regularized'
    TOUGH = 'tough'
    TRUNCATED = 'truncated'


class OptionKey(str, Enum):
    """
    Enumeration of the available options.
    """
    DISTRIBUTION = 'distribution'
    MODIFIER = 'modifier'
    N_RUNS = 'n_runs'
    ORDER = 'order'
    PARAMETER = 'parameter'
    RATE_ERROR = 'rate_error'
    RATE_NAN = 'rate_nan'
    SIGNIFICANT_DIGITS = 'significant_digits'
    TYPE = 'type'


class NoiseType(str, Enum):
    """
    Enumeration of the available noise types.
    """
    ABSOLUTE = 'absolute'
    RELATIVE = 'relative'


class Feature:
    """
    Feature used to modify the objective function.
    """

    def __init__(self, feature_name, **kwargs):
        """
        Initialize a feature.

        Parameters
        ----------
        feature_name : str
            Name of the feature.

        Other Parameters
        ----------------
        distribution : callable, optional
            Distribution used in the noisy feature.
        modifier : callable, optional
            Modifier used in the custom feature.
        n_runs : int, optional
            Number of runs for all features.
        order : int or float, optional
            Order of the regularized feature.
        parameter : int or float, optional
            Regularization parameter of the regularized feature.
        rate_error : int or float, optional
            Rate of errors of the tough feature.
        rate_nan : int or float, optional
            Rate of NaNs of the tough feature.
        significant_digits : int, optional
            Number of significant digits of the truncated feature.
        type : str, optional
            Type of the noisy feature.
        """
        self._name = feature_name.lower()
        self._options = {k.lower(): v for k, v in kwargs.items()}

        # Check whether the feature is valid.
        if self._name not in FeatureName.__members__.values():
            raise ValueError(f'Unknown feature: {self._name}.')
        for key in self._options:
            # Check whether the option is known.
            if key not in OptionKey.__members__.values():
                raise ValueError(f'Unknown option: {key}.')

            # Check whether the option is valid for the feature.
            known_options = [OptionKey.N_RUNS]
            if self._name == FeatureName.CUSTOM:
                known_options.extend([OptionKey.MODIFIER])
            elif self._name == FeatureName.NOISY:
                known_options.extend([OptionKey.DISTRIBUTION, OptionKey.TYPE])
            elif self._name == FeatureName.REGULARIZED:
                known_options.extend([OptionKey.ORDER, OptionKey.PARAMETER])
            elif self._name == FeatureName.TOUGH:
                known_options.extend([OptionKey.RATE_ERROR, OptionKey.RATE_NAN])
            elif self._name == FeatureName.TRUNCATED:
                known_options.extend([OptionKey.SIGNIFICANT_DIGITS])
            elif self._name != FeatureName.PLAIN:
                raise NotImplementedError(f'Unknown feature: {self._name}.')
            if key not in known_options:
                raise ValueError(f'Option {key} is not valid for feature {self._name}.')

            # Check whether the option type is valid.
            if key == OptionKey.N_RUNS:
                if isinstance(self._options[key], float) and self._options[key].is_integer():
                    self._options[key] = int(self._options[key])
                if not isinstance(self._options[key], int) or self._options[key] <= 0:
                    raise TypeError(f'Option {key} must be a positive integer.')
            elif key == OptionKey.MODIFIER:
                if not callable(self._options[key]):
                    raise TypeError(f'Option {key} must be callable.')
            elif key == OptionKey.DISTRIBUTION:
                if not callable(self._options[key]):
                    raise TypeError(f'Option {key} must be callable.')
            elif key == OptionKey.ORDER:
                if not isinstance(self._options[key], (int, float)):
                    raise TypeError(f'Option {key} must be a number.')
            elif key == OptionKey.PARAMETER:
                if not isinstance(self._options[key], (int, float)) or self._options[key] < 0.0:
                    raise TypeError(f'Option {key} must be a nonnegative number.')
            elif key in [OptionKey.RATE_ERROR, OptionKey.RATE_NAN]:
                if not isinstance(self._options[key], (int, float)) or not (0.0 <= self._options[key] <= 1.0):
                    raise TypeError(f'Option {key} must be a number between 0 and 1.')
            elif key == OptionKey.SIGNIFICANT_DIGITS:
                if isinstance(self._options[key], float) and self._options[key].is_integer():
                    self._options[key] = int(self._options[key])
                if not isinstance(self._options[key], int) or self._options[key] <= 0:
                    raise TypeError(f'Option {key} must be a positive integer.')
            elif key == OptionKey.TYPE:
                if not isinstance(self._options[key], str) or self._options[key].lower() not in NoiseType.__members__.values():
                    raise TypeError(f'Option {key} must be either "{NoiseType.ABSOLUTE.value}" or "{NoiseType.RELATIVE.value}".')

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
        return dict(self._options)

    def modifier(self, x, f, seed=None):
        """
        Modify the objective function value.

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
        if self._name == FeatureName.CUSTOM:
            f = self._options[OptionKey.MODIFIER](x, f, seed)
        elif self._name == FeatureName.NOISY:
            rng = self._default_rng(seed, f, sum(ord(letter) for letter in self._options[OptionKey.TYPE]), *x)
            if self._options[OptionKey.TYPE] == NoiseType.ABSOLUTE:
                f += self._options[OptionKey.DISTRIBUTION](rng)
            else:
                f *= 1.0 + self._options[OptionKey.DISTRIBUTION](rng)
        elif self._name == FeatureName.REGULARIZED:
            f += self._options[OptionKey.PARAMETER] * np.linalg.norm(x, self._options[OptionKey.ORDER])
        elif self._name == FeatureName.TOUGH:
            rng = self._default_rng(seed, f, self._options[OptionKey.RATE_ERROR], self._options[OptionKey.RATE_NAN], *x)
            if rng.uniform() < self._options[OptionKey.RATE_ERROR]:
                raise RuntimeError
            elif rng.uniform() < self._options[OptionKey.RATE_NAN]:
                f = np.nan
        elif self._name == FeatureName.TRUNCATED:
            rng = self._default_rng(seed, f, self._options[OptionKey.SIGNIFICANT_DIGITS], *x)
            if f == 0.0:
                digits = self._options[OptionKey.SIGNIFICANT_DIGITS] - 1
            else:
                digits = self._options[OptionKey.SIGNIFICANT_DIGITS] - int(np.floor(np.log10(np.abs(f)))) - 1
            if f >= 0.0:
                f = round(f, digits) + rng.uniform(0.0, 10.0 ** (-digits))
            else:
                f = round(f, digits) - rng.uniform(0.0, 10.0 ** (-digits))
        elif self._name != FeatureName.PLAIN:
            raise NotImplementedError(f'Unknown feature: {self._name}.')
        return f

    def _set_default_options(self):
        """
        Set default options.
        """
        if self._name == FeatureName.PLAIN:
            self._options.setdefault(OptionKey.N_RUNS.value, 1)
        elif self._name == FeatureName.CUSTOM:
            if OptionKey.MODIFIER not in self._options:
                raise TypeError(f'When using a custom feature, you must specify the {OptionKey.MODIFIER} option.')
            self._options.setdefault(OptionKey.N_RUNS.value, 1)
        elif self._name == FeatureName.NOISY:
            self._options.setdefault(OptionKey.DISTRIBUTION.value, lambda rng: 1e-3 * rng.standard_normal())
            self._options.setdefault(OptionKey.N_RUNS.value, 10)
            self._options.setdefault(OptionKey.TYPE.value, NoiseType.RELATIVE.value)
        elif self._name == FeatureName.REGULARIZED:
            self._options.setdefault(OptionKey.N_RUNS.value, 1)
            self._options.setdefault(OptionKey.ORDER.value, 2)
            self._options.setdefault(OptionKey.PARAMETER.value, 1.0)
        elif self._name == FeatureName.TOUGH:
            self._options.setdefault(OptionKey.N_RUNS.value, 10)
            self._options.setdefault(OptionKey.RATE_ERROR.value, 0.0)
            self._options.setdefault(OptionKey.RATE_NAN.value, 0.05)
        elif self._name == FeatureName.TRUNCATED:
            self._options.setdefault(OptionKey.N_RUNS.value, 10)
            self._options.setdefault(OptionKey.SIGNIFICANT_DIGITS.value, 6)
        else:
            raise NotImplementedError(f'Unknown feature: {self._name}.')

    @staticmethod
    def _default_rng(seed, *args):
        """
        Generate a random number generator.

        Parameters
        ----------
        seed : int
            Seed used to generate an initial random number generator.
        args : tuple of int or float
            Arguments used to generate the returned random number generator.

        Returns
        -------
        `numpy.random.Generator`
            Random number generator.
        """
        rng = np.random.default_rng(seed)
        return np.random.default_rng(int(1e9 * abs(np.sin(1e5 * rng.standard_normal()) + sum(np.sin(1e5 * arg) for arg in args))))
