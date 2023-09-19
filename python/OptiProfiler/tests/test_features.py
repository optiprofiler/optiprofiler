import numpy as np
import pytest

from OptiProfiler.features import Feature


class TestFeature:

    @staticmethod
    def rosen(x):
        return np.sum(1e2 * (x[1:] - x[:-1] ** 2) ** 2 + (1.0 - x[:-1]) ** 2)

    @pytest.mark.parametrize('n', [1, 10, 100])
    @pytest.mark.parametrize('seed', [0, 1, 2])
    def test_plain(self, n, seed):
        # Generate random data.
        rng = np.random.default_rng(seed)
        x = rng.standard_normal(n)
        f = self.rosen(x)

        # Test the plain feature.
        feature = Feature('plain')
        assert feature.name == 'plain'
        assert feature.options == {}
        assert feature.modifier(x, f) == f

    @pytest.mark.parametrize('n', [1, 10, 100])
    @pytest.mark.parametrize('seed', [0, 1, 2])
    def test_regularized(self, n, seed):
        # Generate random data.
        rng = np.random.default_rng(seed)
        x = rng.standard_normal(n)
        f = self.rosen(x)

        # Test the regularized feature.
        feature = Feature('regularized')
        assert feature.name == 'regularized'
        assert feature.options == {'parameter': 1.0, 'order': 2}
        np.testing.assert_allclose(feature.modifier(x, f),  f + np.linalg.norm(x))

        # Add custom options.
        feature = Feature('regularized', parameter=2.0, order=3)
        assert feature.name == 'regularized'
        assert feature.options == {'parameter': 2.0, 'order': 3}
        np.testing.assert_allclose(feature.modifier(x, f),  f + 2.0 * np.linalg.norm(x, 3))

    @pytest.mark.parametrize('n', [1, 10, 100])
    @pytest.mark.parametrize('seed', [0, 1, 2])
    def test_noisy(self, n, seed):
        # Generate random data.
        rng = np.random.default_rng(seed)
        x = rng.standard_normal(n)
        f = self.rosen(x)

        # Test the noisy feature.
        feature = Feature('noisy')
        assert feature.name == 'noisy'
        assert 'distribution' in feature.options and feature.options['type'] == 'relative'
        feature.modifier(x, f)

        # Add custom options.
        feature = Feature('noisy', distribution=lambda rng: 1.0, type='absolute')
        assert feature.name == 'noisy'
        assert 'distribution' in feature.options and feature.options['type'] == 'absolute'
        np.testing.assert_allclose(feature.modifier(x, f),  f + 1.0)

    @pytest.mark.parametrize('n', [1, 10, 100])
    @pytest.mark.parametrize('seed', [0, 1, 2])
    def test_truncated(self, n, seed):
        # Generate random data.
        rng = np.random.default_rng(seed)
        x = rng.standard_normal(n)
        f = self.rosen(x)

        # Test the truncated feature.
        feature = Feature('truncated')
        assert feature.name == 'truncated'
        assert feature.options == {'significant_digits': 6}
        np.testing.assert_allclose(feature.modifier(x, f),  f, 1e-5, 1e-5)
        np.testing.assert_allclose(feature.modifier(x, -f),  -f, 1e-5, 1e-5)

        # Add custom options.
        feature = Feature('truncated', significant_digits=4)
        assert feature.name == 'truncated'
        assert feature.options == {'significant_digits': 4}
        np.testing.assert_allclose(feature.modifier(x, f),  f, 1e-3, 1e-3)

    @pytest.mark.parametrize('n', [1, 10, 100])
    @pytest.mark.parametrize('seed', [0, 1, 2])
    def test_tough(self, n, seed):
        # Generate random data.
        rng = np.random.default_rng(seed)
        x = rng.standard_normal(n)
        f = self.rosen(x)

        # Test the tough feature.
        feature = Feature('tough')
        assert feature.name == 'tough'
        assert feature.options == {'rate_error': 0.0, 'rate_nan': 0.05}
        f_tough = feature.modifier(x, f)
        assert f_tough == f or np.isnan(f_tough)

        # Add custom options.
        feature = Feature('tough', rate_error=0.5, rate_nan=0.5)
        assert feature.name == 'tough'
        assert feature.options == {'rate_error': 0.5, 'rate_nan': 0.5}
        try:
            f_tough = feature.modifier(x, f)
            assert f_tough == f or np.isnan(f_tough)
        except RuntimeError:
            pass

    @pytest.mark.parametrize('n', [1, 10, 100])
    @pytest.mark.parametrize('seed', [0, 1, 2])
    def test_custom(self, n, seed):
        # Generate random data.
        rng = np.random.default_rng(seed)
        x = rng.standard_normal(n)
        f = self.rosen(x)

        # Test the custom feature.
        feature = Feature('custom', modifier=lambda x, f, seed: f + 1.0)
        assert feature.name == 'custom'
        assert 'modifier' in feature.options
        np.testing.assert_allclose(feature.modifier(x, f),  f + 1.0)

    def test_exceptions(self):
        with pytest.raises(ValueError):
            Feature('unknown')
        with pytest.raises(ValueError):
            Feature('regularized', unknown=1.0)
        with pytest.raises(ValueError):
            Feature('plain', parameter=1.0)
        with pytest.raises(TypeError):
            Feature('custom', modifier=1.0)
        with pytest.raises(TypeError):
            Feature('noisy', distribution=1.0)
        with pytest.raises(TypeError):
            Feature('regularized', order='1.0')
        with pytest.raises(TypeError):
            Feature('regularized', parameter='1.0')
        with pytest.raises(TypeError):
            Feature('regularized', parameter=-1.0)
        with pytest.raises(TypeError):
            Feature('tough', rate_error='1.0')
        with pytest.raises(TypeError):
            Feature('tough', rate_error=-1.0)
        with pytest.raises(TypeError):
            Feature('tough', rate_error=2.0)
        with pytest.raises(TypeError):
            Feature('tough', rate_nan='1.0')
        with pytest.raises(TypeError):
            Feature('tough', rate_nan=-1.0)
        with pytest.raises(TypeError):
            Feature('tough', rate_nan=2.0)
        with pytest.raises(TypeError):
            Feature('truncated', significant_digits=2.5)
        with pytest.raises(TypeError):
            Feature('noisy', type='+')
        with pytest.raises(TypeError):
            Feature('custom')

    def test_catch(self):
        # The value significant_digits can be a float.
        Feature('truncated', significant_digits=2.0)
