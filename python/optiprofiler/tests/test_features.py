import numpy as np
import pytest

from optiprofiler.opclasses import Feature, Problem, FeaturedProblem


class TestFeature:
    """Test the Feature class."""

    @staticmethod
    def rosen(x):
        """Rosenbrock function."""
        return np.sum(1e2 * (x[1:] - x[:-1] ** 2) ** 2 + (1.0 - x[:-1]) ** 2)

    @pytest.mark.parametrize('n', [2, 10])
    def test_plain(self, n):
        """Test the plain feature."""
        feature = Feature('plain')
        assert feature.name == 'plain'
        assert feature.options == {'n_runs': 1}
        assert not feature.is_stochastic

        # Test with a problem
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)
        featured_problem = FeaturedProblem(problem, feature, 500)
        f = featured_problem.fun(x0)
        assert f == self.rosen(x0)

    @pytest.mark.parametrize('n', [2, 10])
    def test_noisy(self, n):
        """Test the noisy feature."""
        # Default noisy feature
        feature = Feature('noisy')
        assert feature.name == 'noisy'
        assert 'noise_level' in feature.options
        assert feature.options['n_runs'] == 5  # Default is 5
        assert feature.is_stochastic

        # Noisy feature with custom options
        feature = Feature('noisy', noise_level=0.01, noise_type='relative', n_runs=5)
        assert feature.name == 'noisy'
        assert feature.options['noise_level'] == 0.01
        assert feature.options['noise_type'] == 'relative'
        assert feature.options['n_runs'] == 5

        # Test with a problem
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)
        featured_problem = FeaturedProblem(problem, feature, 500)
        f_noisy = featured_problem.fun(x0)
        # The noisy value should be different from the true value (with high probability)
        # but we can't guarantee it, so we just check that it's a number
        assert isinstance(f_noisy, (int, float))

    @pytest.mark.parametrize('n', [2, 10])
    def test_truncated(self, n):
        """Test the truncated feature."""
        feature = Feature('truncated')
        assert feature.name == 'truncated'
        assert 'significant_digits' in feature.options
        
        # Test with custom options
        feature = Feature('truncated', significant_digits=4)
        assert feature.options['significant_digits'] == 4

        # Test with a problem
        x0 = np.ones(n)
        problem = Problem(self.rosen, x0)
        featured_problem = FeaturedProblem(problem, feature, 500)
        f_truncated = featured_problem.fun(x0)
        f_true = self.rosen(x0)
        # Truncated value should be close to true value
        np.testing.assert_allclose(f_truncated, f_true, rtol=1e-3)

    @pytest.mark.parametrize('n', [2, 10])
    def test_random_nan(self, n):
        """Test the random_nan feature."""
        feature = Feature('random_nan')
        assert feature.name == 'random_nan'
        assert 'nan_rate' in feature.options
        assert feature.is_stochastic

        # Test with custom nan rate
        feature = Feature('random_nan', nan_rate=0.5)
        assert feature.options['nan_rate'] == 0.5

    @pytest.mark.parametrize('n', [2, 10])
    def test_perturbed_x0(self, n):
        """Test the perturbed_x0 feature."""
        feature = Feature('perturbed_x0')
        assert feature.name == 'perturbed_x0'
        assert feature.is_stochastic

        # Test with a problem
        x0 = np.zeros(n)
        problem = Problem(self.rosen, x0)
        featured_problem = FeaturedProblem(problem, feature, 500, seed=42)
        # The perturbed x0 should be different from the original
        assert not np.allclose(featured_problem.x0, x0)

    def test_linearly_transformed(self):
        """Test the linearly_transformed feature."""
        feature = Feature('linearly_transformed')
        assert feature.name == 'linearly_transformed'
        
        # Test with rotation
        feature = Feature('linearly_transformed', rotated=True)
        assert feature.options['rotated'] is True
        assert feature.is_stochastic

    def test_permuted(self):
        """Test the permuted feature."""
        feature = Feature('permuted')
        assert feature.name == 'permuted'
        assert feature.is_stochastic

    def test_quantized(self):
        """Test the quantized feature."""
        feature = Feature('quantized')
        assert feature.name == 'quantized'
        
        # Test with custom mesh size
        feature = Feature('quantized', mesh_size=0.01)
        assert feature.options['mesh_size'] == 0.01

    def test_custom(self):
        """Test the custom feature."""
        # Custom feature requires at least one modifier
        # The signature is: (x, rng, problem) -> float
        def custom_mod_fun(x, rng, problem):
            f = problem.fun(x)
            return f + 1.0
        
        feature = Feature('custom', mod_fun=custom_mod_fun)
        assert feature.name == 'custom'
        assert feature.is_stochastic
        
        # Test with a problem
        x0 = np.zeros(2)
        problem = Problem(self.rosen, x0)
        featured_problem = FeaturedProblem(problem, feature, 500)
        f = featured_problem.fun(x0)
        expected = self.rosen(x0) + 1.0
        np.testing.assert_allclose(f, expected)

    def test_exceptions(self):
        """Test that appropriate exceptions are raised."""
        # Invalid feature name type
        with pytest.raises(TypeError):
            Feature(1)
        
        # Unknown feature name
        with pytest.raises(ValueError):
            Feature('unknown')
        
        # Invalid option for feature
        with pytest.raises(ValueError):
            Feature('plain', noise_level=0.1)
        
        # Invalid n_runs type
        with pytest.raises(TypeError):
            Feature('noisy', n_runs=1.5)
        
        # Invalid n_runs value
        with pytest.raises(ValueError):
            Feature('noisy', n_runs=-1)
        
        # Invalid nan_rate type
        with pytest.raises(TypeError):
            Feature('random_nan', nan_rate='1.0')
        
        # Invalid nan_rate value (out of range)
        with pytest.raises(ValueError):
            Feature('random_nan', nan_rate=-1.0)
        with pytest.raises(ValueError):
            Feature('random_nan', nan_rate=2.0)
        
        # Invalid significant_digits type
        with pytest.raises(TypeError):
            Feature('truncated', significant_digits=2.5)
        
        # Invalid significant_digits value
        with pytest.raises(ValueError):
            Feature('truncated', significant_digits=-1)
        
        # Invalid noise_type
        with pytest.raises(TypeError):
            Feature('noisy', noise_type=1)
        with pytest.raises(ValueError):
            Feature('noisy', noise_type='invalid')

    def test_catch(self):
        """Test edge cases that should work."""
        # n_runs can be a float if it's an integer value
        Feature('noisy', n_runs=2.0)

        # significant_digits can be a float if it's an integer value
        Feature('truncated', significant_digits=2.0)
