"""
Outputs of custom callbacks inside a composition are normalized once, at the
custom stage, before any later stage sees them.

Policy: an objective output follows the ``Problem.fun`` scalar policy (a
real scalar or a one-element real array converts to ``float``; anything
else, including complex values, is logged and recorded as NaN). A nonlinear
constraint output must be a real one-dimensional array of the predecessor's
channel size; anything else raises ``ValueError`` naming the stage,
occurrence and callback. Construction callbacks are validated the same way.
The single custom feature keeps its unchanged compatibility behaviour.
"""

import logging

import numpy as np
import pytest

from optiprofiler.opclasses import Feature, FeaturedProblem, Problem

X0 = np.array([0.3, -0.2, 0.5])


def quadratic(x):
    x = np.asarray(x, dtype=float)
    return float(x @ x)


def cub_two(x):
    return np.array([x[0] - 1.0, x[1] + 2.0])


def ceq_one(x):
    return np.array([x[2] - 0.5])


def problem():
    return Problem(quadratic, X0, xl=-np.ones(3), xu=np.ones(3), cub=cub_two, ceq=ceq_one, name='boundary')


SHAPES = {
    'list': lambda values: list(values),
    'tuple': lambda values: tuple(values),
    'ndarray': lambda values: np.array(values),
    'row_matrix': lambda values: np.array(values).reshape(1, -1),
    'float32': lambda values: np.array(values, dtype=np.float32),
}

DOWNSTREAM = {
    'noisy': dict(noise_level=0.0),
    'truncated': dict(significant_digits=12),
    'random_nan': dict(nan_rate=0.0),
    'nonquantifiable_constraints': {},
    'quantized': dict(mesh_size=1e-9, ground_truth=False),
    'unrelaxable_constraints': dict(unrelaxable_bounds=False, unrelaxable_nonlinear_constraints=True),
    'custom': {},
}


def make_mods(shape):
    def mod_cub(x, rng, prob):
        return SHAPES[shape](prob.cub(x) + 0.25)

    def mod_ceq(x, rng, prob):
        return SHAPES[shape](prob.ceq(x) - 0.25)

    return mod_cub, mod_ceq


class TestConstraintOutputs:

    @pytest.mark.parametrize('stage', sorted(DOWNSTREAM))
    @pytest.mark.parametrize('shape', sorted(SHAPES))
    def test_shapes_are_normalized_before_the_next_stage(self, shape, stage):
        mod_cub, mod_ceq = make_mods(shape)
        options = dict(mod_cub=mod_cub, mod_ceq=mod_ceq, **DOWNSTREAM[stage])
        featured = FeaturedProblem(problem(), Feature(f'custom+{stage}', **options), 10, 0)
        cub = featured.cub(X0)
        ceq = featured.ceq(X0)
        assert isinstance(cub, np.ndarray) and cub.dtype == float and cub.shape == (2,)
        assert isinstance(ceq, np.ndarray) and ceq.dtype == float and ceq.shape == (1,)
        if stage == 'nonquantifiable_constraints':
            assert set(cub.tolist()) <= {0.0, 1.0} and set(ceq.tolist()) <= {0.0, 1.0}
        elif stage == 'custom':
            # The second custom stage applies the same callbacks to the first stage's observed values.
            np.testing.assert_allclose(cub, cub_two(X0) + 0.5)
            np.testing.assert_allclose(ceq, ceq_one(X0) - 0.5)
        else:
            # A float32 callback output is exact only to single precision.
            tolerance = 1e-6 if shape == 'float32' else 1e-9
            np.testing.assert_allclose(cub, cub_two(X0) + 0.25, rtol=tolerance)
            np.testing.assert_allclose(ceq, ceq_one(X0) - 0.25, rtol=tolerance)
        # The gate reads the observed constraints of the custom stage during an objective query.
        assert np.isfinite(featured.fun(X0)) or stage == 'unrelaxable_constraints'

    def test_gate_reads_normalized_custom_constraints(self):
        def violated_list(x, rng, prob):
            return [0.5, -1.0]

        feature = Feature('custom+unrelaxable_constraints', mod_cub=violated_list, unrelaxable_bounds=False,
                          unrelaxable_nonlinear_constraints=True)
        assert FeaturedProblem(problem(), feature, 10, 0).fun(X0) == np.inf

    @pytest.mark.parametrize('bad, fragment', [
        (lambda x, rng, prob: [0.5], 'size 2'),
        (lambda x, rng, prob: np.zeros((2, 2)), 'one-dimensional'),
        (lambda x, rng, prob: np.array([1.0 + 1.0j, 2.0]), 'complex'),
        (lambda x, rng, prob: ['a', 'b'], 'real'),
        (lambda x, rng, prob: None, 'size 2'),
    ])
    def test_malformed_outputs_name_stage_occurrence_and_callback(self, bad, fragment):
        featured = FeaturedProblem(problem(), Feature('noisy+custom+noisy', mod_cub=bad, noise_level=0.0), 10, 0)
        with pytest.raises(ValueError) as excinfo:
            featured.cub(X0)
        message = str(excinfo.value)
        assert "stage 2 'custom' (occurrence 1)" in message and '`mod_cub`' in message and fragment in message

    def test_user_buffers_are_not_modified(self):
        buffer = np.array([0.5, 0.5])

        def shared_buffer(x, rng, prob):
            return buffer

        featured = FeaturedProblem(problem(), Feature('custom+random_nan', mod_cub=shared_buffer, nan_rate=1.0), 10, 0)
        assert np.all(np.isnan(featured.cub(X0)))
        np.testing.assert_array_equal(buffer, [0.5, 0.5])


class TestObjectiveOutputs:

    @pytest.mark.parametrize('value, expected', [
        (3, 3.0), (np.float32(1.5), 1.5), (np.array(2.5), 2.5), (np.array([4.0]), 4.0), (np.int64(7), 7.0),
    ])
    def test_real_scalars_convert_to_float(self, value, expected):
        featured = FeaturedProblem(problem(), Feature('custom+noisy', mod_fun=lambda x, rng, prob: value, noise_level=0.0), 10, 0)
        result = featured.fun(X0)
        assert isinstance(result, float) and result == expected

    @pytest.mark.parametrize('value', [1.0 + 2.0j, np.complex128(1.0 + 0.0j), np.array([1.0, 2.0]), [1.0], 'text', None])
    def test_non_real_scalars_are_recorded_as_nan_with_a_warning(self, value, caplog):
        featured = FeaturedProblem(problem(), Feature('custom+noisy', mod_fun=lambda x, rng, prob: value, noise_level=0.0), 10, 0)
        with caplog.at_level(logging.WARNING, logger='optiprofiler'):
            result = featured.fun(X0)
        assert np.isnan(result)
        assert any("stage 1 'custom' (occurrence 1)" in record.getMessage() and '`mod_fun`' in record.getMessage()
                   for record in caplog.records)
        # The scoring reference is untouched by the custom observation.
        assert featured.fun_hist[0] == quadratic(X0)


class TestConstructionOutputs:

    def test_list_outputs_become_float_arrays(self):
        def mod_x0(rng, prob):
            return [0.1, 0.2, 0.3]

        def mod_bounds(rng, prob):
            return [-2, -2, -2], (2, 2, 2)

        def mod_affine(rng, prob):
            return [[2.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]], [0, 0, 0], [[0.5, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]

        feature = Feature('custom+truncated', mod_x0=mod_x0, mod_bounds=mod_bounds, mod_affine=mod_affine, significant_digits=12)
        featured = FeaturedProblem(problem(), feature, 10, 0)
        np.testing.assert_array_equal(featured.x0, [0.1, 0.2, 0.3])
        assert featured.xl.dtype == float and featured.xu.dtype == float
        np.testing.assert_array_equal(featured.xl, [-2.0, -2.0, -2.0])
        assert featured.fun(np.array([0.5, 0.0, 0.0])) == quadratic([1.0, 0.0, 0.0])

    def test_array_pairs_and_triples_are_accepted_as_before(self):
        # The legacy modifiers tuple-unpack whatever a callback returns, so a
        # (2, n) array of bounds or a list of linear constraint parts is valid.
        def mod_bounds(rng, prob):
            return np.array([[-2.0, -2.0, -2.0], [2.0, 2.0, 2.0]])

        def mod_linear_ub(rng, prob):
            return [np.array([[1.0, 1.0, 0.0]]), np.array([5.0])]

        def mod_affine(rng, prob):
            return (np.diag([2.0, 1.0, 1.0]), np.zeros(3), np.diag([0.5, 1.0, 1.0]))

        feature = Feature('custom+truncated', mod_bounds=mod_bounds, mod_linear_ub=mod_linear_ub, mod_affine=mod_affine,
                          significant_digits=12)
        featured = FeaturedProblem(problem(), feature, 10, 0)
        np.testing.assert_array_equal(featured.xl, [-2.0, -2.0, -2.0])
        np.testing.assert_array_equal(featured.xu, [2.0, 2.0, 2.0])
        # A custom linear-constraint callback replaces the constraints verbatim, as in the single feature.
        np.testing.assert_array_equal(featured.aub, [[1.0, 1.0, 0.0]])
        np.testing.assert_array_equal(featured.bub, [5.0])
        assert featured.fun(np.array([0.5, 0.0, 0.0])) == quadratic([1.0, 0.0, 0.0])

    @pytest.mark.parametrize('key', ['mod_bounds', 'mod_linear_ub', 'mod_linear_eq', 'mod_affine'])
    def test_ndarray_containers_are_accepted_for_every_construction_callback(self, key):
        # Valid components packed in ndarray containers (a (2, n) float array for the
        # bounds, object arrays for the others), as the legacy modifiers unpack them.
        def container(*parts):
            packed = np.empty(len(parts), dtype=object)
            for i, part in enumerate(parts):
                packed[i] = part
            return packed

        callbacks = {
            'mod_bounds': lambda rng, prob: np.array([[-2.0, -2.0, -2.0], [2.0, 2.0, 2.0]]),
            'mod_linear_ub': lambda rng, prob: container(np.array([[1.0, 1.0, 0.0]]), np.array([5.0])),
            'mod_linear_eq': lambda rng, prob: container(np.array([[0.0, 1.0, -1.0]]), np.array([-0.5])),
            'mod_affine': lambda rng, prob: container(np.diag([2.0, 1.0, 1.0]), np.zeros(3), np.diag([0.5, 1.0, 1.0])),
        }
        featured = FeaturedProblem(problem(), Feature('custom+noisy', noise_level=0.0, **{key: callbacks[key]}), 10, 0)
        if key == 'mod_bounds':
            np.testing.assert_array_equal(featured.xl, [-2.0, -2.0, -2.0])
            np.testing.assert_array_equal(featured.xu, [2.0, 2.0, 2.0])
        elif key == 'mod_linear_ub':
            np.testing.assert_array_equal(featured.aub, [[1.0, 1.0, 0.0]])
            np.testing.assert_array_equal(featured.bub, [5.0])
        elif key == 'mod_linear_eq':
            np.testing.assert_array_equal(featured.aeq, [[0.0, 1.0, -1.0]])
            np.testing.assert_array_equal(featured.beq, [-0.5])
        else:
            assert featured.fun(np.array([0.5, 0.0, 0.0])) == quadratic([1.0, 0.0, 0.0])
            np.testing.assert_array_equal(featured.x0, X0 / np.array([2.0, 1.0, 1.0]))

    def test_generators_unpack_like_the_legacy_modifiers(self):
        def mod_bounds(rng, prob):
            return (np.full(3, value) for value in (-2.0, 2.0))

        def mod_affine(rng, prob):
            return iter((np.eye(3), np.zeros(3), np.eye(3)))

        featured = FeaturedProblem(problem(), Feature('custom+truncated', mod_bounds=mod_bounds, mod_affine=mod_affine), 10, 0)
        np.testing.assert_array_equal(featured.xl, [-2.0, -2.0, -2.0])
        np.testing.assert_array_equal(featured.xu, [2.0, 2.0, 2.0])

        def three_parts(rng, prob):
            return (np.zeros(3) for _ in range(3))

        with pytest.raises(ValueError, match="`mod_bounds` of stage 1 'custom' \\(occurrence 1\\) must return a pair"):
            FeaturedProblem(problem(), Feature('custom+truncated', mod_bounds=three_parts), 10, 0)

    @pytest.mark.parametrize('key, bad, fragment', [
        ('mod_x0', lambda rng, prob: [0.1, 0.2], 'size 3'),
        ('mod_x0', lambda rng, prob: np.array([1.0 + 1.0j, 0.0, 0.0]), 'complex'),
        ('mod_bounds', lambda rng, prob: ([0.0, 0.0], [1.0, 1.0, 1.0]), 'size 3'),
        ('mod_linear_ub', lambda rng, prob: (np.ones((2, 3)), np.ones(3)), 'size 2'),
        ('mod_linear_ub', lambda rng, prob: ([[1.0, 2.0, 3.0], [1.0]], np.ones(2)), 'real'),
        ('mod_affine', lambda rng, prob: (np.eye(2), np.zeros(3), np.eye(3)), 'shape (3, 3)'),
    ])
    def test_malformed_construction_outputs_name_stage_and_callback(self, key, bad, fragment):
        with pytest.raises(ValueError) as excinfo:
            FeaturedProblem(problem(), Feature(f'custom+truncated', **{key: bad}), 10, 0)
        message = str(excinfo.value)
        assert "stage 1 'custom' (occurrence 1)" in message and f'`{key}`' in message and fragment in message


def top_level_mod_fun(x, rng, problem):
    return problem.fun(x) + 0.5


def top_level_mod_cub(x, rng, problem):
    return list(problem.cub(x) + 0.25)


def top_level_mod_x0(rng, problem):
    return problem.x0 + 0.1


class TestPicklability:

    def test_composed_problem_with_custom_callbacks_survives_pickling(self):
        import pickle

        feature = Feature('custom+noisy', mod_fun=top_level_mod_fun, mod_cub=top_level_mod_cub, mod_x0=top_level_mod_x0,
                          noise_type='absolute', noise_level=1.0)
        original = FeaturedProblem(problem(), feature, 10, 5)
        original.fun(X0)
        clone = pickle.loads(pickle.dumps(original))
        np.testing.assert_array_equal(clone.x0, original.x0)
        np.testing.assert_array_equal(clone.fun_hist, original.fun_hist)
        # Both continue from the same served-query state, so their observations agree.
        for x in (X0, X0 + 0.1, X0):
            assert clone.fun(x) == original.fun(x)
            np.testing.assert_array_equal(clone.cub(x), original.cub(x))
        np.testing.assert_array_equal(clone.fun_hist, original.fun_hist)
        assert clone.n_eval_fun == original.n_eval_fun == 4


class TestLegacySingleCustomUnchanged:

    def test_single_custom_returns_callback_output_verbatim(self):
        def list_cub(x, rng, prob):
            return list(prob.cub(x))

        featured = FeaturedProblem(problem(), Feature('custom', mod_cub=list_cub), 10, 0)
        result = featured.cub(X0)
        assert isinstance(result, list)
