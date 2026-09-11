"""
Ownership of container and array option values.

A ``Feature`` is a pure specification: it never aliases data the caller still
owns. Ordinary containers and arrays given as option values are copied at
ingress (arrays as read-only copies), the inspection views hand out read-only
arrays and fresh container copies, and the native transport carries owned
copies. Callables and other objects keep their identity and no copy or
reduction hook of a user object is invoked by construction, inspection or
provenance.
"""

import pickle

import numpy as np
import pytest

from optiprofiler.opclasses import Feature, FeaturedProblem, Problem
from optiprofiler.provenance import describe_feature, effective_specification


def sphere(x):
    return float(np.dot(x, x))


class Hostile:
    """A callable whose copy and reduction hooks must never run during normalization or inspection."""

    def __call__(self, x, rng, problem):
        return problem.fun(x) + 1.0

    def __copy__(self):
        raise AssertionError('copy hook invoked')

    def __deepcopy__(self, memo):
        raise AssertionError('deepcopy hook invoked')

    def __reduce_ex__(self, protocol):
        raise AssertionError('reduction hook invoked')

    def __getstate__(self):
        raise AssertionError('state hook invoked')


class TestArrayAndListValues:

    def test_array_option_is_copied_at_ingress_and_read_only_in_views(self):
        level = np.array([1e-3, 2e-3])
        feature = Feature('perturbed_x0', perturbation_level=level)
        stored = feature.stages[0].options['perturbation_level']
        assert stored is not level and not stored.flags.writeable
        np.testing.assert_array_equal(stored, [1e-3, 2e-3])
        level[0] = 5.0   # the caller's array is theirs; the specification does not follow it
        np.testing.assert_array_equal(feature.stages[0].options['perturbation_level'], [1e-3, 2e-3])
        with pytest.raises(ValueError, match='read-only'):
            feature.stages[0].options['perturbation_level'][0] = 7.0
        declared = feature.declared.entries[0][1]['perturbation_level']
        assert not declared.flags.writeable
        np.testing.assert_array_equal(declared, [1e-3, 2e-3])
        # The exposed array is an isolated copy: re-enabling writes on it and
        # mutating it cannot reach the specification (no shared base either).
        for exposed in (feature.stages[0].options['perturbation_level'], feature.declared.entries[0][1]['perturbation_level'],
                        effective_specification(feature)[0]['options']['perturbation_level']):
            assert exposed.base is None
            exposed.setflags(write=True)
            exposed[0] = 9.0
        np.testing.assert_array_equal(feature.stages[0].options['perturbation_level'], [1e-3, 2e-3])
        np.testing.assert_array_equal(feature.declared.entries[0][1]['perturbation_level'], [1e-3, 2e-3])
        np.testing.assert_array_equal(effective_specification(feature)[0]['options']['perturbation_level'], [1e-3, 2e-3])

    def test_array_values_survive_transport_and_replay_as_independent_copies(self):
        feature = Feature('perturbed_x0', perturbation_level=np.array([1e-3, 2e-3]))
        clone = pickle.loads(pickle.dumps(feature))
        cloned = clone.stages[0].options['perturbation_level']
        assert not cloned.flags.writeable and cloned is not feature.stages[0].options['perturbation_level']
        np.testing.assert_array_equal(cloned, [1e-3, 2e-3])
        spec = effective_specification(feature)
        replay = spec[0]['options']['perturbation_level']
        assert replay is not feature.stages[0].options['perturbation_level']
        np.testing.assert_array_equal(replay, [1e-3, 2e-3])
        # The replayed specification declares the effective stages structurally; the stages agree.
        assert describe_feature(Feature(spec))['stages'] == describe_feature(feature)['stages']

    def test_numbers_do_not_change_when_the_callers_array_changes(self):
        problem = Problem(sphere, np.array([1.0, 2.0]))
        level = np.array([1e-2, 3e-2])
        feature = Feature('perturbed_x0', perturbation_level=level)
        before = FeaturedProblem(problem, feature, 5, 1).x0.copy()
        expected = FeaturedProblem(problem, Feature('perturbed_x0', perturbation_level=np.array([1e-2, 3e-2])), 5, 1).x0
        np.testing.assert_array_equal(before, expected)
        level[:] = 100.0
        np.testing.assert_array_equal(FeaturedProblem(problem, feature, 5, 1).x0, before)
        # The scalar path is unchanged by the ownership policy.
        scalar = FeaturedProblem(problem, Feature('perturbed_x0', perturbation_level=1e-2), 5, 1).x0
        np.testing.assert_array_equal(scalar, FeaturedProblem(problem, Feature('perturbed_x0', perturbation_level=1e-2), 5, 1).x0)

    def test_list_option_is_copied_and_views_are_fresh_copies(self):
        level = [1e-3, 2e-3]
        feature = Feature([{'name': 'perturbed_x0', 'options': {'perturbation_level': level}}])
        level.append(9.0)
        assert feature.stages[0].options['perturbation_level'] == [1e-3, 2e-3]
        seen = feature.stages[0].options['perturbation_level']
        seen.append(4.0)
        assert feature.stages[0].options['perturbation_level'] == [1e-3, 2e-3]
        entry = feature.declared.entries[0][1]
        assert entry['perturbation_level'] == [1e-3, 2e-3]
        with pytest.raises(TypeError):
            entry['perturbation_level'] = None   # the mapping view itself is read-only

    def test_own_and_view_policy_on_plain_values(self):
        from optiprofiler.feature_definitions import own_value, view_value
        array = np.arange(3.0)
        owned = own_value({'a': [array, (1, 2)], 'b': 'text', 'c': 3.5})
        assert owned['a'][0] is not array and not owned['a'][0].flags.writeable
        assert owned['a'][1] == (1, 2) and owned['b'] == 'text' and owned['c'] == 3.5
        viewed = view_value(owned)
        with pytest.raises(TypeError):
            viewed['b'] = 'other'
        # The viewed array is an isolated read-only copy of the owned one, never the owned object.
        assert viewed['a'][0] is not owned['a'][0] and viewed['a'][0].base is None
        assert not viewed['a'][0].flags.writeable
        np.testing.assert_array_equal(viewed['a'][0], owned['a'][0])
        hostile = Hostile()
        assert own_value(hostile) is hostile and view_value(hostile) is hostile
        assert own_value([hostile])[0] is hostile
        # Exact built-in containers only: subclasses keep their identity.
        class Tagged(list):
            pass
        tagged = Tagged([1, 2])
        assert own_value(tagged) is tagged and view_value(tagged) is tagged
        viewed_array = view_value(owned['a'][0])
        assert viewed_array is not owned['a'][0] and viewed_array.base is None and not viewed_array.flags.writeable


class TestOpaqueCallables:

    def test_hostile_callable_keeps_identity_and_no_hook_runs(self):
        hostile = Hostile()
        feature = Feature('custom', mod_fun=hostile)
        assert feature.stages[0].options['mod_fun'] is hostile
        assert feature.declared.entries[0][1]['mod_fun'] is hostile
        assert effective_specification(feature)[0]['options']['mod_fun'] is hostile
        described = describe_feature(feature)['stages'][0]['options']['mod_fun']
        # Provenance describes the callable by name: plain data, never the object itself.
        assert described is not hostile and isinstance(described, (dict, str)) and 'Hostile' in str(described)
        again = Feature(effective_specification(feature))
        assert again.stages[0].options['mod_fun'] is hostile
        reused = Feature(feature)
        assert reused.stages[0].options['mod_fun'] is hostile
        problem = Problem(sphere, np.array([1.0, 2.0]))
        assert FeaturedProblem(problem, feature, 5, 1).fun(np.array([1.0, 1.0])) == 3.0

    def test_callable_container_subclass_stays_the_opaque_callback(self):
        class CallableSamples(list):
            """A callable that is also a list: it must stay the callback, never be copied as a list."""

            def __call__(self, rng, size):
                return rng.standard_normal(size)

            def __iter__(self):
                raise AssertionError('iteration hook invoked')

        callback = CallableSamples()
        feature = Feature('noisy', distribution=callback)
        assert feature.stages[0].options['distribution'] is callback
        assert feature.declared.entries[0][1]['distribution'] is callback
        assert effective_specification(feature)[0]['options']['distribution'] is callback
        assert Feature(effective_specification(feature)).stages[0].options['distribution'] is callback
        problem = Problem(sphere, np.array([1.0, 2.0]))
        assert np.isfinite(FeaturedProblem(problem, feature, 5, 1).fun(np.array([1.0, 1.0])))

    def test_hostile_element_inside_a_container_keeps_identity(self):
        hostile = Hostile()
        feature = Feature('custom', mod_fun=hostile, mod_x0=lambda rng, problem: problem.x0)
        stored = feature.stages[0].options
        assert stored['mod_fun'] is hostile and callable(stored['mod_x0'])
