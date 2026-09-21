"""Serialization owns runtime caches, not a user's unrelated object graph."""

import copy
import copyreg
import io
import pickle
import weakref
from multiprocessing.reduction import ForkingPickler

import numpy as np
import pytest

from optiprofiler import Feature, FeaturedProblem, Problem
from optiprofiler.composition import ComposedFeaturedProblem
from optiprofiler.opclasses import _restore_featured_problem


def objective(x):
    return float(x @ x)


class Metadata:
    pass


class AffineCallback:
    def __init__(self):
        self.calls = 0

    def __call__(self, rng, problem):
        self.calls += 1
        scale = float(self.calls + 1)
        return np.eye(problem.n) * scale, np.zeros(problem.n), np.eye(problem.n) / scale


class PrivateSlotMixin:
    __slots__ = ('__private', 'parent_metadata')

    def set_parent_private(self, value):
        self.__private = value

    def get_parent_private(self):
        return self.__private


class SlottedSingle(FeaturedProblem, PrivateSlotMixin):
    __slots__ = ('__private', 'child_metadata', 'unset_metadata')

    def set_child_private(self, value):
        self.__private = value

    def get_child_private(self):
        return self.__private


class SlottedComposed(ComposedFeaturedProblem, PrivateSlotMixin):
    __slots__ = ('__private', 'child_metadata', 'unset_metadata')

    def set_child_private(self, value):
        self.__private = value

    def get_child_private(self):
        return self.__private


class FirstReducerPickler(pickle.Pickler):
    """The interim reduction saved state among the reconstructor arguments."""

    def reducer_override(self, obj):
        if isinstance(obj, FeaturedProblem):
            return _restore_featured_problem, (type(obj), obj.__dict__)
        return NotImplemented


class ConstructorStatePickler(pickle.Pickler):
    """Compositions before the explicit reducer used __new__ arguments and state."""

    def reducer_override(self, obj):
        if isinstance(obj, ComposedFeaturedProblem):
            return copyreg.__newobj__, (type(obj), obj._problem, obj._runtime, obj._max_eval, obj._seed), obj.__dict__
        return NotImplemented


MODES = ['copy', 'deepcopy', 'protocol2', 'protocol3', 'protocol4', 'protocol5', 'multiprocessing']


def roundtrip(value, mode):
    if mode == 'copy':
        return copy.copy(value)
    if mode == 'deepcopy':
        return copy.deepcopy(value)
    if mode == 'multiprocessing':
        return pickle.loads(ForkingPickler.dumps(value))
    return pickle.loads(pickle.dumps(value, protocol=int(mode[-1])))


def build(composed=False, cls=None, callback=None):
    feature = Feature('plain') if callback is None else Feature('custom', mod_affine=callback)
    if composed:
        feature = Feature([{'name': 'custom', 'options': {'mod_affine': callback}}, 'noisy']) if callback else \
            Feature(['perturbed_x0', 'noisy'])
    return (cls or FeaturedProblem)(Problem(objective, [1.0, 2.0]), feature, 5, 3)


def affine_runtime(featured):
    return featured._views[0]._runtime if isinstance(featured, ComposedFeaturedProblem) else featured._runtime


@pytest.mark.parametrize('composed', [False, True], ids=['single', 'composed'])
@pytest.mark.parametrize('mode', MODES)
def test_restoring_a_runtime_does_not_visit_global_function_metadata(composed, mode):
    # Functions pickle by global name. Their metadata is neither serialized nor
    # owned by the featured problem, and must not enter its restoration walk.
    metadata = Metadata()
    objective.diagnostic = weakref.proxy(metadata)
    del metadata
    try:
        featured = build(composed)
        restored = roundtrip(featured, mode)
        assert type(restored) is type(featured)
        np.testing.assert_array_equal(restored.x0, featured.x0)
    finally:
        del objective.diagnostic


@pytest.mark.parametrize('composed', [False, True], ids=['single', 'composed'])
@pytest.mark.parametrize('mode', MODES)
def test_slots_from_the_class_and_its_other_base_are_not_discarded(composed, mode):
    cls = SlottedComposed if composed else SlottedSingle
    featured = build(composed, cls)
    shared = {'value': 3}
    featured.parent_metadata = shared
    featured.child_metadata = shared
    featured.set_parent_private('parent')
    featured.set_child_private('child')
    restored = roundtrip(featured, mode)
    assert type(restored) is cls
    assert restored.parent_metadata == shared
    assert restored.parent_metadata is restored.child_metadata
    assert restored.get_parent_private() == 'parent'
    assert restored.get_child_private() == 'child'
    assert not hasattr(restored, 'unset_metadata')
    if mode == 'copy':
        assert restored.parent_metadata is shared
    else:
        assert restored.parent_metadata is not shared


@pytest.mark.parametrize('composed', [False, True], ids=['single', 'composed'])
@pytest.mark.parametrize('mode', MODES[1:])
def test_runtime_restores_its_own_cache_without_resampling_or_breaking_cycles(composed, mode):
    callback = AffineCallback()
    featured = build(composed, callback=callback)
    callback.owner = featured
    restored = roundtrip(featured, mode)
    runtime = affine_runtime(restored)
    saved_callback = runtime._options['mod_affine']
    saved_callback = getattr(saved_callback, 'user', saved_callback)
    assert saved_callback.owner is restored
    assert saved_callback.calls == callback.calls == 1
    assert all(not value.flags.writeable for value in runtime._kept_affine[2])
    assert runtime._kept_affine[0] is (restored._views[0]._predecessor if composed else restored._problem)
    assert restored._problem._x0.flags.writeable  # user problem data is not frozen with the internal cache


@pytest.mark.parametrize('composed', [False, True], ids=['single', 'composed'])
@pytest.mark.parametrize('cyclic', [False, True], ids=['acyclic', 'cyclic'])
def test_first_reducer_preserves_the_state_dictionary_while_pickle_finishes_it(composed, cyclic):
    callback = AffineCallback()
    featured = build(composed, callback=callback)
    if cyclic:
        callback.owner = featured
    stream = io.BytesIO()
    FirstReducerPickler(stream, protocol=4).dump(featured)
    restored = pickle.loads(stream.getvalue())
    np.testing.assert_array_equal(restored.x0, featured.x0)
    runtime = affine_runtime(restored)
    saved_callback = runtime._options['mod_affine']
    saved_callback = getattr(saved_callback, 'user', saved_callback)
    assert saved_callback.calls == callback.calls == 1
    if cyclic:
        assert saved_callback.owner is restored
    assert all(not value.flags.writeable for value in runtime._kept_affine[2])


def test_shallow_copy_does_not_share_the_instance_dictionary():
    featured = build()
    restored = copy.copy(featured)
    restored.independent_attribute = 3
    assert not hasattr(featured, 'independent_attribute')


def test_earlier_composition_with_constructor_arguments_keeps_its_built_cache():
    callback = AffineCallback()
    featured = build(True, callback=callback)
    stream = io.BytesIO()
    ConstructorStatePickler(stream, protocol=4).dump(featured)
    restored = pickle.loads(stream.getvalue())
    np.testing.assert_array_equal(restored.x0, featured.x0)
    runtime = affine_runtime(restored)
    assert runtime._options['mod_affine'].user.calls == callback.calls == 1
    assert all(not array.flags.writeable for array in runtime._kept_affine[2])
