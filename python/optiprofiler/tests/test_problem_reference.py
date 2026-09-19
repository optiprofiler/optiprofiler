"""
Feasible reference fact of a problem (``Problem.reference``).

The record has exactly four fields: a finite real scalar ``merit``, a ``kind``
(``lower_bound``, ``optimum``, ``best_known`` or ``target``; every kind is a
claim over the *feasible* points of the problem), a non-empty ``source`` and a
``mapping`` token of a closed registry (only ``feasible_objective/1``). It
holds no point, no callable and no constraint violation.

What these tests pin:

- the four-field shape, and the rejection of the superseded fields (``fun``,
  ``maxcv``, ``point``), of unknown mappings and of non-finite scalars;
- that a record this version does not accept deserializes as unknown or as an
  error and is never reinterpreted;
- the propagation table of the features, single and composed, including every
  subset of the ``custom`` option names, uniformly over the four kinds;
- that a reference changes nothing else: callback counts, histories and whole
  benchmarks are identical with and without it;
- what the scalar means: the default merit function returns the objective at
  every feasible point whatever the initial violation is (the *feasible
  identity*), and on a constrained problem the merit of a run may legitimately
  be below the reference, which is never used as a clamp.

The reference fact must not be confused with two things the code base also
calls "reference": the *reference channel* of a featured problem (the scoring
truth of a run) and the *plain reference* of the profiles (a run-history
minimum). The fact is author/provider metadata and is neither.
"""

import copy
import copyreg
import inspect
import itertools
import logging
import pickle
import warnings
from pathlib import Path

import numpy as np
import pytest

import optiprofiler
import optiprofiler.opclasses as opclasses
from optiprofiler import Feature, FeaturedProblem, Problem, ProblemReference, benchmark
from optiprofiler.composition import AffineView, ComposedFeaturedProblem
from optiprofiler.feature_definitions import LOCAL_OPTIONS, REFERENCE_SAFE_CUSTOM_OPTIONS, STAGE_NAMES, retains_reference
from optiprofiler.loader import load_results_from_h5
from optiprofiler.problem_libraries import load_problem_library, resolve_problem_library
from optiprofiler.profile_utils import _default_merit, compute_merit_values


MAPPING = 'feasible_objective/1'
KINDS = ('lower_bound', 'optimum', 'best_known', 'target')
FIELDS = ('merit', 'kind', 'source', 'mapping')
CENTER = np.array([1.0, 2.0])


def record(merit=0.0, kind='optimum', source='author', mapping=MAPPING):
    return {'merit': merit, 'kind': kind, 'source': source, 'mapping': mapping}


def shifted_sphere(x):
    return float(np.sum((np.asarray(x, dtype=float) - CENTER) ** 2))


def cub_first_coordinate(x):
    return np.array([x[0] - 3.0])


def constrained_problem(reference=None, fun=shifted_sphere, cub=cub_first_coordinate):
    """Bounds, one linear and one nonlinear constraint, so that every feature has something to act on."""
    return Problem(fun, [0.0, 0.0], xl=[-5.0, -5.0], xu=[5.0, 5.0], aub=[[1.0, 1.0]], bub=[10.0], cub=cub,
                   reference=reference)


# One valid callback per ``custom`` option. They are deliberately harmless (the
# value callbacks return what the predecessor returns, the structure callbacks
# return the predecessor's structure): the propagation rule reads the option
# NAMES only, because the framework cannot know what user code does.
def mod_x0(rng, problem):
    return problem.x0 + 0.25


def mod_affine(rng, problem):
    A = np.array([[2.0, 0.0], [1.0, 1.0]])
    return A, np.array([0.5, -0.5]), np.linalg.inv(A)


def mod_bounds(rng, problem):
    return problem.xl, problem.xu


def mod_linear_ub(rng, problem):
    return problem.aub, problem.bub


def mod_linear_eq(rng, problem):
    return problem.aeq, problem.beq


def mod_fun(x, rng, problem):
    return problem.fun(x)


def mod_cub(x, rng, problem):
    return problem.cub(x)


def mod_ceq(x, rng, problem):
    return problem.ceq(x)


CUSTOM_CALLBACKS = {'mod_x0': mod_x0, 'mod_affine': mod_affine, 'mod_bounds': mod_bounds,
                    'mod_linear_ub': mod_linear_ub, 'mod_linear_eq': mod_linear_eq, 'mod_fun': mod_fun,
                    'mod_cub': mod_cub, 'mod_ceq': mod_ceq}
# Written out here, independently of the implementation's constant.
SAFE_CUSTOM_KEYS = {'mod_x0', 'mod_affine'}


def custom_stage(keys):
    return {'name': 'custom', 'options': {key: CUSTOM_CALLBACKS[key] for key in keys}}


# The stage kinds of the propagation table with the expected outcome written as
# a literal (True: the reference is retained). Non-default options are used on
# purpose: for a retaining kind no option value matters.
STAGES = {
    'noisy': ({'name': 'noisy', 'options': {'noise_level': 10.0}}, True),
    'noisy_deterministic': ({'name': 'noisy', 'options': {'noise_mode': 'deterministic'}}, True),
    'truncated': ({'name': 'truncated', 'options': {'significant_digits': 1}}, True),
    'random_nan': ({'name': 'random_nan', 'options': {'nan_rate': 0.9}}, True),
    'nonquantifiable_constraints': ('nonquantifiable_constraints', True),
    'unrelaxable_constraints': ({'name': 'unrelaxable_constraints',
                                 'options': {'unrelaxable_bounds': True, 'unrelaxable_linear_constraints': True,
                                             'unrelaxable_nonlinear_constraints': True}}, True),
    'perturbed_x0': ({'name': 'perturbed_x0', 'options': {'perturbation_level': 10.0}}, True),
    'permuted': ('permuted', True),
    'linearly_transformed': ({'name': 'linearly_transformed', 'options': {'rotated': True, 'condition_factor': 5}}, True),
    'linearly_transformed_unrotated': ({'name': 'linearly_transformed', 'options': {'rotated': False,
                                                                                   'condition_factor': 5}}, True),
    'quantized_default': ('quantized', False),  # ground_truth defaults to True
    'quantized_ground_truth': ({'name': 'quantized', 'options': {'ground_truth': True}}, False),
    'quantized_observed': ({'name': 'quantized', 'options': {'ground_truth': False}}, True),
    'custom_empty': (custom_stage([]), True),
    'custom_x0': (custom_stage(['mod_x0']), True),
    'custom_affine': (custom_stage(['mod_affine']), True),
    'custom_x0_affine': (custom_stage(['mod_x0', 'mod_affine']), True),
    'custom_fun': (custom_stage(['mod_fun']), False),
    'custom_cub': (custom_stage(['mod_cub']), False),
    'custom_ceq': (custom_stage(['mod_ceq']), False),
    'custom_bounds': (custom_stage(['mod_bounds']), False),
    'custom_linear_ub': (custom_stage(['mod_linear_ub']), False),
    'custom_linear_eq': (custom_stage(['mod_linear_eq']), False),
    'custom_affine_fun': (custom_stage(['mod_affine', 'mod_fun']), False),
}


def solver_coordinates(featured, x):
    """
    Test-local: the solver coordinates of the original point ``x``. Every
    coordinate-changing stage maps ``x_predecessor = A @ y + b``; the maps are
    inverted here, stage by stage, without any help from the reference record
    (which holds no point).
    """
    x = np.asarray(x, dtype=float)
    if isinstance(featured, ComposedFeaturedProblem):
        for view in featured._views:
            if isinstance(view, AffineView):
                x = view._inv @ (x - view._b)
        return x
    _, b, inv = featured._runtime.modifier_affine(featured._seed, featured._problem)
    return inv @ (x - b)


class TestRecordShape:

    def test_record_has_exactly_four_fields(self):
        reference = ProblemReference(2.5, 'best_known', 'DOI 10.0/example', MAPPING)
        assert ProblemReference.FIELDS == FIELDS
        assert ProblemReference.__slots__ == ('_merit', '_kind', '_source', '_mapping')
        assert not hasattr(reference, '__dict__')  # nothing can be attached to a record
        assert list(reference.as_dict()) == list(FIELDS)
        assert reference.as_dict() == {'merit': 2.5, 'kind': 'best_known', 'source': 'DOI 10.0/example',
                                       'mapping': MAPPING}
        assert (reference.merit, reference.kind, reference.source, reference.mapping) == \
            (2.5, 'best_known', 'DOI 10.0/example', MAPPING)
        assert type(reference.merit) is float
        assert ProblemReference.KINDS == KINDS
        assert 'ProblemReference' in optiprofiler.__all__
        assert repr(reference) == ("ProblemReference(merit=2.5, kind='best_known', source='DOI 10.0/example', "
                                   "mapping='feasible_objective/1')")

    def test_superseded_record_and_transport_api_is_gone(self):
        # The superseded layout stored an objective value, a violation and a
        # point, and transported the point through an inverse affine map. None
        # of this exists: the record has no such attribute, the constructor has
        # no such parameter and nothing transports anything.
        for name in ('point', 'with_point', 'fun', 'maxcv', 'affine_inverse'):
            assert not hasattr(ProblemReference, name), name
        assert list(inspect.signature(ProblemReference.__init__).parameters) == ['self', 'merit', 'kind', 'source',
                                                                                 'mapping']
        for name in ('_propagate_reference', '_reference_scalar'):
            assert not hasattr(opclasses, name), name
        assert 'with_point' not in Path(optiprofiler.composition.__file__).read_text(encoding='utf-8')
        with pytest.raises(TypeError):
            ProblemReference(0.0, 'optimum', 'author', MAPPING, point=[1.0, 2.0])
        with pytest.raises(TypeError):
            ProblemReference(fun=0.0, kind='optimum', source='author')

    def test_from_mapping_and_from_instance(self):
        problem = Problem(shifted_sphere, [0.0, 0.0], reference=record())
        assert isinstance(problem.reference, ProblemReference)
        same = Problem(shifted_sphere, [0.0, 0.0], reference=ProblemReference(0.0, 'optimum', 'author', MAPPING))
        assert same.reference == problem.reference and hash(same.reference) == hash(problem.reference)
        assert problem.reference != ProblemReference(0.0, 'lower_bound', 'author', MAPPING)
        assert problem.reference != ProblemReference(0.0, 'optimum', 'someone else', MAPPING)
        assert problem.reference != ProblemReference(1.0, 'optimum', 'author', MAPPING)
        assert problem.reference != record()  # a dict is not a record
        # The field order of the input does not matter; the stored order is fixed.
        shuffled = {'mapping': MAPPING, 'source': 'author', 'kind': 'optimum', 'merit': 0.0}
        assert list(ProblemReference.from_record(shuffled).as_dict()) == list(FIELDS)

    @pytest.mark.parametrize('kind', KINDS)
    def test_every_kind_is_accepted_with_the_same_shape(self, kind):
        reference = Problem(shifted_sphere, [0.0, 0.0], reference=record(kind=kind)).reference
        assert reference.kind == kind and list(reference.as_dict()) == list(FIELDS)

    def test_record_is_immutable(self):
        reference = ProblemReference(1.0, 'target', 'campaign', MAPPING)
        for name in FIELDS + ('_merit', 'point'):
            with pytest.raises(AttributeError):
                setattr(reference, name, 2.0)
        with pytest.raises(AttributeError):
            del reference.kind
        problem = Problem(shifted_sphere, [0.0, 0.0], reference=reference)
        with pytest.raises(AttributeError):
            problem.reference = None  # set at construction only
        assert copy.copy(reference) == reference and copy.deepcopy(reference) == reference

    def test_omitted_reference_is_unknown(self):
        problem = Problem(shifted_sphere, [0.0, 0.0])
        assert problem.reference is None
        assert Problem(shifted_sphere, [0.0, 0.0], reference=None).reference is None
        assert FeaturedProblem(problem, Feature('noisy'), 5, 0).reference is None
        assert FeaturedProblem(problem, Feature('noisy+permuted'), 5, 0).reference is None

    @pytest.mark.parametrize('legacy', [
        {'fun': 0.0, 'kind': 'optimum', 'source': 'author'},
        {'fun': 0.0, 'maxcv': 0.0, 'kind': 'optimum', 'source': 'author', 'point': [1.0, 2.0]},
        {'fun': 6.25, 'maxcv': 0.5, 'kind': 'best_known', 'source': 'test run', 'point': None},
        dict(record(), fun=0.0),
        dict(record(), maxcv=0.0),
        dict(record(), point=[1.0, 2.0]),
        dict(record(), point=None),
    ])
    def test_superseded_fields_are_rejected_not_reinterpreted(self, legacy):
        # A record of the superseded layout carries an objective value in
        # ``fun``. Reading it as ``merit`` would be a guess about another
        # contract, so the whole record is refused, even when the four current
        # fields are present as well.
        with pytest.raises(ValueError, match='superseded record layout'):
            Problem(shifted_sphere, [0.0, 0.0], reference=legacy)
        with pytest.raises(ValueError, match='never reinterpreted'):
            ProblemReference.from_record(legacy)

    @pytest.mark.parametrize('missing', FIELDS)
    def test_every_field_is_required(self, missing):
        # No field has a default. In particular an omitted mapping is not
        # read as 'feasible_objective/1': the author states how the scalar is read.
        incomplete = {key: value for key, value in record().items() if key != missing}
        with pytest.raises(ValueError, match='missing the required field'):
            Problem(shifted_sphere, [0.0, 0.0], reference=incomplete)
        arguments = [value for key, value in record().items() if key != missing]
        with pytest.raises(TypeError):
            ProblemReference(*arguments)

    @pytest.mark.parametrize('extra', ['x', 'optimum', 'value', 'merit_fun', 'Merit', 3])
    def test_unknown_fields_are_rejected(self, extra):
        with pytest.raises(ValueError, match='Unknown problem reference field'):
            ProblemReference.from_record(dict(record(), **{extra: 1}) if isinstance(extra, str)
                                         else {**record(), extra: 1})

    @pytest.mark.parametrize('naked', [1.0, 0, np.float64(2.0), np.array(1.0), [0.0, 'optimum', 'author', MAPPING],
                                       (0.0, 'optimum', 'author', MAPPING), 'optimum'])
    def test_naked_scalars_and_sequences_are_rejected(self, naked):
        # The roadmap forbids a public scalar property: a number alone has no
        # kind, no provenance and no stated reading.
        with pytest.raises(TypeError, match='naked scalar'):
            Problem(shifted_sphere, [0.0, 0.0], reference=naked)

    @pytest.mark.parametrize('kind, error', [('bogus', ValueError), ('Optimum', ValueError), ('optimum ', ValueError),
                                             ('', ValueError), (3, TypeError), (None, TypeError), (b'optimum', TypeError)])
    def test_kind_is_one_of_the_four(self, kind, error):
        with pytest.raises(error):
            ProblemReference(0.0, kind, 'author', MAPPING)

    @pytest.mark.parametrize('source, error', [('', ValueError), ('   ', ValueError), ('\t\n', ValueError),
                                               (None, TypeError), (5, TypeError), (b'author', TypeError)])
    def test_source_is_non_empty_text(self, source, error):
        with pytest.raises(error):
            ProblemReference(0.0, 'optimum', source, MAPPING)


class TestMappingRegistry:

    def test_registry_is_closed(self):
        assert ProblemReference.MAPPINGS == ('feasible_objective/1',)
        assert isinstance(ProblemReference.MAPPINGS, tuple)  # not a container that code can append to
        # There is no registration entry point, anywhere in the public package.
        for owner in (ProblemReference, opclasses, optiprofiler):
            assert not [name for name in dir(owner) if 'register' in name.lower() and 'mapping' in name.lower()]

    def test_registry_cannot_be_extended(self, monkeypatch):
        # A subclass could redefine the registry, a property or the
        # validation: a user-defined mapping by another name. It is refused.
        with pytest.raises(TypeError, match='cannot be subclassed'):
            class Extended(ProblemReference):  # noqa: F841
                MAPPINGS = ProblemReference.MAPPINGS + ('user_mapping/1',)
        # The class attributes are informational. Validation reads private
        # literals, so assigning to the attributes adds neither a mapping nor a kind.
        monkeypatch.setattr(ProblemReference, 'MAPPINGS', ProblemReference.MAPPINGS + ('user_mapping/1',))
        monkeypatch.setattr(ProblemReference, 'KINDS', ProblemReference.KINDS + ('user_kind',))
        with pytest.raises(ValueError, match='closed registry'):
            ProblemReference(0.0, 'optimum', 'author', 'user_mapping/1')
        with pytest.raises(ValueError, match='must be one of'):
            ProblemReference(0.0, 'user_kind', 'author', MAPPING)
        with pytest.warns(RuntimeWarning, match='read as unknown'):
            assert opclasses._restore_problem_reference(record(mapping='user_mapping/1')) is None
        assert ProblemReference(0.0, 'optimum', 'author', MAPPING).mapping == MAPPING

    def test_text_subclasses_cannot_pass_the_registry_or_change_what_is_stored(self):
        class Chameleon(str):
            """Text that claims to equal everything and prints as something else."""

            def __eq__(self, other):
                return True

            def __hash__(self):
                return hash('feasible_objective/1')

            def __str__(self):
                return 'user_mapping/1'

        # Judged by its content, 'user_mapping/1' is unknown, whatever it claims to equal.
        with pytest.raises(ValueError, match='closed registry'):
            ProblemReference(0.0, 'optimum', 'author', Chameleon('user_mapping/1'))
        with pytest.raises(ValueError, match='must be one of'):
            ProblemReference(0.0, Chameleon('user_kind'), 'author', MAPPING)
        # A subclass with valid content is accepted and stored as exact text.
        reference = ProblemReference(0.0, Chameleon('optimum'), Chameleon('author'), Chameleon(MAPPING))
        assert [type(value) for value in (reference.kind, reference.source, reference.mapping)] == [str, str, str]
        assert reference.as_dict() == record()
        assert reference == ProblemReference.from_record(record())
        numpy_text = ProblemReference(0.0, np.str_('optimum'), np.str_('author'), np.str_(MAPPING))
        assert type(numpy_text.kind) is str and numpy_text == reference

    @pytest.mark.parametrize('mapping', ['feasible_objective/2', 'feasible_objective/0', 'feasible_objective',
                                         'feasible_objective/1 ', ' feasible_objective/1', 'Feasible_Objective/1',
                                         'FEASIBLE_OBJECTIVE/1', 'feasible_objective/1.0', 'feasible_objective/01',
                                         'objective/1', 'merit/1', 'default_merit/1', 'identity', ''])
    def test_unknown_tokens_are_rejected(self, mapping):
        # A token is matched exactly. Near misses (another version, another
        # case, surrounding blanks) are different tokens, not spellings of the
        # known one: accepting them would be a guess.
        with pytest.raises(ValueError, match='closed registry'):
            ProblemReference(0.0, 'optimum', 'author', mapping)
        with pytest.raises(ValueError, match='never guessed'):
            Problem(shifted_sphere, [0.0, 0.0], reference=record(mapping=mapping))

    @pytest.mark.parametrize('mapping', [lambda f, v, v0: f, _default_merit, shifted_sphere, float, ProblemReference])
    def test_callables_are_never_mappings(self, mapping):
        # No arbitrary callbacks and no user mappings: even the default merit
        # function itself is refused as a mapping.
        with pytest.raises(TypeError, match='callables and user-defined mappings are not accepted'):
            ProblemReference(0.0, 'optimum', 'author', mapping)

    @pytest.mark.parametrize('mapping', [None, 1, 1.0, b'feasible_objective/1', ['feasible_objective/1'],
                                         ('feasible_objective', 1), {'name': 'feasible_objective', 'version': 1}])
    def test_non_text_tokens_are_rejected(self, mapping):
        with pytest.raises(TypeError, match='string token'):
            ProblemReference(0.0, 'optimum', 'author', mapping)


class TestFiniteValidation:

    @pytest.mark.parametrize('merit', [np.nan, np.inf, -np.inf, float('nan'), float('inf'), np.float64('nan'),
                                       np.float32('inf'), np.array(np.nan), np.float16('inf')])
    def test_non_finite_merits_are_rejected(self, merit):
        with pytest.raises(ValueError, match='must be finite'):
            ProblemReference(merit, 'lower_bound', 'author', MAPPING)
        with pytest.raises(ValueError, match='must be finite'):
            Problem(shifted_sphere, [0.0, 0.0], reference=record(merit=merit, kind='lower_bound'))

    @pytest.mark.parametrize('merit', ['1', b'1', True, False, np.True_, None, [1.0], (1.0,), np.array([1.0]),
                                       np.array([[1.0]]), np.array([1.0, 2.0]), 1.0 + 0.0j, np.complex128(1.0),
                                       {'value': 1.0}, 10 ** 400, object()])
    def test_non_scalar_and_non_real_merits_are_rejected(self, merit):
        # A one-element array is not a scalar, a boolean is not a magnitude,
        # and a complex number with a zero imaginary part is still complex.
        with pytest.raises(TypeError, match='must be a real scalar'):
            ProblemReference(merit, 'optimum', 'author', MAPPING)

    @pytest.mark.parametrize('merit', [2 ** 53 + 1, -(2 ** 53) - 1, 2 ** 63 - 1, np.int64(2 ** 53 + 1),
                                       np.uint64(2 ** 64 - 1)])
    def test_integers_that_a_float_cannot_hold_are_rejected(self, merit):
        # Storing 2**53 + 1 as 2.0**53 would silently change the claim.
        with pytest.raises(ValueError, match='represented exactly'):
            ProblemReference(merit, 'optimum', 'author', MAPPING)

    @pytest.mark.parametrize('merit, stored', [(0, 0.0), (-3, -3.0), (2 ** 53, 2.0 ** 53), (2 ** 60, 2.0 ** 60),
                                               (-0.0, 0.0), (1e308, 1e308), (-1e308, -1e308), (5e-324, 5e-324),
                                               (np.float32(0.5), 0.5), (np.float16(-2.0), -2.0), (np.int8(-7), -7.0),
                                               (np.array(2.25), 2.25), (np.float64(1e-300), 1e-300)])
    def test_finite_real_scalars_are_stored_exactly_as_float(self, merit, stored):
        reference = ProblemReference(merit, 'target', 'author', MAPPING)
        assert type(reference.merit) is float and reference.merit == stored

    def test_validation_never_evaluates_the_problem(self):
        calls = []

        def counted(x):
            calls.append('fun')
            return shifted_sphere(x)

        def counted_cub(x):
            calls.append('cub')
            return cub_first_coordinate(x)

        without = constrained_problem(fun=counted, cub=counted_cub)
        baseline = list(calls)
        with_reference = constrained_problem(reference=record(), fun=counted, cub=counted_cub)
        # Construction makes exactly the calls it makes without a reference,
        # and none of them is an objective evaluation.
        assert calls == baseline + baseline and 'fun' not in calls
        assert with_reference.reference.kind == 'optimum' and without.reference is None
        for _ in range(50):
            with_reference.reference
        assert calls == baseline + baseline


    def test_malformed_record_is_rejected_before_any_callback_is_touched(self):
        # The Python constructor probes the constraint callbacks at x0 for
        # their dimensions. The reference is validated before that, so a
        # malformed record never causes a callback of the problem to run.
        calls = []

        def counted(x):
            calls.append('fun')
            return shifted_sphere(x)

        def counted_cub(x):
            calls.append('cub')
            return cub_first_coordinate(x)

        malformed = [record(merit=np.nan), record(mapping='feasible_objective/2'), record(mapping=counted),
                     record(kind='bogus'), record(source=''), dict(record(), point=[1.0, 2.0]),
                     {'fun': 0.0, 'kind': 'optimum', 'source': 'author'}, 0.0]
        for reference in malformed:
            with pytest.raises((TypeError, ValueError)):
                constrained_problem(reference=reference, fun=counted, cub=counted_cub)
        assert calls == []
        # Control: a valid record lets the constructor go on to its usual probe.
        constrained_problem(reference=record(), fun=counted, cub=counted_cub)
        assert 'cub' in calls and 'fun' not in calls


class TestPropagationRule:
    """The pure rule ``retains_reference(name, options)``; nothing is built or evaluated."""

    @pytest.mark.parametrize('label', sorted(STAGES))
    def test_rule_matches_the_table(self, label):
        stage, expected = STAGES[label]
        name, options = (stage, {}) if isinstance(stage, str) else (stage['name'], stage['options'])
        if name == 'quantized':
            options = dict({'ground_truth': True}, **options)  # the rule sees the defaulted options
        assert retains_reference(name, options) is expected

    def test_every_stage_kind_has_an_entry(self):
        # A stage kind added later must be classified here before it can
        # retain anything: the rule fails closed for a name it does not know.
        assert {(stage if isinstance(stage, str) else stage['name']) for stage, _ in STAGES.values()} == \
            set(STAGE_NAMES) - {'plain'}
        assert retains_reference('plain', {}) is True
        for unknown in ('', 'Noisy', 'rotated', 'mod_affine', 'future_feature', None, 3):
            assert retains_reference(unknown, {}) is False

    def test_custom_whitelist_is_exactly_x0_and_affine(self):
        assert REFERENCE_SAFE_CUSTOM_OPTIONS == frozenset(SAFE_CUSTOM_KEYS)
        assert set(LOCAL_OPTIONS['custom']) == set(CUSTOM_CALLBACKS)  # the eight custom options, all covered below
        # A whitelist, not a blacklist: an option name nobody has defined yet
        # makes the reference unknown.
        assert retains_reference('custom', {'mod_x0': mod_x0, 'mod_future': mod_fun}) is False
        assert retains_reference('custom', {'n_runs': 3}) is False

    @pytest.mark.parametrize('options', [{}, {'ground_truth': True}, {'ground_truth': 1}, {'ground_truth': 0},
                                         {'ground_truth': None}, {'ground_truth': 'false'}, {'ground_truth': np.False_}])
    def test_quantized_is_retained_only_for_an_explicit_false(self, options):
        # Missing or non-boolean values count as ground truth (fail closed).
        assert retains_reference('quantized', options) is False
        assert retains_reference('quantized', {'ground_truth': False}) is True

    def test_rule_reads_names_only(self):
        def explode(*args):
            raise AssertionError('a custom callback was called to decide the propagation')

        assert retains_reference('custom', {'mod_x0': explode, 'mod_affine': explode}) is True
        assert retains_reference('custom', {'mod_fun': explode}) is False


class TestPropagationTable:

    @pytest.mark.parametrize('kind', KINDS)
    @pytest.mark.parametrize('label', ['plain'] + sorted(STAGES))
    def test_single_stage(self, label, kind):
        # The rule is uniform over the kinds: each is a claim over feasible
        # points, so each survives exactly the same stages.
        problem = constrained_problem(record(merit=0.5, kind=kind, source=f'table:{kind}'))
        stage, expected = ('plain', True) if label == 'plain' else STAGES[label]
        featured = FeaturedProblem(problem, Feature([stage]), 10, 3)
        assert not isinstance(featured, ComposedFeaturedProblem)
        if expected:
            # Retained means unchanged: same four fields, nothing transported.
            assert featured.reference == problem.reference
            assert featured.reference.as_dict() == record(merit=0.5, kind=kind, source=f'table:{kind}')
        else:
            assert featured.reference is None
        assert problem.reference == ProblemReference(0.5, kind, f'table:{kind}', MAPPING)  # the original is untouched

    @pytest.mark.parametrize('first, second', list(itertools.product(sorted(STAGES), repeat=2)),
                             ids=lambda value: value)
    def test_every_ordered_pair_is_safe_only_if_both_stages_are(self, first, second):
        problem = constrained_problem(record(kind='lower_bound', merit=-1.0))
        featured = FeaturedProblem(problem, Feature([STAGES[first][0], STAGES[second][0]]), 10, 3)
        assert isinstance(featured, ComposedFeaturedProblem)
        expected = STAGES[first][1] and STAGES[second][1]
        assert (featured.reference == problem.reference) if expected else (featured.reference is None)

    @pytest.mark.parametrize('kind', KINDS)
    @pytest.mark.parametrize('labels, expected', [
        (('perturbed_x0', 'noisy', 'permuted', 'linearly_transformed', 'truncated'), True),
        (('custom_x0_affine', 'quantized_observed', 'unrelaxable_constraints', 'random_nan'), True),
        (('quantized_ground_truth', 'noisy', 'permuted'), False),   # unsafe first
        (('noisy', 'custom_bounds', 'permuted'), False),            # unsafe in the middle
        (('noisy', 'permuted', 'custom_fun'), False),               # unsafe last
        # A later safe stage never restores what an earlier stage made unknown.
        (('quantized_default', 'quantized_observed', 'custom_affine'), False),
    ])
    def test_longer_compositions(self, labels, expected, kind):
        problem = constrained_problem(record(kind=kind))
        featured = FeaturedProblem(problem, Feature([STAGES[label][0] for label in labels]), 10, 3)
        assert isinstance(featured, ComposedFeaturedProblem)
        assert (featured.reference == problem.reference) if expected else (featured.reference is None)
        # Every view reports the fact of its own prefix of the pipeline.
        prefix = True
        for label, view in zip(labels, featured._views):
            prefix = prefix and STAGES[label][1]
            assert (view.reference == problem.reference) if prefix else (view.reference is None)

    def test_plain_tokens_do_not_change_the_outcome(self):
        problem = constrained_problem(record())
        assert FeaturedProblem(problem, Feature('plain+noisy+plain'), 10, 3).reference == problem.reference
        assert FeaturedProblem(problem, Feature('plain+quantized'), 10, 3).reference is None

    def test_shorthand_options_reach_the_rule(self):
        problem = constrained_problem(record())
        assert FeaturedProblem(problem, Feature('quantized', ground_truth=False), 10, 3).reference == problem.reference
        assert FeaturedProblem(problem, Feature('quantized+noisy', ground_truth=False), 10, 3).reference == \
            problem.reference
        assert FeaturedProblem(problem, Feature('custom', mod_affine=mod_affine), 10, 3).reference == problem.reference
        assert FeaturedProblem(problem, Feature('custom', mod_fun=mod_fun), 10, 3).reference is None

    def test_composed_recorder_pickles_with_its_reference(self):
        problem = constrained_problem(record())
        composed = FeaturedProblem(problem, Feature('noisy+permuted'), 10, 3)
        assert pickle.loads(pickle.dumps(composed)).reference == problem.reference


class TestCustomKeySubsets:
    """All 2**8 subsets of the custom option names: retained iff the subset is within {mod_x0, mod_affine}."""

    SUBSETS = [keys for size in range(len(CUSTOM_CALLBACKS) + 1)
               for keys in itertools.combinations(sorted(CUSTOM_CALLBACKS), size)]

    @staticmethod
    def refused(keys):
        # ``mod_affine`` here is dense, so the bounds of the problem are posed
        # as linear rows, and ``mod_linear_ub`` replaces the linear rows: unless
        # ``mod_bounds`` takes over the bounds, they would leave the posed
        # problem silently. The affine safeguard refuses to build that.
        return {'mod_affine', 'mod_linear_ub'} <= set(keys) and 'mod_bounds' not in keys

    def check(self, build):
        assert len(self.SUBSETS) == 256
        problem = constrained_problem(record(kind='best_known', merit=0.5))
        mismatches, retained = [], []
        for keys in self.SUBSETS:
            expected = set(keys) <= SAFE_CUSTOM_KEYS
            if self.refused(keys):
                # No featured problem, so nothing that could carry a claim; none was expected either.
                assert not expected
                with pytest.raises(ValueError, match='would be dropped silently'):
                    build(problem, keys)
                continue
            reference = build(problem, keys).reference
            if (reference == problem.reference) != expected or (reference is None) == expected:
                mismatches.append(keys)
            if reference is not None:
                retained.append(set(keys))
        assert not mismatches
        assert sum(map(self.refused, self.SUBSETS)) == 32
        # Exactly the four subsets of {mod_x0, mod_affine} retain the record.
        assert sorted(retained, key=sorted) == sorted([set(), {'mod_x0'}, {'mod_affine'}, {'mod_x0', 'mod_affine'}],
                                                      key=sorted)

    def test_single_custom_stage(self):
        self.check(lambda problem, keys: FeaturedProblem(problem, Feature([custom_stage(keys)]), 10, 3))

    def test_custom_stage_first_in_a_composition(self):
        self.check(lambda problem, keys: FeaturedProblem(problem, Feature([custom_stage(keys), 'noisy']), 10, 3))

    def test_custom_stage_last_in_a_composition(self):
        self.check(lambda problem, keys: FeaturedProblem(problem, Feature(['permuted', custom_stage(keys)]), 10, 3))

    def test_pure_rule(self):
        for keys in self.SUBSETS:
            assert retains_reference('custom', {key: CUSTOM_CALLBACKS[key] for key in keys}) is \
                (set(keys) <= SAFE_CUSTOM_KEYS), keys


class CallCounter:
    """Counts every call of the problem callbacks and of the custom callbacks, by name."""

    def __init__(self):
        self.counts = {}

    def wrap(self, name, function):
        def counted(*args):
            self.counts[name] = self.counts.get(name, 0) + 1
            return function(*args)
        return counted

    def problem(self, reference):
        return constrained_problem(reference, fun=self.wrap('fun', shifted_sphere),
                                   cub=self.wrap('cub', cub_first_coordinate))

    def stage(self, keys):
        return {'name': 'custom', 'options': {key: self.wrap(key, CUSTOM_CALLBACKS[key]) for key in keys}}


class TestCallbackCountInvariance:

    FEATURES = ['plain', 'noisy', 'permuted', 'linearly_transformed', 'perturbed_x0', 'quantized',
                'unrelaxable_constraints', 'noisy+permuted', 'linearly_transformed+quantized+noisy',
                ('custom', ('mod_x0', 'mod_affine')), ('custom', tuple(sorted(CUSTOM_CALLBACKS))),
                ('custom+noisy', ('mod_affine',)), ('custom+noisy', ('mod_fun', 'mod_cub'))]

    @staticmethod
    def exercise(counter, reference, feature):
        problem = counter.problem(reference)
        if isinstance(feature, tuple):
            name, keys = feature
            stages = [counter.stage(keys)] + name.split('+')[1:]
        else:
            stages = feature.split('+')
        featured = FeaturedProblem(problem, Feature(stages), 20, 11)
        observed = []
        for shift in (0.0, 0.25, -0.5):
            x = featured.x0 + shift
            observed.append((featured.fun(x), tuple(featured.cub(x)), featured.maxcv(x)))
        return problem, featured, observed

    @pytest.mark.parametrize('feature', FEATURES, ids=str)
    def test_a_reference_costs_no_callback_call(self, feature):
        # The superseded design called the affine modifier once more to
        # transport a point. Nothing is transported now, so building and using
        # a featured problem makes exactly the same calls with a reference as
        # without one, user callbacks included.
        with_counter, without_counter = CallCounter(), CallCounter()
        problem, featured, observed = self.exercise(with_counter, record(), feature)
        _, plain_featured, plain_observed = self.exercise(without_counter, None, feature)
        assert with_counter.counts == without_counter.counts
        assert with_counter.counts.get('fun', 0) > 0  # the comparison is not vacuous
        np.testing.assert_array_equal(featured.x0, plain_featured.x0)
        assert repr(observed) == repr(plain_observed)  # NaN-safe comparison of the observed values
        np.testing.assert_array_equal(featured.fun_hist, plain_featured.fun_hist)
        np.testing.assert_array_equal(featured.maxcv_hist, plain_featured.maxcv_hist)
        # Reading the fact, from the problem or from the featured problem, calls nothing.
        before = dict(with_counter.counts)
        for _ in range(100):
            problem.reference
            featured.reference
        assert with_counter.counts == before
        assert plain_featured.reference is None


class TestSafeCoordinateChanges:
    """
    A permutation and a valid affine change of variables keep the reference
    although the record holds no point. The optimum is known to the test only:

        minimize (x1 - 1)^2 + (x2 - 2)^2 + 3  subject to  x1 + x2 <= 2.5,  -5 <= x <= 5.

    The unconstrained minimizer (1, 2) violates the linear constraint, so the
    solution is its projection onto the half-plane, x* = (0.75, 1.75), with
    value 2 * 0.25^2 + 3 = 3.125.
    """

    OPTIMUM = 3.125
    MINIMIZER = np.array([0.75, 1.75])

    @staticmethod
    def problem():
        def objective(x):
            return shifted_sphere(x) + 3.0

        return Problem(objective, [0.0, 0.0], xl=[-5.0, -5.0], xu=[5.0, 5.0], aub=[[1.0, 1.0]], bub=[2.5],
                       reference=record(merit=TestSafeCoordinateChanges.OPTIMUM, source='test-local derivation'))

    @pytest.mark.parametrize('stages', [
        ['permuted'],
        [{'name': 'linearly_transformed', 'options': {'rotated': True, 'condition_factor': 4}}],
        [{'name': 'linearly_transformed', 'options': {'rotated': False, 'condition_factor': 4}}],
        [custom_stage(['mod_affine'])],
        [custom_stage(['mod_x0', 'mod_affine'])],
        ['permuted', 'noisy'],
        ['linearly_transformed', 'permuted'],
        [custom_stage(['mod_affine']), 'linearly_transformed', 'truncated'],
        ['perturbed_x0', 'permuted', custom_stage(['mod_affine']), {'name': 'quantized', 'options': {'ground_truth': False}}],
    ], ids=lambda stages: '+'.join(stage if isinstance(stage, str) else stage['name'] for stage in stages))
    @pytest.mark.parametrize('seed', [0, 3, 17])
    def test_retained_claim_is_true_for_the_featured_problem(self, stages, seed):
        problem = self.problem()
        featured = FeaturedProblem(problem, Feature(stages), 50, seed)
        assert featured.reference == problem.reference
        assert featured.reference.as_dict() == record(merit=self.OPTIMUM, source='test-local derivation')

        # 1. The claimed value is attained at a feasible point of the featured
        #    problem: the image of x*, found here without any stored point.
        minimizer = solver_coordinates(featured, self.MINIMIZER)
        value, violation = featured._evaluate_truth(minimizer)
        assert value == pytest.approx(self.OPTIMUM, abs=1e-9) and violation == pytest.approx(0.0, abs=1e-9)

        # 2. The change of variables is a bijection of the feasible set that
        #    preserves the truth: at the image of any x, the featured truth is
        #    the original (objective, violation). Hence no feasible point of
        #    the featured problem is below the reference either.
        rng = np.random.default_rng(seed)
        feasible = 0
        for x in rng.uniform(-6.0, 6.0, size=(300, 2)):
            value, violation = featured._evaluate_truth(solver_coordinates(featured, x))
            assert value == pytest.approx(problem.fun(x), rel=1e-9, abs=1e-9)
            assert violation == pytest.approx(problem.maxcv(x), rel=1e-9, abs=1e-9)
            if problem.maxcv(x) == 0.0:
                feasible += 1
                assert value >= self.OPTIMUM - 1e-9
        assert feasible > 50  # the sample really covers the feasible set

    def test_coordinates_really_change(self):
        # Guard against a vacuous test: the solver coordinates of x* differ
        # from x*, so a stored point could not have been kept as it was.
        problem = self.problem()
        for stages in (['linearly_transformed'], [custom_stage(['mod_affine'])]):
            featured = FeaturedProblem(problem, Feature(stages), 50, 3)
            assert not np.allclose(solver_coordinates(featured, self.MINIMIZER), self.MINIMIZER)

    def test_invalid_affine_never_yields_a_featured_problem(self):
        # "Valid affine coordinate change" is enforced where the problem is
        # built: a custom map whose inverse is wrong raises, so it never
        # carries a reference anywhere.
        def broken(rng, problem):
            return np.array([[2.0, 0.0], [1.0, 1.0]]), np.zeros(2), np.eye(2)

        with pytest.raises(ValueError, match='not an identity matrix'):
            FeaturedProblem(self.problem(), Feature('custom', mod_affine=broken), 10, 3)
        with pytest.raises(ValueError, match='not an identity matrix'):
            FeaturedProblem(self.problem(), Feature([{'name': 'custom', 'options': {'mod_affine': broken}}, 'noisy']), 10, 3)


class TestDefaultMeritFeasibleIdentity:
    """
    ``feasible_objective/1`` reads the scalar as an objective value over
    feasible points. It is comparable with run merits exactly when the merit
    function maps a feasible point to its objective value, whatever the
    violation of the run's initial point is.
    """

    INITIAL_VIOLATIONS = [0.0, 0.2, 50.0, np.nan]
    VALUES = [-1e30, -3.5, -1e-300, -0.0, 0.0, 5e-324, 2.5, 1e30, np.inf, -np.inf]

    @pytest.mark.parametrize('maxcv_init', INITIAL_VIOLATIONS, ids=str)
    def test_default_merit_is_the_objective_at_feasible_points(self, maxcv_init):
        for value in self.VALUES:
            merit = _default_merit(value, 0.0, maxcv_init)
            assert merit == value and type(merit) is type(value)  # exactly the value, not an approximation
        # The same holds through the vectorized path used for the histories.
        values = np.array(self.VALUES)
        merits = compute_merit_values(_default_merit, values, np.zeros_like(values), maxcv_init)
        np.testing.assert_array_equal(merits, values)

    @pytest.mark.parametrize('maxcv_init', INITIAL_VIOLATIONS, ids=str)
    @pytest.mark.parametrize('kind', KINDS)
    def test_reference_merit_is_a_fixed_point_of_the_default_merit(self, kind, maxcv_init):
        reference = ProblemReference(-7.25, kind, 'author', MAPPING)
        assert _default_merit(reference.merit, 0.0, maxcv_init) == reference.merit

    @pytest.mark.parametrize('maxcv_init, tolerated', [(0.0, 1e-10), (0.2, 1e-10), (50.0, 5e-9), (np.nan, 1e-10)])
    def test_identity_is_a_statement_about_feasible_points_only(self, maxcv_init, tolerated):
        # Away from feasibility the default merit tolerates (below v1),
        # penalizes (up to v2) or discards (beyond v2) a point. This is the
        # violation handling that lets run merits fall below the reference.
        assert _default_merit(2.0, tolerated, maxcv_init) == 2.0
        penalized = _default_merit(2.0, 0.05, maxcv_init)
        assert penalized == pytest.approx(2.0 + 1e5 * (0.05 - tolerated)) and penalized > 2.0
        assert _default_merit(2.0, 1e3, maxcv_init) == np.inf

    def test_custom_merit_must_preserve_the_identity_before_any_consumer_uses_the_reference(self):
        # Test-local check that a future consumer has to make before it
        # compares run merits with a record: the probes are the four initial
        # violations above. A merit function without the identity puts run
        # merits in another space than the record.
        def preserves_feasible_identity(merit_fun):
            return all(merit_fun(value, 0.0, maxcv_init) == value
                       for value in (-3.5, 0.0, 2.5) for maxcv_init in self.INITIAL_VIOLATIONS)

        def shifted(fun_value, maxcv_value, maxcv_init):
            return fun_value + 1.0

        def scaled_by_initial_violation(fun_value, maxcv_value, maxcv_init):
            return fun_value * (1.0 + maxcv_init)

        def quadratic_penalty(fun_value, maxcv_value, maxcv_init):
            return fun_value + 1e3 * maxcv_value ** 2

        assert preserves_feasible_identity(_default_merit)
        assert preserves_feasible_identity(quadratic_penalty)  # a custom merit may well have the identity
        assert not preserves_feasible_identity(shifted)
        assert not preserves_feasible_identity(scaled_by_initial_violation)
        reference = ProblemReference(0.0, 'optimum', 'author', MAPPING)
        # With the shifted merit a feasible optimal run has merit 1, not 0:
        # comparing it with the record would be meaningless, not "close".
        assert shifted(reference.merit, 0.0, 0.0) != reference.merit


def steep_line(x):
    return 1e7 * float(np.asarray(x, dtype=float)[0])


def nonnegative(x):
    return np.array([-np.asarray(x, dtype=float)[0]])


class TestConstrainedMeritBelowReference:
    """
    minimize 1e7 * x  subject to  -x <= 0: the feasible optimum is 0, at x = 0.
    An infeasible point has a negative objective, and the default merit does
    not push it back above 0: below the tolerance v1 it is the objective, and
    up to v2 the penalty 1e5 * (v - v1) is smaller than the gain 1e7 * v.
    """

    TOLERATED = np.array([-5e-11])   # violation 5e-11 <= v1 = 1e-10: merit = objective = -5e-4
    PENALIZED = np.array([-1e-4])    # violation 1e-4 in (v1, v2]: merit = -1e3 + 1e5 * (1e-4 - 1e-10)

    @staticmethod
    def problem(reference=True):
        return Problem(steep_line, [1.0], cub=nonnegative,
                       reference=record(merit=0.0, kind='optimum', source='analytic') if reference else None)

    def test_run_merits_fall_below_the_reference_and_are_not_clamped(self):
        problem = self.problem()
        featured = FeaturedProblem(problem, Feature('plain'), 10, 0)
        assert featured.maxcv_init == 0.0 and featured.reference == problem.reference
        for x in (np.array([1.0]), np.array([0.0]), self.TOLERATED, self.PENALIZED):
            featured.fun(x)
        np.testing.assert_array_equal(featured.fun_hist, [1e7, 0.0, -5e-4, -1e3])
        np.testing.assert_array_equal(featured.maxcv_hist, [0.0, 0.0, 5e-11, 1e-4])
        merits = compute_merit_values(_default_merit, np.array(featured.fun_hist), np.array(featured.maxcv_hist),
                                      featured.maxcv_init)
        np.testing.assert_allclose(merits, [1e7, 0.0, -5e-4, -1e3 + 1e5 * (1e-4 - 1e-10)], rtol=1e-12)
        reference = featured.reference.merit
        # The feasible points respect the claim; the infeasible ones are below
        # it, in the tolerance zone and in the penalty zone alike.
        assert merits[0] >= reference and merits[1] == reference
        assert merits[2] < reference and merits[3] < reference
        assert merits.min() == pytest.approx(-989.99999, rel=1e-6)
        # Nothing was clamped to the reference, and the run did not move the
        # reference either: it is a fact about the problem, not a minimum of histories.
        assert featured.reference == ProblemReference(0.0, 'optimum', 'analytic', MAPPING)
        assert problem.reference == ProblemReference(0.0, 'optimum', 'analytic', MAPPING)

    def test_histories_are_identical_without_a_reference(self):
        histories = []
        for reference in (True, False):
            featured = FeaturedProblem(self.problem(reference), Feature('plain'), 10, 0)
            for x in (np.array([1.0]), self.TOLERATED, self.PENALIZED):
                featured.fun(x)
            histories.append((featured.fun_hist, featured.maxcv_hist, featured.fun_init, featured.maxcv_init))
        for with_reference, without_reference in zip(*histories):
            np.testing.assert_array_equal(with_reference, without_reference)
        np.testing.assert_array_equal(histories[0][0], [1e7, -5e-4, -1e3])  # the comparison is not vacuous

    @pytest.mark.parametrize('kind', KINDS)
    def test_no_kind_is_a_floor(self, kind):
        # Not even a rigorous lower bound bounds the merits of infeasible points.
        featured = FeaturedProblem(Problem(steep_line, [1.0], cub=nonnegative, reference=record(kind=kind)),
                                   Feature('plain'), 10, 0)
        featured.fun(self.PENALIZED)
        merit = _default_merit(featured.fun_hist[-1], featured.maxcv_hist[-1], featured.maxcv_init)
        assert merit < featured.reference.merit


def dive_then_descend(fun, x0, *constraints):
    """Visits one slightly infeasible point, then does a short coordinate search."""
    fun(np.asarray(x0, dtype=float) - 1.0001)
    best_x, best_f = np.array(x0, dtype=float), fun(x0)
    for i in range(len(x0)):
        for delta in (0.5, -0.5, 0.25):
            trial = best_x.copy()
            trial[i] += delta
            value = fun(trial)
            if value < best_f:
                best_x, best_f = trial, value
    return best_x


def stay(fun, x0, *constraints):
    fun(x0)
    return x0


LIBRARY_SOURCE = '''
import numpy as np
from optiprofiler import Problem

WITH_REFERENCE = {with_reference!r}
MAPPING = 'feasible_objective/1'

def steep_line(x):
    return 1e7 * float(x[0])

def nonnegative(x):
    return np.array([-x[0]])

def sphere(x):
    return float(np.sum((x - 1.0) ** 2))

def reference_invariance_select(options):
    return ['STEEP', 'SPHERE']

def reference_invariance_load(name):
    if name == 'STEEP':
        reference = dict(merit=0.0, kind='optimum', source='analytic', mapping=MAPPING)
        return Problem(steep_line, [1.0], cub=nonnegative, name=name, reference=reference if WITH_REFERENCE else None)
    reference = dict(merit=-1.0, kind='lower_bound', source='analytic', mapping=MAPPING)
    return Problem(sphere, [0.0, 2.0], name=name, reference=reference if WITH_REFERENCE else None)
'''

ARRAYS = ('fun_histories', 'maxcv_histories', 'fun_outs', 'maxcv_outs', 'fun_inits', 'maxcv_inits', 'n_evals',
          'merit_histories', 'merit_outs', 'merit_inits')


def same_tree(first, second):
    """NaN-safe deep equality of the nested tuples, lists, dicts and arrays that ``benchmark`` returns."""
    if isinstance(first, dict):
        return isinstance(second, dict) and first.keys() == second.keys() and \
            all(same_tree(first[key], second[key]) for key in first)
    if isinstance(first, (list, tuple)):
        return isinstance(second, (list, tuple)) and len(first) == len(second) and \
            all(same_tree(a, b) for a, b in zip(first, second))
    if first is None or second is None:
        return first is second
    return np.array_equal(np.asarray(first), np.asarray(second), equal_nan=True)


class TestBenchmarkInvariance:

    @pytest.mark.parametrize('feature', ['plain', 'noisy', 'permuted+noisy'])
    def test_benchmark_is_identical_with_and_without_references(self, tmp_path, monkeypatch, feature):
        # Nothing consumes the record: scores, curves and archived histories
        # of a whole benchmark are the same whether the provider states
        # references or not. The two providers have the same name and the same
        # problems; only the ``reference`` argument differs.
        monkeypatch.chdir(tmp_path)
        outputs = {}
        for with_reference in (True, False):
            root = tmp_path / ('with' if with_reference else 'without')
            (root / 'libraries' / 'reference_invariance').mkdir(parents=True)
            (root / 'libraries' / 'reference_invariance' / 'reference_invariance_tools.py').write_text(
                LIBRARY_SOURCE.format(with_reference=with_reference), encoding='utf-8')
            result = benchmark([dive_then_descend, stay], plibs=['reference_invariance'],
                               custom_problem_libs_path=str(root / 'libraries'), ptype='un', mindim=1, maxdim=2,
                               feature_name=feature, n_runs=2, seed=5, n_jobs=1, max_eval_factor=8, max_tol_order=2,
                               score_only=False, draw_hist_plots='none', silent=True, solver_names=['dive', 'stay'],
                               benchmark_id='invariance', savepath=str(root))
            archives = sorted(root.rglob('data_for_loading.h5'))
            assert len(archives) == 1
            outputs[with_reference] = (result, load_results_from_h5(str(archives[0]))[0])

        (with_result, with_archive), (without_result, without_archive) = outputs[True], outputs[False]
        assert same_tree(with_result, without_result)
        assert list(with_archive['problem_names']) == list(without_archive['problem_names']) == ['STEEP', 'SPHERE'] \
            or sorted(with_archive['problem_names']) == ['SPHERE', 'STEEP']
        for key in ARRAYS:
            np.testing.assert_array_equal(with_archive[key], without_archive[key], err_msg=key)

        # The constrained problem was really run below its reference (0), and
        # the archive keeps those merits as they are: never clamped.
        steep = list(with_archive['problem_names']).index('STEEP')
        steep_merits = np.asarray(with_archive['merit_histories'][steep], dtype=float)
        assert np.nanmin(steep_merits) < 0.0
        assert np.nanmin(np.asarray(with_archive['fun_histories'][steep], dtype=float)) < 0.0

    def test_benchmark_has_no_reference_option(self):
        # The public signature is unchanged: no option selects, overrides or
        # consumes a reference.
        assert not [name for name in inspect.signature(benchmark).parameters if 'reference' in name.lower()]
        assert not [option.value for option in optiprofiler.utils.ProfileOption if 'reference' in option.value] \
            and not [option.value for option in optiprofiler.utils.ProblemOption if 'reference' in option.value]


def legacy_newobj_stream(state):
    """
    The protocol-2 byte stream that a ``ProblemReference`` of the superseded
    layout produced: ``GLOBAL`` of the class, ``NEWOBJ`` without arguments,
    then the state dictionary handed to ``__setstate__`` by ``BUILD``. It is
    assembled by hand because the pickler refuses to emit ``NEWOBJ`` for an
    object of another class.
    """
    header = pickle.PROTO + bytes([2])
    body = pickle.dumps(state, 2)
    assert body.startswith(header) and body.endswith(pickle.STOP)
    return (header + pickle.GLOBAL + b'optiprofiler.opclasses\nProblemReference\n' + pickle.EMPTY_TUPLE + pickle.NEWOBJ
            + body[len(header):-len(pickle.STOP)] + pickle.BUILD + pickle.STOP)


class LegacyLayoutStream:
    """
    Stands in for a ``ProblemReference`` of the superseded layout inside a
    pickled object: the instance is created without ``__init__`` and the
    five-field state is handed to ``__setstate__`` (the protocol 0/1 shape of
    the same stream, which the pickler accepts for a foreign class).
    """

    def __init__(self, state):
        self.state = state

    def __reduce_ex__(self, protocol):
        return copyreg._reconstructor, (ProblemReference, object, None), self.state


class RestoredRecordStream:
    """Pickles as the current layout does, with an arbitrary record: what another version may have written."""

    def __init__(self, stored):
        self.stored = stored

    def __reduce__(self):
        return opclasses._restore_problem_reference, (self.stored,)


LEGACY_STATE = {'fun': 0.0, 'maxcv': 0.0, 'kind': 'optimum', 'source': 'author', 'point': [1.0, 2.0]}


class TestLegacyAndProviderLoad:

    def test_pickle_round_trip(self):
        problem = constrained_problem(record(merit=0.5, kind='best_known', source='catalog'))
        clone = pickle.loads(pickle.dumps(problem))
        assert clone.reference == problem.reference and clone.reference is not problem.reference
        assert clone.fun(CENTER) == 0.0
        for protocol in range(2, pickle.HIGHEST_PROTOCOL + 1):
            assert pickle.loads(pickle.dumps(problem.reference, protocol)) == problem.reference

    def test_objects_from_before_the_record_load_as_unknown(self):
        problem = Problem(shifted_sphere, [0.0, 0.0])
        del problem.__dict__['_reference']  # the attribute layout of a Problem pickled before the record existed
        restored = pickle.loads(pickle.dumps(problem))
        assert '_reference' not in restored.__dict__ and restored.reference is None
        featured = FeaturedProblem(restored, Feature('noisy'), 5, 0)
        assert featured.reference is None and np.isfinite(featured.fun(restored.x0))
        assert FeaturedProblem(restored, Feature('noisy+permuted'), 5, 0).reference is None

    def test_superseded_layout_stream_is_an_error(self):
        # Such a stream restores instance state directly. Its ``fun`` is not
        # read as ``merit``: the stream is refused as a whole, alone or inside
        # a pickled problem.
        with pytest.raises(TypeError, match='superseded record layout'):
            pickle.loads(legacy_newobj_stream(LEGACY_STATE))
        with pytest.raises(TypeError, match='superseded record layout'):
            pickle.loads(pickle.dumps(LegacyLayoutStream(LEGACY_STATE)))
        problem = Problem(shifted_sphere, [0.0, 0.0])
        problem.__dict__['_reference'] = LegacyLayoutStream(LEGACY_STATE)
        with pytest.raises(TypeError, match='never reinterpreted'):
            pickle.loads(pickle.dumps(problem))
        # Even a state made of the four current fields is refused on this
        # path: a current record never arrives through ``__setstate__``.
        with pytest.raises(TypeError, match='superseded record layout'):
            pickle.loads(legacy_newobj_stream(record()))

    @pytest.mark.parametrize('stored', [
        record(mapping='feasible_objective/2'),   # a token of a later registry
        record(mapping='penalized_merit/1'),
        record(merit=float('nan')),
        record(merit=float('inf')),
        record(kind='certified_bound'),
        record(source=''),
        LEGACY_STATE,
        dict(record(), point=[1.0, 2.0]),
        {key: value for key, value in record().items() if key != 'mapping'},
        0.0,
        None,
        ['merit', 0.0],
    ], ids=repr)
    def test_rejected_records_deserialize_as_unknown(self, stored):
        # Unknown, never reinterpreted: the problem stays loadable and usable,
        # and no field is guessed, defaulted or borrowed from another layout.
        with pytest.warns(RuntimeWarning, match='read as unknown and is not reinterpreted'):
            assert pickle.loads(pickle.dumps(RestoredRecordStream(stored))) is None
        problem = Problem(shifted_sphere, [0.0, 0.0])
        problem.__dict__['_reference'] = RestoredRecordStream(stored)
        with pytest.warns(RuntimeWarning, match='read as unknown'):
            restored = pickle.loads(pickle.dumps(problem))
        assert restored.reference is None and restored.fun(CENTER) == 0.0
        assert FeaturedProblem(restored, Feature('permuted'), 5, 0).reference is None

    def test_valid_record_restores_without_a_warning(self):
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            assert pickle.loads(pickle.dumps(RestoredRecordStream(record()))) == ProblemReference.from_record(record())

    @pytest.mark.parametrize('stored', [record(), LEGACY_STATE, 0.0, 'optimum', object()], ids=repr)
    def test_only_a_validated_record_is_ever_handed_out(self, stored):
        # Whatever ends up in the instance dictionary, the property returns a
        # validated record or None.
        problem = Problem(shifted_sphere, [0.0, 0.0])
        problem.__dict__['_reference'] = stored
        assert problem.reference is None
        assert FeaturedProblem(problem, Feature('noisy'), 5, 0).reference is None
        assert FeaturedProblem(problem, Feature('noisy+permuted'), 5, 0).reference is None

    def test_project_x0_keeps_the_reference(self):
        problem = Problem(shifted_sphere, [-1.0, 5.0], xl=[0.0, 0.0], xu=[3.0, 3.0], reference=record())
        problem.project_x0()
        np.testing.assert_array_equal(problem.x0, [0.0, 3.0])
        assert problem.reference == ProblemReference(0.0, 'optimum', 'author', MAPPING)

    def test_provider_load_preserves_the_reference_without_executing_anything(self, tmp_path, monkeypatch, caplog):
        library_dir = tmp_path / 'reference_toy'
        library_dir.mkdir()
        (library_dir / 'reference_toy_tools.py').write_text('\n'.join([
            'from optiprofiler import Problem',
            '',
            'def never_evaluated(x):',
            "    raise RuntimeError('the objective was evaluated while loading')",
            '',
            'def reference_toy_select(options):',
            "    return ['REF', 'PLAIN', 'NAN', 'LEGACY', 'FUTURE', 'CALLBACK', 'NAKED']",
            '',
            'def reference_toy_load(problem_name):',
            "    good = {'merit': 0.0, 'kind': 'optimum', 'source': 'catalog:reference_toy', 'mapping': 'feasible_objective/1'}",
            "    reference = {'REF': good,",
            "                 'PLAIN': None,",
            "                 'NAN': dict(good, merit=float('nan')),",
            "                 'LEGACY': {'fun': 0.0, 'kind': 'optimum', 'source': 'catalog:reference_toy', 'point': [1.0, 2.0]},",
            "                 'FUTURE': dict(good, mapping='feasible_objective/2'),",
            "                 'CALLBACK': dict(good, mapping=never_evaluated),",
            "                 'NAKED': 0.0}[problem_name]",
            '    return Problem(never_evaluated, [0.0, 0.0], name=problem_name, reference=reference)',
        ]), encoding='utf-8')

        def no_solver(*args, **kwargs):
            raise AssertionError('a solver was executed while loading a problem')

        monkeypatch.setattr(optiprofiler.profiles, '_solve_one_problem', no_solver)
        plugin = load_problem_library(resolve_problem_library('reference_toy', tmp_path))
        problem = plugin.load('REF', {})
        assert type(problem) is Problem
        assert problem.reference == ProblemReference(0.0, 'optimum', 'catalog:reference_toy', MAPPING)
        assert problem.reference.source == 'catalog:reference_toy'  # provenance survives the load
        assert plugin.load('PLAIN', {}).reference is None
        # A malformed record makes the load fail; the provider is never
        # "helped" by a repaired or reinterpreted record.
        for name, error, message in [('NAN', ValueError, 'must be finite'),
                                     ('LEGACY', ValueError, 'superseded record layout'),
                                     ('FUTURE', ValueError, 'closed registry'),
                                     ('CALLBACK', TypeError, 'callables and user-defined mappings'),
                                     ('NAKED', TypeError, 'naked scalar')]:
            with pytest.raises(error, match=message):
                plugin.load(name, {})
        # The loaded objective is still the raising callable (Problem.fun
        # turns its exception into a logged warning and NaN), so nothing
        # evaluated it during the load.
        with caplog.at_level(logging.WARNING, logger='optiprofiler.opclasses'):
            assert np.isnan(problem.fun(problem.x0))
        assert 'evaluated while loading' in caplog.text
