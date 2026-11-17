"""
One-way provenance of a feature specification and its experiment plans, and the
versioned readers of the archived payloads.

The writer consumes the explicit specification and plan objects (never private
dictionaries), encodes option values with the shared safe encoder (callables are
described by name, nothing is executed) and produces the ``feature_pipeline-v3``
payload: a feature block with local-only stage options and an experiment block
with the authoritative run settings. Payloads written by earlier versions
(``feature_pipeline-v1``, ``feature_pipeline-v2``) are returned verbatim by the
reader with their own schema; nothing is rewritten and unknown facts stay unknown.
"""

import json

from .feature_definitions import SEED_POLICY_COMPOSED, SEED_POLICY_SINGLE, stamped_name
from .metadata import safe_metadata

FEATURE_PIPELINE_SCHEMA = 'feature_pipeline-v3'


def _describe(value):
    """Encode a value with the shared metadata encoder; never raise."""
    try:
        return safe_metadata(value)
    except Exception as exc:  # defensive: the encoder itself must not abort a benchmark
        return {'value': None, 'reason': 'options_not_described', 'error_type': type(exc).__name__}


def seed_policy(feature):
    """The seed policy of the execution of ``feature`` (recorded, not inferred from private state)."""
    return SEED_POLICY_COMPOSED if len(feature.stages) > 1 else SEED_POLICY_SINGLE


def full_feature_stamp(feature):
    """
    The complete stamp of ``feature``: the stamps of its effective stages
    joined with ``__`` in order (a one-stage feature keeps its established
    stamp; the identity is ``'plain'``). Never shortened: the bounded stamp
    used in paths is derived from it by ``profile_utils``.
    """
    if not feature.stages:
        return 'plain'
    return '__'.join(stage.stamp for stage in feature.stages)


def effective_specification(feature):
    """
    The ordered effective specification of ``feature`` as native data: one
    ``{'name', 'options'}`` mapping per effective stage with the validated
    local options after defaults, callables kept as the objects they are. Valid
    ``feature`` input; with the experiment count it reproduces the experiment.
    An identity pipeline is exported as one ``plain`` entry, which is valid
    input (an empty specification is not).
    """
    if not feature.stages:
        return [{'name': 'plain', 'options': {}}]
    return [{'name': stage.name, 'options': stage.native_options()} for stage in feature.stages]


def describe_feature(feature, feature_stamp=None, full_feature_stamp=None, input_route=None,
                     feature_stamp_origin=None):
    """
    The feature block of the provenance payload (plain, JSON-serializable data).

    ``route`` is the benchmark keyword that carried the specification in the
    described invocation (``'feature_name'``, ``'feature'``, or ``None`` when
    neither was given or no invocation is described). ``declaration_route`` is
    the declaration route of the specification itself: ``'feature_name'`` when
    it was declared by the shorthand string, ``'feature'`` when it was declared
    by structured entries, ``None`` (with ``declared`` empty and
    ``declared_name`` ``None``) when the specification was imported from a
    historical object that recorded no declaration. The declaration is a
    property of the specification and does not change with the keyword that
    carried it: a ``Feature`` declared once and passed as
    ``benchmark(feature=obj)`` keeps its declaration and its receipt.
    ``feature_stamp_origin`` says whether ``feature_stamp`` was given explicitly
    or generated (``None`` when no invocation is described).
    """
    return {
        'route': input_route,
        'declaration_route': feature.declared.route,
        'declared_name': feature.declared_name,
        'declared': [{'name': name, 'options': _describe(options)} for name, options in feature.declared.native_entries()],
        'effective_name': feature.name,
        'seed_policy': seed_policy(feature),
        'feature_stamp': feature_stamp,
        'feature_stamp_origin': feature_stamp_origin,
        'full_feature_stamp': full_feature_stamp,
        'stages': [{
            'position': stage.position,
            'name': stage.name,
            'code': stage.code,
            'occurrence': stage.occurrence,
            'identity': stage.identity,
            'options': _describe(stage.native_options()),
        } for stage in feature.stages],
    }


def describe_plan(plan):
    """The experiment block of the provenance payload."""
    return None if plan is None else plan.describe()


def feature_pipeline_payload(feature, plan, feature_stamp=None, full_feature_stamp=None, input_route=None,
                             feature_stamp_origin=None):
    """The ``feature_pipeline-v3`` payload written into archives and reports."""
    return {
        'schema': FEATURE_PIPELINE_SCHEMA,
        'feature': describe_feature(feature, feature_stamp, full_feature_stamp, input_route, feature_stamp_origin),
        'experiment': describe_plan(plan),
    }


def feature_pipeline_text(feature, plan, feature_stamp=None, full_feature_stamp=None, input_route=None,
                          feature_stamp_origin=None):
    """The payload as sorted JSON text (the archive stores text)."""
    return json.dumps(feature_pipeline_payload(feature, plan, feature_stamp, full_feature_stamp, input_route,
                                               feature_stamp_origin), sort_keys=True)


def read_feature_pipeline(value):
    """
    Read an archived ``feature_pipeline`` entry verbatim.

    Returns ``(payload, schema)``: the parsed payload exactly as stored (any
    known version, never migrated or relabelled) and its schema identifier, or
    ``(None, None)`` when the archive carries no entry, or an explicit reason
    record when the text cannot be parsed. Unknown schema identifiers are
    returned as stored with ``schema`` reporting the unknown identifier.
    """
    if value is None:
        return None, None
    if isinstance(value, (bytes, bytearray)):
        value = value.decode('utf-8', errors='replace')
    if isinstance(value, str):
        try:
            payload = json.loads(value)
        except ValueError:
            return {'value': None, 'reason': 'unparsable_feature_pipeline'}, None
    else:
        payload = value
    schema = payload.get('schema') if isinstance(payload, dict) else None
    return payload, schema


def archived_feature_name(result):
    """
    The name of the feature of an archived library record, for display.

    It is the effective name recorded in the ``feature_pipeline`` entry (any
    known version). An archive without that entry (OptiProfiler 1.x) only has
    its ``feature_stamp``: a generated stamp is the feature name followed by
    option values (``linearly_transformed_rotated``; the stages of a
    composition are joined with a double underscore), so the name is read from
    it. A stamp of another form (it can be set by the user) is returned as it
    is, since it is the only name the archive has. Returns ``None`` when the
    record has neither.
    """
    payload, _ = read_feature_pipeline(result.get('feature_pipeline'))
    if isinstance(payload, dict):
        # The feature block of the current schema; the first schema was flat.
        block = payload['feature'] if isinstance(payload.get('feature'), dict) else payload
        name = block.get('effective_name')
        if isinstance(name, str) and name:
            return name
    stamp = _archived_stamp(result)
    if not stamp:
        return None
    names = [stamped_name(part) for part in stamp.split('__')]
    return '+'.join(names) if all(names) else stamp


def archived_feature_names(result):
    """
    The feature names that an archived library record vouches for: the
    effective and declared names of its ``feature_pipeline`` entry, otherwise
    the name read from a generated ``feature_stamp``. Empty when the record
    has neither, or only a stamp of another form, so that nothing can be
    compared with it.
    """
    payload, _ = read_feature_pipeline(result.get('feature_pipeline'))
    if isinstance(payload, dict):
        block = payload['feature'] if isinstance(payload.get('feature'), dict) else payload
        names = {block.get(key) for key in ('effective_name', 'declared_name')}
        names = {name for name in names if isinstance(name, str) and name}
        if names:
            return names
    stamp = _archived_stamp(result)
    if stamp:
        names = [stamped_name(part) for part in stamp.split('__')]
        if all(names):
            return {'+'.join(names)}
    return set()


def _archived_stamp(result):
    stamp = result.get('feature_stamp')
    if isinstance(stamp, bytes):
        stamp = stamp.decode('utf-8', errors='replace')
    return stamp if isinstance(stamp, str) else None
