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

from .feature_definitions import SEED_POLICY_COMPOSED, SEED_POLICY_SINGLE
from .metadata import safe_metadata

FEATURE_PIPELINE_SCHEMA = 'feature_pipeline-v3'
KNOWN_PIPELINE_SCHEMAS = ('feature_pipeline-v1', 'feature_pipeline-v2', FEATURE_PIPELINE_SCHEMA)


def _describe(value):
    """Encode a value with the shared metadata encoder; never raise."""
    try:
        return safe_metadata(value)
    except Exception as exc:  # defensive: the encoder itself must not abort a benchmark
        return {'value': None, 'reason': 'options_not_described', 'error_type': type(exc).__name__}


def seed_policy(feature):
    """The seed policy of the execution of ``feature`` (recorded, not inferred from private state)."""
    return SEED_POLICY_COMPOSED if len(feature.stages) > 1 else SEED_POLICY_SINGLE


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
    return [{'name': stage.name, 'options': dict(stage.options)} for stage in feature.stages]


def describe_feature(feature, feature_stamp=None, full_feature_stamp=None):
    """The feature block of the provenance payload (plain, JSON-serializable data)."""
    return {
        'route': feature.declared.route,
        'declared_name': feature.declared_name,
        'declared': [{'name': name, 'options': _describe(dict(options))} for name, options in feature.declared.entries],
        'effective_name': feature.name,
        'seed_policy': seed_policy(feature),
        'feature_stamp': feature_stamp,
        'full_feature_stamp': full_feature_stamp,
        'stages': [{
            'position': position,
            'name': stage.name,
            'code': stage.code,
            'occurrence': stage.occurrence,
            'identity': stage.identity,
            'options': _describe(dict(stage.options)),
        } for position, stage in enumerate(feature.stages)],
    }


def describe_plan(plan):
    """The experiment block of the provenance payload."""
    return None if plan is None else plan.describe()


def feature_pipeline_payload(feature, plan, feature_stamp=None, full_feature_stamp=None):
    """The ``feature_pipeline-v3`` payload written into archives and reports."""
    return {
        'schema': FEATURE_PIPELINE_SCHEMA,
        'feature': describe_feature(feature, feature_stamp, full_feature_stamp),
        'experiment': describe_plan(plan),
    }


def feature_pipeline_text(feature, plan, feature_stamp=None, full_feature_stamp=None):
    """The payload as sorted JSON text (the archive stores text)."""
    return json.dumps(feature_pipeline_payload(feature, plan, feature_stamp, full_feature_stamp), sort_keys=True)


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
