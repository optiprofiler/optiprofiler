"""Dependency-free checker for the EvalReport JSON Schemas used by the tests.

The two schemas in ``doc/source/_static`` are the single Python/MATLAB reader
contract. The MATLAB test fixture ``evalReportSchemaCheck.m`` implements the
same subset, so a field spelled differently by either emitter fails that
language's own suite. The subset covers exactly the keywords the schemas use:
type, enum, const, required, properties, additionalProperties, items,
prefixItems, minItems, maxItems, uniqueItems, minimum, maximum, minLength,
maxLength, pattern, anyOf, oneOf, allOf, if/then and local ``$ref``.

When ``jsonschema`` is importable the tests additionally validate with it;
this module is the check that always runs.
"""
import json
import math
import re
from pathlib import Path


SCHEMA_DIR = Path(__file__).resolve().parents[3] / 'doc' / 'source' / '_static'


def load_schema(name):
    path = SCHEMA_DIR / name
    if not path.is_file():
        # The schemas are documentation assets of the repository, not package
        # data: these tests run from a checkout, never from an installed wheel.
        raise FileNotFoundError(f'{path} is missing; run the EvalReport tests from a repository checkout')
    with path.open(encoding='utf-8') as stream:
        return json.load(stream, parse_constant=lambda token: (_ for _ in ()).throw(ValueError(token)))


def _type_of(value):
    if value is None:
        return 'null'
    if isinstance(value, bool):
        return 'boolean'
    if isinstance(value, int):
        return 'integer'
    if isinstance(value, float):
        return 'number'
    if isinstance(value, str):
        return 'string'
    if isinstance(value, list):
        return 'array'
    if isinstance(value, dict):
        return 'object'
    raise TypeError(type(value))


def _same(value, expected):
    """JSON equality: 1 and 1.0 are the same number, but True is not 1."""
    if isinstance(value, bool) or isinstance(expected, bool):
        return type(value) is type(expected) and value == expected
    if isinstance(value, (int, float)) and isinstance(expected, (int, float)):
        return value == expected
    if isinstance(value, list) and isinstance(expected, list):
        return len(value) == len(expected) and all(_same(a, b) for a, b in zip(value, expected))
    if isinstance(value, dict) and isinstance(expected, dict):
        return value.keys() == expected.keys() and all(_same(value[k], expected[k]) for k in value)
    return type(value) is type(expected) and value == expected


def _matches_type(value, expected):
    kind = _type_of(value)
    if expected == 'number':
        return kind in ('number', 'integer')
    if expected == 'integer':
        return kind == 'integer' or (kind == 'number' and math.isfinite(value) and value == int(value))
    return kind == expected


def _errors(value, schema, root, path):
    """Yield ``(path, message)`` pairs; an empty result means valid."""
    if '$ref' in schema:
        target = root
        for part in schema['$ref'].lstrip('#/').split('/'):
            target = target[part]
        yield from _errors(value, target, root, path)
        return
    if 'type' in schema:
        expected = schema['type'] if isinstance(schema['type'], list) else [schema['type']]
        if not any(_matches_type(value, item) for item in expected):
            yield path, f'type {_type_of(value)} is not {expected}'
            return
    if 'const' in schema and not _same(value, schema['const']):
        yield path, f'{value!r} is not the constant {schema["const"]!r}'
    if 'enum' in schema and not any(_same(value, item) for item in schema['enum']):
        yield path, f'{value!r} is not one of {schema["enum"]!r}'
    for key in ('anyOf', 'oneOf'):
        if key in schema:
            outcomes = [list(_errors(value, option, root, path)) for option in schema[key]]
            passing = [i for i, errors in enumerate(outcomes) if not errors]
            if not passing or (key == 'oneOf' and len(passing) != 1):
                closest = min(outcomes, key=len) if outcomes else []
                detail = '; '.join(f'{"/".join(str(p) for p in p_)}: {m}' for p_, m in closest[:4])
                yield path, f'{key} matched {len(passing)} alternatives (closest: {detail})'
    for option in schema.get('allOf', []):
        yield from _errors(value, option, root, path)
    if 'if' in schema and not list(_errors(value, schema['if'], root, path)):
        yield from _errors(value, schema.get('then', {}), root, path)
    if isinstance(value, dict):
        for key in schema.get('required', []):
            if key not in value:
                yield path, f'missing required key {key!r}'
        properties = schema.get('properties', {})
        for key, item in value.items():
            if key in properties:
                yield from _errors(item, properties[key], root, path + [key])
            elif 'additionalProperties' in schema:
                extra = schema['additionalProperties']
                if extra is False:
                    yield path, f'unexpected key {key!r}'
                elif isinstance(extra, dict):
                    yield from _errors(item, extra, root, path + [key])
        if 'maxProperties' in schema and len(value) > schema['maxProperties']:
            yield path, 'too many properties'
    if isinstance(value, list):
        if 'minItems' in schema and len(value) < schema['minItems']:
            yield path, f'fewer than {schema["minItems"]} items'
        if 'maxItems' in schema and len(value) > schema['maxItems']:
            yield path, f'more than {schema["maxItems"]} items'
        if schema.get('uniqueItems') and len({json.dumps(v, sort_keys=True) for v in value}) != len(value):
            yield path, 'items are not unique'
        prefix = schema.get('prefixItems', [])
        for index, item in enumerate(value):
            if index < len(prefix):
                yield from _errors(item, prefix[index], root, path + [index])
            elif 'items' in schema:
                yield from _errors(item, schema['items'], root, path + [index])
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if 'minimum' in schema and value < schema['minimum']:
            yield path, f'{value} is below {schema["minimum"]}'
        if 'maximum' in schema and value > schema['maximum']:
            yield path, f'{value} is above {schema["maximum"]}'
    if isinstance(value, str):
        if 'minLength' in schema and len(value) < schema['minLength']:
            yield path, 'string too short'
        if 'maxLength' in schema and len(value) > schema['maxLength']:
            yield path, 'string too long'
        if 'pattern' in schema and not re.search(schema['pattern'], value):
            yield path, f'{value!r} does not match {schema["pattern"]!r}'


def validation_errors(document, schema):
    return [f'{"/".join(str(p) for p in path) or "<root>"}: {message}'
            for path, message in _errors(document, schema, schema, [])]


def assert_valid(document, schema_name):
    """Validate with the built-in subset and, when available, with jsonschema."""
    schema = load_schema(schema_name)
    errors = validation_errors(document, schema)
    assert not errors, '\n'.join(errors[:20])
    try:
        import jsonschema
    except ImportError:  # pragma: no cover - optional cross-check only
        return
    validator = jsonschema.Draft202012Validator(schema)
    problems = [f'{list(error.absolute_path)}: {error.message}' for error in validator.iter_errors(document)]
    assert not problems, '\n'.join(problems[:20])
