function errors = evalReportSchemaCheck(json_path, schema_path)
%EVALREPORTSCHEMACHECK Validate a report file against a JSON Schema (test fixture).
%   ERRORS = evalReportSchemaCheck(JSON_PATH, SCHEMA_PATH) returns a cell
%   array of 'path: message' strings; empty means valid.
%
%   The two schemas in doc/source/_static are the single Python/MATLAB
%   reader contract. This checker implements the same keyword subset as the
%   Python test helper eval_report_contract.py (type, enum, const, required,
%   properties, additionalProperties, items, prefixItems, min/maxItems,
%   uniqueItems, minimum, maximum, min/maxLength, pattern, anyOf, oneOf,
%   allOf, if/then and local $ref), so a differently spelled field fails
%   the MATLAB suite exactly as it fails the Python suite.
%
%   jsondecode cannot be used for this purpose: it turns a one-element
%   array into a scalar and null into [], which destroys the distinctions
%   the schema pins. A small tokenizer therefore builds a tagged tree:
%   struct('kind', 'object'|'array'|'string'|'number'|'boolean'|'null',
%   'value', ...), objects keep 'keys' (cellstr) and 'values' (cell).
    % Both files are UTF-8 by contract; decode explicitly rather than with
    % the platform default encoding (a Windows code page on Windows).
    document = parseJson(readUtf8(json_path));
    schema = parseJson(readUtf8(schema_path));
    errors = validate(document, schema, schema, '');
end

function errors = validate(node, schema, root, path)
    errors = {};
    if ~strcmp(schema.kind, 'object'), return; end
    ref = getKey(schema, '$ref');
    if ~isempty(ref)
        target = root;
        parts = strsplit(ref.value, '/');
        for k = 2:numel(parts), target = getKey(target, parts{k}); end
        errors = validate(node, target, root, path);
        return;
    end
    type_schema = getKey(schema, 'type');
    if ~isempty(type_schema)
        if strcmp(type_schema.kind, 'array'), expected = cellfun(@(t) t.value, type_schema.value, 'UniformOutput', false);
        else, expected = {type_schema.value}; end
        if ~any(cellfun(@(t) matchesType(node, t), expected))
            errors{end+1} = sprintf('%s: type %s is not %s', where(path), node.kind, strjoin(expected, '|'));
            return;
        end
    end
    constant = getKey(schema, 'const');
    if ~isempty(constant) && ~same(node, constant)
        errors{end+1} = sprintf('%s: %s is not the constant %s', where(path), describe(node), describe(constant));
    end
    enumeration = getKey(schema, 'enum');
    if ~isempty(enumeration) && ~any(cellfun(@(option) same(node, option), enumeration.value))
        errors{end+1} = sprintf('%s: %s is not one of the enumerated values', where(path), describe(node));
    end
    for keyword = {'anyOf', 'oneOf'}
        options = getKey(schema, keyword{1});
        if isempty(options), continue; end
        outcomes = cellfun(@(option) validate(node, option, root, path), options.value, 'UniformOutput', false);
        passing = cellfun(@isempty, outcomes);
        if ~any(passing) || (strcmp(keyword{1}, 'oneOf') && sum(passing) ~= 1)
            [~, closest] = min(cellfun(@numel, outcomes));
            errors{end+1} = sprintf('%s: %s matched %d alternatives (closest: %s)', where(path), keyword{1}, sum(passing), strjoin(outcomes{closest}(1:min(4, numel(outcomes{closest}))), '; '));
        end
    end
    all_of = getKey(schema, 'allOf');
    if ~isempty(all_of)
        for k = 1:numel(all_of.value), errors = [errors, validate(node, all_of.value{k}, root, path)]; end
    end
    condition = getKey(schema, 'if');
    if ~isempty(condition) && isempty(validate(node, condition, root, path))
        consequence = getKey(schema, 'then');
        if ~isempty(consequence), errors = [errors, validate(node, consequence, root, path)]; end
    end
    if strcmp(node.kind, 'object')
        required = getKey(schema, 'required');
        if ~isempty(required)
            for k = 1:numel(required.value)
                if ~ismember(required.value{k}.value, node.keys)
                    errors{end+1} = sprintf('%s: missing required key %s', where(path), required.value{k}.value);
                end
            end
        end
        properties = getKey(schema, 'properties');
        additional = getKey(schema, 'additionalProperties');
        for k = 1:numel(node.keys)
            key = node.keys{k}; child = node.values{k}; child_path = [path, '/', key];
            property = [];
            if ~isempty(properties), property = getKey(properties, key); end
            if ~isempty(property)
                errors = [errors, validate(child, property, root, child_path)];
            elseif ~isempty(additional)
                if strcmp(additional.kind, 'boolean') && ~additional.value
                    errors{end+1} = sprintf('%s: unexpected key %s', where(path), key);
                elseif strcmp(additional.kind, 'object')
                    errors = [errors, validate(child, additional, root, child_path)];
                end
            end
        end
        limit = getKey(schema, 'maxProperties');
        if ~isempty(limit) && numel(node.keys) > limit.value
            errors{end+1} = sprintf('%s: too many properties', where(path));
        end
    end
    if strcmp(node.kind, 'array')
        n = numel(node.value);
        limit = getKey(schema, 'minItems');
        if ~isempty(limit) && n < limit.value, errors{end+1} = sprintf('%s: fewer than %d items', where(path), limit.value); end
        limit = getKey(schema, 'maxItems');
        if ~isempty(limit) && n > limit.value, errors{end+1} = sprintf('%s: more than %d items', where(path), limit.value); end
        unique_items = getKey(schema, 'uniqueItems');
        if ~isempty(unique_items) && unique_items.value
            rendered = cellfun(@describe, node.value, 'UniformOutput', false);
            if numel(unique(rendered)) ~= n, errors{end+1} = sprintf('%s: items are not unique', where(path)); end
        end
        prefix = getKey(schema, 'prefixItems');
        items = getKey(schema, 'items');
        for k = 1:n
            child_path = sprintf('%s/%d', path, k - 1);
            if ~isempty(prefix) && k <= numel(prefix.value)
                errors = [errors, validate(node.value{k}, prefix.value{k}, root, child_path)];
            elseif ~isempty(items)
                errors = [errors, validate(node.value{k}, items, root, child_path)];
            end
        end
    end
    if strcmp(node.kind, 'number')
        limit = getKey(schema, 'minimum');
        if ~isempty(limit) && node.value < limit.value, errors{end+1} = sprintf('%s: %g is below %g', where(path), node.value, limit.value); end
        limit = getKey(schema, 'maximum');
        if ~isempty(limit) && node.value > limit.value, errors{end+1} = sprintf('%s: %g is above %g', where(path), node.value, limit.value); end
    end
    if strcmp(node.kind, 'string')
        limit = getKey(schema, 'minLength');
        if ~isempty(limit) && numel(node.value) < limit.value, errors{end+1} = sprintf('%s: string too short', where(path)); end
        limit = getKey(schema, 'maxLength');
        if ~isempty(limit) && numel(node.value) > limit.value, errors{end+1} = sprintf('%s: string too long', where(path)); end
        pattern = getKey(schema, 'pattern');
        if ~isempty(pattern) && isempty(regexp(node.value, pattern.value, 'once'))
            errors{end+1} = sprintf('%s: string does not match the pattern', where(path));
        end
    end
end

function yes = matchesType(node, expected)
    switch expected
        case 'number', yes = strcmp(node.kind, 'number');
        case 'integer', yes = strcmp(node.kind, 'number') && node.value == floor(node.value) && isfinite(node.value);
        otherwise, yes = strcmp(node.kind, expected);
    end
end

function yes = same(a, b)
    yes = false;
    if strcmp(a.kind, 'number') && strcmp(b.kind, 'number'), yes = a.value == b.value; return; end
    if ~strcmp(a.kind, b.kind), return; end
    switch a.kind
        case {'string'}, yes = strcmp(a.value, b.value);
        case {'boolean'}, yes = a.value == b.value;
        case {'null'}, yes = true;
        case {'array'}, yes = numel(a.value) == numel(b.value) && all(cellfun(@(x, y) same(x, y), a.value, b.value));
        case {'object'}
            yes = numel(a.keys) == numel(b.keys) && all(ismember(a.keys, b.keys));
            for k = 1:numel(a.keys)
                if ~yes, return; end
                yes = same(a.values{k}, getKey(b, a.keys{k}));
            end
    end
end

function value = getKey(node, key)
    value = [];
    if ~strcmp(node.kind, 'object'), return; end
    index = find(strcmp(node.keys, key), 1);
    if ~isempty(index), value = node.values{index}; end
end

function text = where(path)
    text = path; if isempty(text), text = '<root>'; end
end

function text = describe(node)
    switch node.kind
        case 'string', text = ['"', node.value, '"'];
        case 'number', text = sprintf('%.17g', node.value);
        case 'boolean', if node.value, text = 'true'; else, text = 'false'; end
        case 'null', text = 'null';
        case 'array', text = ['[', strjoin(cellfun(@describe, node.value, 'UniformOutput', false), ','), ']'];
        otherwise, text = ['{', strjoin(node.keys, ','), '}'];
    end
end

function node = parseJson(text)
    % Tokenize once with a regular expression, then build the tree with an
    % explicit stack so deeply nested documents do not recurse in MATLAB.
    tokens = regexp(text, '"(?:[^"\\]|\\.)*"|-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?|true|false|null|[{}\[\]:,]', 'match');
    stack = {}; keys = {}; pending_key = {};
    node = [];
    for t = 1:numel(tokens)
        token = tokens{t};
        switch token
            case '{'
                stack{end+1} = struct('kind', 'object', 'keys', {{}}, 'values', {{}}); keys{end+1} = ''; pending_key{end+1} = false;
                continue;
            case '['
                stack{end+1} = struct('kind', 'array', 'value', {{}}); keys{end+1} = ''; pending_key{end+1} = false;
                continue;
            case {'}', ']'}
                value = stack{end}; stack(end) = []; keys(end) = []; pending_key(end) = [];
            case {':', ','}
                continue;
            otherwise
                value = leaf(token);
        end
        if isempty(stack)
            node = value;
            continue;
        end
        parent = stack{end};
        if strcmp(parent.kind, 'object')
            if ~pending_key{end}
                keys{end} = value.value; pending_key{end} = true;
            else
                parent.keys{end+1} = keys{end}; parent.values{end+1} = value; pending_key{end} = false;
                stack{end} = parent;
            end
        else
            parent.value{end+1} = value; stack{end} = parent;
        end
    end
    if isempty(node), error('OptiProfiler:EvalReportSchemaCheck', 'Empty JSON document.'); end
end

function value = leaf(token)
    switch token
        case 'true', value = struct('kind', 'boolean', 'value', true);
        case 'false', value = struct('kind', 'boolean', 'value', false);
        case 'null', value = struct('kind', 'null', 'value', []);
        otherwise
            if token(1) == '"'
                value = struct('kind', 'string', 'value', unescape(token(2:end-1)));
            else
                value = struct('kind', 'number', 'value', str2double(token));
            end
    end
end

function text = unescape(text)
    if ~any(text == '\'), return; end
    text = regexprep(text, '\\u([0-9a-fA-F]{4})', '${char(hex2dec($1))}');
    text = regexprep(text, '\\([\\"/])', '$1');
    text = strrep(strrep(strrep(text, '\n', newline), '\t', sprintf('\t')), '\r', sprintf('\r'));
    text = strrep(strrep(text, '\b', sprintf('\b')), '\f', sprintf('\f'));
end

function text = readUtf8(path)
    fid = fopen(path, 'rb');
    assert(fid >= 0, 'OptiProfiler:EvalReportSchemaCheck', 'Cannot open %s', path);
    guard = onCleanup(@() fclose(fid));
    text = native2unicode(reshape(fread(fid, '*uint8'), 1, []), 'UTF-8');
end
