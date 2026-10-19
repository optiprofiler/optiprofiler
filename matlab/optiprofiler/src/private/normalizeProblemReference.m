function reference = normalizeProblemReference(reference)
%NORMALIZEPROBLEMREFERENCE Validate raw objective/violation reference facts.
%   Required fields: fun, maxcv, kind, source. Optional x_ref is source
%   metadata in the original problem's coordinates, never an oracle callback.
%   feasible requires maxcv=0; merit may have finite nonnegative violation.
%   Old kinds and mapping tokens are read for provider/archive compatibility.
%   Empty means unknown. Validation evaluates no objective or constraint.
%   The exact earlier four-field draft with feasible_objective/1 is accepted
%   as (fun=merit,maxcv=0); no arbitrary precomputed merit is reinterpreted.
    if isempty(reference)
        reference = [];
        return
    end
    required = {'fun', 'maxcv', 'kind', 'source'};
    if ~(isstruct(reference) && isscalar(reference))
        error("MATLAB:Problem:reference_NotStruct", "The reference must be a scalar struct with fields %s.", strjoin(required, ', '));
    end
    fields = reshape(fieldnames(reference), 1, []);
    if numel(fields) == 4 && all(ismember(fields, {'merit', 'kind', 'source', 'mapping'}))
        if isa(reference.mapping, 'function_handle')
            error("MATLAB:Problem:reference_mapping_NotToken", "Reference mapping must be a token, not a function handle.");
        end
        if ~strcmp(referenceText(reference.mapping, 'mapping'), 'feasible_objective/1')
            error("MATLAB:Problem:reference_mapping_Unknown", "A legacy merit record requires feasible_objective/1; its meaning cannot be guessed.");
        end
        reference.fun = reference.merit;
        reference.maxcv = 0;
        reference = rmfield(reference, 'merit');
        fields = reshape(fieldnames(reference), 1, []);
    end
    if ismember('point', fields)
        error("MATLAB:Problem:reference_LegacyField", "The superseded point field is not accepted; use original-coordinate x_ref.");
    end
    unknown = fields(~ismember(fields, [required, {'x_ref','mapping'}]));
    if ~isempty(unknown)
        error("MATLAB:Problem:reference_UnknownField", "Unknown problem reference fields: %s.", strjoin(unknown, ', '));
    end
    missing = required(~ismember(required, fields));
    if ~isempty(missing)
        error("MATLAB:Problem:reference_MissingField", "The problem reference is missing required fields: %s.", strjoin(missing, ', '));
    end
    fun = referenceScalar(reference.fun, 'fun');
    maxcv = referenceScalar(reference.maxcv, 'maxcv');
    if maxcv < 0
        error("MATLAB:Problem:reference_maxcv_Negative", "The reference maxcv must be nonnegative.");
    end
    kind = referenceText(reference.kind, 'kind');
    if ismember(kind, {'lower_bound','optimum'}), kind = 'feasible'; end
    if ismember(kind, {'best_known','target'}), kind = 'merit'; end
    if ~ismember(kind, {'feasible', 'merit'})
        error("MATLAB:Problem:reference_kind_Unknown", "Unknown reference kind '%s'.", kind);
    end
    source = referenceText(reference.source, 'source');
    if isempty(strtrim(source))
        error("MATLAB:Problem:reference_source_Empty", "The reference source must be non-empty provenance text.");
    end
    mapping = 'objective_violation/1';
    if isfield(reference,'mapping')
        if isa(reference.mapping, 'function_handle')
            error("MATLAB:Problem:reference_mapping_NotToken", "Reference mapping must be a token, not a function handle.");
        end
        mapping = referenceText(reference.mapping, 'mapping');
    end
    if ~ismember(mapping, {'feasible_objective/1', 'objective_violation/1'})
        error("MATLAB:Problem:reference_mapping_Unknown", "Unknown reference mapping '%s'.", mapping);
    end
    if maxcv ~= 0 && (strcmp(kind,'feasible') || strcmp(mapping, 'feasible_objective/1'))
        error("MATLAB:Problem:reference_maxcv_NotFeasible", "Reference maxcv must be zero for feasible and legacy feasible_objective/1 records.");
    end
    point = [];
    if isfield(reference,'x_ref') && ~isempty(reference.x_ref)
        point = reference.x_ref;
        if ~(isnumeric(point) && isreal(point) && isvector(point) && all(isfinite(point)))
            error('MATLAB:Problem:reference_x_ref_Invalid', 'The reference x_ref must be a finite real vector.');
        end
        point = full(double(point(:)));
    end
    reference = struct('fun', fun, 'maxcv', maxcv, 'kind', kind, 'source', source);
    if ~isempty(point), reference.x_ref = point; end
end

function value = referenceScalar(input, field)
    if ~(isnumeric(input) && isreal(input) && isscalar(input))
        error("MATLAB:Problem:reference_" + field + "_NotRealScalar", "The reference %s must be a real scalar.", field);
    end
    value = full(double(input));
    if ~isfinite(value)
        error("MATLAB:Problem:reference_" + field + "_NotFinite", "The reference %s must be finite.", field);
    end
    % Cast saturation makes a round trip alone insufficient at intmax.
    if isinteger(input) && abs(value) >= flintmax ...
            && ~(value < double(intmax(class(input))) && cast(value, class(input)) == input)
        error("MATLAB:Problem:reference_" + field + "_NotExact", "The reference %s cannot be represented exactly as a double.", field);
    end
end

function value = referenceText(value, field)
    if isstring(value) && isscalar(value) && ~ismissing(value)
        value = char(value);
    end
    if ~(ischar(value) && (isrow(value) || isempty(value)))
        error("MATLAB:Problem:reference_" + field + "_NotText", "The reference %s must be a character vector or string scalar.", field);
    end
    value = reshape(value, 1, []);
end
