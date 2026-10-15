function reference = normalizeProblemReference(reference)
%NORMALIZEPROBLEMREFERENCE Validate raw objective/violation reference facts.
%   Canonical fields: fun, maxcv, kind, source, mapping. No point or callback
%   is stored. lower_bound/optimum require maxcv=0; best_known/target may
%   have a finite nonnegative violation. Scoring uses the current merit_fun.
%   Empty means unknown. Validation evaluates no objective or constraint.
%   The exact earlier four-field draft with feasible_objective/1 is accepted
%   as (fun=merit,maxcv=0); no arbitrary precomputed merit is reinterpreted.
    if isempty(reference)
        reference = [];
        return
    end
    required = {'fun', 'maxcv', 'kind', 'source', 'mapping'};
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
        error("MATLAB:Problem:reference_LegacyField", "Reference points are not accepted or transported.");
    end
    unknown = fields(~ismember(fields, required));
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
    if ~ismember(kind, {'lower_bound', 'optimum', 'best_known', 'target'})
        error("MATLAB:Problem:reference_kind_Unknown", "Unknown reference kind '%s'.", kind);
    end
    source = referenceText(reference.source, 'source');
    if isempty(strtrim(source))
        error("MATLAB:Problem:reference_source_Empty", "The reference source must be non-empty provenance text.");
    end
    if isa(reference.mapping, 'function_handle')
        error("MATLAB:Problem:reference_mapping_NotToken", "Reference mapping must be a token, not a function handle.");
    end
    mapping = referenceText(reference.mapping, 'mapping');
    if ~ismember(mapping, {'feasible_objective/1', 'objective_violation/1'})
        error("MATLAB:Problem:reference_mapping_Unknown", "Unknown reference mapping '%s'.", mapping);
    end
    if maxcv ~= 0 && (ismember(kind, {'lower_bound', 'optimum'}) || strcmp(mapping, 'feasible_objective/1'))
        error("MATLAB:Problem:reference_maxcv_NotFeasible", "Reference maxcv must be zero for lower_bound/optimum and feasible_objective/1.");
    end
    reference = struct('fun', fun, 'maxcv', maxcv, 'kind', kind, 'source', source, 'mapping', mapping);
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
