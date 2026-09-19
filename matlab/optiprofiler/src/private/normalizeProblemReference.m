function reference = normalizeProblemReference(reference)
%NORMALIZEPROBLEMREFERENCE Validate a feasible reference fact of a problem.
%   REFERENCE = normalizeProblemReference(REFERENCE) returns [] for an empty
%   input (an unknown reference) or a scalar struct with exactly the four
%   fields `merit`, `kind`, `source` and `mapping`, in this order. Anything
%   else raises an error; nothing is repaired, defaulted or guessed.
%
%   The record is one scalar stated by the author or provider of the problem,
%   with its kind, its provenance and the mapping that fixes how the scalar is
%   read. It holds no point, no function handle and no constraint violation.
%
%   - merit: a finite real scalar (stored as double). Logical, text, complex,
%     non-scalar, NaN and infinite values are rejected, and so is an integer
%     that a double cannot hold exactly: the stored scalar is exactly the
%     scalar that was given.
%   - kind: 'lower_bound', 'optimum', 'best_known' or 'target'. Every kind is
%     a claim over the FEASIBLE points of the problem: a bound on the
%     objective over them, its exact optimal value over them, the objective
%     value of a known feasible point, or a target level chosen for them.
%   - source: non-empty provenance text.
%   - mapping: a token of the closed registry {'feasible_objective/1'}. With
%     this token, `merit` is an objective value over feasible points, so it
%     equals the merit of a feasible point under every merit function with
%     the feasible identity merit_fun(f, 0, maxcv_init) == f for every
%     maxcv_init. The default merit function has this identity. A mapping is
%     never a function handle and users cannot register one; an unknown token
%     is rejected, never guessed.
%
%   The fields `fun`, `maxcv` and `point` belong to a superseded record layout.
%   A record that has them is rejected and never reinterpreted as a merit.
%
%   What the record is not. It is not a run-history minimum and not the
%   dynamic cohort minimum of a benchmark: the profile baseline is the least
%   merit observed over the selected solvers, runs and evaluations, changes
%   with the cohort and is never stored in a Problem. It is not a floor for
%   run merits either: on a constrained problem the merit of a run may be
%   BELOW the reference, because a merit function tolerates or penalizes small
%   violations, and such values are never clamped. Before a consumer compares
%   run merits with the record under a custom merit function, that function
%   must be known to preserve the feasible identity above.
%
%   Validation is structural and evaluates nothing, so a reference adds no
%   callback call to building or loading a problem.

    if isempty(reference)
        reference = [];
        return
    end
    fields_required = {'merit', 'kind', 'source', 'mapping'};
    if ~(isstruct(reference) && isscalar(reference))
        error("MATLAB:Problem:reference_NotStruct", ...
            "The field `reference` for `Problem` must be a scalar struct with the fields %s; a naked scalar has no kind, provenance or mapping and is rejected.", ...
            strjoin(fields_required, ', '));
    end
    fields = reshape(fieldnames(reference), 1, []);
    legacy = fields(ismember(fields, {'fun', 'maxcv', 'point'}));
    if ~isempty(legacy)
        error("MATLAB:Problem:reference_LegacyField", ...
            "The problem reference field(s) %s belong to a superseded record layout. Such a record is rejected and never reinterpreted; state the fields %s.", ...
            strjoin(legacy, ', '), strjoin(fields_required, ', '));
    end
    unknown = fields(~ismember(fields, fields_required));
    if ~isempty(unknown)
        error("MATLAB:Problem:reference_UnknownField", ...
            "Unknown problem reference field(s): %s. The fields are %s.", strjoin(unknown, ', '), strjoin(fields_required, ', '));
    end
    missing = fields_required(~ismember(fields_required, fields));
    if ~isempty(missing)
        error("MATLAB:Problem:reference_MissingField", ...
            "The problem reference is missing the required field(s): %s. All of %s are required.", ...
            strjoin(missing, ', '), strjoin(fields_required, ', '));
    end

    merit = referenceMerit(reference.merit);

    kind = referenceText(reference.kind, 'kind');
    kinds = {'lower_bound', 'optimum', 'best_known', 'target'};
    if ~ismember(kind, kinds)
        error("MATLAB:Problem:reference_kind_Unknown", ...
            "The field `kind` of a problem reference must be one of %s, not '%s'.", strjoin(kinds, ', '), kind);
    end

    source = referenceText(reference.source, 'source');
    if isempty(strtrim(source))
        error("MATLAB:Problem:reference_source_Empty", ...
            "The field `source` of a problem reference must be a non-empty provenance string.");
    end

    % The registry is closed: one literal list, no registration function.
    mappings = {'feasible_objective/1'};
    if isa(reference.mapping, 'function_handle')
        error("MATLAB:Problem:reference_mapping_NotToken", ...
            "The field `mapping` of a problem reference must be a token of the closed registry {%s}; function handles and user-defined mappings are not accepted.", ...
            strjoin(mappings, ', '));
    end
    mapping = referenceText(reference.mapping, 'mapping');
    if ~ismember(mapping, mappings)
        error("MATLAB:Problem:reference_mapping_Unknown", ...
            "Unknown problem reference mapping '%s'; the closed registry is {%s}. An unknown mapping is rejected, never guessed.", ...
            mapping, strjoin(mappings, ', '));
    end

    % One canonical field order, whatever order the input used.
    reference = struct('merit', merit, 'kind', kind, 'source', source, 'mapping', mapping);
end

function merit = referenceMerit(value)
    % Logical values are not numeric here, so a logical merit is rejected as in Python.
    if ~(isnumeric(value) && isreal(value) && isscalar(value))
        error("MATLAB:Problem:reference_merit_NotRealScalar", ...
            "The field `merit` of a problem reference must be a real scalar.");
    end
    merit = double(value);
    if ~isfinite(merit)
        error("MATLAB:Problem:reference_merit_NotFinite", ...
            "The field `merit` of a problem reference must be finite.");
    end
    % An integer beyond 2^53 may be stored as a different number; that is a
    % silent reinterpretation, so it is rejected. The test starts AT flintmax,
    % because 2^53 + 1 rounds to exactly 2^53. The round trip alone is not a
    % proof, because cast saturates: double(intmax('int64')) is 2^63, which
    % casts back to intmax. Hence the strict upper limit before the cast.
    if isinteger(value) && abs(merit) >= flintmax ...
            && ~(merit < double(intmax(class(value))) && cast(merit, class(value)) == value)
        error("MATLAB:Problem:reference_merit_NotExact", ...
            "The field `merit` of a problem reference cannot be represented exactly as a double.");
    end
end

function value = referenceText(value, name)
    if isstring(value) && isscalar(value) && ~ismissing(value)
        value = char(value);
    end
    if ~(ischar(value) && (isrow(value) || isempty(value)))
        error("MATLAB:Problem:reference_" + name + "_NotText", ...
            "The field `%s` of a problem reference must be a character vector or a string scalar.", name);
    end
    value = reshape(value, 1, []);
end
