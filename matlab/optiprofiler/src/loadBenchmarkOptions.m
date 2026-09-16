function [options, receipt] = loadBenchmarkOptions(source, legacy_feature_name)
%LOADBENCHMARKOPTIONS Read trusted native settings for a fresh benchmark.
% [OPTIONS, RECEIPT] = LOADBENCHMARKOPTIONS(PATH) reads options_refined from a
% trusted MAT file, or accepts the already-loaded scalar native options struct.
% OPTIONS uses exactly one public feature route with n_runs separate. A replay
% runs new solvers; it is not result replotting or resumption of runtime state.
%
% Versioned options_refined-v2 stores effective stages and native callbacks.
% Old flat files must retain feature_name/feature or supply LEGACY_FEATURE_NAME
% explicitly; missing identity is never guessed from paths or option values.
% Options saved by a load (re-plot) invocation (a nonempty load field) are
% rejected: they describe the re-plot, not the archived execution.
% Native MAT loading is trusted input and can execute user code. Callback names
% in JSON reports are descriptions, not recipes for reconstructing functions.
% The inspection duplicate is checked for ordinary values/structure and callback
% positions/types only; callback identity or semantics is not certified. Only
% the canonical native Feature supplies callbacks to the returned replay.
    if nargin < 2, legacy_feature_name = []; end
    if ischar(source) || (isstring(source) && isscalar(source))
        try
            loaded = load(char(source), 'options_refined');
        catch cause
            failure = MException('OptiProfiler:NativeOptionsLoadFailed', ...
                'Unable to load trusted native options. Restore the saved MATLAB callback/class dependencies: %s', cause.message);
            throwAsCaller(addCause(failure, cause));
        end
        if ~isfield(loaded, 'options_refined')
            error('OptiProfiler:MissingNativeOptions', 'The MAT file does not contain options_refined.');
        end
        source = loaded.options_refined;
    end
    if ~isstruct(source) || ~isscalar(source)
        error('OptiProfiler:InvalidNativeOptions', 'Expected a scalar native options struct or a trusted options_refined MAT file.');
    end
    if isfield(source, 'load') && ~isempty(source.load) && ...
            ~(isstring(source.load) && isscalar(source.load) && strlength(source.load) == 0)
        % Options recorded by benchmark(load=...) describe that re-plot, not
        % the archived execution: the label and its defaults are not evidence
        % of the archived experiment, so they are never turned into a replay.
        error('OptiProfiler:LoadInvocationNotReplayable', ...
            ['These options were written by a load (re-plot) invocation and describe that load, not the archived ', ...
             'execution. Replay the source experiment from its own test_log/options_refined.mat.']);
    end
    receipt = struct('schema', 'matlab-benchmark-options-import-v1', ...
        'source_schema', 'legacy-flat-options', 'original_request_origin', 'unknown', ...
        'scope', 'fresh_replay_not_original_execution_or_live_resume');
    if isfield(source, 'schema')
        if ~(ischar(source.schema) || (isstring(source.schema) && isscalar(source.schema))) ...
                || ~strcmp(source.schema, 'options_refined-v2')
            error('OptiProfiler:UnsupportedNativeOptionsVersion', 'Only options_refined-v2 and supported unversioned native options are understood.');
        end
        if ~isfield(source, 'feature_specification') || ~isfield(source, 'n_runs')
            error('OptiProfiler:IncompleteNativeOptions', 'options_refined-v2 requires feature_specification and separate n_runs.');
        end
        if isfield(source, 'feature')
            if ~isa(source.feature, 'Feature') || ~isscalar(source.feature)
                error('OptiProfiler:InvalidNativeFeature', 'The retained v2 feature must be a genuine scalar Feature object.');
            end
            feature = source.feature;
            try
                stages = feature.stages;
            catch cause
                % MATLAB can warn on loadobj failure and return an object of
                % the right class with unusable state. Class identity alone
                % therefore cannot certify a canonical replay object.
                failure = MException('OptiProfiler:InvalidNativeFeature', ...
                    'The retained Feature could not restore its canonical state. Check the MAT load warnings and restore compatible Feature/callback dependencies: %s', cause.message);
                throwAsCaller(addCause(failure, cause));
            end
            inspection = cell(1, numel(stages));
            for i = 1:numel(stages)
                inspection{i} = struct('name', stages{i}.name, 'options', stages{i}.options);
            end
            if isempty(inspection), inspection = {struct('name', 'plain', 'options', struct())}; end
            if ~sameInspection(source.feature_specification, inspection)
                error('OptiProfiler:InconsistentNativeFeature', ...
                    'feature_specification is inspection-only and disagrees with the canonical native Feature.');
            end
            receipt.feature_specification_comparison = 'ordinary_structure_and_values_callback_positions_only';
            receipt.callback_comparison = 'opaque_native_handles_not_identity_or_semantics';
        else
            % A specification-only v2 input is an explicit new replay request;
            % it does not recover a historical declaration or invocation route.
            feature = Feature(source.feature_specification);
        end
        receipt.source_schema = 'options_refined-v2';
        if isfield(source, 'feature_provenance'), receipt.retained_feature_provenance = source.feature_provenance; end
    elseif isfield(source, 'feature')
        if isa(source.feature, 'optiprofiler_internal.LegacyFeatureEnvelope')
            [feature, request] = optiprofiler_internal.importLegacyFeature(source.feature);
            if ~isfield(source, 'n_runs') && request.n_runs_present, source.n_runs = request.n_runs; end
            receipt.feature_import = request;
        else
            feature = Feature(source.feature);
        end
    else
        if isfield(source, 'feature_name')
            name = source.feature_name;
        elseif ~isempty(legacy_feature_name)
            name = legacy_feature_name;
        else
            error('OptiProfiler:LegacyFeatureIdentityMissing', ...
                'This old flat options file did not retain feature identity. Supply legacy_feature_name explicitly; stamps and folders are not evidence.');
        end
        definitions = optiprofiler_internal.featureDefinitions();
        local_keys = {};
        for i = 1:numel(definitions), local_keys = [local_keys, definitions(i).local_keys]; end %#ok<AGROW>
        local_options = struct();
        keys = intersect(fieldnames(source), unique(local_keys));
        for i = 1:numel(keys), local_options.(keys{i}) = source.(keys{i}); end
        saved = struct('name', name, 'options', local_options);
        if isfield(source, 'n_runs'), saved.options.n_runs = source.n_runs; end
        [feature, request] = optiprofiler_internal.importLegacyFeature(optiprofiler_internal.LegacyFeatureEnvelope(saved));
        receipt.feature_import = request;
        if ~isempty(legacy_feature_name) && ~isfield(source, 'feature_name')
            receipt.current_legacy_feature_name_override = legacy_feature_name;
        end
    end
    % Never forward envelope facts or saved output targets to benchmark.
    % The declared solver list is still supplied by the caller separately.
    profile_keys = arrayfun(@(x) x.value, enumeration('ProfileOptionKey'), 'UniformOutput', false);
    problem_keys = arrayfun(@(x) x.value, enumeration('ProblemOptionKey'), 'UniformOutput', false);
    allowed = setdiff([profile_keys(:); problem_keys(:)], ...
        {'load'; 'solvers_to_load'; 'report_path'; 'savepath'; 'benchmark_id'; 'feature'; 'feature_name'});
    options = struct();
    keys = intersect(fieldnames(source), allowed);
    for i = 1:numel(keys), options.(keys{i}) = source.(keys{i}); end
    options.feature = feature;
    if isfield(source, 'n_runs')
        options.n_runs = source.n_runs;
        receipt.retained_n_runs = source.n_runs;
    end
end

function same = sameInspection(left, right)
% Anonymous/closure handles restored from the two native copies need not be
% isequal even when saved together. Compare ordinary fields exactly, but only
% require callback positions/types to agree. Never inspect a callback closure
% or execute it to guess semantic equality; only canonical Feature is replayed.
    if isa(left, 'function_handle') || isa(right, 'function_handle')
        same = isa(left, 'function_handle') && isa(right, 'function_handle') ...
            && isequal(size(left), size(right));
        return;
    end
    if ~strcmp(class(left), class(right)) || ~isequal(size(left), size(right))
        same = false; return;
    end
    if iscell(left)
        same = true;
        for k = 1:numel(left)
            if ~sameInspection(left{k}, right{k}), same = false; return; end
        end
    elseif isstruct(left)
        keys = fieldnames(left);
        same = isequal(sort(keys), sort(fieldnames(right)));
        if ~same, return; end
        for k = 1:numel(left)
            for j = 1:numel(keys)
                if ~sameInspection(left(k).(keys{j}), right(k).(keys{j}))
                    same = false; return;
                end
            end
        end
    elseif isnumeric(left) || islogical(left) || ischar(left) || isstring(left)
        same = isequaln(left, right);
    else
        % Feature options have ordinary data or callback leaves. Unsupported
        % objects must not run overloaded equality during metadata checking.
        same = false;
    end
end
