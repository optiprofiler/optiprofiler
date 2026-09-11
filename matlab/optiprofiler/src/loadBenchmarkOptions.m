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
% Native MAT loading is trusted input and can execute user code. Callback names
% in JSON reports are descriptions, not recipes for reconstructing functions.
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
            stages = feature.stages;
            inspection = cell(1, numel(stages));
            for i = 1:numel(stages)
                inspection{i} = struct('name', stages{i}.name, 'options', stages{i}.options);
            end
            if isempty(inspection), inspection = {struct('name', 'plain', 'options', struct())}; end
            if ~isequaln(source.feature_specification, inspection)
                error('OptiProfiler:InconsistentNativeFeature', ...
                    'feature_specification is inspection-only and disagrees with the canonical native Feature.');
            end
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
