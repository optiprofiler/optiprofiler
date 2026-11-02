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
% Unversioned options with a nonempty load field (the options_user of a load,
% or a flat options_refined that older versions wrote for a load) are rejected:
% they describe the re-plot and their feature label is not evidence of the
% archived execution. A versioned options_refined-v2 record carries its own
% effective feature identity (benchmark never writes one for a load), and its
% obsolete load/solvers_to_load selectors are simply not forwarded.
% options_refined does not store a directly supplied Problem; its problem_route
% field ('direct_problem' or 'library_selection') says whether the experiment
% used one. A direct-problem record raises OptiProfiler:DirectProblemNotSaved
% instead of becoming a library replay, unless the original Problem is set as
% the problem field of the struct passed here. A source with a problem field
% (such as options_user of a direct-problem run) replays that Problem, without
% library selectors. RECEIPT.problem_route records the route: 'direct_problem',
% 'library_selection', or 'not_recorded' for a record older than the field.
% Such a record has default library selectors whatever the experiment used. For
% a file, the options_user.mat next to it settles the route. Otherwise a replay
% warns with OptiProfiler:ProblemRouteNotRecorded, because a direct-problem
% experiment would run on those selectors.
% Native MAT loading is trusted input and can execute user code. Callback names
% in JSON reports are descriptions, not recipes for reconstructing functions.
% The inspection duplicate is checked for ordinary values/structure and callback
% positions/types only; callback identity or semantics is not certified. Only
% the canonical native Feature supplies callbacks to the returned replay.
    if nargin < 2, legacy_feature_name = []; end
    source_file = '';
    if ischar(source) || (isstring(source) && isscalar(source))
        source_file = char(source);
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
    has_load = isfield(source, 'load') && ~isempty(source.load) && ...
        ~(isstring(source.load) && isscalar(source.load) && strlength(source.load) == 0);
    if has_load && ~isfield(source, 'schema')
        % Unversioned options recorded by benchmark(load=...) describe that
        % re-plot, not the archived execution: the label and its defaults are
        % not evidence of the archived experiment, so they are never turned
        % into a replay. (A versioned record is validated below; it carries
        % its own effective feature identity and drops the selector.)
        error('OptiProfiler:LoadInvocationNotReplayable', ...
            ['These options were written by a load (re-plot) invocation and describe that load, not the archived ', ...
             'execution. Replay the source experiment from its own test_log/options_refined.mat.']);
    end
    receipt = struct('schema', 'matlab-benchmark-options-import-v1', ...
        'source_schema', 'legacy-flat-options', 'original_request_origin', 'unknown', ...
        'scope', 'fresh_replay_not_original_execution_or_live_resume');
    has_problem = isfield(source, 'problem');
    if has_problem && ~(isa(source.problem, 'Problem') && isscalar(source.problem))
        error('OptiProfiler:InvalidNativeProblem', 'The problem field of the options must be a scalar Problem.');
    end
    receipt.problem_route = 'not_recorded';
    if isfield(source, 'problem_route')
        if ~(ischar(source.problem_route) || (isstring(source.problem_route) && isscalar(source.problem_route))) ...
                || ~ismember(char(source.problem_route), {'direct_problem', 'library_selection'})
            error('OptiProfiler:InvalidNativeOptions', 'problem_route must be ''direct_problem'' or ''library_selection''.');
        end
        receipt.problem_route = char(source.problem_route);
    elseif ~has_problem
        receipt.problem_route = routeOfUserOptions(source_file);
        if strcmp(receipt.problem_route, 'not_recorded')
            warning('OptiProfiler:ProblemRouteNotRecorded', ...
                ['These options do not record whether the experiment selected its problems from a library or used a ', ...
                 'directly supplied Problem. If it used a Problem, this replay would benchmark the recorded library ', ...
                 'selectors instead: set the original Problem as the problem field of the options.']);
        end
    end
    if strcmp(receipt.problem_route, 'direct_problem') && ~has_problem
        error('OptiProfiler:DirectProblemNotSaved', ...
            ['These options come from a benchmark of a directly supplied Problem, which options_refined does not store. ', ...
             'A replay without it would benchmark a problem library instead. Set the original Problem as the problem ', ...
             'field of the loaded struct, for example S = load(file); S.options_refined.problem = p; ', ...
             'options = loadBenchmarkOptions(S.options_refined).']);
    end
    if has_problem, receipt.problem_route = 'direct_problem'; end
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
    else
        [feature, legacy_receipt] = optiprofiler_internal.importLegacyBenchmarkOptions(source, legacy_feature_name);
        keys = fieldnames(legacy_receipt);
        for i = 1:numel(keys), receipt.(keys{i}) = legacy_receipt.(keys{i}); end
        if isfield(legacy_receipt, 'feature_import')
            request = legacy_receipt.feature_import;
            if ~isfield(source, 'n_runs') && request.n_runs_present, source.n_runs = request.n_runs; end
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
    if has_problem, keys = setdiff(keys, problem_keys); end
    for i = 1:numel(keys), options.(keys{i}) = source.(keys{i}); end
    % A replay of a direct-problem experiment benchmarks the same Problem.
    if has_problem, options.problem = source.problem; end
    options.feature = feature;
    if isfield(source, 'n_runs')
        options.n_runs = source.n_runs;
        receipt.retained_n_runs = source.n_runs;
    end
end

function route = routeOfUserOptions(source_file)
% options_user.mat beside a record holds what the caller passed to benchmark,
% with a Problem if there was one, in every version that wrote such a record.
    route = 'not_recorded';
    % load also takes the path without its extension.
    if isempty(source_file) || ~(isfile(source_file) || isfile([source_file, '.mat'])), return; end
    sibling = fullfile(fileparts(source_file), 'options_user.mat');
    if ~isfile(sibling), return; end
    try
        loaded = load(sibling, 'options_user');
    catch
        return
    end
    if ~isfield(loaded, 'options_user') || ~isstruct(loaded.options_user) || ~isscalar(loaded.options_user)
        return
    end
    if ~isfield(loaded.options_user, 'problem')
        route = 'library_selection';
    elseif isa(loaded.options_user.problem, 'Problem')
        route = 'direct_problem';
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
