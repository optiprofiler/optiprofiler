classdef EvalReport < handle
%EVALREPORT Controller-private, observational JSON collector. Not a public API.
% Numeric completion is not convergence; recorded history excludes padding.
% The emitted vocabulary is shared with the Python collector
% (python/optiprofiler/eval_report.py) and pinned by
% python/optiprofiler/schemas/eval_report.schema.json and plot_data.schema.json
% (Python package resources, also downloadable from the documentation):
% identical concepts use identical keys; genuine MATLAB conventions (sample
% standard deviation, 1-based seed rule, retained log-ratio bar identity,
% native/portable renderer variants) are explicit fields, not renamed ones.
    properties (Access = private)
        path
        document
        started
        outputDirectory = ''
        maxEvalFactor = NaN
        replaceFile
        renderingFailure = false
        renderingFailureCode = ''
        renderingSuccess = false
        sourcePath = ''
        convergence = {}
        ownedIdentity
        outputIdentity = ''
        plotPath
        plotIdentity
        plotDocument
        historyPreparation
        profileOptions
        capturedRenderIds = {}
        % id -> position maps keep problem/history/plot upserts O(1) while
        % the documents themselves stay ordered cell arrays (insertion order
        % is the selection/merge order readers rely on).
        problemPositions
        plotPositions
        historyPositions
        permissionsApplied
    end
    methods
        function self = EvalReport(path, request, replaceFile)
            if ~(ischar(path) && isrow(path)) && ~(isstring(path) && isscalar(path))
                error('OptiProfiler:EvalReportPath', 'report_path must be a scalar char/string path.');
            end
            self.path = char(path);
            if ~optiprofiler_internal.EvalReport.isAbsolute(self.path), self.path = fullfile(pwd, self.path); end
            self.replaceFile = replaceFile;
            parent = fileparts(self.path);
            if isempty(parent), self.path = fullfile(pwd, self.path); parent = pwd; end
            [~, stem] = fileparts(self.path);
            self.plotPath = fullfile(parent, [stem, '.plot_data.json']);
            if ~isfolder(parent), mkdir(parent); end
            if ~optiprofiler_internal.EvalReport.reserve(self.path)
                error('OptiProfiler:EvalReportExists', 'The report target already exists or cannot be reserved.');
            end
            self.ownedIdentity = optiprofiler_internal.EvalReport.fileIdentity(self.path);
            % Both paths are owned before any solver call. On a companion
            % collision, remove only our still-identical empty main shell.
            if ~optiprofiler_internal.EvalReport.reserve(self.plotPath)
                self.checkOwnership(); delete(self.path);
                error('OptiProfiler:EvalReportExists', 'The plot-data target already exists or cannot be reserved.');
            end
            self.plotIdentity = optiprofiler_internal.EvalReport.fileIdentity(self.plotPath);
            self.problemPositions = containers.Map('KeyType', 'char', 'ValueType', 'double');
            self.plotPositions = containers.Map('KeyType', 'char', 'ValueType', 'double');
            self.historyPositions = containers.Map('KeyType', 'char', 'ValueType', 'double');
            self.permissionsApplied = optiprofiler_internal.EvalReport.null();
            self.started = tic;
            null = optiprofiler_internal.EvalReport.null();
            stage = struct('status', 'unknown');
            self.document = struct('schema', 'optiprofiler.eval_report/1', ...
                'evaluation_id', optiprofiler_internal.EvalReport.identifier(), ...
                'operation', 'benchmark', 'status', 'running', ...
                'producer', struct('language', 'matlab', 'version', null, 'revision', null, ...
                    'revision_reason', 'not_available', 'scope', 'current_invocation_not_original_archive_producer'), ...
                'configuration', struct(), 'stages', struct('numerical', stage, 'scoring', stage, 'persistence', stage, 'rendering', stage), ...
                'coverage', struct('selected', 0, 'loaded', 0, 'completed', 0, 'load_failed', 0, ...
                    'scope', 'live_selection', 'original_selection_known', true), ...
                'problems', {{}}, 'scores', null, ...
                'profiles', struct('plot_refs', {{}}, 'convergence', {{}}, ...
                    'work_summary', null, 'work_summary_reason', 'not_supplied'), ...
                'artifact_root', null, 'artifacts', {{}}, 'diagnostics', {{}}, ...
                'plot_data', null, 'source', null, 'report_files', null, ...
                'timing', struct('started_at', optiprofiler_internal.EvalReport.timestamp(), 'finished_at', null, 'elapsed_seconds', null));
            % Shared main-report semantics keys (see eval_report.schema.json);
            % the prose describes THIS producer.
            self.document.semantics = struct('index_base', 1, ...
                'metric_best', 'componentwise_minimum_ignoring_nan;first_tie_index;not_necessarily_a_jointly_attained_point', ...
                'budget', 'problems[].budget applies to every run unless runs[].budget overrides it; budget_reached is a per-run comparison and reaching the cap does not identify the termination cause.', ...
                'convergence', 'Never inferred from solver return values; there is no per-run convergence field. target_work in the companion observes the existing profile construction.', ...
                'run_defaults', 'evaluations, budget_reached, abnormal_termination, output_fallback, execution, oracle_seed and elapsed_seconds are per-run facts. An absent *_reason key means the observation was available; an absent first_invalid_evaluation_index means no invalid evaluation was observed and invalid_evaluations is then the integer 0.', ...
                'configuration', 'configuration.request lists user-supplied options; configuration.effective lists the resolved options of this invocation. In a load operation they describe reanalysis/rendering, not the archived execution.', ...
                'paths', struct('artifacts', 'relative_to_artifact_root', 'artifact_root', 'relative_to_main_report_parent', ...
                    'plot_data', 'relative_to_main_report_parent', 'source', 'relative_to_main_report_parent'), ...
                'privacy', 'controller_private;not_an_allowlisted_agent_prompt;consumers_build_an_allowlisted_feedback_view');
            self.plotDocument = struct('schema', 'optiprofiler.plot_data/1', ...
                'evaluation_id', self.document.evaluation_id, ...
                'semantics', struct('index_base', 1, ...
                    'history_bins', 'exact_samples_up_to_64;otherwise_lossy_max_32_contiguous_bins_with_exact_endpoints_finite_extrema_first_tie_indices_and_nonfinite_counts;internal_order_not_retained_in_bins;excludes_padding;not_a_scoring_input', ...
                    'history_plots', 'display_copy_only;across_run_band_then_shift_then_cummin_then_block_aggregation;block_aggregation_applies_when_padded_length_exceeds_1002_and_keeps_about_1000_points;array_position_is_not_an_evaluation_index_use_evaluation_indices;resolve_renderer_variants_ref_before_interpreting_native_versus_portable_styling', ...
                    'profile_plots', 'full_step_vertices_and_across_run_bands;band_clamped_to_0_1;unreached_work_shown_at_failure_placeholder_not_counted_as_hit;log_ratio_zero_height_ties_retained_in_numeric_data_but_not_drawn;log_ratio_bar_sources_retain_problem_and_run_identity', ...
                    'error_bands', 'matlab_sample_std_ddof_1_for_more_than_one_run_and_0_for_a_single_run', ...
                    'target_work', 'existing_profile_work_arrays;history_first_actual_hit;output_total_evaluations_at_passing_returned_output;nan_is_nonhit_not_solver_diagnosis', ...
                    'privacy', 'controller_private;not_an_allowlisted_agent_prompt'), ...
                'histories', {{}}, 'plots', {{}}, 'target_work', {{}}, 'diagnostics', {{}});
            % Shared variant descriptions avoid repeating the same display
            % policy on every panel. They describe existing renderers; they
            % do not claim that either renderer succeeded in this invocation.
            self.plotDocument.semantics.renderer_variants = struct( ...
                'history', struct( ...
                    'native', struct('geometry','mean line or singleton point; filled lower/upper band when series.band_visible', ...
                        'coordinates','provided x/mean/lower/upper; logarithmic Y axis when panel.y_scale=log'), ...
                    'portable_svg', struct('geometry','mean line, plus separate lower/upper lines when n_runs>1; no filled band', ...
                        'coordinates','provided x; when panel.y_scale=log, nonpositive Y is omitted and log10(Y) is drawn on linear axes')), ...
                'performance_data', struct( ...
                    'native', struct('geometry','provided step mean with filled lower/upper band when series.band_visible'), ...
                    'portable_svg', struct('geometry','provided step mean only; lower/upper bands are not drawn')), ...
                'log_ratio', struct( ...
                    'native', struct('geometry','provided visible bars, including half-opacity both-failed placeholders'), ...
                    'portable_svg', struct('geometry','two line groups connecting negative and positive y separately; zero ties omitted; no bars or half-opacity styling')));
            if isfield(request, 'load') && ~isempty(request.load)
                self.document.operation = 'load';
                self.document.coverage.scope = 'retained_archive';
                self.document.coverage.original_selection_known = false;
            end
            self.document.configuration.request = optiprofiler_internal.EvalReport.configuration(request);
            self.setStage('numerical', 'running');
            self.write();
        end

        function configure(self, problem_options, profile_options, feature_value, historyPreparation)
            self.historyPreparation = historyPreparation;
            self.profileOptions = profile_options;
            % request = what the caller supplied (kept from the constructor);
            % effective = resolved options, stated once, never per run.
            feature = struct('name', feature_value.name, 'options', optiprofiler_internal.EvalReport.configuration(feature_value.options));
            self.document.configuration.effective = struct( ...
                'problem_options', optiprofiler_internal.EvalReport.configuration(problem_options), ...
                'profile_options', optiprofiler_internal.EvalReport.configuration(profile_options), ...
                'feature', feature);
            self.maxEvalFactor = profile_options.max_eval_factor;
            if strcmp(self.document.operation, 'load')
                self.document.configuration.scope = 'current_load_selection_reanalysis_and_rendering';
                self.document.configuration.solver_execution_requested = false;
                self.document.configuration.original_execution_configuration = optiprofiler_internal.EvalReport.null();
                self.document.configuration.original_execution_configuration_reason = 'not_fully_retained_by_archive';
                self.document.configuration.effective.feature = struct('name', optiprofiler_internal.EvalReport.null(), ...
                    'options', optiprofiler_internal.EvalReport.null(), 'scope', 'current_load_context_not_original_execution_feature', ...
                    'reason', 'default_load_feature_is_not_original_execution_provenance');
            else
                self.document.configuration.scope = 'current_benchmark_execution';
                self.document.configuration.solver_execution_requested = true;
                self.document.configuration.effective.feature.scope = 'current_execution_feature';
            end
            if profile_options.score_only
                self.setStage('rendering', 'not_requested', 'score_only');
                self.setStage('persistence', 'not_requested', 'score_only');
            else
                self.setStage('rendering', 'unknown', 'requested_not_yet_observed');
                self.setStage('persistence', 'unknown', 'requested_not_yet_observed');
            end
            self.write();
        end

        function setOutputDirectory(self, path)
            self.outputDirectory = path;
            self.outputIdentity = optiprofiler_internal.EvalReport.fileIdentity(path);
            self.document.artifact_root = self.relative(path);
        end

        function setStage(self, stage, status, reason)
            value = struct('status', status);
            if nargin > 3, value.reason = reason; end
            self.document.stages.(stage) = value;
        end

        function addDiagnostic(self, code, stage, scope)
            if nargin < 4, scope = struct(); end
            if numel(self.document.diagnostics) < 128
                self.document.diagnostics{end+1} = struct('code', code, 'stage', stage, 'scope', scope);
            end
        end

        function addProblem(self, result, library, role)
            if nargin < 4, role = 'primary'; end
            id = jsonencode({library, result.problem_name, role});
            ordinal = self.problemIndex(id);
            if isempty(ordinal), ordinal = numel(self.document.problems) + 1; end
            unknown = optiprofiler_internal.EvalReport.null();
            metadata = struct();
            if isfield(result, 'eval_report_metadata'), metadata = result.eval_report_metadata; end
            % The cap is a per-problem fact of this invocation. Load cannot
            % infer the archived budget from today's max_eval_factor.
            if strcmp(self.document.operation, 'load')
                budget = unknown;
                budget_record = struct('evaluations', unknown, 'reason', 'original_execution_budget_not_retained');
            else
                budget = ceil(self.maxEvalFactor * result.problem_dim);
                budget_record = struct('evaluations', budget, 'rule', 'ceil(max_eval_factor*dimension)');
            end
            entry = struct('id', id, ...
                'library', library, 'name', result.problem_name, 'role', role, ...
                'dimension', result.problem_dim, 'type', result.problem_type, ...
                'selection_status', 'selected', 'load_status', 'loaded', 'status', 'completed', ...
                'budget', budget_record, 'runs', {{}}, 'plot_refs', {{}});
            for solver = 1:size(result.n_eval, 1)
                for run = 1:size(result.n_eval, 2)
                    count = result.n_eval(solver, run);
                    retained = min(count, size(result.fun_history, 3));
                    values = reshape(result.fun_history(solver, run, 1:retained), 1, []);
                    constraints = reshape(result.maxcv_history(solver, run, 1:min(count,size(result.maxcv_history,3))), 1, []);
                    merits = [];
                    if isfield(result, 'merit_history'), merits = reshape(result.merit_history(solver, run, 1:min(count,size(result.merit_history,3))), 1, []); end
                    budget_reached = unknown;
                    if ~isstruct(budget), budget_reached = count >= budget; end
                    elapsed = unknown;
                    if isfield(result, 'computation_time'), elapsed = result.computation_time(solver, run); end
                    item = struct('solver_index', solver, 'run_index', run, 'evaluations', count, ...
                        'budget_reached', budget_reached, ...
                        'objective', optiprofiler_internal.EvalReport.metric(values, count, result.fun_out(solver, run), result.fun_inits(run)), ...
                        'constraint', optiprofiler_internal.EvalReport.metric(constraints, count, result.maxcv_out(solver, run), result.maxcv_inits(run)), ...
                        'merit', optiprofiler_internal.EvalReport.metric([], unknown, unknown, unknown), ...
                        'abnormal_termination', result.solver_abnormal_termination(solver, run), ...
                        'output_fallback', result.solver_output_fallback(solver, run), ...
                        'execution', struct('kind', 'unknown', 'source_run_index', unknown, 'reason', 'execution_metadata_not_retained'), ...
                        'history_ref', sprintf('history-%d-%d-%d', ordinal, solver, run), ...
                        'oracle_seed', unknown, 'elapsed_seconds', elapsed);
                    if isfield(metadata, 'real_n_runs')
                        item.execution = struct('kind', 'actual', 'source_run_index', unknown);
                        if run > metadata.real_n_runs(solver)
                            item.execution = struct('kind', 'repeated', 'source_run_index', 1);
                        end
                    end
                    if isfield(result, 'merit_history')
                        merit_out = unknown;
                        if isfield(result, 'merit_out'), merit_out = result.merit_out(solver, run); end
                        item.merit = optiprofiler_internal.EvalReport.metric(merits, count, merit_out, result.merit_inits(run));
                    end
                    if isstruct(item.abnormal_termination) || isstruct(item.output_fallback)
                        item.termination_metadata_reason = 'solver_termination_metadata_not_retained';
                    end
                    % Actual seed of the featured problem under MATLAB's own
                    % 1-based rule; a repeated slot copies run 1's seed.
                    if isfield(metadata, 'oracle_seeds')
                        if strcmp(item.execution.kind, 'repeated')
                            item.oracle_seed = metadata.oracle_seeds(1);
                            item.oracle_seed_reason = 'copied_from_source_run';
                        else
                            item.oracle_seed = metadata.oracle_seeds(run);
                        end
                    else
                        item.oracle_seed_reason = 'execution_metadata_not_retained';
                    end
                    entry.runs{end+1} = item;
                    channels = struct('objective', optiprofiler_internal.EvalReport.binned(values, count, true), ...
                        'constraint', optiprofiler_internal.EvalReport.binned(constraints, count, true), ...
                        'merit', optiprofiler_internal.EvalReport.binned(merits, count, isfield(result, 'merit_history')));
                    self.upsert('histories', struct('id', item.history_ref, 'problem_id', id, ...
                        'solver_index', solver, 'run_index', run, 'channels', channels));
                end
            end
            index = self.problemIndex(entry.id);
            if strcmp(self.document.operation, 'load'), entry.selection_status = 'retained'; end
            if isfield(metadata, 'render_status')
                entry.rendering = struct('status', metadata.render_status);
                if isempty(index) || isempty(self.document.problems{index}.runs)
                    self.recordRendering(metadata.render_status, struct('library', library, 'problem', result.problem_name));
                end
            end
            if isempty(index)
                self.document.problems{end+1} = entry;
                self.problemPositions(entry.id) = numel(self.document.problems);
            else
                old = self.document.problems{index};
                if isfield(old, 'plot_refs'), entry.plot_refs = old.plot_refs; end
                if isfield(old, 'provider'), entry.provider = old.provider; end
                if isfield(old, 'rendering') && ~isfield(entry, 'rendering'), entry.rendering = old.rendering; end
                if ~isfield(result, 'eval_report_metadata') && numel(old.runs) == numel(entry.runs)
                    for k = 1:numel(entry.runs)
                        entry.runs{k}.execution = old.runs{k}.execution;
                        entry.runs{k}.oracle_seed = old.runs{k}.oracle_seed;
                        if isfield(old.runs{k}, 'oracle_seed_reason')
                            entry.runs{k}.oracle_seed_reason = old.runs{k}.oracle_seed_reason;
                        elseif isfield(entry.runs{k}, 'oracle_seed_reason')
                            entry.runs{k} = rmfield(entry.runs{k}, 'oracle_seed_reason');
                        end
                    end
                end
                self.document.problems{index} = entry;
            end
            self.recount();
            % Shared preparation consumes retained arrays only. It never
            % evaluates an oracle, feature, solver or custom merit callback.
            if isfield(metadata, 'plot_presentation')
                self.captureHistoryPresentation(metadata.plot_presentation, library, result.problem_name, role);
            end
            if ~isempty(self.historyPreparation) && ~ismember(id, self.capturedRenderIds)
                try
                    panels = self.historyPreparation(result, self.profileOptions);
                    self.storeHistoryPanels(panels, id, ordinal, 'retained_scoring_observations');
                catch cause
                    self.plotDocument.diagnostics{end+1} = struct('code', 'history_presentation_unavailable', ...
                        'problem_id', id, 'exception_type', cause.identifier);
                end
            end
        end

        function captureHistoryPresentation(self, panels, library, name, role)
            if isempty(panels), return; end
            id = jsonencode({library,name,role}); ordinal = self.problemIndex(id);
            if isempty(ordinal), return; end
            if isfield(panels{1}, 'preparation_error')
                self.plotDocument.diagnostics{end+1} = struct('code', 'render_input_presentation_unavailable', ...
                    'problem_id', id, 'exception_type', panels{1}.preparation_error);
                return;
            end
            self.storeHistoryPanels(panels,id,ordinal,'rendering_inputs');
            if ~ismember(id,self.capturedRenderIds), self.capturedRenderIds{end+1} = id; end
        end

        function selection(self, library, names, role, provider)
            for k = 1:numel(names)
                id = jsonencode({library, names{k}, role});
                if isempty(self.problemIndex(id))
                    entry = struct('id', id, 'library', library, 'name', names{k}, 'role', role, ...
                        'dimension', optiprofiler_internal.EvalReport.null(), 'type', optiprofiler_internal.EvalReport.null(), ...
                        'selection_status', 'selected', 'load_status', 'pending', 'status', 'pending', ...
                        'budget', optiprofiler_internal.EvalReport.null(), 'runs', {{}}, 'plot_refs', {{}});
                    if nargin > 4 && ~isempty(provider), entry.provider = optiprofiler_internal.EvalReport.configuration(provider); end
                    self.document.problems{end+1} = entry;
                    self.problemPositions(id) = numel(self.document.problems);
                end
            end
            self.recount();
        end

        function loadFailed(self, library, name, role, exception_type)
            index = self.problemIndex(jsonencode({library, name, role}));
            self.document.problems{index}.load_status = 'failed';
            self.document.problems{index}.status = 'failed';
            scope = struct('library', library, 'problem', name, 'role', role);
            if nargin > 4 && ~isempty(exception_type), scope.exception_type = exception_type; end
            self.addDiagnostic('problem_load_failed', 'numerical', scope);
            self.recount();
        end

        function completeNumerical(self)
            if self.document.coverage.load_failed > 0 && self.document.coverage.loaded == 0
                self.setStage('numerical', 'failed', 'all_selected_problem_loads_failed');
            elseif self.document.coverage.load_failed > 0
                self.setStage('numerical', 'partial', 'selected_problem_load_failures');
            elseif self.document.coverage.loaded > 0
                self.setStage('numerical', 'completed');
            else
                self.setStage('numerical', 'not_applicable', 'empty_selection');
            end
            self.write();
        end

        function addResults(self, groups, role)
            if nargin < 3, role = 'primary'; end
            for group = 1:numel(groups)
                value = groups{group};
                for p = 1:numel(value.problem_names)
                    ns = size(value.fun_histories, 2); nr = size(value.fun_histories, 3);
                    item = struct('problem_name', value.problem_names{p}, 'problem_dim', value.problem_dims(p), ...
                        'problem_type', value.problem_types{p});
                    mappings = {'fun_histories','fun_history'; 'maxcv_histories','maxcv_history'; ...
                        'merit_histories','merit_history'; 'fun_outs','fun_out'; 'maxcv_outs','maxcv_out'; ...
                        'merit_outs','merit_out'; 'n_evals','n_eval'; 'computation_times','computation_time'; ...
                        'solver_abnormal_terminations','solver_abnormal_termination'; ...
                        'solver_output_fallbacks','solver_output_fallback'};
                    for k = 1:size(mappings, 1)
                        from = mappings{k, 1}; to = mappings{k, 2};
                        if isfield(value, from)
                            if contains(from, 'histories'), item.(to) = reshape(value.(from)(p,:,:,:), ns, nr, []);
                            else, item.(to) = reshape(value.(from)(p,:,:), ns, nr); end
                        end
                    end
                    for key = {'fun_inits', 'maxcv_inits', 'merit_inits'}
                        if isfield(value, key{1}), item.(key{1}) = value.(key{1})(p,:); end
                    end
                    for key = {'solver_abnormal_termination', 'solver_output_fallback'}
                        if ~isfield(item, key{1}), item.(key{1}) = repmat(optiprofiler_internal.EvalReport.null(), ns, nr); end
                    end
                    self.addProblem(item, value.plib, role);
                end
                if isfield(value, 'results_plib_plain'), self.addResults({value.results_plib_plain}, 'plain_reference'); end
            end
        end

        function setProfiles(self, curves, solver_scores, profile_scores, solver_names, single)
            semantics = 'cohort_relative';
            if single, semantics = 'single_problem_relative_decrease_averaged_over_runs'; end
            types = {'performance', 'data', 'log_ratio'};
            n_types = 0;
            if ~isempty(profile_scores), n_types = size(profile_scores, 4); end
            self.document.scores = struct('solver_names', {optiprofiler_internal.EvalReport.configuration(solver_names)}, 'solver_scores', {num2cell(solver_scores)}, ...
                'profile_scores', {optiprofiler_internal.EvalReport.tensor(profile_scores)}, 'profile_axes', {{'solver', 'tolerance', 'history_or_output', 'profile_type'}}, ...
                'axis_values', struct('history_or_output', {{'history', 'output'}}, 'profile_type', {types(1:n_types)}, ...
                    'tolerance', {num2cell(10.^(-(1:numel(curves))))}), ...
                'semantics', semantics, 'direction', 'higher_is_better_for_default_scoring; custom_callbacks_define_their_own_direction', ...
                'comparability_note', 'Scores depend on the retained cohort, tolerances, options, and scoring callbacks.');
            self.setStage('scoring', 'completed');
        end

        function addProfilePresentation(self, presentation, tolerance_index, channel)
            for k = 1:numel(presentation)
                item = presentation{k};
                if isfield(item, 'preparation_error')
                    self.plotDocument.diagnostics{end+1} = struct('code', 'profile_presentation_unavailable', ...
                        'tolerance_index', tolerance_index, 'history_or_output', channel, 'exception_type', item.preparation_error);
                    continue;
                end
                item.id = sprintf('profile-%d-%s-%s', tolerance_index, channel, item.kind);
                item.tolerance_index = tolerance_index;
                item.history_or_output = channel;
                item.observation_scope = 'scoring_profile_work';
                item.target_work_ref = sprintf('target-work-%d', tolerance_index);
                item.renderer_variants_ref = 'performance_data';
                if strcmp(item.kind,'log_ratio'), item.renderer_variants_ref = 'log_ratio'; end
                self.upsert('plots', item);
                if ~ismember(item.id, self.document.profiles.plot_refs)
                    self.document.profiles.plot_refs{end+1} = item.id;
                end
            end
        end

        function addConvergence(self, tolerance, history, output, problem_ids)
            entry = struct('tolerance_index', numel(self.convergence)+1, 'tolerance', tolerance, ...
                'semantics', 'observed_work_used_by_existing_profile_construction', ...
                'denominator', 'all_retained_problem_run_pairs_per_solver', ...
                'invalid_initial_separation', optiprofiler_internal.EvalReport.null(), ...
                'invalid_initial_reason', 'not_separately_supplied', 'problem_count', size(history, 1), ...
                'history_work_semantics', 'first_actual_evaluation_meeting_existing_target', ...
                'output_work_semantics', 'total_evaluations_at_passing_returned_output_not_first_hit', ...
                'history', {{}}, 'output', {{}});
            arrays = {history, output}; names = {'history', 'output'};
            for k = 1:2
                work = arrays{k};
                for solver = 1:size(work,2)
                    observed = work(:,solver,:); finite = observed(isfinite(observed));
                    low = optiprofiler_internal.EvalReport.null(); high = low;
                    if ~isempty(finite), low = min(finite); high = max(finite); end
                    entry.(names{k}){end+1} = struct('solver_index', solver, 'hits', numel(finite), ...
                        'total', numel(observed), 'not_hit', sum(isnan(observed(:))), ...
                        'nonfinite_other', sum(isinf(observed(:))), 'work_evaluations_min', low, 'work_evaluations_max', high);
                end
            end
            self.convergence{end+1} = entry;
            self.document.profiles.convergence = self.convergence;
            self.document.profiles.work_summary = 'see_convergence';
            self.document.profiles.work_summary_reason = optiprofiler_internal.EvalReport.null();
            item = struct('id', sprintf('target-work-%d', numel(self.convergence)), ...
                'tolerance_index', numel(self.convergence), 'tolerance', tolerance, ...
                'axes', {{'problem', 'solver', 'run'}}, ...
                'axis_values', struct('problem', {problem_ids}, 'solver', {num2cell(1:size(history,2))}, ...
                    'run', {num2cell(1:size(history,3))}), ...
                'history', {optiprofiler_internal.EvalReport.tensor3(history)}, ...
                'output', {optiprofiler_internal.EvalReport.tensor3(output)});
            self.plotDocument.target_work{end+1} = item;
        end

        function setSource(self, path)
            self.sourcePath = path;
            [~, name, ext] = fileparts(path);
            record = struct('kind', 'archive', 'name', [name, ext], ...
                'path', self.relative(path), 'status', 'captured_before_filtering', ...
                'coverage_note', 'Only retained/replayed coverage is certified.');
            try
                [record.bytes, record.sha256] = optiprofiler_internal.EvalReport.digest(path);
            catch
                record.bytes = optiprofiler_internal.EvalReport.null(); record.sha256 = record.bytes;
                record.status = 'unavailable'; record.reason = 'archive_hash_unavailable';
                self.addDiagnostic('source_hash_unavailable', 'persistence');
            end
            self.document.source = record;
            self.write();
        end

        function recordRendering(self, status, scope, code)
            % CODE names the failed export the same way Python's diagnostics
            % do (history_render_failed, summary_pdf_merge_failed,
            % history_merge_failed). The stage reason is the most recent
            % failure code, as in Python; a later 'completed' call must not
            % replace it with the default code.
            if nargin < 4, code = 'history_render_failed'; end
            if strcmp(status, 'failed')
                self.renderingFailure = true;
                self.renderingFailureCode = code;
                self.addDiagnostic(code, 'rendering', scope);
            elseif strcmp(status, 'completed')
                self.renderingSuccess = true;
            end
            if self.renderingFailure
                state = 'failed';
                if self.renderingSuccess, state = 'partial'; end
                self.setStage('rendering', state, self.renderingFailureCode);
            elseif self.renderingSuccess
                self.setStage('rendering', 'completed');
            end
        end

        function finish(self, cause)
            % Reporting must never replace the original exception.
            try
                failed = nargin > 1;
                if strcmp(self.document.stages.rendering.status, 'running') && ~failed
                    if self.renderingFailure
                        state = 'failed';
                        if self.renderingSuccess, state = 'partial'; end
                        self.setStage('rendering', state, 'history_render_failed');
                    elseif self.renderingSuccess
                        self.setStage('rendering', 'completed');
                    elseif ~failed
                        self.setStage('rendering', 'unknown', 'no_render_receipt');
                    end
                end
                if failed
                    self.document.status = 'failed';
                    stages = fieldnames(self.document.stages);
                    for k = 1:numel(stages)
                        if strcmp(self.document.stages.(stages{k}).status, 'running')
                            self.setStage(stages{k}, 'failed', 'interrupted_by_exception');
                        end
                    end
                    self.addDiagnostic('benchmark_exception', 'runtime', struct('exception_type', cause.identifier));
                elseif self.document.coverage.load_failed > 0 && self.document.coverage.loaded == 0
                    self.document.status = 'failed';
                    self.completeNumerical();
                    self.setStage('scoring', 'not_applicable', 'no_loaded_problems');
                elseif self.document.coverage.load_failed > 0
                    self.document.status = 'partial';
                    self.completeNumerical();
                elseif self.document.coverage.loaded == 0
                    self.document.status = 'empty';
                    self.completeNumerical();
                    self.setStage('scoring', 'not_applicable', 'no_loaded_problems');
                else
                    self.document.status = 'completed';
                end
                if ~failed && self.renderingFailure, self.document.status = 'partial'; end
                if ~failed && strcmp(self.document.operation, 'benchmark') && ...
                        any(cellfun(@(p) strcmp(p.role, 'plain_reference') && strcmp(p.load_status, 'failed'), self.document.problems))
                    if ~strcmp(self.document.stages.numerical.status, 'failed')
                        self.setStage('numerical', 'partial', 'requested_plain_reference_incomplete');
                        self.document.status = 'partial';
                    end
                end
                self.harvest();
                if ~failed && isstruct(self.document.scores) && isfield(self.document.scores, 'semantics') && ...
                        startsWith(self.document.scores.semantics, 'single_problem') && ...
                        ~ismember(self.document.stages.persistence.status, {'failed','partial'})
                    % A direct problem never produces a reloadable archive,
                    % whether or not score_only was requested (same as Python).
                    self.setStage('persistence', 'not_applicable', 'single_problem_has_no_reload_archive');
                end
                if ~failed && strcmp(self.document.status, 'completed')
                    stages = fieldnames(self.document.stages);
                    for k = 1:numel(stages)
                        if ismember(self.document.stages.(stages{k}).status, {'partial','failed'})
                            self.document.status = 'partial';
                        end
                    end
                end
                self.document.timing.elapsed_seconds = toc(self.started);
                self.document.timing.finished_at = optiprofiler_internal.EvalReport.timestamp();
                self.write();
            catch write_error
                if nargin < 2, rethrow(write_error); end
            end
        end
    end
    methods (Access = private)
        function path = relative(self, path, base)
            % Relative paths may include '..'; never expose machine roots.
            if ~optiprofiler_internal.EvalReport.isAbsolute(path), path = fullfile(pwd, path); end
            if nargin < 3, base = fileparts(self.path); end
            if usejava('jvm')
                try
                    from = java.io.File(base); to = java.io.File(path);
                    path = char(from.toPath().normalize().relativize(to.toPath().normalize()).toString());
                    path = strrep(path, filesep, '/');
                catch
                    path = optiprofiler_internal.EvalReport.null();
                end
                return;
            end
            a = strsplit(base, filesep); b = strsplit(path, filesep); common = 0;
            while common < min(numel(a),numel(b)) && strcmp(a{common+1}, b{common+1}), common = common+1; end
            path = strjoin([repmat({'..'},1,numel(a)-common), b(common+1:end)], '/');
        end

        function harvest(self)
            if isempty(self.outputDirectory) || ~isfolder(self.outputDirectory), return; end
            if optiprofiler_internal.EvalReport.isLink(self.outputDirectory) || ...
                    ~strcmp(self.outputIdentity, optiprofiler_internal.EvalReport.fileIdentity(self.outputDirectory))
                self.addDiagnostic('artifact_directory_ownership_changed', 'persistence');
                self.setStage('persistence', 'partial', 'artifact_directory_ownership_changed');
                return;
            end
            pending = {self.outputDirectory}; artifacts = {};
            while ~isempty(pending) && numel(artifacts) < 2048
                directory = pending{1}; pending(1) = [];
                if optiprofiler_internal.EvalReport.isLink(directory), continue; end
                files = dir(directory);
                for k = 1:numel(files)
                    if numel(artifacts) >= 2048, break; end
                    if ismember(files(k).name, {'.', '..'}), continue; end
                    path = fullfile(files(k).folder, files(k).name);
                    if optiprofiler_internal.EvalReport.isLink(path), continue; end
                    if files(k).isdir, pending{end+1} = path; continue; end
                    if strcmp(path,self.path) || strcmp(path,self.plotPath) || strcmp(path,self.sourcePath), continue; end
                    [~,~,ext] = fileparts(path); kind = 'supporting_file'; media = 'application/octet-stream';
                    if ismember(ext,{'.pdf','.svg','.png'}), kind = 'plot'; end
                    if strcmp(ext,'.mat'), kind = 'raw_data'; end
                    if strcmp(ext,'.pdf'), media = 'application/pdf'; elseif strcmp(ext,'.svg'), media = 'image/svg+xml';
                    elseif strcmp(ext,'.png'), media = 'image/png'; elseif strcmp(ext,'.html'), media = 'text/html';
                    elseif ismember(ext,{'.txt','.m'}), media = 'text/plain'; end
                    try
                        [bytes, hash] = optiprofiler_internal.EvalReport.digest(path);
                        artifacts{end+1} = struct('path', self.relative(path, self.outputDirectory), 'media_type', media, ...
                            'kind', kind, 'status', 'present', 'bytes', bytes, 'sha256', hash);
                    catch
                        self.addDiagnostic('artifact_hash_unavailable', 'persistence');
                        self.setStage('persistence', 'partial', 'artifact_hash_unavailable');
                    end
                end
            end
            self.document.artifacts = artifacts;
            if numel(artifacts) >= 2048
                self.addDiagnostic('artifact_manifest_limit_reached', 'persistence', struct('limit', 2048));
            end
        end

        function index = problemIndex(self, id)
            index = [];
            if isKey(self.problemPositions, id), index = self.problemPositions(id); end
        end

        function upsert(self, collection, item)
            % Insert or replace by id; positions are tracked in a map so
            % large whole-library reports do not rescan every record.
            if strcmp(collection, 'plots'), positions = self.plotPositions; else, positions = self.historyPositions; end
            entries = self.plotDocument.(collection);
            if isKey(positions, item.id)
                entries{positions(item.id)} = item;
            else
                entries{end+1} = item;
                positions(item.id) = numel(entries);
            end
            self.plotDocument.(collection) = entries;
        end

        function storeHistoryPanels(self, panels, id, ordinal, source)
            for k = 1:numel(panels)
                panel = panels{k};
                panel.id = sprintf('history-plot-%d-%s-%s', ordinal, panel.channel, panel.mode);
                panel.problem_id = id; panel.observation_scope = source;
                panel.renderer_variants_ref = 'history';
                self.upsert('plots', panel);
                self.document.problems{ordinal}.plot_refs{k} = panel.id;
            end
        end

        function recount(self)
            coverage = self.document.coverage;
            coverage.selected = 0; coverage.loaded = 0; coverage.completed = 0; coverage.load_failed = 0;
            for k = 1:numel(self.document.problems)
                value = self.document.problems{k};
                if ~strcmp(value.role, 'primary'), continue; end
                coverage.selected = coverage.selected + 1;
                coverage.loaded = coverage.loaded + strcmp(value.load_status, 'loaded');
                coverage.completed = coverage.completed + strcmp(value.status, 'completed');
                coverage.load_failed = coverage.load_failed + strcmp(value.load_status, 'failed');
            end
            self.document.coverage = coverage;
        end

        function write(self)
            self.checkOwnership();
            self.checkPlotOwnership();
            self.publish(self.plotDocument, self.plotPath, true);
            receipt = struct('schema', 'optiprofiler.plot_data/1', 'path', self.relative(self.plotPath));
            try
                [receipt.bytes, receipt.sha256] = optiprofiler_internal.EvalReport.digest(self.plotPath);
                receipt.status = 'completed';
            catch cause
                if ~strcmp(cause.identifier,'OptiProfiler:EvalReportHash'), rethrow(cause); end
                % Missing system hashing tools must not suppress numerical
                % facts in -nojvm. This is explicitly NOT a verified pair;
                % consumers must reject the missing hash as such.
                receipt.bytes = optiprofiler_internal.EvalReport.null(); receipt.sha256 = receipt.bytes;
                receipt.status = 'partial'; receipt.reason = 'companion_sha256_unavailable';
                if ~any(cellfun(@(d) strcmp(d.code,'companion_sha256_unavailable'),self.document.diagnostics))
                    self.addDiagnostic('companion_sha256_unavailable','reporting');
                    warning('OptiProfiler:EvalReportHash', '%s', cause.message);
                end
                if strcmp(self.document.status,'completed'), self.document.status = 'partial'; end
            end
            self.checkPlotOwnership();
            receipt.history_count = numel(self.plotDocument.histories);
            receipt.plot_count = numel(self.plotDocument.plots);
            self.document.plot_data = receipt;
            self.document.report_files = struct('permission_policy', 'owner_read_write_only_best_effort', ...
                'permissions_applied', self.permissionsApplied, ...
                'platform_note', 'unix_chmod_600_after_each_publish;not_enforced_on_windows;windows_file_identity_is_creation_time_size_and_mtime_not_a_file_index;windows_directory_identity_is_creation_time_only;directory_privacy_is_the_caller_responsibility');
            snapshot = self.document;
            if strcmp(snapshot.operation, 'load')
                snapshot.coverage.load_failed = optiprofiler_internal.EvalReport.null();
                snapshot.coverage.reason = 'original_selection_and_load_failures_not_retained';
            end
            self.publish(snapshot, self.path, false);
        end

        function publish(self, snapshot, target, is_plot)
            encoded = unicode2native(jsonencode(optiprofiler_internal.EvalReport.safe(snapshot)), 'UTF-8');
            % A fresh, exclusively created same-directory staging file avoids
            % overwriting unrelated adjacent files and permits atomic rename.
            if usejava('jvm')
                parent = java.io.File(fileparts(self.path));
                file = java.io.File.createTempFile('op-eval-report-', '.json', parent);
                stage = char(file.getPath());
            else
                [status, stage] = system(['/usr/bin/mktemp ', optiprofiler_internal.EvalReport.quote(fullfile(fileparts(self.path), '.op-eval-report.XXXXXXXX'))]);
                if status ~= 0, error('OptiProfiler:EvalReportWrite', 'Could not reserve report stage.'); end
                stage = strtrim(stage);
            end
            fid = fopen(stage, 'w');
            if fid < 0, error('OptiProfiler:EvalReportWrite', 'Could not write report stage.'); end
            guard = onCleanup(@() fclose(fid));
            written = fwrite(fid, encoded, 'uint8');
            clear guard;
            if written ~= numel(encoded), error('OptiProfiler:EvalReportWrite', 'Incomplete report write.'); end
            if is_plot, self.checkPlotOwnership(); else, self.checkOwnership(); end
            self.replaceFile(stage, target);
            identity = optiprofiler_internal.EvalReport.fileIdentity(target);
            if is_plot, self.plotIdentity = identity; else, self.ownedIdentity = identity; end
            self.restrictPermissions(target);
        end

        function restrictPermissions(self, target)
            % Best effort owner-only mode, mirroring Python's 0600 files. The
            % staging file inherits the umask, so restrict after every
            % publish; a failure is recorded, never fatal for the numbers.
            applied = optiprofiler_internal.EvalReport.null();
            if isunix
                applied = false;
                try
                    % chmod works identically with and without the JVM and on
                    % macOS; the outcome is verified, never assumed.
                    [status, ~] = system(['chmod 600 ', optiprofiler_internal.EvalReport.quote(target)]);
                    [ok, info] = fileattrib(target);
                    applied = status == 0 && ok && ~info.GroupRead && ~info.GroupWrite && ~info.OtherRead && ~info.OtherWrite;
                catch
                end
            end
            if isstruct(self.permissionsApplied) || isstruct(applied)
                self.permissionsApplied = applied;
            else
                self.permissionsApplied = self.permissionsApplied && applied;
            end
        end

        function checkOwnership(self)
            if optiprofiler_internal.EvalReport.isLink(self.path) || ...
                    ~strcmp(self.ownedIdentity, optiprofiler_internal.EvalReport.fileIdentity(self.path))
                error('OptiProfiler:EvalReportOwnership', 'The report target no longer belongs to this invocation.');
            end
        end

        function checkPlotOwnership(self)
            if optiprofiler_internal.EvalReport.isLink(self.plotPath) || ...
                    ~strcmp(self.plotIdentity, optiprofiler_internal.EvalReport.fileIdentity(self.plotPath))
                error('OptiProfiler:EvalReportOwnership', 'The plot-data target no longer belongs to this invocation.');
            end
        end
    end
    methods (Static)
        function yes = isLink(path)
        %ISLINK True for a symbolic link and, on Windows, for a junction or any
        % other reparse point that redirects the name. Reparse points that keep
        % the name (cloud-file placeholders, compressed or deduplicated files)
        % are ordinary entries. Public so the platform tests can exercise it.
            if usejava('jvm')
                target = java.io.File(path).toPath();
                yes = java.nio.file.Files.isSymbolicLink(target);
                if ~yes && ispc
                    yes = optiprofiler_internal.EvalReport.isNameSurrogate(target);
                end
            else
                [status, ~] = system(['test -L ', optiprofiler_internal.EvalReport.quote(path)]); yes = status == 0;
            end
        end
    end
    methods (Static, Access = private)
        function yes = isAbsolute(path)
            yes = startsWith(path, filesep) || ~isempty(regexp(path, '^[A-Za-z]:[\\/]', 'once')) || startsWith(path, '\\');
        end

        function text = timestamp()
            text = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS''+00:00'''));
        end

        function value = fileIdentity(path)
            if usejava('jvm')
                file = java.io.File(path);
                options = javaArray('java.nio.file.LinkOption',1);
                options(1) = java.nio.file.LinkOption.NOFOLLOW_LINKS;
                attributes = java.nio.file.Files.readAttributes(file.toPath(), ...
                    'basic:fileKey,creationTime,lastModifiedTime,size,isDirectory', options);
                key = attributes.get('fileKey');
                if ~isempty(key)
                    value = char(key.toString());
                elseif logical(attributes.get('isDirectory'))
                    % Windows directory: its modification time changes whenever
                    % an entry is added or removed (the benchmark keeps writing
                    % into its output directory) and its size is not
                    % meaningful, so only the creation time identifies it: a
                    % replaced directory is detected, a modified one is normal.
                    value = sprintf('windows-directory:%s', char(attributes.get('creationTime').toString()));
                else
                    % Windows: Java exposes no file index (fileKey is null).
                    % Creation time (100 ns NTFS resolution), size and
                    % modification time identify a replaced or rewritten file.
                    % Every publish creates a fresh file and records its new
                    % identity, so a foreign move or in-place write between
                    % publishes is detected. This is weaker than an inode or
                    % fileKey (forged timestamps are outside this model) and is
                    % stated in report_files.platform_note.
                    value = sprintf('windows:%s|%s|%d', char(attributes.get('creationTime').toString()), ...
                        char(attributes.get('lastModifiedTime').toString()), double(attributes.get('size')));
                end
            else
                format = 'stat -c ''%d:%i'' ';
                if ismac, format = 'stat -f ''%d:%i'' '; end
                [status,value] = system([format,optiprofiler_internal.EvalReport.quote(path)]);
                if status ~= 0, error('OptiProfiler:EvalReportOwnership', 'File identity is unavailable.'); end
                value = strtrim(value);
            end
        end

        function result = tensor(value)
            result = {};
            if isempty(value), return; end
            for i = 1:size(value,1)
                rows = {};
                for j = 1:size(value,2)
                    channels = {};
                    for k = 1:size(value,3), channels{k} = num2cell(reshape(value(i,j,k,:),1,[])); end
                    rows{j} = channels;
                end
                result{i} = rows;
            end
        end

        function result = tensor3(value)
            result = cell(1,size(value,1));
            for i = 1:size(value,1)
                rows = cell(1,size(value,2));
                for j = 1:size(value,2), rows{j} = num2cell(reshape(value(i,j,:),1,[])); end
                result{i} = rows;
            end
        end

        function value = identifier()
            if usejava('jvm'), value = char(java.util.UUID.randomUUID()); return; end
            persistent sequence;
            if isempty(sequence), sequence = 0; end
            sequence = sequence+1;
            value = sprintf('matlab-%s-%d-%d', char(datetime('now','Format','yyyyMMddHHmmssSSS')), feature('getpid'), sequence);
        end

        function yes = isNameSurrogate(target)
            % Java exposes no reparse tag. An entry redirects the name when its
            % resolved real path differs from the resolved path of its parent
            % joined with the entry's own canonical name; an entry that exists
            % but cannot be resolved (a dangling junction) is a link as well.
            % java.nio methods taking LinkOption varargs bind from MATLAB only
            % with an explicit, possibly empty, LinkOption[] argument.
            follow = javaArray('java.nio.file.LinkOption', 0);
            nofollow = javaArray('java.nio.file.LinkOption', 1);
            nofollow(1) = java.nio.file.LinkOption.NOFOLLOW_LINKS;
            yes = false;
            if ~java.nio.file.Files.exists(target, nofollow), return; end
            parent = target.toAbsolutePath().getParent();
            if isempty(parent), return; end
            try
                resolved = target.toRealPath(follow);
                canonical = target.toRealPath(nofollow);
                yes = ~parent.toRealPath(follow).resolve(canonical.getFileName()).equals(resolved);
            catch
                attributes = java.nio.file.Files.readAttributes(target, 'basic:isOther', nofollow);
                yes = logical(attributes.get('isOther'));
            end
        end

        function [bytes, hash] = digest(path)
            % Best-effort checks for a directory assumed private to the caller:
            % links are refused before the open and size/modification time are
            % compared after hashing; the opened handle's identity is not
            % compared, so this detects ordinary concurrent changes, not a
            % hostile redirection timed between the checks.
            if optiprofiler_internal.EvalReport.isLink(path), error('OptiProfiler:EvalReportSymlink', 'Symlinks are not report artifacts.'); end
            before = dir(path);
            if numel(before) ~= 1 || before.isdir, error('OptiProfiler:EvalReportHash', 'Expected a regular file.'); end
            bytes = before.bytes;
            if usejava('jvm')
                digest = java.security.MessageDigest.getInstance('SHA-256');
                fid = fopen(path, 'rb');
                if fid < 0, error('OptiProfiler:EvalReportHash', 'File is unavailable.'); end
                guard = onCleanup(@() fclose(fid));
                while ~feof(fid)
                    block = fread(fid, 1048576, '*uint8');
                    % feof becomes true only after a read: empty files and a
                    % final empty block must not enter Java's overloaded
                    % update method. digest() still hashes empty input normally.
                    if ~isempty(block), digest.update(block); end
                end
                hash = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2).', 1, []));
            else
                % Minimal Linux systems commonly provide GNU sha256sum,
                % whereas macOS provides shasum. Neither plotting nor a JVM
                % is required to certify the requested numeric companion.
                [status, ~] = system('command -v sha256sum >/dev/null 2>&1');
                if status == 0, command = 'sha256sum ';
                else
                    [status, ~] = system('command -v shasum >/dev/null 2>&1');
                    if status ~= 0
                        error('OptiProfiler:EvalReportHash', ...
                            'EvalReport SHA256 without JVM requires sha256sum or shasum on PATH.');
                    end
                    command = 'shasum -a 256 ';
                end
                % Hash stdin so a backslash/newline in a valid filename does
                % not trigger the command's escaped-filename output prefix.
                [status, hash] = system([command, '< ', optiprofiler_internal.EvalReport.quote(path)]);
                if status ~= 0, error('OptiProfiler:EvalReportHash', 'SHA256 tool could not read the report artifact.'); end
                hash = hash(1:64);
            end
            after = dir(path);
            if numel(after) ~= 1 || before.bytes ~= after.bytes || before.datenum ~= after.datenum
                error('OptiProfiler:EvalReportHash', 'File changed while hashing.');
            end
        end

        function value = null()
            value = struct('eval_report_null', true);
        end
        function ok = reserve(path)
            if usejava('jvm')
                file = java.io.File(path);
                ok = file.createNewFile();
            elseif isunix
                [status, ~] = system(['(set -C; : > ', optiprofiler_internal.EvalReport.quote(path), ') 2>/dev/null']);
                ok = status == 0;
            else
                error('OptiProfiler:EvalReportJVM', 'Exclusive report creation requires the JVM on Windows.');
            end
        end

        function text = quote(value)
            text = ['''', strrep(value, '''', '''"''"'''), ''''];
        end

        function value = metric(history, count, output, initial)
            % Same record as Python eval_report._metric: absent
            % availability_reason means every observation was available;
            % absent first_invalid_evaluation_index means none was observed
            % and invalid_evaluations is then the integer 0.
            null = optiprofiler_internal.EvalReport.null();
            value = struct('output', output, 'initial', initial, 'best', null, ...
                'best_evaluation_index', null, 'invalid_evaluations', null);
            if isstruct(count)
                value.availability_reason = 'history_or_evaluation_count_unavailable';
                return;
            end
            if ~isempty(history)
                value.best = min(history, [], 'omitnan');
                if ~isnan(value.best), value.best_evaluation_index = find(history == value.best, 1); end
            end
            invalid_index = find(~isfinite(history), 1);
            if isempty(invalid_index) && numel(history) == count
                value.invalid_evaluations = 0;
            else
                value.invalid_evaluations = struct('nan', sum(isnan(history)), 'positive_infinity', sum(history == Inf), ...
                    'negative_infinity', sum(history == -Inf), 'observed_evaluations', numel(history));
            end
            if ~isempty(invalid_index), value.first_invalid_evaluation_index = invalid_index; end
            if numel(history) < count
                value.availability_reason = 'history_shorter_than_evaluation_count';
            elseif count == 0
                value.availability_reason = 'no_evaluations';
            elseif isstruct(output) || isstruct(initial)
                value.availability_reason = 'output_or_initial_unavailable';
            end
        end

        function value = binned(history, count, available)
            value = struct('status', 'available', 'reason', optiprofiler_internal.EvalReport.null(), ...
                'count', numel(history), 'total_evaluations', count, ...
                'representation', 'exact_samples', 'values', {num2cell(history)}, 'bins', {{}});
            if ~available
                value.status = 'unavailable'; value.reason = 'history_or_evaluation_count_unavailable';
                value.representation = optiprofiler_internal.EvalReport.null();
                value.values = optiprofiler_internal.EvalReport.null(); return;
            end
            if numel(history) < count, value.reason = 'history_shorter_than_evaluation_count'; end
            if isempty(history)
                if count == 0, value.reason = 'no_evaluations'; end
                return;
            end
            if numel(history) <= 64, return; end
            value.representation = 'bin_extrema'; value.values = optiprofiler_internal.EvalReport.null();
            % n_bins is always 32 here (numel > 64), so linspace steps are
            % multiples of n/32, an exact binary fraction: floor() yields the
            % same integer edges as Python's i*n//32 for every n.
            edges = floor(linspace(0, numel(history), min(32,numel(history))+1));
            for k = 1:numel(edges)-1
                start = edges(k)+1; stop = edges(k+1); part = history(start:stop);
                finite_indices = find(isfinite(part));
                low = optiprofiler_internal.EvalReport.null(); high = low;
                if ~isempty(finite_indices)
                    [minimum, lo] = min(part(finite_indices)); [maximum, hi] = max(part(finite_indices));
                    low = struct('value', minimum, 'evaluation_index', start+finite_indices(lo)-1);
                    high = struct('value', maximum, 'evaluation_index', start+finite_indices(hi)-1);
                end
                value.bins{end+1} = struct('start_index', start, 'end_index', stop, ...
                    'first', part(1), 'last', part(end), 'finite_min', low, 'finite_max', high, ...
                    'nonfinite', struct('nan', sum(isnan(part)), 'positive_infinity', sum(part == Inf), 'negative_infinity', sum(part == -Inf)));
            end
        end

        function value = configuration(value)
            if isstruct(value)
                names = fieldnames(value);
                for k = 1:numel(names)
                    key = names{k};
                    if ~isempty(regexpi(key, 'path|token|secret|password|credential|api_key')) || strcmp(key, 'problem')
                        value = rmfield(value, key);
                    else
                        value.(key) = optiprofiler_internal.EvalReport.configuration(value.(key));
                    end
                end
            elseif isa(value, 'function_handle')
                % Anonymous source/captured workspaces are deliberately omitted.
                value = struct('kind', 'callback', 'description', 'function_handle');
            elseif iscell(value)
                value = cellfun(@optiprofiler_internal.EvalReport.configuration, value, 'UniformOutput', false);
            elseif ischar(value) || isstring(value)
                value = char(value);
                if startsWith(value, '/') || ~isempty(regexp(value, '^[A-Za-z]:[\\/]', 'once'))
                    value = '[absolute_path_omitted]';
                elseif startsWith(value, '@(')
                    value = '[anonymous_callback_source_omitted]';
                end
                value = value(1:min(numel(value),256));
            elseif ~isnumeric(value) && ~islogical(value)
                value = struct('kind', 'object', 'description', 'not_serialized');
            end
        end

        function value = safe(value)
            if isstruct(value) && isscalar(value) && isfield(value, 'eval_report_null')
                value = NaN; % jsonencode emits JSON null; not a scientific NaN.
            elseif isstruct(value)
                names = fieldnames(value);
                for i = 1:numel(value)
                    for k = 1:numel(names), value(i).(names{k}) = optiprofiler_internal.EvalReport.safe(value(i).(names{k})); end
                end
            elseif iscell(value)
                value = cellfun(@optiprofiler_internal.EvalReport.safe, value, 'UniformOutput', false);
            elseif isnumeric(value) && any(~isfinite(value(:)))
                if isscalar(value)
                    reason = 'nan';
                    if value == Inf, reason = 'positive_infinity'; elseif value == -Inf, reason = 'negative_infinity'; end
                    value = struct('value', NaN, 'reason', reason);
                else
                    value = cellfun(@optiprofiler_internal.EvalReport.safe, num2cell(value), 'UniformOutput', false);
                end
            end
        end
    end
end
