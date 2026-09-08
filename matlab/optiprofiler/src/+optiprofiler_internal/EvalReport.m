classdef EvalReport < handle
%EVALREPORT Controller-private, observational JSON collector. Not a public API.
% Numeric completion is not convergence; recorded history excludes padding.
    properties (Access = private)
        path
        document
        started
        outputDirectory = ''
        maxEvalFactor = NaN
        replaceFile
        renderingFailure = false
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
            self.started = tic;
            stage = struct('status', 'unknown');
            self.document = struct('schema', 'optiprofiler.eval_report/1', ...
                'evaluation_id', optiprofiler_internal.EvalReport.identifier(), ...
                'operation', 'benchmark', 'status', 'running', ...
                'producer', struct('language', 'matlab', 'version', 'unknown', ...
                    'revision', optiprofiler_internal.EvalReport.null(), 'revision_reason', 'not_available'), ...
                'configuration', struct(), 'stages', struct('numerical', stage, 'scoring', stage, 'persistence', stage, 'rendering', stage), ...
                'coverage', struct('selected', 0, 'loaded', 0, 'completed', 0, 'load_failed', 0, ...
                    'scope', 'live_selection', 'original_selection_known', true), ...
                'problems', {{}}, 'scores', optiprofiler_internal.EvalReport.null(), ...
                'profiles', struct('plot_refs', {{}}, 'convergence', {{}}, ...
                    'work_summary', optiprofiler_internal.EvalReport.null(), 'work_summary_reason', 'not_supplied'), ...
                'artifact_root', '.', 'artifacts', {{}}, 'diagnostics', {{}}, ...
                'plot_data', optiprofiler_internal.EvalReport.null(), ...
                'source', optiprofiler_internal.EvalReport.null(), 'timing', struct());
            self.document.semantics = struct('index_base', 1, ...
                'best', 'componentwise_minimum_ignoring_nan; earliest retained evaluation attaining it', ...
                'budget', 'Reaching the cap does not identify the termination cause.', ...
                'convergence', 'Not inferred from solver return; target_work observes the existing profile construction.', ...
                'paths', 'artifacts.path is relative to artifact_root; artifact_root, plot_data.path and source.path are relative to the report parent.', ...
                'privacy', 'Controller-private facts; consumers must construct an allowlisted public feedback view.');
            self.plotDocument = struct('schema', 'optiprofiler.plot_data/1', ...
                'evaluation_id', self.document.evaluation_id, ...
                'semantics', struct('index_base', 1, 'histories', ...
                    'Up to 64 observations are exact_samples. Longer histories use lossy bins: at most 32 contiguous bins with exact endpoints, finite extrema/earliest indices and nonfinite counts, but no internal order. Never use bins as scoring inputs.', ...
                    'plots', 'Exact shared prepared numeric representation, not raw history or proof of a rendered file. score_only creates no figures. Resolve each renderer_variants_ref below before interpreting native PDF versus portable SVG styling.', ...
                    'error_bands', 'MATLAB meanstd uses sample standard deviation (ddof=1 for N>1, zero for N=1).', ...
                    'target_work', 'Existing complete work arrays; history is first qualifying evaluation, output is total evaluations at a qualifying returned point.'), ...
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
            self.document.configuration = struct('problem_options', optiprofiler_internal.EvalReport.configuration(problem_options), ...
                'profile_options', optiprofiler_internal.EvalReport.configuration(profile_options), ...
                'feature', struct('name', feature_value.name, 'options', optiprofiler_internal.EvalReport.configuration(feature_value.options)));
            self.maxEvalFactor = profile_options.max_eval_factor;
            if strcmp(self.document.operation, 'load')
                self.document.configuration.scope = 'current_load_selection_reanalysis_and_rendering';
                self.document.configuration.solver_execution_requested = false;
                self.document.configuration.original_execution_configuration = optiprofiler_internal.EvalReport.null();
                self.document.configuration.original_execution_reason = 'not_retained_in_mat_archive';
                self.document.configuration.feature = struct('name', optiprofiler_internal.EvalReport.null(), ...
                    'options', optiprofiler_internal.EvalReport.null(), 'scope', 'original_execution_unknown', ...
                    'reason', 'default_load_feature_is_not_original_execution_provenance');
            end
            if profile_options.score_only
                self.setStage('rendering', 'not_requested', 'score_only');
                self.setStage('persistence', 'not_requested', 'score_only');
            else
                self.setStage('rendering', 'unknown', 'requested_not_started');
                self.setStage('persistence', 'unknown', 'requested_not_started');
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
            entry = struct('id', id, ...
                'library', library, 'name', result.problem_name, 'role', role, ...
                'dimension', result.problem_dim, 'type', result.problem_type, ...
                'selection_status', 'selected', 'load_status', 'loaded', 'status', 'completed', 'availability', 'retained', 'runs', {{}}, 'plot_refs', {{}});
            for solver = 1:size(result.n_eval, 1)
                for run = 1:size(result.n_eval, 2)
                    count = result.n_eval(solver, run);
                    retained = min(count, size(result.fun_history, 3));
                    values = reshape(result.fun_history(solver, run, 1:retained), 1, []);
                    constraints = reshape(result.maxcv_history(solver, run, 1:min(count,size(result.maxcv_history,3))), 1, []);
                    merits = [];
                    if isfield(result, 'merit_history'), merits = reshape(result.merit_history(solver, run, 1:min(count,size(result.merit_history,3))), 1, []); end
                    budget = ceil(self.maxEvalFactor * result.problem_dim);
                    unknown = optiprofiler_internal.EvalReport.null();
                    item = struct('solver_index', solver, 'run_index', run, 'evaluations', count, ...
                        'budget', budget, 'budget_reached', count >= budget, ...
                        'objective', optiprofiler_internal.EvalReport.metric(values, result.fun_out(solver, run), result.fun_inits(run)), ...
                        'constraint', optiprofiler_internal.EvalReport.metric(constraints, result.maxcv_out(solver, run), result.maxcv_inits(run)), ...
                        'merit', struct('output', unknown, 'initial', unknown, 'best', unknown, ...
                            'best_evaluation_index', unknown, 'first_invalid_evaluation_index', unknown, ...
                            'invalid_evaluations', unknown, 'availability', 'not_retained'), ...
                        'abnormal_termination', result.solver_abnormal_termination(solver, run), ...
                        'output_fallback', result.solver_output_fallback(solver, run), ...
                        'execution', struct('kind', 'unknown', 'source_run_index', unknown), ...
                        'history_ref', sprintf('history-%d-%d-%d', ordinal, solver, run));
                    if isfield(result, 'eval_report_metadata') && isfield(result.eval_report_metadata, 'real_n_runs')
                        item.execution.kind = 'actual';
                        item.execution.source_run_index = run;
                        if run > result.eval_report_metadata.real_n_runs(solver)
                            item.execution.kind = 'repeated';
                            item.execution.source_run_index = 1;
                        end
                    end
                    if isfield(result, 'merit_history')
                        item.merit = optiprofiler_internal.EvalReport.metric(merits, unknown, result.merit_inits(run));
                        item.merit.output_availability = 'not_retained';
                        if isfield(result, 'merit_out'), item.merit.output = result.merit_out(solver, run); item.merit.output_availability = 'retained'; end
                    end
                    item.budget_reason = unknown;
                    item.convergence = unknown;
                    item.convergence_reason = 'not_inferred_from_solver_return';
                    item.termination_metadata_reason = unknown;
                    if isstruct(item.abnormal_termination) || isstruct(item.output_fallback)
                        item.termination_metadata_reason = 'solver_termination_metadata_not_retained';
                    end
                    if strcmp(self.document.operation, 'load')
                        item.budget = unknown; item.budget_reached = unknown;
                        item.budget_reason = 'original_execution_budget_not_retained';
                    end
                    if strcmp(item.execution.kind, 'actual'), item.execution.source_run_index = unknown; end
                    if strcmp(item.execution.kind, 'unknown'), item.execution.reason = 'execution_metadata_not_retained'; end
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
            if isfield(result, 'eval_report_metadata') && isfield(result.eval_report_metadata, 'render_status') && ...
                    (isempty(index) || isempty(self.document.problems{index}.runs))
                self.recordRendering(result.eval_report_metadata.render_status, struct('library', library, 'name', result.problem_name));
            end
            if isempty(index)
                self.document.problems{end+1} = entry;
            else
                old = self.document.problems{index};
                if isfield(old,'plot_refs'), entry.plot_refs = old.plot_refs; end
                if ~isfield(result, 'eval_report_metadata') && numel(old.runs) == numel(entry.runs)
                    for k = 1:numel(entry.runs), entry.runs{k}.execution = old.runs{k}.execution; end
                end
                self.document.problems{index} = entry;
            end
            self.recount();
            % Shared preparation consumes retained arrays only. It never
            % evaluates an oracle, feature, solver or custom merit callback.
            if isfield(result, 'eval_report_metadata') && isfield(result.eval_report_metadata, 'plot_presentation')
                self.captureHistoryPresentation(result.eval_report_metadata.plot_presentation, library, result.problem_name, role);
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

        function selection(self, library, names, role)
            for k = 1:numel(names)
                id = jsonencode({library, names{k}, role});
                if isempty(self.problemIndex(id))
                    self.document.problems{end+1} = struct('id', id, 'library', library, 'name', names{k}, 'role', role, ...
                        'dimension', optiprofiler_internal.EvalReport.null(), 'type', optiprofiler_internal.EvalReport.null(), ...
                        'selection_status', 'selected', 'load_status', 'pending', 'status', 'running', ...
                        'availability', 'not_loaded', 'runs', {{}}, 'plot_refs', {{}});
                end
            end
            self.recount();
        end

        function loadFailed(self, library, name, role)
            index = self.problemIndex(jsonencode({library, name, role}));
            self.document.problems{index}.load_status = 'failed';
            self.document.problems{index}.status = 'failed';
            self.addDiagnostic('problem_load_failed', 'numerical', struct('library', library, 'name', name, 'role', role));
            self.recount();
        end

        function completeNumerical(self)
            if self.document.coverage.load_failed > 0 && self.document.coverage.loaded == 0
                self.setStage('numerical', 'failed', 'all_selected_problem_loads_failed');
            elseif self.document.coverage.load_failed > 0
                self.setStage('numerical', 'partial', 'problem_load_failed');
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
                        'merit_outs','merit_out'; 'n_evals','n_eval'; ...
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
                self.document.profiles.plot_refs{end+1} = item.id;
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

        function recordRendering(self, status, scope)
            if strcmp(status, 'failed')
                self.renderingFailure = true;
                self.addDiagnostic('render_failed', 'rendering', scope);
            elseif strcmp(status, 'completed')
                self.renderingSuccess = true;
            end
            if self.renderingFailure
                state = 'failed';
                if self.renderingSuccess, state = 'partial'; end
                self.setStage('rendering', state, 'requested_render_failed');
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
                        self.setStage('rendering', state, 'requested_render_failed');
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
                            self.setStage(stages{k}, 'failed', 'benchmark_exception');
                        end
                    end
                    self.addDiagnostic('benchmark_exception', 'controller', struct('exception_type', cause.identifier));
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
                        ~ismember(self.document.stages.persistence.status, {'failed','partial','not_requested'})
                    self.setStage('persistence', 'not_applicable', 'single_problem_raw_archive_not_produced');
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
            index = find(cellfun(@(value) strcmp(value.id, id), self.document.problems), 1);
        end

        function upsert(self, collection, item)
            entries = self.plotDocument.(collection);
            index = find(cellfun(@(entry) strcmp(entry.id, item.id), entries), 1);
            if isempty(index), entries{end+1} = item; else, entries{index} = item; end
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
            receipt_status = 'completed'; reason = optiprofiler_internal.EvalReport.null();
            try
                [bytes, hash] = optiprofiler_internal.EvalReport.digest(self.plotPath);
            catch cause
                if ~strcmp(cause.identifier,'OptiProfiler:EvalReportHash'), rethrow(cause); end
                % Missing system hashing tools must not suppress numerical
                % facts in -nojvm. This is explicitly NOT a verified pair;
                % consumers must reject the missing hash as such.
                bytes = optiprofiler_internal.EvalReport.null(); hash = bytes;
                receipt_status = 'partial'; reason = 'companion_sha256_unavailable';
                if ~any(cellfun(@(d) strcmp(d.code,'companion_sha256_unavailable'),self.document.diagnostics))
                    self.addDiagnostic('companion_sha256_unavailable','report');
                    warning('OptiProfiler:EvalReportHash', '%s', cause.message);
                end
                if strcmp(self.document.status,'completed'), self.document.status = 'partial'; end
            end
            self.checkPlotOwnership();
            self.document.plot_data = struct('schema', 'optiprofiler.plot_data/1', ...
                'path', self.relative(self.plotPath), 'bytes', bytes, 'sha256', hash, ...
                'status', receipt_status, 'reason', reason, 'history_count', numel(self.plotDocument.histories), ...
                'plot_count', numel(self.plotDocument.plots));
            snapshot = self.document;
            if strcmp(snapshot.operation, 'load')
                snapshot.coverage.load_failed = optiprofiler_internal.EvalReport.null();
                snapshot.coverage.load_failed_reason = 'original_selection_failures_not_retained';
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
    methods (Static, Access = private)
        function yes = isAbsolute(path)
            yes = startsWith(path, filesep) || ~isempty(regexp(path, '^[A-Za-z]:[\\/]', 'once')) || startsWith(path, '\\');
        end

        function value = fileIdentity(path)
            if usejava('jvm')
                file = java.io.File(path);
                options = javaArray('java.nio.file.LinkOption',1);
                options(1) = java.nio.file.LinkOption.NOFOLLOW_LINKS;
                attributes = java.nio.file.Files.readAttributes(file.toPath(), 'basic:fileKey', options);
                key = attributes.get('fileKey'); value = char(key.toString());
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

        function yes = isLink(path)
            if usejava('jvm')
                file = java.io.File(path); yes = java.nio.file.Files.isSymbolicLink(file.toPath());
            else
                [status, ~] = system(['test -L ', optiprofiler_internal.EvalReport.quote(path)]); yes = status == 0;
            end
        end

        function [bytes, hash] = digest(path)
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

        function value = metric(history, output, initial)
            best = optiprofiler_internal.EvalReport.null();
            best_index = optiprofiler_internal.EvalReport.null();
            if ~isempty(history)
                best = min(history, [], 'omitnan');
                if ~isnan(best), best_index = find(history == best, 1); end
            end
            invalid_index = find(~isfinite(history), 1);
            if isempty(invalid_index), invalid_index = optiprofiler_internal.EvalReport.null(); end
            value = struct('output', output, 'initial', initial, 'best', best, ...
                'best_evaluation_index', best_index, 'first_invalid_evaluation_index', invalid_index, ...
                'invalid_evaluations', struct('nan', sum(isnan(history)), 'positive_infinity', sum(history == Inf), ...
                    'negative_infinity', sum(history == -Inf), 'observed_evaluations', numel(history)), ...
                'history_reason', optiprofiler_internal.EvalReport.null());
            if isempty(history), value.history_reason = 'no_evaluations'; end
        end

        function value = binned(history, count, available)
            value = struct('status', 'available', 'reason', optiprofiler_internal.EvalReport.null(), ...
                'count', numel(history), 'total_evaluations', count, ...
                'representation', 'exact_samples', 'values', {num2cell(history)}, 'bins', {{}});
            if ~available
                value.status = 'unavailable'; value.reason = 'not_retained';
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
