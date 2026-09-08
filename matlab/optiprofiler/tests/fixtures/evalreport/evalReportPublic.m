function evalReportPublic(source_root, output_root, slice)
% Public benchmark/report JSON regression; run only on syu-ubuntu.
    if nargin < 3, slice = 'identity'; end
    addpath(fullfile(source_root, 'matlab', 'optiprofiler', 'src'));
    if ~isfolder(output_root), mkdir(output_root); end
    problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [1; 2], 'name', 'report_quadratic'));
    options = struct('problem', problem, 'score_only', true, 'silent', true, ...
        'n_jobs', 1, 'max_eval_factor', 4, 'solver_names', {{'stay', 'zero'}});
    solvers = {@stay, @zero};
    schema_root = fullfile(source_root, 'doc', 'source', '_static');
    % Every report a slice reads is validated against the shared schemas
    % first: the schemas are the Python/MATLAB reader contract, so a field
    % spelled differently by this emitter fails here.
    readReport = @(path) readValidated(path, fullfile(schema_root, 'eval_report.schema.json'));
    readDetail = @(path) readValidated(path, fullfile(schema_root, 'plot_data.schema.json'));
    if strcmp(slice, 'identity')
        options.solver_names = {'保持', '归零'};
        options.problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [1; 2], 'name', '中文二次'));
        rng(71); state_before = rng;
        [s0, p0, c0] = benchmark(solvers, options);
        state_without = rng;
        rng(state_before);
        options.report_path = fullfile(output_root, 'identity.json');
        [s1, p1, c1] = benchmark(solvers, options);
        assert(isequaln(s0, s1) && isequaln(p0, p1) && isequaln(c0, c1), 'Report changed benchmark outputs.');
        assert(isequal(state_without, rng), 'Report consumed scientific RNG state.');
        report = readReport(options.report_path);
        assert(strcmp(report.schema, 'optiprofiler.eval_report/1'));
        assert(strcmp(report.status, 'completed'));
        assert(strcmp(report.stages.numerical.status, 'completed'));
        assert(strcmp(report.stages.rendering.status, 'not_requested'));
        assert(report.coverage.selected == 1 && report.coverage.completed == 1);
        assert(isempty(report.artifacts), 'score_only report enabled artifacts.');
        detail = readDetail(fullfile(output_root, report.plot_data.path));
        assert(strcmp(detail.schema, 'optiprofiler.plot_data/1'));
        assert(strcmp(detail.evaluation_id, report.evaluation_id));
        assert(strcmp(report.plot_data.sha256, hashFile(fullfile(output_root, report.plot_data.path))));
        assert(numel(detail.histories) == 2 && ~isfield(report.problems.runs, 'history_preview'));
        assert(strcmp(report.problems.runs(1).history_ref, detail.histories(1).id));
        assert(strcmp(detail.histories(1).channels.objective.representation, 'exact_samples'));
        assert(isequal(detail.histories(1).channels.objective.values, 5) && isempty(detail.histories(1).channels.objective.bins));
        assert(~isempty(detail.plots), 'score_only must retain numeric history presentations without figures.');
        assert(report.problems.runs(1).evaluations == 1);
        assert(report.problems.runs(2).objective.output == 0);
        assert(strcmp(report.problems.name, '中文二次') && strcmp(report.scores.solver_names{1}, '保持'), 'UTF-8 labels were corrupted.');
        % Shared vocabulary with Python (pinned by the schemas).
        assert(strcmp(report.producer.language, 'matlab') && strcmp(report.problems.library, 'user'));
        assert(strcmp(report.problems.id, '["user","中文二次","primary"]'), 'Direct problems use the user library label.');
        assert(isequal(sort(fieldnames(report.semantics)), sort({'index_base'; 'metric_best'; 'budget'; 'convergence'; 'run_defaults'; 'configuration'; 'paths'; 'privacy'})));
        assert(isequal(sort(fieldnames(detail.semantics)), sort({'index_base'; 'history_bins'; 'history_plots'; 'profile_plots'; 'error_bands'; 'target_work'; 'privacy'; 'renderer_variants'})));
        run = report.problems.runs(1);
        assert(isequal(sort(fieldnames(run)), sort({'solver_index'; 'run_index'; 'evaluations'; 'budget_reached'; 'objective'; 'constraint'; 'merit'; 'abnormal_termination'; 'output_fallback'; 'execution'; 'history_ref'; 'oracle_seed'; 'elapsed_seconds'})), 'Run record keys differ from the shared contract.');
        seed = report.configuration.effective.profile_options.seed;
        assert(run.oracle_seed == mod(23333 * seed + 211 * 1, 2^32), 'Actual MATLAB seed rule (1-based) must be recorded, not Python''s.');
        assert(isequal(run.objective.invalid_evaluations, 0) && ~isfield(run.objective, 'first_invalid_evaluation_index'));
        assert(strcmp(run.merit.availability_reason, 'output_or_initial_unavailable'));
        assert(strcmp(report.stages.persistence.reason, 'single_problem_has_no_reload_archive'));
        assert(strcmp(report.stages.rendering.reason, 'score_only'));
        assert(strcmp(report.configuration.scope, 'current_benchmark_execution') && isfield(report.configuration, 'request') && isfield(report.configuration, 'effective'));
        panel = detail.plots(1);
        assert(strcmp(panel.x_transform, 'evaluation_index/(dimension+1)') && panel.display_limit == 1e100 && panel.padded_length == 8);
        assert(strcmp(panel.aggregation_trigger, 'padded_history_length_above_1002') && panel.std_ddof == 0 && panel.n_runs == 1);
        % Best-effort owner-only file modes on Unix, and the report says so.
        assert(strcmp(report.report_files.permission_policy, 'owner_read_write_only_best_effort'));
        if isunix
            for target = {options.report_path, fullfile(output_root, report.plot_data.path)}
                [~, info] = fileattrib(target{1});
                assert(~info.GroupRead && ~info.OtherRead && ~info.GroupWrite && ~info.OtherWrite, 'Report files must be owner-only on Unix.');
            end
            assert(report.report_files.permissions_applied, 'The report must record that owner-only modes were applied.');
        end
        fprintf('PASS identity: public outputs, RNG, schema, coverage, score_only, shared vocabulary, permissions\n');
    elseif strcmp(slice, 'coverage')
        fixture_root = fileparts(mfilename('fullpath'));
        setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(output_root, 'registry.mat'));
        for name = {'report_a', 'report_b'}
            registerProblemLibrary(struct('name', name{1}, 'root', fixture_root, ...
                'select_function', 'eval_fixture_select', 'load_function', 'eval_fixture_load'));
        end
        options = rmfield(options, 'problem');
        options.plibs = {'report_a', 'report_b'};
        options.n_runs = 3;
        options.max_tol_order = 2;
        options.report_path = fullfile(output_root, 'coverage.json');
        benchmark(solvers, options);
        report = readReport(options.report_path);
        assert(strcmp(report.status, 'partial'));
        assert(report.coverage.selected == 4 && report.coverage.loaded == 2 && report.coverage.load_failed == 2);
        entries = report.problems;
        assert(numel(unique({entries.id})) == 4, 'Same-name cross-library identity collapsed.');
        assert(isequal({entries.id}, {'["report_a","shared","primary"]', '["report_a","unloadable","primary"]', ...
            '["report_b","shared","primary"]', '["report_b","unloadable","primary"]'}), 'Problem order must follow selection order.');
        assert(strcmp(entries(1).provider.name, 'report_a') && ~isfield(entries(1).provider, 'root'), 'Provider metadata must not expose machine paths.');
        failed_entry = entries(strcmp({entries.name}, 'unloadable'));
        assert(all(strcmp({failed_entry.load_status}, 'failed')) && all(strcmp({failed_entry.status}, 'failed')));
        load_diagnostics = report.diagnostics(strcmp({report.diagnostics.code}, 'problem_load_failed'));
        assert(numel(load_diagnostics) == 2 && isfield(load_diagnostics(1).scope, 'problem') && ~isempty(load_diagnostics(1).scope.exception_type), ...
            'Load failures must carry the problem name and exception identifier.');
        good = entries(strcmp({entries.name}, 'shared'));
        assert(numel(good) == 2 && numel(good(1).runs) == 6);
        assert(strcmp(good(1).runs(2).execution.kind, 'repeated'));
        assert(good(1).runs(2).execution.source_run_index == 1);
        before = fileread(options.report_path);
        try
            benchmark(solvers, options);
            error('Expected existing-report rejection.');
        catch cause
            assert(strcmp(cause.identifier, 'OptiProfiler:EvalReportExists'));
        end
        assert(strcmp(before, fileread(options.report_path)), 'Collision changed existing report.');
        options.report_path = fullfile(output_root, 'companion_collision.json');
        companion_path = fullfile(output_root, 'companion_collision.plot_data.json');
        fid = fopen(companion_path, 'w'); fprintf(fid, 'foreign companion'); fclose(fid);
        try
            benchmark({@mustNotRun, @mustNotRun}, options);
            error('Expected existing-companion rejection.');
        catch cause
            assert(strcmp(cause.identifier, 'OptiProfiler:EvalReportExists'));
        end
        assert(strcmp(fileread(companion_path), 'foreign companion'));
        assert(~isfile(options.report_path), 'Failed pair reservation left an owned main-file shell.');
        options.problem_names = {'absent'};
        options.report_path = fullfile(output_root, 'empty.json');
        benchmark(solvers, options);
        report = readReport(options.report_path);
        assert(strcmp(report.status, 'empty') && report.coverage.selected == 0);
        fprintf('PASS coverage: live failed loads, cross-library names, repeats, collision, empty\n');
    elseif strcmp(slice, 'semantics')
        options.problem = Problem(struct('fun', @undefinedObjective, 'x0', 1, 'name', 'nonfinite'));
        options.max_eval_factor = 128;
        options.report_path = fullfile(output_root, 'nonfinite.json');
        benchmark({@many, @zero}, options);
        report = readReport(options.report_path);
        assert(isempty(report.producer.revision), 'Unknown revision must be null, not scientific NaN.');
        assert(strcmp(report.stages.numerical.status, 'completed'));
        first = report.problems.runs(1);
        assert(first.evaluations == 100 && ~first.budget_reached);
        assert(report.problems.budget.evaluations == 128 && ~isfield(first, 'budget'), 'Budget is a per-problem fact inherited by runs.');
        assert(~isfield(first, 'convergence'), 'Convergence must never be inferred per run.');
        assert(first.objective.invalid_evaluations.nan == 100 && first.objective.invalid_evaluations.observed_evaluations == 100, 'Padded history counted as evaluations.');
        detail = readDetail(fullfile(output_root, report.plot_data.path));
        bins = detail.histories(1).channels.objective.bins;
        assert(numel(bins) == 32);
        assert(sum(arrayfun(@(bin) bin.nonfinite.nan, bins)) == 100);
        assert(first.objective.first_invalid_evaluation_index == 1);
        assert(strcmp(first.objective.output.reason, 'nan'));
        assert(bins(1).start_index == 1 && bins(end).end_index == 100);
        options.report_path = fullfile(output_root, 'fatal.json');
        options.merit_fun = @fatalMerit;
        try
            benchmark(solvers, options);
            error('Expected merit failure.');
        catch cause
            assert(strcmp(cause.identifier, 'MATLAB:benchmark:merit_fun_error'));
        end
        report = readReport(options.report_path);
        assert(strcmp(report.status, 'failed'));
        assert(report.coverage.completed == 1 && strcmp(report.stages.numerical.status, 'completed'));
        assert(strcmp(report.stages.scoring.status, 'failed'));
        fprintf('PASS semantics: explicit nonfinite, bounded actual history, fatal score separation\n');
    elseif strcmp(slice, 'compact')
        options.problem = Problem(struct('fun', @spikeObjective, 'x0', 1, 'name', 'isolated_spike'));
        options.max_eval_factor = 3000;
        options.report_path = fullfile(output_root, 'compact.json');
        [before_s,before_p,before_c] = benchmark({@longWalk, @zero}, rmfield(options,'report_path'));
        [after_s,after_p,after_c] = benchmark({@longWalk, @zero}, options);
        assert(isequaln(before_s,after_s) && isequaln(before_p,after_p) && isequaln(before_c,after_c));
        report = readReport(options.report_path);
        detail = readDetail(fullfile(output_root,report.plot_data.path));
        first = report.problems.runs(1);
        assert(first.objective.best == -100 && first.objective.best_evaluation_index == 51);
        assert(first.objective.first_invalid_evaluation_index == 73);
        bins = detail.histories(1).channels.objective.bins;
        assert(numel(bins) == 32 && bins(end).end_index == 2500);
        assert(bins(1).finite_max.value == 1e6 && bins(1).finite_max.evaluation_index == 2, ...
            'An isolated spike missed by old uniform previews must survive bin extrema.');
        assert(bins(1).finite_min.value == -100 && bins(1).finite_min.evaluation_index == 51);
        assert(bins(1).nonfinite.nan == 1);
        assert(~isfield(first,'history_preview') && ~isfield(first.objective,'best_semantics'));
        assert(numel(detail.plots(1).series(1).x) <= 1002 && strcmp(detail.plots(1).fidelity,'exact_rendered_data'));
        assert(detail.plots(1).padded_length == 3000, 'Padded length (max_eval) drives block aggregation.');
        edges = [[bins.start_index]; [bins.end_index]];
        assert(isequal(edges, [floor((0:31) * 2500 / 32) + 1; floor((1:32) * 2500 / 32)]), 'Bin edges must be the exact integer edges shared with Python.');
        assert(numel(dir(fullfile(output_root,'*.json'))) == 2, 'score_only must produce only the two requested JSON files.');
        fprintf('PASS compact: exact extrema and indices, bounded bins/shared plots, unchanged scores, two files only\n');
    elseif strcmp(slice, 'aggregation')
        % Block aggregation is keyed on the PADDED history length
        % max_eval = ceil(max_eval_factor*dimension). Above 1002 the shared
        % renderer preparation keeps about 1000 interior points, so an array
        % position is no longer an evaluation number: readers must locate a
        % sample through evaluation_indices. MATLAB's rule is verified here
        % without altering it; the Python suite covers the same boundary.
        options.problem = Problem(struct('fun', @spikeAt17, 'x0', 1, 'name', 'padded_spike'));
        options.max_eval_factor = 1002;
        options.report_path = fullfile(output_root, 'boundary_1002.json');
        benchmark({@walk1003, @zero}, options);
        report = readReport(options.report_path);
        detail = readDetail(fullfile(output_root, report.plot_data.path));
        raw = detail.plots(strcmp({detail.plots.mode}, 'raw'));
        series = raw.series(1);
        assert(report.problems.budget.evaluations == 1002 && report.problems.runs(1).evaluations == 1002);
        assert(isequal(series.evaluation_indices(:)', 1:1002), 'At the 1002 boundary every evaluation keeps its own vertex.');
        assert(series.mean(17) == 1e6, 'Without aggregation the position equals the evaluation index.');
        options.max_eval_factor = 1003;
        options.report_path = fullfile(output_root, 'boundary_1003.json');
        benchmark({@walk1003, @zero}, options);
        report = readReport(options.report_path);
        detail = readDetail(fullfile(output_root, report.plot_data.path));
        raw = detail.plots(strcmp({detail.plots.mode}, 'raw'));
        series = raw.series(1);
        assert(raw.padded_length == 1003 && strcmp(raw.aggregation_trigger, 'padded_history_length_above_1002'));
        % 1001 interior positions in 1000 blocks: exactly one block holds two
        % points and its 'min' representative drops the other evaluation.
        % (MATLAB linspace(1, n, 1) places that block last; NumPy's
        % linspace(0, n-1, 1) places it first. Neither renderer is altered.)
        indices = series.evaluation_indices(:)';
        assert(numel(indices) < 1003 && indices(1) == 1 && indices(end) == 1003 && any(diff(indices) > 1), ...
            'One interior evaluation must be dropped by the two-point min block.');
        position = find(indices == 17, 1);
        assert(~isempty(position) && series.mean(position) == 1e6, 'The spike must be located through its evaluation index.');
        assert(find(indices == 1003, 1) ~= 1003, 'Above 1002 an array position is not an evaluation number.');
        bins = detail.histories(1).channels.objective.bins;
        assert(bins(1).finite_max.value == 1e6 && bins(1).finite_max.evaluation_index == 17, 'Bin extrema keep the spike independently of display aggregation.');
        fprintf('PASS aggregation: padded-length trigger at 1002/1003, spike located by evaluation index, bins unchanged\n');
    elseif strcmp(slice, 'stateful_render')
        options.problem = Problem(struct('fun', @(x) x.^2, 'x0', 1, 'xl', 0, 'name', 'stateful_merit'));
        options.max_eval_factor = 1; options.score_only = false;
        options.savepath = output_root; options.draw_hist_plots = 'parallel';
        options.merit_fun = @statefulMerit;
        options.benchmark_id = 'without_report';
        statefulMerit('reset');
        [s0,p0,c0] = benchmark({@stayConstrained,@zeroConstrained},options);
        calls0 = statefulMerit('count');
        options.benchmark_id = 'with_report'; options.report_path = fullfile(output_root,'stateful-render.json');
        statefulMerit('reset');
        [s1,p1,c1] = benchmark({@stayConstrained,@zeroConstrained},options);
        assert(isequaln(s0,s1) && isequaln(p0,p1) && isequaln(c0,c1));
        assert(statefulMerit('count') == calls0 && calls0 == 6, 'Reporting added custom-merit calls.');
        report = readReport(options.report_path);
        detail = readDetail(fullfile(output_root,report.plot_data.path));
        merit_panels = detail.plots(arrayfun(@(p) isfield(p, 'channel') && ~isempty(p.channel), detail.plots));
        panel = merit_panels(strcmp({merit_panels.channel},'merit') & strcmp({merit_panels.mode},'raw'));
        assert(numel(report.problems.plot_refs) == 6);
        assert(all(ismember(report.problems.plot_refs,{detail.plots.id})), 'Late scoring overwrote captured history references.');
        assert(strcmp(panel.observation_scope,'rendering_inputs'));
        assert(panel.series(1).mean == 1 && panel.series(2).mean == 2, ...
            'Report did not retain the independently evaluated renderer merit.');
        assert(report.problems.runs(1).merit.best == 4 && report.problems.runs(2).merit.best == 5, ...
            'Scoring facts must stay distinct from stateful renderer values.');
        fprintf('PASS stateful_render: actual renderer merit distinct from scoring, no extra callback, unchanged outputs\n');
    elseif strcmp(slice, 'render')
        options.score_only = false;
        options.savepath = output_root;
        options.benchmark_id = 'render_artifacts';
        options.draw_hist_plots = 'parallel';
        options.report_path = fullfile(output_root, 'render.json');
        setenv('EVAL_REPORT_PUBLIC_OUTPUT', fullfile(output_root, options.benchmark_id));
        [scores, ~, ~] = benchmark({@breakRendering, @zero}, options);
        report = readReport(options.report_path);
        assert(all(isfinite(scores)));
        assert(strcmp(report.stages.numerical.status, 'completed'));
        assert(strcmp(report.stages.scoring.status, 'completed'));
        assert(strcmp(report.stages.rendering.status, 'failed'));
        assert(strcmp(report.status, 'partial'));
        fprintf('PASS render: missing requested history does not invalidate numerical scores\n');
    elseif strcmp(slice, 'archive')
        fixture_root = fileparts(mfilename('fullpath'));
        setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(output_root, 'registry.mat'));
        registerProblemLibrary(struct('name', 'report_archive', 'root', fixture_root, ...
            'select_function', 'eval_fixture_select', 'load_function', 'eval_fixture_load'));
        previous = pwd; cleanup = onCleanup(@() cd(previous)); cd(output_root);
        options = rmfield(options, 'problem');
        options.plibs = {'report_archive'}; options.problem_names = {'shared'};
        options.score_only = false; options.savepath = output_root; options.benchmark_id = 'archive_data';
        options.draw_hist_plots = 'none'; options.max_tol_order = 1;
        options.report_path = fullfile(output_root, 'fresh.json');
        prior_hash_root = getenv('EVAL_REPORT_HASH_FIXTURE_ROOT');
        hash_guard = onCleanup(@() setenv('EVAL_REPORT_HASH_FIXTURE_ROOT', prior_hash_root));
        setenv('EVAL_REPORT_HASH_FIXTURE_ROOT', fullfile(output_root, options.benchmark_id));
        [fresh_scores, ~, ~] = benchmark({@stayWithHashBoundary, @zero}, options);
        fresh = readReport(options.report_path);
        assert(strcmp(fresh.status, 'completed'), 'Ordinary empty log files must not make persistence partial.');
        assert(~isempty(fresh.artifacts), 'Existing produced artifacts must be enumerated.');
        empty_logs = fresh.artifacts(endsWith({fresh.artifacts.path}, 'test_log/log.txt'));
        assert(numel(empty_logs) == 1 && empty_logs.bytes == 0, 'Empty log artifact was omitted.');
        assert(strcmp(empty_logs.sha256, 'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855'), ...
            'Empty file must have the standard SHA256 digest.');
        blocks = fresh.artifacts(endsWith({fresh.artifacts.path}, 'hash-block.bin'));
        assert(numel(blocks) == 1 && blocks.bytes == 1048576);
        assert(strcmp(blocks.sha256, '30e14955ebf1352266dc2ff8067e68104607e750abb9d3b36582b8af909fcb58'), ...
            'A full read block must retain its standard SHA256 digest.');
        [~, caller_name] = fileparts(mfilename('fullpath'));
        assert(any(endsWith({fresh.artifacts.path}, [caller_name, '.m'])), 'Wrapper displaced immediate caller provenance.');
        assert(~isempty(fresh.profiles.plot_refs), 'Actual computed curves must reference shared numeric presentations.');
        assert(any(endsWith({fresh.artifacts.path}, 'data_for_loading.mat')), 'The raw archive must be an enumerated artifact.');
        detail = readDetail(fullfile(output_root, fresh.plot_data.path));
        assert(~isempty(detail.target_work) && numel(detail.target_work(1).axis_values.problem) == 1);
        profile_plots = detail.plots(arrayfun(@(p) ~strcmp(p.kind, 'history'), detail.plots));
        kinds = {profile_plots.kind};
        assert(all(ismember({'performance', 'data', 'log_ratio'}, kinds)), 'All three profile presentations must be captured.');
        channels = {profile_plots.history_or_output};
        assert(isequal(sort(unique(channels)), {'history', 'output'}), 'History- and output-based presentations must both be captured.');
        perf = profile_plots(strcmp(kinds, 'performance') & strcmp(channels, 'history'));
        assert(strcmp(perf.x_transform, 'log2(work/best_work)') && strcmp(perf.observation_scope, 'scoring_profile_work') && strcmp(perf.target_work_ref, 'target-work-1'));
        assert(abs(perf.failure_placeholder - 1.1 * perf.ratio_max) < 1e-12 && perf.std_ddof == 0 && perf.n_runs == 1);
        ratio = profile_plots(strcmp(kinds, 'log_ratio') & strcmp(channels, 'history'));
        assert(strcmp(ratio.problem_mapping, 'bar_sources') && strcmp(ratio.x_transform, 'sorted_bar_position'));
        assert(numel(ratio.bar_sources.problem_index) == numel(ratio.series.x), 'Every sorted bar keeps its problem/run identity.');
        assert(ratio.tie_pairs == sum(ratio.bar_sources.tie) && ratio.both_failed_pairs == sum(ratio.bar_sources.solver1_failed & ratio.bar_sources.solver2_failed), ...
            'Tie and both-failed counts must agree with the retained bar identities.');
        source_files = dir(fullfile(output_root, 'archive_data', '**', 'data_for_loading.mat'));
        source_path = fullfile(source_files(1).folder, source_files(1).name);
        before = hashFile(source_path);
        load_options = struct('load', 'latest', 'benchmark_id', 'archive_data', 'score_only', true, ...
            'silent', true, 'report_path', fullfile(output_root, 'loaded.json'), 'max_tol_order', 1);
        [loaded_scores, ~, ~] = benchmark({@mustNotRun, @mustNotRun}, load_options);
        loaded = readReport(load_options.report_path);
        assert(isequaln(fresh_scores, loaded_scores), 'Load changed score semantics.');
        assert(strcmp(loaded.operation, 'load') && strcmp(loaded.source.sha256, before));
        assert(strcmp(hashFile(source_path), before), 'Loading changed the original archive.');
        assert(~loaded.coverage.original_selection_known && strcmp(loaded.coverage.scope, 'retained_archive'));
        assert(isempty(loaded.coverage.load_failed), 'Archive coverage must not certify zero historical failed loads.');
        assert(isempty(loaded.problems.budget.evaluations) && strcmp(loaded.problems.budget.reason, 'original_execution_budget_not_retained'), 'Current defaults are not original archive budgets.');
        assert(isempty(loaded.problems.runs(1).budget_reached) && isempty(loaded.problems.runs(1).oracle_seed), 'Load cannot infer budgets or seeds.');
        assert(strcmp(loaded.problems.runs(1).oracle_seed_reason, 'execution_metadata_not_retained'));
        assert(isempty(loaded.configuration.effective.feature.name), 'Default load feature is not original execution provenance.');
        assert(strcmp(loaded.coverage.reason, 'original_selection_and_load_failures_not_retained'));
        options.score_only = true; options.problem_names = {'unloadable'};
        options.report_path = fullfile(output_root, 'all_failed.json');
        benchmark(solvers, options);
        failed = readReport(options.report_path);
        assert(strcmp(failed.status, 'failed') && failed.coverage.load_failed == 1);
        fprintf('PASS archive: real artifacts and curves, SHA256, no execution, honest loaded provenance, all loads failed\n');
    elseif strcmp(slice, 'parallel')
        fixture_root = fileparts(mfilename('fullpath'));
        setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(output_root, 'registry.mat'));
        registerProblemLibrary(struct('name', 'report_parallel', 'root', fixture_root, ...
            'select_function', 'eval_fixture_select', 'load_function', 'eval_fixture_load'));
        options = rmfield(options, 'problem'); options.plibs = {'report_parallel'};
        options.max_tol_order = 1;
        [s0,p0,c0] = benchmark(solvers, options);
        options.n_jobs = 2; options.report_path = fullfile(output_root, 'parallel.json');
        [s1,p1,c1] = benchmark({@stayOnWorker, @zeroOnWorker}, options);
        assert(isequaln(s0,s1) && isequaln(p0,p1) && isequaln(c0,c1));
        report = readReport(options.report_path);
        assert(report.coverage.selected == 2 && report.coverage.load_failed == 1);
        assert(strcmp(report.status, 'partial'));
        completed = report.problems(strcmp({report.problems.name}, 'shared'));
        assert(all(~[completed.runs.abnormal_termination]), 'Parallel request fell back to controller execution.');
        fprintf('PASS parallel: parent-owned report, unchanged public numerical results\n');
    elseif strcmp(slice, 'ownership')
        options.report_path = fullfile(output_root, 'replaced-success.json');
        setenv('EVAL_REPORT_PUBLIC_REPLACE_TARGET', options.report_path);
        options.merit_fun = @replaceAndReturn;
        try
            benchmark(solvers, options);
            error('Expected report ownership failure.');
        catch cause
            assert(strcmp(cause.identifier, 'OptiProfiler:EvalReportOwnership'));
        end
        assert(strcmp(fileread(options.report_path), 'external replacement'));
        options.report_path = fullfile(output_root, 'companion-replaced.json');
        companion_path = fullfile(output_root, 'companion-replaced.plot_data.json');
        setenv('EVAL_REPORT_PUBLIC_REPLACE_TARGET', companion_path);
        options.merit_fun = @replaceAndReturn;
        try
            benchmark(solvers, options);
            error('Expected companion ownership failure.');
        catch cause
            assert(strcmp(cause.identifier, 'OptiProfiler:EvalReportOwnership'));
        end
        assert(strcmp(fileread(companion_path), 'external replacement'), 'Foreign companion was overwritten.');
        options.report_path = fullfile(output_root, 'replaced-error.json');
        setenv('EVAL_REPORT_PUBLIC_REPLACE_TARGET', options.report_path);
        options.merit_fun = @replaceAndFail;
        try
            benchmark(solvers, options);
            error('Expected original merit failure.');
        catch cause
            assert(strcmp(cause.identifier, 'MATLAB:benchmark:merit_fun_error'), 'Secondary report failure masked original exception.');
        end
        assert(strcmp(fileread(options.report_path), 'external replacement'));
        fprintf('PASS ownership: foreign replacement untouched, success write error visible, original failure preserved\n');
    elseif strcmp(slice, 'plain_reference')
        fixture_root = fileparts(mfilename('fullpath'));
        setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(output_root, 'registry.mat'));
        registerProblemLibrary(struct('name', 'report_plain', 'root', fixture_root, ...
            'select_function', 'eval_fixture_select', 'load_function', 'eval_fixture_load'));
        prior = getenv('EVAL_REPORT_PLAIN_FAILURE_COUNTER');
        cleanup = onCleanup(@() setenv('EVAL_REPORT_PLAIN_FAILURE_COUNTER', prior));
        setenv('EVAL_REPORT_PLAIN_FAILURE_COUNTER', fullfile(output_root, 'first_load_marker'));
        options = rmfield(options, 'problem'); options.plibs = {'report_plain'};
        options.problem_names = {'shared'}; options.run_plain = true; options.max_tol_order = 1;
        original_error = '';
        try, benchmark(solvers, options); catch cause, original_error = cause.identifier; end
        assert(~isempty(original_error), 'Expected legacy all-plain-loads-failed exception.');
        setenv('EVAL_REPORT_PLAIN_FAILURE_COUNTER', fullfile(output_root, 'second_first_load_marker'));
        options.report_path = fullfile(output_root, 'plain-reference.json');
        reported_error = '';
        try, benchmark(solvers, options); catch cause, reported_error = cause.identifier; end
        assert(strcmp(original_error, reported_error), 'Report changed the legacy exception.');
        report = readReport(options.report_path);
        assert(report.coverage.selected == 1 && report.coverage.completed == 1 && report.coverage.load_failed == 0);
        assert(strcmp(report.status, 'failed'));
        assert(any(strcmp({report.problems.role}, 'plain_reference')));
        fprintf('PASS plain reference: existing fatal behavior preserved, primary rows and failed baseline recorded\n');
        old_mode = getenv('EVAL_REPORT_PLAIN_FAILURE_MODE');
        mode_cleanup = onCleanup(@() setenv('EVAL_REPORT_PLAIN_FAILURE_MODE', old_mode));
        setenv('EVAL_REPORT_PLAIN_FAILURE_MODE', 'partial');
        setenv('EVAL_REPORT_PLAIN_FAILURE_COUNTER', fullfile(output_root, 'partial_first_load_marker'));
        options.problem_names = {'first', 'second'};
        options.report_path = fullfile(output_root, 'plain-reference-partial.json');
        benchmark(solvers, options);
        report = readReport(options.report_path);
        assert(report.coverage.selected == 2 && report.coverage.completed == 2 && report.coverage.load_failed == 0);
        assert(strcmp(report.status, 'partial') && strcmp(report.stages.numerical.status, 'partial'));
        assert(strcmp(report.stages.numerical.reason, 'requested_plain_reference_incomplete'));
        assert(strcmp(report.stages.scoring.status, 'completed'));
        fprintf('PASS plain reference: returned scoring preserved with truthful partial protocol status\n');
    else
        error('Unknown slice.');
    end
end

function x = stay(fun, x0)
    fun(x0);
    x = x0;
end

function x = zero(fun, x0)
    x = zeros(size(x0));
    fun(x);
end

function x = stayWithHashBoundary(fun, x0)
    x = stay(fun, x0);
    folders = dir(fullfile(getenv('EVAL_REPORT_HASH_FIXTURE_ROOT'), '*'));
    folders = folders([folders.isdir] & ~ismember({folders.name}, {'.', '..'}));
    assert(numel(folders) == 1, 'Expected the current benchmark output directory.');
    target = fullfile(folders.folder, folders.name, 'hash-block.bin');
    fid = fopen(target, 'wb'); guard = onCleanup(@() fclose(fid));
    fwrite(fid, zeros(1048576, 1, 'uint8'), 'uint8');
end

function value = undefinedObjective(x)
    value = NaN;
    if x == 0, value = -Inf; end
end

function x = many(fun, x0)
    x = x0;
    for k = 1:100, fun(x); end
end

function value = spikeObjective(x)
    value = x;
    if x == 2, value = 1e6; elseif x == 51, value = -100; elseif x == 73, value = NaN; end
end

function x = longWalk(fun, ~)
    for k = 1:2500, fun(k); end
    x = 2500;
end

function x = stayConstrained(fun,x0,varargin)
    x = stay(fun,x0);
end

function x = zeroConstrained(fun,x0,varargin)
    x = zero(fun,x0);
end

function value = statefulMerit(fun_value,varargin)
    persistent count;
    if ischar(fun_value)
        if strcmp(fun_value,'reset'), count = 0; end
        value = count; return;
    end
    count = count+1; value = count;
end

function value = fatalMerit(varargin)
    error('Audit:ExpectedMeritFailure', 'Deliberate public merit callback failure.');
    value = 0; %#ok<UNRCH>
end

function x = breakRendering(fun, x0)
    x = x0; fun(x);
    folders = dir(fullfile(getenv('EVAL_REPORT_PUBLIC_OUTPUT'), '*'));
    for k = 1:numel(folders)
        if ismember(folders(k).name, {'.', '..'}), continue; end
        target = fullfile(folders(k).folder, folders(k).name, 'history_plots');
        if isfolder(target)
            rmdir(target, 's');
            fid = fopen(target, 'w'); fclose(fid);
        end
    end
end

function value = hashFile(path)
    if ~usejava('jvm')
        [status,~] = system('command -v sha256sum >/dev/null 2>&1');
        command = 'shasum -a 256 '; if status == 0, command = 'sha256sum '; end
        [status, text] = system([command, '< ', '''', strrep(path, '''', '''"''"'''), '''']);
        assert(status == 0, 'Test prerequisite SHA256 command is unavailable.'); value = text(1:64); return;
    end
    digest = java.security.MessageDigest.getInstance('SHA-256');
    fid = fopen(path, 'rb'); guard = onCleanup(@() fclose(fid));
    while ~feof(fid)
        block = fread(fid, 1048576, '*uint8');
        if ~isempty(block), digest.update(block); end
    end
    value = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2).', 1, []));
end

function x = mustNotRun(varargin)
    error('Audit:UnexpectedSolverCall', 'Load must not execute a solver.');
    x = []; %#ok<UNRCH>
end

function value = replaceAndReturn(fun_value, varargin)
    target = getenv('EVAL_REPORT_PUBLIC_REPLACE_TARGET');
    replacement = [target, '.replacement'];
    fid = fopen(replacement, 'w'); fwrite(fid, 'external replacement', 'char'); fclose(fid);
    movefile(replacement, target, 'f');
    value = fun_value;
end

function value = replaceAndFail(fun_value, varargin)
    replaceAndReturn(fun_value, varargin{:});
    error('Audit:ExpectedMeritFailure', 'Deliberate public callback failure.');
    value = []; %#ok<UNRCH>
end

function x = stayOnWorker(fun, x0)
    assert(~isempty(getCurrentTask()), 'Expected actual MATLAB worker execution.');
    x = stay(fun, x0);
end

function x = zeroOnWorker(fun, x0)
    assert(~isempty(getCurrentTask()), 'Expected actual MATLAB worker execution.');
    x = zero(fun, x0);
end

function document = readValidated(path, schema_path)
    errors = evalReportSchemaCheck(path, schema_path);
    assert(isempty(errors), 'Schema violations in %s:\n%s', path, strjoin(errors(1:min(20, numel(errors))), newline));
    document = normalizeCells(jsondecode(fileread(path)));
end

function value = normalizeCells(value)
% jsondecode returns a cell array for a JSON array of objects whose key sets
% differ (optional keys such as provider, rendering or *_reason). Test code
% indexes records as struct arrays, so unify the fields (absent -> []).
    if iscell(value) && ~isempty(value) && all(cellfun(@(v) isstruct(v) && isscalar(v), value))
        names = {};
        for k = 1:numel(value), names = union(names, fieldnames(value{k}), 'stable'); end
        array = repmat(struct(), numel(value), 1);
        for k = 1:numel(value)
            for n = 1:numel(names)
                if isfield(value{k}, names{n}), array(k).(names{n}) = normalizeCells(value{k}.(names{n})); else, array(k).(names{n}) = []; end
            end
        end
        value = array;
    elseif iscell(value)
        value = cellfun(@normalizeCells, value, 'UniformOutput', false);
    elseif isstruct(value)
        names = fieldnames(value);
        for k = 1:numel(value)
            for n = 1:numel(names), value(k).(names{n}) = normalizeCells(value(k).(names{n})); end
        end
    end
end

function value = spikeAt17(x)
    value = x;
    if x == 17, value = 1e6; end
end

function x = walk1003(fun, ~)
    for k = 1:1003, fun(k); end
    x = 1003;
end
