classdef TestFeatureReportV2 < matlab.unittest.TestCase
% Public benchmark receipt contract; runtime tests execute only on syu-ubuntu.
    methods (Test)
        function identityReportSeparatesSpecificationAndExperiment(testCase)
            output = tempname(getenv('OP_ARTIFACTS'));
            mkdir(output);
            feature = Feature('plain+plain');
            problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [2; 1], 'name', 'REPORT_V2_IDENTITY'));
            options = struct('feature', feature, 'n_runs', 3, 'problem', problem, ...
                'solver_names', {{'stay', 'half'}}, 'solver_isrand', [false, true], ...
                'n_jobs', 1, 'max_eval_factor', 2, 'seed', 17, 'score_only', true, ...
                'draw_hist_plots', 'none', 'silent', true, 'savepath', output, ...
                'report_path', fullfile(output, 'identity.json'));
            benchmark({@stay, @half}, options);
            report = jsondecode(fileread(options.report_path));
            testCase.verifyEqual(report.schema, 'optiprofiler.eval_report/2');
            testCase.verifyEqual(report.status, 'completed');
            testCase.verifyEmpty(feature.stages);
            testCase.verifyTrue(isfield(report.configuration.effective, 'experiment'));
            testCase.verifyTrue(isfield(report.configuration.effective.feature, 'stages'));
            if isfield(report.configuration.effective, 'experiment')
                testCase.verifyEqual(report.configuration.effective.experiment.primary.n_runs, 3);
                testCase.verifyEqual(report.configuration.effective.experiment.primary.execution_strategy, 'identity');
                testCase.verifyEmpty(report.configuration.effective.feature.stages);
                testCase.verifyEqual(report.configuration.effective.feature.route, 'feature');
                testCase.verifyEqual(report.configuration.effective.feature.declaration_route, 'feature_name');
            end
        end

        function trustedRefinedIdentityIsAValidReplayInput(testCase)
            native = struct('schema', 'options_refined-v2', 'feature_route', 'feature_name', ...
                'feature_name', 'plain', 'feature_specification', {{struct('name', 'plain', 'options', struct())}}, ...
                'n_runs', 3, 'seed', 17, 'load', 'latest', 'solvers_to_load', [2 1], ...
                'report_path', 'old-report.json', 'savepath', '/obsolete/source');
            [options, receipt] = loadBenchmarkOptions(native);
            testCase.verifyTrue(isa(options.feature, 'Feature'));
            testCase.verifyTrue(options.feature.is_identity);
            testCase.verifyEqual(options.n_runs, 3);
            testCase.verifyEqual(options.seed, 17);
            testCase.verifyFalse(any(isfield(options, {'schema', 'feature_name', 'feature_specification', 'load', 'solvers_to_load', 'report_path', 'savepath'})));
            testCase.verifyEqual(receipt.source_schema, 'options_refined-v2');
            native.feature = Feature('plain+plain');
            native.feature_specification = {struct('name', 'noisy', 'options', struct())};
            testCase.verifyError(@() loadBenchmarkOptions(native), 'OptiProfiler:InconsistentNativeFeature');
            native.feature = struct('name', 'plain');
            testCase.verifyError(@() loadBenchmarkOptions(native), 'OptiProfiler:InvalidNativeFeature');
        end

        function defaultPlainDoesNotInventAnInputKeyword(testCase)
            output = tempname(getenv('OP_ARTIFACTS')); mkdir(output);
            problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [2; 1], 'name', 'REPORT_V2_DEFAULT'));
            options = struct('problem', problem, 'n_runs', 1, 'n_jobs', 1, ...
                'solver_names', {{'stay', 'half'}}, 'max_eval_factor', 2, 'seed', 17, ...
                'score_only', true, 'silent', true, 'report_path', fullfile(output, 'default.json'));
            benchmark({@stay, @half}, options);
            report = jsondecode(fileread(options.report_path));
            testCase.verifyEmpty(report.configuration.effective.feature.route);
            testCase.verifyEqual(report.configuration.effective.feature.effective_name, 'plain');
            testCase.verifyEqual(report.configuration.effective.experiment.primary.n_runs, 1);
        end

        function legacyReplayRequiresIdentityAndPreservesNativeValues(testCase)
            old = struct('n_runs', 3, 'noise_level', 0, 'noise_map', @sin, 'noise_mode', 'deterministic');
            testCase.verifyError(@() loadBenchmarkOptions(old), 'OptiProfiler:LegacyFeatureIdentityMissing');
            [options, receipt] = loadBenchmarkOptions(old, 'noisy');
            testCase.verifyEqual(options.n_runs, 3);
            testCase.verifyEqual(options.feature.stages{1}.options.noise_map, @sin);
            testCase.verifyEmpty(options.feature.declared);
            testCase.verifyEqual(receipt.original_request_origin, 'unknown');
            testCase.verifyEqual(receipt.current_legacy_feature_name_override, 'noisy');
            bad = old; bad.schema = 'options_refined-v999';
            testCase.verifyError(@() loadBenchmarkOptions(bad), 'OptiProfiler:UnsupportedNativeOptionsVersion');
        end

        function nativeProvenanceAndBoundedProjectionNeverRunCallbacks(testCase)
            fixture = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', ...
                'feature-v2', 'native-legacy-enum-feature.mat');
            before = readBytes(fixture); loaded = load(fixture);
            [options, receipt] = loadBenchmarkOptions(struct('feature', loaded.legacy_feature, 'n_runs', 4));
            testCase.verifyEqual(options.n_runs, 4);
            testCase.verifyEqual(receipt.feature_import.n_runs, 3);
            testCase.verifyEmpty(options.feature.declared);
            testCase.verifyEqual(readBytes(fixture), before);
            native_feature = options.feature;
            output = tempname(getenv('OP_ARTIFACTS')); mkdir(output);
            file = fullfile(output, 'canonical.mat');
            save(file, 'native_feature', '-v7'); restored = load(file);
            testCase.verifyEqual(restored.native_feature.stages, native_feature.stages);
            testCase.verifyEmpty(restored.native_feature.declared);

            calls = 0; text = repmat(char(955), 1, 300);
            value = struct('effective_name', repmat('noisy+', 1, 60), 'name', '/private/name', ...
                'full_feature_stamp', text, 'stages', {{struct('options', struct('noise_map', @forbidden))}});
            projected = jsondecode(optiprofiler_internal.EvalReport.encodeMetadata(value));
            testCase.verifyEqual(calls, 0);
            testCase.verifyEmpty(projected.full_feature_stamp);
            testCase.verifyEqual(projected.full_feature_stamp_bytes, 600);
            testCase.verifyEqual(projected.full_feature_stamp_reason, 'omitted_from_bounded_metadata_projection');
            testCase.verifyEmpty(projected.effective_name);
            testCase.verifyEmpty(projected.name);
            testCase.verifyEqual(projected.stages.options.noise_map.kind, 'callback');
            function out = forbidden(varargin) %#ok<INUSD,STOUT>
                calls = calls + 1;
                error('OptiProfiler:MetadataExecutedCallback', 'Metadata must never call a callback.');
            end
        end

        function wholeLibraryRolesNativeReplayAndFilteredLoad(testCase)
            output = tempname(getenv('OP_ARTIFACTS')); mkdir(output);
            fixture = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'feature-report-v2');
            old_path = path; old_cwd = pwd;
            old_registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            old_forbid = getenv('OP_FEATURE_D_FORBID_PROVIDER');
            testCase.addTeardown(@() path(old_path));
            testCase.addTeardown(@() cd(old_cwd));
            testCase.addTeardown(@() setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', old_registry));
            testCase.addTeardown(@() setenv('OP_FEATURE_D_FORBID_PROVIDER', old_forbid));
            addpath(fixture); cd(output);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(output, 'registry.mat'));
            setenv('OP_FEATURE_D_FORBID_PROVIDER', '0');
            registerProblemLibrary(struct('name', 'report_v2', 'root', fixture, ...
                'select_function', 'feature_report_select', 'load_function', 'feature_report_load'));
            calls = [0 0 0];
            solvers = {@(fun, x0) probe(1, fun, x0), @(fun, x0) probe(2, fun, x0), @(fun, x0) probe(3, fun, x0)};
            specification = {struct('name', 'noisy', 'options', struct('noise_level', 0, 'noise_mode', 'deterministic')), ...
                struct('name', 'truncated', 'options', struct('significant_digits', 5))};
            options = struct('feature', Feature(specification), 'plibs', {{'report_v2'}}, ...
                'ptype', 'u', 'mindim', 1, 'maxdim', 3, 'n_runs', 3, 'run_plain', true, ...
                'solver_isrand', [false true false], 'solver_names', {{'stay', 'half', 'stay_again'}}, ...
                'n_jobs', 1, 'max_eval_factor', 2, 'max_tol_order', 1, 'seed', 17, ...
                'score_only', false, 'draw_hist_plots', 'none', 'silent', true, ...
                'benchmark_id', 'fresh', 'savepath', output, ...
                'report_path', fullfile(output, 'fresh.json'));
            [scores, profiles, curves] = benchmark(solvers, options);
            testCase.verifyEqual(calls, [4 8 4]);
            report = jsondecode(fileread(options.report_path));
            testCase.verifyEqual(report.configuration.effective.experiment.primary.n_runs, 3);
            testCase.verifyEqual(report.configuration.effective.experiment.plain_reference.n_runs, 1);
            files = dir(fullfile(output, '**', 'data_for_loading.mat'));
            testCase.assertNumElements(files, 1);
            source_path = fullfile(files(1).folder, files(1).name);
            bytes_before = readBytes(source_path);
            loaded = load(source_path, 'results_plibs');
            group = loaded.results_plibs{1};
            testCase.verifyEqual(size(group.fun_histories, 3), 3);
            testCase.verifyEqual(size(group.results_plib_plain.fun_histories, 3), 1);
            testCase.verifyEqual(group.execution_metadata{1}.real_n_runs(:)', [1 3 1]);
            testCase.verifyEqual(group.results_plib_plain.execution_metadata{1}.real_n_runs(:)', [1 1 1]);
            pipeline = jsondecode(group.feature_pipeline);
            testCase.verifyEqual(pipeline.schema, 'feature_pipeline-v3');
            testCase.verifyEqual({pipeline.feature.stages.identity}, {'noisy#0', 'truncated#0'});
            testCase.verifyFalse(isfield(pipeline.feature.stages(1).options, 'n_runs'));
            [replay, import_receipt] = loadBenchmarkOptions(fullfile(files(1).folder, 'options_refined.mat'));
            testCase.verifyEqual(import_receipt.source_schema, 'options_refined-v2');
            replay.score_only = false; replay.savepath = output; replay.benchmark_id = 'replay';
            replay.report_path = fullfile(output, 'replay.json'); calls = [0 0 0];
            [scores2, profiles2, curves2] = benchmark(solvers, replay);
            testCase.verifyEqual(calls, [4 8 4]);
            testCase.verifyEqual(scores2, scores); testCase.verifyEqual(profiles2, profiles);
            testCase.verifyTrue(isequaln(curves2, curves));
            markers = dir(fullfile(files(1).folder, 'time_stamp_*.txt'));
            testCase.assertNumElements(markers, 1);
            stamp = markers(1).name(12:end-4);
            cd(fileparts(files(1).folder));
            setenv('OP_FEATURE_D_FORBID_PROVIDER', '1'); calls = [0 0 0];
            load_options = struct('load', stamp, 'benchmark_id', '.', 'score_only', true, ...
                'silent', true, 'n_jobs', 1, 'mindim', 3, 'maxdim', 3, 'max_tol_order', 1, ...
                'solvers_to_load', [3 2], 'report_path', fullfile(output, 'load.json'));
            benchmark({@forbidden, @forbidden}, load_options);
            testCase.verifyEqual(calls, [0 0 0]);
            testCase.verifyEqual(readBytes(source_path), bytes_before);
            reloaded = jsondecode(fileread(load_options.report_path));
            testCase.verifyEmpty(fieldnames(reloaded.configuration.effective.experiment));
            testCase.verifyEmpty(reloaded.configuration.effective.feature.stages);
            testCase.verifyEmpty(reloaded.configuration.effective.feature.name);
            testCase.verifyEqual(reloaded.configuration.effective.feature.seed_policy, 'not_applicable_no_solver_execution');
            retained = reloaded.configuration.retained_result_metadata;
            primary = retained(strcmp({retained.role}, 'primary'));
            testCase.verifyEqual(primary.feature_pipeline, pipeline);
            testCase.verifyEqual(primary.retained_problem_indices, 2);
            % Existing MATLAB validation sorts unique solver indices. Check a
            % genuine subset and align retained metadata with that public rule.
            testCase.verifyEqual(primary.retained_solver_indices(:)', [2 3]);
            record = reloaded.problems(strcmp({reloaded.problems.role}, 'primary'));
            testCase.verifyEqual(record.name, 'WIDE');
            runs = record.runs;
            if ~iscell(runs), runs = num2cell(runs); end
            kinds = cellfun(@(x) x.execution.kind, runs, 'UniformOutput', false);
            kinds = kinds(:);
            testCase.verifyEqual(kinds(1:3), {'actual'; 'actual'; 'actual'});
            testCase.verifyEqual(kinds(4:6), {'actual'; 'repeated'; 'repeated'});
            for j = 1:numel(runs)
                if strcmp(runs{j}.execution.kind, 'actual')
                    testCase.verifyEqual(sort(fieldnames(runs{j}.runtime)), ...
                        sort({'language'; 'execution_strategy'; 'runtime_policy'; 'seed_policy'}));
                else
                    testCase.verifyFalse(isfield(runs{j}, 'runtime'));
                end
            end
            % These copies are explicit metadata-negative fixtures, not claimed
            % to be genuine historical producers. The genuine source remains
            % unchanged and the numerical channels are never regenerated.
            payloads = {[], '{not-json', struct('schema', 'feature_pipeline-v999', 'retained_marker', 71), ...
                struct('schema', 'feature_pipeline-v1', 'retained_marker', 72, 'options', struct('n_runs', 5))};
            statuses = {'absent', 'malformed', 'unsupported', 'known'};
            for variant = 1:numel(payloads)
                variant_root = fullfile(output, sprintf('metadata-fixture-%d', variant));
                mkdir(variant_root);
                copyfile(files(1).folder, fullfile(variant_root, 'test_log'));
                results_plibs = loaded.results_plibs;
                if variant == 1
                    results_plibs{1} = rmfield(results_plibs{1}, 'feature_pipeline');
                elseif variant == 2
                    results_plibs{1}.feature_pipeline = payloads{variant};
                else
                    results_plibs{1}.feature_pipeline = jsonencode(payloads{variant});
                end
                variant_source = fullfile(variant_root, 'test_log', 'data_for_loading.mat');
                save(variant_source, 'results_plibs', '-v7.3');
                variant_before = readBytes(variant_source);
                cd(variant_root); calls = [0 0 0];
                load_options.report_path = fullfile(output, sprintf('metadata-load-%d.json', variant));
                benchmark({@forbidden, @forbidden}, load_options);
                testCase.verifyEqual(calls, [0 0 0]);
                testCase.verifyEqual(readBytes(variant_source), variant_before);
                variant_report = jsondecode(fileread(load_options.report_path));
                records = variant_report.configuration.retained_result_metadata;
                record = records(strcmp({records.role}, 'primary'));
                testCase.verifyEqual(record.feature_pipeline_interpretation.status, statuses{variant});
                if variant > 2, testCase.verifyEqual(record.feature_pipeline, payloads{variant}); end
                testCase.verifyEmpty(fieldnames(variant_report.configuration.effective.experiment));
            end
            testCase.verifyEqual(readBytes(source_path), bytes_before);

            function x = probe(index, fun, x0)
                calls(index) = calls(index) + 1; fun(x0);
                x = x0;
                if index == 2, x = x0 / 2; end
                fun(x);
            end
            function x = forbidden(varargin) %#ok<STOUT,INUSD>
                calls(1) = calls(1) + 1;
                error('OptiProfiler:ForbiddenSolver', 'Saved-result load executed a solver.');
            end
        end
    end
end

function x = stay(fun, x0)
    fun(x0); x = x0;
end

function x = half(fun, x0)
    fun(x0); x = x0 / 2; fun(x);
end

function value = readBytes(path)
    fid = fopen(path, 'rb'); cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    value = fread(fid, Inf, '*uint8');
end
