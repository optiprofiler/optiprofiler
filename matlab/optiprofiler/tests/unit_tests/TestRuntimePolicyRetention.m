classdef TestRuntimePolicyRetention < matlab.unittest.TestCase
% New identity and single-stage executions record matlab-legacy-single-v2
% (per-query constraint counters). A reloaded archive reports the runtime
% receipts it retained, verbatim: an archive written under version 1 keeps
% 'matlab-legacy-single-v1'. The version-1 copy below is an explicit metadata
% fixture derived from a genuine run, not a claimed historical producer; its
% numerical channels are never regenerated.

    methods (Test)
        function reloadedArchiveKeepsRecordedPolicy(testCase)
            output = tempname;
            mkdir(output);
            fixture = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'feature-report-v2');
            old_path = path;
            old_cwd = pwd;
            old_registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            old_forbid = getenv('OP_FEATURE_D_FORBID_PROVIDER');
            testCase.addTeardown(@() rmdir(output, 's'));
            testCase.addTeardown(@() path(old_path));
            testCase.addTeardown(@() cd(old_cwd));
            testCase.addTeardown(@() setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', old_registry));
            testCase.addTeardown(@() setenv('OP_FEATURE_D_FORBID_PROVIDER', old_forbid));
            addpath(fixture);
            cd(output);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(output, 'registry.mat'));
            setenv('OP_FEATURE_D_FORBID_PROVIDER', '0');
            registerProblemLibrary(struct('name', 'report_v2', 'root', fixture, ...
                'select_function', 'feature_report_select', 'load_function', 'feature_report_load'));
            feature = Feature({struct('name', 'noisy', 'options', struct('noise_level', 0, 'noise_mode', 'deterministic'))});
            options = struct('feature', feature, 'plibs', {{'report_v2'}}, 'ptype', 'u', 'mindim', 1, 'maxdim', 3, ...
                'n_runs', 1, 'run_plain', true, 'solver_names', {{'stay', 'half'}}, 'solver_isrand', [false false], ...
                'n_jobs', 1, 'max_eval_factor', 2, 'max_tol_order', 1, 'seed', 17, 'score_only', false, ...
                'draw_hist_plots', 'none', 'silent', true, 'benchmark_id', 'fresh', 'savepath', output, ...
                'report_path', fullfile(output, 'fresh.json'));
            benchmark({@stayAtStart, @halfStep}, options);
            fresh = jsondecode(fileread(options.report_path));
            testCase.verifyEqual(unique(TestRuntimePolicyRetention.policies(fresh)), {'matlab-legacy-single-v2'}, ...
                'New identity and single-stage executions record version 2.');

            files = dir(fullfile(output, '**', 'data_for_loading.mat'));
            testCase.assertNumElements(files, 1);
            source_path = fullfile(files(1).folder, files(1).name);
            source_bytes = TestRuntimePolicyRetention.readBytes(source_path);
            loaded = load(source_path, 'results_plibs');
            markers = dir(fullfile(files(1).folder, 'time_stamp_*.txt'));
            testCase.assertNumElements(markers, 1);
            stamp = markers(1).name(12:end-4);

            variant_root = fullfile(output, 'version-1-fixture');
            mkdir(variant_root);
            copyfile(files(1).folder, fullfile(variant_root, 'test_log'));
            results_plibs = loaded.results_plibs;
            results_plibs{1} = TestRuntimePolicyRetention.markVersionOne(results_plibs{1});
            if isfield(results_plibs{1}, 'results_plib_plain')
                results_plibs{1}.results_plib_plain = TestRuntimePolicyRetention.markVersionOne(results_plibs{1}.results_plib_plain);
            end
            save(fullfile(variant_root, 'test_log', 'data_for_loading.mat'), 'results_plibs', '-v7.3');

            setenv('OP_FEATURE_D_FORBID_PROVIDER', '1');
            cd(variant_root);
            load_options = struct('load', stamp, 'benchmark_id', '.', 'score_only', true, 'silent', true, ...
                'n_jobs', 1, 'max_tol_order', 1, 'report_path', fullfile(output, 'load-version-1.json'));
            benchmark({@forbiddenSolver, @forbiddenSolver}, load_options);
            reloaded = jsondecode(fileread(load_options.report_path));
            testCase.verifyEqual(unique(TestRuntimePolicyRetention.policies(reloaded)), {'matlab-legacy-single-v1'}, ...
                'A reloaded archive must keep the policy string it recorded.');

            % Control: the unmodified archive reloads with the policy it recorded.
            cd(fileparts(files(1).folder));
            load_options.report_path = fullfile(output, 'load-version-2.json');
            benchmark({@forbiddenSolver, @forbiddenSolver}, load_options);
            control = jsondecode(fileread(load_options.report_path));
            testCase.verifyEqual(unique(TestRuntimePolicyRetention.policies(control)), {'matlab-legacy-single-v2'});
            testCase.verifyEqual(TestRuntimePolicyRetention.readBytes(source_path), source_bytes);
        end
    end

    methods (Static)
        function group = markVersionOne(group)
            for p = 1:numel(group.execution_metadata)
                record = group.execution_metadata{p};
                if ~isstruct(record) || ~isfield(record, 'runtime_receipts')
                    continue;
                end
                for k = 1:numel(record.runtime_receipts)
                    if ~isempty(record.runtime_receipts{k})
                        record.runtime_receipts{k}.runtime_policy = 'matlab-legacy-single-v1';
                    end
                end
                group.execution_metadata{p} = record;
            end
        end

        function values = policies(report)
            values = {};
            problems = report.problems;
            if ~iscell(problems), problems = num2cell(problems); end
            for i = 1:numel(problems)
                if ~isfield(problems{i}, 'runs'), continue; end
                runs = problems{i}.runs;
                if ~iscell(runs), runs = num2cell(runs); end
                for j = 1:numel(runs)
                    if isfield(runs{j}, 'runtime')
                        values{end + 1} = runs{j}.runtime.runtime_policy; %#ok<AGROW>
                    end
                end
            end
        end

        function value = readBytes(path)
            fid = fopen(path, 'r');
            value = fread(fid, Inf, '*uint8');
            fclose(fid);
        end
    end
end

function x = stayAtStart(fun, x0)
    fun(x0);
    x = x0;
end

function x = halfStep(fun, x0)
    fun(x0);
    x = x0 / 2;
    fun(x);
end

function x = forbiddenSolver(varargin) %#ok<STOUT>
    error('OptiProfiler:ForbiddenSolver', 'Saved-result load executed a solver.');
end
