classdef TestMeritReloadConsistency < matlab.unittest.TestCase
%TESTMERITRELOADCONSISTENCY Cached merits, reference targets and raw archives stay in one merit space.
%   Public save/load contracts mirrored by Python's test_merit_reload_consistency.py:
%   - a load without merit_fun never compares cached custom merits with
%     reference targets converted by the default merit;
%   - archives without explicit references keep the V1 behaviour (cached
%     merits are reused);
%   - a merit failure after solving keeps a reloadable raw archive, and a
%     reference that the selected merit maps to a non-finite target fails
%     before any solver call (primary role only).
    properties
        OldCwd
        OldPath
        OldRegistry
        WorkDir
    end
    methods (TestMethodSetup)
        function setupState(testCase)
            testCase.OldCwd = pwd;
            testCase.OldPath = path;
            testCase.OldRegistry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            testCase.WorkDir = tempname;
            mkdir(testCase.WorkDir);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(testCase.WorkDir, 'registry.mat'));
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
            cd(testCase.WorkDir);
            setappdata(0, 'OptiProfilerMeritReloadSolverCalls', 0);
        end
    end
    methods (TestMethodTeardown)
        function restoreState(testCase)
            cd(testCase.OldCwd);
            path(testCase.OldPath);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', testCase.OldRegistry);
            if isappdata(0, 'OptiProfilerMeritReloadSolverCalls')
                rmappdata(0, 'OptiProfilerMeritReloadSolverCalls');
            end
            if isfolder(testCase.WorkDir), rmdir(testCase.WorkDir, 's'); end
        end
    end
    methods (Static)
        function x = finish(fun, x0, varargin)
            TestMeritReloadConsistency.countCall();
            x = 0.5 * ones(size(x0));
            fun(x);
        end
        function x = stay(fun, x0, varargin)
            TestMeritReloadConsistency.countCall();
            fun(x0);
            x = x0;
        end
        function countCall()
            setappdata(0, 'OptiProfilerMeritReloadSolverCalls', ...
                getappdata(0, 'OptiProfilerMeritReloadSolverCalls') + 1);
        end
        function value = squareMerit(f, v, ~)
            value = f * f + v;
        end
        function value = affineMerit(f, v, ~)
            value = 2 * f + 100 + v;
        end
        function value = failingMerit(~, ~, ~) %#ok<STOUT>
            error('TestMeritReloadConsistency:deliberate', 'deliberate merit failure');
        end
        function writeFile(file, code)
            fid = fopen(file, 'w');
            assert(fid >= 0);
            cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
            fprintf(fid, '%s', code);
        end
        function library(root, name, problem_code)
            % problem_code builds the struct `s` for problem `name` (P1/P2 or B1).
            mkdir(root);
            TestMeritReloadConsistency.writeFile(fullfile(root, [name, '_select.m']), sprintf( ...
                'function names = %s_select(options)\n names = {''P1'', ''P2''};\n if contains(''%s'', ''bounded''), names = {''B1''}; end\nend\n', ...
                name, name));
            TestMeritReloadConsistency.writeFile(fullfile(root, [name, '_load.m']), sprintf( ...
                'function problem = %s_load(name)\n%s\nproblem = Problem(s);\nend\n', name, problem_code));
            registerProblemLibrary(struct('name', name, 'root', root, ...
                'select_function', [name, '_select'], 'load_function', [name, '_load']));
        end
        function options = runOptions(name, id, ptype)
            options = struct('plibs', {{name}}, 'ptype', ptype, 'mindim', 1, 'maxdim', 1, ...
                'n_jobs', 1, 'max_eval_factor', 2, 'max_tol_order', 3, 'draw_hist_plots', 'none', ...
                'silent', true, 'benchmark_id', id, 'seed', 7);
        end
        function options = loadOptions(id)
            options = struct('load', 'latest', 'benchmark_id', id, 'score_only', true, ...
                'silent', true, 'max_tol_order', 3, 'max_eval_factor', 2);
        end
        function bytes = archiveBytes(folder)
            files = dir(fullfile(folder, '**', 'data_for_loading.mat'));
            assert(numel(files) == 1, 'Expected exactly one saved archive.');
            fid = fopen(fullfile(files(1).folder, files(1).name), 'r');
            cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
            bytes = fread(fid, Inf, '*uint8');
        end
        function codes = diagnosticCodes(report_path)
            report = jsondecode(fileread(report_path));
            codes = {};
            if ~isfield(report, 'diagnostics') || isempty(report.diagnostics), return; end
            entries = report.diagnostics;
            for i = 1:numel(entries)
                if iscell(entries), entry = entries{i}; else, entry = entries(i); end
                codes{end+1} = entry.code; %#ok<AGROW>
            end
        end
        function code = feasibleProblems()
            code = ['s = struct(''fun'', @(x) sum(x.^2), ''x0'', 2 + strcmp(name, ''P2''), ''name'', name);', newline, ...
                's.reference = struct(''fun'', 0, ''maxcv'', 0, ''kind'', ''feasible'', ''source'', ''analytic test'');'];
        end
    end
    methods (Test)
        function reloadWithoutMeritFunRescoresReferenceArchives(testCase)
            TestMeritReloadConsistency.library(fullfile(testCase.WorkDir, 'lib'), 'meritreload', ...
                TestMeritReloadConsistency.feasibleProblems());
            solvers = {@TestMeritReloadConsistency.finish, @TestMeritReloadConsistency.stay};
            for feature_case = {{'plain', 1}, {'perturbed_x0+noisy', 2}}
                feature_name = feature_case{1}{1};
                custom_id = ['custom-', strrep(feature_name, '+', '-')];
                default_id = ['default-', strrep(feature_name, '+', '-')];
                options = TestMeritReloadConsistency.runOptions('meritreload', custom_id, 'u');
                options.feature_name = feature_name;
                options.n_runs = feature_case{1}{2};
                options.run_plain = true;
                options.merit_fun = @TestMeritReloadConsistency.squareMerit;
                [custom_scores, custom_profiles] = benchmark(solvers, options);
                options = rmfield(options, 'merit_fun');
                options.benchmark_id = default_id;
                [default_scores, default_profiles] = benchmark(solvers, options);
                source = TestMeritReloadConsistency.archiveBytes(fullfile(testCase.WorkDir, custom_id));

                load_options = TestMeritReloadConsistency.loadOptions(custom_id);
                load_options.run_plain = true;
                load_options.report_path = fullfile(testCase.WorkDir, [custom_id, '-reload.json']);
                [scores, profiles] = benchmark({}, load_options);
                % Explicit references are converted with the default merit, so
                % every cached quantity is rescored with it too.
                testCase.verifyEqual(scores, default_scores);
                testCase.verifyEqual(profiles, default_profiles);
                testCase.verifyTrue(any(strcmp(TestMeritReloadConsistency.diagnosticCodes( ...
                    load_options.report_path), 'cached_merits_recomputed')));

                load_options = rmfield(load_options, 'report_path');
                load_options.merit_fun = @TestMeritReloadConsistency.squareMerit;
                [scores, profiles] = benchmark({}, load_options);
                testCase.verifyEqual(scores, custom_scores);
                testCase.verifyEqual(profiles, custom_profiles);
                testCase.verifyEqual(TestMeritReloadConsistency.archiveBytes( ...
                    fullfile(testCase.WorkDir, custom_id)), source);
            end
        end
        function defaultMeritReferenceArchiveReloadsUnchanged(testCase)
            TestMeritReloadConsistency.library(fullfile(testCase.WorkDir, 'lib'), 'meritreload', ...
                TestMeritReloadConsistency.feasibleProblems());
            solvers = {@TestMeritReloadConsistency.finish, @TestMeritReloadConsistency.stay};
            [scores, profiles] = benchmark(solvers, TestMeritReloadConsistency.runOptions('meritreload', 'default', 'u'));
            load_options = TestMeritReloadConsistency.loadOptions('default');
            load_options.report_path = fullfile(testCase.WorkDir, 'reload.json');
            [loaded_scores, loaded_profiles] = benchmark({}, load_options);
            testCase.verifyEqual(loaded_scores, scores);
            testCase.verifyEqual(loaded_profiles, profiles);
            testCase.verifyFalse(any(strcmp(TestMeritReloadConsistency.diagnosticCodes( ...
                load_options.report_path), 'cached_merits_recomputed')));
        end
        function reloadWithoutReferencesKeepsCachedMerits(testCase)
            TestMeritReloadConsistency.library(fullfile(testCase.WorkDir, 'lib'), 'meritreload', ...
                's = struct(''fun'', @(x) sum(x.^2), ''x0'', 2 + strcmp(name, ''P2''), ''name'', name);');
            solvers = {@TestMeritReloadConsistency.finish, @TestMeritReloadConsistency.stay};
            options = TestMeritReloadConsistency.runOptions('meritreload', 'custom', 'u');
            options.run_plain = true;
            options.merit_fun = @TestMeritReloadConsistency.squareMerit;
            [scores, profiles] = benchmark(solvers, options);
            load_options = TestMeritReloadConsistency.loadOptions('custom');
            load_options.run_plain = true;
            [loaded_scores, loaded_profiles] = benchmark({}, load_options);
            testCase.verifyEqual(loaded_scores, scores);
            testCase.verifyEqual(loaded_profiles, profiles);
        end
        function failingMeritKeepsReloadableRawArchive(testCase)
            TestMeritReloadConsistency.library(fullfile(testCase.WorkDir, 'lib'), 'meritreload', ...
                's = struct(''fun'', @(x) sum(x.^2), ''x0'', 2 + strcmp(name, ''P2''), ''name'', name);');
            solvers = {@TestMeritReloadConsistency.finish, @TestMeritReloadConsistency.stay};
            options = TestMeritReloadConsistency.runOptions('meritreload', 'failed', 'u');
            options.merit_fun = @TestMeritReloadConsistency.failingMerit;
            testCase.verifyError(@() benchmark(solvers, options), 'MATLAB:benchmark:merit_fun_error');
            testCase.verifyGreaterThan(getappdata(0, 'OptiProfilerMeritReloadSolverCalls'), 0);
            testCase.verifyNotEmpty(dir(fullfile(testCase.WorkDir, 'failed', '**', 'time_stamp_*.txt')));
            source = TestMeritReloadConsistency.archiveBytes(fullfile(testCase.WorkDir, 'failed'));

            options.merit_fun = @TestMeritReloadConsistency.affineMerit;
            options.benchmark_id = 'fresh';
            [fresh_scores, fresh_profiles] = benchmark(solvers, options);
            load_options = TestMeritReloadConsistency.loadOptions('failed');
            load_options.merit_fun = @TestMeritReloadConsistency.affineMerit;
            [scores, profiles] = benchmark({}, load_options);
            testCase.verifyEqual(scores, fresh_scores);
            testCase.verifyEqual(profiles, fresh_profiles);

            options = rmfield(options, 'merit_fun');
            options.benchmark_id = 'fresh-default';
            [default_scores, default_profiles] = benchmark(solvers, options);
            [scores, profiles] = benchmark({}, TestMeritReloadConsistency.loadOptions('failed'));
            testCase.verifyEqual(scores, default_scores);
            testCase.verifyEqual(profiles, default_profiles);
            testCase.verifyEqual(TestMeritReloadConsistency.archiveBytes(fullfile(testCase.WorkDir, 'failed')), source);
        end
        function incompatibleReferenceFailsBeforeSolverCalls(testCase)
            TestMeritReloadConsistency.library(fullfile(testCase.WorkDir, 'lib'), 'meritbounded', ...
                ['s = struct(''fun'', @(x) sum(x.^2), ''x0'', 0.5, ''xl'', 0, ''xu'', 1, ''name'', name);', newline, ...
                's.reference = struct(''fun'', -3, ''maxcv'', 0.5, ''kind'', ''merit'', ''source'', ''best known infeasible point'');']);
            solvers = {@TestMeritReloadConsistency.finish, @TestMeritReloadConsistency.stay};
            options = TestMeritReloadConsistency.runOptions('meritbounded', 'preflight', 'b');
            try
                benchmark(solvers, options);
                testCase.verifyFail('An incompatible reference must stop the benchmark.');
            catch error_info
                testCase.verifyEqual(error_info.identifier, 'MATLAB:applyReferenceTargets:invalidMerit');
                testCase.verifySubstring(error_info.message, 'incompatible with the selected merit');
                testCase.verifySubstring(error_info.message, 'finite real scalar');
                testCase.verifySubstring(error_info.message, 'B1');
            end
            testCase.verifyEqual(getappdata(0, 'OptiProfilerMeritReloadSolverCalls'), 0);
            testCase.verifyEmpty(dir(fullfile(testCase.WorkDir, 'preflight', '**', 'data_for_loading.mat')));

            % A merit that is finite at the reference scores the same problem.
            options.benchmark_id = 'compatible';
            options.merit_fun = @TestMeritReloadConsistency.affineMerit;
            scores = benchmark(solvers, options);
            testCase.verifyTrue(all(isfinite(scores)));
        end
        function plainBaselineNeverTriggersPreflight(testCase)
            % Primary runs start from a strongly infeasible perturbed point, where
            % the default merit is finite at the reference; the plain baseline
            % starts from the feasible x0, where it is not, but plain records are
            % never targets.
            TestMeritReloadConsistency.library(fullfile(testCase.WorkDir, 'lib'), 'meritbounded', ...
                ['s = struct(''fun'', @(x) sum(x.^2), ''x0'', 0.5, ''xl'', 0, ''xu'', 1, ''name'', name);', newline, ...
                's.reference = struct(''fun'', -3, ''maxcv'', 0.5, ''kind'', ''merit'', ''source'', ''best known infeasible point'');']);
            solvers = {@TestMeritReloadConsistency.finish, @TestMeritReloadConsistency.stay};
            options = TestMeritReloadConsistency.runOptions('meritbounded', 'plain-baseline', 'b');
            options.feature_name = 'perturbed_x0';
            options.perturbation_level = 50;
            options.n_runs = 1;
            options.run_plain = true;
            scores = benchmark(solvers, options);
            testCase.verifyTrue(all(isfinite(scores)));
        end
    end
end
