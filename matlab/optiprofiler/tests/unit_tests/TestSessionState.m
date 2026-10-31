classdef TestSessionState < matlab.unittest.TestCase
%TESTSESSIONSTATE benchmark restores the diary, root defaults and library path it changes.
    properties
        WorkDir
        Root
    end
    methods (TestMethodSetup)
        function isolate(testCase)
            old_directory = pwd;
            old_path = path;
            old_registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            old_diary = get(0, 'Diary');
            old_diary_file = get(0, 'DiaryFile');
            testCase.WorkDir = tempname;
            mkdir(testCase.WorkDir);
            % Teardowns run in reverse order: leave and close the folder first.
            testCase.addTeardown(@() rmdir(testCase.WorkDir, 's'));
            testCase.addTeardown(@() TestSessionState.restore(old_directory, old_path, old_registry, old_diary, old_diary_file));
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(testCase.WorkDir, 'registry.mat'));
            testCase.Root = fullfile(testCase.WorkDir, 'lib');
            mkdir(testCase.Root);
            TestSessionState.writeFile(fullfile(testCase.Root, 'sessionlib_select.m'), ...
                sprintf('function names = sessionlib_select(~)\n names = {''Q1'', ''Q2''};\nend\n'));
            TestSessionState.writeFile(fullfile(testCase.Root, 'sessionlib_load.m'), sprintf(['function p = sessionlib_load(name)\n', ...
                ' p = Problem(struct(''fun'', @(x) sum((x - strcmp(name, ''Q2'')).^2), ''x0'', 2, ''name'', name));\nend\n']));
            registerProblemLibrary(struct('name', 'sessionlib', 'root', testCase.Root, ...
                'select_function', 'sessionlib_select', 'load_function', 'sessionlib_load'));
            cd(testCase.WorkDir);
            set(groot, 'DefaultAxesFontSize', 'remove');
            set(groot, 'DefaultAxesFontName', 'remove');
            testCase.addTeardown(@() set(groot, 'DefaultAxesFontSize', 'remove'));
            testCase.addTeardown(@() set(groot, 'DefaultAxesFontName', 'remove'));
        end
    end
    methods (Static)
        function restore(directory, old_path, registry, diary_state, diary_file)
            diary('off');
            set(0, 'DiaryFile', diary_file);
            if strcmp(diary_state, 'on'), diary('on'); end
            cd(directory);
            path(old_path);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', registry);
        end
        function writeFile(file, text)
            fid = fopen(file, 'w');
            fwrite(fid, text, 'char');
            fclose(fid);
        end
        function options = runOptions(id, varargin)
            options = struct('plibs', {{'sessionlib'}}, 'mindim', 1, 'maxdim', 1, 'n_jobs', 1, ...
                'max_eval_factor', 2, 'max_tol_order', 1, 'draw_hist_plots', 'none', 'silent', true, ...
                'benchmark_id', id, 'summarize_performance_profiles', false, ...
                'summarize_data_profiles', false, 'summarize_output_based_profiles', false);
            for k = 1:2:numel(varargin)
                options.(varargin{k}) = varargin{k + 1};
            end
        end
        function x = stay(fun, x0)
            fun(x0);
            x = x0;
        end
        function x = step(fun, x0)
            fun(x0);
            x = x0 - 1;
            fun(x);
        end
        function x = extendPath(fun, x0)
            % A solver may change the path itself; that entry is not benchmark's.
            folder = getappdata(0, 'OptiProfilerSessionSolverFolder');
            if ~isempty(folder), addpath(folder); end
            fun(x0);
            x = x0;
        end
        function value = failingMerit(~, ~, ~) %#ok<STOUT>
            error('TestSessionState:merit', 'deliberate merit failure');
        end
        function yes = onPath(folder)
            yes = any(strcmp(strsplit(path, pathsep), folder));
        end
        function runEither(solvers, options)
            % Success, or the merit failure raised after the solves.
            try
                benchmark(solvers, options);
            catch cause
                assert(strcmp(cause.identifier, 'MATLAB:benchmark:merit_fun_error'), cause.message);
            end
        end
    end
    methods (Test)
        function diaryIsRestoredOnSuccessAndError(testCase)
            solvers = {@TestSessionState.stay, @TestSessionState.step};
            for failing = [false, true]
                options = TestSessionState.runOptions(sprintf('diary%d', failing));
                if failing, options.merit_fun = @TestSessionState.failingMerit; end
                user_diary = fullfile(testCase.WorkDir, sprintf('user%d.txt', failing));
                diary(user_diary);
                diary('on');
                fprintf('before benchmark\n');
                TestSessionState.runEither(solvers, options);
                testCase.verifyEqual(get(0, 'Diary'), 'on');
                testCase.verifyEqual(get(0, 'DiaryFile'), user_diary);
                fprintf('after benchmark\n');
                diary('off');
                text = fileread(user_diary);
                testCase.verifySubstring(text, 'before benchmark');
                testCase.verifySubstring(text, 'after benchmark');
                logs = dir(fullfile(testCase.WorkDir, options.benchmark_id, '**', 'log.txt'));
                testCase.verifyNotEmpty(logs);
                testCase.verifyFalse(contains(fileread(fullfile(logs(1).folder, logs(1).name)), 'after benchmark'));
                % A caller whose diary is off keeps it off, with its file name.
                other = fullfile(testCase.WorkDir, sprintf('off%d.txt', failing));
                set(0, 'DiaryFile', other);
                options.benchmark_id = sprintf('diary-off%d', failing);
                TestSessionState.runEither(solvers, options);
                testCase.verifyEqual(get(0, 'Diary'), 'off');
                testCase.verifyEqual(get(0, 'DiaryFile'), other);
            end
        end

        function rootDefaultsAreRestoredOnSuccessAndError(testCase)
            solvers = {@TestSessionState.stay, @TestSessionState.step};
            for failing = [false, true]
                options = TestSessionState.runOptions(sprintf('defaults%d', failing), 'score_only', true);
                if failing, options.merit_fun = @TestSessionState.failingMerit; end
                set(groot, 'DefaultAxesFontSize', 'remove');
                set(groot, 'DefaultAxesFontName', 'remove');
                TestSessionState.runEither(solvers, options);
                defaults = get(groot, 'default');
                testCase.verifyFalse(any(isfield(defaults, {'defaultAxesFontSize', 'defaultAxesFontName'})));
                set(groot, 'DefaultAxesFontSize', 17);
                set(groot, 'DefaultAxesFontName', 'Courier');
                TestSessionState.runEither(solvers, options);
                testCase.verifyEqual(get(groot, 'DefaultAxesFontSize'), 17);
                testCase.verifyEqual(get(groot, 'DefaultAxesFontName'), 'Courier');
            end
        end

        function libraryRootsAddedByBenchmarkAreRemoved(testCase)
            solvers = {@TestSessionState.extendPath, @TestSessionState.step};
            solver_folder = fullfile(testCase.WorkDir, 'solverpath');
            mkdir(solver_folder);
            setappdata(0, 'OptiProfilerSessionSolverFolder', solver_folder);
            testCase.addTeardown(@() rmappdata(0, 'OptiProfilerSessionSolverFolder'));
            for failing = [false, true]
                options = TestSessionState.runOptions(sprintf('path%d', failing), 'score_only', true);
                if failing, options.merit_fun = @TestSessionState.failingMerit; end
                rmpath(testCase.Root);
                if TestSessionState.onPath(solver_folder), rmpath(solver_folder); end
                TestSessionState.runEither(solvers, options);
                testCase.verifyFalse(TestSessionState.onPath(testCase.Root));
                testCase.verifyTrue(TestSessionState.onPath(solver_folder));
                % A root already on the caller's path stays there.
                addpath(testCase.Root);
                TestSessionState.runEither(solvers, options);
                testCase.verifyTrue(TestSessionState.onPath(testCase.Root));
            end
        end
    end
end
