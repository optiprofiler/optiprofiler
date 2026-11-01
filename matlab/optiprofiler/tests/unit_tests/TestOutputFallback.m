classdef TestOutputFallback < matlab.unittest.TestCase
%TESTOUTPUTFALLBACK A run without an output never passes the output-based test.
%   Mirrors Python's test_output_fallback.py.
    properties (Constant)
        History = 1
        Output = 2
    end
    properties (TestParameter)
        kind = {'loose_reference', 'infinite_initial_merit'};
        noOutput = struct('crashBefore', @TestOutputFallback.crashBefore, ...
            'crashAfter', @TestOutputFallback.crashAfter, ...
            'returnsEmpty', @TestOutputFallback.returnsEmpty);
    end
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
            rmappdata_if_present();
        end
    end
    methods (TestMethodTeardown)
        function restoreState(testCase)
            cd(testCase.OldCwd);
            path(testCase.OldPath);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', testCase.OldRegistry);
            rmappdata_if_present();
            if isfolder(testCase.WorkDir), rmdir(testCase.WorkDir, 's'); end
        end
    end
    methods (Static)
        function x = reachesZero(fun, x0, varargin)
            fun(x0);
            x = zeros(size(x0));
            fun(x);
        end
        function x = crashBefore(~, ~, varargin) %#ok<STOUT>
            error('TestOutputFallback:Crash', 'crashes before its first evaluation');
        end
        function x = crashAfter(fun, x0, varargin) %#ok<STOUT>
            fun(x0);
            fun(zeros(size(x0)));
            error('TestOutputFallback:Crash', 'crashes after reaching the best point');
        end
        function x = returnsEmpty(fun, x0, varargin)
            fun(x0);
            x = [];
        end
        function x = returnsX0(fun, x0, varargin)
            fun(x0);
            x = x0;
        end
        function x = returnsX0WithoutEvaluating(~, x0, varargin)
            x = x0;
        end
        function x = crashesOnProblemB(fun, x0, varargin)
            if x0(1) == 3
                error('TestOutputFallback:Crash', 'crashes on one problem');
            end
            x = TestOutputFallback.reachesZero(fun, x0);
        end
        function x = crashesOnEachFirstRun(fun, x0, varargin)
            seen = getappdata(0, 'TestOutputFallbackSeen');
            if ~ismember(x0(1), seen)
                setappdata(0, 'TestOutputFallbackSeen', [seen, x0(1)]);
                error('TestOutputFallback:Crash', 'crashes in the first run of a problem');
            end
            x = TestOutputFallback.reachesZero(fun, x0);
        end
        function value = infiniteAwayFromZero(x)
            value = sum(x.^2);
            if x(1) > 1, value = Inf; end
        end
        function writeFile(file, code)
            fid = fopen(file, 'w');
            assert(fid >= 0);
            cleanup = onCleanup(@() fclose(fid));
            fprintf(fid, '%s', code);
        end
    end
    methods (Test)
        function runWithoutAnOutputNeverPassesTheOutputTest(testCase, kind, noOutput)
            options = testCase.library(kind);
            options.score_only = true;
            [~, profiles] = benchmark({@TestOutputFallback.reachesZero, noOutput}, options);
            testCase.verifyEqual(squeeze(profiles(1, 1, testCase.Output, :)), ones(3, 1));
            testCase.verifyEqual(squeeze(profiles(2, 1, testCase.Output, :)), zeros(3, 1));
        end

        function solverThatReturnsX0KeepsPassing(testCase, kind)
            options = testCase.library(kind);
            options.score_only = true;
            [~, profiles] = benchmark({@TestOutputFallback.reachesZero, @TestOutputFallback.returnsX0}, options);
            testCase.verifyEqual(squeeze(profiles(2, 1, testCase.Output, :)), ones(3, 1));
        end

        function outputWithoutEvaluationsIsStillChargedOneUnit(testCase, kind)
            options = testCase.library(kind);
            options.score_only = true;
            [~, profiles] = benchmark({@TestOutputFallback.reachesZero, @TestOutputFallback.returnsX0WithoutEvaluating}, options);
            % One evaluation against two: the genuine output is the fastest.
            testCase.verifyEqual(squeeze(profiles(2, 1, testCase.Output, :)), ones(3, 1));
            testCase.verifyTrue(all(profiles(1, 1, testCase.Output, :) < 1));
        end

        function onlyTheFlaggedRunLosesItsOutput(testCase)
            options = testCase.library('loose_reference');
            options.score_only = true;
            options.feature_name = 'noisy';
            options.n_runs = 2;
            options.report_path = fullfile(testCase.WorkDir, 'report.json');
            benchmark({@TestOutputFallback.reachesZero, @TestOutputFallback.crashesOnEachFirstRun}, options);
            report = jsondecode(fileread(options.report_path));
            convergence = report.profiles.convergence(1);
            testCase.verifyEqual([convergence.output.hits], [4, 2]);
            testCase.verifyEqual([convergence.output.total], [4, 4]);
            testCase.verifyEqual([convergence.history.hits], [4, 2]);
            flags = false(1, 0);
            for p = 1:numel(report.problems)
                runs = report.problems(p).runs;
                runs = runs([runs.solver_index] == 2);
                flags = [flags, [runs.output_fallback]]; %#ok<AGROW>
            end
            testCase.verifyEqual(flags, [true, false, true, false]);
        end

        function solverFreeLoadReproducesTheScores(testCase)
            options = testCase.library('loose_reference');
            options.benchmark_id = 'saved';
            options.draw_hist_plots = 'none';
            [scores, profiles] = benchmark({@TestOutputFallback.reachesZero, @TestOutputFallback.crashAfter}, options);
            testCase.verifyEqual(squeeze(profiles(2, 1, testCase.Output, :)), zeros(3, 1));
            [loaded_scores, loaded_profiles] = benchmark({}, testCase.loadOptions());
            testCase.verifyEqual(loaded_scores, scores);
            testCase.verifyEqual(loaded_profiles, profiles);
        end

        function loadSlicesTheFlagsWithTheSolversAndProblems(testCase)
            options = testCase.library('loose_reference');
            options.benchmark_id = 'saved';
            options.draw_hist_plots = 'none';
            benchmark({@TestOutputFallback.reachesZero, @TestOutputFallback.crashAfter, ...
                @TestOutputFallback.crashesOnProblemB}, options);
            load_options = testCase.loadOptions();
            load_options.solvers_to_load = [1, 2];
            [~, selected] = benchmark({}, load_options);
            testCase.verifyEqual(squeeze(selected(2, 1, testCase.Output, :)), zeros(3, 1));
            load_options.solvers_to_load = [1, 3];
            [~, selected] = benchmark({}, load_options);
            testCase.verifyTrue(all(selected(2, 1, testCase.Output, 1:2) > 0));
            testCase.verifyTrue(all(selected(2, 1, testCase.Output, 1:2) < 1));
            load_options.problem_names = {'A'};
            [~, selected] = benchmark({}, load_options);
            testCase.verifyEqual(squeeze(selected(2, 1, testCase.Output, 1:2)), ones(2, 1));
            load_options.problem_names = {'B'};
            [~, selected] = benchmark({}, load_options);
            testCase.verifyEqual(squeeze(selected(2, 1, testCase.Output, 1:2)), zeros(2, 1));
        end

        function archiveWithoutTheFlagsKeepsItsBehaviour(testCase)
            options = testCase.library('loose_reference');
            options.benchmark_id = 'saved';
            options.draw_hist_plots = 'none';
            [~, profiles] = benchmark({@TestOutputFallback.reachesZero, @TestOutputFallback.returnsEmpty}, options);
            files = dir(fullfile(testCase.WorkDir, 'saved', '**', 'data_for_loading.mat'));
            testCase.assertEqual(numel(files), 1);
            file = fullfile(files(1).folder, files(1).name);
            stored = load(file);
            testCase.assertTrue(isfield(stored.results_plibs{1}, 'solver_output_fallbacks'));
            stored.results_plibs{1} = rmfield(stored.results_plibs{1}, 'solver_output_fallbacks');
            save(file, '-struct', 'stored');
            [~, legacy] = benchmark({}, testCase.loadOptions());
            % Without the flags the recorded value decides, as in an old archive ...
            testCase.verifyEqual(squeeze(legacy(2, 1, testCase.Output, :)), ones(3, 1));
            % ... and nothing that is history-based depends on them.
            testCase.verifyEqual(legacy(:, :, testCase.History, :), profiles(:, :, testCase.History, :));
        end

        function historyBasedResultsAreTheSameWithAndWithoutAnOutput(testCase)
            options = testCase.library('loose_reference');
            options.score_only = true;
            [scores, profiles] = benchmark({@TestOutputFallback.reachesZero, @TestOutputFallback.crashAfter}, options);
            [other_scores, other_profiles] = benchmark({@TestOutputFallback.reachesZero, @TestOutputFallback.reachesZero}, options);
            testCase.verifyEqual(scores, other_scores);
            testCase.verifyEqual(profiles(:, :, testCase.History, :), other_profiles(:, :, testCase.History, :));
        end

        function reportTextStatesTheRule(testCase)
            options = testCase.library('loose_reference');
            options.benchmark_id = 'saved';
            options.draw_hist_plots = 'none';
            benchmark({@TestOutputFallback.reachesZero, @TestOutputFallback.crashBefore}, options);
            files = dir(fullfile(testCase.WorkDir, 'saved', '**', 'report.txt'));
            testCase.assertEqual(numel(files), 1);
            testCase.verifySubstring(fileread(fullfile(files(1).folder, files(1).name)), ...
                'output fallback (never passes the output-based test)');
        end

        function singleProblemScoreReadsHistoriesOnly(testCase)
            reference = struct('fun', 4, 'maxcv', 0, 'kind', 'merit', 'source', 'f(x0)');
            options = struct('score_only', true, 'silent', true, 'solver_verbose', 0, 'n_jobs', 1, 'max_eval_factor', 2);
            options.problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', 2, 'reference', reference));
            scores = benchmark({@TestOutputFallback.reachesZero, @TestOutputFallback.crashAfter}, options);
            testCase.verifyEqual(scores(1), scores(2));
        end
    end
    methods
        function options = library(testCase, kind)
            root = fullfile(testCase.WorkDir, 'library');
            mkdir(root);
            TestOutputFallback.writeFile(fullfile(root, 'fallbacklib_select.m'), sprintf( ...
                'function names = fallbacklib_select(options)\n names = {''A'', ''B''};\n end\n'));
            if strcmp(kind, 'loose_reference')
                % The initial point meets the target of a reference at f(x0).
                body = ['s = struct(''fun'', @(x) sum(x.^2), ''x0'', x0, ''name'', name);\n' ...
                    's.reference = struct(''fun'', x0^2, ''maxcv'', 0, ''kind'', ''merit'', ''source'', ''f(x0)'');\n'];
            else
                % The convention for merit_init = +Inf passes its merit of +Inf.
                body = 's = struct(''fun'', @TestOutputFallback.infiniteAwayFromZero, ''x0'', x0, ''name'', name);\n';
            end
            TestOutputFallback.writeFile(fullfile(root, 'fallbacklib_load.m'), sprintf([ ...
                'function problem = fallbacklib_load(name)\n' ...
                'x0 = 2;\nif strcmp(name, ''B''), x0 = 3; end\n', body, ...
                'problem = Problem(s);\nend\n']));
            registerProblemLibrary(struct('name', 'fallbacklib', 'root', root, ...
                'select_function', 'fallbacklib_select', 'load_function', 'fallbacklib_load'));
            options = struct('plibs', {{'fallbacklib'}}, 'ptype', 'u', 'mindim', 1, 'maxdim', 1, ...
                'n_jobs', 1, 'max_eval_factor', 2, 'max_tol_order', 1, 'silent', true, 'solver_verbose', 0);
        end

        function options = loadOptions(~)
            options = struct('load', 'latest', 'benchmark_id', 'saved', 'score_only', true, ...
                'silent', true, 'max_eval_factor', 2, 'max_tol_order', 1);
        end
    end
end

function rmappdata_if_present()
    if isappdata(0, 'TestOutputFallbackSeen')
        rmappdata(0, 'TestOutputFallbackSeen');
    end
end
