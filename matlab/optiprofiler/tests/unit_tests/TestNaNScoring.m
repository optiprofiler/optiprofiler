classdef TestNaNScoring < matlab.unittest.TestCase
    % Public benchmark and MAT reload must not turn undefined values into passes.
    properties
        Options
        Directory
    end
    methods (TestMethodSetup)
        function setUpLibrary(testCase)
            testCase.Directory = tempname;
            mkdir(testCase.Directory);
            old_cwd = pwd;
            old_path = path;
            old_registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            testCase.addTeardown(@() rmdir(testCase.Directory, 's'));
            testCase.addTeardown(@() cd(old_cwd));
            testCase.addTeardown(@() path(old_path));
            testCase.addTeardown(@() setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', old_registry));
            cd(testCase.Directory);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(testCase.Directory, 'registry.mat'));
            fixtures = fullfile(fileparts(mfilename('fullpath')), '../fixtures/nanvalidity');
            if exist('registerProblemLibrary', 'file') == 2
                registerProblemLibrary(struct('name', 'nanvalidity', 'root', fixtures, ...
                    'select_function', 'nanvalidity_select', 'load_function', 'nanvalidity_load'));
            else
                % The paper maintenance line uses the older folder catalog.
                % Add only a temporary empty catalog marker; fixture functions
                % remain on the test path, and teardown removes the marker.
                marker = fullfile(fileparts(which('benchmark')), '../problem_libs/nanvalidity');
                assert(exist(marker, 'dir') ~= 7, 'The test catalog marker already exists.');
                mkdir(marker);
                testCase.addTeardown(@() rmdir(marker));
                addpath(fixtures);
            end
            testCase.Options = struct('plibs', {{'nanvalidity'}}, 'ptype', 'n', ...
                'mindim', 1, 'maxdim', 1, 'feature_name', 'plain', 'n_runs', 2, ...
                'n_jobs', 1, 'max_eval_factor', 3, 'max_tol_order', 1, ...
                'draw_hist_plots', 'none', 'solver_names', {{'one', 'other'}}, ...
                'benchmark_id', 'nan-scoring', 'savepath', testCase.Directory, ...
                'score_only', true, 'silent', true);
        end
    end
    methods (Test)
        function testUndefinedEvaluationsCannotPass(testCase)
            for name = {'NAN_INIT_CV', 'NAN_INIT_FUN', 'NAN_LATER_CV', 'NAN_LATER_FUN', 'INF_INIT_NAN_LATER'}
                for custom_merit = [false, true]
                    options = testCase.Options;
                    options.problem_names = name;
                    if custom_merit
                        options.merit_fun = @(f, cv, cv0) 0;
                    end
                    [~, ~, curves] = benchmark({@TestNaNScoring.visitOne, @TestNaNScoring.visitOneAgain}, options);
                    testCase.verifyEqual(TestNaNScoring.fraction(curves, 'hist', 1), 0, name{1});
                    testCase.verifyEqual(TestNaNScoring.fraction(curves, 'out', 1), 0, name{1});
                end
            end
        end

        function testValidInitialConventions(testCase)
            for name = {'FINITE', 'VALID_INF'}
                options = testCase.Options;
                options.problem_names = name;
                [~, ~, curves] = benchmark({@TestNaNScoring.visitOne, @TestNaNScoring.visitOneAgain}, options);
                testCase.verifyEqual(TestNaNScoring.fraction(curves, 'hist', 1), 1);
                testCase.verifyEqual(TestNaNScoring.fraction(curves, 'out', 1), 1);
            end
        end

        function testValidPointsAfterNaN(testCase)
            options = testCase.Options;
            options.problem_names = {'INF_INIT_NAN_LATER'};
            [~, ~, curves] = benchmark({@TestNaNScoring.invalidThenValid, @TestNaNScoring.visitOne}, options);
            testCase.verifyEqual(TestNaNScoring.fraction(curves, 'hist', 1), 1);
            testCase.verifyEqual(TestNaNScoring.fraction(curves, 'out', 1), 1);
            testCase.verifyEqual(TestNaNScoring.fraction(curves, 'hist', 2), 0);
            testCase.verifyEqual(TestNaNScoring.fraction(curves, 'out', 2), 0);
            curve = curves{1}.hist.data{1, end};
            testCase.verifyEqual(curve(1, find(curve(2, :) > 0, 1)), 1);
        end

        function testSavedValiditySurvivesLoad(testCase)
            options = testCase.Options;
            options.problem_names = {'FINITE', 'NAN_INIT_CV', 'INF_INIT_NAN_LATER'};
            options.score_only = false;
            [~, ~, original] = benchmark({@TestNaNScoring.visitOne, @TestNaNScoring.visitOneAgain}, options);
            files = dir(fullfile(testCase.Directory, '**', 'data_for_loading.mat'));
            source = fullfile(files(1).folder, files(1).name);
            before = TestNaNScoring.bytes(source);
            options.load = 'latest';
            options.score_only = true;
            [~, ~, loaded] = benchmark({}, options);
            for channel = {'hist', 'out'}
                testCase.verifyEqual(TestNaNScoring.fraction(original, channel{1}, 1), 1/3, 'AbsTol', 1e-12);
                testCase.verifyEqual(TestNaNScoring.fraction(loaded, channel{1}, 1), 1/3, 'AbsTol', 1e-12);
            end
            testCase.verifyEqual(TestNaNScoring.bytes(source), before);
            report = fileread(fullfile(files(1).folder, 'report.txt'));
            testCase.verifyTrue(contains(report, 'undefined initial objective or constraint'));
        end
    end
    methods (Static)
        function x = visitOne(fun, x0, varargin)
            x = ones(size(x0));
            fun(x);
        end
        function x = visitOneAgain(fun, x0, varargin)
            x = TestNaNScoring.visitOne(fun, x0, varargin{:});
        end
        function x = invalidThenValid(fun, x0, varargin)
            fun(ones(size(x0)));
            fun(x0);
            x = x0;
        end
        function value = fraction(curves, channel, solver)
            curve = curves{1}.(channel).perf{solver, end};
            if isempty(curve)
                value = 0;
            else
                value = curve(2, end);
            end
        end
        function values = bytes(filename)
            fid = fopen(filename, 'rb');
            cleanup = onCleanup(@() fclose(fid));
            values = fread(fid, Inf, '*uint8');
        end
    end
end
