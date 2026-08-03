classdef TestMeritFunCompute < matlab.unittest.TestCase
%TESTMERITFUNCOMPUTE Regression tests for meritFunCompute input shapes.
%   The 'single' case receives per-problem data of shape
%   (n_solvers, n_runs, n_evals) with maxcv_inits of length n_runs. Loaded
%   histories slice a row out of the saved (n_problems, n_runs) matrix while
%   fresh runs pass a column; both orientations must produce the same result.

    properties (Access = private)
        CurrentDirectory
    end

    methods (TestMethodSetup)

        function goToPrivateFolder(testCase)
            testCase.CurrentDirectory = pwd;
            test_dir = fileparts(mfilename('fullpath'));
            source_dir = fullfile(test_dir, '../../src');
            addpath(source_dir);
            cd(fullfile(source_dir, 'private'));
        end

    end

    methods (TestMethodTeardown)

        function restoreFolder(testCase)
            cd(testCase.CurrentDirectory);
        end

    end

    methods (Test)

        function singleColumnInits(testCase)
            fun_values = zeros(2, 5, 6);
            maxcv_values = zeros(size(fun_values));
            inits = (1:5).';
            merit_fun = @(f, cv, cv_init) f + cv + cv_init;
            merits = meritFunCompute(merit_fun, fun_values, maxcv_values, inits, 'single');
            testCase.verifyEqual(size(merits), [2, 5, 6]);
            testCase.verifyEqual(squeeze(merits(1, :, 1)), 1:5);
        end

        function singleRowInitsMatchesColumn(testCase)
            fun_values = zeros(2, 5, 6);
            maxcv_values = zeros(size(fun_values));
            merit_fun = @(f, cv, cv_init) f + cv + cv_init;
            merits_column = meritFunCompute(merit_fun, fun_values, maxcv_values, (1:5).', 'single');
            merits_row = meritFunCompute(merit_fun, fun_values, maxcv_values, 1:5, 'single');
            testCase.verifyEqual(merits_row, merits_column);
        end

        function singleScalarInits(testCase)
            fun_values = zeros(2, 1, 6);
            maxcv_values = zeros(size(fun_values));
            merit_fun = @(f, cv, cv_init) f + cv + cv_init;
            merits = meritFunCompute(merit_fun, fun_values, maxcv_values, 3, 'single');
            testCase.verifyEqual(size(merits), [2, 1, 6]);
            testCase.verifyEqual(merits(2, 1, 4), 3);
        end

        function singleTwoDimensional(testCase)
            % Case 5: (n_solvers, n_runs) values with n_runs inits.
            fun_values = zeros(3, 4);
            maxcv_values = zeros(size(fun_values));
            merit_fun = @(f, cv, cv_init) f + cv + cv_init;
            merits_column = meritFunCompute(merit_fun, fun_values, maxcv_values, (1:4).', 'single');
            merits_row = meritFunCompute(merit_fun, fun_values, maxcv_values, 1:4, 'single');
            testCase.verifyEqual(size(merits_column), [3, 4]);
            testCase.verifyEqual(merits_row, merits_column);
        end

        function multipleCaseUnchanged(testCase)
            % Case 1: (n_problems, n_solvers, n_runs, n_evals) with
            % (n_problems, n_runs) inits must keep working unchanged.
            fun_values = zeros(2, 3, 4, 5);
            maxcv_values = zeros(size(fun_values));
            inits = repmat((1:4), 2, 1);
            merit_fun = @(f, cv, cv_init) f + cv + cv_init;
            merits = meritFunCompute(merit_fun, fun_values, maxcv_values, inits, 'multiple');
            testCase.verifyEqual(size(merits), [2, 3, 4, 5]);
            testCase.verifyEqual(squeeze(merits(1, 1, :, 1)).', 1:4);
        end

        function mismatchedInitsStillError(testCase)
            fun_values = zeros(2, 5, 6);
            maxcv_values = zeros(size(fun_values));
            merit_fun = @(f, cv, cv_init) f + cv + cv_init;
            testCase.verifyError(@() meritFunCompute(merit_fun, fun_values, maxcv_values, (1:4).', 'single'), ...
                'MATLAB:meritFunCompute:size_mismatch');
        end

    end

end
