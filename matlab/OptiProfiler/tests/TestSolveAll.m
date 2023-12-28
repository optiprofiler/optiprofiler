classdef TestSolveAll < matlab.unittest.TestCase
    methods (Test)
        function testWithValidInput(testCase)
            % Define valid problem names, options, solvers, labels, feature, max_eval_factor, and profile_options
            problem_names = {'AKIVA', 'BEALE'};
            problem_options = struct('maxdim', 2);
            solvers = {@fminsearch_test};
            labels = {'fminsearch'};
            feature = Feature('plain');
            max_eval_factor = 500;
            profile_options = struct('n_jobs', 1);

            % Call the function under test
            [fun_values, maxcv_values, fun_inits, maxcv_inits, n_evals, problem_names_output, problem_dimensions] = solveAll(problem_names, problem_options, solvers, labels, feature, max_eval_factor, profile_options);

            % Verify the function outputs
            testCase.verifyNotEmpty(fun_values);
            testCase.verifyNotEmpty(maxcv_values);
            testCase.verifyNotEmpty(fun_inits);
            testCase.verifyNotEmpty(maxcv_inits);
            testCase.verifyNotEmpty(n_evals);
            testCase.verifyEqual(problem_names, problem_names_output);
            testCase.verifyGreaterThan(problem_dimensions, 0);
        end

        function testParallelComputing(testCase)
            problem_names = {'AKIVA', 'BEALE'};
            problem_options = struct('maxdim', 2);
            solvers = {@fminsearch_test, @fminunc_test};
            labels = {'fminsearch', 'fminunc'};
            feature = Feature('plain');
            max_eval_factor = 500;
            profile_options = struct('n_jobs', 2);  % Parallel computing option set to 2 jobs

            % Call the function under test
            [fun_values, maxcv_values, fun_inits, maxcv_inits, n_evals, problem_names_output, problem_dimensions] = solveAll(problem_names, problem_options, solvers, labels, feature, max_eval_factor, profile_options);

            % Verify the function outputs
            testCase.verifyNotEmpty(fun_values);
            testCase.verifyNotEmpty(maxcv_values);
            testCase.verifyNotEmpty(fun_inits);
            testCase.verifyNotEmpty(maxcv_inits);
            testCase.verifyNotEmpty(n_evals);
            testCase.verifyEqual(problem_names, problem_names_output);
            testCase.verifyGreaterThan(problem_dimensions, 0);
        end

        function testFevalFewerThanMaxEval(testCase)
            problem_names = {'AKIVA', 'BEALE'};
            problem_options = struct('maxdim', 2);
            solvers = {@fewerfevalsolver};
            labels = {'fewerfevalsolver'};
            feature = Feature('plain');
            max_eval_factor = 500;  % Set max_eval_factor to 0.1
            profile_options = struct('n_jobs', 1);
            
            % Define a solver function that will be called less than max_eval_factor times
            function [x, fval] = fewerfevalsolver(fun, x0, varargin)
                x = x0;
                fval = fun(x);
            end

            % Call the function under test
            [fun_values, maxcv_values, fun_inits, maxcv_inits, n_evals, problem_names_output, problem_dimensions] = solveAll(problem_names, problem_options, solvers, labels, feature, max_eval_factor, profile_options);

            % Verify the function outputs are in the correct size
            testCase.verifyEqual(size(fun_values), [2 1 1 1000]);
            testCase.verifyEqual(size(maxcv_values), [2 1 1 1000]);
        end
    end
end