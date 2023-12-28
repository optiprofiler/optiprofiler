classdef TestSolveOne < matlab.unittest.TestCase
    methods (Test)
        function testWithValidInput(testCase)
            % Define problem name, options, solvers, labels, feature, max_eval_factor
            problem_name = 'AKIVA';
            problem_options = struct('maxdim', 2);
            solvers = {@fminsearch_test};
            labels = {'fminsearch'};
            feature = Feature('plain');
            max_eval_factor = 500;

            % Call the function under test
            [fun_values, maxcv_values, fun_init, maxcv_init, n_eval, problem_name_output, problem_n] = solveOne(problem_name, problem_options, solvers, labels, feature, max_eval_factor);
            
            % Verify the function outputs
            testCase.verifyNotEmpty(fun_values);
            testCase.verifyNotEmpty(maxcv_values);
            testCase.verifyNotEmpty(fun_init);
            testCase.verifyNotEmpty(maxcv_init);
            testCase.verifyNotEmpty(n_eval);
            testCase.verifyEqual(problem_name, problem_name_output);
            testCase.verifyGreaterThan(problem_n, 0);

        end

        function testWithInvalidProblemName(testCase)

            problem_name = 'InvalidProblemName';
            problem_options = struct('maxdim', 2);
            solvers = {@fminsearch_test};
            labels = {'fminsearch'};
            feature = Feature('plain');
            max_eval_factor = 500;

            % Call the function under test and verify it behaves as expected when problem load fails
            [fun_values, maxcv_values, fun_init, maxcv_init, n_eval, problem_name_output, problem_n] = solveOne(problem_name, problem_options, solvers, labels, feature, max_eval_factor);

            % Verify the function output
            testCase.verifyEmpty(fun_values);
            testCase.verifyEmpty(maxcv_values);
            testCase.verifyEmpty(problem_n);
        end

        function testWithErrorSolver(testCase)
            
            problem_name = 'AKIVA';
            problem_options = struct('maxdim', 2);
            solvers = {@errorsolver};
            labels = {'errorsolver'};
            feature = Feature('plain');
            max_eval_factor = 500;

            % This function handle represents a solver that throws an error
            function [x, fval] = errorsolver(fun, x0, varargin)
                error("This solver only throws error!")
            end

            [fun_values, maxcv_values, fun_init, maxcv_init, n_eval, problem_name_output, problem_n] = solveOne(problem_name, problem_options, solvers, labels, feature, max_eval_factor);

        end

    end
end