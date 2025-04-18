classdef TestSolveOne < matlab.unittest.TestCase
    methods (Test)
        function testWithValidInput(testCase)
            % Test whether the function returns the correct outputs when given valid input.

            problem_name = 'ALLINITU';
            solvers = {@fminsearch_test1, @fminsearch_test2};
            feature = Feature('permuted', 'n_runs', 1);
            len_problem_names = length('ALLINITU');
            other_options.solver_names = {'fminsearch', 'fminunc'};
            other_options.solver_isrand = [false, false];
            other_options.custom_problem_loader = [];
            profile_options.project_x0 = true;
            profile_options.max_eval_factor = 500;
            profile_options.keep_pool = false;
            profile_options.silent = false;
            profile_options.seed = 1;
            profile_options.solver_verbose = 1;
            [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_name, problem_type, problem_dim, problem_con, computation_time, solvers_success] = solveOneProblem(problem_name, solvers, feature, len_problem_names, profile_options, other_options, false, '.');
            testCase.verifyNotEmpty(fun_histories);
            testCase.verifyNotEmpty(maxcv_histories);
            testCase.verifyNotEmpty(fun_out);
            testCase.verifyNotEmpty(maxcv_out);
            testCase.verifyNotEmpty(fun_init);
            testCase.verifyNotEmpty(maxcv_init);
            testCase.verifyNotEmpty(n_eval);
            testCase.verifyNotEmpty(problem_name);
            testCase.verifyNotEmpty(problem_type);
            testCase.verifyNotEmpty(problem_dim);
            testCase.verifyNotEmpty(problem_con);
            testCase.verifyNotEmpty(computation_time);
            testCase.verifyNotEmpty(solvers_success);
        end
    end
end

function x = fminsearch_test1(fun, x0)

x = fminsearch(fun, x0);
end

function x = fminsearch_test2(fun, x0)

options = optimset('MaxFunEvals', 1000);
x = fminunc(fun, x0, options);
end