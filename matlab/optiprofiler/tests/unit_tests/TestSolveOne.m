classdef TestSolveOne < matlab.unittest.TestCase
    methods (Test)
        function testWithValidInput(testCase)
            % Test whether the function returns the correct outputs when given valid input.

            problem_name = 'ALLINITU';
            problem = s2mpj_load(problem_name);
            solvers = {@fminsearch_test1, @fminsearch_test2};
            feature = Feature('permuted', 'n_runs', 1);
            len_problem_names = length('ALLINITU');
            profile_options.solver_names = {'fminsearch', 'fminunc'};
            profile_options.solver_isrand = [false, false];
            profile_options.project_x0 = true;
            profile_options.max_eval_factor = 500;
            profile_options.keep_pool = false;
            profile_options.silent = false;
            profile_options.seed = 1;
            profile_options.solver_verbose = 1;
            profile_options.score_only = false;
            result = solveOneProblem(solvers, problem, feature, problem_name, len_problem_names, profile_options, true, '');
            testCase.verifyNotEmpty(result);
            testCase.verifyNotEmpty(result.fun_history);
            testCase.verifyNotEmpty(result.maxcv_history);
            testCase.verifyNotEmpty(result.fun_out);
            testCase.verifyNotEmpty(result.maxcv_out);
            testCase.verifyNotEmpty(result.fun_init);
            testCase.verifyNotEmpty(result.maxcv_init);
            testCase.verifyNotEmpty(result.n_eval);
            testCase.verifyNotEmpty(result.problem_name);
            testCase.verifyNotEmpty(result.problem_type);
            testCase.verifyNotEmpty(result.problem_dim);
            testCase.verifyNotEmpty(result.problem_mb);
            testCase.verifyNotEmpty(result.problem_con);
            testCase.verifyNotEmpty(result.computation_time);
            testCase.verifyNotEmpty(result.solvers_success);
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