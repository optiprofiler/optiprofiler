classdef TestSolveAll < matlab.unittest.TestCase
    methods (Test)
        function testWithValidInput(testCase)
            % Test whether the function returns the correct outputs when given valid input.

            solvers = {@fminsearch_test1, @fminsearch_test2};
            other_options.cutest_problem_names = {'ALLINITU', 'BARD'};
            other_options.custom_problem_loader = [];
            other_options.custom_problem_names = [];
            other_options.solver_names = {'fminsearch1', 'fminsearch2'};
            other_options.solver_isrand = [false, false];
            feature = Feature('plain');
            profile_options.max_eval_factor = 500;
            profile_options.n_jobs = 1;
            profile_options.project_x0 = false;
            profile_options.keep_pool = false;
            profile_options.silent = false;
            profile_options.seed = 1;
            profile_options.solver_verbose = 1;
            [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_evals, problem_names, problem_types, problem_dims, problem_cons, computation_times, solvers_successes] = solveAllProblems(solvers, feature, profile_options, other_options, false, '.');
            testCase.verifyNotEmpty(fun_histories);
            testCase.verifyNotEmpty(maxcv_histories);
            testCase.verifyNotEmpty(fun_out);
            testCase.verifyNotEmpty(maxcv_out);
            testCase.verifyNotEmpty(fun_init);
            testCase.verifyNotEmpty(maxcv_init);
            testCase.verifyNotEmpty(n_evals);
            testCase.verifyNotEmpty(problem_names);
            testCase.verifyNotEmpty(problem_types);
            testCase.verifyNotEmpty(problem_dims);
            testCase.verifyNotEmpty(problem_cons);
            testCase.verifyNotEmpty(computation_times);
            testCase.verifyNotEmpty(solvers_successes);
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