classdef TestSolveAll < matlab.unittest.TestCase
    methods (Test)
        function testWithValidInput(testCase)
            % Test whether the function returns the correct outputs when given valid input.

            cutest_problem_names = {'ALLINITU', 'BARD'};
            custom_problem_loader = [];
            custom_problem_names = [];
            solvers = {@fminsearch_test, @fminunc_test};
            labels = {'fminsearch', 'fminunc'};
            feature = Feature('plain');
            profile_options.max_eval_factor = 500;
            profile_options.n_jobs = 1;
            profile_options.project_x0 = false;
            [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_names, problem_dimensions, computation_times] = solveAllProblems(cutest_problem_names, custom_problem_loader, custom_problem_names, solvers, labels, feature, profile_options);
            testCase.verifyNotEmpty(fun_histories);
            testCase.verifyNotEmpty(maxcv_histories);
            testCase.verifyNotEmpty(fun_out);
            testCase.verifyNotEmpty(maxcv_out);
            testCase.verifyNotEmpty(fun_init);
            testCase.verifyNotEmpty(maxcv_init);
            testCase.verifyNotEmpty(n_eval);
            testCase.verifyNotEmpty(problem_names);
            testCase.verifyNotEmpty(problem_dimensions);
            testCase.verifyNotEmpty(computation_times);

            profile_options.n_jobs = 2;
            [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_names, problem_dimensions, computation_times] = solveAllProblems(cutest_problem_names, custom_problem_loader, custom_problem_names, solvers, labels, feature, profile_options);
            testCase.verifyNotEmpty(fun_histories);
            testCase.verifyNotEmpty(maxcv_histories);
            testCase.verifyNotEmpty(fun_out);
            testCase.verifyNotEmpty(maxcv_out);
            testCase.verifyNotEmpty(fun_init);
            testCase.verifyNotEmpty(maxcv_init);
            testCase.verifyNotEmpty(n_eval);
            testCase.verifyNotEmpty(problem_names);
            testCase.verifyNotEmpty(problem_dimensions);
            testCase.verifyNotEmpty(computation_times);
        end

        function testAllProblemFailedToLoad(testCase)
            % Test whether the function returns empty outputs when all problems failed to load.

            cutest_problem_names = {'A', 'B'};
            custom_problem_loader = [];
            custom_problem_names = [];
            solvers = {@fminsearch_test, @fminunc_test};
            labels = {'fminsearch', 'fminunc'};
            feature = Feature('plain');
            profile_options.max_eval_factor = 500;
            profile_options.n_jobs = 1;
            profile_options.project_x0 = false;
            [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_names, problem_dimensions, computation_times] = solveAllProblems(cutest_problem_names, custom_problem_loader, custom_problem_names, solvers, labels, feature, profile_options);
            testCase.verifyEmpty(fun_histories);
            testCase.verifyEmpty(maxcv_histories);
            testCase.verifyEmpty(fun_out);
            testCase.verifyEmpty(maxcv_out);
            testCase.verifyEmpty(fun_init);
            testCase.verifyEmpty(maxcv_init);
            testCase.verifyEmpty(n_eval);
            testCase.verifyEmpty(problem_names);
            testCase.verifyEmpty(problem_dimensions);
            testCase.verifyEmpty(computation_times);
        end
    end
end

function x = fminsearch_test(fun, x0)

    x = fminsearch(fun, x0);

end

function x = fminunc_test(fun, x0)

    x = fminunc(fun, x0);

end