classdef TestSolveOne < matlab.unittest.TestCase
    methods (Test)
        function testWithValidInput(testCase)
            % Test whether the function returns the correct outputs when given valid input.

            problem_name = 'ALLINITU';
            solvers = {@fminsearch_test, @fminunc_test};
            labels = {'fminsearch', 'fminunc'};
            feature = Feature('permuted', 'n_runs', 1);
            custom_problem_loader = [];
            profile_options.project_x0 = true;
            profile_options.max_eval_factor = 500;
            [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_name, problem_n, computation_time] = solveOneProblem(problem_name, solvers, labels, feature, custom_problem_loader, profile_options);
            testCase.verifyNotEmpty(fun_histories);
            testCase.verifyNotEmpty(maxcv_histories);
            testCase.verifyNotEmpty(fun_out);
            testCase.verifyNotEmpty(maxcv_out);
            testCase.verifyNotEmpty(fun_init);
            testCase.verifyNotEmpty(maxcv_init);
            testCase.verifyNotEmpty(n_eval);
            testCase.verifyNotEmpty(problem_name);
            testCase.verifyNotEmpty(problem_n);
            testCase.verifyNotEmpty(computation_time);

            problem_name = {'EXTRA', 'A'};
            custom_problem_loader = @custom_loader;
            feature = Feature('affine_transformed', 'n_runs', 1);
            [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_name, problem_n, computation_time] = solveOneProblem(problem_name, solvers, labels, feature, custom_problem_loader, profile_options);
            testCase.verifyNotEmpty(fun_histories);
            testCase.verifyNotEmpty(maxcv_histories);
            testCase.verifyNotEmpty(fun_out);
            testCase.verifyNotEmpty(maxcv_out);
            testCase.verifyNotEmpty(fun_init);
            testCase.verifyNotEmpty(maxcv_init);
            testCase.verifyNotEmpty(n_eval);
            testCase.verifyNotEmpty(problem_name);
            testCase.verifyNotEmpty(problem_n);
            testCase.verifyNotEmpty(computation_time);
        end

        function testSolverSignature(testCase)
            % Test whether the function returns correct outputs or throws an error when given a solver with corresponding signatures.

            problem_name = 'ALLINITU';
            solvers = {@fminsearch_test2, @fminunc_test2};
            labels = {'fminsearch', 'fminunc'};
            feature = Feature('plain');
            custom_problem_loader = [];
            profile_options.project_x0 = false;
            profile_options.max_eval_factor = 500;
            [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_name, problem_n, computation_time] = solveOneProblem(problem_name, solvers, labels, feature, custom_problem_loader, profile_options);

            solvers = {@fminsearch_test3, @fminunc_test3};
            solveOneProblem(problem_name, solvers, labels, feature, custom_problem_loader, profile_options);

            solvers = {@fminsearch_test4, @fminunc_test4};
            solveOneProblem(problem_name, solvers, labels, feature, custom_problem_loader, profile_options);

            solvers = {@fminsearch_test5, @fminunc_test5};
            testCase.verifyError(@() solveOneProblem(problem_name, solvers, labels, feature, custom_problem_loader, profile_options), "MATLAB:solveOneProblem:InvalidSolverSignature");
        end


    end
end

function x = fminsearch_test(fun, x0)

    x = fminsearch(fun, x0);

end

function x = fminunc_test(fun, x0)

    x = fminunc(fun, x0);

end

function x = fminsearch_test2(fun, x0, xl, xu)

    x = fminsearch(fun, x0);

end

function x = fminunc_test2(fun, x0, xl, xu)

    x = fminunc(fun, x0);

end

function x = fminsearch_test3(fun, x0, xl, xu, aub, bub, aeq, beq)

    x = fminsearch(fun, x0);

end

function x = fminunc_test3(fun, x0, xl, xu, aub, bub, aeq, beq)

    x = fminunc(fun, x0);
    
end

function x = fminsearch_test4(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq)

    x = fminsearch(fun, x0);

end

function x = fminunc_test4(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq)

    x = fminunc(fun, x0);

end

function x = fminsearch_test5(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, options)

    x = fminsearch(fun, x0);

end

function x = fminunc_test5(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, options)

    x = fminunc(fun, x0);

end

function p = custom_loader(problem_name)

    p = Problem(struct('fun', @(x) x' * x, 'x0', [2; 2]));

end