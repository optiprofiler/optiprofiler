classdef TestSolveAll < matlab.unittest.TestCase
    methods (Test)
        function testWithValidInput(testCase)
            % Test whether the function returns the correct outputs when given valid input.

            solvers = {@fminsearch_test1, @fminsearch_test2};
            plib = 's2mpj';
            feature = Feature('plain');

            problem_options.ptype = 'u';
            problem_options.maxdim = 11;
            problem_options.mindim = 11;

            profile_options.n_jobs = 1;
            profile_options = getDefaultProfileOptions(solvers, feature, profile_options);
            
            results = solveAllProblems(solvers, plib, feature, problem_options, profile_options, false, '');
            testCase.verifyNotEmpty(results);
            testCase.verifyNotEmpty(results.plib);
            testCase.verifyNotEmpty(results.solver_names);
            testCase.verifyNotEmpty(results.feature_stamp);
            testCase.verifyNotEmpty(results.fun_histories);
            testCase.verifyNotEmpty(results.maxcv_histories);
            testCase.verifyNotEmpty(results.fun_outs);
            testCase.verifyNotEmpty(results.maxcv_outs);
            testCase.verifyNotEmpty(results.fun_inits);
            testCase.verifyNotEmpty(results.maxcv_inits);
            testCase.verifyNotEmpty(results.n_evals);
            testCase.verifyNotEmpty(results.problem_names);
            testCase.verifyNotEmpty(results.problem_types);
            testCase.verifyNotEmpty(results.problem_dims);
            testCase.verifyNotEmpty(results.problem_mbs);
            testCase.verifyNotEmpty(results.problem_mlcons);
            testCase.verifyNotEmpty(results.problem_mnlcons);
            testCase.verifyNotEmpty(results.problem_mcons);
            testCase.verifyNotEmpty(results.computation_times);
            testCase.verifyNotEmpty(results.solvers_successes);
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