classdef TestBenchmark < matlab.unittest.TestCase
    methods (Test)
        
        function testWithValidInput(testCase)

            solvers = {@fminsearch_test1, @fminsearch_test2};
            options.benchmark_id = 'unit-test';

            % Test the case where problem set only contains one problem.
            options.problem = s_load('ROSENBR');
            benchmark(solvers, options);

            % Test plain feature.
            options.maxdim = 2;
            options.max_tol_order = 2;
            options = rmfield(options, 'problem');
            options.feature_name = 'plain';
            benchmark(solvers, options);

            % Test perturbed_x0 feature.
            options.feature_name = 'perturbed_x0';
            options.n_runs = 2;
            options.run_plain = true;
            options.feature_name = 'perturbed_x0';
            benchmark(solvers, options);

            % Test noisy feature.
            options.feature_name = 'noisy';
            options.noise_level = 1e-4;
            options.run_plain = false;
            benchmark(solvers, options);
            options = rmfield(options, 'noise_level');

            % Test truncated feature.
            options.feature_name = 'truncated';
            options.significant_digits = 4;
            options.perturbed_trailing_zeros = false;
            benchmark(solvers, options);
            options = rmfield(options, 'significant_digits');
            options = rmfield(options, 'perturbed_trailing_zeros');

            % Test permuted feature.
            options.feature_name = 'permuted';
            benchmark(solvers, options);

            % Test linearly_transformed feature.
            options.feature_name = 'linearly_transformed';
            options.rotated = true;
            options.condition_factor = 1;
            benchmark(solvers, options);
            options = rmfield(options, 'rotated');
            options = rmfield(options, 'condition_factor');

            % Test random_nan feature.
            options.feature_name = 'random_nan';
            options.rate_nan = 0.1;
            benchmark(solvers, options);

            rmdir('out', 's');

        end

        function testErrors(testCase)
            % Test whether the function throws errors as expected.

            testCase.verifyError(@() benchmark(), "MATLAB:benchmark:solverMustBeProvided")

            solvers = {@fminsearch_test1, @fminsearch_test2};
            options.problem = 'ROSENBR';
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:problemNotProblem")
            options = rmfield(options, 'problem');

            options.custom_problem_loader = 1;
            options.custom_problem_names = {};
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:customnamesEmpty")
            options = rmfield(options, 'custom_problem_loader');

            options.custom_problem_names = 'example';
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:LoaderAndNamesNotSameTime")
            options = rmfield(options, 'custom_problem_names');

            testCase.verifyError(@() benchmark(solvers, {1}), "MATLAB:benchmark:SecondArgumentWrongType")

            testCase.verifyError(@() benchmark(solvers, 'plain', options, 1), "MATLAB:benchmark:TooMuchInput")

            testCase.verifyError(@() benchmark({1, 2}, options), "MATLAB:benchmark:solversWrongType")

            testCase.verifyError(@() benchmark({@fminsearch_test}), "MATLAB:benchmark:solversAtLeastTwo")

            options.feature_name = 1;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:feature_nameNotcharstr")
            options = rmfield(options, 'feature_name');

            options.feature_name = 'a';
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:feature_nameNotValid")
            options = rmfield(options, 'feature_name');

            options.labels = {1, 2};
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:labelsNotCellOfcharstr")

            options.labels = {'a', 'b', 'c'};
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:labelsAndsolversLengthNotSame")
            options = rmfield(options, 'labels');

            options.custom_problem_names = {'A', 'B'};
            options.custom_problem_loader = 1;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:customloaderNotFunctionHandle")

            options.custom_problem_loader = @(x) x;
            options.custom_problem_names = 1;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:customnamesNotcharstrOrCellOfcharstr")

            options.custom_problem_names = {'A', 'B'};
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:customloaderNotAcceptcustomnames")
            options = rmfield(options, 'custom_problem_names');
            options = rmfield(options, 'custom_problem_loader');

            options.a = {'a'};
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:UnknownOptions")

        end
    end
end

function x = fminsearch_test1(fun, x0)

    x = fminsearch(fun, x0);
end

function x = fminsearch_test2(fun, x0)

    options = optimset('MaxFunEvals', 200);
    x = fminunc(fun, x0, options);
end