classdef TestBenchmark < matlab.unittest.TestCase
    methods (Test)
        
        function testWithValidInput(testCase)
            % Test whether the function works well with valid input.

            solvers = {@fminsearch_test, @fminunc_test};
            benchmark(solvers);

            benchmark(solvers, {'plain', 'noisy'});
            
            options.feature_names = 'noisy';
            options.n_jobs = 0;
            options.run_plain = false;
            options.n_runs = 5;
            options.max_tol_order = 3;
            options.max_eval_factor = 10;
            options.problem_type = 'u';
            options.maxdim = 3;
            options.benchmark_id = 'unit-test';
            options.summarize_log_ratio_profiles = true;
            options.labels = {'fminsearch', 'fminunc'};
            benchmark(solvers, options);

            options.feature_names = 'plain';
            options.labels = {};
            options.custom_problem_names = {"A"};
            options.custom_problem_loader = @custom_loader;
            options.excludelist = {"SISSER2"};
            options.n_jobs = 1000;
            options.summarize_data_profiles = false;
            options.summarize_performance_profiles = false;
            options.cutest_problem_names = {'S308'};
            options.benchmark_id = 'test-special-case';
            savepath = fullfile(pwd, 'out', 'test-special-case');
            mkdir(savepath);
            % Create a txt file named 'a.txt' in the savepath.
            fid = fopen(fullfile(savepath, 'a.txt'), 'w');
            fclose(fid);
            mkdir(fullfile(savepath, 'time-unknown'));
            benchmark(solvers, options);
            mkdir(fullfile(savepath, 'time-unknown-abc'));
            benchmark(solvers, options);
            mkdir(fullfile(savepath, 'time-unknown-100'));
            benchmark(solvers, options);
        end

        function testErrors(testCase)
            % Test whether the function throws errors as expected.

            testCase.verifyError(@() benchmark(), "MATLAB:benchmark:solverMustBeProvided")

            solvers = {@fminsearch_test, @fminunc_test};
            options.custom_problem_names = {'A', 'B'};
            options.labels = {'fminsearch', 'fminunc'};
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:LoaderAndNamesNotSameTime")

            testCase.verifyError(@() benchmark(solvers, {1}), "MATLAB:benchmark:SecondArgumentWrongType")

            testCase.verifyError(@() benchmark(solvers, 'plain', options, 1), "MATLAB:benchmark:TooMuchInput")

            options = rmfield(options, 'custom_problem_names');
            options.feature_names = 'plain';
            testCase.verifyError(@() benchmark({1, 2}, options), "MATLAB:benchmark:solversWrongType")

            testCase.verifyError(@() benchmark({@fminsearch_test}), "MATLAB:benchmark:solversAtLeastTwo")

            options.labels = {1};
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:labelsNotCellOfcharstr")

            options.labels = {'fminsearch'};
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:labelsAndsolversLengthNotSame")

            options.custom_problem_names = {'A', 'B'};
            options.custom_problem_loader = 1;
            options = rmfield(options, 'labels');
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:customloaderNotFunctionHandle")

            options.custom_problem_names = {};
            options.custom_problem_loader = @(x) [1; 1] * [1; 1];
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:customnamesCanNotBeEmptyWhenHavingcustomloader")

            options.custom_problem_names = {'A', 'B'};
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:customloaderNotAcceptcustomnames")

            options.custom_problem_loader = {};
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:customloaderCanNotBeEmptyWhenHavingcustomnames")

            options.custom_problem_names = {1};
            options.custom_problem_loader = @custom_loader;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:customnamesNotcharstrOrCellOfcharstr")

            options = rmfield(options, 'custom_problem_names');
            options = rmfield(options, 'custom_problem_loader');
            options.a = {};
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:UnknownOptions")
            options = rmfield(options, 'a');

            options.problem_type = 1;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:problem_typeNotcharstr")

            options.problem_type = 'd';
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:problem_typeNotubln")
            options.problem_type = 'u';

            options.mindim = 'a';
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:mindimNotValid")
            options.mindim = 2;

            options.maxdim = 'a';
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:maxdimNotValid")

            options.maxdim = 1;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:maxdimSmallerThanmindim")
            options.maxdim = 2;

            options.mincon = 'a';
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:minconNotValid")
            options.mincon = 2;

            options.maxcon = 'a';
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:maxconNotValid")

            options.maxcon = 1;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:maxconSmallerThanmincon")
            options.maxcon = 2;

            options.excludelist = 1;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:excludelistNotCellOfcharstr")
            options = rmfield(options, 'excludelist');

            options.n_jobs = 1.1;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:n_jobsNotValid")
            options.n_jobs = 1;

            options.benchmark_id = 1;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:benchmark_idNotValid")
            options.benchmark_id = 'unit-test';

            options.range_type = 1;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:range_typeNotValid")
            options.range_type = 'meanstd';

            options.std_factor = 'a';
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:std_factorNotValid")
            options.std_factor = 1;

            options.savepath = 1;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:savepathNotValid")
            options = rmfield(options, 'savepath');

            options.max_tol_order = 'a';
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:max_tol_orderNotValid")
            options.max_tol_order = 1;

            options.max_eval_factor = 'a';
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:max_eval_factorNotValid")
            options.max_eval_factor = 1;

            options.project_x0 = 2;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:project_x0NotValid")
            options.project_x0 = true;

            options.run_plain = 2;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:run_plainNotValid")
            options.run_plain = true;

            options.summarize_performance_profiles = 2;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:summarize_performance_profilesNotValid")
            options.summarize_performance_profiles = true;

            options.summarize_data_profiles = 2;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:summarize_data_profilesNotValid")
            options.summarize_data_profiles = true;

            options.summarize_log_ratio_profiles = 2;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:summarize_log_ratio_profilesNotValid")
            options.summarize_log_ratio_profiles = false;

            options.cutest_problem_names = 1;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:cutest_problem_namesNotValid")

            options.cutest_problem_names = {};
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:AtLeastOneProblem")
            options.cutest_problem_names = {'A', 'B'};

            options.feature_names = 'all';
            options.n_runs = 2;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:OnlyOneFeatureWhenHavingfeature_options")
            options.feature_names = 'plain';
        end
    end
end

function x = fminsearch_test(fun, x0)

    x = fminsearch(fun, x0);

end

function x = fminunc_test(fun, x0)

    x = fminunc(fun, x0);

end

function p = custom_loader(problem_name)

    p = Problem(struct('fun', @(x) x' * x, 'x0', [1; 1]));

end