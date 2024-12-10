classdef TestCheckValidityProfileOptions < matlab.unittest.TestCase
    methods (Test)

        function testErrors(testCase)

            options = struct();

            options.n_jobs = 1.1;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:n_jobsNotValid");
            options = rmfield(options, 'n_jobs');

            options.keep_pool = 2;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:keep_poolNotValid");
            options = rmfield(options, 'keep_pool');

            options.seed = 0;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:seedNotValid");
            options = rmfield(options, 'seed');

            options.benchmark_id = 1;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:benchmark_idNotValid");
            options = rmfield(options, 'benchmark_id');

            options.range_type = 'other';
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:range_typeNotValid");
            options = rmfield(options, 'range_type');

            options.savepath = 1;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:savepathNotValid");
            options = rmfield(options, 'savepath');

            options.max_tol_order = 20;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:max_tol_orderNotValid");
            options = rmfield(options, 'max_tol_order');

            options.max_eval_factor = 0;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:max_eval_factorNotValid");
            options = rmfield(options, 'max_eval_factor');

            options.project_x0 = 2;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:project_x0NotValid");
            options = rmfield(options, 'project_x0');

            options.run_plain = 2;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:run_plainNotValid");
            options = rmfield(options, 'run_plain');

            options.summarize_performance_profiles = 2;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:summarize_performance_profilesNotValid");
            options = rmfield(options, 'summarize_performance_profiles');

            options.summarize_data_profiles = 2;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:summarize_data_profilesNotValid");
            options = rmfield(options, 'summarize_data_profiles');

            options.summarize_log_ratio_profiles = 2;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:summarize_log_ratio_profilesNotValid");
            options = rmfield(options, 'summarize_log_ratio_profiles');

            options.summarize_output_based_profiles = 2;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:summarize_output_based_profilesNotValid");
            options = rmfield(options, 'summarize_output_based_profiles');

            options.silent = 2;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:silentNotValid");
            options = rmfield(options, 'silent');

            options.solver_verbose = 3;
            testCase.verifyError(@() checkValidityProfileOptions(options), "MATLAB:checkValidityProfileOptions:solver_verboseNotValid");
        end

    end

end