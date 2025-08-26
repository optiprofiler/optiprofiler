classdef TestCheckValidityProfileOptions < matlab.unittest.TestCase
    methods (Test)

        function testErrors(testCase)

            solvers = {@fminsearch, @fminunc};
            options = struct();


            options.n_jobs = 1.1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:n_jobsNotValid");
            options = rmfield(options, 'n_jobs');


            options.seed = -1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:seedNotValid");
            options = rmfield(options, 'seed');


            options.benchmark_id = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:benchmark_idNotValid");
            options = rmfield(options, 'benchmark_id');


            options.solver_names = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:solver_namesNotCellOfcharstr");
            options = rmfield(options, 'solver_names');
            options.solver_names = {'a', 'b', 'c'};
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:solver_namesAndsolversLengthNotSame");
            options = rmfield(options, 'solver_names');


            options.solver_isrand = 2;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:solver_israndNotLogical");
            options = rmfield(options, 'solver_isrand');


            options.feature_stamp = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:feature_stampNotcharstr");
            options = rmfield(options, 'feature_stamp');


            options.feature_stamp = 'abc$';
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:feature_stampNotValid");
            options = rmfield(options, 'feature_stamp');


            options.errorbar_type = 'other';
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:errorbar_typeNotValid");
            options = rmfield(options, 'errorbar_type');


            options.savepath = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:savepathNotValid");
            options = rmfield(options, 'savepath');


            options.max_tol_order = 20;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:max_tol_orderNotValid");
            options = rmfield(options, 'max_tol_order');


            options.max_eval_factor = 0;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:max_eval_factorNotValid");
            options = rmfield(options, 'max_eval_factor');


            options.merit_fun = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:merit_funNotValid");
            options = rmfield(options, 'merit_fun');


            options.project_x0 = 2;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:project_x0NotValid");
            options = rmfield(options, 'project_x0');


            options.run_plain = 2;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:run_plainNotValid");
            options = rmfield(options, 'run_plain');


            options.score_only = 2;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:score_onlyNotValid");
            options = rmfield(options, 'score_only');


            options.summarize_performance_profiles = 2;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:summarize_performance_profilesNotValid");
            options = rmfield(options, 'summarize_performance_profiles');


            options.summarize_data_profiles = 2;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:summarize_data_profilesNotValid");
            options = rmfield(options, 'summarize_data_profiles');


            options.summarize_log_ratio_profiles = 2;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:summarize_log_ratio_profilesNotValid");
            options = rmfield(options, 'summarize_log_ratio_profiles');


            options.summarize_output_based_profiles = 2;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:summarize_output_based_profilesNotValid");
            options = rmfield(options, 'summarize_output_based_profiles');


            options.silent = 2;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:silentNotValid");
            options = rmfield(options, 'silent');


            options.solver_verbose = 3;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:solver_verboseNotValid");
            options = rmfield(options, 'solver_verbose');


            options.semilogx = 2;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:semilogxNotValid");
            options = rmfield(options, 'semilogx');


            options.normalized_scores = 2;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:normalized_scoresNotValid");
            options = rmfield(options, 'normalized_scores');


            options.score_weight_fun = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:score_weight_funNotValid");
            options = rmfield(options, 'score_weight_fun');


            options.score_fun = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:score_funNotValid");
            options = rmfield(options, 'score_fun');


            options.load = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:loadNotValid");
            options = rmfield(options, 'load');


            options.solvers_to_load = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:solvers_to_loadNoLoad");
            options.load = 'latest';
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:solvers_to_loadNotValid");
            options.solvers_to_load = [1 2];
            options.solver_names = {'a' 'b' 'c'};
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:solver_namesAndsolversLengthNotSame");
            options = rmfield(options, 'solvers_to_load');
            options = rmfield(options, 'load');
            options = rmfield(options, 'solver_names');


            options.line_colors = {'w'};
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:line_colorsCellNotValid");
            options.line_colors = [0 0 0 0; 0 0 0 0];
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:line_colorsMatrixNotValid");
            options.line_colors = 'w';
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:line_colorsNotValid");
            options = rmfield(options, 'line_colors');


            options.line_styles = {'--.p'};
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:line_stylesNotValid");
            options = rmfield(options, 'line_styles');


            options.line_widths = 0;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:line_widthsNotValid");
            options = rmfield(options, 'line_widths');


            options.bar_colors = {'w'};
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:bar_colorsCellNotValid");
            options.bar_colors = [0 0 0 0; 0 0 0 0];
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:bar_colorsMatrixNotValid");
            options.bar_colors = 'w';
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:bar_colorsNotValid");
            options = rmfield(options, 'bar_colors');


            options.xlabel_performance_profile = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:xlabel_performance_profileNotValid");
            options = rmfield(options, 'xlabel_performance_profile');
            options.ylabel_performance_profile = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:ylabel_performance_profileNotValid");
            options = rmfield(options, 'ylabel_performance_profile');
            options.xlabel_data_profile = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:xlabel_data_profileNotValid");
            options = rmfield(options, 'xlabel_data_profile');
            options.ylabel_data_profile = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:ylabel_data_profileNotValid");
            options = rmfield(options, 'ylabel_data_profile');
            options.xlabel_log_ratio_profile = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:xlabel_log_ratio_profileNotValid");
            options = rmfield(options, 'xlabel_log_ratio_profile');
            options.ylabel_log_ratio_profile = 1;
            testCase.verifyError(@() checkValidityProfileOptions(solvers, options), "MATLAB:checkValidityProfileOptions:ylabel_log_ratio_profileNotValid");
            options = rmfield(options, 'ylabel_log_ratio_profile');
        end

    end

end