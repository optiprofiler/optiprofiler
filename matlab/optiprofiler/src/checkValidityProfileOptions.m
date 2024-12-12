function profile_options = checkValidityProfileOptions(profile_options, solvers)
%CHECKVALIDITYPROFILEOPTIONS Check the validity of the options in profile_options

    if exist('parcluster', 'file') == 2
        myCluster = parcluster('local');
        nb_cores = myCluster.NumWorkers;
    else
        nb_cores = 1;
    end

    % Judge whether profile_options.n_jobs is a integer between 1 and nb_cores.
    if isfield(profile_options, ProfileOptionKey.N_JOBS.value)
        if ~isintegerscalar(profile_options.(ProfileOptionKey.N_JOBS.value))
            error("MATLAB:checkValidityProfileOptions:n_jobsNotValid", "The field 'n_jobs' of options should be a integer.");
        elseif profile_options.(ProfileOptionKey.N_JOBS.value) < 1
            profile_options.(ProfileOptionKey.N_JOBS.value) = 1;
        elseif profile_options.(ProfileOptionKey.N_JOBS.value) > nb_cores
            profile_options.(ProfileOptionKey.N_JOBS.value) = nb_cores;
        else
            profile_options.(ProfileOptionKey.N_JOBS.value) = round(profile_options.(ProfileOptionKey.N_JOBS.value));
        end
    end
    % Judge whether profile_options.keep_pool is a boolean.
    if isfield(profile_options, ProfileOptionKey.KEEP_POOL.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.KEEP_POOL.value))
            error("MATLAB:checkValidityProfileOptions:keep_poolNotValid", "The field 'keep_pool' of options should be a boolean.");
        end
    end
    % Judge whether profile_options.seed is a positive integer.
    if isfield(profile_options, ProfileOptionKey.SEED.value)
        if ~isintegerscalar(profile_options.(ProfileOptionKey.SEED.value)) || profile_options.(ProfileOptionKey.SEED.value) <= 0
            error("MATLAB:checkValidityProfileOptions:seedNotValid", "The field 'seed' of options should be a positive integer.");
        end
    end
    % Judge whether profile_options.benchmark_id is a char or a string and satisfies the file name requirements.
    if isfield(profile_options, ProfileOptionKey.BENCHMARK_ID.value)
        is_valid_foldername = @(x) ischarstr(x) && ~isempty(x) && all(ismember(char(x), ['a':'z', 'A':'Z', '0':'9', '_', '-', '.']));
        if ~ischarstr(profile_options.(ProfileOptionKey.BENCHMARK_ID.value)) || ~is_valid_foldername(profile_options.(ProfileOptionKey.BENCHMARK_ID.value))
            error("MATLAB:checkValidityProfileOptions:benchmark_idNotValid", "The field 'benchmark_id' of options should be a char or a string satisfying the strict file name requirements (only containing letters, numbers, underscores, hyphens, and dots).");
        end
    end
    % Judge whether profile_options.range_type is among 'minmax' and 'meanstd'.
    if isfield(profile_options, ProfileOptionKey.RANGE_TYPE.value)
        if ~ischarstr(profile_options.(ProfileOptionKey.RANGE_TYPE.value)) || ~ismember(char(profile_options.(ProfileOptionKey.RANGE_TYPE.value)), {'minmax', 'meanstd'})
            error("MATLAB:checkValidityProfileOptions:range_typeNotValid", "The field 'range_type' of options should be either 'minmax' or 'meanstd'.");
        end
    end
    % Judge whether profile_options.savepath is a string and exists. If not exists, create it.
    if isfield(profile_options, ProfileOptionKey.SAVEPATH.value)
        if ~ischarstr(profile_options.(ProfileOptionKey.SAVEPATH.value))
            error("MATLAB:checkValidityProfileOptions:savepathNotValid", "The field 'savepath' of options should be a char or a string.");
        elseif ~exist(profile_options.(ProfileOptionKey.SAVEPATH.value), 'dir')
            status = mkdir(profile_options.(ProfileOptionKey.SAVEPATH.value));
            if ~status
                error("MATLAB:checkValidityProfileOptions:savepathNotExist", "The field 'savepath' of options does not exist and cannot be created.");
            end
        end
    end
    % Judge whether profile_options.max_tol_order is a positive integer smaller than or equal to 15.
    if isfield(profile_options, ProfileOptionKey.MAX_TOL_ORDER.value)
        if ~isintegerscalar(profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value)) || profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) <= 0 || profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) > 15
            error("MATLAB:checkValidityProfileOptions:max_tol_orderNotValid", "The field 'max_tol_order' of options should be a positive integer smaller than or equal to 15.");
        end
    end
    % Judge whether profile_options.max_eval_factor is a positive integer.
    if isfield(profile_options, ProfileOptionKey.MAX_EVAL_FACTOR.value)
        if ~isintegerscalar(profile_options.(ProfileOptionKey.MAX_EVAL_FACTOR.value)) || profile_options.(ProfileOptionKey.MAX_EVAL_FACTOR.value) <= 0
            error("MATLAB:checkValidityProfileOptions:max_eval_factorNotValid", "The field 'max_eval_factor' of options should be a positive integer.");
        end
    end
    % Judge whether profile_options.project_x0 is a boolean.
    if isfield(profile_options, ProfileOptionKey.PROJECT_X0.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.PROJECT_X0.value))
            error("MATLAB:checkValidityProfileOptions:project_x0NotValid", "The filed 'project_x0' of options should be a boolean.");
        end
    end
    % Judge whether profile_options.run_plain is a boolean.
    if isfield(profile_options, ProfileOptionKey.RUN_PLAIN.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.RUN_PLAIN.value))
            error("MATLAB:checkValidityProfileOptions:run_plainNotValid", "The field 'run_plain' of options should be a boolean.");
        end
    end
    % Judge whether profile_options.summarize_performance_profiles is a boolean.
    if isfield(profile_options, ProfileOptionKey.SUMMARIZE_PERFORMANCE_PROFILES.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.SUMMARIZE_PERFORMANCE_PROFILES.value))
            error("MATLAB:checkValidityProfileOptions:summarize_performance_profilesNotValid", "The field 'summarize_performance_profiles' of options should be a boolean.");
        end
    end
    % Judge whether profile_options.summarize_data_profiles is a boolean.
    if isfield(profile_options, ProfileOptionKey.SUMMARIZE_DATA_PROFILES.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.SUMMARIZE_DATA_PROFILES.value))
            error("MATLAB:checkValidityProfileOptions:summarize_data_profilesNotValid", "The field of 'summarize_data_profiles' of options should be a boolean.");
        end
    end
    % Judge whether profile_options.summarize_log_ratio_profiles is a boolean.
    if isfield(profile_options, ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value))
            error("MATLAB:checkValidityProfileOptions:summarize_log_ratio_profilesNotValid", "The field 'summarize_log_ratio_profiles' of options should be a boolean.");
        end
        if profile_options.(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value) && numel(solvers) > 2
            warning("MATLAB:checkValidityProfileOptions:summarize_log_ratio_profilesOnlyWhenTwoSolvers", "The log-ratio profiles are available only when there are exactly two solvers. We will not generate the log-ratio profiles.");
            profile_options.(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value) = false;
        end
    end
    % Judge whether profile_options.summarize_output_based_profiles is a boolean.
    if isfield(profile_options, ProfileOptionKey.SUMMARIZE_OUTPUT_BASED_PROFILES.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.SUMMARIZE_OUTPUT_BASED_PROFILES.value))
            error("MATLAB:checkValidityProfileOptions:summarize_output_based_profilesNotValid", "The field 'summarize_output_based_profiles' of options should be a boolean.");
        end
    end
    % Judge whether profile_options.silent is a boolean.
    if isfield(profile_options, ProfileOptionKey.SILENT.value)
        if ~islogicalscalar(profile_options.(ProfileOptionKey.SILENT.value))
            error("MATLAB:checkValidityProfileOptions:silentNotValid", "The field 'silent' of options should be a boolean.");
        end
    end
    % Judge whether profile_options.solver_verbose is among 0, 1, and 2. If silent is true, solver_verbose should be 0 or 1 (only print errors). If it is 2, print a message saying that solver_verbose will be set to 1.
    if isfield(profile_options, ProfileOptionKey.SOLVER_VERBOSE.value)
        if ~isintegerscalar(profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value)) || ~ismember(profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value), [0, 1, 2])
            error("MATLAB:checkValidityProfileOptions:solver_verboseNotValid", "The field 'solver_verbose' of options should be either 0, 1, or 2.");
        end
        if profile_options.(ProfileOptionKey.SILENT.value) && profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) == 2
            warning("MATLAB:checkValidityProfileOptions:solver_verboseSetToZero", "solver_verbose will be set to 1 because silent is true.");
            profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) = 1;
        end
    end

end