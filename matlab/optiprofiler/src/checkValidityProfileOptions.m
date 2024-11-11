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
            error("MATLAB:benchmark:n_jobsNotValid", "profile_options.n_jobs should be a integer.");
        elseif profile_options.(ProfileOptionKey.N_JOBS.value) < 1
            profile_options.(ProfileOptionKey.N_JOBS.value) = 1;
        elseif profile_options.(ProfileOptionKey.N_JOBS.value) > nb_cores
            profile_options.(ProfileOptionKey.N_JOBS.value) = nb_cores;
        else
            profile_options.(ProfileOptionKey.N_JOBS.value) = round(profile_options.(ProfileOptionKey.N_JOBS.value));
        end
    end
    % Judge whether profile_options.keep_pool is a boolean.
    if ~islogicalscalar(profile_options.(ProfileOptionKey.KEEP_POOL.value))
        error("MATLAB:benchmark:keep_poolNotValid", "keep_pool should be a boolean.");
    end
    % Judge whether profile_options.seed is a positive integer.
    if ~isintegerscalar(profile_options.(ProfileOptionKey.SEED.value)) || profile_options.(ProfileOptionKey.SEED.value) <= 0
        error("MATLAB:benchmark:seedNotValid", "seed should be a positive integer.");
    end
    % Judge whether profile_options.benchmark_id is a char or a string and satisfies the file name requirements (but it can be '.').
    is_valid_foldername = @(x) ischarstr(x) && ~isempty(x) && all(ismember(char(x), ['a':'z', 'A':'Z', '0':'9', '_', '-']));
    if ~ischarstr(profile_options.(ProfileOptionKey.BENCHMARK_ID.value)) || ~is_valid_foldername(profile_options.(ProfileOptionKey.BENCHMARK_ID.value)) && ~strcmp(profile_options.(ProfileOptionKey.BENCHMARK_ID.value), '.')
        error("MATLAB:benchmark:benchmark_idNotValid", "profile_options.benchmark_id should be a char or a string satisfying the strict file name requirements (only containing letters, numbers, '_', and '-').");
    end
    % Judge whether profile_options.range_type is among 'minmax' and 'meanstd'.
    if ~ischarstr(ProfileOptionKey.RANGE_TYPE.value) || ~ismember(char(profile_options.(ProfileOptionKey.RANGE_TYPE.value)), {'minmax', 'meanstd'})
        error("MATLAB:benchmark:range_typeNotValid", "range_type should be either 'minmax' or 'meanstd'.");
    end
    % Judge whether profile_options.std_factor is a positive real number.
    if ~isrealscalar(profile_options.(ProfileOptionKey.STD_FACTOR.value)) || profile_options.(ProfileOptionKey.STD_FACTOR.value) <= 0
        error("MATLAB:benchmark:std_factorNotValid", "std_factor should be a positive real number.");
    end
    % Judge whether profile_options.savepath is a string and exists. If not exists, create it.
    if ~ischarstr(profile_options.(ProfileOptionKey.SAVEPATH.value))
        error("MATLAB:benchmark:savepathNotValid", "savepath should be a char or a string.");
    elseif ~exist(profile_options.(ProfileOptionKey.SAVEPATH.value), 'dir')
        status = mkdir(profile_options.(ProfileOptionKey.SAVEPATH.value));
        if ~status
            error("MATLAB:benchmark:savepathNotExist", "profile_options.savepath does not exist and cannot be created.");
        end
    end
    % Judge whether profile_options.max_tol_order is a positive integer.
    if ~isintegerscalar(profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value)) || profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) <= 0
        error("MATLAB:benchmark:max_tol_orderNotValid", "max_tol_order should be a positive integer.");
    end
    % Judge whether profile_options.max_eval_factor is a positive integer.
    if ~isintegerscalar(profile_options.(ProfileOptionKey.MAX_EVAL_FACTOR.value)) || profile_options.(ProfileOptionKey.MAX_EVAL_FACTOR.value) <= 0
        error("MATLAB:benchmark:max_eval_factorNotValid", "max_eval_factor should be a positive integer.");
    end
    % Judge whether profile_options.project_x0 is a boolean.
    if ~islogicalscalar(profile_options.(ProfileOptionKey.PROJECT_X0.value))
        error("MATLAB:benchmark:project_x0NotValid", "project_x0 should be a boolean.");
    end
    % Judge whether profile_options.run_plain is a boolean.
    if ~islogicalscalar(profile_options.(ProfileOptionKey.RUN_PLAIN.value))
        error("MATLAB:benchmark:run_plainNotValid", "run_plain should be a boolean.");
    end
    % Judge whether profile_options.summarize_performance_profiles is a boolean.
    if ~islogicalscalar(profile_options.(ProfileOptionKey.SUMMARIZE_PERFORMANCE_PROFILES.value))
        error("MATLAB:benchmark:summarize_performance_profilesNotValid", "summarize_performance_profiles should be a boolean.");
    end
    % Judge whether profile_options.summarize_data_profiles is a boolean.
    if ~islogicalscalar(profile_options.(ProfileOptionKey.SUMMARIZE_DATA_PROFILES.value))
        error("MATLAB:benchmark:summarize_data_profilesNotValid", "summarize_data_profiles should be a boolean.");
    end
    % Judge whether profile_options.summarize_log_ratio_profiles is a boolean.
    if ~islogicalscalar(profile_options.(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value))
        error("MATLAB:benchmark:summarize_log_ratio_profilesNotValid", "summarize_log_ratio_profiles should be a boolean.");
    end
    if profile_options.(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value) && numel(solvers) > 2
        warning("MATLAB:benchmark:summarize_log_ratio_profilesOnlyWhenTwoSolvers", "The log-ratio profiles are available only when there are exactly two solvers. We will not generate the log-ratio profiles.");
        profile_options.(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value) = false;
    end
    % Judge whether profile_options.summarize_output_based_profiles is a boolean.
    if ~islogicalscalar(profile_options.(ProfileOptionKey.SUMMARIZE_OUTPUT_BASED_PROFILES.value))
        error("MATLAB:benchmark:summarize_output_based_profilesNotValid", "summarize_output_based_profiles should be a boolean.");
    end
    % Judge whether profile_options.silent is a boolean.
    if ~islogicalscalar(profile_options.(ProfileOptionKey.SILENT.value))
        error("MATLAB:benchmark:silentNotValid", "silent should be a boolean.");
    end
    % Judge whether profile_options.solver_verbose is among 0, 1, and 2. If silent is true, solver_verbose should be 0 or 1 (only print errors). If it is 2, print a message saying that solver_verbose will be set to 1.
    if ~isintegerscalar(profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value)) || ~ismember(profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value), [0, 1, 2])
        error("MATLAB:benchmark:solver_verboseNotValid", "solver_verbose should be either 0, 1, or 2.");
    end
    if profile_options.(ProfileOptionKey.SILENT.value) && profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) == 2
        warning("MATLAB:benchmark:solver_verboseSetToZero", "solver_verbose will be set to 1 because silent is true.");
        profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) = 1;
    end

end