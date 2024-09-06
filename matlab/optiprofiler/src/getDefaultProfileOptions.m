function profile_options = getDefaultProfileOptions()

    if exist('parcluster', 'file') == 2
        myCluster = parcluster('local');
        nb_cores = myCluster.NumWorkers;
    else
        nb_cores = 1;
    end
    profile_options.(ProfileOptionKey.N_JOBS.value) = nb_cores;
    profile_options.(ProfileOptionKey.BENCHMARK_ID.value) = '.';
    profile_options.(ProfileOptionKey.RANGE_TYPE.value) = 'minmax';
    profile_options.(ProfileOptionKey.STD_FACTOR.value) = 1;
    profile_options.(ProfileOptionKey.SAVEPATH.value) = pwd;
    profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) = 10;
    profile_options.(ProfileOptionKey.MAX_EVAL_FACTOR.value) = 500;
    profile_options.(ProfileOptionKey.PROJECT_X0.value) = false;
    profile_options.(ProfileOptionKey.RUN_PLAIN.value) = true;
    profile_options.(ProfileOptionKey.SUMMARIZE_PERFORMANCE_PROFILES.value) = true;
    profile_options.(ProfileOptionKey.SUMMARIZE_DATA_PROFILES.value) = true;
    profile_options.(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value) = false;
    profile_options.(ProfileOptionKey.SUMMARIZE_OUTPUT_BASED_PROFILES.value) = true;
    profile_options.(ProfileOptionKey.SUMMARIZE_FUNHIST.value) = true;
    profile_options.(ProfileOptionKey.SUMMARIZE_MAXCVHIST.value) = true;
    profile_options.(ProfileOptionKey.SUMMARIZE_MERITHIST.value) = true;
    profile_options.(ProfileOptionKey.SUMMARIZE_CUMMINHIST.value) = true;
end