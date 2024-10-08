classdef ProfileOptionKey
%PROFILEOPTIONKEY enumerates options for creating profiles
    
    enumeration
        N_JOBS ('n_jobs')
        BENCHMARK_ID ('benchmark_id')
        RANGE_TYPE ('range_type')
        STD_FACTOR ('std_factor')
        SAVEPATH ('savepath')
        MAX_TOL_ORDER ('max_tol_order')
        MAX_EVAL_FACTOR ('max_eval_factor')
        PROJECT_X0 ('project_x0')
        RUN_PLAIN ('run_plain')
        SUMMARIZE_PERFORMANCE_PROFILES ('summarize_performance_profiles')
        SUMMARIZE_DATA_PROFILES ('summarize_data_profiles')
        SUMMARIZE_LOG_RATIO_PROFILES ('summarize_log_ratio_profiles')
        SUMMARIZE_OUTPUT_BASED_PROFILES ('summarize_output_based_profiles')
    end
    properties
        value
    end
    methods
        function obj = ProfileOptionKey(inputValue)
            obj.value = inputValue;
        end
    end
end