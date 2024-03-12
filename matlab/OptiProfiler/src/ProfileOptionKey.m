classdef ProfileOptionKey
%PROFILEOPTIONKEY enumerates options for creating profiles
    
    enumeration
        N_JOBS ('n_jobs')
        BENCHMARK_ID ('benchmark_id')
        MAX_TOL_ORDER ('max_tol_order')
        MAX_EVAL_FACTOR ('max_eval_factor')
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