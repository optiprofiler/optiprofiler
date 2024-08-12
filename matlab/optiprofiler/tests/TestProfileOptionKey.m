classdef TestProfileOptionKey < matlab.unittest.TestCase
    methods (Test)
        
        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned
            
            enumValues = enumeration('ProfileOptionKey');
            enumValues = cellstr(arrayfun(@char, enumValues, 'UniformOutput', false));
            expectedValues = {'N_JOBS'; 'BENCHMARK_ID'; 'RANGE_TYPE'; 'STD_FACTOR'; 'SAVEPATH'; 'MAX_TOL_ORDER'; 'MAX_EVAL_FACTOR'; 'PROJECT_X0'; 'RUN_PLAIN'; 'SUMMARIZE_PERFORMANCE_PROFILES'; 'SUMMARIZE_DATA_PROFILES'; 'SUMMARIZE_LOG_RATIO_PROFILES'};
            testCase.verifyEqual(enumValues, expectedValues);
        end

        function testConstructor(testCase)
            % Test if the constructor works as expected

            clear obj;
            clear ProfileOptionKey;
            testCase.verifyEqual(ProfileOptionKey.N_JOBS.value, 'n_jobs');
            testCase.verifyEqual(ProfileOptionKey.BENCHMARK_ID.value, 'benchmark_id');
            testCase.verifyEqual(ProfileOptionKey.RANGE_TYPE.value, 'range_type');
            testCase.verifyEqual(ProfileOptionKey.STD_FACTOR.value, 'std_factor');
            testCase.verifyEqual(ProfileOptionKey.SAVEPATH.value, 'savepath');
            testCase.verifyEqual(ProfileOptionKey.MAX_TOL_ORDER.value, 'max_tol_order');
            testCase.verifyEqual(ProfileOptionKey.MAX_EVAL_FACTOR.value, 'max_eval_factor');
            testCase.verifyEqual(ProfileOptionKey.PROJECT_X0.value, 'project_x0');
            testCase.verifyEqual(ProfileOptionKey.RUN_PLAIN.value, 'run_plain');
            testCase.verifyEqual(ProfileOptionKey.SUMMARIZE_PERFORMANCE_PROFILES.value, 'summarize_performance_profiles');
            testCase.verifyEqual(ProfileOptionKey.SUMMARIZE_DATA_PROFILES.value, 'summarize_data_profiles');
            testCase.verifyEqual(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value, 'summarize_log_ratio_profiles');
        end
        
    end

end