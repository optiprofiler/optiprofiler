classdef TestFeatureOptionKey < matlab.unittest.TestCase
    methods (Test)

        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned

            enumValues = enumeration('FeatureOptionKey');
            enumValues = cellstr(arrayfun(@char, enumValues, 'UniformOutput', false));
            expectedValues = {'N_RUNS'; 'DISTRIBUTION'; 'NOISE_LEVEL'; 'NOISE_TYPE'; 'SIGNIFICANT_DIGITS'; 'PERTURBED_TRAILING_ZEROS'; 'ROTATED'; 'CONDITION_FACTOR'; 'RATE_NAN'; 'UNRELAXABLE_BOUNDS'; 'UNRELAXABLE_LINEAR_CONSTRAINTS'; 'UNRELAXABLE_NONLINEAR_CONSTRAINTS'; 'MESH_SIZE'; 'IS_TRUTH'; 'MOD_X0'; 'MOD_AFFINE'; 'MOD_BOUNDS'; 'MOD_LINEAR_UB'; 'MOD_LINEAR_EQ'; 'MOD_FUN'; 'MOD_CUB'; 'MOD_CEQ'};
            testCase.verifyEqual(enumValues, expectedValues);
        end

        function testConstructor(testCase)
            % Test if the constructor works as expected

            clear obj;
            clear FeatureOptionKey;
            testCase.verifyEqual(FeatureOptionKey.N_RUNS.value, 'n_runs');
            testCase.verifyEqual(FeatureOptionKey.DISTRIBUTION.value, 'distribution');
            testCase.verifyEqual(FeatureOptionKey.NOISE_LEVEL.value, 'noise_level');
            testCase.verifyEqual(FeatureOptionKey.NOISE_TYPE.value, 'noise_type');
            testCase.verifyEqual(FeatureOptionKey.SIGNIFICANT_DIGITS.value, 'significant_digits');
            testCase.verifyEqual(FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value, 'perturbed_trailing_zeros');
            testCase.verifyEqual(FeatureOptionKey.ROTATED.value, 'rotated');
            testCase.verifyEqual(FeatureOptionKey.CONDITION_FACTOR.value, 'condition_factor');
            testCase.verifyEqual(FeatureOptionKey.RATE_NAN.value, 'rate_nan');
            testCase.verifyEqual(FeatureOptionKey.UNRELAXABLE_BOUNDS.value, 'unrelaxable_bounds');
            testCase.verifyEqual(FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value, 'unrelaxable_linear_constraints');
            testCase.verifyEqual(FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value, 'unrelaxable_nonlinear_constraints');
            testCase.verifyEqual(FeatureOptionKey.MESH_SIZE.value, 'mesh_size');
            testCase.verifyEqual(FeatureOptionKey.IS_TRUTH.value, 'is_truth');
            testCase.verifyEqual(FeatureOptionKey.MOD_X0.value, 'mod_x0');
            testCase.verifyEqual(FeatureOptionKey.MOD_AFFINE.value, 'mod_affine');
            testCase.verifyEqual(FeatureOptionKey.MOD_BOUNDS.value, 'mod_bounds');
            testCase.verifyEqual(FeatureOptionKey.MOD_LINEAR_UB.value, 'mod_linear_ub');
            testCase.verifyEqual(FeatureOptionKey.MOD_LINEAR_EQ.value, 'mod_linear_eq');
            testCase.verifyEqual(FeatureOptionKey.MOD_FUN.value, 'mod_fun');
            testCase.verifyEqual(FeatureOptionKey.MOD_CUB.value, 'mod_cub');
            testCase.verifyEqual(FeatureOptionKey.MOD_CEQ.value, 'mod_ceq');
        end

    end

end