classdef TestFeatureOptionKey < matlab.unittest.TestCase
    methods (Test)

        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned

            enumValues = enumeration('FeatureOptionKey');
            enumValues = cellstr(arrayfun(@char, enumValues, 'UniformOutput', false));
            expectedValues = {'N_RUNS'; 'DISTRIBUTION'; 'PERTURBATION_LEVEL'; 'NOISE_LEVEL'; 'NOISE_TYPE'; 'SIGNIFICANT_DIGITS'; 'PERTURBED_TRAILING_DIGITS'; 'ROTATED'; 'CONDITION_FACTOR'; 'NAN_RATE'; 'UNRELAXABLE_BOUNDS'; 'UNRELAXABLE_LINEAR_CONSTRAINTS'; 'UNRELAXABLE_NONLINEAR_CONSTRAINTS'; 'MESH_SIZE'; 'MESH_TYPE'; 'GROUND_TRUTH'; 'MOD_X0'; 'MOD_AFFINE'; 'MOD_BOUNDS'; 'MOD_LINEAR_UB'; 'MOD_LINEAR_EQ'; 'MOD_FUN'; 'MOD_CUB'; 'MOD_CEQ'};
            testCase.verifyEqual(enumValues, expectedValues);
        end

        function testConstructor(testCase)
            % Test if the constructor works as expected

            clear obj;
            clear FeatureOptionKey;
            testCase.verifyEqual(FeatureOptionKey.N_RUNS.value, 'n_runs');
            testCase.verifyEqual(FeatureOptionKey.DISTRIBUTION.value, 'distribution');
            testCase.verifyEqual(FeatureOptionKey.PERTURBATION_LEVEL.value, 'perturbation_level');
            testCase.verifyEqual(FeatureOptionKey.NOISE_LEVEL.value, 'noise_level');
            testCase.verifyEqual(FeatureOptionKey.NOISE_TYPE.value, 'noise_type');
            testCase.verifyEqual(FeatureOptionKey.SIGNIFICANT_DIGITS.value, 'significant_digits');
            testCase.verifyEqual(FeatureOptionKey.PERTURBED_TRAILING_DIGITS.value, 'perturbed_trailing_digits');
            testCase.verifyEqual(FeatureOptionKey.ROTATED.value, 'rotated');
            testCase.verifyEqual(FeatureOptionKey.CONDITION_FACTOR.value, 'condition_factor');
            testCase.verifyEqual(FeatureOptionKey.NAN_RATE.value, 'nan_rate');
            testCase.verifyEqual(FeatureOptionKey.UNRELAXABLE_BOUNDS.value, 'unrelaxable_bounds');
            testCase.verifyEqual(FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value, 'unrelaxable_linear_constraints');
            testCase.verifyEqual(FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value, 'unrelaxable_nonlinear_constraints');
            testCase.verifyEqual(FeatureOptionKey.MESH_SIZE.value, 'mesh_size');
            testCase.verifyEqual(FeatureOptionKey.MESH_TYPE.value, 'mesh_type');
            testCase.verifyEqual(FeatureOptionKey.GROUND_TRUTH.value, 'ground_truth');
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