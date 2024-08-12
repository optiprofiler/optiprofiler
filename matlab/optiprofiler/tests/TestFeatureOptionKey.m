classdef TestFeatureOptionKey < matlab.unittest.TestCase
    methods (Test)

        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned

            enumValues = enumeration('FeatureOptionKey');
            enumValues = cellstr(arrayfun(@char, enumValues, 'UniformOutput', false));
            expectedValues = {'N_RUNS'; 'DISTRIBUTION'; 'NOISE_LEVEL'; 'NOISE_TYPE'; 'SIGNIFICANT_DIGITS'; 'ROTATED'; 'CONDITION_NUMBER'; 'UNRELAXABLE_BOUNDS'; 'UNRELAXABLE_LINEAR_CONSTRAINTS'; 'UNRELAXABLE_NONLINEAR_CONSTRAINTS'; 'PERTURBED_TRAILING_ZEROS'; 'RATE_NAN'; 'MODIFIER'};
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
            testCase.verifyEqual(FeatureOptionKey.ROTATED.value, 'rotated');
            testCase.verifyEqual(FeatureOptionKey.CONDITION_NUMBER.value, 'condition_number');
            testCase.verifyEqual(FeatureOptionKey.UNRELAXABLE_BOUNDS.value, 'unrelaxable_bounds');
            testCase.verifyEqual(FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value, 'unrelaxable_linear_constraints');
            testCase.verifyEqual(FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value, 'unrelaxable_nonlinear_constraints');
            testCase.verifyEqual(FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value, 'perturbed_trailing_zeros');
            testCase.verifyEqual(FeatureOptionKey.RATE_NAN.value, 'rate_nan');
            testCase.verifyEqual(FeatureOptionKey.MODIFIER.value, 'modifier');
        end

    end

end