classdef TestFeatureOptionKey < matlab.unittest.TestCase
    methods (Test)

        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned
            enumValues = enumeration('FeatureOptionKey');
            expectedValues = {'DISTRIBUTION', 'MODIFIER', 'N_RUNS', 'ORDER', 'PARAMETER', 'RATE_ERROR', 'RATE_NAN', 'SIGNIFICANT_DIGITS', 'TYPE'};

            diff1 = setdiff(enumValues, expectedValues);
            diff2 = setdiff(expectedValues, enumValues);
            testCase.verifyEmpty(diff1, "The enumeration values are not correctly assigned");
            testCase.verifyEmpty(diff2, "The enumeration values are not correctly assigned");
        end

        function testConstructor(testCase)
            % Test if the constructor works as expected
            testCase.verifyEqual(FeatureOptionKey.DISTRIBUTION.value, 'distribution');
            testCase.verifyEqual(FeatureOptionKey.MODIFIER.value, 'modifier');
            testCase.verifyEqual(FeatureOptionKey.N_RUNS.value, 'n_runs');
            testCase.verifyEqual(FeatureOptionKey.ORDER.value, 'order');
            testCase.verifyEqual(FeatureOptionKey.PARAMETER.value, 'parameter');
            testCase.verifyEqual(FeatureOptionKey.RATE_ERROR.value, 'rate_error');
            testCase.verifyEqual(FeatureOptionKey.RATE_NAN.value, 'rate_nan');
            testCase.verifyEqual(FeatureOptionKey.SIGNIFICANT_DIGITS.value, 'significant_digits');
            testCase.verifyEqual(FeatureOptionKey.TYPE.value, 'type');
        end

    end

end