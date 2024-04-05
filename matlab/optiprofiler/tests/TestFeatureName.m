classdef TestFeatureName < matlab.unittest.TestCase
    methods (Test)
        
        function testEnumerationValues(testCase)

            % Test if the enumeration values are correctly assigned
            enumValues = enumeration('FeatureName');
            expectedValues = {'CUSTOM', 'NOISY', 'PLAIN', 'RANDOMIZE_X0', 'REGULARIZED', 'TOUGH', 'TRUNCATED'};

            diff1 = setdiff(enumValues, expectedValues);
            diff2 = setdiff(expectedValues, enumValues);
            testCase.verifyEmpty(diff1, "The enumeration values are not correctly assigned");
            testCase.verifyEmpty(diff2, "The enumeration values are not correctly assigned");
        end

        function testConstructor(testCase)
            % Test if the constructor works as expected
            testCase.verifyEqual(FeatureName.CUSTOM.value, 'custom');
            testCase.verifyEqual(FeatureName.NOISY.value, 'noisy');
            testCase.verifyEqual(FeatureName.PLAIN.value, 'plain');
            testCase.verifyEqual(FeatureName.RANDOMIZE_X0.value, 'randomize_x0');
            testCase.verifyEqual(FeatureName.REGULARIZED.value, 'regularized');
            testCase.verifyEqual(FeatureName.TOUGH.value, 'tough');
            testCase.verifyEqual(FeatureName.TRUNCATED.value, 'truncated');
        end
        
    end

end