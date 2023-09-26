classdef TestFeatureName < matlab.unittest.TestCase
    methods (Test)
        
        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned
            enumValues = enumeration('FeatureName');
            expectedValues = {'custom', 'noisy', 'plain', 'regularized', 'tough', 'truncated'};

            for i = 1:numel(enumValues)
                testCase.verifyEqual(enumValues(i).value, expectedValues{i});
            end
        end
        
    end

end