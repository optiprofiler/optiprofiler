classdef TestProfileOptionKey < matlab.unittest.TestCase
    methods (Test)
        
        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned
            enumValues = enumeration('ProfileOptionKey');
            expectedValues = {'n_jobs'};

            for i = 1:numel(enumValues)
                testCase.verifyEqual(enumValues(i).value, expectedValues{i});
            end
        end
        
    end

end