classdef TestProblemOptionKey < matlab.unittest.TestCase
    methods (Test)
        
        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned
            enumValues = enumeration('ProblemOptionKey');
            expectedValues = {'mindim', 'maxdim', 'mincon', 'maxcon'};

            for i = 1:numel(enumValues)
                testCase.verifyEqual(enumValues(i).value, expectedValues{i});
            end
        end
        
    end

end