classdef TestProblemOptionKey < matlab.unittest.TestCase
    methods (Test)
        
        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned
            enumValues = enumeration('ProblemOptionKey');
            expectedValues = {'N_MIN', 'N_MAX', 'M_MIN', 'M_MAX'};

            diff1 = setdiff(enumValues, expectedValues);
            diff2 = setdiff(expectedValues, enumValues);
            testCase.verifyEmpty(diff1, "The enumeration values are not correctly assigned");
            testCase.verifyEmpty(diff2, "The enumeration values are not correctly assigned");
        end

        function testConstructor(testCase)
            % Test if the constructor works as expected
            testCase.verifyEqual(ProblemOptionKey.N_MIN.value, 'mindim');
            testCase.verifyEqual(ProblemOptionKey.N_MAX.value, 'maxdim');
            testCase.verifyEqual(ProblemOptionKey.M_MIN.value, 'mincon');
            testCase.verifyEqual(ProblemOptionKey.M_MAX.value, 'maxcon');
        end
        
    end

end