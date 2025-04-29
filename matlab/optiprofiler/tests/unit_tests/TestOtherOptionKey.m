classdef TestOtherOptionKey < matlab.unittest.TestCase
    methods (Test)
        
        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned
            
            enumValues = enumeration('OtherOptionKey');
            enumValues = cellstr(arrayfun(@char, enumValues, 'UniformOutput', false));
            expectedValues = {'PROBLEM'; 'CUTEST_PROBLEM_NAMES'; 'CUSTOM_PROBLEM_LOADER'; 'CUSTOM_PROBLEM_NAMES'};
            testCase.verifyEqual(enumValues, expectedValues);
        end

        function testConstructor(testCase)
            % Test if the constructor works as expected
            
            clear obj;
            clear OtherOptionKey;
            testCase.verifyEqual(OtherOptionKey.PROBLEM.value, 'problem');
            testCase.verifyEqual(OtherOptionKey.CUTEST_PROBLEM_NAMES.value, 'cutest_problem_names');
            testCase.verifyEqual(OtherOptionKey.CUSTOM_PROBLEM_LOADER.value, 'custom_problem_loader');
            testCase.verifyEqual(OtherOptionKey.CUSTOM_PROBLEM_NAMES.value, 'custom_problem_names');
        end
        
    end

end