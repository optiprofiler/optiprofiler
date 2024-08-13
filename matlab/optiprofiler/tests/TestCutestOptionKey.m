classdef TestCutestOptionKey < matlab.unittest.TestCase
    methods (Test)

        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned

            enumValues = enumeration('CutestOptionKey');
            enumValues = cellstr(arrayfun(@char, enumValues, 'UniformOutput', false));
            expectedValues = {'PROBLEM_TYPE'; 'MINDIM'; 'MAXDIM'; 'MINCON'; 'MAXCON'; 'EXCLUDELIST'};
            testCase.verifyEqual(enumValues, expectedValues);
        end

        function testConstructor(testCase)
            % Test if the constructor works as expected

            clear obj;
            clear CutestOptionKey;
            testCase.verifyEqual(CutestOptionKey.PROBLEM_TYPE.value, 'problem_type');
            testCase.verifyEqual(CutestOptionKey.MINDIM.value, 'mindim');
            testCase.verifyEqual(CutestOptionKey.MAXDIM.value, 'maxdim');
            testCase.verifyEqual(CutestOptionKey.MINCON.value, 'mincon');
            testCase.verifyEqual(CutestOptionKey.MAXCON.value, 'maxcon');
            testCase.verifyEqual(CutestOptionKey.EXCLUDELIST.value, 'excludelist');
        end
    end
end