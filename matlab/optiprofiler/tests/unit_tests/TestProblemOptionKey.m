classdef TestProblemOptionKey < matlab.unittest.TestCase
    methods (Test)

        function testEnumerationValues(testCase)
            % Test if the enumeration values are correctly assigned

            enumValues = enumeration('ProblemOptionKey');
            enumValues = cellstr(arrayfun(@char, enumValues, 'UniformOutput', false));
            expectedValues = {'PLIBS'; 'PTYPE'; 'MINDIM'; 'MAXDIM'; 'MINB'; 'MAXB'; 'MINLCON'; 'MAXLCON'; 'MINNLCON'; 'MAXNLCON'; 'MINCON'; 'MAXCON'; 'EXCLUDELIST'; 'PROBLEM_NAMES'};
            testCase.verifyEqual(enumValues, expectedValues);
        end

        function testConstructor(testCase)
            % Test if the constructor works as expected

            clear obj;
            clear ProblemOptionKey;
            testCase.verifyEqual(ProblemOptionKey.PLIBS.value, 'plibs');
            testCase.verifyEqual(ProblemOptionKey.PTYPE.value, 'ptype');
            testCase.verifyEqual(ProblemOptionKey.MINDIM.value, 'mindim');
            testCase.verifyEqual(ProblemOptionKey.MAXDIM.value, 'maxdim');
            testCase.verifyEqual(ProblemOptionKey.MINB.value, 'minb');
            testCase.verifyEqual(ProblemOptionKey.MAXB.value, 'maxb');
            testCase.verifyEqual(ProblemOptionKey.MINLCON.value, 'minlcon');
            testCase.verifyEqual(ProblemOptionKey.MAXLCON.value, 'maxlcon');
            testCase.verifyEqual(ProblemOptionKey.MINNLCON.value, 'minnlcon');
            testCase.verifyEqual(ProblemOptionKey.MAXNLCON.value, 'maxnlcon');
            testCase.verifyEqual(ProblemOptionKey.MINCON.value, 'mincon');
            testCase.verifyEqual(ProblemOptionKey.MAXCON.value, 'maxcon');
            testCase.verifyEqual(ProblemOptionKey.EXCLUDELIST.value, 'excludelist');
            testCase.verifyEqual(ProblemOptionKey.PROBLEM_NAMES.value, 'problem_names');
        end
    end
end