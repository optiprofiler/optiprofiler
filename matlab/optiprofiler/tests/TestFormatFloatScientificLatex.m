classdef TestFormatFloatScientificLatex < matlab.unittest.TestCase
    methods (Test)
        
        function testValidInput(testCase)
            % Test whether the function works well.

            x = 0;
            [raw, formatted] = formatFloatScientificLatex(x);
            testCase.verifyEqual(raw, '0');
            testCase.verifyEqual(formatted, '0');

            x = 5;
            [raw, formatted] = formatFloatScientificLatex(x);
            testCase.verifyEqual(raw, '5');
            testCase.verifyEqual(formatted, '5');

            x = 10;
            [raw, formatted] = formatFloatScientificLatex(x);
            testCase.verifyEqual(raw, '10');
            testCase.verifyEqual(formatted, sprintf('10^{1}'));

            x = 20;
            [raw, formatted] = formatFloatScientificLatex(x);
            testCase.verifyEqual(raw, '2e1');
            testCase.verifyEqual(formatted, sprintf('2 \\times 10^{1}'));


        end

    end
end