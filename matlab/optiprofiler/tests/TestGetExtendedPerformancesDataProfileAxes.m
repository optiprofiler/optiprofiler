classdef TestGetExtendedPerformancesDataProfileAxes < matlab.unittest.TestCase
    methods (Test)
        
        function testWithValidInput(testCase)
            % Test whether the function works well with valid input.

            work = rand(3, 2, 5);
            problem_dimensions = [1, 2, 3];
            [x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data] = getExtendedPerformancesDataProfileAxes(work, problem_dimensions);
            testCase.verifyEqual(size(x_perf), [17, 2]);
            testCase.verifyEqual(size(y_perf), [17, 2, 5]);
            testCase.verifyEqual(size(ratio_max_perf), [1, 1]);
            testCase.verifyEqual(size(x_data), [17, 2]);
            testCase.verifyEqual(size(y_data), [17, 2, 5]);
            testCase.verifyEqual(size(ratio_max_data), [1, 1]);
        end
    end
end