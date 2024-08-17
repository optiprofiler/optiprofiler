classdef TestGetPerformanceDataProfileAxes < matlab.unittest.TestCase
    methods (Test)
        
        function testWithValidInput(testCase)
            % Test whether the function works well with valid input.

            work = randn(3, 2, 5);
            problem_dimensions = [1, 2, 3];
            denominator_perf = @(i_problem, i_run) min(work(i_problem, :, i_run), [], 'omitnan');
            perf_or_data = 'perf';
            [x_perf, y_perf, ratio_max_perf] = getPerformanceDataProfileAxes(work, denominator_perf, perf_or_data);
            testCase.verifyEqual(size(x_perf), [15, 2]);
            testCase.verifyEqual(size(y_perf), [15, 2, 5]);
            testCase.verifyEqual(size(ratio_max_perf), [1, 1]);

            denominator_data = @(i_problem, i_run) problem_dimensions(i_problem) + 1;
            perf_or_data = 'data';
            [x_data, y_data, ratio_max_data] = getPerformanceDataProfileAxes(work, denominator_data, perf_or_data);
            testCase.verifyEqual(size(x_data), [15, 2]);
            testCase.verifyEqual(size(y_data), [15, 2, 5]);
            testCase.verifyEqual(size(ratio_max_data), [1, 1]);

            work = NaN(3, 2, 5);
            [~, ~, ratio_max] = getPerformanceDataProfileAxes(work, denominator_perf, perf_or_data);
            testCase.verifyEqual(ratio_max, eps);
        end

        function testErrors(testCase)
            % Test whether the function throws errors as expected.

            work = randn(3, 2, 5);
            problem_dimensions = [1, 2, 3];
            denominator_perf = @(i_problem, i_run) min(work(i_problem, :, i_run), [], 'omitnan');
            perf_or_data = 'unknown';
            testCase.verifyError(@() getPerformanceDataProfileAxes(work, denominator_perf, perf_or_data), "MATLAB:getPerformanceDataProfileAxes:UnknownNameForgetPerformanceDataProfileAxes")
            
        end
    end
end