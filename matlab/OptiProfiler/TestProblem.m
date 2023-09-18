classdef TestProblem < matlab.unittest.TestCase
    methods(Test)

        function testConstructorWithRequiredFields(testCase)
            % Test that the constructor works correctly when only the required fields are provided.
            s = struct('fun', @(x) x.^2, 'x0', 0);
            problemInstance = Problem(s);
            testCase.verifyEqual(problemInstance.fun(0), 0);  % Verify that the `fun` field was correctly set
        end
        
        function testConstructorWithAllFields(testCase)
            % Test that the constructor works correctly when all fields are provided.
            s = struct('fun', @(x) x^2, 'x0', 0, 'cub', @(x) x-1, 'ceq', @(x) x);
            problemInstance = Problem(s);
            testCase.verifyEqual(problemInstance.fun(0), 0);  % Verify that the `fun` field was correctly set
            testCase.verifyEqual(problemInstance.cub(0), -1);  % Verify that the `cub` field was correctly set
            testCase.verifyEqual(problemInstance.ceq(0), 0);  % Verify that the `ceq` field was correctly set
        end

    end
end

