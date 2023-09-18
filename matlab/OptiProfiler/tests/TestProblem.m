classdef TestProblem < matlab.unittest.TestCase
    methods(Test)

        % Test the constructor.

        function testConstructorWithEmpty(testCase)
            % Test that the constructor works correctly when no argument is provided.
            testCase.verifyError(@() Problem(), "MATLAB:Problem:MissingArguments");
        end

        function testConstructorWithEmptyStruct(testCase)
            % Test that the constructor works correctly when an empty struct is provided.
            testCase.verifyError(@() Problem(struct()), "MATLAB:Problem:MissingFields");
        end

        function testConstructorWithRequiredFields(testCase)
            % Test that the constructor works correctly when only the two required fields are provided.
            fun = @(x) x.' * x;
            x0 = [0, 0];
            s = struct('fun', fun, 'x0', x0);
            problemInstance = Problem(s);
            testCase.verifyEqual(problemInstance.fun([100; 100]), fun([100; 100]));
            testCase.verifyEqual(problemInstance.x0, [0; 0]);
        end
        
        function testConstructorWithAllFields(testCase)
            % Test that the constructor works correctly when all fields are provided.
            fun = @(x) x.' * x;
            x0 = [0, 0];
            xl = [-1, -1];
            xu = [1, 1];
            aub = [1 1; 1 1];
            bub = [1; 1];
            aeq = [1 1; 1 1];
            beq = [1; 1];
            cub = @(x) x.' * x;
            ceq = @(x) x.' * x;
            m_nonlinear_ub = 1;
            m_nonlinear_eq = 1;
            s = struct('fun', fun, 'x0', x0, 'xl', xl, 'xu', xu, 'aub', aub, 'bub', bub, 'aeq', aeq, 'beq', beq, 'cub', cub, 'ceq', ceq, 'm_nonlinear_ub', m_nonlinear_ub, 'm_nonlinear_eq', m_nonlinear_eq);
            problemInstance = Problem(s);
            testCase.verifyEqual(problemInstance.fun([0; 0]), 0);
            testCase.verifyEqual(problemInstance.cub([0; 0]), 0);
            testCase.verifyEqual(problemInstance.ceq([0; 0]), 0);
            testCase.verifyEqual(problemInstance.x0, [0; 0]);
            testCase.verifyEqual(problemInstance.xl, [-1; -1]);
            testCase.verifyEqual(problemInstance.xu, [1; 1]);
            testCase.verifyEqual(problemInstance.aub, [1 1; 1 1]);
            testCase.verifyEqual(problemInstance.bub, [1; 1]);
            testCase.verifyEqual(problemInstance.aeq, [1 1; 1 1]);
            testCase.verifyEqual(problemInstance.beq, [1; 1]);
            testCase.verifyEqual(problemInstance.m_nonlinear_ub, 1);
            testCase.verifyEqual(problemInstance.m_nonlinear_eq, 1);
        end

        function testConstructorWithNotStruct(testCase)
            % Test that the constructor works correctly when a non-struct is provided.
            testCase.verifyError(@() Problem(1), "MATLAB:Problem:NotStruct");
        end

        function testConstructorConsistent(testCase)
            % Test that the constructor works correctly in the case that some fields are not consistent.
            s1 = struct('fun', @(x) x.' * x, 'x0', [0; 0], 'xl', [1; 1; 1]);
            testCase.verifyError(@() Problem(s1), "MATLAB:Problem:xl_x0_NotConsistent");
            s2 = struct('fun', @(x) x.' * x, 'x0', [0; 0], 'xu', [1; 1; 1]);
            testCase.verifyError(@() Problem(s2), "MATLAB:Problem:xu_x0_NotConsistent");
            s3 = struct('fun', @(x) x.' * x, 'x0', [0; 0], 'aub', [1 1; 1 1], 'bub', [1; 1; 1]);
            testCase.verifyError(@() Problem(s3), "MATLAB:Problem:aub_m_linear_ub_n_NotConsistent");
            s4 = struct('fun', @(x) x.' * x, 'x0', [0; 0], 'aeq', [1 1; 1 1], 'beq', [1; 1; 1]);
            testCase.verifyError(@() Problem(s4), "MATLAB:Problem:aeq_m_linear_eq_n_NotConsistent");
            s5 = struct('fun', @(x) x.' * x, 'x0', [0; 0], 'm_nonlinear_ub', 2);
            testCase.verifyError(@() Problem(s5), "MATLAB:Problem:m_nonlinear_ub_cub_NotConsistent");
            s6 = struct('fun', @(x) x.' * x, 'x0', [0; 0], 'm_nonlinear_eq', 2);
            testCase.verifyError(@() Problem(s6), "MATLAB:Problem:m_nonlinear_eq_ceq_NotConsistent");
        end

    end
end

