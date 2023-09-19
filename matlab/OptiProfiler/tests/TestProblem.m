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


        % Test the setters.

        function testSetter_x0(testCase)
            % Test that the error is thrown when the x0 is not a vector.
            fun = @(x) x.' * x;
            x0 = [0 0; 0 0];
            s = struct('fun', fun, 'x0', x0);
            testCase.verifyError(@() Problem(s), "MATLAB:Problem:x0_NotVector");
        end

        function testSetter_fun_(testCase)
            % Test that the errors are thrown when the fun is not provided correctly.
            s1 = struct('fun', 1, 'x0', 0);
            testCase.verifyError(@() Problem(s1), "MATLAB:Problem:fun_NotFunctionHandle");
            s2 = struct('fun', @(x, y) x.' * y, 'x0', 0);
            testCase.verifyError(@() Problem(s2), "MATLAB:Problem:fun_NotOneArguementFun");
        end

        function testSetter_xl(testCase)
            % Test that the error is thrown when the xl is not a vector.
            fun = @(x) x.' * x;
            x0 = [0; 0];
            xl = [0 0; 0 0];
            s = struct('fun', fun, 'x0', x0, 'xl', xl);
            testCase.verifyError(@() Problem(s), "MATLAB:Problem:xl_NotVector");
        end

        function testSetter_xu(testCase)
            % Test that the error is thrown when the xu is not a vector.
            fun = @(x) x.' * x;
            x0 = [0; 0];
            xu = [0 0; 0 0];
            s = struct('fun', fun, 'x0', x0, 'xu', xu);
            testCase.verifyError(@() Problem(s), "MATLAB:Problem:xu_NotVector");
        end

        function testSetter_aub(testCase)
            % Test that the error is thrown when the aub is not a matrix.
            fun = @(x) x.' * x;
            x0 = [0; 0];
            aub = ones(2, 2, 3);
            s = struct('fun', fun, 'x0', x0, 'aub', aub);
            testCase.verifyError(@() Problem(s), "MATLAB:Problem:aub_NotMatrix");
        end

        function testSetter_bub(testCase)
            % Test that the error is thrown when the bub is not a vector.
            fun = @(x) x.' * x;
            x0 = [0; 0];
            bub = ones(2, 2, 3);
            s = struct('fun', fun, 'x0', x0, 'bub', bub);
            testCase.verifyError(@() Problem(s), "MATLAB:Problem:bub_NotVector");
        end

        function testSetter_aeq(testCase)
            % Test that the error is thrown when the aeq is not a matrix.
            fun = @(x) x.' * x;
            x0 = [0; 0];
            aeq = ones(2, 2, 3);
            s = struct('fun', fun, 'x0', x0, 'aeq', aeq);
            testCase.verifyError(@() Problem(s), "MATLAB:Problem:aeq_NotMatrix");
        end

        function testSetter_beq(testCase)
            % Test that the error is thrown when the beq is not a vector.
            fun = @(x) x.' * x;
            x0 = [0; 0];
            beq = ones(2, 2, 3);
            s = struct('fun', fun, 'x0', x0, 'beq', beq);
            testCase.verifyError(@() Problem(s), "MATLAB:Problem:beq_NotVector");
        end

        function testSetter_cub_(testCase)
            % Test that the error is thrown when cub is not a function handle.
            s = struct('fun', @(x) x.' * x, 'x0', 0, 'cub', 1);
            testCase.verifyError(@() Problem(s), "MATLAB:Problem:cub_NotFunctionHandle");
        end

        function testSetter_ceq_(testCase)
            % Test that the error is thrown when ceq is not a function handle.
            s = struct('fun', @(x) x.' * x, 'x0', 0, 'ceq', 1);
            testCase.verifyError(@() Problem(s), "MATLAB:Problem:ceq_NotFunctionHandle");
        end

        function testSetter_m_nonlinear_ub(testCase)
            % Test that the error is thrown when m_nonlinear_ub is not a scalar.
            fun = @(x) x.' * x;
            x0 = [0; 0];
            m_nonlinear_ub = -2;
            s = struct('fun', fun, 'x0', x0, 'm_nonlinear_ub', m_nonlinear_ub);
            testCase.verifyError(@() Problem(s), "MATLAB:Problem:m_nonlinear_ub_NotPositiveScalar");
        end

        function testSetter_m_nonlinear_eq(testCase)
            % Test that the error is thrown when m_nonlinear_eq is not a scalar.
            fun = @(x) x.' * x;
            x0 = [0; 0];
            m_nonlinear_eq = -2;
            s = struct('fun', fun, 'x0', x0, 'm_nonlinear_eq', m_nonlinear_eq);
            testCase.verifyError(@() Problem(s), "MATLAB:Problem:m_nonlinear_eq_NotPositiveScalar");
        end

    end
end

