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


        % Test the getters.

        function testGetter_n(testCase)
            % Test that obj.n returns the correct value.
            s = struct('fun', @(x) x.' * x, 'x0', zeros(10, 1));
            problemInstance = Problem(s);
            testCase.verifyEqual(problemInstance.n, 10);
        end

        function testGetter_m_linear_ub(testCase)
            % Test that obj.m_linear_ub returns the correct value.
            s = struct('fun', @(x) x.' * x, 'x0', zeros(10, 1), 'aub', zeros(5, 10), 'bub', zeros(5, 1));
            problemInstance = Problem(s);
            testCase.verifyEqual(problemInstance.m_linear_ub, 5);
        end

        function testGetter_m_linear_eq(testCase)
            % Test that obj.m_linear_eq returns the correct value.
            s = struct('fun', @(x) x.' * x, 'x0', zeros(10, 1), 'aeq', zeros(5, 10), 'beq', zeros(5, 1));
            problemInstance = Problem(s);
            testCase.verifyEqual(problemInstance.m_linear_eq, 5);
        end

        function testGetter_fun(testCase)
            % Test that obj.fun(x) returns the correct value.
            fun = @(x) x.' * x;
            s = struct('fun', fun, 'x0', [0 0]);
            problemInstance = Problem(s);
            testCase.verifyEqual(problemInstance.fun([1; 1]), fun([1; 1]));
        end

        function testGetter_cub(testCase)
            % Test that obj.cub(x) returns the correct value.
            cub = @(x) x.' * x;
            s = struct('fun', @(x) x.' * x, 'x0', [0 0], 'cub', cub);
            problemInstance = Problem(s);
            testCase.verifyEqual(problemInstance.cub([1; 1]), cub([1; 1]));
        end

        function testGetter_ceq(testCase)
            % Test that obj.ceq(x) returns the correct value.
            ceq = @(x) x.' * x;
            s = struct('fun', @(x) x.' * x, 'x0', [0 0], 'ceq', ceq);
            problemInstance = Problem(s);
            testCase.verifyEqual(problemInstance.ceq([1; 1]), ceq([1; 1]));
        end

        function testGetter_m_nonlinear_ub(testCase)
            % Test that obj.m_nonlinear_ub returns the correct value.
            s1 = struct('fun', @(x) x.' * x, 'x0', [0 0]);
            problemInstance1 = Problem(s1);
            testCase.verifyEqual(problemInstance1.m_nonlinear_ub, 0);
            s2 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'cub', @(x) x.' * x);
            problemInstance2 = Problem(s2);
            testCase.verifyEqual(problemInstance2.m_nonlinear_ub, []);
            s3 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'cub', @(x) x.' * x, 'm_nonlinear_ub', 1);
            problemInstance3 = Problem(s3);
            testCase.verifyEqual(problemInstance3.m_nonlinear_ub, 1);
        end
        
        function testGetter_m_nonlinear_eq(testCase)
            % Test that obj.m_nonlinear_eq returns the correct value.
            s1 = struct('fun', @(x) x.' * x, 'x0', [0 0]);
            problemInstance1 = Problem(s1);
            testCase.verifyEqual(problemInstance1.m_nonlinear_eq, 0);
            s2 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'ceq', @(x) x.' * x);
            problemInstance2 = Problem(s2);
            testCase.verifyEqual(problemInstance2.m_nonlinear_eq, []);
            s3 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'ceq', @(x) x.' * x, 'm_nonlinear_eq', 1);
            problemInstance3 = Problem(s3);
            testCase.verifyEqual(problemInstance3.m_nonlinear_eq, 1);
        end

        function testGetter_xl(testCase)
            % Test that obj.xl returns the correct value.
            s1 = struct('fun', @(x) x.' * x, 'x0', [0 0]);
            problemInstance1 = Problem(s1);
            testCase.verifyEqual(problemInstance1.xl, -inf(2, 1));
            s2 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'xl', [1 1]);
            problemInstance2 = Problem(s2);
            testCase.verifyEqual(problemInstance2.xl, [1; 1]);
        end

        function testGetter_xu(testCase)
            % Test that obj.xu returns the correct value.
            s1 = struct('fun', @(x) x.' * x, 'x0', [0 0]);
            problemInstance1 = Problem(s1);
            testCase.verifyEqual(problemInstance1.xu, inf(2, 1));
            s2 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'xu', [1 1]);
            problemInstance2 = Problem(s2);
            testCase.verifyEqual(problemInstance2.xu, [1; 1]);
        end

        function testGetter_aub(testCase)
            % Test that obj.aub returns the correct value.
            s1 = struct('fun', @(x) x.' * x, 'x0', [0 0]);
            problemInstance1 = Problem(s1);
            testCase.verifyEqual(problemInstance1.aub, NaN(0, 2));
            s2 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'aub', [1 1; 1 1], 'bub', [1; 1]);
            problemInstance2 = Problem(s2);
            testCase.verifyEqual(problemInstance2.aub, [1 1; 1 1]);
        end

        function testGetter_bub(testCase)
            % Test that obj.bub returns the correct value.
            s1 = struct('fun', @(x) x.' * x, 'x0', [0 0]);
            problemInstance1 = Problem(s1);
            testCase.verifyEqual(problemInstance1.bub, NaN(0, 1));
            s2 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'aub', [1 1; 1 1], 'bub', [1; 1]);
            problemInstance2 = Problem(s2);
            testCase.verifyEqual(problemInstance2.bub, [1; 1]);
        end

        function testGetter_aeq(testCase)
            % Test that obj.aeq returns the correct value.
            s1 = struct('fun', @(x) x.' * x, 'x0', [0 0]);
            problemInstance1 = Problem(s1);
            testCase.verifyEqual(problemInstance1.aeq, NaN(0, 2));
            s2 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'aeq', [1 1; 1 1], 'beq', [1; 1]);
            problemInstance2 = Problem(s2);
            testCase.verifyEqual(problemInstance2.aeq, [1 1; 1 1]);
        end

        function testGetter_beq(testCase)
            % Test that obj.beq returns the correct value.
            s1 = struct('fun', @(x) x.' * x, 'x0', [0 0]);
            problemInstance1 = Problem(s1);
            testCase.verifyEqual(problemInstance1.beq, NaN(0, 1));
            s2 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'aeq', [1 1; 1 1], 'beq', [1; 1]);
            problemInstance2 = Problem(s2);
            testCase.verifyEqual(problemInstance2.beq, [1; 1]);
        end


        % Test other functions.

        function testMaxCV(testCase)
            % Test that obj.maxCV(x) returns the correct value.
            s1 = struct('fun', @(x) x.' * x, 'x0', [0 0]);
            problemInstance1 = Problem(s1);
            testCase.verifyEqual(problemInstance1.maxcv([10; 10]), 0);
            s2 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'xl', [10; 10]);
            problemInstance2 = Problem(s2);
            testCase.verifyEqual(problemInstance2.maxcv([0; 0]), 10);
            s3 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'xu', [-10; -10]);
            problemInstance3 = Problem(s3);
            testCase.verifyEqual(problemInstance3.maxcv([0; 0]), 10);
            s4 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'aub', [1 1; 1 1], 'bub', [-1; -1]);
            problemInstance4 = Problem(s4);
            testCase.verifyEqual(problemInstance4.maxcv([0; 0]), 1);
            s5 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'aeq', [1 1; 1 1], 'beq', [1; 1]);
            problemInstance5 = Problem(s5);
            testCase.verifyEqual(problemInstance5.maxcv([0; 0]), 1);
            s6 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'cub', @(x) x.' * x, 'm_nonlinear_ub', 1);
            problemInstance6 = Problem(s6);
            testCase.verifyEqual(problemInstance6.maxcv([1; 1]), 2);
            s7 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'ceq', @(x) x.' * x, 'm_nonlinear_eq', 1);
            problemInstance7 = Problem(s7);
            testCase.verifyEqual(problemInstance7.maxcv([1; 1]), 2);
        end


        % Test Private methods

        function testPrivateFUN(testCase)
            % Test that FUN works as expected.
            s = struct('fun', @(x) x.' * x, 'x0', [0 0]);
            problemInstance = Problem(s);
            testCase.verifyError(@() problemInstance.fun(ones(2, 2, 3)), "MATLAB:Problem:InvalidInputForFUN");
            testCase.verifyError(@() problemInstance.fun(ones(3, 1)), "MATLAB:Problem:WrongSizeInputForFUN");
        end

        function testPrivateCUB(testCase)
            % Test that CUB works as expected.
            s1 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'cub', @(x) x.' * x, 'm_nonlinear_ub', 1);
            problemInstance1 = Problem(s1);
            testCase.verifyError(@() problemInstance1.cub(ones(2, 2, 3)), "MATLAB:Problem:InvalidInputForCUB");
            testCase.verifyError(@() problemInstance1.cub(ones(3, 1)), "MATLAB:Problem:WrongSizeInputForCUB");
            testCase.verifyError(@() problemInstance1.cub([0 0]), "MATLAB:Problem:InvalidOutputForCUB");

            s2 = struct('fun', @(x) x.' * x, 'x0', [0 0]);
            problemInstance2 = Problem(s2);
            testCase.verifyEqual(problemInstance2.cub(ones(2, 1)), NaN(0, 1));

            s3 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'cub', @(x) x.' * x);
            problemInstance3 = Problem(s3);
            problemInstance3.cub([0; 0]);
            testCase.verifyEqual(problemInstance3.m_nonlinear_ub, 1);

            s4 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'cub', @(x) x.' * x, 'm_nonlinear_ub', 2);
            problemInstance4 = Problem(s4);
            testCase.verifyError(@() problemInstance4.cub([0; 0]), "MATLAB:Problem:cubx_m_nonlinear_ub_NotConsistent");
        end

        function testPrivateCEQ(testCase)
            % Test that CEQ works as expected.
            s1 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'ceq', @(x) x.' * x, 'm_nonlinear_eq', 1);
            problemInstance1 = Problem(s1);
            testCase.verifyError(@() problemInstance1.ceq(ones(2, 2, 3)), "MATLAB:Problem:InvalidInputForCEQ");
            testCase.verifyError(@() problemInstance1.ceq(ones(3, 1)), "MATLAB:Problem:WrongSizeInputForCEQ");
            testCase.verifyError(@() problemInstance1.ceq([0 0]), "MATLAB:Problem:InvalidOutputForCEQ");

            s2 = struct('fun', @(x) x.' * x, 'x0', [0 0]);
            problemInstance2 = Problem(s2);
            testCase.verifyEqual(problemInstance2.ceq(ones(2, 1)), NaN(0, 1));

            s3 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'ceq', @(x) x.' * x);
            problemInstance3 = Problem(s3);
            problemInstance3.ceq([0; 0]);
            testCase.verifyEqual(problemInstance3.m_nonlinear_eq, 1);

            s4 = struct('fun', @(x) x.' * x, 'x0', [0 0], 'ceq', @(x) x.' * x, 'm_nonlinear_eq', 2);
            problemInstance4 = Problem(s4);
            testCase.verifyError(@() problemInstance4.ceq([0; 0]), "MATLAB:Problem:ceqx_m_nonlinear_eq_NotConsistent");
        end


    end
end

