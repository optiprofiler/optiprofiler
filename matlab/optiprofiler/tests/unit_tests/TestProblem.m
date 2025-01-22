classdef TestProblem < matlab.unittest.TestCase
    methods (Static)

        function f = rosen(x)
            f = sum(100*(x(2:end) - x(1:end-1).^2).^2 + (1 - x(1:end-1)).^2);
        end

        function f = sum_cos(x)
            f = sum(cos(x));
        end

        function f = sum_sin(x)
            f = sum(sin(x));
        end

        function f = errortest(x)
            f = 'x';
        end

    end

    methods (Test)

        function testNormalCase(testCase)
            % Test the normal case.

            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'xl', -ones(10, 1), 'xu', ones(10, 1), 'aub', diag(ones(10, 1)), 'bub', ones(10, 1), 'aeq', diag(ones(10, 1)), 'beq', ones(10, 1), 'cub', @TestProblem.sum_cos, 'ceq', @TestProblem.sum_sin));
            testCase.verifyEqual(p.fun(p.x0), TestProblem.rosen(p.x0));
            testCase.verifyEqual(p.x0, ones(10, 1));
            testCase.verifyEqual(p.xl, -ones(10, 1));
            testCase.verifyEqual(p.xu, ones(10, 1));
            testCase.verifyEqual(p.aub, diag(ones(10, 1)));
            testCase.verifyEqual(p.bub, ones(10, 1));
            testCase.verifyEqual(p.aeq, diag(ones(10, 1)));
            testCase.verifyEqual(p.beq, ones(10, 1));
            testCase.verifyEqual(p.cub(p.x0), TestProblem.sum_cos(p.x0));
            testCase.verifyEqual(p.ceq(p.x0), TestProblem.sum_sin(p.x0));
            testCase.verifyEqual(p.name, 'Unnamed Problem');
            testCase.verifyEqual(p.x_type, 'real');
            testCase.verifyEqual(p.n, 10);
            testCase.verifyEqual(p.m_linear_ub, 10);
            testCase.verifyEqual(p.m_linear_eq, 10);
            testCase.verifyEqual(p.m_nonlinear_ub, 1);
            testCase.verifyEqual(p.m_nonlinear_eq, 1);
            testCase.verifyEqual(p.p_type, 'nonlinearly constrained');
        end

        function testErrors(testCase)
            % Test whether the function throws errors as expected.

            testCase.verifyError(@() Problem(), "MATLAB:Problem:MissingArguments")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen)), "MATLAB:Problem:MissingFields")

            testCase.verifyError(@() Problem('fun'), "MATLAB:Problem:NotStruct")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'x_type', 'mixed')), "MATLAB:Problem:x_type_NotValid")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', 'a')), "MATLAB:Problem:x0_NotRealVector")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'xl', 'a')), "MATLAB:Problem:xl_x0_NotConsistent")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'xl', -ones(10, 1), 'xu', 'a')), "MATLAB:Problem:xu_x0_NotConsistent")
            
            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'xl', -ones(10, 1), 'xu', ones(10, 1), 'aub', ones(10, 10), 'bub', 'a')), "MATLAB:Problem:bub_NotRealVector")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'xl', -ones(10, 1), 'xu', ones(10, 1), 'aub', 'a')), "MATLAB:Problem:aub_m_linear_ub_n_NotConsistent")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'xl', -ones(10, 1), 'xu', ones(10, 1), 'aeq', ones(10, 10), 'beq', 'a')), "MATLAB:Problem:beq_NotRealVector")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'xl', -ones(10, 1), 'xu', ones(10, 1), 'aeq', 'a')), "MATLAB:Problem:aeq_m_linear_eq_n_NotConsistent")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'xl', -ones(10, 1), 'xu', ones(10, 1), 'm_nonlinear_ub', 1)), "MATLAB:Problem:m_nonlinear_ub_cub_NotConsistent")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'xl', -ones(10, 1), 'xu', ones(10, 1), 'm_nonlinear_eq', 1)), "MATLAB:Problem:m_nonlinear_eq_ceq_NotConsistent")

            testCase.verifyError(@() Problem(struct('fun', 1, 'x0', ones(10, 1))), "MATLAB:Problem:fun_NotFunctionHandle")

            testCase.verifyError(@() Problem(struct('fun', @(x,y) x.^2 + y.^2, 'x0', ones(10, 1))), "MATLAB:Problem:fun_NotOneArguementFun")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', 1, 'cub', 1)), "MATLAB:Problem:cub_NotFunctionHandle")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', 1, 'ceq', 1)), "MATLAB:Problem:ceq_NotFunctionHandle")

            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1)));
            testCase.verifyError(@() p.maxcv('a'), "MATLAB:Problem:InvalidInputForMaxCV")

            testCase.verifyError(@() p.maxcv(1), "MATLAB:Problem:WrongSizeInputForMaxCV")

            testCase.verifyError(@() p.fun('a'), "MATLAB:Problem:InvalidInputForFUN")

            testCase.verifyError(@() p.fun(1), "MATLAB:Problem:WrongSizeInputForFUN")

            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'cub', @TestProblem.sum_cos, 'ceq', @TestProblem.sum_sin));
            testCase.verifyError(@() p.cub('a'), "MATLAB:Problem:InvalidInputForCUB")

            testCase.verifyError(@() p.cub(1), "MATLAB:Problem:WrongSizeInputForCUB")

            testCase.verifyError(@() p.ceq('a'), "MATLAB:Problem:InvalidInputForCEQ")

            testCase.verifyError(@() p.ceq(1), "MATLAB:Problem:WrongSizeInputForCEQ")

            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'cub', @TestProblem.errortest, 'ceq', @TestProblem.errortest));
            testCase.verifyError(@() p.cub(ones(10, 1)), "MATLAB:Problem:InvalidOutputForCUB")

            testCase.verifyError(@() p.ceq(ones(10, 1)), "MATLAB:Problem:InvalidOutputForCEQ")
        end
    end
end

