classdef TestProblem < matlab.unittest.TestCase
    methods (Static)

        function f = rosen(x)
            f = sum(100*(x(2:end) - x(1:end-1).^2).^2 + (1 - x(1:end-1)).^2);
        end

        function g = rosen_grad(x)
            n = length(x);
            g = zeros(n, 1);
            g(1) = -400 * x(1) * (x(2) - x(1)^2) - 2 * (1 - x(1));
            for i = 2:n-1
                g(i) = 200 * (x(i) - x(i-1)^2) - 400 * x(i) * (x(i+1) - x(i)^2) - 2 * (1 - x(i));
            end
            g(n) = 200 * (x(n) - x(n-1)^2);
        end

        function h = rosen_hess(x)
            n = length(x);
            h = zeros(n, n);
            h(1, 1) = 1200 * x(1)^2 - 400 * x(2) + 2;
            h(1, 2) = -400 * x(1);
            for i = 2:n-1
                h(i, i-1) = -400 * x(i);
                h(i, i) = 202 + 1200 * x(i)^2 - 400 * x(i+1);
                h(i, i+1) = -400 * x(i);
            end
            h(n, n-1) = -400 * x(n);
            h(n, n) = 200;
        end

        function f = sum_cos(x)
            f = sum(cos(x));
        end

        function f = sum_sin(x)
            f = sum(sin(x));
        end

        function g = sum_cos_grad(x)
            g = -sin(x);
        end

        function g = sum_sin_grad(x)
            g = cos(x);
        end

        function h = sum_cos_hess(x)
            h = -diag(cos(x));
        end

        function h = sum_sin_hess(x)
            h = -diag(sin(x));
        end

        function f = errortest(x)
            f = 'x';
        end

    end

    methods (Test)

        function testNormalCase(testCase)
            % Test the normal case.

            option.fun = @TestProblem.rosen;
            option.x0 = -2 * ones(10, 1);
            option.xl = -ones(10, 1);
            option.xu = ones(10, 1);
            option.aub = diag(ones(10, 1));
            option.bub = ones(10, 1);
            option.aeq = diag(ones(10, 1));
            option.beq = ones(10, 1);
            option.cub = @TestProblem.sum_cos;
            option.ceq = @TestProblem.sum_sin;
            option.grad = @TestProblem.rosen_grad;
            option.hess = @TestProblem.rosen_hess;
            option.Jcub = @(x) TestProblem.sum_cos_grad(x)';
            option.Jceq = @(x) TestProblem.sum_sin_grad(x)';
            option.Hcub = @(x) {TestProblem.sum_cos_hess(x)};
            option.Hceq = @(x) {TestProblem.sum_sin_hess(x)};
            p = Problem(option);

            testCase.verifyEqual(p.fun(p.x0), TestProblem.rosen(p.x0));
            testCase.verifyEqual(p.x0, -2 * ones(10, 1));
            testCase.verifyEqual(p.xl, -ones(10, 1));
            testCase.verifyEqual(p.xu, ones(10, 1));
            testCase.verifyEqual(p.aub, diag(ones(10, 1)));
            testCase.verifyEqual(p.bub, ones(10, 1));
            testCase.verifyEqual(p.aeq, diag(ones(10, 1)));
            testCase.verifyEqual(p.beq, ones(10, 1));
            testCase.verifyEqual(p.cub(p.x0), TestProblem.sum_cos(p.x0));
            testCase.verifyEqual(p.ceq(p.x0), TestProblem.sum_sin(p.x0));
            testCase.verifyEqual(p.grad(p.x0), TestProblem.rosen_grad(p.x0));
            testCase.verifyEqual(p.hess(p.x0), TestProblem.rosen_hess(p.x0));
            testCase.verifyEqual(p.Jcub(p.x0), TestProblem.sum_cos_grad(p.x0)');
            testCase.verifyEqual(p.Jceq(p.x0), TestProblem.sum_sin_grad(p.x0)');
            testCase.verifyEqual(p.Hcub(p.x0), {TestProblem.sum_cos_hess(p.x0)});
            testCase.verifyEqual(p.Hceq(p.x0), {TestProblem.sum_sin_hess(p.x0)});
            testCase.verifyEqual(p.name, 'Unnamed Problem');
            testCase.verifyEqual(p.x_type, 'real');
            testCase.verifyEqual(p.n, 10);
            testCase.verifyEqual(p.m_linear_ub, 10);
            testCase.verifyEqual(p.m_linear_eq, 10);
            testCase.verifyEqual(p.m_nonlinear_ub, 1);
            testCase.verifyEqual(p.m_nonlinear_eq, 1);
            testCase.verifyEqual(p.p_type, 'nonlinearly constrained');

            p.project_x0;
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

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'grad', 1)), "MATLAB:Problem:grad_NotFunctionHandle")
            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'grad', @(x) 1));
            testCase.verifyError(@() p.grad('a'), "MATLAB:Problem:InvalidInputForGRAD")
            testCase.verifyError(@() p.grad(1), "MATLAB:Problem:WrongSizeInputForGRAD")
            testCase.verifyError(@() p.grad(ones(10, 1)), "MATLAB:Problem:InvalidOutputForGRAD")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'hess', 1)), "MATLAB:Problem:hess_NotFunctionHandle")
            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'hess', @(x) 1));
            testCase.verifyError(@() p.hess('a'), "MATLAB:Problem:InvalidInputForHESS")
            testCase.verifyError(@() p.hess(1), "MATLAB:Problem:WrongSizeInputForHESS")
            testCase.verifyError(@() p.hess(ones(10, 1)), "MATLAB:Problem:InvalidOutputForHESS")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'cub', @TestProblem.sum_cos, 'Jcub', 1)), "MATLAB:Problem:Jcub_NotFunctionHandle")
            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'cub', @TestProblem.sum_cos, 'Jcub', @(x) 1j));
            testCase.verifyError(@() p.Jcub('a'), "MATLAB:Problem:InvalidInputForJCUB")
            testCase.verifyError(@() p.Jcub(1), "MATLAB:Problem:WrongSizeInputForJCUB")
            testCase.verifyError(@() p.Jcub(ones(10, 1)), "MATLAB:Problem:InvalidOutputForJCUB")
            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'cub', @TestProblem.sum_cos, 'Jcub', @(x) 1));
            testCase.verifyError(@() p.Jcub(ones(10, 1)), "MATLAB:Problem:Jcubx_m_nonlinear_ub_n_NotConsistent")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'ceq', @TestProblem.sum_sin, 'Jceq', 1)), "MATLAB:Problem:Jceq_NotFunctionHandle")
            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'ceq', @TestProblem.sum_sin, 'Jceq', @(x) 1j));
            testCase.verifyError(@() p.Jceq('a'), "MATLAB:Problem:InvalidInputForJCEQ")
            testCase.verifyError(@() p.Jceq(1), "MATLAB:Problem:WrongSizeInputForJCEQ")
            testCase.verifyError(@() p.Jceq(ones(10, 1)), "MATLAB:Problem:InvalidOutputForJCEQ")
            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'ceq', @TestProblem.sum_sin, 'Jceq', @(x) 1));
            testCase.verifyError(@() p.Jceq(ones(10, 1)), "MATLAB:Problem:Jceqx_m_nonlinear_eq_n_NotConsistent")

            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'cub', @TestProblem.sum_cos, 'Hcub', 1)), "MATLAB:Problem:Hcub_NotFunctionHandle")
            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'cub', @TestProblem.sum_cos, 'Hcub', @(x) 1));
            testCase.verifyError(@() p.Hcub('a'), "MATLAB:Problem:InvalidInputForHCUB")
            testCase.verifyError(@() p.Hcub(1), "MATLAB:Problem:WrongSizeInputForHCUB")
            testCase.verifyError(@() p.Hcub(ones(10, 1)), "MATLAB:Problem:InvalidOutputForHCUB")
            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'cub', @TestProblem.sum_cos, 'Hcub', @(x) {1j}));
            testCase.verifyError(@() p.Hcub(ones(10, 1)), "MATLAB:Problem:InvalidOutputForHCUB")
            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'cub', @TestProblem.sum_cos, 'Hcub', @(x) {1}));
            testCase.verifyError(@() p.Hcub(ones(10, 1)), "MATLAB:Problem:Hcubx_m_nonlinear_ub_n_NotConsistent")
            
            testCase.verifyError(@() Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'ceq', @TestProblem.sum_sin, 'Hceq', 1)), "MATLAB:Problem:Hceq_NotFunctionHandle")
            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'ceq', @TestProblem.sum_sin, 'Hceq', @(x) 1));
            testCase.verifyError(@() p.Hceq('a'), "MATLAB:Problem:InvalidInputForHCEQ")
            testCase.verifyError(@() p.Hceq(1), "MATLAB:Problem:WrongSizeInputForHCEQ")
            testCase.verifyError(@() p.Hceq(ones(10, 1)), "MATLAB:Problem:InvalidOutputForHCEQ")
            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'ceq', @TestProblem.sum_sin, 'Hceq', @(x) {1j}));
            testCase.verifyError(@() p.Hceq(ones(10, 1)), "MATLAB:Problem:InvalidOutputForHCEQ")
            p = Problem(struct('fun', @TestProblem.rosen, 'x0', ones(10, 1), 'ceq', @TestProblem.sum_sin, 'Hceq', @(x) {1}));
            testCase.verifyError(@() p.Hceq(ones(10, 1)), "MATLAB:Problem:Hceqx_m_nonlinear_eq_n_NotConsistent")
            
        end
    end
end

