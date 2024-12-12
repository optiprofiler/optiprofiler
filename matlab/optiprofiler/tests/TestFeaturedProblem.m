classdef TestFeaturedProblem < matlab.unittest.TestCase
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

        function x0 = custom_mod_x0(rand_stream, p)
            x0 = p.x0 + rand_stream.rand(size(p.x0));
        end

        function [A, b, inv] = custom_mod_affine(rand_stream, p)
            n = size(p.x0, 1);
            A = diag(1 + rand_stream.rand(n, 1));
            b = rand_stream.rand(n, 1);
            inv = diag(1 ./ diag(A));
        end

        function [xl, xu] = custom_mod_bounds(rand_stream, p)
            n = size(p.x0, 1);
            xl = -rand_stream.rand(n, 1);
            xu = rand_stream.rand(n, 1);
        end

        function [aub, bub] = custom_mod_linear_ub(rand_stream, p)
            n = size(p.x0, 1);
            aub = rand_stream.rand(1, n);
            bub = rand_stream.rand(1);
        end

        function [aeq, beq] = custom_mod_linear_eq(rand_stream, p)
            n = size(p.x0, 1);
            aeq = rand_stream.rand(1, n);
            beq = rand_stream.rand(1);
        end

        function f = custom_mod_fun(x, rand_stream, p)
            f = p.fun(x) + rand_stream.rand(1);
        end

        function cub = custom_mod_cub(x, rand_stream, p)
            m = size(p.cub(x), 1);
            cub = p.cub(x) + rand_stream.rand(m, 1);
        end

        function ceq = custom_mod_ceq(x, rand_stream, p)
            m = size(p.ceq(x), 1);
            ceq = p.ceq(x) + rand_stream.rand(m, 1);
        end
    end

    methods (Test)
        
        function testNormalCase(testCase)
            % Test the normal case.

            p = Problem(struct('fun', @TestFeaturedProblem.rosen, 'x0', ones(10, 1), 'xl', -ones(10, 1), 'xu', ones(10, 1), 'aub', diag(ones(10, 1)), 'bub', ones(10, 1), 'aeq', diag(ones(10, 1)), 'beq', ones(10, 1), 'cub', @TestFeaturedProblem.sum_cos, 'ceq', @TestFeaturedProblem.sum_sin));
            ft = Feature('plain');
            fp = FeaturedProblem(p, ft, 500, 1);
            testCase.verifyEqual(fp.problem, p);
            testCase.verifyEqual(fp.feature, ft);
            testCase.verifyEqual(fp.max_eval, 500);
            testCase.verifyEqual(fp.seed, 1);
            testCase.verifyEqual(fp.fun_hist, []);
            testCase.verifyEqual(fp.maxcv_hist, []);
            testCase.verifyEqual(fp.last_cub, []);
            testCase.verifyEqual(fp.last_ceq, []);
            testCase.verifyEqual(fp.n_eval, 0);

            p = s_load('ALLINITA');
            options.mod_x0 = @TestFeaturedProblem.custom_mod_x0;
            options.mod_affine = @TestFeaturedProblem.custom_mod_affine;
            options.mod_bounds = @TestFeaturedProblem.custom_mod_bounds;
            options.mod_linear_ub = @TestFeaturedProblem.custom_mod_linear_ub;
            options.mod_linear_eq = @TestFeaturedProblem.custom_mod_linear_eq;
            options.mod_fun = @TestFeaturedProblem.custom_mod_fun;
            options.mod_cub = @TestFeaturedProblem.custom_mod_cub;
            options.mod_ceq = @TestFeaturedProblem.custom_mod_ceq;
            ft = Feature('custom', options);
            fp = FeaturedProblem(p, ft, 500, 1);
            fp.x0;
            fp.xl;
            fp.xu;
            fp.fun(fp.x0);
            fp.aub;
            fp.bub;
            fp.aeq;
            fp.beq;
            fp.cub(fp.x0);
            fp.ceq(fp.x0);
        end

        function testError(testCase)
            % Test whether the function throws errors as expected.

            testCase.verifyError(@() FeaturedProblem('problem', Feature('plain'), 500, 1), "MATLAB:FeaturedProblem:NotProblemClass");

            p = Problem(struct('fun', @TestFeaturedProblem.rosen, 'x0', ones(10, 1), 'xl', -ones(10, 1), 'xu', ones(10, 1), 'aub', diag(ones(10, 1)), 'bub', ones(10, 1), 'aeq', diag(ones(10, 1)), 'beq', ones(10, 1), 'cub', @TestFeaturedProblem.sum_cos, 'ceq', @TestFeaturedProblem.sum_sin));
            testCase.verifyError(@() FeaturedProblem(p, 'feature', 500, 1), "MATLAB:FeaturedProblem:NotFeatureClass");

            testCase.verifyError(@() FeaturedProblem(p, Feature('plain'), -1, 1), "MATLAB:FeaturedProblem:max_evalNotPositiveInteger");

            testCase.verifyError(@() FeaturedProblem(p, Feature('plain'), 500, -1), "MATLAB:FeaturedProblem:seedNotNonnegativeInteger");
        end
    end
end