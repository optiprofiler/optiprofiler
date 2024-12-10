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