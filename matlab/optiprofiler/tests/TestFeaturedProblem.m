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
        
        function testPlainProblem(testCase)
            % Test whether the featured problem (feature plain) works correctly.

            p = Problem(struct('fun', @TestFeaturedProblem.rosen, 'x0', ones(10, 1)));
            ft = Feature('plain');
            ftp = FeaturedProblem(p, ft, 500);
            testCase.verifyEqual(ftp.fun(ones(10, 1)), p.fun(ones(10, 1)));
        end

        function testPermutedProblem(testCase)
            % Test whether the featured problem (feature permuted) works correctly.

            p = Problem(struct('fun', @TestFeaturedProblem.rosen, 'x0', 0.5 * linspace(1, 10, 10)', 'xl', linspace(-10, -1, 10)', 'xu', linspace(1, 10, 10)', 'aub', diag(linspace(1, 10, 10)), 'bub', linspace(1, 10, 10), 'aeq', diag(linspace(1, 10, 10)), 'beq', linspace(1, 10, 10), 'cub', @TestFeaturedProblem.sum_cos, 'ceq', @TestFeaturedProblem.sum_sin));
            ft = Feature('permuted');
            ftp = FeaturedProblem(p, ft, 500, 1);
            testCase.verifyEqual(ftp.fun(ftp.x0), p.fun(p.x0));
            [~, reverse_permutation] = sort(ftp.permutation);
            testCase.verifyEqual(ftp.x0, p.x0(reverse_permutation));
            testCase.verifyEqual(ftp.xl, p.xl(reverse_permutation));
            testCase.verifyEqual(ftp.xu, p.xu(reverse_permutation));
            testCase.verifyEqual(ftp.aub, p.aub(:, reverse_permutation));
            testCase.verifyEqual(ftp.bub, p.bub);
            testCase.verifyEqual(ftp.aeq, p.aeq(:, reverse_permutation));
            testCase.verifyEqual(ftp.beq, p.beq);
            testCase.verifyEqual(ftp.fun(ftp.x0), p.fun(p.x0));
            testCase.verifyEqual(ftp.cub(ftp.x0), p.cub(p.x0));
            testCase.verifyEqual(ftp.ceq(ftp.x0), p.ceq(p.x0));
        end

        function testAffine_transformedProblem(testCase)
            % Test whether the featured problem (feature affine_transformed) works correctly.

            p = Problem(struct('fun', @TestFeaturedProblem.rosen, 'x0', ones(10, 1), 'xl', -ones(10, 1), 'xu', ones(10, 1), 'aub', diag(ones(10, 1)), 'bub', ones(10, 1), 'aeq', diag(ones(10, 1)), 'beq', ones(10, 1), 'cub', @TestFeaturedProblem.sum_cos, 'ceq', @TestFeaturedProblem.sum_sin));
            ft = Feature('affine_transformed', 'condition_number', 'dimension_dependent');
            ftp = FeaturedProblem(p, ft, 500, 1);
            scaler = ftp.scaler;
            rotation = ftp.rotation;
            testCase.verifyEqual(ftp.x0, (1 ./ scaler) .* (rotation' * p.x0));
            testCase.verifyEqual(ftp.cub(ones(10, 1)), p.cub(rotation * (scaler .* ones(10, 1))));
            testCase.verifyEqual(ftp.ceq(ones(10, 1)), p.ceq(rotation * (scaler .* ones(10, 1))));
        end

        function testReach_maxfun(testCase)
            % Test whether the featured problem (feature plain) works correctly.

            p = Problem(struct('fun', @TestFeaturedProblem.rosen, 'x0', ones(10, 1), 'xl', -ones(10, 1), 'xu', ones(10, 1), 'aub', diag(ones(10, 1)), 'bub', ones(10, 1), 'aeq', diag(ones(10, 1)), 'beq', ones(10, 1), 'cub', @TestFeaturedProblem.sum_cos, 'ceq', @TestFeaturedProblem.sum_sin));
            ft = Feature('plain');
            ftp = FeaturedProblem(p, ft, 1);
            ftp.fun(ones(10, 1));
            testCase.verifyEqual(ftp.fun(randn(10, 1)), ftp.fun(randn(10, 1)));
            testCase.verifyEqual(ftp.cub(ones(10, 1)), ftp.cub(ones(10, 1)));
            testCase.verifyEqual(ftp.ceq(ones(10, 1)), ftp.ceq(ones(10, 1)));
        end

        function testError(testCase)
            % Test whether all the errors are thrown correctly.

            p = Problem(struct('fun', @TestFeaturedProblem.rosen, 'x0', ones(10, 1)));
            ft = Feature('plain');
            testCase.verifyError(@() FeaturedProblem(p, p, 500), "MATLAB:FeaturedProblem:NotFeatureClass");
            testCase.verifyError(@() FeaturedProblem(p, ft, 0), "MATLAB:FeaturedProblem:max_evalNotPositiveInteger");
            testCase.verifyError(@() FeaturedProblem(p, ft, 500, -1), "MATLAB:FeaturedProblem:seedNotNonnegativeInteger");
            testCase.verifyError(@() FeaturedProblem(ft, ft, 500, 1), "MATLAB:FeaturedProblem:NotProblemClass");
        end

    end
end