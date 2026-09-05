classdef TestConstraintValidity < matlab.unittest.TestCase
    % Undefined constraint values must survive initialization and histories.
    methods (Test)
        function testNaNNonlinearConstraints(testCase)
            for kind = {'cub', 'ceq'}
                options = struct('fun', @(x) 7, 'x0', 0, 'xl', -1, 'xu', 1);
                options.(kind{1}) = @(x) [NaN; -2];
                problem = Problem(options);
                testCase.verifyTrue(isnan(problem.maxcv(0)));
            end
            options = struct('fun', @(x) 7, 'x0', 0, ...
                'cub', @(x) NaN, 'ceq', @(x) 3);
            problem = Problem(options);
            testCase.verifyTrue(isnan(problem.maxcv(0)));
        end

        function testRowConstraintValues(testCase)
            for kind = {'cub', 'ceq'}
                for values = {[-1, 2], [-1, NaN]}
                    row = values{1};
                    problem = Problem(struct('fun', @(x) 7, 'x0', 0, ...
                        kind{1}, @(x) row));
                    expected = 2;
                    if any(isnan(row))
                        expected = NaN;
                    end
                    testCase.verifyEqual(problem.maxcv(0), expected);
                    for name = {'plain', 'quantized'}
                        featured = FeaturedProblem(problem, Feature(name{1}), 10, 0);
                        testCase.verifyEqual(featured.maxcv_init, expected);
                        testCase.verifyEqual(featured.maxcv(0), expected);
                    end
                end
            end
        end

        function testNaNLinearConstraint(testCase)
            problem = Problem(struct('fun', @(x) 7, 'x0', 0, ...
                'xl', -1, 'xu', 1, 'aub', 1, 'bub', NaN));
            testCase.verifyTrue(isnan(problem.maxcv(0)));
        end

        function testInitialValueAndHistory(testCase)
            problem = Problem(struct('fun', @(x) 7, 'x0', 0, 'cub', @(x) NaN));
            for name = {'plain', 'quantized'}
                featured = FeaturedProblem(problem, Feature(name{1}), 10, 0);
                testCase.verifyTrue(isnan(featured.maxcv_init));
                testCase.verifyTrue(isnan(featured.maxcv(0)));
                testCase.verifyEqual(featured.fun(0), 7);
                testCase.verifyTrue(all(isnan(featured.maxcv_hist)));
                testCase.verifyEqual(featured.n_eval_fun, 1);
            end
        end

        function testValidPointsSurviveFailure(testCase)
            problem = Problem(struct('fun', @(x) x(1)^2, 'x0', 0, ...
                'cub', @TestConstraintValidity.constraint));
            featured = FeaturedProblem(problem, Feature('plain'), 10, 0);
            featured.fun(0);
            featured.fun(1);
            featured.fun(0);
            testCase.verifyEqual(featured.maxcv_init, 0);
            testCase.verifyEqual(featured.maxcv_hist([1, 3]), [0, 0]);
            testCase.verifyTrue(isnan(featured.maxcv_hist(2)));
        end

        function testInfiniteBoundsAreSupported(testCase)
            problem = Problem(struct('fun', @(x) 0, 'x0', 0, ...
                'xl', -Inf, 'xu', Inf, 'cub', @(x) -Inf));
            testCase.verifyEqual(problem.maxcv(0), 0);
        end

        function testAffineFeasibleRegionBaseline(testCase)
            original = Problem(struct('fun', @(x) x'*x, 'x0', [0.5; 0.5; 0], ...
                'xl', [0; 0.5; -Inf], 'xu', [1; 0.5; Inf], ...
                'aub', [1, 0, 1], 'bub', 1.5, 'aeq', [0, 0, 1], 'beq', 0));
            matrices = {diag([-2, 1, 0.5]), [1, 1, 0; 0, 1, 1; 0, 0, 1], ...
                [0, 1, 0; 0, 0, 1; 1, 0, 0]};
            shift = [0.25; -0.5; 1];
            features = {Feature('permuted'), Feature('plain'), Feature('noisy'), ...
                Feature('perturbed_x0'), Feature('linearly_transformed'), ...
                Feature('linearly_transformed', struct('rotated', false))};
            for i = 1:numel(matrices)
                matrix = matrices{i};
                features{end+1} = Feature('custom', struct('mod_affine', ...
                    @(rng, p) deal(matrix, shift, inv(matrix))));
            end
            points = [0.5, -1, 0.5, 0.5; 0.5, 0.5, 1, 0.5; 0, 0, 0, 1];
            for i = 1:numel(features)
                feature = features{i};
                featured = FeaturedProblem(original, feature, 20, 7);
                [~, b, inverse] = feature.modifier_affine(7, original);
                solver_view = Problem(struct('fun', @(x) x'*x, 'x0', featured.x0, ...
                    'xl', featured.xl, 'xu', featured.xu, 'aub', featured.aub, ...
                    'bub', featured.bub, 'aeq', featured.aeq, 'beq', featured.beq));
                for j = 1:size(points, 2)
                    x = points(:, j);
                    y = inverse*(x - b);
                    testCase.verifyEqual(solver_view.maxcv(y) < 1e-12, original.maxcv(x) < 1e-12);
                end
            end
        end
    end

    methods (Static)
        function value = constraint(x)
            if x(1) == 0
                value = -1;
            else
                value = NaN;
            end
        end
    end
end
