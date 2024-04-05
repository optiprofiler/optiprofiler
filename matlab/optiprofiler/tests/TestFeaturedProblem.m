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
        
        function testSimple(testCase)
            nValues = [1, 10, 100];

            for n = nValues
                % Construct a simple problem.
                x0 = zeros(n, 1);
                pb_struct = struct('fun', @TestFeaturedProblem.rosen, 'x0', x0);
                problem = Problem(pb_struct);

                % Construct a featured problem.
                feature = Feature("PLAIN");
                featured_problem = FeaturedProblem(problem, feature, 10);
                

                % Check that the initialization of featured_problem is correct.
                testCase.verifyEqual(featured_problem.n_eval, 0);
                testCase.verifyEqual(size(featured_problem.fun_hist), [0 0]);
                testCase.verifyEqual(size(featured_problem.maxcv_hist), [0 0]);

                % Evaluate the objective function at x0.
                f = featured_problem.fun(x0);
                testCase.verifyEqual(featured_problem.n_eval, 1);
                testCase.verifyEqual(featured_problem.fun_hist, f);
                testCase.verifyEqual(featured_problem.maxcv_hist, 0);

                % Construct a featured problem with a different feature.
                feature = Feature("custom", 'modifier', @(x, f, seed) f + 1.0);
                featured_problem = FeaturedProblem(problem, feature, 10);

                % Evaluate the objective function at x0.
                f = featured_problem.fun(x0);
                testCase.verifyEqual(featured_problem.n_eval, 1);
                testCase.verifyEqual(featured_problem.fun_hist, f - 1);
                testCase.verifyEqual(featured_problem.maxcv_hist, 0);

                % % Construct a featured problem with randomized x0.
                feature = Feature("randomize_x0", 'distribution', @(rand_stream, n) ones(n, 1));
                featured_problem = FeaturedProblem(problem, feature, 10);

                % Verify that the x0 is changed.
                testCase.verifyEqual(featured_problem.x0, problem.x0 + 1);

            end
        end

        function testCatch(testCase)
            % Construct a nonlinearly constrained problem.
            x0 = zeros(2,1);
            pb_struct = struct('fun', @TestFeaturedProblem.rosen, 'x0', x0, 'cub', @TestFeaturedProblem.sum_cos, 'ceq', @TestFeaturedProblem.sum_sin);
            problem = Problem(pb_struct);

            % Construct a featured problem.
            feature = Feature("PLAIN");
            featured_problem = FeaturedProblem(problem, feature, 10);
        end

    end
end