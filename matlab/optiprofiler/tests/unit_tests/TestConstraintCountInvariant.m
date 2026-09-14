classdef TestConstraintCountInvariant < matlab.unittest.TestCase
% Constraint counters count recorded queries on every execution strategy.
% Literal oracles only. The identity and legacy-single rows failed on
% e7d5341, whose counters returned length() of m-by-k histories: the
% component count m before m queries, stale cached values once m reached
% MAX_EVAL, m*k for row results and 0 without constraints.

    properties (TestParameter)
        variant = {'m0_none', 'm0_empty_callback', 'm1', 'm2', 'm3', 'm2_row'};
        strategy = {'identity', 'legacy_single', 'composed'};
        budget = {2, 3};
    end

    methods (Test)
        function countsAreRecordedQueries(testCase, variant, strategy, budget)
            [p, m, cubTruth, ceqTruth] = TestConstraintCountInvariant.problem(variant);
            fp = FeaturedProblem(p, TestConstraintCountInvariant.feature(strategy), budget, 17);
            testCase.verifyEqual(fp.execution_strategy, strrep(strrep(strategy, '_', '-'), 'composed', 'composed-views'));
            for k = 1:(budget + 1)
                x = k - 1;
                c = fp.cub(x);
                e = fp.ceq(x);
                recorded = min(k, budget);
                testCase.verifyEqual([fp.n_eval_cub, fp.n_eval_ceq], [recorded, recorded], sprintf('counts after query %d', k));
                testCase.verifySize(fp.cub_hist, [m, recorded], sprintf('cub_hist after query %d', k));
                testCase.verifySize(fp.ceq_hist, [m, recorded], sprintf('ceq_hist after query %d', k));
                if k <= budget
                    testCase.verifyEqual(c(:), cubTruth(x), sprintf('fresh cub at query %d', k));
                    testCase.verifyEqual(e(:), ceqTruth(x), sprintf('fresh ceq at query %d', k));
                else
                    testCase.verifyEqual(c(:), cubTruth(budget - 1), 'cached cub after the budget');
                    testCase.verifyEqual(e(:), ceqTruth(budget - 1), 'cached ceq after the budget');
                end
            end
            % The hard stop is unchanged: 2*MAX_EVAL accepted public calls per channel.
            for k = (budget + 2):(2 * budget)
                fp.cub(0);
                fp.ceq(0);
            end
            testCase.verifyError(@() fp.cub(0), 'MATLAB:FeaturedProblem:cubExceedTerminationEval');
            testCase.verifyError(@() fp.ceq(0), 'MATLAB:FeaturedProblem:ceqExceedTerminationEval');
        end

        function unrecordedQueriesAreNotCounted(testCase, strategy)
            p = TestConstraintCountInvariant.problem('m2');
            fp = FeaturedProblem(p, TestConstraintCountInvariant.feature(strategy), 3, 17);
            fp.cub(0, false);
            fp.ceq(0, false);
            testCase.verifyEqual([fp.n_eval_cub, fp.n_eval_ceq], [0, 0]);
            fp.cub(1);
            fp.ceq(1);
            testCase.verifyEqual([fp.n_eval_cub, fp.n_eval_ceq], [1, 1]);
            testCase.verifySize(fp.cub_hist, [2, 1]);
            % Unrecorded calls still count toward the hard stop at 2*MAX_EVAL.
            for k = 1:4
                fp.cub(2);
            end
            testCase.verifyEqual(fp.n_eval_cub, 3);
            testCase.verifyError(@() fp.cub(2), 'MATLAB:FeaturedProblem:cubExceedTerminationEval');
        end

        function stochasticLegacyIndexIsPerQuery(testCase)
            % The served index of stochastic legacy constraint modifiers is the
            % number of earlier recorded queries, so repeated queries at one
            % point receive distinct streams. e7d5341 passed 0,3,3,3,4 for m=3.
            [p, m] = TestConstraintCountInvariant.problem('m3');
            testCase.assertEqual(m, 3);
            stages = {struct('name', 'noisy', 'options', struct('noise_level', 1e-3)), ...
                struct('name', 'custom', 'options', struct('mod_cub', @(x, s, q) q.cub(x) + s.rand(size(q.cub(x)))))};
            for i = 1:numel(stages)
                feature = Feature(stages(i));
                fp = FeaturedProblem(p, feature, 10, 17);
                kernel = optiprofiler_internal.FeatureKernel(feature.stages{1}.name, feature.stages{1}.options);
                previous = [];
                for k = 1:5
                    value = fp.cub(0.5);
                    testCase.verifyEqual(value, kernel.modifier_cub(0.5, 17, p, k - 1), ...
                        sprintf('%s query %d uses index %d', feature.name, k, k - 1));
                    if k > 1
                        testCase.verifyNotEqual(value, previous, sprintf('%s query %d repeats the previous value', feature.name, k));
                    end
                    previous = value;
                end
            end
        end
    end

    methods (Static)
        function [p, m, cubTruth, ceqTruth] = problem(variant)
            s = struct('name', 'count', 'x0', 0, 'fun', @(x) x(1)^2);
            switch variant
                case 'm0_none'
                    m = 0;
                case 'm0_empty_callback'
                    m = 0;
                    s.cub = @(x) [];
                    s.ceq = @(x) [];
                case 'm1'
                    m = 1;
                    s.cub = @(x) x(1);
                    s.ceq = @(x) x(1) - 2;
                case 'm2'
                    m = 2;
                    s.cub = @(x) [x(1); x(1) + 1];
                    s.ceq = @(x) [x(1) - 2; x(1) - 3];
                case 'm3'
                    m = 3;
                    s.cub = @(x) x(1) + (0:2)';
                    s.ceq = @(x) x(1) - (2:4)';
                case 'm2_row'
                    m = 2;
                    s.cub = @(x) [x(1), x(1) + 1];
                    s.ceq = @(x) [x(1) - 2, x(1) - 3];
            end
            p = Problem(s);
            cubTruth = @(x) reshape(x + (0:m-1)', [], 1);
            ceqTruth = @(x) reshape(x - (2:m+1)', [], 1);
        end

        function f = feature(strategy)
            zero = struct('name', 'noisy', 'options', struct('noise_level', 0));
            switch strategy
                case 'identity'
                    f = Feature('plain');
                case 'legacy_single'
                    f = Feature({zero});
                case 'composed'
                    f = Feature({zero, zero});
            end
        end
    end
end
