classdef TestAffineAcceptance < matlab.unittest.TestCase
% Independent acceptance regressions for the affine coordinate contract.
% These cases exercise public construction, so acceptance must mean that the
% solver receives an invertible, representable change of the original problem.

    methods (Static)
        function feature = makeFeature(A, b, inverse, composed)
            stage = struct('name', 'custom', 'options', struct( ...
                'mod_affine', @(~, ~) deal(A, b, inverse)));
            if composed
                feature = Feature({stage, 'noisy'});
            else
                feature = Feature({stage});
            end
        end

        function problem = makeProblem(x0, varargin)
            % Keep the objective finite even for the large-point controls;
            % these tests concern construction, not objective overflow.
            data = struct('fun', @(x) sum(x / 1e308), 'x0', x0);
            for index = 1:2:numel(varargin)
                data.(varargin{index}) = varargin{index + 1};
            end
            problem = Problem(data);
        end
    end

    methods (Test)
        function singularMapCannotHideBehindHugeInverseTerms(testCase)
            % Both products vanish, but a large fabricated inverse used to
            % make the componentwise rounding allowance larger than I.
            A = ones(2);
            inverse = 1e14 * [1, -1; -1, 1];
            problem = TestAffineAcceptance.makeProblem([0; 0]);
            for composed = [false, true]
                feature = TestAffineAcceptance.makeFeature(A, [0; 0], inverse, composed);
                testCase.verifyError(@() FeaturedProblem(problem, feature, 10, 3), ...
                    'MATLAB:Feature:AffineTransformationNotInvertible');
            end
        end

        function finiteInitialPointCannotBecomeInfinite(testCase)
            cases = {0.5, 0, 2; 2, -1e308, 0.5};
            problem = TestAffineAcceptance.makeProblem(1e308);
            for composed = [false, true]
                for index = 1:size(cases, 1)
                    feature = TestAffineAcceptance.makeFeature( ...
                        cases{index, 1}, cases{index, 2}, cases{index, 3}, composed);
                    testCase.verifyError(@() FeaturedProblem(problem, feature, 10, 3), ...
                        'MATLAB:Feature:AffineInitialPointNotRepresentable');
                end
            end
        end

        function largeValidInitialPointAndOrdinaryScalingRemainUsable(testCase)
            cases = {1, 1e308, 1e308; 0.5, 1, 2};
            for composed = [false, true]
                for index = 1:size(cases, 1)
                    scale = cases{index, 1};
                    problem = TestAffineAcceptance.makeProblem(cases{index, 2});
                    feature = TestAffineAcceptance.makeFeature(scale, 0, 1 / scale, composed);
                    featured = FeaturedProblem(problem, feature, 10, 3);
                    testCase.verifyEqual(featured.x0, cases{index, 3});
                    testCase.verifyTrue(all(isfinite(featured.x0)));
                end
            end
        end

        function overflowingAllowanceCannotAcceptAnInaccurateInitialPoint(testCase)
            % The supplied inverse passes the matrix-level tolerance, but its
            % point misses x0 by about 1e299. An infinite allowance used to
            % accept that miss; the equation solve must recover the identity.
            problem = TestAffineAcceptance.makeProblem(1e308);
            for composed = [false, true]
                feature = TestAffineAcceptance.makeFeature(1, 0, 1 + 1e-9, composed);
                featured = FeaturedProblem(problem, feature, 10, 3);
                testCase.verifyEqual(featured.x0, 1e308);
                testCase.verifyTrue(all(isfinite(featured.x0)));
            end
        end

        function unboundedForwardErrorCannotCertifyTheInitialPoint(testCase)
            % The signed products may be finite on some BLAS kernels, but
            % their absolute sum overflows. That cannot provide a finite
            % componentwise rounding bound, so construction must fail closed.
            A = [1, 1, -1; 0, 1, 0; 0, 0, 1];
            inverse = [1, -1, 1; 0, 1, 0; 0, 0, 1];
            problem = TestAffineAcceptance.makeProblem(1e308 * ones(3, 1));
            for composed = [false, true]
                feature = TestAffineAcceptance.makeFeature(A, zeros(3, 1), inverse, composed);
                testCase.verifyError(@() FeaturedProblem(problem, feature, 10, 3), ...
                    'MATLAB:Feature:AffineInitialPointNotRepresentable');
            end
        end

        function translationCannotCollapseStrictBounds(testCase)
            % x0=0 isolates the box collapse from initial-point cancellation.
            % The generic path must enforce the same restriction as the box
            % shortcut, even though it stores bounds as linear rows.
            problems = {TestAffineAcceptance.makeProblem(0, 'xl', 0, 'xu', 1), ...
                TestAffineAcceptance.makeProblem([0; 0], 'xl', [0; -Inf], 'xu', [1; Inf])};
            matrices = {1, [1, 1; 0, 1]};
            inverses = {1, [1, -1; 0, 1]};
            shifts = {1e16, [1e16; 0]};
            for composed = [false, true]
                for index = 1:numel(problems)
                    feature = TestAffineAcceptance.makeFeature( ...
                        matrices{index}, shifts{index}, inverses{index}, composed);
                    testCase.verifyError(@() FeaturedProblem(problems{index}, feature, 10, 3), ...
                        'MATLAB:Feature:AffineBoundsNotRepresentable');
                end
            end
        end

        function fixedAndRepresentableNarrowBoundsRemainUsable(testCase)
            % Exact fixed variables are deliberate. Positive intervals must
            % not be rejected merely because their width is very small.
            cases = { ...
                0, 0, 1e16, 0; ...
                1, 1 + eps, 1, 1; ...
                1e-200, 2e-200, 0, 1e-200; ...
                0, 4, 1e16, 2};
            for composed = [false, true]
                for index = 1:size(cases, 1)
                    lower = cases{index, 1}; upper = cases{index, 2};
                    shift = cases{index, 3}; initial = cases{index, 4};
                    problem = TestAffineAcceptance.makeProblem(initial, 'xl', lower, 'xu', upper);
                    feature = TestAffineAcceptance.makeFeature(1, shift, 1, composed);
                    featured = FeaturedProblem(problem, feature, 10, 3);
                    testCase.verifyEqual(featured.xl, lower - shift);
                    testCase.verifyEqual(featured.xu, upper - shift);
                    if lower < upper
                        testCase.verifyLessThan(featured.xl, featured.xu);
                    else
                        testCase.verifyEqual(featured.xl, featured.xu);
                    end
                end
            end
        end

        function scalingCannotCollapseStrictBounds(testCase)
            % Both products remain normal and finite. Their rounding alone
            % can nevertheless remove the width of a genuine interval.
            lower = 1.5 + 2 * eps;
            upper = 1.5 + 3 * eps;
            testCase.assertLessThan(lower, upper);
            testCase.assertEqual(0.75 * lower, 0.75 * upper);
            problem = TestAffineAcceptance.makeProblem(lower, 'xl', lower, 'xu', upper);
            for composed = [false, true]
                feature = TestAffineAcceptance.makeFeature(4 / 3, 0, 0.75, composed);
                testCase.verifyError(@() FeaturedProblem(problem, feature, 10, 3), ...
                    'MATLAB:Feature:AffineBoundsNotRepresentable');
            end
        end

        function replacedBoundsDoNotValidateTheDiscardedInterval(testCase)
            % mod_bounds replaces the original box; its old interval must
            % not be transformed or rejected after replacement.
            problems = {TestAffineAcceptance.makeProblem(0, 'xl', 0, 'xu', 1), ...
                TestAffineAcceptance.makeProblem([0; 0], 'xl', [0; -Inf], 'xu', [1; Inf])};
            matrices = {1, [1, 1; 0, 1]};
            inverses = {1, [1, -1; 0, 1]};
            shifts = {1e16, [1e16; 0]};
            for composed = [false, true]
                for index = 1:numel(problems)
                    A = matrices{index}; b = shifts{index}; inverse = inverses{index};
                    stage = struct('name', 'custom', 'options', struct( ...
                        'mod_affine', @(~, ~) deal(A, b, inverse), ...
                        'mod_bounds', @(~, p) deal(-ones(p.n, 1), ones(p.n, 1))));
                    if composed, feature = Feature({stage, 'noisy'});
                    else, feature = Feature({stage}); end
                    featured = FeaturedProblem(problems{index}, feature, 10, 3);
                    testCase.verifyEqual(featured.xl, -ones(featured.n, 1));
                    testCase.verifyEqual(featured.xu, ones(featured.n, 1));
                    testCase.verifyEmpty(featured.aub);
                    testCase.verifyEmpty(featured.aeq);
                end
            end
        end
    end
end
