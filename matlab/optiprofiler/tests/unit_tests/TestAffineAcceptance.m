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

        function illConditionedInverseCannotMoveAFeasibleStart(testCase)
            % This explicit inverse omits a unit even at factor=1. Its large
            % entries used to inflate both verification allowances enough to
            % accept images (0,1) or (0,10), instead of the feasible (1,2).
            for composed = [false, true]
                for dimension = [2, 3]
                    x0 = [1; 2]; lower = [0.5; 1.5]; upper = [1.5; 2.5];
                    if dimension == 3
                        % An unrelated huge coordinate must not inflate the
                        % cap on the two coordinates of the ill-conditioned
                        % block. A normwise x0 cap would hide their corruption.
                        x0 = [x0; 1e14]; lower = [lower; 0]; upper = [upper; 2e14]; %#ok<AGROW>
                    end
                    problem = TestAffineAcceptance.makeProblem(x0, ...
                        'fun', @(x) sum((x - x0).^2), 'xl', lower, 'xu', upper);
                    for delta = [3e-14, 1e-14]
                        A = [1, 1; 1, 1 + delta];
                        for factor = [1, 10]
                            inverse = (factor / delta) * [1, -1; -1, 1];
                            if dimension == 3
                                matrix = blkdiag(A, 1); inverse = blkdiag(inverse, 1);
                            else
                                matrix = A;
                            end
                            feature = TestAffineAcceptance.makeFeature( ...
                                matrix, zeros(dimension, 1), inverse, composed);
                            try
                                featured = FeaturedProblem(problem, feature, 10, 3);
                            catch err
                                % Fail-closed rejection is safe; a validated
                                % independent solve may also recover this case.
                                testCase.verifyTrue(any(strcmp(err.identifier, { ...
                                    'MATLAB:Feature:AffineTransformationNotInvertible', ...
                                    'MATLAB:Feature:AffineInitialPointNotRepresentable'})), ...
                                    err.message);
                                continue;
                            end
                            testCase.verifyEqual(matrix * featured.x0, x0, 'AbsTol', 1e-8);
                            testCase.verifyEqual(featured.maxcv_init, 0);
                        end
                    end
                end
            end
        end

        function ordinaryRotationsKeepRepresentableLargeDynamicRangeStarts(testCase)
            % STREG/STREGNE have this initial point. Mixing its large and
            % small coordinates naturally rounds the small image coordinates
            % by more than sqrt(eps). Built-in construction must retain its
            % established point. The same user-supplied map may trigger an
            % independent solve but must still remain representable.
            x0 = [-1.2; 1; 1e10; 1e10];
            problem = TestAffineAcceptance.makeProblem(x0);
            options = struct('rotated', true, 'condition_factor', 0);
            kernel = optiprofiler_internal.FeatureKernel('linearly_transformed', options);
            exceeded_component_cap = false;
            for seed = 0:9
                for composed = [false, true]
                    stage = struct('name', 'linearly_transformed', 'options', options);
                    if composed, feature = Feature({stage, 'noisy'});
                    else, feature = Feature({stage}); end
                    construction_seed = seed;
                    if composed
                        normalized_stage = feature.stages{1};
                        construction_seed = optiprofiler_internal.deriveFeatureStageSeed( ...
                            seed, normalized_stage.code, normalized_stage.occurrence, 3);
                    end
                    [A, ~, inverse] = kernel.modifier_affine(construction_seed, problem);
                    expected = inverse * x0;
                    amplification = 64 * numel(x0) * eps * norm(A, inf) * norm(inverse, inf);
                    testCase.assertLessThan(amplification, 1);
                    exceeded_component_cap = exceeded_component_cap || ...
                        any(abs(A * expected - x0) > sqrt(eps) * abs(x0));
                    featured = FeaturedProblem(problem, feature, 10, seed);
                    testCase.verifyEqual(featured.x0, expected);
                    testCase.verifyTrue(all(isfinite(featured.x0)));
                    allowance = 64 * numel(x0) * eps * (abs(A) * abs(expected) + abs(x0));
                    testCase.verifyTrue(all(abs(A * featured.x0 - x0) <= allowance));
                    custom_expected = expected;
                    if any(abs(A * expected - x0) > sqrt(eps) * abs(x0))
                        custom_expected = A \ x0;
                    end
                    custom_feature = TestAffineAcceptance.makeFeature(A, zeros(4, 1), inverse, composed);
                    custom_featured = FeaturedProblem(problem, custom_feature, 10, seed);
                    testCase.verifyEqual(custom_featured.x0, custom_expected);
                    custom_allowance = 64 * numel(x0) * eps * (abs(A) * abs(custom_expected) + abs(x0));
                    testCase.verifyTrue(all(abs(A * custom_featured.x0 - x0) <= custom_allowance));
                end
            end
            % Ensure this regression would expose an unconditional cap.
            testCase.verifyTrue(exceeded_component_cap);
        end

        function honestCoarseShearRetainsTheSolvedPointPolicy(testCase)
            % This exact affine pair can represent x0 only coarsely in double
            % precision. The local threshold must replace the inverse's
            % candidate by a solve, not become a universal final-error cap.
            A = [1, 1e15; 0, 1]; inverse = [1, -1e15; 0, 1];
            x0 = [1 / 3; 1 / 7];
            expected = A \ x0;
            testCase.assertTrue(any(abs(A * expected - x0) > sqrt(eps) * abs(x0)));
            problem = TestAffineAcceptance.makeProblem(x0);
            for composed = [false, true]
                feature = TestAffineAcceptance.makeFeature(A, [0; 0], inverse, composed);
                featured = FeaturedProblem(problem, feature, 10, 3);
                testCase.verifyEqual(featured.x0, expected);
                allowance = 64 * numel(x0) * eps * (abs(A) * abs(expected) + abs(x0));
                testCase.verifyTrue(all(abs(A * featured.x0 - x0) <= allowance));
            end
        end

        function partialInverseErrorCannotHideBelowAConditionTrigger(testCase)
            % A coarse norm-product trigger at one misses partial inverse
            % errors: these budgets are below one, yet the initial image used
            % to become (1,2.1) or to drift by 5e-4. Verification must detect
            % the actual point error, not only catastrophic inverse inflation.
            cases = [5e-13, 0.1; 1e-11, 0.001];
            for composed = [false, true]
                for dimension = [2, 3]
                    x0 = [1; 2];
                    if dimension == 3, x0 = [x0; 1e14]; end %#ok<AGROW>
                    problem = TestAffineAcceptance.makeProblem(x0);
                    for index = 1:size(cases, 1)
                        delta = cases(index, 1); distortion = cases(index, 2);
                        A = [1, 1; 1, 1 + delta];
                        actual_delta = A(2, 2) - 1;
                        inverse = [A(2, 2), -1; -1, 1] / actual_delta + ...
                            (distortion / delta) * [1, -1; -1, 1];
                        if dimension == 3
                            A = blkdiag(A, 1); inverse = blkdiag(inverse, 1);
                        end
                        feature = TestAffineAcceptance.makeFeature( ...
                            A, zeros(dimension, 1), inverse, composed);
                        try
                            featured = FeaturedProblem(problem, feature, 10, 3);
                        catch err
                            testCase.verifyTrue(any(strcmp(err.identifier, { ...
                                'MATLAB:Feature:AffineTransformationNotInvertible', ...
                                'MATLAB:Feature:AffineInitialPointNotRepresentable'})), ...
                                err.message);
                            continue;
                        end
                        testCase.verifyEqual(A * featured.x0, x0, 'AbsTol', 1e-8);
                    end
                end
            end
        end

        function inverseAccuracyTriggerHasNoAbsoluteUnitFloor(testCase)
            % A fixed absolute tolerance can hide the same wrong inverse by
            % expressing the problem in smaller units. Compare to an actual
            % independent solve, not an absolute test tolerance of our own.
            for scale = [1e-9, 1e-100]
                x0 = scale * [1; 2];
                problem = TestAffineAcceptance.makeProblem(x0);
                for delta = [3e-14, 1e-14]
                    A = [1, 1; 1, 1 + delta];
                    expected = A \ x0;
                    for factor = [1, 10]
                        inverse = (factor / delta) * [1, -1; -1, 1];
                        for composed = [false, true]
                            feature = TestAffineAcceptance.makeFeature(A, [0; 0], inverse, composed);
                            try
                                featured = FeaturedProblem(problem, feature, 10, 3);
                            catch err
                                testCase.verifyTrue(any(strcmp(err.identifier, { ...
                                    'MATLAB:Feature:AffineTransformationNotInvertible', ...
                                    'MATLAB:Feature:AffineInitialPointNotRepresentable'})), err.message);
                                continue;
                            end
                            testCase.verifyEqual(featured.x0, expected);
                        end
                    end
                end
            end
        end

        function zeroCoordinateTriggerCannotBorrowALargeCoordinateScale(testCase)
            % A nonzero floor inferred from another coordinate would hide
            % errors in the zero coordinates, especially after translation.
            % The unrelated last coordinate intentionally spans 1e14.
            x0 = [0; 0; 1e14];
            problem = TestAffineAcceptance.makeProblem(x0);
            for scale = [1, 1e-9, 1e-100]
                b = [scale; 2 * scale; 0];
                A = blkdiag([1, 1; 1, 1 + 3e-14], 1);
                inverse = blkdiag((10 / 3e-14) * [1, -1; -1, 1], 1);
                expected = A \ (x0 - b);
                for composed = [false, true]
                    feature = TestAffineAcceptance.makeFeature(A, b, inverse, composed);
                    try
                        featured = FeaturedProblem(problem, feature, 10, 3);
                    catch err
                        testCase.verifyTrue(any(strcmp(err.identifier, { ...
                            'MATLAB:Feature:AffineTransformationNotInvertible', ...
                            'MATLAB:Feature:AffineInitialPointNotRepresentable'})), err.message);
                        continue;
                    end
                    testCase.verifyEqual(featured.x0, expected);
                end
            end
        end

        function centeringCannotHideAnInaccurateInverseCandidate(testCase)
            % The original coordinates are large, but centering leaves a
            % small right-hand side. A tolerance scaled only by x0 used to
            % accept an error of eight while x0's original box was feasible.
            x0 = [1e14; 1e14];
            b = x0 + [1; 2];
            problem = TestAffineAcceptance.makeProblem(x0, ...
                'xl', x0 - 0.5, 'xu', x0 + 0.5);
            for delta = [3e-14, 1e-14]
                A = [1, 1; 1, 1 + delta];
                inverse = (10 / delta) * [1, -1; -1, 1];
                expected = A \ (x0 - b);
                for composed = [false, true]
                    feature = TestAffineAcceptance.makeFeature(A, b, inverse, composed);
                    featured = FeaturedProblem(problem, feature, 10, 3);
                    testCase.verifyEqual(featured.x0, expected);
                    testCase.verifyEqual(A * featured.x0 + b, x0);
                    testCase.verifyEqual(featured.maxcv_init, 0);
                end
            end
        end

        function largeTranslationCannotMaskASmallOriginalCoordinate(testCase)
            % The symmetric centering case needs the original-point scale:
            % a large right-hand side must not hide a nonzero image of x0=0.
            % These binary-exact shifts distinguish inverse multiplication
            % from independent division without trigonometric fixtures.
            x0 = 0;
            b = 2^46 + 2^-6;
            problem = TestAffineAcceptance.makeProblem(x0);
            exercised_solve = false;
            for A = [3, 10]
                inverse = 1 / A;
                candidate = inverse * (x0 - b);
                expected = A \ (x0 - b);
                image_error = abs(A * candidate + b - x0);
                testCase.assertLessThanOrEqual(image_error, sqrt(eps) * abs(x0 - b));
                exercised_solve = exercised_solve || ...
                    (image_error > 0 && candidate ~= expected);
                for composed = [false, true]
                    feature = TestAffineAcceptance.makeFeature(A, b, inverse, composed);
                    featured = FeaturedProblem(problem, feature, 10, 3);
                    if image_error > 0
                        testCase.verifyEqual(featured.x0, expected);
                    else
                        testCase.verifyEqual(featured.x0, candidate);
                    end
                end
            end
            % The regression must kill an RHS-only tolerance mutation, not
            % merely pass on a fixture whose two construction paths coincide.
            testCase.verifyTrue(exercised_solve);
        end

        function nonfiniteAllowanceAloneCannotCertifyAFiniteCandidate(testCase)
            % Both the candidate and its signed image stay finite. Only the
            % sum of absolute products overflows, so the allowance's own
            % finiteness check (not either point check) must reject it.
            A = [1, -1; 0, 1];
            x0 = [0; 1e308];
            problem = TestAffineAcceptance.makeProblem(x0);
            for composed = [false, true]
                for inverse_error = [0, 1e-9]
                    inverse = [1, 1 + inverse_error; 0, 1];
                    candidate = inverse * x0;
                    testCase.assertTrue(all(isfinite(candidate)));
                    testCase.assertTrue(all(isfinite(A * candidate)));
                    testCase.assertTrue(any(~isfinite(abs(A) * abs(candidate))));
                    if inverse_error == 0
                        % The forward error is exactly zero here. The new
                        % accuracy cap cannot hide removal of the separate
                        % finite-allowance premise.
                        testCase.assertEqual(A * candidate, x0);
                    end
                    feature = TestAffineAcceptance.makeFeature(A, [0; 0], inverse, composed);
                    testCase.verifyError(@() FeaturedProblem(problem, feature, 10, 3), ...
                        'MATLAB:Feature:AffineInitialPointNotRepresentable');
                end
            end
        end

        function recoverableIllConditionedMapRemainsUsable(testCase)
            % Do not replace point verification by a stricter global
            % condition cutoff: this ill-conditioned map has an accurate
            % image of its initial point and remains useful.
            x0 = [1; 2];
            problem = TestAffineAcceptance.makeProblem(x0);
            for composed = [false, true]
                for delta = [3e-14, 1e-14]
                    A = [1, 1; 1, 1 + delta];
                    stored_delta = A(2, 2) - 1;
                    inverse = [A(2, 2), -1; -1, 1] / stored_delta;
                    feature = TestAffineAcceptance.makeFeature(A, [0; 0], inverse, composed);
                    featured = FeaturedProblem(problem, feature, 10, 3);
                    testCase.verifyEqual(A * featured.x0, x0, 'AbsTol', 1e-8);
                end
            end
        end

        function equilibrationCannotSilentlyDiscardANonzeroEntry(testCase)
            % Dividing the first row by its maximum underflows the off-
            % diagonal term. The round trip must detect that information
            % loss, not use the altered diagonal matrix to certify A.
            A = [1e300, 1e-300; 0, 1];
            inverse = [1e-300, 0; 0, 1];
            problem = TestAffineAcceptance.makeProblem([0; 0]);
            for composed = [false, true]
                feature = TestAffineAcceptance.makeFeature(A, [0; 0], inverse, composed);
                try
                    FeaturedProblem(problem, feature, 10, 3);
                    testCase.verifyFail('An equilibration that loses data was accepted.');
                catch err
                    testCase.verifyEqual(err.identifier, ...
                        'MATLAB:Feature:AffineTransformationNotInvertible');
                    testCase.verifyTrue(contains(err.message, 'losing data'), err.message);
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
