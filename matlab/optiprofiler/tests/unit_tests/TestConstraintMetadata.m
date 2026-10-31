classdef TestConstraintMetadata < matlab.unittest.TestCase
    % Static constraint data: what construction refuses and what the legal
    % sentinels mean.
    properties (Constant)
        Sphere = @(x) sum(x.^2)
    end

    methods (Test)
        function malformedDataIsRefusedBeforeAnyCallback(testCase)
            cases = { ...
                'xl', struct('xl', NaN); ...
                'xl', struct('xl', Inf); ...
                'xl', struct('xl', [0; NaN]); ...
                'xu', struct('xu', NaN); ...
                'xu', struct('xu', -Inf); ...
                'xu', struct('xu', [1; -Inf]); ...
                'bub', struct('aub', 1, 'bub', -Inf); ...
                'bub', struct('aub', 1, 'bub', NaN); ...
                'bub', struct('aub', [1; 1], 'bub', [Inf; -Inf]); ...
                'beq', struct('aeq', 1, 'beq', NaN); ...
                'beq', struct('aeq', 1, 'beq', Inf); ...
                'beq', struct('aeq', 1, 'beq', -Inf); ...
                'aub', struct('aub', NaN, 'bub', 1); ...
                'aub', struct('aub', -Inf, 'bub', 1); ...
                'aub', struct('aub', NaN, 'bub', Inf); ...
                'aeq', struct('aeq', Inf, 'beq', 1); ...
                'aeq', struct('aeq', NaN, 'beq', 1)};
            for k = 1:size(cases, 1)
                options = cases{k, 2};
                options.fun = testCase.Sphere;
                options.x0 = 0;
                if isfield(options, 'xl')
                    options.x0 = zeros(numel(options.xl), 1);
                elseif isfield(options, 'xu')
                    options.x0 = zeros(numel(options.xu), 1);
                end
                options.cub = @TestConstraintMetadata.forbiddenCallback;
                testCase.verifyError(@() Problem(options), ...
                    "MATLAB:Problem:" + cases{k, 1} + "_InvalidValue", sprintf('case %d', k));
            end
        end

        function laterAssignmentIsRefusedWhenATrialIsBuilt(testCase)
            % A rotation poses only finite bounds, so the posed problem would look fine.
            for feature = {Feature('plain'), Feature('noisy'), Feature('linearly_transformed'), ...
                    Feature('permuted'), Feature({'noisy', 'permuted'}), Feature({'perturbed_x0', 'linearly_transformed'})}
                problem = Problem(struct('fun', testCase.Sphere, 'x0', [0; 0], 'xu', [1; 1]));
                problem.xl = [Inf; Inf];
                testCase.verifyError(@() FeaturedProblem(problem, feature{1}, 2, 7), "MATLAB:Problem:xl_InvalidValue");
                problem = Problem(struct('fun', testCase.Sphere, 'x0', [0; 0], 'aub', [1, 1], 'bub', 1));
                problem.bub = -Inf;
                testCase.verifyError(@() FeaturedProblem(problem, feature{1}, 2, 7), "MATLAB:Problem:bub_InvalidValue");
            end
        end

        function malformedProblemInAFileLoadsAsItWasAndIsRefusedWhenUsed(testCase)
            % A setter that refused it would drop the property on load, and
            % the problem would come back without its bound.
            problem = Problem(struct('fun', @(x) x^2, 'x0', 0));
            problem.xl = Inf;
            file = [tempname, '.mat'];
            cleanup = onCleanup(@() delete(file));
            save(file, 'problem');
            loaded = load(file);
            testCase.verifyEqual(loaded.problem.xl, Inf);
            testCase.verifyError(@() FeaturedProblem(loaded.problem, Feature('plain'), 2, 7), "MATLAB:Problem:xl_InvalidValue");
        end

        function largeSparseRowsAreCheckedWithoutFillingThem(testCase)
            n = 100000;
            aub = sparse(1:n, 1:n, 1, n, n);
            problem = Problem(struct('fun', testCase.Sphere, 'x0', zeros(n, 1), 'aub', aub, 'bub', ones(n, 1)));
            testCase.verifyEqual(problem.m_linear_ub, n);
            aub(5, 7) = NaN;
            testCase.verifyError(@() Problem(struct('fun', testCase.Sphere, 'x0', zeros(n, 1), ...
                'aub', aub, 'bub', ones(n, 1))), "MATLAB:Problem:aub_InvalidValue");
        end

        function messageNamesTheEntryAndTheRule(testCase)
            try
                Problem(struct('fun', testCase.Sphere, 'x0', 0, 'xl', Inf));
                testCase.verifyFail('The bound was accepted.');
            catch err
                testCase.verifySubstring(err.message, '`xl`');
                testCase.verifySubstring(err.message, 'Inf');
                testCase.verifySubstring(err.message, '-Inf for no bound');
            end
        end

        function oneSidedAndAbsentBoundsAreLegal(testCase)
            problem = Problem(struct('fun', testCase.Sphere, 'x0', [0; 0], ...
                'xl', [-Inf; 0], 'xu', [1; Inf]));
            testCase.verifyEqual(problem.ptype, 'b');
            testCase.verifyEqual(problem.mb, 2);
            testCase.verifyEqual(problem.maxcv([2; -3]), 3);
            problem = Problem(struct('fun', testCase.Sphere, 'x0', 0, 'xl', -Inf, 'xu', Inf));
            testCase.verifyEqual(problem.ptype, 'u');
        end

        function rowWithInfiniteRightHandSideImposesNoConstraint(testCase)
            only = Problem(struct('fun', testCase.Sphere, 'x0', 5, 'aub', 1, 'bub', Inf));
            testCase.verifyEqual(only.ptype, 'u');
            testCase.verifyEqual(only.m_linear_ub, 0);
            testCase.verifyEqual(only.maxcv(5), 0);
            mixed = Problem(struct('fun', testCase.Sphere, 'x0', 5, 'aub', [1; 2], 'bub', [Inf; 1]));
            testCase.verifyEqual(mixed.ptype, 'l');
            testCase.verifyEqual(mixed.m_linear_ub, 1);
            testCase.verifyEqual(size(mixed.aub), [2, 1]);
            testCase.verifyEqual(mixed.bub, [Inf; 1]);
            testCase.verifyEqual(mixed.maxcv(5), 9);
        end

        function inconsistentShapeNamesTheRequiredShape(testCase)
            try
                Problem(struct('fun', testCase.Sphere, 'x0', [0; 0], 'aub', 1, 'bub', Inf));
                testCase.verifyFail('The shape was accepted.');
            catch err
                testCase.verifyEqual(err.identifier, 'MATLAB:Problem:aub_m_linear_ub_n_NotConsistent');
                testCase.verifySubstring(err.message, '(1, 2)');
            end
        end

        function reversedFiniteBoundsAreMeasuredNotRefused(testCase)
            problem = Problem(struct('fun', testCase.Sphere, 'x0', 0, 'xl', 1, 'xu', 0));
            testCase.verifyEqual(problem.ptype, 'b');
            testCase.verifyEqual(problem.maxcv(0), 1);
            featured = FeaturedProblem(problem, Feature('plain'), 2, 7);
            testCase.verifyEqual(featured.maxcv_init, 1);
            featured.fun(0);
            testCase.verifyEqual(featured.maxcv_hist, 1);
        end

        function sentinelsKeepInitialAndHistoryTruthUnderFeatures(testCase)
            features = {Feature('plain'), Feature('noisy'), Feature('perturbed_x0'), ...
                Feature('permuted'), Feature('linearly_transformed'), ...
                Feature('linearly_transformed', struct('rotated', false)), ...
                Feature('quantized', struct('ground_truth', false)), ...
                Feature({'perturbed_x0', 'linearly_transformed'}), Feature({'permuted', 'noisy'})};
            for k = 1:numel(features)
                problem = testCase.withSentinels();
                featured = FeaturedProblem(problem, features{k}, 10, 7);
                testCase.verifyEqual(problem.m_linear_ub, 1);
                testCase.verifyTrue(isfinite(featured.maxcv_init), sprintf('feature %d', k));
                featured.fun(featured.x0);
                featured.fun(featured.x0 + 3);
                history = featured.maxcv_hist;
                testCase.verifyTrue(all(isfinite(history)), sprintf('feature %d', k));
                % The finite rows and bounds stay in force at the far point.
                testCase.verifyGreaterThan(history(2), 0, sprintf('feature %d', k));
                testCase.verifyEqual(history(1), featured.maxcv_init, 'RelTol', 1e-12, sprintf('feature %d', k));
            end
        end

        function sentinelRowSurvivesTransportAsARowWithoutForce(testCase)
            featured = FeaturedProblem(testCase.withSentinels(), ...
                Feature('linearly_transformed', struct('rotated', true)), 10, 7);
            testCase.verifyEqual(sum(featured.bub == Inf), 1);
            testCase.verifyEqual(size(featured.aub, 1), numel(featured.bub));
            % Two finite bounds become rows, as does the finite original row.
            testCase.verifyEqual(featured.m_linear_ub, 3);
            testCase.verifyEqual(featured.ptype, 'l');
        end

        function sentinelRowKeepsNoForceUnderACustomAffineMap(testCase)
            shift = [0.25; -0.5];
            for matrix = {diag([2, -1]), [1, 1; 0, 1]}
                stage = struct('name', 'custom', 'options', struct('mod_affine', ...
                    @(s, p) deal(matrix{1}, shift, inv(matrix{1}))));
                for stages = {{stage}, {stage, 'noisy'}}
                    problem = testCase.withSentinels();
                    featured = FeaturedProblem(problem, Feature(stages{1}), 10, 7);
                    testCase.verifyEqual(sum(featured.bub == Inf), 1);
                    testCase.verifyEqual(size(featured.aub, 1), numel(featured.bub));
                    testCase.verifyEqual(featured.maxcv_init, problem.maxcv(problem.x0), 'AbsTol', 1e-14);
                    point = featured.x0 + 3;
                    featured.fun(point);
                    testCase.verifyEqual(featured.maxcv_hist, problem.maxcv(matrix{1} * point + shift), 'RelTol', 1e-12);
                end
            end
        end

        function quantizedGroundTruthMeasuresOnlyRowsWithForce(testCase)
            problem = Problem(struct('fun', testCase.Sphere, 'x0', 5, 'aub', [1; 2], 'bub', [Inf; 1]));
            feature = Feature('quantized', struct('ground_truth', true));
            testCase.verifyEqual(FeaturedProblem(problem, feature, 10, 7).maxcv_init, 9);
            problem = Problem(struct('fun', testCase.Sphere, 'x0', 5, 'aub', 1, 'bub', Inf));
            featured = FeaturedProblem(problem, feature, 10, 7);
            testCase.verifyEqual(featured.ptype, 'u');
            testCase.verifyEqual(featured.maxcv_init, 0);
        end

        function customStageCannotPoseMalformedData(testCase)
            outputs = { ...
                'mod_bounds', @(s, p) deal(NaN, p.xu), 'xl'; ...
                'mod_bounds', @(s, p) deal(Inf, p.xu), 'xl'; ...
                'mod_bounds', @(s, p) deal(p.xl, NaN), 'xu'; ...
                'mod_bounds', @(s, p) deal(p.xl, -Inf), 'xu'; ...
                'mod_linear_ub', @(s, p) deal(1, -Inf), 'bub'; ...
                'mod_linear_ub', @(s, p) deal(1, NaN), 'bub'; ...
                'mod_linear_ub', @(s, p) deal(-Inf, 1), 'aub'; ...
                'mod_linear_eq', @(s, p) deal(1, Inf), 'beq'; ...
                'mod_linear_eq', @(s, p) deal(NaN, 0), 'aeq'};
            problem = Problem(struct('fun', testCase.Sphere, 'x0', 0));
            for k = 1:size(outputs, 1)
                custom = struct('name', 'custom', 'options', struct(outputs{k, 1}, outputs{k, 2}));
                testCase.verifyError(@() FeaturedProblem(problem, Feature({custom}), 10, 7), ...
                    "MATLAB:Problem:" + outputs{k, 3} + "_InvalidValue", sprintf('single %d', k));
                testCase.verifyError(@() FeaturedProblem(problem, Feature({custom, 'noisy'}), 10, 7), ...
                    'MATLAB:Feature:InvalidCustomOutput', sprintf('first %d', k));
                testCase.verifyError(@() FeaturedProblem(problem, Feature({'noisy', custom}), 10, 7), ...
                    'MATLAB:Feature:InvalidCustomOutput', sprintf('last %d', k));
            end
        end

        function customStageMayPoseTheLegalSentinels(testCase)
            custom = struct('name', 'custom', 'options', struct( ...
                'mod_bounds', @(s, p) deal(-Inf, Inf), 'mod_linear_ub', @(s, p) deal(1, Inf)));
            problem = Problem(struct('fun', testCase.Sphere, 'x0', 0));
            for stages = {{custom}, {custom, 'noisy'}}
                featured = FeaturedProblem(problem, Feature(stages{1}), 10, 7);
                testCase.verifyEqual(featured.m_linear_ub, 0);
                testCase.verifyEqual(featured.mb, 0);
                testCase.verifyEqual(featured.ptype, 'u');
            end
        end
    end

    methods
        function problem = withSentinels(testCase)
            problem = Problem(struct('fun', testCase.Sphere, 'x0', [0.5; 0.5], ...
                'xl', [0; -Inf], 'xu', [Inf; 1], 'aub', [1, 1; 1, 0], 'bub', [Inf; 1], ...
                'aeq', [1, -1], 'beq', 0));
        end
    end

    methods (Static)
        function value = forbiddenCallback(~)
            error('TestConstraintMetadata:CallbackCalled', 'A callback ran before the data were checked.');
            value = []; %#ok<UNRCH>
        end
    end
end
