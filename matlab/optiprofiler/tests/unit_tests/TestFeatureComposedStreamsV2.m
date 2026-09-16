classdef TestFeatureComposedStreamsV2 < matlab.unittest.TestCase
% Composed-view payload streams (seed_policy matlab-stage-horner32-v2).
%
% The per-query random stream of a composed stage must depend on the whole
% observed payload (values, point, served index) even when a coordinate, a
% value or the index is zero: version 1 handed the payload to the legacy
% product mixer, which collapsed to the bare seed whenever any element was
% zero. The identity/single-feature strategies keep the legacy mixer and its
% established streams. The fold is a finite 32-bit hash: these tests pin
% distinct draws for specific payloads, not independence or absence of
% collisions in general.
    properties (Access = private)
        Registry
    end
    methods (TestMethodSetup)
        function isolateRegistry(testCase)
            folder = fileparts(mfilename('fullpath'));
            registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            original_path = path;
            testCase.addTeardown(@() path(original_path));
            testCase.addTeardown(@() setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', registry));
            testCase.Registry = [tempname, '.mat'];
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', testCase.Registry);
            addpath(fullfile(folder, '..', 'fixtures', 'composedstreams'));
        end
    end
    methods (Test)
        function zeroCoordinatesValuesAndIndexKeepPointDependence(testCase)
            feature = TestFeatureComposedStreamsV2.composedNoise();
            zero_valued = Problem(struct('fun', @(x) 0, 'x0', [0; 0; 0], ...
                'cub', @(x) [0; 0], 'ceq', @(x) 0, 'name', 'zeros'));
            points = {[0; 0; 0], [1; 0; 0], [0; 1; 0], [0; 0; 1], [0.5; 0; 0], [0; 0; 0]};
            draws = zeros(1, numel(points));
            cub_draws = zeros(2, numel(points));
            ceq_draws = zeros(1, numel(points));
            for k = 1:numel(points)
                % A fresh runtime per point: every query is the first one served
                % (index 0), so only the point differs between the draws.
                fp = FeaturedProblem(zero_valued, feature, 10, 17);
                draws(k) = fp.fun(points{k});
                cub_draws(:, k) = fp.cub(points{k});
                ceq_draws(k) = fp.ceq(points{k});
            end
            testCase.verifyEqual(draws(6), draws(1), 'Equal seed, payload and index must give equal draws.');
            testCase.verifyEqual(cub_draws(:, 6), cub_draws(:, 1));
            testCase.verifyEqual(ceq_draws(6), ceq_draws(1));
            testCase.verifyEqual(numel(unique(draws(1:5))), 5, ...
                'Objective draws collapsed across points with zero coordinates and a zero value.');
            testCase.verifyEqual(size(unique(cub_draws(:, 1:5)', 'rows'), 1), 5, ...
                'Inequality-constraint draws collapsed across zero-coordinate points.');
            testCase.verifyEqual(numel(unique(ceq_draws(1:5))), 5, ...
                'Equality-constraint draws collapsed across zero-coordinate points.');
            % The draws are exactly the documented policy arithmetic: stage
            % seed (identity fold) then the word fold over (values, point, index).
            x = points{2};
            stage_seed = optiprofiler_internal.deriveFeatureStageSeed(17, 2, 0, 0);
            stream = optiprofiler_internal.FeatureKernel.horner32_payload_rng(stage_seed, 0, x(1), x(2), x(3), 0);
            testCase.verifyEqual(draws(2), randn(stream), 'AbsTol', 0);
            % Repeated queries at one point advance the served index and draw again.
            fp = FeaturedProblem(zero_valued, feature, 10, 17);
            first = fp.fun([0; 0; 0]); second = fp.fun([0; 0; 0]);
            testCase.verifyNotEqual(first, second, 'A repeated query must not reuse the first draw.');
            testCase.verifyEqual(first, draws(1));
            receipt = fp.runtimeReceipt();
            testCase.verifyEqual(receipt.seed_policy, 'matlab-stage-horner32-v2');
            testCase.verifyEqual(receipt.stages{1}.served.fun, 2);
        end

        function randomNanIsNotAllOrNothingAcrossZeroCoordinatePoints(testCase)
            feature = Feature({struct('name', 'random_nan', 'options', struct('nan_rate', 0.5)), ...
                struct('name', 'noisy', 'options', struct('noise_level', 0, 'noise_type', 'absolute'))});
            problem = Problem(struct('fun', @(x) 0, 'x0', [0; 0; 0], 'name', 'zeros'));
            points = {[0; 0; 0], [1; 0; 0], [0; 1; 0], [0; 0; 1]};
            mixed = 0;
            for seed = 1:20
                flags = false(1, numel(points));
                for k = 1:numel(points)
                    fp = FeaturedProblem(problem, feature, 10, seed);
                    flags(k) = isnan(fp.fun(points{k}));
                end
                mixed = mixed + (any(flags) && ~all(flags));
            end
            % Version 1 gave 0 mixed seeds here (the draw depended on the seed
            % only). A pinned count, not a statistical claim.
            testCase.verifyGreaterThan(mixed, 0, 'random_nan was all-or-nothing on every seed.');
        end

        function singleFeatureStreamsKeepTheLegacyMixer(testCase)
            % The identity/single strategies are unchanged: their draws equal
            % the legacy formula, including its documented collapse for zero
            % payload elements (retained for reproducibility of old archives).
            feature = Feature('noisy', struct('noise_type', 'absolute', 'noise_level', 1));
            problem = Problem(struct('fun', @(x) 0, 'x0', [0; 0; 0], 'name', 'zeros'));
            points = {[0; 0; 0], [1; 0; 0], [0; 1; 0], [0.5; 0; 0]};
            draws = zeros(1, numel(points));
            for k = 1:numel(points)
                fp = FeaturedProblem(problem, feature, 10, 17);
                testCase.verifyEqual(fp.seed_policy, 'legacy-run-seed');
                draws(k) = fp.fun(points{k});
                x = points{k};
                legacy = randn(Feature.default_rng(17, 0, x(1), x(2), x(3), 0));
                testCase.verifyEqual(draws(k), legacy, 'AbsTol', 0);
            end
            testCase.verifyEqual(numel(unique(draws)), 1, 'The legacy single-feature stream must be unchanged.');
        end

        function composedStreamsAreIdenticalOnActualWorkers(testCase)
            testCase.assumeFalse(isempty(ver('parallel')), 'Parallel Computing Toolbox is not installed.');
            fixture_root = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'composedstreams');
            registerProblemLibrary(struct('name', 'composed_streams', 'root', fixture_root, ...
                'select_function', 'composed_fixture_select', 'load_function', 'composed_fixture_load'));
            output = tempname; mkdir(output);
            testCase.addTeardown(@() rmdir(output, 's'));
            options = struct('plibs', {{'composed_streams'}}, 'feature_name', 'noisy+truncated', ...
                'noise_level', 0.5, 'significant_digits', 4, 'n_runs', 2, 'seed', 3, 'score_only', true, ...
                'silent', true, 'max_eval_factor', 4, 'max_tol_order', 1, 'savepath', output, 'n_jobs', 1);
            [s1, p1, c1] = benchmark({@composedStreamsStay, @composedStreamsZero}, options);
            options.n_jobs = 2;
            [s2, p2, c2] = benchmark({@composedStreamsStayOnWorker, @composedStreamsZeroOnWorker}, options);
            testCase.verifyTrue(isequaln(s1, s2) && isequaln(p1, p2) && isequaln(c1, c2), ...
                'Composed streams differ between the controller and actual workers.');
        end
    end
    methods (Static)
        function feature = composedNoise()
            % Absolute noise of level one exposes the standardized draw of the
            % first stage; the second noisy stage (level zero) keeps two
            % effective stages without changing the observed value.
            feature = Feature({struct('name', 'noisy', 'options', struct('noise_type', 'absolute', 'noise_level', 1)), ...
                struct('name', 'noisy', 'options', struct('noise_type', 'absolute', 'noise_level', 0))});
        end
    end
end

function x = composedStreamsStay(fun, x0)
    fun(x0); fun(x0 + 1); x = x0;
end

function x = composedStreamsZero(fun, x0)
    x = zeros(size(x0)); fun(x); fun(x0);
end

function x = composedStreamsStayOnWorker(fun, x0)
    assert(~isempty(getCurrentTask()), 'Expected actual MATLAB worker execution.');
    x = composedStreamsStay(fun, x0);
end

function x = composedStreamsZeroOnWorker(fun, x0)
    assert(~isempty(getCurrentTask()), 'Expected actual MATLAB worker execution.');
    x = composedStreamsZero(fun, x0);
end
