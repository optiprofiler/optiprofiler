classdef TestFeatureExperimentPlan < matlab.unittest.TestCase
% Public experiment ownership: a reusable feature never owns repetitions.
    methods (Test)
        function canonicalFeatureUsesTopLevelCount(testCase)
            output = tempname;
            mkdir(output);
            testCase.addTeardown(@() rmdir(output, 's'));
            calls = [0, 0];
            feature = Feature('plain');
            problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [2; 1], 'name', 'PLAN_PROBE'));
            options = struct('feature', feature, 'n_runs', 3, ...
                'problem', problem, 'solver_isrand', [false, true], ...
                'solver_names', {{'deterministic', 'randomized'}}, ...
                'n_jobs', 1, 'max_eval_factor', 2, 'seed', 17, ...
                'score_only', true, 'draw_hist_plots', 'none', ...
                'silent', true, 'savepath', output);
            scores = benchmark({@(fun, x0) probe(1, fun, x0), ...
                @(fun, x0) probe(2, fun, x0)}, options);
            % These are actual user solver invocations, not planned run counts.
            testCase.verifyEqual(calls, [1, 3]);
            testCase.verifySize(scores, [2, 1]);
            testCase.verifyEmpty(feature.stages);
            testCase.verifyTrue(feature.is_identity);
            testCase.verifyEmpty(dir(fullfile(output, '**', 'data_for_loading.mat')));

            function x = probe(index, fun, x0)
                calls(index) = calls(index) + 1;
                fun(x0);
                x = 0.5 * x0;
                fun(x);
            end
        end

        function reusableRepeatedStagesKeepIndependentExperimentCounts(testCase)
            feature = Feature({struct('name', 'noisy', 'options', ...
                struct('noise_level', 0.01, 'distribution', 'gaussian')), ...
                struct('name', 'noisy', 'options', ...
                struct('noise_level', 0.1, 'distribution', 'uniform'))});
            before_stages = feature.stages;
            before_declared = feature.declared;
            [three_calls, three_traces] = runPublicProbe(feature, struct('n_runs', 3));
            [one_calls, one_traces] = runPublicProbe(feature, struct('n_runs', 1));
            testCase.verifyEqual(three_calls, [3, 3]);
            testCase.verifyEqual(one_calls, [1, 1]);
            % Fresh runtimes must restart at the same legacy run seed. Reusing
            % immutable configuration must not continue a previous run stream.
            for solver = 1:2
                testCase.verifyEqual(three_traces{solver}{1}, one_traces{solver}{1});
            end
            testCase.verifyEqual(feature.stages, before_stages);
            testCase.verifyEqual(feature.declared, before_declared);
            for stage = feature.stages
                testCase.verifyFalse(isfield(stage{1}.options, 'n_runs'));
            end
        end

        function reportedPlanMatchesObservedRuntimePolicy(testCase)
            % A plan declares implementation policy; the runtime independently
            % records what actually ran. Verify agreement at the public report
            % boundary for identity, effective single and composed execution.
            features = {Feature('plain+plain'), ...
                Feature('noisy', struct('noise_level', 0, 'noise_mode', 'deterministic')), ...
                Feature('noisy+truncated', struct('noise_level', 0, 'noise_mode', 'deterministic'))};
            expected = {'matlab-legacy-single-v1', 'matlab-legacy-single-v1', ...
                'matlab-composed-views-v1'};
            problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [2; 1], 'name', 'POLICY_PROBE'));
            output = tempname;
            mkdir(output);
            testCase.addTeardown(@() rmdir(output, 's'));
            for i = 1:numel(features)
                report_path = fullfile(output, sprintf('policy-%d.json', i));
                options = struct('feature', features{i}, 'n_runs', 1, 'problem', problem, ...
                    'solver_names', {{'first', 'second'}}, 'solver_isrand', [false, false], ...
                    'n_jobs', 1, 'max_eval_factor', 2, 'seed', 17, 'score_only', true, ...
                    'draw_hist_plots', 'none', 'silent', true, 'savepath', output, ...
                    'report_path', report_path);
                benchmark({@axisProbe, @axisProbe}, options);
                report = jsondecode(fileread(report_path));
                testCase.verifyEqual(report.status, 'completed');
                plan = report.configuration.effective.experiment.primary;
                testCase.assertTrue(isfield(plan, 'runtime_policy'), ...
                    'The report must retain the controller-owned implementation policy.');
                testCase.verifyEqual(plan.runtime_policy, expected{i});
                runs = report.problems(1).runs;
                testCase.assertNumElements(runs, 2);
                for run = 1:numel(runs)
                    item = runs(run);
                    if iscell(runs), item = runs{run}; end
                    testCase.verifyEqual(item.execution.kind, 'executed');
                    testCase.verifyEqual(item.runtime.runtime_policy, expected{i});
                    testCase.verifyEqual(item.runtime.runtime_policy, plan.runtime_policy);
                end
            end
        end

        function literalHintsAreNotStochasticPredicates(testCase)
            % Affine without rotation is deterministic but keeps literal hint
            % 5; custom is classified stochastic but its historical hint is 1.
            cases = {Feature('plain'), Feature('custom'), ...
                Feature('linearly_transformed', struct('rotated', false)), ...
                Feature('noisy', struct('noise_mode', 'deterministic')), ...
                Feature('noisy'), Feature('truncated'), ...
                Feature('truncated', struct('perturbed_trailing_digits', true)), ...
                Feature('perturbed_x0'), Feature('permuted'), Feature('random_nan'), ...
                Feature('unrelaxable_constraints'), Feature('nonquantifiable_constraints'), ...
                Feature('quantized'), Feature('linearly_transformed')};
            hints = [1, 1, 5, 1, 5, 1, 5, 5, 5, 5, 1, 1, 1, 5];
            actual = [1, 1, 1, 1, 5, 1, 5, 5, 5, 5, 1, 1, 1, 5];
            for i = 1:numel(cases)
                plan = optiprofiler_internal.resolveFeatureExperiment(cases{i}, ...
                    struct('solver_isrand', [false, false]), 'primary');
                testCase.verifyEqual(plan.n_runs, hints(i));
                testCase.verifyFalse(plan.request_present);
                testCase.verifyEqual(plan.origin, 'stage_hints');
                testCase.verifyEqual(runPublicProbe(cases{i}, struct()), [actual(i), actual(i)]);
            end
            % The solver rule has priority over stage hints, but an explicit
            % experiment count still wins. Assertions use real solver calls.
            testCase.verifyEqual(runPublicProbe(Feature('custom'), ...
                struct('solver_isrand', [false, true])), [5, 5]);
            testCase.verifyEqual(runPublicProbe(Feature('custom'), ...
                struct('solver_isrand', [false, true], 'n_runs', 2)), [2, 2]);
        end

        function invalidCountsFailBeforeOutputsOrSolvers(testCase)
            routes = {struct('feature_name', 'plain'), ...
                struct('feature_name', 'noisy+truncated'), ...
                struct('feature', struct('name', 'noisy')), ...
                struct('feature', {{'noisy', 'truncated'}}), ...
                struct('feature', Feature('plain'))};
            invalid = {[], NaN, Inf, 0, -1, 1.5, true, '3'};
            problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [2; 1]));
            for route = 1:numel(routes)
                for randomized = [false, true]
                    for value = 1:numel(invalid)
                        output = tempname;  % Must not be created by validation.
                        options = routes{route};
                        options.n_runs = invalid{value};
                        options.problem = problem;
                        options.solver_isrand = [false, randomized];
                        options.savepath = output;
                        options.score_only = true;
                        options.silent = true;
                        testCase.verifyError(@() benchmark({@forbiddenSolver, @forbiddenSolver}, options), ...
                            'MATLAB:checkValidityProfileOptions:n_runsNotValid');
                        testCase.verifyFalse(isfolder(output));
                    end
                end
            end
        end

        function entryAmbiguitiesFailBeforeOutputs(testCase)
            cases = {struct('feature', Feature('plain'), 'feature_name', 'plain'), ...
                struct('feature', 'noisy'), ...
                struct('feature', Feature('plain'), 'noise_level', 0.1), ...
                struct('feature', Feature('plain'), 'load', 'NO_ARCHIVE')};
            errors = {'MATLAB:benchmark:ConflictingFeatureInputs', ...
                'MATLAB:benchmark:InvalidFeatureInput', ...
                'MATLAB:benchmark:FeatureLocalOverride', 'MATLAB:benchmark:FeatureWithLoad'};
            for i = 1:numel(cases)
                output = tempname;
                options = cases{i};
                options.savepath = output;
                options.score_only = true;
                options.silent = true;
                testCase.verifyError(@() benchmark({@forbiddenSolver, @forbiddenSolver}, options), errors{i});
                testCase.verifyFalse(isfolder(output));
            end
        end

        function libraryReferenceHasItsOwnOneRunPlan(testCase)
            output = tempname;
            mkdir(output);
            original_path = path;
            original_registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            cleanup = onCleanup(@() restoreRegistry(original_path, original_registry, output)); %#ok<NASGU>
            registry = fullfile(output, 'registry.mat');
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', registry);
            fixture = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'featureplan');
            library = registerProblemLibrary(struct('name', 'featureplan', 'root', fixture, ...
                'select_function', 'plan_select', 'load_function', 'plan_load'));
            testCase.verifyEqual(library.name, 'featureplan');
            testCase.verifyTrue(isfile(registry));
            calls = [0, 0];
            feature = Feature('plain');
            options = struct('feature', feature, 'n_runs', 3, 'run_plain', true, ...
                'plibs', {{'featureplan'}}, 'ptype', 'u', 'mindim', 2, 'maxdim', 2, ...
                'solver_isrand', [false, true], 'solver_names', {{'deterministic', 'randomized'}}, ...
                'n_jobs', 1, 'max_eval_factor', 2, 'seed', 17, 'score_only', true, ...
                'draw_hist_plots', 'none', 'silent', true, 'savepath', output);
            scores = benchmark({@(fun, x0) probe(1, fun, x0), ...
                @(fun, x0) probe(2, fun, x0)}, options);
            % One selected problem: primary actual [1,3], plus independent
            % reference [1,1]. Inheriting primary count for randomized plain
            % would wrongly yield six calls for the second solver.
            testCase.verifyEqual(calls, [2, 4]);
            testCase.verifySize(scores, [2, 1]);
            testCase.verifyTrue(feature.is_identity);

            function x = probe(index, fun, x0)
                calls(index) = calls(index) + 1;
                fun(x0); x = 0.5 * x0; fun(x);
            end
        end

        function retainedAxesSeparateActualAndCopiedSlots(testCase)
            original = pwd;
            cleanup = onCleanup(@() cd(original)); %#ok<NASGU>
            cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'private'));
            get_defaults = @getDefaultProfileOptions;
            solve = @solveOneProblem;
            cd(original);
            feature = Feature('plain');
            problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [2; 1], 'name', 'AXIS_PROBE'));
            options = struct('n_runs', 3, 'solver_isrand', [false, true], ...
                'solver_names', {{'deterministic', 'randomized'}}, 'n_jobs', 1, ...
                'max_eval_factor', 2, 'seed', 17, 'score_only', true, 'silent', true);
            solvers = {@axisProbe, @axisProbe};
            options = get_defaults(solvers, feature, options);
            plan = optiprofiler_internal.resolveFeatureExperiment(feature, options, 'primary');
            result = solve(solvers, problem, feature, problem.name, length(problem.name), ...
                options, false, '', false, plan);
            testCase.verifySize(result.fun_history, [2, 3, 4]);
            testCase.verifyEqual(result.eval_report_metadata.real_n_runs, [1; 3]);
            testCase.verifyEqual(result.n_eval, 2 * ones(2, 3));
            testCase.verifyEqual(result.fun_out, 1.25 * ones(2, 3));
            testCase.verifyEqual(result.fun_history(1, 1, :), result.fun_history(1, 3, :));
            receipts = result.eval_report_metadata.runtime_receipts;
            testCase.verifySize(receipts, [2, 3]);
            testCase.verifyTrue(all(cellfun(@isempty, receipts(1, 2:3))));
            for i = 1:3
                testCase.verifyEqual(receipts{2, i}.run_seed, mod(23333 * 17 + 211 * i, 2^32));
            end
            reference = optiprofiler_internal.resolveFeatureExperiment(feature, options, 'plain_reference');
            result = solve(solvers, problem, feature, problem.name, length(problem.name), ...
                options, false, '', false, reference);
            testCase.verifySize(result.fun_history, [2, 1, 4]);
            testCase.verifyEqual(result.eval_report_metadata.real_n_runs, [1; 1]);
            testCase.verifyEqual(size(result.fun_inits, 1), 1);
        end
    end
end

function x = forbiddenSolver(varargin) %#ok<STOUT,INUSD>
    error('TestFeatureExperimentPlan:UnexpectedSolver', 'Validation must not invoke a solver.');
end

function x = axisProbe(fun, x0)
    fun(x0); x = 0.5 * x0; fun(x);
end

function restoreRegistry(original_path, original_registry, output)
    % The registry is wholly test-owned; avoid unregister's persisted startup
    % cleanup path. Restore the path/env and delete only our temporary folder.
    path(original_path);
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', original_registry);
    if isfolder(output), rmdir(output, 's'); end
end

function [calls, traces] = runPublicProbe(feature, extra)
    output = tempname;
    mkdir(output);
    cleanup = onCleanup(@() rmdir(output, 's')); %#ok<NASGU>
    calls = [0, 0];
    traces = {{}, {}};
    problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [2; 1], 'name', 'PLAN_PROBE'));
    options = struct('feature', feature, 'problem', problem, ...
        'solver_isrand', [false, false], 'solver_names', {{'first', 'second'}}, ...
        'n_jobs', 1, 'max_eval_factor', 2, 'seed', 17, 'score_only', true, ...
        'draw_hist_plots', 'none', 'silent', true, 'savepath', output);
    for key = fieldnames(extra)'
        options.(key{1}) = extra.(key{1});
    end
    benchmark({@(fun, x0) probe(1, fun, x0), @(fun, x0) probe(2, fun, x0)}, options);

    function x = probe(index, fun, x0)
        calls(index) = calls(index) + 1;
        first = fun(x0);
        x = 0.5 * x0;
        second = fun(x);
        traces{index}{end + 1} = [first, second];
    end
end
