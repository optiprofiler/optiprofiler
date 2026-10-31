classdef TestNativeBenchmarkReplay < matlab.unittest.TestCase
% Public trusted-MAT replay: callback identity is not a serializable equality.
    properties (TestParameter)
        callback_kind = {'named', 'anonymous', 'closure', 'custom'}
    end
    methods (Test)
        function nativeCallbackRoundTrip(testCase, callback_kind)
            output = tempname(artifactRoot()); mkdir(output);
            global OP_NATIVE_REPLAY_CALLBACK_CALLS
            OP_NATIVE_REPLAY_CALLBACK_CALLS = 0;
            testCase.addTeardown(@() clearCount());
            name = 'noisy'; key = 'noise_map';
            switch callback_kind
                case 'named', callback = @sum;
                case 'anonymous', callback = @(x) observeValue(x, 3);
                case 'closure', callback = capturedCallback(3);
                case 'custom'
                    name = 'custom'; key = 'mod_fun'; callback = capturedCallback(3);
            end
            feature = Feature({struct('name', name, 'options', struct(key, callback))});
            options_refined = nativeOptions(feature);
            file = fullfile(output, 'options.mat');
            save(file, 'options_refined', '-v7');
            before = fileBytes(file);
            [options, receipt] = loadBenchmarkOptions(file);
            testCase.verifyEqual(OP_NATIVE_REPLAY_CALLBACK_CALLS, 0, 'Import executed a user callback.');
            testCase.verifyEqual(fileBytes(file), before);
            testCase.verifyEqual(receipt.feature_specification_comparison, ...
                'ordinary_structure_and_values_callback_positions_only');
            testCase.verifyEqual(receipt.callback_comparison, 'opaque_native_handles_not_identity_or_semantics');
            restored = options.feature.stages{1}.options.(key);
            if strcmp(callback_kind, 'named'), expected = 3; else, expected = 6; end
            testCase.verifyEqual(restored([1, 2]), expected);
        end

        function duplicateIsInspectionOnlyNotCallbackAuthority(testCase)
            testCase.addTeardown(@() clearCount());
            feature = Feature({struct('name', 'noisy', 'options', struct('noise_map', capturedCallback(3)))});
            source = nativeOptions(feature);
            % A native callback cannot be certified semantically equal to its
            % duplicate. The canonical Feature remains the sole replay source.
            source.feature_specification{1}.options.noise_map = @(x) 999;
            [options, receipt] = loadBenchmarkOptions(source);
            testCase.verifyEqual(receipt.callback_comparison, 'opaque_native_handles_not_identity_or_semantics');
            callback = options.feature.stages{1}.options.noise_map;
            testCase.verifyEqual(callback([1, 2]), 6);
            source.feature_specification{1}.options.noise_level = 42;
            testCase.verifyError(@() loadBenchmarkOptions(source), 'OptiProfiler:InconsistentNativeFeature');
            source = nativeOptions(feature); source.feature_specification{1}.name = 'truncated';
            testCase.verifyError(@() loadBenchmarkOptions(source), 'OptiProfiler:InconsistentNativeFeature');
            source = nativeOptions(feature); source.feature_specification{1}.options.noise_map = 'chebyshev';
            testCase.verifyError(@() loadBenchmarkOptions(source), 'OptiProfiler:InconsistentNativeFeature');
            source = nativeOptions(feature); source.feature_specification = source.feature_specification';
            source.feature_specification{2} = source.feature_specification{1};
            testCase.verifyError(@() loadBenchmarkOptions(source), 'OptiProfiler:InconsistentNativeFeature');
        end

        function directProblemReplayRequiresTheProblem(testCase)
            % options_refined of a direct-Problem run does not store the
            % Problem; it must not replay as a default-library experiment.
            output = tempname(artifactRoot()); mkdir(output);
            testCase.addTeardown(@() rmdir(output, 's'));
            problem = Problem(struct('fun', @(x) sum((x - 1).^2), 'x0', [0; 0], 'name', 'replayquad'));
            solvers = {@replayStay, @replayStep};
            options = struct('problem', problem, 'feature_name', 'noisy', 'n_runs', 2, 'seed', 3, ...
                'savepath', output, 'benchmark_id', 'direct', 'silent', true, 'draw_hist_plots', 'none', ...
                'max_eval_factor', 5, 'max_tol_order', 2, 'n_jobs', 1);
            scores = benchmark(solvers, options);
            log = dir(fullfile(output, 'direct', '**', 'options_refined.mat'));
            file = fullfile(log(1).folder, 'options_refined.mat');
            S = load(file);
            testCase.verifyEqual(S.options_refined.problem_route, 'direct_problem');
            testCase.verifyFalse(any(isfield(S.options_refined, {'plibs', 'ptype', 'mindim', 'maxdim', 'problem_names'})));
            testCase.verifyError(@() loadBenchmarkOptions(file), 'OptiProfiler:DirectProblemNotSaved');
            S.options_refined.problem = problem;
            [replay, receipt] = loadBenchmarkOptions(S.options_refined);
            testCase.verifyEqual(receipt.problem_route, 'direct_problem');
            testCase.verifySameHandle(replay.problem, problem);
            replay.savepath = output; replay.benchmark_id = 'replayed';
            testCase.verifyEqual(benchmark(solvers, replay), scores);
            % The raw input of the same run records the Problem itself.
            U = load(fullfile(log(1).folder, 'options_user.mat'));
            [replay, receipt] = loadBenchmarkOptions(U.options_user);
            testCase.verifyEqual(receipt.problem_route, 'direct_problem');
            replay.savepath = output; replay.benchmark_id = 'user-replayed';
            testCase.verifyEqual(benchmark(solvers, replay), scores);
            % A library experiment records its route and keeps its selectors.
            source = nativeOptions(Feature('noisy'));
            source.problem_route = 'library_selection'; source.plibs = {'s2mpj'};
            [replay, receipt] = loadBenchmarkOptions(source);
            testCase.verifyEqual(receipt.problem_route, 'library_selection');
            testCase.verifyEqual(replay.plibs, {'s2mpj'});
            testCase.verifyFalse(isfield(replay, 'problem'));
            source.problem_route = 'elsewhere';
            testCase.verifyError(@() loadBenchmarkOptions(source), 'OptiProfiler:InvalidNativeOptions');
            source.problem_route = 'direct_problem'; source.problem = 3;
            testCase.verifyError(@() loadBenchmarkOptions(source), 'OptiProfiler:InvalidNativeProblem');
        end

        function malformedRetainedNativeFeatureHasActionableError(testCase)
            output = tempname(artifactRoot()); mkdir(output);
            options_refined = nativeOptions(Feature('noisy'));
            options_refined.feature = Feature.empty;
            file = fullfile(output, 'malformed-options.mat');
            save(file, 'options_refined', '-v7');
            testCase.verifyError(@() loadBenchmarkOptions(file), 'OptiProfiler:InvalidNativeFeature');
        end
    end
end

function source = nativeOptions(feature)
    stages = feature.stages; inspection = cell(1, numel(stages));
    for k = 1:numel(stages)
        inspection{k} = struct('name', stages{k}.name, 'options', stages{k}.options);
    end
    source = struct('schema', 'options_refined-v2', 'feature', feature, ...
        'feature_specification', {inspection}, 'n_runs', 3);
end

function x = replayStay(fun, x0)
    fun(x0);
    x = x0;
end

function x = replayStep(fun, x0)
    fun(x0);
    x = x0 + 0.5;
    fun(x);
end

function callback = capturedCallback(offset)
    callback = @(x) observeValue(x, offset);
end

function value = observeValue(x, offset)
    global OP_NATIVE_REPLAY_CALLBACK_CALLS
    if isempty(OP_NATIVE_REPLAY_CALLBACK_CALLS), OP_NATIVE_REPLAY_CALLBACK_CALLS = 0; end
    OP_NATIVE_REPLAY_CALLBACK_CALLS = OP_NATIVE_REPLAY_CALLBACK_CALLS + 1;
    value = sum(x) + offset;
end

function clearCount()
    clear global OP_NATIVE_REPLAY_CALLBACK_CALLS
end

function value = fileBytes(file)
    fid = fopen(file, 'rb'); cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    value = fread(fid, Inf, '*uint8');
end

function root = artifactRoot()
% CI may collect artifacts from OP_ARTIFACTS. Without it, use the system temp
% folder: tempname('') is an error in R2026a (MATLAB:tempname:MustBeString).
    root = getenv('OP_ARTIFACTS');
    if isempty(root), root = tempdir; end
end
