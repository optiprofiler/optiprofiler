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
