classdef TestHistoricalNativeFeatures < matlab.unittest.TestCase
% Native Features written by e7d5341 itself, read by the repaired validators.
% A stored value that no longer describes a valid stage fails explicitly when
% loaded: MATLAB reports the loadobj error as a warning and returns an object
% without usable state, so it can never execute. Valid stored values load and
% are canonicalized (numeric 0/1 to logical, integer classes to double).
% See fixtures/feature-v2/e7-native-features.README.md for the provenance.

    properties (TestParameter)
        invalid = struct( ...
            'nan_perturbation', {{'nan_perturbation', 'MATLAB:Feature:perturbation_level_InvalidInput'}}, ...
            'negative_perturbation', {{'negative_perturbation', 'MATLAB:Feature:perturbation_level_InvalidInput'}}, ...
            'row_perturbation', {{'row_perturbation', 'MATLAB:Feature:perturbation_level_InvalidInput'}}, ...
            'inf_noise', {{'inf_noise', 'MATLAB:Feature:noise_level_NotPositive'}}, ...
            'nan_mesh', {{'nan_mesh', 'MATLAB:Feature:mesh_size_NotPositive'}});
    end

    methods (Test)
        function invalidHistoricalValuesFailExplicitly(testCase, invalid)
            [name, identifier] = invalid{:};
            testCase.verifyWarning(@() TestHistoricalNativeFeatures.loadVariable(name), identifier);
            state = warning('off', identifier);
            testCase.addTeardown(@() warning(state));
            loaded = TestHistoricalNativeFeatures.loadVariable(name);
            feature = loaded.(name);
            testCase.verifyError(@() feature.stages, ?MException);
        end

        function validHistoricalValuesLoadAndCanonicalize(testCase)
            testCase.verifyWarningFree(@() TestHistoricalNativeFeatures.loadVariable('valid_control'));
            control = TestHistoricalNativeFeatures.loadVariable('valid_control');
            testCase.verifyEqual(control.valid_control.stages{1}.options.noise_level, 1e-3);
            testCase.verifyWarningFree(@() TestHistoricalNativeFeatures.loadVariable('numeric_logical'));
            logical_case = TestHistoricalNativeFeatures.loadVariable('numeric_logical');
            testCase.verifyClass(logical_case.numeric_logical.stages{1}.options.rotated, 'logical');
            testCase.verifyTrue(logical_case.numeric_logical.stages{1}.options.rotated);
            testCase.verifyWarningFree(@() TestHistoricalNativeFeatures.loadVariable('int32_digits'));
            digits = TestHistoricalNativeFeatures.loadVariable('int32_digits');
            testCase.verifyClass(digits.int32_digits.stages{1}.options.significant_digits, 'double');
            testCase.verifyEqual(digits.int32_digits.stages{1}.options.significant_digits, 3);
        end
    end

    methods (Static)
        function file = fixture()
            file = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'feature-v2', 'e7-native-features.mat');
        end

        function loaded = loadVariable(name)
            loaded = load(TestHistoricalNativeFeatures.fixture(), name);
        end
    end
end
