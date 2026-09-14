classdef TestBenchmarkFeatureNameErrors < matlab.unittest.TestCase
% benchmark keeps MATLAB:benchmark:feature_nameNotValid for an unknown or
% malformed feature_name, with the Feature error as its cause. e7d5341 raised
% the Feature identifier instead, because the 2.0 redesign moved name checks
% into the Feature constructor (commit 77dfbea); ac4a67e raised the benchmark
% identifier. The structured feature route and option errors keep their
% Feature identifiers. Every case fails before any solver runs or any output
% folder is created.

    methods (Test)
        function unknownShorthandNamesKeepBenchmarkIdentifier(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            for name = {'a', 'noisy+a', 'noisy+'}
                testCase.verifyError(@() benchmark({@stayAtStart, @stayAtStart}, struct('feature_name', name{1})), ...
                    'MATLAB:benchmark:feature_nameNotValid', sprintf('feature_name ''%s''', name{1}));
            end
            raised = [];
            try
                benchmark({@stayAtStart, @stayAtStart}, struct('feature_name', 'a'));
            catch failure
                raised = failure;
            end
            testCase.assertNotEmpty(raised, 'benchmark must reject an unknown feature_name.');
            testCase.verifyEqual(raised.identifier, 'MATLAB:benchmark:feature_nameNotValid');
            testCase.assertNotEmpty(raised.cause, 'The Feature error must be kept as the cause.');
            testCase.verifyEqual(raised.cause{1}.identifier, 'MATLAB:Feature:UnknownFeature');
            listing = dir(pwd);
            testCase.verifyEmpty(setdiff({listing.name}, {'.', '..'}), 'No output may be created.');
        end

        function otherRoutesKeepFeatureIdentifiers(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            testCase.verifyError(@() benchmark({@stayAtStart, @stayAtStart}, struct('feature', struct('name', 'a'))), ...
                'MATLAB:Feature:UnknownFeature');
            testCase.verifyError(@() benchmark({@stayAtStart, @stayAtStart}, struct('feature_name', 'noisy', 'noise_level', -1)), ...
                'MATLAB:Feature:noise_level_NotPositive');
            listing = dir(pwd);
            testCase.verifyEmpty(setdiff({listing.name}, {'.', '..'}), 'No output may be created.');
        end
    end
end

function x = stayAtStart(fun, x0)
    fun(x0);
    x = x0;
end
