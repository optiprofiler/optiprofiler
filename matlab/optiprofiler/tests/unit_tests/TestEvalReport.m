classdef TestEvalReport < matlab.unittest.TestCase
% Public benchmark/report regression: no private collector API assertions.
    properties (Access = private)
        SourceRoot
        OutputRoot
    end
    methods (TestMethodSetup)
        function isolatePublicFixture(testCase)
            folder = fileparts(mfilename('fullpath'));
            [~, root] = fileattrib(fullfile(folder, '..', '..', '..', '..'));
            testCase.SourceRoot = root.Name;
            testCase.OutputRoot = tempname;
            mkdir(testCase.OutputRoot);
            original_path = path;
            original_directory = pwd;
            registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            testCase.addTeardown(@() path(original_path));
            testCase.addTeardown(@() cd(original_directory));
            testCase.addTeardown(@() setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', registry));
            testCase.addTeardown(@() removeOutput(testCase.OutputRoot));
            addpath(fullfile(folder, '..', 'fixtures', 'evalreport'));
        end
    end
    methods (Test)
        function numericalIdentityAndUtf8(testCase), testCase.runCase('identity'); end
        function coverageAndNoOverwrite(testCase), testCase.runCase('coverage'); end
        function nonfiniteAndScoringFailure(testCase), testCase.runCase('semantics'); end
        function compactExactExtremaAndCompanion(testCase), testCase.runCase('compact'); end
        function paddedLengthAggregationBoundary(testCase), testCase.runCase('aggregation'); end
        function statefulRendererMeritObserved(testCase), testCase.runCase('stateful_render'); end
        function renderingFailureSeparate(testCase), testCase.runCase('render'); end
        function archiveProvenanceAndArtifacts(testCase), testCase.runCase('archive'); end
        function actualParallelEquivalence(testCase)
            % License preflight can return false even when pool checkout
            % succeeds. Skip only absent installations; the fixture requires
            % real worker execution and fails if the runtime cannot supply it.
            testCase.assumeFalse(isempty(ver('parallel')), 'Parallel Computing Toolbox is not installed.');
            testCase.runCase('parallel');
        end
        function writeFailureAndOwnership(testCase), testCase.runCase('ownership'); end
        function optionalPlainReferenceStatus(testCase), testCase.runCase('plain_reference'); end
    end
    methods (Access = private)
        function runCase(testCase, name)
            evalReportPublic(testCase.SourceRoot, testCase.OutputRoot, name);
        end
    end
end

function removeOutput(path)
    if isfolder(path), rmdir(path, 's'); end
end
