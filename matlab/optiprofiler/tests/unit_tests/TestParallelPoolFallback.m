classdef TestParallelPoolFallback < matlab.unittest.TestCase
%TESTPARALLELPOOLFALLBACK A run that cannot open a parallel pool says so in its report.
    properties
        OldCwd
        OldPath
        OldRegistry
        WorkDir
    end
    methods (TestMethodSetup)
        function setupState(testCase)
            testCase.OldCwd = pwd;
            testCase.OldPath = path;
            testCase.OldRegistry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            testCase.WorkDir = tempname;
            mkdir(testCase.WorkDir);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(testCase.WorkDir, 'registry.mat'));
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
            cd(testCase.WorkDir);
        end
    end
    methods (TestMethodTeardown)
        function restoreState(testCase)
            cd(testCase.OldCwd);
            path(testCase.OldPath);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', testCase.OldRegistry);
            if isfolder(testCase.WorkDir), rmdir(testCase.WorkDir, 's'); end
        end
    end
    methods (Static)
        function x = halve(fun, x0, varargin)
            fun(x0);
            x = x0 / 2;
            fun(x);
        end
        function writeFile(file, code)
            fid = fopen(file, 'w');
            assert(fid >= 0);
            cleanup = onCleanup(@() fclose(fid));
            fprintf(fid, '%s', code);
        end
    end
    methods (Test)
        function poolThatCannotStartIsRecordedAndTheRunIsSequential(testCase)
            % No pool exists, and opening one fails, as without a licence seat for a worker.
            stubs = fullfile(testCase.WorkDir, 'stubs');
            mkdir(stubs);
            TestParallelPoolFallback.writeFile(fullfile(stubs, 'gcp.m'), strjoin({ ...
                'function pool = gcp(varargin)', 'pool = [];', 'end', ''}, newline));
            TestParallelPoolFallback.writeFile(fullfile(stubs, 'parpool.m'), strjoin({ ...
                'function varargout = parpool(varargin)', ...
                'error(''TestParallelPoolFallback:NoSeat'', ''no licence seat'');', 'end', ''}, newline));
            addpath(stubs);
            root = fullfile(testCase.WorkDir, 'library');
            mkdir(root);
            TestParallelPoolFallback.writeFile(fullfile(root, 'pool_fallback_library_select.m'), sprintf( ...
                'function names = pool_fallback_library_select(options)\n names = {''A'', ''B''};\n end\n'));
            TestParallelPoolFallback.writeFile(fullfile(root, 'pool_fallback_library_load.m'), sprintf([ ...
                'function problem = pool_fallback_library_load(name)\n' ...
                'x0 = 2;\nif strcmp(name, ''B''), x0 = 3; end\n' ...
                'problem = Problem(struct(''fun'', @(x) sum(x.^2), ''x0'', x0, ''name'', name));\nend\n']));
            registerProblemLibrary(struct('name', 'pool_fallback_library', 'root', root, ...
                'select_function', 'pool_fallback_library_select', 'load_function', 'pool_fallback_library_load'));
            report_path = fullfile(testCase.WorkDir, 'report.json');
            scores = benchmark({@TestParallelPoolFallback.halve, @TestParallelPoolFallback.halve}, ...
                struct('plibs', {{'pool_fallback_library'}}, 'ptype', 'u', 'mindim', 1, 'maxdim', 1, ...
                'n_jobs', 2, 'max_tol_order', 1, 'score_only', true, 'silent', true, 'solver_verbose', 0, ...
                'report_path', report_path));
            testCase.verifyNumElements(scores, 2);
            report = jsondecode(fileread(report_path));
            diagnostics = report.diagnostics;
            if isstruct(diagnostics), diagnostics = num2cell(diagnostics); end
            if isempty(diagnostics), diagnostics = {}; end
            codes = cellfun(@(d) d.code, diagnostics, 'UniformOutput', false);
            match = diagnostics(strcmp(codes, 'parallel_pool_unavailable'));
            testCase.assertNumElements(match, 1);
            testCase.verifyEqual(match{1}.stage, 'numerical');
            testCase.verifyEqual(match{1}.scope.library, 'pool_fallback_library');
            testCase.verifyEqual(match{1}.scope.requested_workers, 2);
            testCase.verifyEqual(match{1}.scope.reason, 'parallel_pool_start_failed');
        end
    end
end
