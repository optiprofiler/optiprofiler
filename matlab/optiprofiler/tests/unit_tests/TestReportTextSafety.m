classdef TestReportTextSafety < matlab.unittest.TestCase
%TESTREPORTTEXTSAFETY Report provenance is bounded honestly and keeps no credentials.
%   Mirrors Python's test_report_text_safety.py.
    properties
        OldCwd
        OldPath
        WorkDir
    end
    methods (TestMethodSetup)
        function setupState(testCase)
            testCase.OldCwd = pwd;
            testCase.OldPath = path;
            testCase.WorkDir = tempname;
            mkdir(testCase.WorkDir);
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
            cd(testCase.WorkDir);
        end
    end
    methods (TestMethodTeardown)
        function restoreState(testCase)
            cd(testCase.OldCwd);
            path(testCase.OldPath);
            if isfolder(testCase.WorkDir), rmdir(testCase.WorkDir, 's'); end
        end
    end
    methods (Static)
        function x = finish(fun, x0, varargin)
            x = ones(size(x0));
            fun(x);
        end
        function x = stay(fun, x0, varargin)
            fun(x0);
            x = x0;
        end
        function source = reportedSource(work_dir, source)
            reference = struct('fun', 0, 'maxcv', 0, 'kind', 'feasible', 'source', source);
            options = struct('score_only', true, 'silent', true, 'n_jobs', 1, 'max_eval_factor', 2);
            options.problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', 2, 'reference', reference));
            options.report_path = [tempname(work_dir), '.json'];
            benchmark({@TestReportTextSafety.finish, @TestReportTextSafety.stay}, options);
            report = jsondecode(fileread(options.report_path));
            problem = report.problems;
            if iscell(problem), problem = problem{1}; else, problem = problem(1); end
            records = problem.reference_records;
            if iscell(records), record = records{1}; else, record = records(1); end
            source = record.source;
        end
    end
    methods (Test)
        function referenceSourceKeepsLocationWithoutSecrets(testCase)
            reported = TestReportTextSafety.reportedSource(testCase.WorkDir, ...
                'https://user:pw@doi.org/10.1080/10556788.2025.2490640?utm=mail#top');
            testCase.verifyEqual(reported, ...
                'https://doi.org/10.1080/10556788.2025.2490640 [url parameters removed]');
            reported = TestReportTextSafety.reportedSource(testCase.WorkDir, ...
                'https://doi.org/10.1080/10556788.2025.2490640');
            testCase.verifyEqual(reported, 'https://doi.org/10.1080/10556788.2025.2490640');
        end
        function longReferenceSourceIsMarkedWhenClipped(testCase)
            reported = TestReportTextSafety.reportedSource(testCase.WorkDir, ...
                ['best known value reported by ', repmat('y', 1, 300)]);
            testCase.verifyEqual(numel(reported), 256);
            testCase.verifyTrue(endsWith(reported, '...[truncated]'));
        end
    end
end
