classdef TestUndrawableProfiles < matlab.unittest.TestCase
%TESTUNDRAWABLEPROFILES A profile branch in which no solver meets the convergence test scores 0.
%   Mirrors Python's test_undrawable_profiles.py.
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
        function x = toThreeTenths(fun, x0, varargin)
            fun(x0);
            x = 0.3 * x0;
            fun(x);
        end
        function x = toAQuarter(fun, x0, varargin)
            fun(x0);
            x = x0 / 4;
            fun(x);
        end
        function x = crashBefore(~, ~, varargin) %#ok<STOUT>
            error('TestUndrawableProfiles:Crash', 'crashes before its first evaluation');
        end
        function writeFile(file, code)
            fid = fopen(file, 'w');
            assert(fid >= 0);
            cleanup = onCleanup(@() fclose(fid));
            fprintf(fid, '%s', code);
        end
    end
    methods (Access = private)
        function options = library(testCase)
            root = fullfile(testCase.WorkDir, 'library');
            mkdir(root);
            TestUndrawableProfiles.writeFile(fullfile(root, 'undrawable_library_select.m'), sprintf( ...
                'function names = undrawable_library_select(options)\n names = {''A'', ''B''};\n end\n'));
            TestUndrawableProfiles.writeFile(fullfile(root, 'undrawable_library_load.m'), sprintf([ ...
                'function problem = undrawable_library_load(name)\n' ...
                'x0 = 2;\nif strcmp(name, ''B''), x0 = 3; end\n' ...
                's = struct(''fun'', @(x) sum(x.^2), ''x0'', x0, ''name'', name);\n' ...
                's.reference = struct(''fun'', 0, ''maxcv'', 0, ''kind'', ''feasible'', ''source'', ''analytic'');\n' ...
                'problem = Problem(s);\nend\n']));
            registerProblemLibrary(struct('name', 'undrawable_library', 'root', root, ...
                'select_function', 'undrawable_library_select', 'load_function', 'undrawable_library_load'));
            options = struct('plibs', {{'undrawable_library'}}, 'ptype', 'u', 'mindim', 1, 'maxdim', 1, ...
                'n_jobs', 1, 'max_eval_factor', 2, 'max_tol_order', 3, 'score_only', true, 'silent', true, ...
                'solver_verbose', 0, 'report_path', fullfile(testCase.WorkDir, 'report.json'));
        end
    end
    methods (Test)
        function tolerancesNoSolverReachesScoreZeroInEveryProfile(testCase)
            % Both solvers reach 9 % and 6 % of f(x0) against the reference 0,
            % so only the loosest tolerance, 1e-1, has a solver that meets the test.
            options = testCase.library();
            [~, profiles, curves] = benchmark({@TestUndrawableProfiles.toThreeTenths, @TestUndrawableProfiles.toAQuarter}, options);
            testCase.verifyTrue(all(profiles(:, 1, 1, 1:2) > 0, 'all'));
            testCase.verifyEqual(profiles(:, 2:3, :, :), zeros(2, 2, 2, 3));
            for i_tol = 2:3
                for branch = {'hist', 'out'}
                    curve = curves{i_tol}.(branch{1});
                    testCase.verifyTrue(all(cellfun(@isempty, [curve.perf(:); curve.data(:); curve.log_ratio(:)])));
                end
            end
            refs = plotRefs(options.report_path);
            testCase.verifyNotEmpty(refs);
            testCase.verifyTrue(all(startsWith(refs, 'profile-1-')));
        end
        function solversThatAllCrashScoreZeroInEveryProfile(testCase)
            options = testCase.library();
            [~, profiles] = benchmark({@TestUndrawableProfiles.crashBefore, @TestUndrawableProfiles.crashBefore}, options);
            testCase.verifyEqual(profiles, zeros(2, 3, 2, 3));
            testCase.verifyEmpty(plotRefs(options.report_path));
        end
        function readmeListsTheSummaryPdfs(testCase)
            % Only the first of the three tolerances has profiles.
            options = testCase.library();
            options.score_only = false;
            options.draw_hist_plots = 'none';
            options.savepath = fullfile(testCase.WorkDir, 'out');
            benchmark({@TestUndrawableProfiles.toThreeTenths, @TestUndrawableProfiles.toAQuarter}, options);
            summary = dir(fullfile(options.savepath, '**', 'summary_*.pdf'));
            testCase.assertNumElements(summary, 1);
            % README.txt wraps a long name over the first 33 columns of its following lines.
            lines = splitlines(fileread(fullfile(summary.folder, 'README.txt')));
            readme = strjoin(cellfun(@(line) strtrim(line(1:min(33, end))), lines, 'UniformOutput', false), '');
            testCase.verifySubstring(readme, summary.name);
            for name = {'perf_hist.pdf', 'data_hist.pdf', 'log-ratio_hist.pdf', 'perf_out.pdf', 'data_out.pdf', 'log-ratio_out.pdf'}
                testCase.verifyTrue(isfile(fullfile(summary.folder, name{1})), name{1});
                testCase.verifySubstring(readme, name{1});
            end
        end
    end
end

function refs = plotRefs(report_path)
    report = jsondecode(fileread(report_path));
    refs = report.profiles.plot_refs;
    if isempty(refs)
        refs = {};
    end
    refs = cellstr(refs);
end
