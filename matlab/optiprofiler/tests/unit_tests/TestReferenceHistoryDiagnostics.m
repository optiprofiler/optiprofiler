classdef TestReferenceHistoryDiagnostics < matlab.unittest.TestCase
    methods (Test)
        function historyFailureWarnsOncePerTrial(testCase)
            state = warning;
            cleanup = onCleanup(@() warning(state)); %#ok<NASGU>
            id = 'OptiProfiler:FeaturedProblem:ReferenceHistoryUnavailable';
            warning('on', id);
            for name = {'plain', 'truncated+truncated'}
                p = Problem(struct('fun', @(x) sum(x), 'x0', [0;0], ...
                    'cub', @TestReferenceHistoryDiagnostics.malformed));
                fp = FeaturedProblem(p, Feature(name{1}), 3, 17);
                testCase.verifyWarning(@() fp.fun([1;1]), id);
                lastwarn('');
                testCase.verifyEqual(fp.fun([1;1]), 2);
                testCase.verifyEqual(fp.fun([1;1]), 2);
                testCase.verifyTrue(all(isnan(fp.maxcv_hist)));
                [message, ~] = lastwarn;
                testCase.verifyEmpty(message);
                failure = fp.historyFailure();
                testCase.verifyEqual(failure.identifier, 'MATLAB:Problem:cubx_m_nonlinear_ub_NotConsistent');
                testCase.verifyNotEmpty(failure.message);
            end
        end
        function publicBenchmarkReportsSuppressedFailure(testCase)
            p = Problem(struct('fun', @(x) sum(x), 'x0', [0;0], ...
                'cub', @TestReferenceHistoryDiagnostics.malformed));
            options = struct('problem',p,'score_only',true,'silent',false, ...
                'solver_verbose',0,'n_jobs',1,'n_runs',1,'run_plain',false, ...
                'max_tol_order',1,'max_eval_factor',2,'project_x0',false);
            text = evalc('benchmark({@TestReferenceHistoryDiagnostics.probeSolver, @TestReferenceHistoryDiagnostics.probeSolver}, options);');
            testCase.verifySubstring(text,'Reference history unavailable');
            testCase.verifySubstring(text,'cubx_m_nonlinear_ub_NotConsistent');
        end
        function disabledWarningsKeepFailureStatus(testCase)
            state = warning;
            cleanup = onCleanup(@() warning(state)); %#ok<NASGU>
            warning('off', 'all');
            p = Problem(struct('fun', @(x) sum(x), 'x0', [0;0], ...
                'cub', @TestReferenceHistoryDiagnostics.malformed));
            fp = FeaturedProblem(p, Feature('plain'), 3, 17);
            testCase.verifyEqual(fp.fun([1;1]), 2);
            testCase.verifyNotEmpty(fp.historyFailure());
            testCase.verifyTrue(isnan(fp.maxcv_hist(end)));
        end
        function silentSavedBenchmarkRetainsActualTrialFailures(testCase)
            work = tempname;
            mkdir(work);
            cleanup = onCleanup(@() rmdir(work, 's')); %#ok<NASGU>
            for report_enabled = [false, true]
                destination = fullfile(work, sprintf('report-%d', report_enabled));
                mkdir(destination);
                p = Problem(struct('fun', @(x) sum(x), 'x0', [0;0], ...
                    'cub', @TestReferenceHistoryDiagnostics.malformed));
                options = struct('problem',p,'score_only',false,'silent',true, ...
                    'solver_verbose',0,'n_jobs',1,'n_runs',3,'run_plain',false, ...
                    'max_tol_order',1,'max_eval_factor',2,'project_x0',false, ...
                    'draw_hist_plots','none','savepath',destination);
                if report_enabled
                    options.report_path = fullfile(destination, 'report.json');
                end
                solvers = {@TestReferenceHistoryDiagnostics.repeatedProbeSolver, ...
                    @TestReferenceHistoryDiagnostics.repeatedProbeSolver};
                text = evalc('benchmark(solvers, options);');
                testCase.verifyFalse(contains(text, 'Reference history unavailable'));
                testCase.verifyFalse(contains(text, 'cubx_m_nonlinear_ub_NotConsistent'));
                logs = dir(fullfile(destination, '**', 'log.txt'));
                testCase.assertEqual(numel(logs), 1);
                log_text = fileread(fullfile(logs(1).folder, logs(1).name));
                % Plain repetitions are copied slots, not new actual trials.
                testCase.verifyEqual(count(log_text, 'Reference history unavailable'), 2);
                testCase.verifyEqual(count(log_text, 'cubx_m_nonlinear_ub_NotConsistent'), 2);
                if report_enabled
                    report = jsondecode(fileread(options.report_path));
                    diagnostic = report.diagnostics;
                    testCase.assertFalse(isempty(diagnostic));
                    diagnostic = diagnostic(strcmp({diagnostic.code}, 'reference_history_unavailable'));
                    testCase.assertEqual(numel(diagnostic), 2);
                    testCase.verifyEqual(arrayfun(@(d) d.scope.run_index, diagnostic(:)), [1;1]);
                    testCase.verifyEqual(sort(arrayfun(@(d) d.scope.solver_index, diagnostic(:))), [1;2]);
                    for k = 1:numel(diagnostic)
                        testCase.verifyEqual(diagnostic(k).scope.exception_type, ...
                            'MATLAB:Problem:cubx_m_nonlinear_ub_NotConsistent');
                    end
                    testCase.verifyEqual(report.stages.numerical.status, 'completed');
                end
            end
        end
    end
    methods (Static)
        function x = probeSolver(fun,x0,varargin)
            fun([1;1]); x=x0;
        end
        function c = malformed(x)
            if all(x == 0), c = 1; else, c = [1;2]; end
        end
        function x = repeatedProbeSolver(fun,x0,varargin)
            fun([1;1]); fun([1;1]); fun([1;1]); x=x0;
        end
    end
end
