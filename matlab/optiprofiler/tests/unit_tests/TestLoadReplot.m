classdef TestLoadReplot < matlab.unittest.TestCase
%TESTLOADREPLOT Regression tests for the load/replot path of `benchmark`.
%   The tests synthesize a saved experiment directory (no solver execution and
%   no problem library required) and exercise the load path end to end.

    properties (Access = private)
        CurrentDirectory
        WorkDir
    end

    methods (TestMethodSetup)

        function setupWorkDir(testCase)
            testCase.CurrentDirectory = pwd;
            test_dir = fileparts(mfilename('fullpath'));
            source_dir = fullfile(test_dir, '..', '..', 'src');
            addpath(source_dir);
            testCase.WorkDir = tempname;
            mkdir(testCase.WorkDir);
            cd(testCase.WorkDir);
        end

    end

    methods (TestMethodTeardown)

        function restoreState(testCase)
            cd(testCase.CurrentDirectory);
            if exist(testCase.WorkDir, 'dir')
                try
                    rmdir(testCase.WorkDir, 's');
                catch
                end
            end
        end

    end

    methods (Static)

        function makeExperiment(work, benchmark_id, n_solvers, n_runs)
            %MAKEEXPERIMENT Synthesize a saved experiment that `benchmark` can load.
            %   The saved problem selection is deliberately non-default
            %   ('u', mindim 2, maxdim 3) so that stamp regressions are visible.
            time_stamp = '20200101_000000';
            n_problems = 2;
            n_evals_max = 8;
            fun_histories = zeros(n_problems, n_solvers, n_runs, n_evals_max);
            for i_problem = 1:n_problems
                for i_solver = 1:n_solvers
                    for i_run = 1:n_runs
                        start_value = 10 + i_problem + i_solver + 0.1 * i_run;
                        fun_histories(i_problem, i_solver, i_run, :) = linspace(start_value, 0.5 * i_solver, n_evals_max);
                    end
                end
            end
            maxcv_histories = zeros(size(fun_histories));

            results = struct();
            results.plib = 'synthetic';
            results.solver_names = arrayfun(@(k) sprintf('solver%d', k), 1:n_solvers, 'UniformOutput', false);
            results.ptype = 'u';
            results.mindim = 2;
            results.maxdim = 3;
            results.minb = 0;
            results.maxb = 10;
            results.minlcon = 0;
            results.maxlcon = 10;
            results.minnlcon = 0;
            results.maxnlcon = 10;
            results.mincon = 0;
            results.maxcon = 10;
            results.problem_names_options = {};
            results.excludelist = {};
            results.feature_stamp = 'plain';
            results.fun_histories = fun_histories;
            results.maxcv_histories = maxcv_histories;
            results.fun_outs = fun_histories(:, :, :, end);
            results.maxcv_outs = maxcv_histories(:, :, :, end);
            results.fun_inits = reshape(fun_histories(:, 1, :, 1), n_problems, n_runs);
            results.maxcv_inits = zeros(n_problems, n_runs);
            results.n_evals = repmat(n_evals_max, n_problems, n_solvers, n_runs);
            results.problem_names = {'SYNPROB1', 'SYNPROB2'};
            results.problem_types = {'u', 'u'};
            results.problem_dims = [2; 3];
            results.problem_mbs = zeros(n_problems, 1);
            results.problem_mlcons = zeros(n_problems, 1);
            results.problem_mnlcons = zeros(n_problems, 1);
            results.problem_mcons = zeros(n_problems, 1);
            results.computation_times = 0.01 * ones(n_problems, n_solvers, n_runs);
            results.solvers_successes = true(n_problems, n_solvers, n_runs);
            % With the default merit function and zero constraint violation the
            % merit values equal the objective values.
            results.merit_histories = fun_histories;
            results.merit_outs = results.fun_outs;
            results.merit_inits = results.fun_inits;

            results_plibs = {results};
            path_log = fullfile(work, benchmark_id, ['synthetic_plain_u_2_3_', time_stamp], 'test_log');
            mkdir(path_log);
            save(fullfile(path_log, 'data_for_loading.mat'), 'results_plibs');
            fid = fopen(fullfile(path_log, ['time_stamp_', time_stamp, '.txt']), 'w');
            fclose(fid);
        end

        function options = loadOptions(varargin)
            %LOADOPTIONS Common options for loading the synthetic experiment.
            options = struct('load', '20200101_000000', 'benchmark_id', 'bench', ...
                'max_tol_order', 1, 'silent', true);
            for i = 1:2:numel(varargin)
                options.(varargin{i}) = varargin{i + 1};
            end
        end

    end

    methods (Test)

        function twoSolverLoadGeneratesLogRatio(testCase)
            TestLoadReplot.makeExperiment(pwd, 'bench', 2, 1);
            options = TestLoadReplot.loadOptions('summarize_log_ratio_profiles', true);
            scores = benchmark(options);
            testCase.verifyEqual(numel(scores), 2);
            pdfs = dir(fullfile(pwd, 'bench', '**', 'log-ratio_hist_*.pdf'));
            testCase.verifyNotEmpty(pdfs, 'The log-ratio profile PDFs are missing after a two-solver load.');
        end

        function threeSolverLoadWarnsAndSkipsLogRatio(testCase)
            TestLoadReplot.makeExperiment(pwd, 'bench', 3, 1);
            options = TestLoadReplot.loadOptions('summarize_log_ratio_profiles', true);
            testCase.verifyWarning(@() benchmark(options), 'MATLAB:benchmark:summarize_log_ratio_profilesOnlyWhenTwoSolvers');
            pdfs = dir(fullfile(pwd, 'bench', '**', 'log-ratio_hist_*.pdf'));
            testCase.verifyEmpty(pdfs, 'Log-ratio profile PDFs were generated for a three-solver load.');
        end

        function selectTwoOfThreeGeneratesLogRatio(testCase)
            TestLoadReplot.makeExperiment(pwd, 'bench', 3, 1);
            options = TestLoadReplot.loadOptions('summarize_log_ratio_profiles', true, 'solvers_to_load', [1, 3]);
            scores = benchmark(options);
            testCase.verifyEqual(numel(scores), 2);
            pdfs = dir(fullfile(pwd, 'bench', '**', 'log-ratio_hist_*.pdf'));
            testCase.verifyNotEmpty(pdfs, 'The log-ratio profile PDFs are missing after selecting two of three solvers.');
        end

        function solverNamesWrongLengthErrors(testCase)
            TestLoadReplot.makeExperiment(pwd, 'bench', 2, 1);
            options = TestLoadReplot.loadOptions('solver_names', {'only_one'});
            testCase.verifyError(@() benchmark(options), 'MATLAB:benchmark:solver_namesAndLoadedSolversLengthNotSame');
        end

        function singleSolversToLoadErrors(testCase)
            TestLoadReplot.makeExperiment(pwd, 'bench', 2, 1);
            options = TestLoadReplot.loadOptions('solvers_to_load', 1);
            testCase.verifyError(@() benchmark(options), 'MATLAB:checkValidityProfileOptions:solvers_to_loadNotValid');
        end

    end

end
