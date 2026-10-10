classdef TestReferenceScoring < matlab.unittest.TestCase
%TESTREFERENCESCORING Raw facts define targets, including solver-free reload.
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
            src = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src');
            addpath(src);
            addpath(fullfile(src, '..', 'tests', 'fixtures', 'evalreport'));
            cd(fullfile(src, 'private'));
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
        function record = reference(fun, maxcv, kind)
            record = struct('fun', fun, 'maxcv', maxcv, 'kind', kind, ...
                'source', 'analytic test', 'mapping', 'objective_violation/1');
        end
        function block = results(record)
            block = struct('problem_names', {{'P'}}, 'problem_dims', 1, ...
                'merit_histories', ones(1, 2, 2, 2), 'merit_outs', ones(1, 2, 2), ...
                'merit_inits', [4, 4], 'maxcv_inits', [1, 2], 'n_evals', 2 * ones(1, 2, 2));
            block.reference_records = repmat({record}, 1, 2);
        end
        function x = finish(fun, x0, varargin)
            x = ones(size(x0));
            fun(x);
        end
        function x = stay(fun, x0, varargin)
            fun(x0);
            x = x0;
        end
        function writeFile(file, code)
            fid = fopen(file, 'w');
            assert(fid >= 0);
            cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
            fprintf(fid, '%s', code);
        end
    end
    methods (Test)
        function rawValidation(testCase)
            for kind = {'lower_bound', 'optimum', 'best_known', 'target'}
                record = TestReferenceScoring.reference(0, 0, kind{1});
                problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', 2, 'reference', record));
                testCase.verifyEqual(problem.reference, record);
                testCase.verifyFalse(isfield(problem.reference, 'point'));
            end
            for kind = {'lower_bound', 'optimum'}
                record = TestReferenceScoring.reference(0, 0.1, kind{1});
                testCase.verifyError(@() normalizeProblemReference(record), 'MATLAB:Problem:reference_maxcv_NotFeasible');
            end
        end
        function targetsUseCurrentMeritPerRun(testCase)
            record = TestReferenceScoring.reference(2, 0.5, 'best_known');
            block = TestReferenceScoring.results(record);
            options = struct('merit_fun', @(f, v, v0) 2*f + 10*v + v0);
            [~, ~, ~, targets] = processResults({block}, options);
            testCase.verifyEqual(targets, reshape([10, 11], 1, 1, 2));
            block = rmfield(block, 'reference_records');
            [~, ~, ~, targets] = processResults({block}, options);
            testCase.verifyEqual(targets, ones(1, 1, 2));
        end
        function plainDoesNotOverrideAnExplicitTarget(testCase)
            block = TestReferenceScoring.results(TestReferenceScoring.reference(2, 0, 'target'));
            plain = TestReferenceScoring.results([]);
            plain.merit_histories(:) = -5;
            block.results_plib_plain = plain;
            [~, ~, ~, targets] = processResults({block}, struct('merit_fun', @(f, v, v0) f, 'run_plain', true));
            testCase.verifyEqual(targets, 2 * ones(1, 1, 2));
        end
        function invalidConvertedMeritIsNotIgnored(testCase)
            block = TestReferenceScoring.results(TestReferenceScoring.reference(0, 0, 'optimum'));
            testCase.verifyError(@() processResults({block}, struct('merit_fun', @(f, v, v0) NaN)), ...
                'MATLAB:applyReferenceTargets:invalidMerit');
            block = TestReferenceScoring.results(TestReferenceScoring.reference(0, 1e3, 'target'));
            testCase.verifyError(@() processResults({block}, struct('merit_fun', @defaultMerit)), ...
                'MATLAB:applyReferenceTargets:invalidMerit');
        end
        function nonzeroViolationIsDroppedAcrossGeneralAffine(testCase)
            record = TestReferenceScoring.reference(0, 0.5, 'best_known');
            problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [1; 2], 'reference', record));
            testCase.verifyEqual(FeaturedProblem(problem, Feature('noisy+permuted'), 10, 7).reference, record);
            testCase.verifyEmpty(FeaturedProblem(problem, Feature('linearly_transformed'), 10, 7).reference);
        end
        function directProblemScoreConsumesReference(testCase)
            solvers = {@TestReferenceScoring.finish, @TestReferenceScoring.stay};
            options = struct('score_only', true, 'silent', true, 'n_jobs', 1, 'max_eval_factor', 2);
            options.problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', 2));
            implicit = benchmark(solvers, options);
            options.problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', 2, ...
                'reference', TestReferenceScoring.reference(0, 0, 'optimum')));
            explicit = benchmark(solvers, options);
            testCase.verifyEqual(implicit, [1; 0]);
            testCase.verifyEqual(explicit, [0.75; 0]);
        end
        function publicSaveLoadUsesReferences(testCase)
            root = fullfile(testCase.WorkDir, 'library');
            mkdir(root);
            TestReferenceScoring.writeFile(fullfile(root, 'refscore_select.m'), sprintf( ...
                'function names = refscore_select(options)\nnames = {''REF'', ''UNKNOWN''};\nend\n'));
            TestReferenceScoring.writeFile(fullfile(root, 'refscore_load.m'), sprintf([ ...
                'function problem = refscore_load(name)\n' ...
                's = struct(''fun'', @(x) sum(x.^2), ''x0'', 2, ''name'', name);\n' ...
                'if strcmp(name,''REF''), s.reference = struct(''fun'',0,''maxcv'',0,''kind'',''optimum'',''source'',''analytic test'',''mapping'',''objective_violation/1''); end\n' ...
                'problem = Problem(s);\nend\n']));
            registerProblemLibrary(struct('name', 'refscore', 'root', root, ...
                'select_function', 'refscore_select', 'load_function', 'refscore_load'));
            cd(testCase.WorkDir);
            % The CI lifecycle lane supplies PCT. Other installations still
            % exercise the serial public contract rather than skip the test.
            jobs = 1;
            if ~isempty(ver('parallel')), jobs = [1, 2]; end
            for n_jobs = jobs
                id = sprintf('saved-jobs-%d', n_jobs);
                options = struct('plibs', {{'refscore'}}, 'ptype', 'u', 'mindim', 1, 'maxdim', 1, ...
                    'n_jobs', n_jobs, 'n_runs', 2, 'max_eval_factor', 2, 'max_tol_order', 2, ...
                    'draw_hist_plots', 'none', 'silent', true, 'benchmark_id', id, 'run_plain', true);
                options.report_path = fullfile(testCase.WorkDir, [id, '-report.json']);
                [scores, profiles] = benchmark({@TestReferenceScoring.finish, @TestReferenceScoring.stay}, options);
                if n_jobs == 2, testCase.verifyNotEmpty(gcp('nocreate')); end
                files = dir(fullfile(testCase.WorkDir, id, '**', 'data_for_loading.mat'));
                testCase.assertEqual(numel(files), 1);
                stored = load(fullfile(files(1).folder, files(1).name));
                records = stored.results_plibs{1}.reference_records;
                testCase.verifyEqual(size(records), [2, 2]);
                testCase.verifyEqual(records(1, :), repmat({TestReferenceScoring.reference(0, 0, 'optimum')}, 1, 2));
                testCase.verifyTrue(all(cellfun(@isempty, records(2, :))));
                report = jsondecode(fileread(options.report_path));
                schema = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', '..', ...
                    'python', 'optiprofiler', 'schemas', 'eval_report-v2.schema.json');
                testCase.verifyEmpty(evalReportSchemaCheck(options.report_path, schema));
                for p = 1:numel(report.problems)
                    if iscell(report.problems), entry = report.problems{p};
                    else, entry = report.problems(p); end
                    if strcmp(entry.name, 'REF')
                        testCase.verifyEqual(entry.reference_records(1).fun, 0);
                        testCase.verifyEqual(entry.reference_records(1).maxcv, 0);
                    end
                end
                load_options = struct('load', 'latest', 'benchmark_id', id, 'score_only', true, ...
                    'silent', true, 'max_tol_order', 2, 'max_eval_factor', 2, 'run_plain', true);
                [loaded_scores, loaded_profiles] = benchmark({}, load_options);
                testCase.verifyEqual(loaded_scores, scores);
                testCase.verifyEqual(loaded_profiles, profiles);
                load_options.merit_fun = @(f, v, v0) f*f + v;
                [~, changed_profiles] = benchmark({}, load_options);
                testCase.verifyFalse(isequaln(changed_profiles, profiles));
                load_options.problem_names = {'REF'};
                benchmark({}, load_options);
            end
        end
    end
end
