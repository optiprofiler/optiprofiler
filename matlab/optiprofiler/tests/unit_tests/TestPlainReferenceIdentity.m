classdef TestPlainReferenceIdentity < matlab.unittest.TestCase
%TESTPLAINREFERENCEIDENTITY Keep plain references within their saved library.

    properties (Access = private)
        CurrentDirectory
        OriginalPath
        WarningState
        WorkDir
    end

    methods (TestMethodSetup)
        function setupPrivateFunctions(testCase)
            testCase.CurrentDirectory = pwd;
            testCase.OriginalPath = path;
            testCase.WarningState = warning;
            test_dir = fileparts(mfilename('fullpath'));
            src = fullfile(test_dir, '..', '..', 'src');
            addpath(src);
            cd(fullfile(src, 'private'));
            testCase.WorkDir = tempname(pwd);
            mkdir(testCase.WorkDir);
        end
    end

    methods (TestMethodTeardown)
        function restoreState(testCase)
            cd(testCase.CurrentDirectory);
            path(testCase.OriginalPath);
            warning(testCase.WarningState);
            if isfolder(testCase.WorkDir)
                rmdir(testCase.WorkDir, 's');
            end
        end
    end

    methods (Test)
        function separateLibrariesWithSameName(testCase)
            blocks = {TestPlainReferenceIdentity.paired('library_a', 5, 1), ...
                TestPlainReferenceIdentity.paired('library_b', 20, 15)};
            testCase.verifyEqual(TestPlainReferenceIdentity.mins(blocks), [1, 1; 15, 15]);
        end

        function separateVariantsOfSameProvider(testCase)
            blocks = {TestPlainReferenceIdentity.paired('provider', 5, 1), ...
                TestPlainReferenceIdentity.paired('provider', 20, 15)};
            blocks{1}.plib_options = struct('n', 2);
            blocks{2}.plib_options = struct('n', 10);
            testCase.verifyEqual(TestPlainReferenceIdentity.mins(blocks), [1, 1; 15, 15]);
        end

        function reorderedAndMissingPlainRows(testCase)
            block = TestPlainReferenceIdentity.block('library_a', {'A', 'B', 'C'}, [10, 20, 30]);
            block.results_plib_plain = TestPlainReferenceIdentity.block('library_a', {'C', 'A'}, [3, 1]);
            testCase.verifyEqual(TestPlainReferenceIdentity.mins({block}), [1, 1; 20, 20; 3, 3]);
        end

        function absentPlainLibrary(testCase)
            blocks = {TestPlainReferenceIdentity.block('library_a', {'DUP'}, 5), ...
                TestPlainReferenceIdentity.paired('library_b', 20, 15)};
            testCase.verifyEqual(TestPlainReferenceIdentity.mins(blocks), [5, 5; 15, 15]);
        end

        function emptyPlainLibrary(testCase)
            first = TestPlainReferenceIdentity.block('library_a', {'DUP'}, 5);
            first.results_plib_plain = [];
            blocks = {first, TestPlainReferenceIdentity.paired('library_b', 20, 15)};
            testCase.verifyEqual(TestPlainReferenceIdentity.mins(blocks), [5, 5; 15, 15]);
        end

        function invalidPlainKeepsValidMain(testCase)
            blocks = {TestPlainReferenceIdentity.paired('library_a', 5, NaN)};
            testCase.verifyEqual(TestPlainReferenceIdentity.mins(blocks), [5, 5]);
        end

        function usesAllPlainRunsAndSolvers(testCase)
            block = TestPlainReferenceIdentity.paired('library_a', 20, 15);
            block.results_plib_plain.merit_histories(1, 3, 2, 3) = 1;
            testCase.verifyEqual(TestPlainReferenceIdentity.mins({block}), [1, 1]);
        end

        function disabledPlain(testCase)
            block = TestPlainReferenceIdentity.paired('library_a', 5, 1);
            [~, ~, ~, mins] = processResults({block}, struct('run_plain', false));
            testCase.verifyEqual(reshape(mins, 1, []), [5, 5]);
        end

        function noPlainOptionOrData(testCase)
            block = TestPlainReferenceIdentity.block('library_a', {'A'}, 5);
            [~, ~, ~, mins] = processResults({block}, struct());
            testCase.verifyEqual(reshape(mins, 1, []), [5, 5]);
        end

        function reportAcceptsPartialPlainResults(testCase)
            block = TestPlainReferenceIdentity.block('library_a', {'A', 'B'}, [10, 20]);
            block.results_plib_plain = TestPlainReferenceIdentity.block('library_a', {'B'}, 15);
            report = fullfile(testCase.WorkDir, 'report.txt');
            readme = fullfile(testCase.WorkDir, 'README.txt');
            fid = fopen(readme, 'w');
            fclose(fid);
            options = struct('score_only', false, 'run_plain', true, 'silent', false);
            writeReport(options, {block}, report, readme);
            text = fileread(report);
            testCase.verifySubstring(text, 'Number of problems selected: 2');
            testCase.verifySubstring(text, 'Wall-clock time spent by all the solvers: 18.00 secs');
        end

        function ambiguousNamesWithinLibrary(testCase)
            block = TestPlainReferenceIdentity.block('library_a', {'A'}, 10);
            block.results_plib_plain = TestPlainReferenceIdentity.block('library_a', {'A', 'A'}, [1, 2]);
            testCase.verifyError(@() TestPlainReferenceIdentity.mins({block}), ...
                'MATLAB:processResults:duplicateProblemNames');
            block = TestPlainReferenceIdentity.block('library_a', {'A', 'A'}, [10, 20]);
            block.results_plib_plain = TestPlainReferenceIdentity.block('library_a', {'A'}, 1);
            testCase.verifyError(@() TestPlainReferenceIdentity.mins({block}), ...
                'MATLAB:processResults:duplicateProblemNames');
        end

        function savedMatReloadFiltersSolversAndKeepsSource(testCase)
            results_plibs = {TestPlainReferenceIdentity.paired('library_a', 5, 1), ...
                TestPlainReferenceIdentity.paired('library_b', 20, 15)};
            results_plibs{2}.results_plib_plain.merit_histories(:, 1, :, :) = -100;
            log_dir = fullfile(testCase.WorkDir, 'saved', 'log');
            mkdir(log_dir);
            mat_file = fullfile(log_dir, 'data_for_loading.mat');
            save(mat_file, 'results_plibs', '-v7.3');
            fid = fopen(fullfile(log_dir, 'time_stamp_20200101_000000.txt'), 'w');
            fclose(fid);
            before = TestPlainReferenceIdentity.readBytes(mat_file);
            [~, work_name] = fileparts(testCase.WorkDir);
            options = struct('load', '20200101_000000', 'benchmark_id', work_name, ...
                'solvers_to_load', [3, 2], 'run_plain', true);
            [loaded, loaded_options] = loadResults(struct('problem_names', {{'DUP'}}), options);
            testCase.verifyEqual(loaded_options.solver_names, {'c', 'b'});
            testCase.verifyEqual(TestPlainReferenceIdentity.mins(loaded), [1, 1; 15, 15]);
            % Exercise the public load entry point as well: no solver runs or
            % figures are needed to recompute the profiles from this MAT file.
            options.score_only = true;
            options.draw_hist_plots = 'none';
            options.max_tol_order = 1;
            options.silent = true;
            scores = benchmark(options);
            testCase.verifyEqual(numel(scores), 2);
            testCase.verifyTrue(all(isfinite(scores)));
            testCase.verifyEqual(TestPlainReferenceIdentity.readBytes(mat_file), before);
        end

        function legacyCachedMeritsCannotHideRawNaN(testCase)
            block = TestPlainReferenceIdentity.paired('library_a', 5, 1);
            block.fun_histories(:, :, :, 1) = NaN;
            block.merit_histories(:, :, :, 1) = -100;
            block.fun_outs(:) = NaN;
            block.merit_outs(:) = 0;
            block.maxcv_inits = NaN(size(block.problem_dims));
            block.merit_inits(:) = 0;
            block.results_plib_plain.maxcv_histories(:, :, :, 1) = NaN;
            block.results_plib_plain.merit_histories(:, :, :, 1) = -200;
            results_plibs = {block};
            mat_file = fullfile(testCase.WorkDir, 'cached-merits.mat');
            save(mat_file, 'results_plibs', '-v7.3');
            before = TestPlainReferenceIdentity.readBytes(mat_file);
            saved = load(mat_file, 'results_plibs');
            [~, outs, inits, mins] = processResults(saved.results_plibs, struct('run_plain', true));
            testCase.verifyEqual(reshape(mins, 1, []), [1, 1]);
            testCase.verifyTrue(all(isinf(outs(:)) & outs(:) > 0));
            testCase.verifyTrue(all(isinf(inits(:)) & inits(:) > 0));
            testCase.verifyEqual(TestPlainReferenceIdentity.readBytes(mat_file), before);
        end

        function meritOnlyHelperDataRemainsSupported(testCase)
            block = TestPlainReferenceIdentity.paired('library_a', 5, 1);
            fields = {'fun_histories', 'maxcv_histories', 'fun_outs', ...
                'maxcv_outs', 'fun_inits', 'maxcv_inits'};
            block = rmfield(block, fields);
            block.results_plib_plain = rmfield(block.results_plib_plain, fields);
            testCase.verifyEqual(TestPlainReferenceIdentity.mins({block}), [1, 1]);
        end

        function legacySingleInitMeritPreservesPerRunValidity(testCase)
            block = TestPlainReferenceIdentity.block('library_a', {'A', 'B'}, [5, 15]);
            block.merit_inits = [5; 15];
            block.maxcv_inits(2, 2) = NaN;
            [~, ~, inits] = processResults({block}, struct('run_plain', false));
            testCase.verifyEqual(inits, [5, 5; 15, Inf]);
        end
    end

    methods (Static, Access = private)
        function block = block(plib, names, values)
            values = values(:);
            histories = repmat(values, 1, 3, 2, 3);
            inits = repmat(values, 1, 2);
            block = struct();
            block.plib = plib;
            block.plib_options = struct('variant', plib);
            block.problem_names = names;
            block.problem_dims = repmat(2, numel(names), 1);
            block.problem_types = repmat({'u'}, 1, numel(names));
            block.solver_names = {'a', 'b', 'c'};
            block.fun_histories = histories;
            block.maxcv_histories = zeros(size(histories));
            block.fun_outs = histories(:, :, :, end);
            block.maxcv_outs = zeros(size(block.fun_outs));
            block.fun_inits = inits;
            block.maxcv_inits = zeros(size(inits));
            block.merit_histories = histories;
            block.merit_outs = block.fun_outs;
            block.merit_inits = inits;
            block.n_evals = 3 * ones(size(block.fun_outs));
            block.computation_times = ones(size(block.fun_outs));
            block.solvers_successes = true(size(block.fun_outs));
            block.feature_stamp = 'noisy';
            block.ptype = 'u';
            block.mindim = 1;
            block.maxdim = 10;
            block.problem_names_options = {};
            block.excludelist = {};
            for field = {'mbs', 'mlcons', 'mnlcons', 'mcons'}
                block.(['problem_', field{1}]) = zeros(numel(names), 1);
            end
        end

        function block = paired(plib, value, plain_value)
            block = TestPlainReferenceIdentity.block(plib, {'DUP'}, value);
            block.results_plib_plain = TestPlainReferenceIdentity.block(plib, {'DUP'}, plain_value);
        end

        function mins = mins(blocks)
            [~, ~, ~, mins] = processResults(blocks, struct('run_plain', true));
            % processResults keeps a singleton solver axis; downstream callers
            % use two-dimensional indexing of its (problem, run) values.
            mins = reshape(mins, size(mins, 1), []);
        end

        function bytes = readBytes(file)
            fid = fopen(file, 'rb');
            cleanup = onCleanup(@() fclose(fid));
            bytes = fread(fid, Inf, '*uint8');
        end
    end
end
