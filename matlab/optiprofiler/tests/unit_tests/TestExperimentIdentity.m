classdef TestExperimentIdentity < matlab.unittest.TestCase
%TESTEXPERIMENTIDENTITY Preserve earlier experiments and load exact run IDs.
    properties
        WorkDir
        OriginalDirectory
    end

    methods (TestMethodSetup)
        function isolate(testCase)
            testCase.OriginalDirectory = pwd;
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
            testCase.WorkDir = tempname;
            mkdir(testCase.WorkDir);
            cd(testCase.WorkDir);
        end
    end

    methods (TestMethodTeardown)
        function restore(testCase)
            cd(testCase.OriginalDirectory);
            rmdir(testCase.WorkDir, 's');
        end
    end

    methods (Static)
        function x = identitySolver(fun, x)
            fun(x);
        end

        function options = options(varargin)
            options = struct('load', '20200101_000000', 'benchmark_id', 'bench', ...
                'max_tol_order', 1, 'silent', true, 'draw_hist_plots', 'none', ...
                'summarize_performance_profiles', false, 'summarize_data_profiles', false);
            for k = 1:2:numel(varargin)
                options.(varargin{k}) = varargin{k + 1};
            end
        end
    end

    methods (Test)
        function timestampParserPreservesNumericSuffix(testCase)
            % Keep a focused guard for MATLAB regexp's optional capture semantics.
            old_directory = pwd;
            cleanup = onCleanup(@() cd(old_directory));
            cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'private'));
            base = parseExperimentTimeStamp('20200101_000000');
            testCase.verifyEqual(base(2), 0);
            for sequence = [1, 999, 1000]
                key = parseExperimentTimeStamp(sprintf('20200101_000000_%03d', sequence));
                testCase.verifyEqual(key, [base(1), sequence]);
            end
            for invalid = {'20200101_000000_000', '20200101_000000_00', ...
                    '20200101_000000_001_extra', 'not_an_id'}
                testCase.verifyTrue(all(isnan(parseExperimentTimeStamp(invalid{1}))));
            end
        end

        function latestOrdersSuffixesNumericallyAndRejectsZero(testCase)
            % Deliberately reverse the path order relative to numeric sequence.
            TestExperimentIdentity.makeSequenceExperiment('z_999', '20200101_000000_999', 2);
            TestExperimentIdentity.makeSequenceExperiment('a_1000', '20200101_000000_1000', 3);
            options = TestExperimentIdentity.options('load', 'latest', 'score_only', true);
            testCase.verifySize(benchmark(options), [3, 1]);
            TestExperimentIdentity.makeSequenceExperiment('invalid_zero', '20990101_000000_000', 4);
            testCase.verifySize(benchmark(options), [3, 1]);
            testCase.verifyError(@() benchmark(TestExperimentIdentity.options( ...
                'load', '20990101_000000_000', 'score_only', true)), ...
                'MATLAB:checkValidityProfileOptions:loadNotValid');
        end

        function aggregateSummaryOrdersSuffixesNumericallyAndRejectsZero(testCase)
            [status, ~] = system('pdftotext -v 2>&1');
            testCase.assumeEqual(status, 0, 'pdftotext is required to inspect aggregate page order.');
            TestLoadReplot.makeExperiment(pwd, 'bench', 2, 1);
            folder = fullfile('bench', 'summary_inputs');
            mkdir(folder);
            TestExperimentIdentity.writeLabelPdf(fullfile(folder, ...
                'summary_a_20990101_000000_999.pdf'), 'SUFFIX999');
            TestExperimentIdentity.writeLabelPdf(fullfile(folder, ...
                'summary_b_20990101_000000_1000.pdf'), 'SUFFIX1000');
            TestExperimentIdentity.writeLabelPdf(fullfile(folder, ...
                'summary_c_20990101_000000_000.pdf'), 'INVALIDZERO');
            benchmark(TestExperimentIdentity.options());
            [status, text] = system(sprintf('pdftotext "%s" -', fullfile(pwd, 'bench', 'summary.pdf')));
            testCase.assertEqual(status, 0, text);
            first = strfind(text, 'SUFFIX1000');
            second = strfind(text, 'SUFFIX999');
            testCase.assertNumElements(first, 1);
            testCase.assertNumElements(second, 1);
            testCase.verifyLessThan(first, second);
            testCase.verifyFalse(contains(text, 'INVALIDZERO'));
            testCase.verifyEmpty(dir(fullfile('bench', '*', '.summary_merge')));
        end

        function existingSameSecondOutputSurvives(testCase)
            % Cover the current second and nearby boundaries without changing
            % MATLAB's datetime implementation or the solver random stream.
            start_time = datetime('now');
            sentinels = cell(1, 63);
            for k = -2:60
                stamp = char(datetime(start_time + seconds(k), 'Format', 'yyyyMMdd_HHmmss'));
                directory = fullfile('collision', ['a_b_u_1_2_plain_', stamp]);
                mkdir(directory);
                sentinels{k + 3} = fullfile(directory, 'keep.txt');
                fid = fopen(sentinels{k + 3}, 'w');
                fprintf(fid, 'original');
                fclose(fid);
            end
            options = struct('problem', Problem(struct('fun', @(x) x(1)^2, 'x0', 1)), ...
                'solver_names', {{'a', 'b'}}, 'benchmark_id', 'collision', ...
                'draw_hist_plots', 'none', 'silent', true);
            benchmark({@TestExperimentIdentity.identitySolver, @TestExperimentIdentity.identitySolver}, options);
            testCase.verifyTrue(all(cellfun(@isfile, sentinels)), 'An earlier experiment was deleted.');
            directories = dir('collision');
            testCase.verifyEqual(sum([directories.isdir]) - 2, 64);
        end

        function oldAndSuffixedIdsLoadWithoutChangingSource(testCase)
            TestLoadReplot.makeExperiment(pwd, 'bench', 2, 1);
            original = dir(fullfile('bench', '*', 'test_log', 'data_for_loading.mat'));
            original_file = fullfile(original.folder, original.name);
            original_bytes = TestExperimentIdentity.readBytes(original_file);
            suffixed = fullfile('bench', 'copy_20200101_000000_001', 'test_log');
            mkdir(suffixed);
            copyfile(original_file, fullfile(suffixed, 'data_for_loading.mat'));
            fid = fopen(fullfile(suffixed, 'time_stamp_20200101_000000_001.txt'), 'w');
            fclose(fid);
            score_old = benchmark(TestExperimentIdentity.options());
            score_suffix = benchmark(TestExperimentIdentity.options('load', '20200101_000000_001'));
            testCase.verifyEqual(score_suffix, score_old);
            testCase.verifyEqual(TestExperimentIdentity.readBytes(original_file), original_bytes);
            testCase.verifyEqual(numel(dir(fullfile('bench', '*', 'test_log', 'data_for_loading.mat'))), 4);
        end

        function latestSkipsIncompleteAndMalformedMarkers(testCase)
            TestLoadReplot.makeExperiment(pwd, 'bench', 2, 1);
            incomplete = fullfile('bench', 'incomplete', 'test_log');
            mkdir(incomplete);
            fid = fopen(fullfile(incomplete, 'time_stamp_20990101_000000.txt'), 'w');
            fclose(fid);
            fid = fopen(fullfile(incomplete, 'time_stamp_not_an_id.txt'), 'w');
            fclose(fid);
            scores = benchmark(TestExperimentIdentity.options('load', 'latest', 'score_only', true));
            testCase.verifySize(scores, [2, 1]);
        end

        function duplicateExactIdRequiresNarrowerDirectory(testCase)
            TestLoadReplot.makeExperiment(pwd, 'bench', 2, 1);
            original = dir(fullfile('bench', '*', 'test_log', 'data_for_loading.mat'));
            duplicate = fullfile('bench', 'duplicate', 'test_log');
            copyfile(original.folder, duplicate);
            testCase.verifyError(@() benchmark(TestExperimentIdentity.options('score_only', true)), ...
                'OptiProfiler:AmbiguousTimeStamp');
            cd(fileparts(original.folder));
            options = TestExperimentIdentity.options('score_only', true, 'benchmark_id', '.');
            testCase.verifySize(benchmark(options), [2, 1]);
        end

        function aggregateSummaryAcceptsSuffixedId(testCase)
            TestLoadReplot.makeExperiment(pwd, 'bench', 2, 1);
            options = TestExperimentIdentity.options('summarize_performance_profiles', true);
            benchmark(options);
            summaries = dir(fullfile('bench', '*', 'summary_*.pdf'));
            testCase.assertNotEmpty(summaries);
            aggregate_before = TestExperimentIdentity.readBytes(fullfile('bench', 'summary.pdf'));
            source = fullfile(summaries(1).folder, summaries(1).name);
            copyfile(source, [source(1:end-4), '_001.pdf']);
            benchmark(options);
            testCase.verifyTrue(isfile(fullfile('bench', 'summary.pdf')));
            testCase.verifyNotEqual(TestExperimentIdentity.readBytes(fullfile('bench', 'summary.pdf')), aggregate_before);
            testCase.verifyEmpty(dir(fullfile('bench', '*', '.summary_merge')));
        end

        function failedAggregateMergePreservesPreviousPdf(testCase)
            TestLoadReplot.makeExperiment(pwd, 'bench', 2, 1);
            options = TestExperimentIdentity.options('summarize_performance_profiles', true);
            benchmark(options);
            aggregate_file = fullfile('bench', 'summary.pdf');
            original_bytes = TestExperimentIdentity.readBytes(aggregate_file);
            mkdir(fullfile('bench', 'damaged'));
            fid = fopen(fullfile('bench', 'damaged', 'summary_broken_20990101_000000.pdf'), 'w');
            fprintf(fid, 'This is not a valid PDF.');
            fclose(fid);
            benchmark(options);
            testCase.verifyEqual(TestExperimentIdentity.readBytes(aggregate_file), original_bytes);
            testCase.verifyEmpty(dir(fullfile('bench', '*', '.summary_merge')));
        end
    end

    methods (Static, Access = private)
        function makeSequenceExperiment(folder, stamp, n_solvers)
            root = fullfile('bench', folder);
            TestLoadReplot.makeExperiment(pwd, root, n_solvers, 1);
            marker = dir(fullfile(root, '*', 'test_log', 'time_stamp_*.txt'));
            movefile(fullfile(marker.folder, marker.name), ...
                fullfile(marker.folder, ['time_stamp_', stamp, '.txt']));
        end

        function writeLabelPdf(file, label)
            fig = figure('Visible', 'off');
            cleanup = onCleanup(@() close(fig));
            axes('Parent', fig, 'Visible', 'off');
            text(0.1, 0.5, label, 'Interpreter', 'none');
            if ispc
                print(fig, file, '-dpdf', '-vector');
            else
                exportgraphics(fig, file, 'ContentType', 'vector');
            end
        end

        function bytes = readBytes(file)
            fid = fopen(file, 'rb');
            cleanup = onCleanup(@() fclose(fid));
            bytes = fread(fid, Inf, '*uint8');
        end
    end
end
