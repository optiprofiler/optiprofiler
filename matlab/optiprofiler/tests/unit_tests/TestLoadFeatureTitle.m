classdef TestLoadFeatureTitle < matlab.unittest.TestCase
%TESTLOADFEATURETITLE A load names the archived feature in the summary title, the README and the log.
%   The saved experiments are synthesized (see TestLoadReplot.makeExperiment), so no solver runs.

    properties (Access = private)
        CurrentDirectory
        WorkDir
        NameOf
    end

    methods (TestMethodSetup)

        function setupWorkDir(testCase)
            testCase.CurrentDirectory = pwd;
            source_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src');
            addpath(source_dir);
            cd(fullfile(source_dir, 'private'));
            testCase.NameOf = @archivedFeatureName;
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

        function archive(benchmark_id, feature_stamp)
            % Synthesize a saved experiment with the given feature stamp.
            TestLoadReplot.makeExperiment(pwd, benchmark_id, 2, 1);
            file = dir(fullfile(pwd, benchmark_id, '**', 'data_for_loading.mat'));
            file = fullfile(file.folder, file.name);
            saved = load(file, 'results_plibs');
            results_plibs = saved.results_plibs;
            results_plibs{1}.feature_stamp = feature_stamp;
            save(file, 'results_plibs');
        end

        function [names, scores, experiment] = namedFeatures(benchmark_id, varargin)
            % Load the latest experiment under benchmark_id into a new folder and return the feature
            % that the summary title, the README and the log of the load name.
            options = struct('load', 'latest', 'benchmark_id', benchmark_id, 'max_tol_order', 1, ...
                'draw_hist_plots', 'none', 'savepath', tempname(pwd));
            for i = 1:2:numel(varargin)
                options.(varargin{i}) = varargin{i + 1};
            end
            [printed, scores] = evalc('benchmark(options)');
            in_log = regexp(regexprep(printed, '\s+', ' '), 'Finished creating profiles with the "([^"]*)" feature\.', 'tokens', 'once');

            saved = dir(fullfile(options.savepath, '**', 'summary_*.fig'));
            assert(isscalar(saved), 'Expected one saved summary figure.');
            figure = openfig(fullfile(saved.folder, saved.name), 'invisible');
            guard = onCleanup(@() close(figure)); %#ok<NASGU>
            layouts = findall(figure, 'Type', 'tiledlayout');
            titles = arrayfun(@(layout) char(layout.Title.String), layouts, 'UniformOutput', false);
            titles = titles(startsWith(titles, 'Profiles with '));
            assert(isscalar(titles), 'Expected one summary title.');
            in_title = regexp(strrep(titles{1}, '\', ''), '^Profiles with ``(.*)'''' feature$', 'tokens', 'once');

            % The README is a table; a long file name wraps below the first column, beside the description.
            pdf = dir(fullfile(options.savepath, '**', 'summary_*.pdf'));
            assert(isscalar(pdf), 'Expected one summary PDF.');
            experiment = pdf.folder;
            rows = strsplit(fileread(fullfile(experiment, 'README.txt')), newline);
            column = strfind(rows{startsWith(rows, 'File/Folder Name')}, 'Description');
            description = '';
            for i = 1:numel(rows)
                if length(rows{i}) >= column
                    description = [description, ' ', rows{i}(column:end)]; %#ok<AGROW>
                end
            end
            in_readme = regexp(regexprep(description, '\s+', ' '), 'for the feature ''([^'']*)''\.', 'tokens', 'once');
            names = [in_title, in_readme, in_log];
        end

    end

    methods (Test)

        function loadWithoutFeatureNameNamesTheArchivedFeature(testCase)
            % A stamp set by the user is the only name that the saved experiment has.
            stamps = {'plain', 'linearly_transformed_rotated', 'my_noise'};
            names = {'plain', 'linearly_transformed', 'my_noise'};
            for i = 1:numel(stamps)
                benchmark_id = sprintf('stamp%d', i);
                TestLoadFeatureTitle.archive(benchmark_id, stamps{i});
                [named, ~, experiment] = TestLoadFeatureTitle.namedFeatures(benchmark_id);
                testCase.verifyEqual(named, repmat(names(i), 1, 3));
                % The new folder keeps the archived stamp, as before.
                testCase.verifySubstring(experiment, ['_u_2_3_', stamps{i}, '_']);
            end
        end

        function featureNameGivenToALoadStaysItsLabel(testCase)
            TestLoadFeatureTitle.archive('bench', 'linearly_transformed_rotated');
            archive = dir(fullfile(pwd, 'bench', '**', 'data_for_loading.mat'));
            bytes = readBytes(fullfile(archive.folder, archive.name));
            [~, scores] = TestLoadFeatureTitle.namedFeatures('bench');
            labels = {'noisy', 'plain'};
            for i = 1:numel(labels)
                [named, labelled_scores, experiment] = TestLoadFeatureTitle.namedFeatures('bench', 'feature_name', labels{i});
                testCase.verifyEqual(named, repmat(labels(i), 1, 3));
                % The label changes neither the results nor the source, nor the stamp of the new folder.
                testCase.verifyEqual(labelled_scores, scores);
                testCase.verifySubstring(experiment, '_u_2_3_linearly_transformed_rotated_');
            end
            testCase.verifyEqual(readBytes(fullfile(archive.folder, archive.name)), bytes);
        end

        function outputOfALoadKeepsTheArchivedName(testCase)
            % The output of a load is an archive of the same experiment.
            TestLoadFeatureTitle.archive('bench', 'linearly_transformed_rotated');
            [~, scores] = TestLoadFeatureTitle.namedFeatures('bench', 'savepath', fullfile(pwd, 'first'));
            [named, same_scores] = TestLoadFeatureTitle.namedFeatures('first');
            testCase.verifyEqual(named, repmat({'linearly_transformed'}, 1, 3));
            testCase.verifyEqual(same_scores, scores);
        end

        function nameOfAGeneratedStamp(testCase)
            % Stamps generated for every feature by both languages, and a composition of OptiProfiler 2.
            cases = {'plain', 'plain'; 'permuted', 'permuted'; 'custom', 'custom'; ...
                'nonquantifiable_constraints', 'nonquantifiable_constraints'; ...
                'perturbed_x0_0.001_spherical', 'perturbed_x0'; 'perturbed_x0_0.5_gaussian', 'perturbed_x0'; ...
                'noisy_0.001_mixed_gaussian', 'noisy'; 'noisy_0.001_mixed_deterministic_chebyshev', 'noisy'; ...
                'truncated_6', 'truncated'; 'truncated_4_perturbed_trailing_digits', 'truncated'; ...
                'linearly_transformed', 'linearly_transformed'; 'linearly_transformed_rotated', 'linearly_transformed'; ...
                'linearly_transformed_rotated_2', 'linearly_transformed'; 'linearly_transformed_cond1.5', 'linearly_transformed'; ...
                'random_nan_0.05', 'random_nan'; 'unrelaxable_constraints', 'unrelaxable_constraints'; ...
                'unrelaxable_constraints_bounds_linear_nonlinear', 'unrelaxable_constraints'; ...
                'quantized_0.001_ground_truth', 'quantized'; 'quantized_0.1', 'quantized'; ...
                'linearly_transformed_rotated__noisy_0.01_mixed_gaussian', 'linearly_transformed+noisy'};
            for i = 1:size(cases, 1)
                testCase.verifyEqual(testCase.NameOf(struct('feature_stamp', cases{i, 1})), cases{i, 2});
                testCase.verifyEqual(testCase.NameOf(struct('feature_stamp', string(cases{i, 1}))), cases{i, 2});
            end
        end

        function aStampOfAnotherFormIsTheName(testCase)
            stamps = {'my_label', 'plain_vs_noisy', 'custom_shift', 'permuted2', 'noisy2', 'Noisy_0.001', 'noisy__mine', ...
                '__noisy_0.001'};
            for i = 1:numel(stamps)
                testCase.verifyEqual(testCase.NameOf(struct('feature_stamp', stamps{i})), stamps{i});
            end
            testCase.verifyEmpty(testCase.NameOf(struct()));
            testCase.verifyEmpty(testCase.NameOf(struct('feature_stamp', '')));
            testCase.verifyEmpty(testCase.NameOf(struct('feature_stamp', 3)));
        end

    end

end

function bytes = readBytes(file)
    fid = fopen(file, 'rb');
    assert(fid >= 0);
    guard = onCleanup(@() fclose(fid)); %#ok<NASGU>
    bytes = fread(fid, Inf, '*uint8');
end
