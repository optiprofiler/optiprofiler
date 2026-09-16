classdef TestFeatureReviewRegressions < matlab.unittest.TestCase
% Regressions repaired after the 2026-09 independent review of the Feature
% and EvalReport candidate: duplicate option names, quantized runtime notes,
% load-invocation replay rejection, report diagnostics/metadata/staging and
% reservation rollback, relative and symlinked report paths, string-array
% options through the public benchmark, and restoring a 1.x FeaturedProblem.
    properties (Access = private)
        Work
    end
    methods (TestMethodSetup)
        function isolate(testCase)
            testCase.Work = tempname;
            mkdir(testCase.Work);
            original_directory = pwd;
            original_path = path;
            registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            testCase.addTeardown(@() cd(original_directory));
            testCase.addTeardown(@() path(original_path));
            testCase.addTeardown(@() setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', registry));
            testCase.addTeardown(@() removeDirectory(testCase.Work));
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(testCase.Work, 'registry.mat'));
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'composedstreams'));
        end
    end
    methods (Test)
        function duplicateCaseVariantOptionNamesAreRejected(testCase)
            duplicate = struct('noise_level', 0.1, 'NOISE_LEVEL', 0.2);
            testCase.verifyError(@() Feature('noisy', duplicate), 'MATLAB:Feature:DuplicateOption');
            testCase.verifyError(@() Feature('noisy', 'noise_level', 0.1, 'Noise_Level', 0.2), 'MATLAB:Feature:DuplicateOption');
            testCase.verifyError(@() Feature(struct('name', 'noisy', 'options', duplicate)), 'MATLAB:Feature:DuplicateOption');
            testCase.verifyError(@() Feature({'plain', struct('name', 'noisy', 'options', duplicate)}), 'MATLAB:Feature:DuplicateOption');
            message = '';
            try
                Feature('noisy', duplicate);
            catch cause
                message = lower(cause.message);
            end
            testCase.verifyTrue(contains(message, 'duplicate') && contains(message, 'ambiguous'), message);
            % One spelling of any case is still accepted and canonicalized.
            feature = Feature('noisy', struct('NOISE_LEVEL', 0.2));
            testCase.verifyEqual(feature.stages{1}.options.noise_level, 0.2);
            problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [1; 2]));
            options = struct('problem', problem, 'silent', true, 'score_only', true, ...
                'feature', {{struct('name', 'noisy', 'options', duplicate)}});
            testCase.verifyError(@() benchmark({@reviewStay, @reviewZero}, options), 'MATLAB:Feature:DuplicateOption');
        end

        function quantizedNotesComeFromStagesWithoutDeprecatedAccessors(testCase)
            registerFixtureLibrary();
            state = warning;
            restore = onCleanup(@() warning(state));
            warning('error', 'MATLAB:Feature:DeprecatedOptions');
            warning('error', 'MATLAB:Feature:DeprecatedModifier');
            options = struct('plibs', {{'composed_streams'}}, 'problem_names', {{'shift2'}}, 'silent', true, ...
                'draw_hist_plots', 'none', 'max_eval_factor', 3, 'max_tol_order', 1, 'n_jobs', 1, ...
                'savepath', testCase.Work, 'seed', 1, 'n_runs', 1);
            options.feature_name = 'quantized';
            options.benchmark_id = 'single';
            benchmark({@reviewStay, @reviewZero}, options);
            report = readReportText(fullfile(testCase.Work, 'single'));
            testCase.verifySubstring(report, '## Quantized truth');
            testCase.verifySubstring(report, 'ground_truth=1 scores the featured problem');
            options.feature_name = 'noisy+quantized';
            options.benchmark_id = 'composed';
            benchmark({@reviewStay, @reviewZero}, options);
            report = readReportText(fullfile(testCase.Work, 'composed'));
            testCase.verifySubstring(report, '## Quantized truth (stage 2, quantized#0)');
            testCase.verifySubstring(report, 'scoring reference (objective and nonlinear constraints) at the snapped point');
        end

        function loadInvocationOptionsAreNotReplayable(testCase)
            source = struct('schema', 'options_refined-v2', 'n_runs', 1, 'load', '20200101_000000', ...
                'feature_specification', {{struct('name', 'plain', 'options', struct())}});
            testCase.verifyError(@() loadBenchmarkOptions(source), 'OptiProfiler:LoadInvocationNotReplayable');
            user = struct('feature_name', 'noisy', 'load', 'latest', 'n_runs', 2);
            testCase.verifyError(@() loadBenchmarkOptions(user), 'OptiProfiler:LoadInvocationNotReplayable');
            options_refined = source;
            file = fullfile(testCase.Work, 'load-options.mat');
            save(file, 'options_refined', '-v7');
            testCase.verifyError(@() loadBenchmarkOptions(file), 'OptiProfiler:LoadInvocationNotReplayable');
            % An empty load selector is an absent one.
            source.load = '';
            options = loadBenchmarkOptions(source);
            testCase.verifyEqual(options.feature.name, 'plain');
            testCase.verifyFalse(isfield(options, 'load'));
        end

        function reportDiagnosticsAreDeduplicatedAndCappedWithACount(testCase)
            target = fullfile(testCase.Work, 'diagnostics.json');
            report = optiprofiler_internal.EvalReport(target, struct(), @reviewReplace);
            for k = 1:140
                report.addDiagnostic('review_probe', 'runtime', struct('index', k));
                report.addDiagnostic('review_probe', 'runtime', struct('index', k));
            end
            report.finish();
            document = jsondecode(fileread(target));
            testCase.verifyEqual(numel(document.diagnostics), 128);
            testCase.verifyEqual(document.diagnostics_omitted, 12);
            testCase.verifyEqual(numel(unique(arrayfun(@(d) d.scope.index, document.diagnostics))), 128);
        end

        function metadataMatricesKeepTheirShape(testCase)
            text = optiprofiler_internal.EvalReport.encodeMetadata(struct( ...
                'finite', [1 2; 3 4], 'with_nan', [1 NaN; 3 4], 'vector', [1 Inf]));
            testCase.verifyEqual(text, ['{"finite":[[1,2],[3,4]],"with_nan":[[1,{"value":null,"reason":"nan"}],[3,4]],', ...
                '"vector":[1,{"value":null,"reason":"positive_infinity"}]}']);
        end

        function reportStagingFilesAreOwnerOnlyBeforeContentIsWritten(testCase)
            if ~isunix, return; end
            setappdata(0, 'OP_REVIEW_STAGE_MODES', {});
            testCase.addTeardown(@() rmappdata(0, 'OP_REVIEW_STAGE_MODES'));
            target = fullfile(testCase.Work, 'staging.json');
            report = optiprofiler_internal.EvalReport(target, struct(), @reviewRecordingReplace);
            report.finish();
            modes = getappdata(0, 'OP_REVIEW_STAGE_MODES');
            testCase.verifyNotEmpty(modes);
            for k = 1:numel(modes)
                info = modes{k};
                testCase.verifyFalse(info.GroupRead || info.OtherRead || info.GroupWrite || info.OtherWrite, ...
                    'A staging file was readable or writable by others before publication.');
            end
        end

        function failedFirstWriteLeavesNoReservationsAndAllowsRetry(testCase)
            target = fullfile(testCase.Work, 'retry.json');
            companion = fullfile(testCase.Work, 'retry.plot_data.json');
            % The companion is published first; the main publish then fails.
            testCase.verifyError(@() optiprofiler_internal.EvalReport(target, struct(), @reviewFailMainReplace), ...
                'Review:InjectedReplaceFailure');
            testCase.verifyFalse(isfile(target), 'The failed main reservation was left behind.');
            testCase.verifyFalse(isfile(companion), 'The published companion of a failed constructor was left behind.');
            leftovers = dir(fullfile(testCase.Work, '*op-eval-report*'));
            testCase.verifyEmpty(leftovers, 'A staging file was left behind.');
            report = optiprofiler_internal.EvalReport(target, struct(), @reviewReplace); %#ok<NASGU>
            document = jsondecode(fileread(target));
            testCase.verifyEqual(document.status, 'running');
            % Negative control: a pre-existing companion is refused and preserved.
            other = fullfile(testCase.Work, 'kept.json');
            other_companion = fullfile(testCase.Work, 'kept.plot_data.json');
            fid = fopen(other_companion, 'w'); fwrite(fid, 'owner-sentinel', 'char'); fclose(fid);
            testCase.verifyError(@() optiprofiler_internal.EvalReport(other, struct(), @reviewReplace), ...
                'OptiProfiler:EvalReportExists');
            testCase.verifyEqual(fileread(other_companion), 'owner-sentinel');
            testCase.verifyFalse(isfile(other));
            % Negative control: a file another writer put in place of the main
            % target during the failed first write is preserved; only the
            % still-owned companion is removed.
            foreign = fullfile(testCase.Work, 'foreign.json');
            testCase.verifyError(@() optiprofiler_internal.EvalReport(foreign, struct(), @reviewForeignThenFail), ...
                'Review:InjectedReplaceFailure');
            testCase.verifyEqual(fileread(foreign), 'foreign writer');
            testCase.verifyFalse(isfile(fullfile(testCase.Work, 'foreign.plot_data.json')));
        end

        function foreignInPlaceRewriteAndInodeReuseAreDetected(testCase)
            % A foreign in-place rewrite keeps the inode, and after two foreign
            % replace-over-target operations ext4 hands the recorded inode
            % back; an identity made of the file key alone treated both as
            % owned (observed once in the EvalReport CI gate on syu-ubuntu).
            % Runs with and without the JVM (Java file key, or stat).
            target = fullfile(testCase.Work, 'inplace.json');
            report = optiprofiler_internal.EvalReport(target, struct(), @reviewReplace);
            fid = fopen(target, 'w'); fwrite(fid, 'external in-place rewrite', 'char'); fclose(fid);
            testCase.verifyError(@() report.finish(), 'OptiProfiler:EvalReportOwnership');
            testCase.verifyEqual(fileread(target), 'external in-place rewrite');
            owner = fullfile(testCase.Work, 'companion.json');
            companion = fullfile(testCase.Work, 'companion.plot_data.json');
            report = optiprofiler_internal.EvalReport(owner, struct(), @reviewReplace);
            fid = fopen(companion, 'w'); fwrite(fid, 'external companion rewrite', 'char'); fclose(fid);
            testCase.verifyError(@() report.finish(), 'OptiProfiler:EvalReportOwnership');
            testCase.verifyEqual(fileread(companion), 'external companion rewrite');
            replaced = fullfile(testCase.Work, 'replaced.json');
            report = optiprofiler_internal.EvalReport(replaced, struct(), @reviewReplace);
            for k = 1:2
                stage = [replaced, '.replacement'];
                fid = fopen(stage, 'w'); fwrite(fid, sprintf('external replacement %d', k), 'char'); fclose(fid);
                [ok, message] = movefile(stage, replaced, 'f'); assert(ok, message);
            end
            testCase.verifyError(@() report.finish(), 'OptiProfiler:EvalReportOwnership');
            testCase.verifyEqual(fileread(replaced), 'external replacement 2');
            % An untouched report still publishes afterwards.
            fresh = fullfile(testCase.Work, 'fresh.json');
            report = optiprofiler_internal.EvalReport(fresh, struct(), @reviewReplace);
            report.finish();
            document = jsondecode(fileread(fresh));
            testCase.verifyEqual(document.status, 'empty');
        end

        function relativeAndSymlinkedReportPathsResolveOnDisk(testCase)
            cd(testCase.Work);
            mkdir(fullfile('relative', 'reports'));
            mkdir(fullfile('relative', 'output'));
            report_path = fullfile('relative', 'reports', '..', 'reports', 'r.json');
            report = optiprofiler_internal.EvalReport(report_path, struct(), @reviewReplace);
            report.setOutputDirectory(fullfile('relative', 'output'));
            fid = fopen(fullfile('relative', 'output', 'README.txt'), 'w'); fprintf(fid, 'artifact'); fclose(fid);
            report.finish();
            data = jsondecode(fileread(report_path));
            testCase.verifyEqual(data.artifact_root, '../output');
            testCase.verifyEqual(data.artifacts(1).path, 'README.txt');
            testCase.verifyFalse(isfield(data.artifacts(1), 'path_reason'));
            actual = fullfile(fileparts(report_path), data.artifact_root, data.artifacts(1).path);
            testCase.verifyTrue(isfile(actual), 'The artifact reference must resolve through the filesystem.');
            testCase.verifyEqual(fileread(actual), 'artifact');
            testCase.verifyFalse(contains(data.artifacts(1).path, testCase.Work) || contains(data.artifact_root, testCase.Work), ...
                'Relative references must not expose the machine path.');
            if ~isunix, return; end
            % An output directory reached through a symlink of a different depth.
            mkdir(fullfile('physical', 'deep', 'out'));
            [status, ~] = system(sprintf('ln -s ''%s'' ''%s''', fullfile(testCase.Work, 'physical', 'deep', 'out'), ...
                fullfile(testCase.Work, 'alias')));
            if status ~= 0, return; end
            linked_report = fullfile(testCase.Work, 'alias_reports', 'r.json');
            linked = optiprofiler_internal.EvalReport(linked_report, struct(), @reviewReplace);
            linked.setOutputDirectory('alias');
            fid = fopen(fullfile('alias', 'data.txt'), 'w'); fprintf(fid, 'linked'); fclose(fid);
            linked.finish();
            data = jsondecode(fileread(linked_report));
            testCase.verifyEqual(data.artifacts(1).path, 'data.txt', 'Artifacts must not escape artifact_root through the link.');
            if usejava('jvm')
                % Java canonicalizes the (then empty) output directory; without
                % the JVM the link may keep its given name, which still resolves.
                testCase.verifyEqual(data.artifact_root, '../physical/deep/out');
            end
            resolved = fullfile(fileparts(linked_report), data.artifact_root, data.artifacts(1).path);
            testCase.verifyTrue(isfile(resolved));
            testCase.verifyEqual(fileread(resolved), 'linked');
        end

        function parentComponentAfterASymlinkFollowsTheFileSystem(testCase)
            % 'alias/../report.json' means '<target parent>/report.json' to the
            % operating system; folding the '..' lexically would silently write
            % to '<root>/report.json' instead (observed without the JVM).
            if ~isunix, return; end
            root = testCase.Work;
            target = fullfile(root, 'physical', 'deep', 'target');
            mkdir(target);
            [status, ~] = system(sprintf('ln -s ''%s'' ''%s''', target, fullfile(root, 'alias')));
            if status ~= 0, return; end
            requested = fullfile(root, 'alias', '..', 'report.json');
            report = optiprofiler_internal.EvalReport(requested, struct(), @reviewReplace);
            report.finish();
            testCase.verifyTrue(isfile(requested), 'The report must exist at the requested location.');
            testCase.verifyTrue(isfile(fullfile(root, 'physical', 'deep', 'report.json')), 'The OS location of alias/.. is the link target parent.');
            testCase.verifyFalse(isfile(fullfile(root, 'report.json')), 'The report must not be written to the lexical location.');
            % Missing final directories after the same construct are created under the physical parent.
            nested = fullfile(root, 'alias', '..', 'missing', 'sub', 'nested.json');
            report = optiprofiler_internal.EvalReport(nested, struct(), @reviewReplace);
            report.finish();
            testCase.verifyTrue(isfile(nested));
            testCase.verifyTrue(isfile(fullfile(root, 'physical', 'deep', 'missing', 'sub', 'nested.json')));
            testCase.verifyFalse(isfolder(fullfile(root, 'missing')));
            % An ordinary path with '.' and '..' components is unaffected.
            plain = fullfile(root, 'ordinary', '.', 'x', '..', 'plain.json');
            report = optiprofiler_internal.EvalReport(plain, struct(), @reviewReplace);
            report.finish();
            testCase.verifyTrue(isfile(fullfile(root, 'ordinary', 'plain.json')));
        end

        function stringArrayOptionsRunWithAReport(testCase)
            problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [1; 2]));
            options = struct('problem', problem, 'score_only', true, 'silent', true, 'n_runs', 1, ...
                'solver_names', {{'stay', 'zero'}}, 'excludelist', ["alpha", "beta"], ...
                'xlabel_performance_profile', ["Performance ratio", "(two-line label)"], ...
                'report_path', fullfile(testCase.Work, 'string-array.json'));
            scores = benchmark({@reviewStay, @reviewZero}, options);
            testCase.verifyTrue(all(isfinite(scores)));
            report = jsondecode(fileread(options.report_path));
            testCase.verifyEqual(report.status, 'completed');
            testCase.verifyEqual(report.configuration.request.excludelist, {'alpha'; 'beta'});
            testCase.verifyEqual(report.configuration.request.xlabel_performance_profile, {'Performance ratio'; '(two-line label)'});
            % The same request without a report still runs (legacy behaviour).
            legacy = rmfield(options, 'report_path');
            testCase.verifyEqual(benchmark({@reviewStay, @reviewZero}, legacy), scores);
        end

        function legacyFeaturedProblemRestoresAndContinuesItsStream(testCase)
            fixture_root = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'feature-v2');
            expected = load(fullfile(fixture_root, 'b0-featured-problem-continuation.mat'));
            expected = expected.continuation;
            lastwarn('');
            loaded = load(fullfile(fixture_root, 'b0-featured-problem.mat'));
            testCase.verifyEmpty(lastwarn(), 'Restoring the 1.x FeaturedProblem must not warn.');
            fp = loaded.fp_noisy;
            testCase.verifyClass(fp, 'FeaturedProblem');
            testCase.verifyClass(fp.feature, 'Feature');
            testCase.verifyEqual(fp.feature.name, 'noisy');
            testCase.verifyEqual(fp.execution_strategy, 'legacy-single');
            testCase.verifyEqual(fp.runtime_policy, 'matlab-legacy-single-v2');
            testCase.verifyEqual(fp.seed_policy, 'legacy-run-seed');
            testCase.verifyEqual(fp.n_eval_fun, 2);
            % The restored trial continues exactly where the 1.x object stopped.
            testCase.verifyEqual(fp.fun([0.25; 0.5]), expected.noisy_f3);
            testCase.verifyEqual(fp.cub([0.25; 0.5]), expected.noisy_c2);
            testCase.verifyEqual(fp.maxcv([0.25; 0.5]), expected.noisy_maxcv);
            testCase.verifyEqual(fp.fun_hist, expected.noisy_fun_hist);
            testCase.verifyEqual(fp.n_eval_fun, expected.noisy_n_eval_fun);
            plain = loaded.fp_plain;
            testCase.verifyEqual(plain.execution_strategy, 'identity');
            testCase.verifyEqual(plain.fun([0.25; 0.5]), expected.plain_f2);
            testCase.verifyEqual(plain.fun_hist, expected.plain_fun_hist);
        end
    end
end

function x = reviewStay(fun, x0)
    fun(x0);
    x = x0;
end

function x = reviewZero(fun, x0)
    x = zeros(size(x0));
    fun(x);
end

function reviewReplace(source, target)
    [ok, message] = movefile(source, target, 'f');
    assert(ok, message);
end

function reviewFailMainReplace(source, target)
    if endsWith(target, '.plot_data.json'), reviewReplace(source, target); return; end
    error('Review:InjectedReplaceFailure', 'Injected main-report publish failure.');
end

function reviewForeignThenFail(source, target)
    if endsWith(target, '.plot_data.json'), reviewReplace(source, target); return; end
    foreign = [target, '.foreign'];
    fid = fopen(foreign, 'w'); fwrite(fid, 'foreign writer', 'char'); fclose(fid);
    movefile(foreign, target, 'f');
    error('Review:InjectedReplaceFailure', 'Injected main-report publish failure after a foreign replacement.');
end

function reviewRecordingReplace(source, target)
    [~, info] = fileattrib(source);
    modes = getappdata(0, 'OP_REVIEW_STAGE_MODES');
    modes{end+1} = info;
    setappdata(0, 'OP_REVIEW_STAGE_MODES', modes);
    reviewReplace(source, target);
end

function text = readReportText(root)
    files = dir(fullfile(root, '**', 'report.txt'));
    assert(isscalar(files), 'Expected exactly one report.txt under %s, found %d.', root, numel(files));
    text = fileread(fullfile(files(1).folder, files(1).name));
end

function registerFixtureLibrary()
    fixture_root = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'composedstreams');
    registerProblemLibrary(struct('name', 'composed_streams', 'root', fixture_root, ...
        'select_function', 'composed_fixture_select', 'load_function', 'composed_fixture_load'));
end

function removeDirectory(directory)
    if isfolder(directory), rmdir(directory, 's'); end
end
