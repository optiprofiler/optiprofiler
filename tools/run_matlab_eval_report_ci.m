function run_matlab_eval_report_ci(repository_root, output_root, parallel_case)
%RUN_MATLAB_EVAL_REPORT_CI Run TestEvalReport and TestProfileBands on a platform with a strict gate.
% Unlike the general full-unit entry point, an Incomplete (filtered) case is a
% failure here: the selected platform cases (UTF-8 labels, owner-only
% permissions where applicable, hash receipts, either-target collision, load
% without solver calls, and real workers unless omitted as below) must
% actually execute on this OS and MATLAB release, and the discovered, selected
% and reported case names must match the expected sets exactly. The suite
% needs only the repository sources: no setup(), no persisted paths, and the
% registry/preferences come from the isolated environment variables exported
% by the workflow.
%
% PARALLEL_CASE (default 'real-workers') keeps the full class. The only other
% value, 'parallel-omitted-batch-licensing', is accepted solely for the
% approved hosted CI partition: the workflow marks that environment with
% OPTIPROFILER_CI_HOSTED_BATCH_LICENSING=1 and the release must be exactly
% R2021b, the approved boundary (batch-licensed parallel pools are
% unsupported before R2023a, so the actualParallelEquivalence case cannot
% start workers there). The omitted case is listed in receipt.json with its
% reason; the real-worker rows (R2023a and latest) are separate required rows
% whose own receipts decide. Parallel computation itself is unchanged for
% normally licensed users.
    if nargin < 3 || isempty(parallel_case), parallel_case = 'real-workers'; end
    expected_cases = {'TestEvalReport/numericalIdentityAndUtf8', 'TestEvalReport/coverageAndNoOverwrite', ...
        'TestEvalReport/nonfiniteAndScoringFailure', 'TestEvalReport/compactExactExtremaAndCompanion', ...
        'TestEvalReport/paddedLengthAggregationBoundary', 'TestEvalReport/statefulRendererMeritObserved', ...
        'TestEvalReport/renderingFailureSeparate', 'TestEvalReport/archiveProvenanceAndArtifacts', ...
        'TestEvalReport/actualParallelEquivalence', 'TestEvalReport/writeFailureAndOwnership', ...
        'TestEvalReport/optionalPlainReferenceStatus'};
    parallel_name = 'TestEvalReport/actualParallelEquivalence';
    old_directory = pwd;
    old_path = path;
    cleanup = onCleanup(@() restoreState(old_directory, old_path));
    if ~isfolder(output_root), mkdir(output_root); end
    [ok, attributes] = fileattrib(output_root);
    assert(ok, 'OptiProfiler:CiOutputUnavailable', 'Cannot access the CI output directory.');
    output_root = attributes.Name;
    log_file = fullfile(output_root, 'eval-report-ci.log');
    appendLog(log_file, sprintf('MATLAB %s on %s (jvm=%d, parallel toolbox installed=%d)', ...
        version, computer, usejava('jvm'), ~isempty(ver('parallel'))));
    stage = 'test discovery';
    try
        omitted = struct('name', {}, 'reason', {});
        switch parallel_case
            case 'real-workers'
                assert(~isempty(ver('parallel')), 'OptiProfiler:ParallelToolboxRequired', ...
                    'Parallel Computing Toolbox must be installed: the worker-equivalence case is mandatory in this row.');
                % Classify the runner before the suite: a pool that cannot start
                % for a reason inside the toolbox itself is an infrastructure
                % limitation of this runner image, distinct from a product
                % failure. The suite still runs and its worker case still
                % fails, so the gate stays honest.
                pool_check = 'ok';
                try
                    pool = parpool(2);
                    delete(pool);
                catch pool_error
                    pool_check = sprintf('parpool failed [%s]: %s', pool_error.identifier, pool_error.message);
                end
                appendLog(log_file, ['parpool precheck: ', pool_check]);
                selected_cases = expected_cases;
            case 'parallel-omitted-batch-licensing'
                assert(strcmp(getenv('OPTIPROFILER_CI_HOSTED_BATCH_LICENSING'), '1'), 'OptiProfiler:ParallelOmissionNotApproved', ...
                    'The serial partition is approved only for the hosted batch-licensed CI environment (OPTIPROFILER_CI_HOSTED_BATCH_LICENSING=1).');
                assert(strcmp(version('-release'), '2021b'), 'OptiProfiler:ParallelOmissionNotAllowed', ...
                    'The parallel case may be omitted only on R2021b, the approved hosted-CI partition boundary (this is %s).', version('-release'));
                pool_check = 'not_requested_for_serial_partition';
                omitted(1) = struct('name', parallel_name, 'reason', ...
                    'batch_licensed_parallel_pools_are_unsupported_before_R2023a_on_hosted_runners; the real-worker case is a separate required row (R2023a and latest)');
                selected_cases = setdiff(expected_cases, {parallel_name}, 'stable');
                appendLog(log_file, ['case omitted in this row: ', omitted(1).name, ' (', omitted(1).reason, ')']);
            otherwise
                error('OptiProfiler:UnknownParallelCasePolicy', 'Unknown parallel_case value: %s', parallel_case);
        end
        cd(fullfile(repository_root, 'matlab', 'optiprofiler'));
        suite = testsuite(fullfile(pwd, 'tests', 'unit_tests', 'TestEvalReport.m'));
        assert(~isempty(suite), 'OptiProfiler:EmptyCiSuite', 'TestEvalReport was not found.');
        discovered = {suite.Name};
        assert(isequal(sort(discovered), sort(expected_cases)), 'OptiProfiler:UnexpectedCiSuite', ...
            'Discovered cases differ from the expected %d cases: %s', numel(expected_cases), strjoin(discovered, ', '));
        assert(nnz(strcmp(discovered, parallel_name)) == 1, 'OptiProfiler:ParallelCaseCount', 'The parallel case must be discovered exactly once.');
        suite = suite(ismember(discovered, selected_cases));
        assert(numel(suite) == numel(selected_cases) && isequal(sort({suite.Name}), sort(selected_cases)), ...
            'OptiProfiler:SelectionMismatch', 'Selected cases do not match the expected selection.');
        runner = matlab.unittest.TestRunner.withTextOutput( ...
            'OutputDetail', matlab.unittest.Verbosity.Detailed);
        runner.addPlugin(matlab.unittest.plugins.XMLPlugin.producingJUnitFormat( ...
            fullfile(output_root, 'junit.xml')));
        stage = 'TestEvalReport';
        results = runner.run(suite);
        nfailed = nnz([results.Failed]);
        nincomplete = nnz([results.Incomplete]);
        assert(isequal(sort({results.Name}), sort(selected_cases)), 'OptiProfiler:ResultSetMismatch', ...
            'Reported results (%d) do not match the selected cases (%d).', numel(results), numel(selected_cases));
        assert(all([results.Passed]) == (nfailed == 0 && nincomplete == 0), 'OptiProfiler:ResultFlagMismatch', 'Result flags are inconsistent.');
        appendLog(log_file, sprintf('Tests=%d Failed=%d Incomplete=%d', numel(results), nfailed, nincomplete));
        for k = 1:numel(results)
            appendLog(log_file, sprintf('  %s passed=%d failed=%d incomplete=%d duration=%.2fs', ...
                results(k).Name, results(k).Passed, results(k).Failed, results(k).Incomplete, results(k).Duration));
        end
        writeReceipt(fullfile(output_root, 'receipt.json'), results, pool_check, parallel_case, selected_cases, omitted);
        assert(nfailed == 0, 'OptiProfiler:EvalReportTestsFailed', '%d test(s) failed.', nfailed);
        assert(nincomplete == 0, 'OptiProfiler:EvalReportTestsFiltered', ...
            '%d test(s) were filtered/incomplete; every selected EvalReport case is mandatory on this platform.', nincomplete);
        assert(all([results.Passed]), 'OptiProfiler:EvalReportTestsNotPassed', 'Every selected case must be Passed.');
        appendLog(log_file, sprintf('TestEvalReport completed with Failed=0 and Incomplete=0 for %d selected of %d expected cases.', numel(selected_cases), numel(expected_cases)));
        % Second, serial gate on every row: the rendered profile-band regression (TestProfileBands runs its fixtures
        % with n_jobs = 1 and never requests a pool). Its exact case set, Failed=0, Incomplete=0 and all-Passed are
        % required on every platform and release; its own receipt and JUnit file are persisted before the gate.
        stage = 'TestProfileBands';
        band_cases = {'TestProfileBands/identicalRunsMeanStdPaintNothing', 'TestProfileBands/identicalRunsMinMaxPaintNothing', ...
            'TestProfileBands/singleRunDrawsNoBand', 'TestProfileBands/unequalRunsKeepTheirBand', ...
            'TestProfileBands/helperGeometryOnStairsArrays', 'TestProfileBands/adjacentFacesShowNoSeam'};
        band_suite = testsuite(fullfile(pwd, 'tests', 'unit_tests', 'TestProfileBands.m'));
        assert(isequal(sort({band_suite.Name}), sort(band_cases)), 'OptiProfiler:UnexpectedBandSuite', ...
            'Discovered TestProfileBands cases differ from the expected %d cases: %s', numel(band_cases), strjoin({band_suite.Name}, ', '));
        band_runner = matlab.unittest.TestRunner.withTextOutput('OutputDetail', matlab.unittest.Verbosity.Detailed);
        band_runner.addPlugin(matlab.unittest.plugins.XMLPlugin.producingJUnitFormat(fullfile(output_root, 'profile-bands-junit.xml')));
        band_results = band_runner.run(band_suite);
        appendLog(log_file, sprintf('TestProfileBands: Tests=%d Failed=%d Incomplete=%d', numel(band_results), nnz([band_results.Failed]), nnz([band_results.Incomplete])));
        for k = 1:numel(band_results)
            appendLog(log_file, sprintf('  %s passed=%d failed=%d incomplete=%d duration=%.2fs', ...
                band_results(k).Name, band_results(k).Passed, band_results(k).Failed, band_results(k).Incomplete, band_results(k).Duration));
        end
        writeBandReceipt(fullfile(output_root, 'profile-bands-receipt.json'), band_results, band_cases);
        assert(isequal(sort({band_results.Name}), sort(band_cases)), 'OptiProfiler:BandResultSetMismatch', 'TestProfileBands results do not match the expected cases.');
        assert(nnz([band_results.Failed]) == 0 && nnz([band_results.Incomplete]) == 0 && all([band_results.Passed]), ...
            'OptiProfiler:BandTestsNotPassed', 'TestProfileBands: %d failed, %d incomplete.', nnz([band_results.Failed]), nnz([band_results.Incomplete]));
        appendLog(log_file, sprintf('TestProfileBands completed with Failed=0 and Incomplete=0 for %d of %d cases (serial, n_jobs=1).', numel(band_results), numel(band_cases)));
    catch cause
        try
            appendLog(log_file, sprintf('Failed during %s [%s]:\n%s', stage, cause.identifier, ...
                getReport(cause, 'extended', 'hyperlinks', 'off')));
        catch write_error
            fprintf(2, 'Could not persist CI diagnostics: %s\n', write_error.message);
        end
        rethrow(cause);
    end
end

function writeReceipt(filename, results, pool_check, parallel_case, selected_cases, omitted)
    names = {results.Name};
    receipt = struct('matlab', version, 'platform', computer, 'jvm', usejava('jvm'), ...
        'tests', numel(results), 'failed', nnz([results.Failed]), 'incomplete', nnz([results.Incomplete]), ...
        'passed', nnz([results.Passed]), 'names', {names}, 'durations', [results.Duration], 'parpool_precheck', pool_check, ...
        'parallel_case', parallel_case, 'selected_cases', {selected_cases}, 'omitted_cases', omitted);
    fid = fopen(filename, 'w');
    assert(fid >= 0, 'OptiProfiler:CiReceiptUnavailable', 'Cannot write %s', filename);
    guard = onCleanup(@() fclose(fid));
    fwrite(fid, unicode2native(jsonencode(receipt), 'UTF-8'), 'uint8');
end

function writeBandReceipt(filename, band_results, band_cases)
% The rendered profile-band regression (serial, n_jobs = 1) is reported in its own receipt, persisted before any of
% its gates can fail, so the TestEvalReport receipt.json and its 10/11-case checks keep their exact meaning.
    receipt = struct('matlab', version, 'platform', computer, 'suite', 'TestProfileBands', 'serial', true, ...
        'expected_cases', {band_cases}, 'tests', numel(band_results), 'passed', nnz([band_results.Passed]), ...
        'failed', nnz([band_results.Failed]), 'incomplete', nnz([band_results.Incomplete]), ...
        'names', {{band_results.Name}}, 'durations', [band_results.Duration]);
    fid = fopen(filename, 'w');
    assert(fid >= 0, 'OptiProfiler:CiReceiptUnavailable', 'Cannot write %s', filename);
    guard = onCleanup(@() fclose(fid));
    fwrite(fid, unicode2native(jsonencode(receipt), 'UTF-8'), 'uint8');
end

function appendLog(filename, message)
    fid = fopen(filename, 'a');
    assert(fid >= 0, 'OptiProfiler:CiLogUnavailable', 'Cannot write CI log: %s', filename);
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', message);
    fprintf('%s\n', message);
end

function restoreState(old_directory, old_path)
    cd(old_directory);
    path(old_path);
end
