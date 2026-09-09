function run_matlab_eval_report_ci(repository_root, output_root)
%RUN_MATLAB_EVAL_REPORT_CI Run TestEvalReport on a platform with a strict gate.
% Unlike the general full-unit entry point, an Incomplete (filtered) case is a
% failure here: the mandatory platform cases (real workers, UTF-8 labels,
% owner-only permissions where applicable, hash receipts, either-target
% collision, load without solver calls) must actually execute on this OS and
% MATLAB release. The suite needs only the repository sources: no setup(),
% no persisted paths, and the registry/preferences come from the isolated
% environment variables exported by the workflow.
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
        assert(~isempty(ver('parallel')), 'OptiProfiler:ParallelToolboxRequired', ...
            'Parallel Computing Toolbox must be installed: the worker-equivalence case is mandatory here.');
        cd(fullfile(repository_root, 'matlab', 'optiprofiler'));
        suite = testsuite(fullfile(pwd, 'tests', 'unit_tests', 'TestEvalReport.m'));
        assert(~isempty(suite), 'OptiProfiler:EmptyCiSuite', 'TestEvalReport was not found.');
        runner = matlab.unittest.TestRunner.withTextOutput( ...
            'OutputDetail', matlab.unittest.Verbosity.Detailed);
        runner.addPlugin(matlab.unittest.plugins.XMLPlugin.producingJUnitFormat( ...
            fullfile(output_root, 'junit.xml')));
        stage = 'TestEvalReport';
        results = runner.run(suite);
        nfailed = nnz([results.Failed]);
        nincomplete = nnz([results.Incomplete]);
        appendLog(log_file, sprintf('Tests=%d Failed=%d Incomplete=%d', numel(results), nfailed, nincomplete));
        for k = 1:numel(results)
            appendLog(log_file, sprintf('  %s passed=%d failed=%d incomplete=%d duration=%.2fs', ...
                results(k).Name, results(k).Passed, results(k).Failed, results(k).Incomplete, results(k).Duration));
        end
        writeReceipt(fullfile(output_root, 'receipt.json'), results);
        assert(nfailed == 0, 'OptiProfiler:EvalReportTestsFailed', '%d test(s) failed.', nfailed);
        assert(nincomplete == 0, 'OptiProfiler:EvalReportTestsFiltered', ...
            '%d test(s) were filtered/incomplete; every EvalReport case is mandatory on this platform.', nincomplete);
        appendLog(log_file, 'TestEvalReport completed with Failed=0 and Incomplete=0.');
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

function writeReceipt(filename, results)
    names = {results.Name};
    receipt = struct('matlab', version, 'platform', computer, 'jvm', usejava('jvm'), ...
        'tests', numel(results), 'failed', nnz([results.Failed]), 'incomplete', nnz([results.Incomplete]), ...
        'names', {names}, 'durations', [results.Duration]);
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
