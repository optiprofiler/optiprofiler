function run_matlab_full_unit_ci(repository_root, output_root)
%RUN_MATLAB_FULL_UNIT_CI Run the ordinary unit suite with durable CI receipts.
% Keep suite selection, verbosity, coverage scope, and the Failed gate aligned
% with matlab/optiprofiler/runTests.m. This CI-only entry does not change that
% script's workspace interface or put receipts inside cleaned test fixtures.
    old_directory = pwd;
    old_path = path;
    cleanup = onCleanup(@() restoreState(old_directory, old_path));
    if ~isfolder(output_root), mkdir(output_root); end
    [ok, attributes] = fileattrib(output_root);
    assert(ok, 'OptiProfiler:CiOutputUnavailable', 'Cannot access CI output directory.');
    output_root = attributes.Name;
    log_file = fullfile(output_root, 'ci.log');
    appendLog(log_file, 'Starting MATLAB full-unit CI.');
    stage = 'PDF prerequisite checks';
    try
        % Check from MATLAB, not just the runner shell: its children can have
        % different loader paths. The workflow's child-only wrappers stay used.
        for tool = {'pdfunite', 'pdftotext'}
            [status, diagnostic] = system([tool{1}, ' -v 2>&1']);
            appendLog(log_file, diagnostic);
            assert(status == 0, 'OptiProfiler:MissingPdfTestTool', ...
                '%s is required for PDF regression tests: %s', tool{1}, diagnostic);
        end
        stage = 'setup';
        cd(repository_root);
        setup(struct('install_matcutest', false, 'install_solar', false));

        stage = 'test discovery and runner configuration';
        cd(fullfile(repository_root, 'matlab', 'optiprofiler'));
        suite = testsuite(fullfile(pwd, 'tests', 'unit_tests'), 'IncludeSubfolders', true);
        assert(~isempty(suite), 'OptiProfiler:EmptyCiSuite', 'The full unit suite is empty.');
        runner = matlab.unittest.TestRunner.withTextOutput( ...
            'OutputDetail', matlab.unittest.Verbosity.Detailed);
        runner.addPlugin(matlab.unittest.plugins.XMLPlugin.producingJUnitFormat( ...
            fullfile(output_root, 'junit.xml')));
        runner.addPlugin(matlab.unittest.plugins.CodeCoveragePlugin.forFolder( ...
            fullfile(pwd, 'src'), 'IncludingSubfolders', true, 'Producing', ...
            matlab.unittest.plugins.codecoverage.CoberturaFormat(fullfile(output_root, 'coverage.xml'))));
        % benchmark can change or disable diary. Framework streams preserve
        % live progress/diagnostics independently, without buffering a full run.
        stream = matlab.unittest.plugins.ToFile(fullfile(output_root, 'runner.log'));
        runner.addPlugin(matlab.unittest.plugins.TestRunProgressPlugin.withVerbosity( ...
            matlab.unittest.Verbosity.Detailed, stream));
        runner.addPlugin(matlab.unittest.plugins.DiagnosticsOutputPlugin(stream));

        stage = 'unit tests';
        results = runner.run(suite);
        nfailed = nnz([results.Failed]);
        appendLog(log_file, sprintf('Tests=%d Failed=%d Incomplete=%d', ...
            numel(results), nfailed, nnz([results.Incomplete])));
        % JUnit records actual names, failures and assumptions, including when
        % the following unchanged Failed gate rejects the completed run.
        for artifact = {'junit.xml', 'coverage.xml', 'runner.log'}
            info = dir(fullfile(output_root, artifact{1}));
            assert(~isempty(info) && info.bytes > 0, 'OptiProfiler:MissingCiArtifact', ...
                'Expected nonempty CI artifact: %s', artifact{1});
        end
        assert(nfailed == 0, 'OptiProfiler:UnitTestsFailed', '%d test(s) failed.', nfailed);
        appendLog(log_file, 'Full unit suite completed without failed tests.');
    catch cause
        % Preserve the primary setup/test error even if diagnostics cannot be
        % written (for example, a full disk). Upload runs independently in CI.
        try
            appendLog(log_file, sprintf('Failed during %s [%s]:\n%s', stage, cause.identifier, ...
                getReport(cause, 'extended', 'hyperlinks', 'off')));
        catch write_error
            fprintf(2, 'Could not persist CI diagnostics: %s\n', write_error.message);
        end
        rethrow(cause);
    end
end

function appendLog(filename, message)
    fid = fopen(filename, 'a');
    assert(fid >= 0, 'OptiProfiler:CiLogUnavailable', 'Cannot write CI log: %s', filename);
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', message);
end

function restoreState(directory, original_path)
    cd(directory);
    path(original_path);
end
