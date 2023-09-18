function runTestsAndGenerateCoverageReport()

    clc
    addpath('../src');

    import matlab.unittest.plugins.CodeCoveragePlugin
    import matlab.unittest.plugins.codecoverage.CoverageReport

    runner = matlab.unittest.TestRunner.withTextOutput;
    sourceCodeFolder = '../src';
    reportFolder = './coverageReport';

    reportFormat = CoverageReport(reportFolder);
    p = CodeCoveragePlugin.forFolder(sourceCodeFolder, 'Producing', reportFormat);
    runner.addPlugin(p);

    suite = matlab.unittest.TestSuite.fromFolder('.');
    results = runner.run(suite);

    rmpath('../src');
end