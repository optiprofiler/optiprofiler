% Copyright 2021 The MathWorks, Inc.

import matlab.unittest.TestRunner
import matlab.unittest.Verbosity
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoberturaFormat
 
% Add paths
addpath(genpath(pwd));

% Create a test suite 
suite = testsuite(pwd, 'IncludeSubfolders', true);

% Create a test runner that displays test run progress at the matlab.unittest.Verbosity.Detailed level
runner = TestRunner.withTextOutput('OutputDetail',Verbosity.Detailed); 

% Create a CodeCoveragePlugin instance and add it to the test runner
sourceFolder = 'src';
reportFile = 'coverage.xml';
reportFormat = CoberturaFormat(reportFile);
p = CodeCoveragePlugin.forFolder(sourceFolder,'IncludingSubfolders', true,'Producing',reportFormat);
runner.addPlugin(p)
 
% Run the tests and fail the build if any of the tests fails
results = runner.run(suite);  
nfailed = nnz([results.Failed]);
assert(nfailed == 0,[num2str(nfailed) ' test(s) failed.'])

% Remove paths
rmpath(genpath(pwd));