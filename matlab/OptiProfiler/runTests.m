% Copyright 2021 The MathWorks, Inc.

import matlab.unittest.TestRunner
import matlab.unittest.Verbosity
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoberturaFormat
 
% Add the source folder to the MATLAB search path 
addpath('src') 

% Create a test suite 
runner = matlab.unittest.TestRunner.withTextOutput;

% Create a test runner
runner = matlab.unittest.TestRunner.withTextOutput;

% Create a CodeCoveragePlugin instance and add it to the test runner
sourceFolder = './src';
reportFile = 'coverage.xml';
reportFormat = CoberturaFormat(reportFile);
p = CodeCoveragePlugin.forFolder(sourceFolder,'IncludingSubfolders', true,'Producing',reportFormat);
runner.addPlugin(p)
 
% Run the tests and fail the build if any of the tests fails
results = runner.run(suite);  
nfailed = nnz([results.Failed]);
assert(nfailed == 0,[num2str(nfailed) ' test(s) failed.'])