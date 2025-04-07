classdef TestCheckValidityCutestOptions < matlab.unittest.TestCase
    methods (Test)

        function testErrors(testCase)

            options = struct();

            options.ptype = 1;
            testCase.verifyError(@() checkValidityCutestOptions(options), "MATLAB:checkValidityCutestOptions:ptypeNotcharstr");

            options.ptype = 'a';
            testCase.verifyError(@() checkValidityCutestOptions(options), "MATLAB:checkValidityCutestOptions:ptypeNotubln");
            options = rmfield(options, 'ptype');

            options.mindim = 0;
            testCase.verifyError(@() checkValidityCutestOptions(options), "MATLAB:checkValidityCutestOptions:mindimNotValid");
            options = rmfield(options, 'mindim');

            options.maxdim = 0;
            testCase.verifyError(@() checkValidityCutestOptions(options), "MATLAB:checkValidityCutestOptions:maxdimNotValid");

            options.mindim = 2;
            options.maxdim = 1;
            testCase.verifyError(@() checkValidityCutestOptions(options), "MATLAB:checkValidityCutestOptions:maxdimSmallerThanmindim");
            options = rmfield(options, 'mindim');
            options = rmfield(options, 'maxdim');

            options.mincon = -1;
            testCase.verifyError(@() checkValidityCutestOptions(options), "MATLAB:checkValidityCutestOptions:minconNotValid");
            options.mincon = 2;

            options.maxcon = -1;
            testCase.verifyError(@() checkValidityCutestOptions(options), "MATLAB:checkValidityCutestOptions:maxconNotValid");

            options.maxcon = 1;
            testCase.verifyError(@() checkValidityCutestOptions(options), "MATLAB:checkValidityCutestOptions:maxconSmallerThanmincon");
            options = rmfield(options, 'mincon');
            options = rmfield(options, 'maxcon');

            options.excludelist = 1;
            testCase.verifyError(@() checkValidityCutestOptions(options), "MATLAB:checkValidityCutestOptions:excludelistNotCellOfcharstr");
        end

    end

end