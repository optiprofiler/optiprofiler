classdef TestCheckValidityCutestOptions < matlab.unittest.TestCase
    methods (Test)

        function testErrors(testCase)

            options = struct();

            options.p_type = 1;
            testCase.verifyError(@() checkValidityCutestOptions(options), "MATLAB:checkValidityCutestOptions:p_typeNotcharstr");

            options.p_type = 'a';
            testCase.verifyError(@() checkValidityCutestOptions(options), "MATLAB:checkValidityCutestOptions:p_typeNotubln");
            options = rmfield(options, 'p_type');

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