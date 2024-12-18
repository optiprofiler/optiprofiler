classdef TestCheckValidityOtherOptions < matlab.unittest.TestCase
    methods (Test)

        function testErrors(testCase)

            solvers = {@fminsearch, @fminunc};
            options = struct();

            options.solver_names = 1;
            testCase.verifyError(@() checkValidityOtherOptions(solvers, options), "MATLAB:checkValidityOtherOptions:solver_namesNotCellOfcharstr");
            options = rmfield(options, 'solver_names');

            options.solver_names = {'a', 'b', 'c'};
            testCase.verifyError(@() checkValidityOtherOptions(solvers, options), "MATLAB:checkValidityOtherOptions:solver_namesAndsolversLengthNotSame");
            options = rmfield(options, 'solver_names');

            options.solver_isrand = 2;
            testCase.verifyError(@() checkValidityOtherOptions(solvers, options), "MATLAB:checkValidityOtherOptions:solver_israndNotLogical");
            options = rmfield(options, 'solver_isrand');

            options.problem = 1;
            testCase.verifyError(@() checkValidityOtherOptions(solvers, options), "MATLAB:checkValidityOtherOptions:problemNotProblem");
            options = rmfield(options, 'problem');

            options.cutest_problem_names = 1;
            testCase.verifyError(@() checkValidityOtherOptions(solvers, options), "MATLAB:checkValidityOtherOptions:cutest_problem_namesNotValid");
            options = rmfield(options, 'cutest_problem_names');

            options.custom_problem_loader = 1;
            testCase.verifyError(@() checkValidityOtherOptions(solvers, options), "MATLAB:checkValidityOtherOptions:customloaderNotFunctionHandle");
            options = rmfield(options, 'custom_problem_loader');

            options.custom_problem_names = 1;
            testCase.verifyError(@() checkValidityOtherOptions(solvers, options), "MATLAB:checkValidityOtherOptions:customnamesNotcharstrOrCellOfcharstr");
            options = rmfield(options, 'custom_problem_names');

            options.custom_problem_loader = @(x) x;
            options.custom_problem_names = {'A', 'B'};
            testCase.verifyError(@() checkValidityOtherOptions(solvers, options), "MATLAB:checkValidityOtherOptions:customloaderNotAcceptcustomnames");
        end
    end

end