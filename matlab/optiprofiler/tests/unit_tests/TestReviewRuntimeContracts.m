classdef TestReviewRuntimeContracts < matlab.unittest.TestCase
    properties (Access = private)
        ValidateOptions
    end
    methods (TestMethodSetup)
        function locateValidator(testCase)
            original = pwd;
            cleanup = onCleanup(@() cd(original)); %#ok<NASGU>
            cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'private'));
            testCase.ValidateOptions = @checkValidityProfileOptions;
        end
    end
    methods (Test)
        function budgetMustBeFinite(testCase)
            for value = [NaN, Inf, -Inf, 0, -1]
                testCase.verifyError(@() testCase.ValidateOptions({}, ...
                    struct('max_eval_factor', value)), ...
                    'MATLAB:checkValidityProfileOptions:max_eval_factorNotValid');
            end
            options = testCase.ValidateOptions({}, struct('max_eval_factor', 2));
            testCase.verifyEqual(options.max_eval_factor, 2);
        end
        function quantizedViolationNormalizesVectors(testCase)
            problems = {struct(), struct('xl', [0;1], 'xu', [2;3]), ...
                struct('cub', @(x) ones(1,1))};
            for k = 1:numel(problems)
                options = problems{k};
                options.fun = @(x) sum(x);
                options.x0 = [0;0];
                fp = FeaturedProblem(Problem(options), Feature('quantized'), 20, 0);
                testCase.verifyEqual(fp.maxcv([1,2]), fp.maxcv([1;2]));
                testCase.verifyError(@() fp.maxcv(1), ...
                    'MATLAB:FeaturedProblem:maxcvInvalidPoint');
            end
        end
        function composedViolationNormalizesVectors(testCase)
            p = Problem(struct('fun', @(x) sum(x), 'x0', [0;0], 'cub', @(x) x(1)));
            fp = FeaturedProblem(p, Feature('truncated+truncated'), 20, 0);
            testCase.verifyEqual(fp.maxcv([1,2]), fp.maxcv([1;2]));
        end
    end
end
