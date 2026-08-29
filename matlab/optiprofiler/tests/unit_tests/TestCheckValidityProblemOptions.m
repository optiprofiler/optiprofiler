classdef TestCheckValidityProblemOptions < matlab.unittest.TestCase
    methods (Test)

        function testInfinityBounds(testCase)

            mydir = fileparts(mfilename('fullpath'));
            original_dir = pwd;
            testCase.addTeardown(@() cd(original_dir));
            cd(fullfile(mydir, '../../src/private'));

            profile_options = struct();
            upper_fields = {'maxdim', 'maxb', 'maxlcon', 'maxnlcon', 'maxcon'};
            lower_fields = {'mindim', 'minb', 'minlcon', 'minnlcon', 'mincon'};

            for i = 1:numel(upper_fields)
                upper = struct(upper_fields{i}, Inf);
                checked = checkValidityProblemOptions(upper, profile_options);
                testCase.verifyEqual(checked.(upper_fields{i}), Inf);

                negative_upper = struct(upper_fields{i}, -Inf);
                testCase.verifyError(@() checkValidityProblemOptions(negative_upper, profile_options), ...
                    ['MATLAB:checkValidityProblemOptions:' upper_fields{i} 'NotValid']);
            end

            for i = 1:numel(lower_fields)
                lower = struct(lower_fields{i}, Inf);
                testCase.verifyError(@() checkValidityProblemOptions(lower, profile_options), ...
                    ['MATLAB:checkValidityProblemOptions:' lower_fields{i} 'NotValid']);
            end
        end

        function testErrors(testCase)

            mydir = fileparts(mfilename('fullpath'));
            cd(fullfile(mydir, '../../src/private'));

            options = struct();
            profile_options = struct();

            options.plibs = 1;
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:plibsNotValid");
            if ~isunix || ismac
                mydir = fileparts(mfilename('fullpath'));
                problem_libs_dir = fullfile(mydir, '../../problem_libs');
                matcutest_dir = fullfile(problem_libs_dir, 'matcutest');
                created_matcutest_dir = exist(matcutest_dir, 'dir') ~= 7;
                if created_matcutest_dir
                    mkdir(matcutest_dir);
                    testCase.addTeardown(@() rmdir(matcutest_dir));
                end
                options.plibs = {'matcutest'};
                testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:plibsNotLinux");
            end
            options = rmfield(options, 'plibs');


            options.ptype = 1;
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:ptypeNotcharstr");
            options.ptype = 'a';
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:ptypeNotubln");
            options = rmfield(options, 'ptype');


            options.mindim = 0;
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:mindimNotValid");
            options = rmfield(options, 'mindim');
            options.maxdim = 0;
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:maxdimNotValid");
            options.mindim = 2;
            options.maxdim = 1;
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:maxdimSmallerThanmindim");
            options = rmfield(options, 'mindim');
            options = rmfield(options, 'maxdim');


            options.minb = -1;
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:minbNotValid");
            options.minb = 2;
            options.maxb = -1;
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:maxbNotValid");
            options.maxb = 1;
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:maxbSmallerThanminb");
            options = rmfield(options, 'minb');
            options = rmfield(options, 'maxb');


            options.mincon = -1;
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:minconNotValid");
            options.mincon = 2;
            options.maxcon = -1;
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:maxconNotValid");
            options.maxcon = 1;
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:maxconSmallerThanmincon");
            options = rmfield(options, 'mincon');
            options = rmfield(options, 'maxcon');


            options.excludelist = 1;
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:excludelistNotCellOfcharstr");
            options = rmfield(options, 'excludelist');


            options.problem_names = 1;
            testCase.verifyError(@() checkValidityProblemOptions(options, profile_options), "MATLAB:checkValidityProblemOptions:problem_namesNotCellOfcharstr");
        end

        function testOptionalSolarPlibName(testCase)

            mydir = fileparts(mfilename('fullpath'));
            original_dir = pwd;
            testCase.addTeardown(@() cd(original_dir));

            problem_libs_dir = fullfile(mydir, '../../problem_libs');
            solar_dir = fullfile(problem_libs_dir, 'solar');
            created_solar_dir = exist(solar_dir, 'dir') ~= 7;
            if created_solar_dir
                mkdir(solar_dir);
                testCase.addTeardown(@() rmdir(solar_dir));
            end

            cd(fullfile(mydir, '../../src/private'));

            options = struct('plibs', 'solar');
            profile_options = struct();
            options = checkValidityProblemOptions(options, profile_options);

            testCase.verifyEqual(options.plibs, {'solar'});
        end

    end

end
