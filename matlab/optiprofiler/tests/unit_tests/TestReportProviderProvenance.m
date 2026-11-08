classdef TestReportProviderProvenance < matlab.unittest.TestCase
% The report names the provider role and the lock's reviewed commit.
    methods (TestMethodSetup)
        function useRepositorySources(testCase)
            root = fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'..','..','src')));
        end
    end
    methods (Test)
        function bundledLibraryReportsLockedRevision(testCase)
            root = tempname; mkdir(root);
            old = pwd;
            % Teardowns run in reverse order: leave the directory before removing it.
            testCase.addTeardown(@() rmdir(root,'s'));
            testCase.addTeardown(@() cd(old));
            cd(root);
            options = struct('plibs',{{'s2mpj'}},'problem_names',{{'ROSENBR'}},'mindim',2,'maxdim',2, ...
                'n_jobs',1,'n_runs',1,'max_eval_factor',2,'max_tol_order',1,'run_plain',false, ...
                'score_only',true,'silent',true,'report_path',fullfile(root,'report.json'));
            benchmark({@TestReferencePoints.stay,@TestReferencePoints.zero},options);
            report = jsondecode(fileread(options.report_path));
            provider = firstProvider(report);
            lock = jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'..','..','..','..','matlab_problem_libraries.lock')));
            entries = lock.libraries;
            if isstruct(entries), entries = num2cell(entries); end
            s2mpj = entries{cellfun(@(entry) strcmp(entry.name,'s2mpj'),entries)};
            testCase.verifyEqual(provider.role,'bundled-default');
            testCase.verifyEqual(provider.locked_revision,s2mpj.commit);
        end
        function unlockedLibraryHasNoLockedRevision(testCase)
            root = tempname; mkdir(root);
            registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            oldpath = path; old = pwd;
            testCase.addTeardown(@() rmdir(root,'s'));
            testCase.addTeardown(@() path(oldpath));
            testCase.addTeardown(@() setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY',registry));
            testCase.addTeardown(@() cd(old));
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY',fullfile(root,'registry.mat'));
            library = fullfile(root,'library'); mkdir(library);
            TestReferencePoints.write(fullfile(library,'prov_probe_select.m'), sprintf( ...
                'function names=prov_probe_select(~)\n names={''P''};\nend\n'));
            TestReferencePoints.write(fullfile(library,'prov_probe_load.m'), sprintf( ...
                'function p=prov_probe_load(name)\n p=Problem(struct(''fun'',@(x)sum(x.^2),''x0'',2,''name'',name));\nend\n'));
            registerProblemLibrary(struct('name','prov_probe','root',library, ...
                'select_function','prov_probe_select','load_function','prov_probe_load'));
            cd(root);
            options = struct('plibs',{{'prov_probe'}},'n_jobs',1,'n_runs',1,'max_eval_factor',2, ...
                'max_tol_order',1,'run_plain',false,'score_only',true,'silent',true, ...
                'report_path',fullfile(root,'report.json'));
            benchmark({@TestReferencePoints.stay,@TestReferencePoints.zero},options);
            report = jsondecode(fileread(options.report_path));
            provider = firstProvider(report);
            testCase.verifyEqual(provider.role,'external');
            testCase.verifyEmpty(provider.locked_revision);
            testCase.verifySubstring(fileread(options.report_path),'"locked_revision":null');
        end
    end
end

function provider = firstProvider(report)
    problems = report.problems;
    if iscell(problems), provider = problems{1}.provider; else, provider = problems(1).provider; end
end
